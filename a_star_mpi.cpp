#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <cassert>
#include <cstddef>
using namespace std;

struct Node {
    long long id;
    double lat, lon;
};

struct Edge {
    long long target;
    double length;
};

double haversine(double lat1, double lon1, double lat2, double lon2) {
    const double R = 6371000.0;
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    lat1 *= M_PI / 180.0;
    lat2 *= M_PI / 180.0;
    double a = sin(dLat/2) * sin(dLat/2) +
               cos(lat1) * cos(lat2) *
               sin(dLon/2) * sin(dLon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return R * c;
}

// Structure to represent an update (for a neighbor that belongs to a remote rank).
struct Update {
    long long neighbor;
    long long parent;
    double tentative_g;
};

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    double time_start = MPI_Wtime();
    MPI_Status status;
    int flag = 0;

    // Load graph data.
    unordered_map<long long, Node> nodes;
    unordered_map<long long, vector<Edge>> graph;

    ifstream nodeFile("nodes_large.csv");
    string line;
    getline(nodeFile, line); // Skip header
    while(getline(nodeFile, line)){
        stringstream ss(line);
        string id_str, lat_str, lon_str;
        getline(ss, id_str, ',');
        getline(ss, lat_str, ',');
        getline(ss, lon_str, ',');
        long long id = stoll(id_str);
        double lat = stod(lat_str);
        double lon = stod(lon_str);
        nodes[id] = {id, lat, lon};
    }
    nodeFile.close();

    ifstream edgeFile("edges_large.csv");
    getline(edgeFile, line); // Skip header
    while(getline(edgeFile, line)){
        stringstream ss(line);
        string u_str, v_str, len_str;
        getline(ss, u_str, ',');
        getline(ss, v_str, ',');
        getline(ss, len_str, ',');
        long long u = stoll(u_str);
        long long v = stoll(v_str);
        double len = stod(len_str);
        graph[u].push_back({v, len});
        graph[v].push_back({u, len});
    }
    edgeFile.close();

    // Simple ownership: node belongs to process if (id % world_size == rank)
    auto owns = [&](long long id) { return (id % world_size) == rank; };

    // Set start and goal (either from command-line or defaults).
    long long start, goal;
    if (argc >= 3) {
        start = stoll(argv[1]);
        goal  = stoll(argv[2]);
    } else {
        start = 1563209699;   // Provided start
        goal  = 7882777921;   // Provided goal
    }

    // A* data structures.
    unordered_map<long long, double> gScore;
    unordered_map<long long, long long> cameFrom;
    // The open list: priority queue storing (f, node).
    priority_queue<pair<double, long long>, vector<pair<double, long long>>, greater<>> open;

    if (owns(start)) {
        gScore[start] = 0.0;
        double h = haversine(nodes[start].lat, nodes[start].lon, nodes[goal].lat, nodes[goal].lon);
        open.push({h, start});
        if(rank==0)
            cout << "🚀 [Rank " << rank << "] Starting A* from " << start << endl;
    }

    // Create an MPI_Datatype for the Update structure.
    MPI_Datatype MPI_Update;
    int block_lengths[3] = {1, 1, 1};
    MPI_Aint offsets[3];
    offsets[0] = offsetof(Update, neighbor);
    offsets[1] = offsetof(Update, parent);
    offsets[2] = offsetof(Update, tentative_g);
    MPI_Datatype types[3] = {MPI_LONG_LONG, MPI_LONG_LONG, MPI_DOUBLE};
    MPI_Type_create_struct(3, block_lengths, offsets, types, &MPI_Update);
    MPI_Type_commit(&MPI_Update);

    // Main A* loop.
    bool progress_made = false;
    while (true) {
        progress_made = false;

        // --- Drain incoming batched messages ---
        while (true) {
            MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
            if (!flag)
                break;
            int count;
            MPI_Recv(&count, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vector<Update> updates(count);
            MPI_Recv(updates.data(), count, MPI_Update, status.MPI_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (auto &upd: updates) {
                // Optionally: check and avoid two-cycles.
                if (cameFrom.count(upd.parent) && cameFrom[upd.parent] == upd.neighbor) {
                    // Skip this update.
                } else if (!gScore.count(upd.neighbor) || upd.tentative_g < gScore[upd.neighbor]) {
                    gScore[upd.neighbor] = upd.tentative_g;
                    cameFrom[upd.neighbor] = upd.parent;
                    double f = upd.tentative_g + haversine(nodes[upd.neighbor].lat, nodes[upd.neighbor].lon, nodes[goal].lat, nodes[goal].lon);
                    open.push({f, upd.neighbor});
                }
            }
            progress_made = true;
        }

        // --- Process local open list ---
        while (!open.empty()) {
            auto [f, current] = open.top();
            open.pop();

            // If goal is reached, break out (we continue until termination detection below).
            if (current == goal) {
                // Optionally, you can break from the loop, but here we let the processing continue.
                // cout << "🎯 [Rank " << rank << "] Expanded goal node " << goal << " (but continuing search)." << endl;
            }

            // Process neighbors of 'current'.
            vector<Update> batch; // Batch for remote updates.
            for (Edge &e : graph[current]) {
                long long neighbor = e.target;
                double tentative_g = gScore[current] + e.length;
                // Avoid a two-cycle: if current's parent is neighbor, skip.
                if (cameFrom.count(current) && cameFrom[current] == neighbor)
                    continue;
                if (!gScore.count(neighbor) || tentative_g < gScore[neighbor]) {
                    gScore[neighbor] = tentative_g;
                    cameFrom[neighbor] = current;
                    double f_val = tentative_g + haversine(nodes[neighbor].lat, nodes[neighbor].lon, nodes[goal].lat, nodes[goal].lon);
                    if (owns(neighbor)) {
                        open.push({f_val, neighbor});
                    } else {
                        // Accumulate update for remote process.
                        batch.push_back({neighbor, current, tentative_g});
                    }
                }
            }
            progress_made = true;

            // Group batched updates by target rank.
            unordered_map<int, vector<Update>> outgoing;
            for (auto &upd : batch) {
                int target_rank = upd.neighbor % world_size;
                outgoing[target_rank].push_back(upd);
            }
            // Send each batch.
            for (auto &entry : outgoing) {
                int trg = entry.first;
                int num = entry.second.size();
                MPI_Send(&num, 1, MPI_INT, trg, 0, MPI_COMM_WORLD);
                MPI_Send(entry.second.data(), num, MPI_Update, trg, 1, MPI_COMM_WORLD);
            }
        }

        // --- Global termination detection ---
        int local_idle = (open.empty() && !progress_made) ? 1 : 0;
        int all_idle = 0;
        MPI_Allreduce(&local_idle, &all_idle, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        if (all_idle)
            break;
    }

    double time_end = MPI_Wtime();
    if (rank == 0) {
        cout << "\n⏱️ Total execution time: " << (time_end - time_start) << " seconds\n";
    }

    // --- Aggregation Phase for Global cameFrom ---
    int aggregator = goal % world_size;
    unordered_map<long long, long long> globalCameFrom;
    if (rank == aggregator) {
        globalCameFrom = cameFrom;
        for (int r = 0; r < world_size; r++) {
            if (r == aggregator)
                continue;
            int map_size = 0;
            MPI_Recv(&map_size, 1, MPI_INT, r, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            vector<long long> buf(2 * map_size);
            MPI_Recv(buf.data(), 2 * map_size, MPI_LONG_LONG, r, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < map_size; i++) {
                long long node = buf[2 * i];
                long long parent = buf[2 * i + 1];
                if (globalCameFrom.find(node) == globalCameFrom.end()) {
                    globalCameFrom[node] = parent;
                }
            }
        }

        // Reconstruct path using aggregated cameFrom.
        cout << "\n🕯️ [Rank " << rank << "] Reconstructing path to " << goal << " using aggregated data:" << endl;
        vector<long long> path;
        long long current = goal;
        int safety_counter = 0;
        while (current != start && globalCameFrom.count(current) && safety_counter < 10000) {
            path.push_back(current);
            current = globalCameFrom[current];
            safety_counter++;
        }
        if (current != start) {
            cout << "❌ Failed to reconstruct path: incomplete chain." << endl;
        } else {
            path.push_back(start);
            reverse(path.begin(), path.end());
            cout << "✅ Reconstructed Path: ";
            for (auto id : path)
                cout << id << " ";
            cout << "\nPath Length (gScore): " << gScore[goal] << " meters" << endl;
        }
    } else {
        int map_size = cameFrom.size();
        MPI_Send(&map_size, 1, MPI_INT, aggregator, 10, MPI_COMM_WORLD);
        vector<long long> buf;
        buf.reserve(2 * map_size);
        for (auto &entry : cameFrom) {
            buf.push_back(entry.first);
            buf.push_back(entry.second);
        }
        MPI_Send(buf.data(), buf.size(), MPI_LONG_LONG, aggregator, 11, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "✅ [Rank " << rank << "] Finalizing." << endl;
    MPI_Type_free(&MPI_Update);
    MPI_Finalize();
    return 0;
}