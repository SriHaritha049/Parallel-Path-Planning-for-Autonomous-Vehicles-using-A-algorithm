#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <cuda_runtime.h>
#include <stack>

using namespace std;

struct Node {
    long long id;
    double lat, lon;
};

struct Edge {
    long long target;
    double length;
};

__device__ double haversine(double lat1, double lon1, double lat2, double lon2) {
    const double R = 6371000.0;
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    lat1 *= M_PI / 180.0;
    lat2 *= M_PI / 180.0;
    double a = sin(dLat/2)*sin(dLat/2) + cos(lat1)*cos(lat2)*sin(dLon/2)*sin(dLon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return R * c;
}

__global__ void expand_kernel(
    Edge* edges,
    Node* nodes,
    int* edge_offsets,
    long long* open_current,
    int open_size_current,
    long long* open_next,
    int* open_size_next,
    double* gScore,
    long long* cameFrom,
    long long goal,
    bool* goal_found,
    int num_nodes
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= open_size_current || *goal_found) return;

    long long current = open_current[idx];

    int start = edge_offsets[current];
    int end = (current + 1 < num_nodes) ? edge_offsets[current + 1] : edge_offsets[num_nodes];

    for (int e = start; e < end; ++e) {
        long long neighbor = edges[e].target;
        double tentative_g = gScore[current] + edges[e].length;

        if (tentative_g < gScore[neighbor]) {
            gScore[neighbor] = tentative_g;
            cameFrom[neighbor] = current;

            int pos = atomicAdd(open_size_next, 1);
            open_next[pos] = neighbor;

            if (neighbor == goal) {
                *goal_found = true;
            }
        }
    }
}

int main() {
    unordered_map<long long, Node> nodes_map;
    unordered_map<long long, vector<Edge>> graph;

    ifstream nodeFile("nodes_large.csv");
    string line;
    getline(nodeFile, line); // skip header
    while (getline(nodeFile, line)) {
        stringstream ss(line);
        string id_str, lat_str, lon_str;
        getline(ss, id_str, ',');
        getline(ss, lat_str, ',');
        getline(ss, lon_str, ',');
        long long id = stoll(id_str);
        double lat = stod(lat_str);
        double lon = stod(lon_str);
        nodes_map[id] = {id, lat, lon};
    }

    ifstream edgeFile("edges_large.csv");
    getline(edgeFile, line); // skip header
    while (getline(edgeFile, line)) {
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

    // Compress nodes into arrays
    vector<Node> nodes_vec;
    unordered_map<long long, int> id_to_idx;
    int idx = 0;
    for (auto& [id, node] : nodes_map) {
        id_to_idx[id] = idx++;
        nodes_vec.push_back(node);
    }

    int num_nodes = nodes_vec.size();
    vector<Edge> edges_vec;
    vector<int> edge_offsets(num_nodes + 1, 0);

    idx = 0;
    for (auto& node : nodes_vec) {
        edge_offsets[idx] = edges_vec.size();
        for (auto& e : graph[node.id]) {
            edges_vec.push_back({id_to_idx[e.target], e.length});
        }
        idx++;
    }
    edge_offsets[num_nodes] = edges_vec.size();

    long long start_id = 195386940;
    long long goal_id = 1616221113;

    if (id_to_idx.find(start_id) == id_to_idx.end() || id_to_idx.find(goal_id) == id_to_idx.end()) {
        cerr << "Start or goal not found!" << endl;
        return 1;
    }

    int start = id_to_idx[start_id];
    int goal = id_to_idx[goal_id];

    Edge* d_edges;
    Node* d_nodes;
    int* d_edge_offsets;
    long long* d_open_current;
    long long* d_open_next;
    int* d_open_size_next;
    double* d_gScore;
    long long* d_cameFrom;
    bool* d_goal_found;

    cudaMalloc(&d_edges, edges_vec.size() * sizeof(Edge));
    cudaMalloc(&d_nodes, nodes_vec.size() * sizeof(Node));
    cudaMalloc(&d_edge_offsets, edge_offsets.size() * sizeof(int));
    cudaMalloc(&d_open_current, num_nodes * sizeof(long long));
    cudaMalloc(&d_open_next, num_nodes * sizeof(long long));
    cudaMalloc(&d_open_size_next, sizeof(int));
    cudaMalloc(&d_gScore, num_nodes * sizeof(double));
    cudaMalloc(&d_cameFrom, num_nodes * sizeof(long long));
    cudaMalloc(&d_goal_found, sizeof(bool));

    cudaMemcpy(d_edges, edges_vec.data(), edges_vec.size() * sizeof(Edge), cudaMemcpyHostToDevice);
    cudaMemcpy(d_nodes, nodes_vec.data(), nodes_vec.size() * sizeof(Node), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edge_offsets, edge_offsets.data(), edge_offsets.size() * sizeof(int), cudaMemcpyHostToDevice);

    double* h_gScore = new double[num_nodes];
    long long* h_cameFrom = new long long[num_nodes];
    fill(h_gScore, h_gScore + num_nodes, 1e9);
    fill(h_cameFrom, h_cameFrom + num_nodes, -1);
    h_gScore[start] = 0.0;

    cudaMemcpy(d_gScore, h_gScore, num_nodes * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cameFrom, h_cameFrom, num_nodes * sizeof(long long), cudaMemcpyHostToDevice);

    long long h_open_current[num_nodes];
    h_open_current[0] = start;
    int h_open_size_current = 1;
    cudaMemcpy(d_open_current, h_open_current, num_nodes * sizeof(long long), cudaMemcpyHostToDevice);

    bool h_goal_found = false;
    cudaMemcpy(d_goal_found, &h_goal_found, sizeof(bool), cudaMemcpyHostToDevice);

    int iterations = 0;

    // Timing
    cudaEvent_t start_event, stop_event;
    cudaEventCreate(&start_event);
    cudaEventCreate(&stop_event);
    cudaEventRecord(start_event);

    while (!h_goal_found && h_open_size_current > 0) {
        int zero = 0;
        cudaMemcpy(d_open_size_next, &zero, sizeof(int), cudaMemcpyHostToDevice);

        expand_kernel<<<(h_open_size_current + 255) / 256, 256>>>(
            d_edges, d_nodes, d_edge_offsets,
            d_open_current, h_open_size_current,
            d_open_next, d_open_size_next,
            d_gScore, d_cameFrom, goal, d_goal_found, num_nodes
        );
        cudaDeviceSynchronize();

        cudaMemcpy(&h_goal_found, d_goal_found, sizeof(bool), cudaMemcpyDeviceToHost);
        cudaMemcpy(&h_open_size_current, d_open_size_next, sizeof(int), cudaMemcpyDeviceToHost);

        swap(d_open_current, d_open_next);
        iterations++;
    }

    cudaEventRecord(stop_event);
    cudaEventSynchronize(stop_event);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start_event, stop_event);

    cudaMemcpy(h_cameFrom, d_cameFrom, num_nodes * sizeof(long long), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_gScore, d_gScore, num_nodes * sizeof(double), cudaMemcpyDeviceToHost);

    cout << "\n🎯 Goal reached after " << iterations << " kernel launches.\n";

    cout << "Path:\n";
    vector<int> path_indices;
    for (int at = goal; at != -1; at = h_cameFrom[at]) {
        path_indices.push_back(at);
    }
    reverse(path_indices.begin(), path_indices.end());
    for (auto idx : path_indices) {
        cout << nodes_vec[idx].id << " ";
    }
    cout << endl;


    cout << "\n📏 Path length: " << h_gScore[goal] << " meters\n";
    cout << "⏱️ CUDA A* Execution time: " << milliseconds / 1000.0 << " seconds\n";

    cudaFree(d_edges);
    cudaFree(d_nodes);
    cudaFree(d_edge_offsets);
    cudaFree(d_open_current);
    cudaFree(d_open_next);
    cudaFree(d_open_size_next);
    cudaFree(d_gScore);
    cudaFree(d_cameFrom);
    cudaFree(d_goal_found);

    delete[] h_gScore;
    delete[] h_cameFrom;

    cudaEventDestroy(start_event);
    cudaEventDestroy(stop_event);

    return 0;
}
