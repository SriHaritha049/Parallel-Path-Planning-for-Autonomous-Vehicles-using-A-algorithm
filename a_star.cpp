
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <stack>
#include <chrono>

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
    lat1 = lat1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;

    double a = sin(dLat/2) * sin(dLat/2) +
               cos(lat1) * cos(lat2) *
               sin(dLon/2) * sin(dLon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return R * c;
}

int main() {
    unordered_map<long long, Node> nodes;
    unordered_map<long long, vector<Edge>> graph;

    ifstream nodeFile("nodes_1.csv");
    string line;
    getline(nodeFile, line);
    while (getline(nodeFile, line)) {
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

    ifstream edgeFile("edges_1.csv");
    getline(edgeFile, line);
    while (getline(edgeFile, line)) {
        stringstream ss(line);
        string u_str, v_str, len_str;
        getline(ss, u_str, ',');
        getline(ss, v_str, ',');
        getline(ss, len_str, ',');
        long long u = stoll(u_str);
        long long v = stoll(v_str);
        double length = stod(len_str);
        graph[u].push_back({v, length});
        graph[v].push_back({u, length});
    }

    long long start = 195386940;
    long long goal  = 1616221113;

    if (nodes.count(start) == 0 || nodes.count(goal) == 0) {
        cerr << "Start or goal node not found!" << endl;
        return 1;
    }

    unordered_map<long long, double> gScore, fScore;
    unordered_map<long long, long long> cameFrom;
    unordered_set<long long> visited;

    auto cmp = [&](long long a, long long b) {
        return fScore[a] > fScore[b];
    };

    priority_queue<long long, vector<long long>, decltype(cmp)> openList(cmp);

    gScore[start] = 0;
    fScore[start] = haversine(nodes[start].lat, nodes[start].lon, nodes[goal].lat, nodes[goal].lon);
    openList.push(start);

    auto t1 = chrono::high_resolution_clock::now();

    while (!openList.empty()) {
        long long current = openList.top();
        openList.pop();

        if (current == goal) {
            auto t2 = chrono::high_resolution_clock::now();
            double duration = chrono::duration<double>(t2 - t1).count();
            cout << "✅ Path found!" << endl;
            stack<long long> path;
            while (cameFrom.count(current)) {
                path.push(current);
                current = cameFrom[current];
            }
            path.push(start);
            while (!path.empty()) {
                cout << path.top() << " ";
                path.pop();
            }
            cout << endl;
            cout << "⏱️  Sequential A* Time: " << duration << " seconds" << endl;
            return 0;
        }

        visited.insert(current);

        for (Edge& e : graph[current]) {
            long long neighbor = e.target;
            if (visited.count(neighbor)) continue;

            double tentative_g = gScore[current] + e.length;

            if (!gScore.count(neighbor) || tentative_g < gScore[neighbor]) {
                cameFrom[neighbor] = current;
                gScore[neighbor] = tentative_g;
                fScore[neighbor] = tentative_g + haversine(
                    nodes[neighbor].lat, nodes[neighbor].lon,
                    nodes[goal].lat, nodes[goal].lon
                );
                openList.push(neighbor);
            }
        }
    }

    cout << "❌ No path found." << endl;
    return 0;
}
