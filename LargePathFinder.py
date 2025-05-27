import csv
import heapq
from collections import defaultdict

# Step 1: Load Graph
graph = defaultdict(list)

with open("edges_large.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        u = int(row["start_node"])
        v = int(row["end_node"])
        d = float(row["length"])
        graph[u].append((v, d))
        graph[v].append((u, d))

# Step 2: Pick arbitrary start node
start = next(iter(graph))
dist = {start: 0}
prev = {}

# Step 3: Dijkstra to find farthest node
pq = [(0, start)]
while pq:
    cost, u = heapq.heappop(pq)
    for v, w in graph[u]:
        alt = cost + w
        if v not in dist or alt < dist[v]:
            dist[v] = alt
            prev[v] = u
            heapq.heappush(pq, (alt, v))

# Step 4: Find farthest node
goal = max(dist.items(), key=lambda x: x[1])[0]
path_length = dist[goal]

# Optional: Reconstruct the path
path = []
cur = goal
while cur != start:
    path.append(cur)
    cur = prev[cur]
path.append(start)
path.reverse()

print(f"✅ Longest connected path found:")
print(f"Start: {start}")
print(f"Goal:  {goal}")
print(f"Length: {path_length:.2f} meters")
print(f"Path length: {len(path)} nodes")
# import csv
# import networkx as nx

# # Paths to your partial files
# nodes_file = "nodes_large.csv"
# edges_file = "edges_large.csv"

# # Load nodes
# nodes = {}
# with open(nodes_file, 'r') as f:
#     reader = csv.DictReader(f)
#     for row in reader:
#         nodes[int(row['node_id'])] = (float(row['latitude']), float(row['longitude']))

# # Load edges
# edges = []
# with open(edges_file, 'r') as f:
#     reader = csv.DictReader(f)
#     for row in reader:
#         u = int(row['start_node'])
#         v = int(row['end_node'])
#         length = float(row['length'])
#         if u in nodes and v in nodes:
#             edges.append((u, v, length))

# print(f"✅ Loaded {len(nodes)} nodes and {len(edges)} edges")

# # Build graph
# G = nx.Graph()
# for u, v, w in edges:
#     G.add_edge(u, v, weight=w)

# print(f"📈 Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

# # Pick two farthest nodes using graph eccentricity
# try:
#     # (This is an approximate method because full eccentricity is expensive)
#     nodes_list = list(G.nodes)
#     if len(nodes_list) < 2:
#         print("❌ Not enough nodes.")
#     else:
#         # Pick two random nodes and run Dijkstra to find path
#         source = nodes_list[0]
#         target = nodes_list[-1]
#         print(f"🔍 Trying to find path from {source} to {target}...")

#         path = nx.dijkstra_path(G, source, target, weight="weight")
#         length = nx.dijkstra_path_length(G, source, target, weight="weight")

#         print(f"✅ Path found with {len(path)} nodes and total length {length:.2f} meters.")
#         print(path)

# except Exception as e:
#     print(f"❌ Failed to find path: {e}")
