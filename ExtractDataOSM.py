import osmium
import math
import csv

# Bigger Bounding Box (Raleigh + Durham + Cary + Chapel Hill)
BBOX = (35.5, -79.0, 36.2, -78.4)  # (south, west, north, east)

pbf_file = "north-carolina-latest.osm.pbf"
nodes_csv = "nodes_large.csv"
edges_csv = "edges_large.csv"

def haversine(lat1, lon1, lat2, lon2):
    R = 6371000  # Earth radius in meters
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dlambda/2)**2
    return R * (2 * math.atan2(math.sqrt(a), math.sqrt(1-a)))

# -- Load and Save nodes first
class NodeCollector(osmium.SimpleHandler):
    def __init__(self):
        super(NodeCollector, self).__init__()
        self.nodes = {}
        self.counter = 0

    def node(self, n):
        self.counter += 1
        if self.counter % 100000 == 0:
            print(f"🔄 Processed {self.counter} nodes...")
        if (BBOX[0] <= n.location.lat <= BBOX[2]) and (BBOX[1] <= n.location.lon <= BBOX[3]):
            self.nodes[n.id] = (n.location.lat, n.location.lon)

print("⏳ Collecting nodes inside bounding box...")
node_handler = NodeCollector()
node_handler.apply_file(pbf_file, locations=True)
print(f"✅ Nodes collected: {len(node_handler.nodes)}")

# Save nodes
with open(nodes_csv, "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["node_id", "latitude", "longitude"])
    for nid, (lat, lon) in node_handler.nodes.items():
        writer.writerow([nid, lat, lon])
print("💾 Saved nodes_large.csv")

# -- Now find edges (but write immediately without storing all)

class EdgeCollector(osmium.SimpleHandler):
    def __init__(self, valid_nodes):
        super(EdgeCollector, self).__init__()
        self.valid_nodes = valid_nodes
        self.counter = 0
        self.edge_writer = csv.writer(open(edges_csv, "w", newline=''))
        self.edge_writer.writerow(["start_node", "end_node", "length"])

    def way(self, w):
        if "highway" in w.tags and len(w.nodes) >= 2:
            for i in range(len(w.nodes) - 1):
                n1 = w.nodes[i].ref
                n2 = w.nodes[i+1].ref
                if n1 in self.valid_nodes and n2 in self.valid_nodes:
                    lat1, lon1 = self.valid_nodes[n1]
                    lat2, lon2 = self.valid_nodes[n2]
                    dist = haversine(lat1, lon1, lat2, lon2)
                    self.edge_writer.writerow([n1, n2, round(dist, 2)])

print("⏳ Collecting edges...")
edge_handler = EdgeCollector(node_handler.nodes)
edge_handler.apply_file(pbf_file, locations=True)
print("📄 Saved edges_large.csv")
