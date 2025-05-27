# Parallel Path Planning for Autonomous Vehicles using A* Algorithm

This project implements parallel versions of the A* path planning algorithm (OpenMP, MPI, and CUDA) using real-world road network data from OpenStreetMap.

## How to Run

1. **Extract graph data** from OSM:
   ```bash
   python ExtractdataOSM.py
   ```

2. **Run pathfinding**:
   ```bash
   python LargePathFinder.py
   ```

## Files
- `ExtractdataOSM.py`: Extracts nodes and edges from OSM into CSV files.
- `LargePathFinder.py`: Runs A* algorithm variants using the preprocessed graph.

## Parallel Implementations
- OpenMP A*
- MPI A*
- CUDA A*

Each implementation uses the same graph format for consistency and performance comparison.
