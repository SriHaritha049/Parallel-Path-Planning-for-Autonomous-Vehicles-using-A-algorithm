# Compiler settings
CXX = g++
NVCC = nvcc
MPICXX = mpic++

CXXFLAGS = -O2 -std=c++11
OMPFLAGS = -fopenmp

# Source files (your filenames)
SEQ_SRC = a_star.cpp
OMP_SRC = a_star_openmp.cpp
MPI_SRC = a_star_mpi.cpp
CUDA_SRC = a_star_cuda.cu

# Executables
SEQ_EXE = a_star
OMP_EXE = a_star_openmp
MPI_EXE = a_star_mpi
CUDA_EXE = a_star_cuda

.PHONY: all clean run_seq run_omp run_mpi run_cuda

all: $(SEQ_EXE) $(OMP_EXE) $(MPI_EXE) $(CUDA_EXE)

$(SEQ_EXE): $(SEQ_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OMP_EXE): $(OMP_SRC)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $^

$(MPI_EXE): $(MPI_SRC)
	$(MPICXX) $(CXXFLAGS) -o $@ $^

$(CUDA_EXE): $(CUDA_SRC)
	$(NVCC) -O2 -o $@ $^

run_seq: $(SEQ_EXE)
	./$(SEQ_EXE)

run_omp: $(OMP_EXE)
	./$(OMP_EXE)

run_mpi: $(MPI_EXE)
	mpirun -np 4 ./$(MPI_EXE)

run_cuda: $(CUDA_EXE)
	./$(CUDA_EXE)

clean:
	rm -f $(SEQ_EXE) $(OMP_EXE) $(MPI_EXE) $(CUDA_EXE)
