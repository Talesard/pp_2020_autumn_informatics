// Copyright 2020 Vagina Valeria
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <stdexcept>
#include "./max_vector_value.h"

std::vector<int> getRandomVector(int size) {
  std::mt19937 gen;
  gen.seed(static_cast<unsigned int>(time(0)));
  std::vector<int> vec(size);
  for (int i = 0; i < size; i++) {
    vec[i] = gen() % 100;
  }

  return vec;
}

int maxVectorValue(std::vector<int> vec) {
  int size = vec.size();
  if (size <= 0) {
    throw std::runtime_error("Vector size must be positive");
  }

  return *std::max_element(vec.begin(), vec.end());
}

int maxVectorValueParallel(std::vector<int> globalVector) {
  int sizeVector = globalVector.size();
  if (sizeVector <= 0) {
    throw std::runtime_error("Vector size must be positive");
  }

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int delta = sizeVector / size;
  const int remainder = sizeVector % size;

  int localMax = 0, globalMax = 0;
  std::vector<int> localVector;

  if ((rank == 0) && (delta != 0)) {
    for (int i = 1; i < size; i++) {
          MPI_Send(&globalVector[0] + i * delta + remainder, delta, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
  }

  if (rank == 0) {
      if (delta != 0) {
        localVector = std::vector<int>(globalVector.begin(), globalVector.begin() + delta + remainder);
        localMax = maxVectorValue(localVector);
      } else {
        localMax = maxVectorValue(globalVector);
      }

  } else if (delta != 0) {
    localVector.resize(delta);
    MPI_Status status;
    MPI_Recv(&localVector[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    localMax = maxVectorValue(localVector);
  }

  MPI_Reduce(&localMax, &globalMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  return globalMax;
}
