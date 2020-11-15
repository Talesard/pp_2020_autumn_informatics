// Copyright 2020 Vagina Valeria
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include "./max_vector_value.h"

TEST(Parallel_maxVectorValue, TestAverageSize) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int sizeVector = 100;
  std::vector<int> globalVector(sizeVector);

  if (rank == 0) {
    globalVector = getRandomVector(sizeVector);
  }

  int globalMax = maxVectorValueParallel(globalVector);

  if (rank == 0) {
    int seqentialMax = maxVectorValue(globalVector);
    ASSERT_EQ(seqentialMax, globalMax);
  }
}

TEST(Parallel_maxVectorValue, TestSmallSize) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int sizeVector = 5;
  std::vector<int> globalVector(sizeVector);

  if (rank == 0) {
    globalVector = getRandomVector(sizeVector);
  }

  int globalMax = maxVectorValueParallel(globalVector);

  if (rank == 0) {
    int seqentialMax = maxVectorValue(globalVector);
    ASSERT_EQ(seqentialMax, globalMax);
  }
}

TEST(Parallel_maxVectorValue, TestLargeSize) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int sizeVector = 1000;
  std::vector<int> globalVector(sizeVector);

  if (rank == 0) {
    globalVector = getRandomVector(sizeVector);
  }

  int globalMax = maxVectorValueParallel(globalVector);

  if (rank == 0) {
    int seqentialMax = maxVectorValue(globalVector);
    ASSERT_EQ(seqentialMax, globalMax);
  }
}

TEST(Parallel_maxVectorValue, TestIncorrectSize) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int sizeVector = 0;
  std::vector<int> globalVector(sizeVector);

  if (rank == 0) {
    globalVector = getRandomVector(sizeVector);
  }

  ASSERT_ANY_THROW(maxVectorValueParallel(globalVector));
}

TEST(Parallel_maxVectorValue, TestCreateVector) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<int> globalVector = {12, 34, 154, 86, 37, 12};

  int globalMax = maxVectorValueParallel(globalVector);

  if (rank == 0) {
    ASSERT_EQ(154, globalMax);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
  ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
  listeners.Release(listeners.default_result_printer());
  listeners.Release(listeners.default_xml_generator());
  listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
  return RUN_ALL_TESTS();
}
