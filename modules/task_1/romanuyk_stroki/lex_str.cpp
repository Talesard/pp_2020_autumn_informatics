// Copyright 2020 Romanuyk Sergey
#include "./../../modules/task_1/romanuyk_stroki/lex_str.h"
#include <mpi.h>
#include <vector>
#include <ctime>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <string>
int getComparison(const char* str1, const char* str2, int len) {
for (int i = 0; i < len; i++) {
if (str1[i] > str2[i]) return 1;
else if (str1[i] < str2[i]) return -1;
}
return 0;
}
int getComparisonLIN(const char* str1, const char* str2) {
int len1 = strlen(str1);
int len2 = strlen(str2);
if (len1 < len2) {
return -1;
} else if (len1 > len2) {
return 1;
}
for (int i = 0; i < len1; i++) {
if (str1[i] > str2[i]) return 1;
else if (str1[i] < str2[i]) return -1;
}
return 0;
}
int getParallelComprasion(const char* str1, const char* str2) {
int rank, size;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
const unsigned int len1 = (unsigned int)strlen(str1);
const unsigned int len2 = (unsigned int)strlen(str2);
std::vector<int>check(size, 0);
int len_size = 0;
if (len1 < len2) {
return -1;
} else if (len1 > len2) {
return 1;
}
len_size = strlen(str1);
int tmp = len_size % size;
int delta = len_size / size;
if (rank == 0) {
for (int proc = 1; proc < size; proc++) {
MPI_Send(&str1[0] + tmp + proc * delta, delta, MPI_CHAR, proc, 0, MPI_COMM_WORLD);
MPI_Send(&str2[0] + tmp + proc * delta, delta, MPI_CHAR, proc, 1, MPI_COMM_WORLD);
}
}
char* str1_local;
char* str2_local;
int size_local_len = 0;
if (rank == 0) {
size_local_len = tmp + delta;
str1_local = new char[size_local_len];
str2_local = new char[size_local_len];
for (int i = 0; i < tmp + delta; i++) {
str1_local[i] = str1[i];
str2_local[i] = str2[i];
}
check[0] = getComparison(str1_local, str2_local, size_local_len);
} else {
size_local_len = delta;
str1_local = new char[size_local_len];
str2_local = new char[size_local_len];
MPI_Status status1;
MPI_Status status2;
int loc_resk = 0;
MPI_Recv(&str1_local[0], delta, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status1);
MPI_Recv(&str2_local[0], delta, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &status2);
loc_resk = getComparison(str1_local, str2_local, size_local_len);
MPI_Send(&loc_resk, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
}
if (rank == 0) {
for (int j = 1; j < size; j++) {
MPI_Status status;
MPI_Recv(&check[j], 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
}
}
delete[] str1_local;
delete[] str2_local;
for (int i = 0; i < size; i++) {
if (check[i] == 1) {
return 1;
}
if (check[i] == -1) {
return -1;
}
}
return 0;
}
