// Copyright 2020 Kulandin Denis
#include <mpi.h>
#include <ctime>
#include <iostream>
#include <random>
#include <utility>
#include <vector>
#include "../../../modules/task_3/kulandin_d_quick_sort_simple/quick_sort.h"

std::mt19937 gen(time(0));

std::vector<int> randomVector(int sz) {
    std::vector<int> ans(sz);
    for (int i = 0; i < sz; ++i) {
        ans[i] = gen();
    }
    return ans;
}

void division(std::vector<int>* a, int l, int r, int* t) {
    int mid = (*a)[l + gen() % (r - l + 1)];
    while (l <= r) {
        while ((*a)[l] < mid) l++;
        while ((*a)[r] > mid) r--;
        if (l >= r) break;
        std::swap((*a)[l++], (*a)[r--]);
    }
    *t = r;
}

void quickSort(std::vector<int>* a, int l, int r) {
    // std::cout << l << ' ' << r << '\n';
    if (l >= r) return;
    int t = 0;
    division(a, l, r, &t);
    quickSort(a, l, t);
    quickSort(a, t + 1, r);
}

std::vector<int> parallelQuickSort(const std::vector<int> &a) {
    int length = a.size();
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    int first               = 0;
    int plength             = length;
    std::vector<int> part   = a;
    int leftProc            = rank;
    int rightProc           = size;

    if (rank != 0) {
        MPI_Recv(&first,        1,          MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&plength,      1,          MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        part = std::vector<int>(plength);
        MPI_Recv(&part[0],      plength,    MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&rightProc,    1,          MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
    }

    int l = 0;
    int r = plength - 1;
    // printf("initial size: rank = %d, l = %d, r = %d\n", rank, l, r);
    int m;
    int leftP = (rightProc + leftProc) / 2;
    int rightP = rightProc;
    while (r > l && leftP > leftProc) {
        division(&part, l, r, &m);
        int start   = 0;
        int count   = 0;
        if (m - l + 1 <= r - m) {
            start   = l;
            count   = m - l + 1;
            l       = m + 1;
        } else {
            start   = m + 1;
            count   = r - m;
            r       = m;
        }
        if (1) {
            // printf("from %d to %d, cnt = %d\n", rank, leftP, count);
            MPI_Send(&start,            1,      MPI_INT, leftP, 0, MPI_COMM_WORLD);
            MPI_Send(&count,            1,      MPI_INT, leftP, 0, MPI_COMM_WORLD);
            MPI_Send(&part[0] + start,  count,  MPI_INT, leftP, 0, MPI_COMM_WORLD);
            MPI_Send(&rightP,           1,      MPI_INT, leftP, 0, MPI_COMM_WORLD);
        }
        rightP = leftP;
        leftP = (rightP + leftProc) / 2;
    }
    // printf("%d %d %d\n", rank, l, r);
    quickSort(&part, l, r);
    // MPI_Barrier(MPI_COMM_WORLD);
    leftP = (rightProc + leftProc) / 2;
    while (leftP > leftProc) {
        int start;
        int count;
        MPI_Status tmp;

        MPI_Recv(&start,            1,      MPI_INT, leftP, 0, MPI_COMM_WORLD, &tmp);
        MPI_Recv(&count,            1,      MPI_INT, leftP, 0, MPI_COMM_WORLD, &tmp);
        MPI_Recv(&part[0] + start,  count,  MPI_INT, leftP, 0, MPI_COMM_WORLD, &tmp);
        // printf("recieve from %d to %d, cnt = %d\n", leftP, rank, count);
        leftP = (leftP + leftProc) / 2;
    }
    if (rank != 0) {
        MPI_Send(&first,    1,          MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        MPI_Send(&plength,  1,          MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
        MPI_Send(&part[0],  plength,    MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
    }

    return part;
}
