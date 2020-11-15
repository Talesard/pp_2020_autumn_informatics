// Copyright 2020 Galindo Buitrago Julio Cesar

#include<mpi.h>
#include<vector>
#include <ctime>
#include<random>
#include<iostream>

std::vector<int> VectorAleatorio(int n) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vectorial(n);
    for (int i = 0; i < n; i++) {
        vectorial[i] = gen()%100;
    }
    return vectorial;
}


double sigma(std::vector<int> vectorial) {
    int size = vectorial.size();
    double AdicionInterna = 0;
    for (int i = 0; i < size; i++) {
        AdicionInterna += vectorial[i];
    }
    return AdicionInterna;
}


double IntermedioSigma(std::vector<int> vectorial) {
    int size = vectorial.size();
    double AdicionInterna = 0;
    for (int i = 0; i < size; i++) {
        AdicionInterna += vectorial[i];
    }
    double PromedioInterno = AdicionInterna / size;
    return PromedioInterno;
}


double VectorMedia(std::vector<int> vectorial, int n) {
    int size, rango;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rango);
    int delta = n / size;
    int s = n - delta * size;
    if (rango == 0) {
        for (int proc = 1; proc < size; proc++) {
            MPI_Send(&vectorial[0] + proc * delta + s, delta, MPI_INT, proc, 0, MPI_COMM_WORLD);
        }
    }

    std::vector<int> VectorialInterno;
    if (rango == 0) {
        VectorialInterno = std::vector<int>(delta + s);
    } else {
         VectorialInterno = std::vector<int>(delta);
    }
    if (rango == 0) {
         VectorialInterno = std::vector<int>(vectorial.begin(), vectorial.begin() + delta+s);
    } else {
        MPI_Status status;
        MPI_Recv(&VectorialInterno[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }


    double PromedioGeneral = 0;
    double AdicionInterna = sigma(VectorialInterno);

    MPI_Reduce(&AdicionInterna, &PromedioGeneral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    PromedioGeneral /= n;
    return PromedioGeneral;
}






