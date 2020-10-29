// Copyright 2020 Volkova Anastasia
#ifndef MODULES_TASK_1_VOLKOVA_A_SENTS_IN_LINE_SENTS_MPI_H_
#define MODULES_TASK_1_VOLKOVA_A_SENTS_IN_LINE_SENTS_MPI_H_

#include <vector>

std::vector<char> getRandomLine(int  size);
int getParallelSumSentences(std::vector<char> global_line, int size_line);
int getSequentialSumSentences(std::vector<char> line);


#endif  // MODULES_TASK_1_VOLKOVA_A_SENTS_IN_LINE_SENTS_MPI_H_
