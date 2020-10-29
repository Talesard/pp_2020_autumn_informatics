// Copyright 2020 Alibekov Murad
#ifndef MODULES_TASK_1_ALIBEKOV_M_COUNT_SENTENCES_COUNT_SENTENCES_H_
#define MODULES_TASK_1_ALIBEKOV_M_COUNT_SENTENCES_COUNT_SENTENCES_H_

#include <string>

constexpr char alphabet[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz .!,?";

std::string getRandomString(int sz);
int getSentencesCountSequential(std::string str);
int getSentencesCountParallel(std::string global_str, int count_chars);

#endif  // MODULES_TASK_1_ALIBEKOV_M_COUNT_SENTENCES_COUNT_SENTENCES_H_
