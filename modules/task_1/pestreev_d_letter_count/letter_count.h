// Copyright 2020 Pestreev Daniil
#ifndef MODULES_TASK_1_PESTREEV_D_LETTER_COUNT_LETTER_COUNT_H_
#define MODULES_TASK_1_PESTREEV_D_LETTER_COUNT_LETTER_COUNT_H_
#include <string>

std::string getRandomString(int  size);
int64_t LetterInLineSequential(std::string str);
int64_t LetterInLineParallel(std::string global_str, int size_str);

#endif  // MODULES_TASK_1_PESTREEV_D_LETTER_COUNT_LETTER_COUNT_H_
