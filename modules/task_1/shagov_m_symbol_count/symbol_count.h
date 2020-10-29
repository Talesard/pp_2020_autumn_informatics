// Copyright 2020 Shagov Maksim
#ifndef MODULES_TASK_1_SHAGOV_M_SYMBOL_COUNT_SYMBOL_COUNT_H_
#define MODULES_TASK_1_SHAGOV_M_SYMBOL_COUNT_SYMBOL_COUNT_H_
#include <string>

std::string createRandomString(int  sz);
int calculateSymbolsCountSequental(const std::string& string, char symbol);
int getParallelSymbolsCount(const std::string& global_string,
                            int count_size_vector, char symbol);
#endif  // MODULES_TASK_1_SHAGOV_M_SYMBOL_COUNT_SYMBOL_COUNT_H_
