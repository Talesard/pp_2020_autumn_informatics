// Copyright 2020 Kumbrasev Mark
#ifndef MODULES_TASK_1_KUMBRASEV_M_DIFF_CHAR_COUNTER_DIFF_CHAR_COUNTER_H_
#define MODULES_TASK_1_KUMBRASEV_M_DIFF_CHAR_COUNTER_DIFF_CHAR_COUNTER_H_

#include <utility>
#include <string>
#include <random>
#include <ctime>
#include <cstring>

template <int policy>
struct execute_policy {};

using seq_policy = execute_policy<0>;
using par_policy = execute_policy<1>;

std::string get_random_string();
std::size_t difference_count(seq_policy, const std::string& str_lhs, const std::string& str_rhs);
std::size_t difference_count(par_policy, const std::string str_lhs, const std::string str_rhs);

#endif  // MODULES_TASK_1_KUMBRASEV_M_DIFF_CHAR_COUNTER_DIFF_CHAR_COUNTER_H_
