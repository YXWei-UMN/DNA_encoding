//
// Created by wyx on 19-6-20.
//

#ifndef CONFIGURABLE_DEDUP_GLOBAL_H
#define CONFIGURABLE_DEDUP_GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <ctime>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <list>
#include <sstream>
#include <time.h>
#include <chrono>
#include <deque>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <bitset>
#include <utility>
#include <math.h>
#include <inttypes.h>
#include <set>
#include <unordered_set>
#include <queue>
#include <cstdlib>
#include <ctime>
using namespace std;

#define PRIMER_CAPACITY 736

extern string g_data_path;
extern string g_payload_path;

extern string g_blast_result_path_1;
extern string g_blast_result_path_2;
extern string g_blast_result_path_3;
extern string g_blast_result_path_4;

extern string g_blast_result_path_varlen;


extern bool g_if_chunk;
extern bool g_if_pre_stranding;
extern bool g_if_ECC;
extern bool g_if_mapping;
extern bool g_if_randomization;
extern int g_encoding_scheme;
extern int g_program;
extern int g_dedup;
extern long g_payload_size;
extern long g_chunk_size;
extern long g_strand_length;
extern int g_swap_granularity;
extern int g_num_bit_per_triplet;
int Parse(string cfgfile);

typedef std::uint64_t hash_t;
constexpr hash_t prime = 0x100000001B3ull;
constexpr hash_t basis = 0xCBF29CE484222325ull;
constexpr hash_t hash_(char const* str, hash_t last_value = basis)
{
    return *str ? hash_(str+1, (*str ^ last_value) * prime) : last_value;
}
#endif //CONFIGURABLE_DEDUP_GLOBAL_H