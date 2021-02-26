//
// Created by eason on 2/22/21.
//

#ifndef DNA_ENCODING_DNA_ENCODER_H
#define DNA_ENCODING_DNA_ENCODER_H
#include "global.h"
#include "rabin.h"
#include <iostream>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <bitset>

class DNA_encoder {
public:
    DNA_encoder();
    ~DNA_encoder(){}


    static bool isDir(string dir);
    void listFiles(string baseDir, bool recursive);
    // if apply chunking, may have dedupe inside it
    void chunking_encode();
    // encode a payload sequence in base 3 with rotate manner
    string base3_rotate_encoding(string digital_data);
    void initial_rotating_encoding_table();
    vector<string> all_files_;

    unordered_set<uint64_t> chunk_index_;
    unsigned int distinct_chunks_ = 0;
    unsigned int total_chunks_ = 0;
    size_t bytes=0;
    vector<unordered_map<string,string>> rotating_encoding_table_;
    string last_bit = "A"; // used for the fist rotating encoding
};


#endif //DNA_ENCODING_DNA_ENCODER_H
