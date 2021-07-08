//
// Created by eason on 2/22/21.
//

#ifndef DNA_ENCODING_DNA_ENCODER_H
#define DNA_ENCODING_DNA_ENCODER_H
#include "global.h"
#include <iostream>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <bitset>
#include <sstream>
#include <cstring>



class DNA_encoder {
public:
    DNA_encoder();
    ~DNA_encoder(){}

    //first translate to nt sequence and then stranding
    void encoding_stranding();
    // translate to nt sequence, without strand
    void encoding_no_strand();
    // translate to nt sequence in file level (1 file 1 strand)
    void encoding_file();

    // encode a payload sequence in base 3 with rotate manner
    string base3_rotate_encoding(string digital_data);
    void initial_rotating_encoding_table();

    string FEC_encoding(string digital_data);
    void initial_FEC_table();
    uint16_t CCITT16(char *ptr, int length);

    string direct_encoding(string digital_data);

    void randomize_XOR(string& digital_data);
    void listFiles(string baseDir, bool recursive);

    string pseudo_random_sequence_;
    vector<string> all_files_;


    vector<unordered_map<string,string>> rotating_encoding_table_;
    string FEC_table1;
    vector<vector<string>> FEC_table2;
    string last_bit = "C"; // used for the fist rotating encoding
};


#endif //DNA_ENCODING_DNA_ENCODER_H
