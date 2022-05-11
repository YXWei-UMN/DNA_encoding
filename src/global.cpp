//
// Created by wyx on 19-6-20.
//

#include "../include/global.h"


 string g_blast_result_path_1;
 string g_blast_result_path_2;
 string g_blast_result_path_3;
 string g_blast_result_path_4;

string g_blast_result_path_varlen;

string g_data_path;
string g_payload_path;
int g_program;
bool g_if_mapping;
bool g_if_randomization;
bool g_if_pre_stranding;
bool g_if_ECC;
int g_encoding_scheme;
int g_swap_granularity;
int g_dedup;
int g_num_bit_per_triplet;
int Parse(string cfgfile){
    ifstream filestream(cfgfile, ios_base::in);
    if (filestream.fail()) {
        cerr << "open cfgfile:" << cfgfile << " fails!\n";
        return -1;
    }
    string line;

    while(getline(filestream, line)) {
        if (line.size()<=1 || line[0]== '#')
            continue;

        stringstream ss(line);
        string key, value;
        getline(ss, key, ' ');
        getline(ss, value, ' ');

        switch(hash_(key.c_str())){
            case hash_("data_path"):
                g_data_path = value;
                break;
            case hash_("blast_result_path_1"):
                g_blast_result_path_1 = value;
                break;
            case hash_("blast_result_path_2"):
                g_blast_result_path_2 = value;
                break;
            case hash_("blast_result_path_3"):
                g_blast_result_path_3 = value;
                break;
            case hash_("blast_result_path_4"):
                g_blast_result_path_4 = value;
                break;
            case hash_("blast_result_path_varlen"):
                g_blast_result_path_varlen = value;
                break;
            case hash_("payload_path"):
                g_payload_path = value;
                break;
            case hash_("if_pre_stranding"):
                g_if_pre_stranding = (value=="true");
                break;
            case hash_("if_ECC"):
                g_if_ECC = (value=="true");
                break;
            case hash_("if_randomization"):
                g_if_randomization = (value=="true");
                break;
            case hash_("encoding_scheme"):
                g_encoding_scheme = stoi(value);
                break;
            case hash_("num_bit_per_triplet"):
                g_num_bit_per_triplet = stoi(value);
                break;
            case hash_("dedup"):
                g_dedup = stoi(value);
                break;
            case hash_("program"):
                g_program = stoi(value);
                break;
            default:
                cout<<"unknown cfg: "<<key<<endl;
                return -1;
        }
    }
    return 0;
}
