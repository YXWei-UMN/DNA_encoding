//
// Created by wyx on 19-6-20.
//

#include "../include/global.h"



string g_data_path;
string g_payload_path;
bool g_if_dedupe;
bool g_if_chunk;
bool g_if_mapping;
bool g_if_randomization;
bool g_if_pre_stranding;
bool g_if_ECC;
int g_encoding_scheme;
long g_payload_size;
long g_chunk_size;
long g_strand_length;
int g_swap_granularity;
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
            case hash_("swap_granularity"):
                g_swap_granularity = stoi(value);
                break;
            case hash_("strand_length"):
                g_strand_length = stol(value);
                break;
            case hash_("payload_size"):
                g_payload_size = stol(value);
                break;
            case hash_("data_path"):
                g_data_path = value;
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
            case hash_("if_mapping"):
                g_if_dedupe = (value=="true");
                break;
            case hash_("encoding_scheme"):
                g_encoding_scheme = stoi(value);
                break;
            default:
                cout<<"unknown cfg: "<<key<<endl;
                return -1;
        }
    }
    return 0;
}
