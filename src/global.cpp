//
// Created by wyx on 19-6-20.
//

#include "../include/global.h"



string g_data_path;
string g_payload_path;
bool g_if_dedupe;
bool g_if_chunk;
bool g_base3_rotate_encoding;
long g_payload_size;
long g_chunk_size;
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
            case hash_("chunk_size"):
                g_chunk_size = stol(value);
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
            case hash_("if_dedupe"):
                g_if_dedupe = (value=="true");
                break;
            case hash_("if_chunk"):
                g_if_chunk = (value=="true");
                break;
            case hash_("base3_rotate_encoding"):
                g_base3_rotate_encoding = (value=="true");
                break;
            default:
                cout<<"unknown cfg: "<<key<<endl;
                return -1;
        }
    }
    return 0;
}
