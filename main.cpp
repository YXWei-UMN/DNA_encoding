#include <iostream>
#include "include/global.h"
#include "include/DNA_encoder.h"
#include "rs.hpp"
#include "transform_selection.h"
#include "variable_length.h"
#define ECC_LENGTH 8

int main(int argc, char** argv) {
    ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);
    /*if (argc != 2) {
        cerr << "argc must be 2" << endl;
        return -1;
    }*/

    string cfgfile = argv[1];

    if (Parse(cfgfile)) {
        cerr<< "parse config file " << cfgfile << " failed!\n";
        return -1;
    }

    if (g_program==1){
        DNA_encoder dnaEncoder;
    }
    else if (g_program == 2){
        string encoding_blast_result_path[4] = {
                g_blast_result_path_1,
                g_blast_result_path_2,
                g_blast_result_path_3,
                g_blast_result_path_4};

        TransformSelection selection(encoding_blast_result_path);
        selection.Select();
        selection.PrintStatistics();
    }
    else if (g_program == 3){
        VariableLength varlen(g_blast_result_path_varlen);
        varlen.Cut();
    }



    return 0;
}
