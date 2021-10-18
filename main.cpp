#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <set>
#include <ctime>

#include "include/global.h"
#include "include/DNA_encoder.h"
#include "rs.hpp"
#include "transform_selection.h"
#include "variable_length.h"
#include "var_transform.h"
#define ECC_LENGTH 8

int main(int argc, char** argv) {
    srand(time(NULL));
    clock_t time_s0, time_t0;
    time_s0 = clock();


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

    if (g_program==1){ // encoding
        DNA_encoder dnaEncoder;
    }
    else if (g_program == 2){ // transform selection
        string encoding_blast_result_path[4] = {
                g_blast_result_path_1,
                g_blast_result_path_2,
                g_blast_result_path_3,
                g_blast_result_path_4};

        TransformSelection selection(encoding_blast_result_path);
        selection.Select();
        selection.PrintStatistics();
    }
    else if (g_program == 3){ // variable length
        VariableLength varlen(g_blast_result_path_varlen);
        varlen.Cut();
    } else if (g_program == 4){ // baseline combination
        string encoding_blast_result_path[4] = {
                g_blast_result_path_1,
                g_blast_result_path_2,
                g_blast_result_path_3,
                g_blast_result_path_4};
        VarTransform var_transform(encoding_blast_result_path);
        var_transform.Run();
        var_transform.PrintStatistics();
    }

    time_t0 = clock();
    double test_time = (double)(time_t0 - time_s0) / CLOCKS_PER_SEC / 60.0;
    cout << "Time: " << test_time << "min" << endl;

    return 0;
}
