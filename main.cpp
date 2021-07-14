#include <iostream>
#include "include/global.h"
#include "include/DNA_encoder.h"
#include "rs.hpp"
#define ECC_LENGTH 8

int main(int argc, char** argv) {
    /*if (argc != 2) {
        cerr << "argc must be 2" << endl;
        return -1;
    }*/

    string cfgfile = argv[1];

    if (Parse(cfgfile)) {
        cerr<< "parse config file " << cfgfile << " failed!\n";
        return -1;
    }

    DNA_encoder dnaEncoder;
    
    // test RS encoding
    // string result = dnaEncoder.ReedSolomon_encoding("Since synthesis and sequencing of very long DNA strands is technically impeded, , data must be stored on several short DNA segments, which cannot be arranged geometrically.");
    // cout << "result = " << result << endl;

    return 0;
}
