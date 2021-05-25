#include <iostream>
#include "include/global.h"
#include "include/DNA_encoder.h"


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

    //


    return 0;
}
