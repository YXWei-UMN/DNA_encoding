#include <iostream>
#include "include/global.h"
#include "include/DNA_encoder.h"
#include "rs.hpp"
#include "transform_selection.h"
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

    DNA_encoder dnaEncoder;

    
    // test RS encoding
    // string result = dnaEncoder.ReedSolomon_encoding("Since synthesis and sequencing of very long DNA strands is technically impeded, , data must be stored on several short DNA segments, which cannot be arranged geometrically.");
    // cout << "result = " << result << endl;

    // test transformation selection baseline
    string encoding_blast_result_path[4] = {
        "/home/umhadmin/huibing/dna/data/test_payload_C/blast_test_payload_C_evalue50_200strand",
        "/home/umhadmin/huibing/dna/data/test_swap_1/blast_evalue50_test_swap_1_200strand",
        "/home/umhadmin/huibing/dna/data/test_swap_3/blast_evalue50_test_swap_3_200strand",
        "/home/umhadmin/huibing/dna/data/test_swap_5/blast_test_swap_5_evalue50_200strand"};
    
    // string encoding_blast_result_path[4] = {
    //     "/home/umhadmin/huibing/dna/data/selection_test/blast0",
    //     "/home/umhadmin/huibing/dna/data/selection_test/blast1",
    //     "/home/umhadmin/huibing/dna/data/selection_test/blast2",
    //     "/home/umhadmin/huibing/dna/data/selection_test/blast3"};
    
    TransformSelection selection(encoding_blast_result_path);
    selection.Select();
    selection.PrintStatistics();


    return 0;
}
