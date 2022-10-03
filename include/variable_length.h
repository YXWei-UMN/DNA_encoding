#ifndef __VARIABLE_LENGTH__H__
#define __VARIABLE_LENGTH__H__

#include<unordered_map>
#include<vector>
#include<set>
#include<string>
#include<unordered_map>
#include<unordered_set>
#include<tuple>
#include "tool.h"

using namespace std;


// Format assumptions:
// 1. primer name is "primer[id]", e.g., "primer233"
// 2. the `path` is a directory containing all the blast result files. Each file is a single long strand, and we assueme it is cut to 200nt strands
// 3. both the primer and strand id are counted from 0, and the max id does not exceed uint range
// 4. For each primer, there is a line beginning with "# Query:", followed by the primer name. This is used to count the number of total primers.
// 5. n_strand only count the strand that appears in the blast result(i.e., strands with collisions)
// 6. If a collision is at [start, end], we cut it at the nearest "multiple of 10 position" with (start+end)/2, i.e., round((start+end)/2)*0.1)*10

class VariableLength {
public:
    VariableLength(string path);
    ~VariableLength() {}

    void Cut();
    void PrintStatistics(int total_collided_primers);
    unordered_set<PrimerID> get_recovered_primers() const{ return recovered_primers;}
    unordered_set<PrimerID> get_discarded_primers() const{ return discarded_primers;}

private:
    long total_collision_num=0;
    vector<bool> blind_spot;
    bool IsBlindSpot(unsigned int distance);
    void ReadCollisions(string path); // read the blast collision result of a single file
    unordered_map<PrimerID, int> n_primer_collisions; // mapping from <primer id> to <the number of collisions>

    vector<string> all_files;

    // data structure for each single file
    unordered_map<StrandID, string> strand_id2name;
    unordered_map<string, StrandID> strand_name2id;
    unordered_map<PrimerID, string> primer_id2name;
    unordered_map<string, PrimerID> primer_name2id;
    vector<Collision> collision_linear_order;

    unordered_set<string> all_primers; // all primer names(w/o collisions)
    unordered_map<PrimerID, unordered_set<PrimerID> > primer_confilct_list; // conflict between primers (primers that cannot exist together)
    unordered_map<PrimerID, int > primer_collision_num_;
    unordered_map<PrimerID, double > primer_capacity_;

    vector<pair<int, PrimerID>> primer_process_order;
    unordered_set<PrimerID> discarded_primers;
    unordered_set<PrimerID> recovered_primers;

    int n_primer = 0;
    int n_strand = 0;

};
#endif // __VARIABLE_LENGTH__H__