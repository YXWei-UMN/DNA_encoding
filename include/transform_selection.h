#ifndef __TRANSFORM_SELECTION__H__
#define __TRANSFORM_SELECTION__H__

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <string>

#include "DNA_encoder.h"

using namespace std;


// Format assumptions:
// 1. primer name is "primer[id]", e.g., "primer233"
// 2. strand name is "payload[id], e.g., "payload12345"
// 3. both the primer and strand id are counted from 0, and the max id does not exceed uint range
// 4. For each primer, there is a line beginning with "# Query:", followed by the primer name. This is used to count the number of total primers.
// 5. n_strand only count the strand that appears in the blast result(i.e., strands with collisions)
// 6. Each encoding scheme has a seperate folder

class TransformSelection{
public:
    TransformSelection(string path[4]); // path: the path to four encoding blast results folders
    ~TransformSelection() {}

    void Select(); // the main function of this class; select an encoding for each payload
    void PrintStatistics();
private:
    void ReadCollisions(string path[4], string strand_id_prefix); // read the blas collision result for four encodings
    vector<string> all_files[4];

    unordered_map<unsigned int, string> primer_id2name;
    unordered_map<string, unsigned int> primer_name2id;
    unordered_set<string> strands;

    set<string> all_primers; // all primer names(w/o collisions)

    unordered_map<unsigned int, vector<unsigned int> > collision_list[4]; // collision_list[i][j]: a set of primer ids of [encoding i, payload j]
    set<unsigned int> collided_primers; // collided_primers after selection
    // vector<bool> primer_is_collided;
    set<unsigned int> collided_primers_default; // collided_primers using the default encoding

    int n_primer = 0;
    int n_strand = 0;
    vector<int> selection_list; // selection_list[i]: the selection of encoding of the i_th payload
    // vector<string> strand_transformed_result;
};

#endif // __TRANSFORM_SELECTION__H__