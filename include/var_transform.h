#ifndef DNA_ENCODING_TRANSFORMATION_H
#define DNA_ENCODING_TRANSFORMATION_H

#include<unordered_map>
#include<vector>
#include<set>
#include<string>
#include<unordered_map>
#include<unordered_set>
#include<tuple>

#include "variable_length.h"

using namespace std;

typedef unsigned int PrimerID;
typedef unsigned int StrandID;
typedef tuple<StrandID, unsigned int, unsigned int, PrimerID> Collision;
typedef int CollisionPositionIndex;

class VarTransform {
public:
    VarTransform(string path[4]);
    ~VarTransform() {}

    void run();
    void PrintStatistics();

private:

    VariableLength *variable_length;

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

    set<string> all_primers; // all primer names(w/o collisions)
    unordered_map<PrimerID, unordered_set<PrimerID> > primer_confilct_list; // conflict between primers (primers that cannot exist together)
    vector<pair<int, PrimerID>> primer_process_order;
    unordered_set<PrimerID> discarded_primers;
    unordered_set<PrimerID> recovered_primers;

    int n_primer = 0;
    int n_strand = 0;

    vector<string> all_files[4];
    set<string> all_primers; // all primer names(w/o collisions)

    vector<Collision> collision_list[4]; // collision_list[i][j]: a set of primer ids of [encoding i, payload j]

    int n_primer = 0;
    int n_strand = 0;
    vector<int> selection_list; 

};
#endif