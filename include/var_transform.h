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
    ~VarTransform() { delete variable_length; }

    void Run();
    void PrintStatistics();

private:

    VariableLength *variable_length;

    void ReadCollisions(string path); // read the blast collision result of a single file

    vector<string> all_files[4];
    unordered_set<string> all_primers; // all primer names(w/o collisions)
    unordered_set<PrimerID> discarded_primers;
    unordered_set<PrimerID> recovered_primers;
    unordered_set<PrimerID> collided_primers_default;

    vector<Collision> collision_list[4]; // collision_list[i][j]: a set of primer ids of [encoding i, payload j]

    int n_primer = 0;
    int n_strand = 0;

    vector<int> selection_list; 

};
#endif