#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include <utility>
#include <climits>

#include "var_transform.h"
#include "tool.h"

using namespace std;

VarTransform::VarTransform(string path[4]) {
    n_primer = 0;

    for (int i = 0; i < 4; i++) {
        all_files[i] = listFiles(path[i], true);
        sort(all_files[i].begin(), all_files[i].end());
    }

    assert(all_files[0].size() == all_files[1].size());
    assert(all_files[0].size() == all_files[2].size());
    assert(all_files[0].size() == all_files[3].size());

    variable_length = new VariableLength(path[0]);

}

void VarTransform::run() {
    variable_length->Cut();
    unordered_set<PrimerID> recovered_primers = variable_length->get_recovered_primers();
    unordered_set<PrimerID> discarded_primers = variable_length->get_discarded_primers();

    for (int file_id = 0; file_id < all_files[0].size(); file_id++){
        cout << "Selecting tranformation in file " << file_id << endl;
        string path[4] = {all_files[0][file_id], all_files[1][file_id], all_files[2][file_id], all_files[3][file_id]};

        for (int i = 0; i < 4; i++) {
            string &cur_path = path[i];

            ifstream myfile;
            myfile.open(cur_path);
            string line;

            while(getline(myfile, line)) {
                istringstream iss(line);
                string current_field;

                // line example:
                // primer4	payload1176045	93.750	16	1	0	1	16	55	40	27	25.6
                if (line[0] == '#') {
                    iss >> current_field;
                    iss >> current_field;
                    if (i == 0 && current_field == "Query:") {
                        iss >> current_field;
                        all_primers.insert(current_field);
                    }
                    continue;
                }

                string primer_name;
                string strand_name;
                unsigned int primer_id;
                unsigned int strand_id;
                iss >> primer_name >> strand_name;
                primer_id = stoul(primer_name.substr(6));
                strand_id = stoul(strand_name.substr(7));

                strand_id2name[strand_id] = strand_name;
                strand_name2id[strand_name] = strand_id;

                primer_id2name[primer_id] = primer_name;
                primer_name2id[primer_name] = primer_id;

                // read collsion position
                for (int i = 0; i < 6; i++) {
                    iss >> current_field;
                }
                iss >> current_field;
                unsigned int strand_start = stoul(current_field);
                iss >> current_field;
                unsigned int strand_end = stoul(current_field);
                if (strand_start > strand_end) swap(strand_start, strand_end);

                collision_list[i][strand_id].push_back(make_tuple(strand_id, strand_start, strand_end, primer_id));
            }
        }

        assert(strand_id2name.size() == strand_name2id.size());
        assert(primer_id2name.size() == primer_name2id.size());

        // For each strand, sort all of its collisions based on their position
        for (int i = 0; i < 4; i ++) {
            unordered_map<unsigned int, vector<Collision> > &cur_collision_list = collision_list[i];
            for (auto it = cur_collision_list.begin(); it != cur_collision_list.end(); ++it) {
                StrandID strand_id = it->first;
                vector<Collision> &strand_collision = it->second;
                sort(strand_collision.begin(), strand_collision.end());
            }
        }

        // Cut the strand accodring to:
        //  1. Recovered primer positions
        //  2. 200-length(not using for now)
        unordered_map<unsigned int, vector<Collision> > &default_collision_list = collision_list[0];
        for (auto it = default_collision_list.begin(); it != default_collision_list.end(); ++it) {
            StrandID strand_id = it->first;
        }
    }
}