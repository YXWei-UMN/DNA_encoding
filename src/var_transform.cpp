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
    
    const int STRAND_LENGTH = 200;
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
                if (strand_start / STRAND_LENGTH != strand_end / STRAND_LENGTH) continue;

                collision_list[i].push_back(make_tuple(strand_id, strand_start, strand_end, primer_id));
            }
        }

        assert(strand_id2name.size() == strand_name2id.size());
        assert(primer_id2name.size() == primer_name2id.size());

        // For each mapping, sort all of its collisions based on their position
        for (int i = 0; i < 4; i ++) {
            vector<Collision> &cur_list = collision_list[i];
            sort(cur_list.begin(), cur_list.end());
        }

        // Cut the strand accodring to:
        //  1. Recovered primer positions
        //  2. 200-length(not using for now)
        int pointer[4] = {0};
        vector<Collision> cur_strand_default_collsions;
        StrandID cur_strand_id = -1;
        // unsigned int start = UINT_MAX;
        // unsigned int end = 0;
        while(pointer[0] < collision_list[0].size() && 
            (pointer[1] < collision_list[1].size() || pointer[2] < collision_list[2].size() || pointer[3] < collision_list[3].size())
        ) {
            Collision &cur_collision = collision_list[0][pointer[0]];
            bool flag_new_strand = false;
            if (cur_strand_default_collsions.size() != 0) {
                int primer_id_cur = get<4>(cur_collision);
                if (discarded_primers.find(primer_id_cur) != discarded_primers.end()) {
                    pointer[0]++; // discarded_primers position is cut point -> skip this collision
                    flag_new_strand = true;
                } else {
                    unsigned int strand_id = get<1>(cur_collision) / STRAND_LENGTH;
                    if (strand_id != cur_strand_id) {
                        assert(strand_id > cur_strand_id);
                        flag_new_strand = true;
                    } 
                }
                
            }
            
            if (flag_new_strand) {
                assert(cur_strand_default_collsions.size() > 0);
                // assert(start < UINT_MAX);
                // assert(end > 0);

                unordered_set<PrimerID> increased_primer[4];
                for (auto it = cur_strand_default_collsions.begin(); it != cur_strand_default_collsions.end(); it++) {
                    PrimerID primer_id = get<3>(*it);
                    if (discarded_primers.find(primer_id) == discarded_primers.end())
                        increased_primer[0].insert(primer_id);
                }

                for (int i = 1; i < 4; i++) {
                    if (pointer[i] >= collision_list[i].size()) continue;
                    while(true) {
                        if (pointer[i] >= collision_list[i].size()) break;
                        Collision &c = collision_list[i][pointer[i]];
                        unsigned int cut = get_collisiton_cut_point(c);
                        StrandID strand_id = cut / STRAND_LENGTH;

                        if (strand_id < cur_strand_id) {
                            pointer[i]++;
                        } else if (strand_id == cur_strand_id) {
                            PrimerID primer_id = get<3>(c);
                            if (discarded_primers.find(primer_id) == discarded_primers.end())
                                increased_primer[i].insert(primer_id);
                        }
                    }
                }

                
            }
            else {
                if (cur_strand_default_collsions.size() == 0) {
                    cur_strand_id = get<3>(cur_collision);
                }
                cur_strand_default_collsions.push_back(cur_collision);
            }
        }

    }
}