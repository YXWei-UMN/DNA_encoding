#include "transform_selection.h"
#include "tool.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <climits>
#include <string>
using namespace std;

TransformSelection::TransformSelection(string path[4]) {
    n_primer = 0;

    for (int i = 0; i < 4; i++) {
        all_files[i] = listFiles(path[i], true);
        sort(all_files[i].begin(), all_files[i].end());
    }

    assert(all_files[0].size() == all_files[1].size());
    assert(all_files[0].size() == all_files[2].size());
    assert(all_files[0].size() == all_files[3].size());

    cout << "The number of seperated files: " << all_files[0].size() << endl;
}

void TransformSelection::ReadCollisions(string path[4], string strand_id_prefix) {
    primer_id2name.clear();
    primer_name2id.clear();
    for (int i = 0; i < 4; i++) {
        collision_list[i].clear();
    }

    const int STRAND_LENGTH = 200;
    
    for (int i = 0; i < 4; i++) {
        string &cur_path = path[i];
        ifstream myfile;
        myfile.open(cur_path);
        string line;

        while(getline(myfile, line)) {
            istringstream iss(line);
            string current_field;
            // line example:
            //      primer4	payload1176045	93.750	16	1	0	1	16	55	40	27	25.6
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

            // read collsion position
            for (int i = 0; i < 6; i++) {
                iss >> current_field;
            }
            iss >> current_field;
            unsigned int strand_start = stoul(current_field);
            iss >> current_field;
            unsigned int strand_end = stoul(current_field);
            if (strand_start > strand_end) swap(strand_start, strand_end);

            strand_id = strand_start / STRAND_LENGTH;
            if (strand_id != strand_end / STRAND_LENGTH) continue;

            // string str_id =  strand_id_prefix + to_string(strand_id);
            // strands.insert(str_id);
            primer_id2name[primer_id] = primer_name;
            primer_name2id[primer_name] = primer_id;

            // cout << "i " << i << "strand" << strand_id << ' ' << "primer" << primer_id << endl;
            collision_list[i][strand_id].push_back(primer_id);
        }
    }

    assert(primer_id2name.size() == primer_name2id.size());
    // n_strand = strands.size();
    n_primer = all_primers.size();
}

void TransformSelection::Select() {
    for (int file_id = 0; file_id < all_files[0].size(); file_id++){
        cout << "Transform selection: processing file #" << file_id << endl;
        string path[4] = {all_files[0][file_id], all_files[1][file_id], all_files[2][file_id], all_files[3][file_id]};
        ReadCollisions(path, to_string(file_id)); 

        unordered_map<unsigned int, vector<unsigned int> > &collision_list_default = collision_list[0];
        for(auto it = collision_list_default.begin(); it != collision_list_default.end(); it++) {
            vector<unsigned int> &current_collision_list = it->second;
            for (auto it2 = current_collision_list.begin(); it2 != current_collision_list.end(); it2++) {
                int primer = *it2;
                collided_primers_default.insert(primer);
            }
        }

        /*
        set< pair<int, unsigned int> > strand_order; // pair: collision_list.size(), strand_id
                                                    // we sort strand ascendingly based on the collision_list size
        for (int i = 0; i < n_strand; i++) {
            int max_size = 0;
            for (int j = 0; j < 4; j++) {
                int cur_size = collision_list[j][i].size();
                if (cur_size > max_size) {
                    max_size = cur_size;
                } 
            }
            strand_order.insert(make_pair(max_size, i));
        }
        

        for (auto it = strand_order.begin(); it != strand_order.end(); it++) {
        */

        for(auto it = collision_list[0].begin(); it != collision_list[0].end(); it++){
            int n_collided_increased_min = 2000000000;
            int selection = -1;

            unsigned int strand_id = it->first;

            for (int j = 0; j < 4; j++) {
                vector<unsigned int> &current_collision_list = collision_list[j][strand_id];
                int n_collided_increased = current_collision_list.size();
                for(auto it2 = current_collision_list.begin(); it2 != current_collision_list.end(); it2++) {
                    int primer = *it2;
                    if (collided_primers.find(primer) != collided_primers.end()) {
                        n_collided_increased--;
                    }
                }

                assert(n_collided_increased >= 0);

                if (n_collided_increased < n_collided_increased_min) {
                    n_collided_increased_min = n_collided_increased;
                    selection = j;
                }
            }

            vector<unsigned int> &current_collision_list = collision_list[selection][strand_id];
            for(auto it = current_collision_list.begin(); it != current_collision_list.end(); it++) {
                int primer = *it;
                collided_primers.insert(primer);
            }
        }
    }


    
}

void TransformSelection::PrintStatistics() {
    cout << "n_primer = " << n_primer << endl << flush;
    int n_available_primer = n_primer - collided_primers_default.size();
    cout << "Default encoding: n_available_primer =\t" << n_available_primer  << "\tratio = " << (double)n_available_primer/n_primer*100 << "%" << endl << flush;
    n_available_primer = n_primer - collided_primers.size();
    cout << "After selection: n_available_primer =\t" << n_available_primer << "\tratio = " << (double)n_available_primer/n_primer*100 << "%" << endl << flush;
}