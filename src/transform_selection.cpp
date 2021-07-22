#include "transform_selection.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <climits>
using namespace std;

TransformSelection::TransformSelection(string path[4]) {
    n_primer = 0;
    
    for (int i = 0; i < 4; i++) {
        string &cur_path = path[i];
        cout << "Current path:" << cur_path << endl << flush;
        ifstream myfile;
        myfile.open(cur_path);
        string line;

        // int times = 20;
        while(getline(myfile, line)) {
            istringstream iss(line);
            string current_field;
            // line example:
            // primer4	payload1176045	93.750	16	1	0	1	16	55	40	27	25.6
            if (line[0] == '#') {
                // cout << line << endl;
                iss >> current_field;
                iss >> current_field;
                if (i == 0 && current_field == "Query:") n_primer++;
                continue;
            }
            // cout << line << endl;

            string primer_name;
            string strand_name;
            unsigned int primer_id;
            unsigned int strand_id;
            iss >> primer_name >> strand_name;
            // cout << primer_name.substr(6) << ' ' << strand_name.substr(7) << endl; 
            primer_id = stoul(primer_name.substr(6));
            strand_id = stoul(strand_name.substr(7));

            /*
            for (int i = 1; i <= 6; i++) {
                iss >> current_field;
            }
            iss >> current_field;
            unsigned int strand_start = stoul(current_field);
            iss >> current_field;
            unsigned int strand_end = stoul(current_field);
            if (strand_start > strand_end) swap(strand_start, strand_end);
            */
            strand_id2name[strand_id] = strand_name;
            strand_name2id[strand_name] = strand_id;

            primer_id2name[primer_id] = primer_name;
            primer_name2id[primer_name] = primer_id;

            // cout << "i " << i << "strand" << strand_id << ' ' << "primer" << primer_id << endl;
            collision_list[i][strand_id].push_back(primer_id);
            // if(times > 0) {
            //     times--;
            //     cout << primer_id << ' ' << strand_id << endl;
            // }
        }
    }

    assert(strand_id2name.size() == strand_name2id.size());
    assert(primer_id2name.size() == primer_name2id.size());
    n_strand = strand_id2name.size();
    cout << "n_primer = " << n_primer << endl;
}

void TransformSelection::Select() {
    for (int i = 0; i < n_strand; i++) {
        vector<unsigned int> &current_collision_list = collision_list[0][i];
        for (auto it = current_collision_list.begin(); it != current_collision_list.end(); it++) {
            int primer = *it;
            collided_primers_default.insert(primer);
        }
    }

    for (int i = 0; i < n_strand; i++) {
        int n_collided_increased_min = 2000000000;
        int selection = -1;

        for (int j = 0; j < 4; j++) {
            vector<unsigned int> &current_collision_list = collision_list[j][i];
            int n_collided_increased = current_collision_list.size();
            for(auto it = current_collision_list.begin(); it != current_collision_list.end(); it++) {
                int primer = *it;
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

        vector<unsigned int> &current_collision_list = collision_list[selection][i];
        for(auto it = current_collision_list.begin(); it != current_collision_list.end(); it++) {
            int primer = *it;
            collided_primers.insert(primer);
        }    

        // selection_list.push_back(selection);
        // payload_transformed_result.push_back()
    }
}

void TransformSelection::PrintStatistics() {
    cout << "n_strand = " << n_strand << endl << flush;
    cout << "n_primer = " << n_primer << endl << flush;
    int n_available_primer = n_primer - collided_primers_default.size();
    cout << "Default encoding: n_available_primer =\t" << n_available_primer  << "\tratio = " << (double)n_available_primer/n_primer*100 << "%" << endl << flush;
    n_available_primer = n_primer - collided_primers.size();
    cout << "After selection: n_available_primer =\t" << n_available_primer << "\tratio = " << (double)n_available_primer/n_primer*100 << "%" << endl << flush;
}