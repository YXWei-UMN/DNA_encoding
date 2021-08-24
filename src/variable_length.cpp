#include<fstream>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <climits>
#include <sstream>
#include <algorithm>

#include "variable_length.h"
#include "tool.h"

using namespace std;

bool CollisionPositionCMP(const Collision &c1,const Collision &c2)
{
    StrandID strand_id1 = get<0>(c1);
    unsigned int start1 = get<1>(c1);
    unsigned int end1 = get<2>(c1);
    PrimerID primer_id1 = get<3>(c1);
	unsigned int cut1 = start1 + (end1 - start1) / 2;
    cut1 = (int)(cut1 * 0.1 + 0.5);

    StrandID strand_id2 = get<0>(c2);
    unsigned int start2 = get<1>(c2);
    unsigned int end2 = get<2>(c2);
    PrimerID primer_id2 = get<3>(c2);
	unsigned int cut2 = start2 + (end2 - start2) / 2;
    cut2 = (int)(cut2 * 0.1 + 0.5);
    
    if (strand_id1 != strand_id2) return strand_id1 < strand_id2;
    if (cut1 != cut2) return cut1 < cut2;
    if (start1 != start2) return start1 < start2;
    if (end1 != end2) return end1 < end2;
    return primer_id1 < primer_id2;
}

VariableLength::VariableLength(string path) {
    n_primer = 0;
    all_files = listFiles(path, true);
    cout << "The number of files: " << all_files.size() << endl;

    //Initialize blind spot
    for (int i = 0;i<=44;i++){
        blind_spot.push_back(false);
    }
    for (int i=1;i<=14;i++){
        blind_spot[i]=true;
    }
    blind_spot[17]=true;
    blind_spot[18]=true;
    for (int i=21;i<=29;i++){
        blind_spot[i]=true;
    }
    blind_spot[33]=true;
    blind_spot[37]=true;
    for (int i=41;i<=44;i++){
        blind_spot[i]=true;
    }

}

bool VariableLength::IsBlindSpot(unsigned int distance) {
    if (distance >= 45) return false;
    return blind_spot[distance];
}

void VariableLength::ReadCollisions(string path) {
    cout << "Processing " << path << endl;
    strand_id2name.clear();
    strand_name2id.clear();
    primer_id2name.clear();
    primer_name2id.clear();
    collision_linear_order.clear();

    ifstream myfile;
    myfile.open(path);
    string line;

    while(getline(myfile, line)) {
        istringstream iss(line);
        string current_field;
        // line example:
        // primer4	payload1176045	93.750	16	1	0	1	16	55	40	27	25.6
        if (line[0] == '#') {
            iss >> current_field;
            iss >> current_field;
            if (current_field == "Query:") {
                iss >> current_field;
                all_primers.insert(current_field);
            }
            continue;
        }

        string primer_name;
        string strand_name;
        PrimerID primer_id;
        StrandID strand_id;
        iss >> primer_name >> strand_name;
        primer_id = stoul(primer_name.substr(6));
        strand_id = stoul(strand_name.substr(7));
        primer_id2name[primer_id] = primer_name;
        primer_name2id[primer_name] = primer_id;
        strand_id2name[strand_id] = strand_name;
        strand_name2id[strand_name] = strand_id;

        // cout << "PrimerID: " << primer_id << endl << flush;
        // cout << "StrandID: " << strand_id << endl << flush;

        // read collsion position
        for (int i = 0; i < 6; i++) {
            iss >> current_field;
        }
        iss >> current_field;
        unsigned int strand_start = stoul(current_field);
        iss >> current_field;
        unsigned int strand_end = stoul(current_field);
        if (strand_start > strand_end) swap(strand_start, strand_end);
        collision_linear_order.push_back(make_tuple(strand_id, strand_start, strand_end, primer_id));

        // cout << "strand_start: " << strand_start << endl << flush;
        // cout << "strand_end: " << strand_end << endl << flush;
    }

    sort(collision_linear_order.begin(), collision_linear_order.end(), CollisionPositionCMP);
    // cout << "Collision linear order: " << endl;
    // for (auto it = collision_linear_order.begin(); it != collision_linear_order.end(); ++it) {
    //     cout << get<0>(*it) << ' ' << get<1>(*it) << ' ' << get<2>(*it) << ' ' << get<3>(*it) << endl;
    // }
    for (int i = 1; i < collision_linear_order.size(); i++) {
        auto &cur_collision = collision_linear_order[i];
        StrandID strand_id = get<0>(cur_collision);
        unsigned int start = get<1>(cur_collision);
        unsigned int end = get<2>(cur_collision);
        PrimerID primer_id = get<3>(cur_collision);

        unsigned int cut = start + (end - start) / 2;
        cut = (int)(cut * 0.1 + 0.5);
        
        auto &last_collision = collision_linear_order[i-1];
        StrandID strand_id_last = get<0>(last_collision);
        if (strand_id_last == strand_id) {
            unsigned int start_last = get<1>(last_collision);
            unsigned int end_last = get<2>(last_collision);
            PrimerID primer_id_last = get<3>(last_collision);

            unsigned int cut_last = start_last + (end_last - start_last)/2;
            cut_last = (int)(cut_last * 0.1 + 0.5);

            // cout << "Primer[" << primer_id_last << "," << primer_id << "]: cut[" << cut_last << "," << cut << "]" << endl << flush;

            assert(cut_last <= cut);
            if (IsBlindSpot(cut-cut_last)) {
                if (primer_id_last == primer_id) {
                    discarded_primers.insert(primer_id);
                } else {
                    // cout << "Conflict pair: primer" << primer_id_last << " primer" << primer_id << endl << flush;
                    primer_confilct_list[primer_id].insert(primer_id_last);
                    primer_confilct_list[primer_id_last].insert(primer_id);
                }
            }

        }
    }
    
    assert(strand_id2name.size() == strand_name2id.size());
    assert(primer_id2name.size() == primer_name2id.size());
    n_strand += strand_id2name.size();
    n_primer = all_primers.size();
}


void VariableLength::Cut() {
    for (int i = 0; i < all_files.size(); i++) {
        ReadCollisions(all_files[i]);
    }

    for (auto it = primer_confilct_list.begin(); it != primer_confilct_list.end(); it++) {
        primer_process_order.push_back(make_pair(it->second.size(), it->first));
    }
    sort(primer_process_order.begin(), primer_process_order.end());

    // for (auto it = primer_process_order.begin(); it != primer_process_order.end(); it++) {
    //     cout << "Primer" << it->second << " has " << it->first <<  " dense collisions" << endl;
    // }

    for (int i = 0; i < primer_process_order.size(); i++) {
        PrimerID primer_id = primer_process_order[i].second;
        if(discarded_primers.find(primer_id) != discarded_primers.end()) {
            continue;
        }
        unordered_set<PrimerID> &conflicts = primer_confilct_list[primer_id];


        for (auto it = conflicts.begin(); it != conflicts.end(); it++){
            discarded_primers.insert(*it);
        }
        
        recovered_primers.insert(primer_id);
    }

    // cout << "Discarded primers: ";
    // for (auto it = discarded_primers.begin(); it != discarded_primers.end(); it++) {
    //     cout << (*it) << " ";
    // }
    // cout << endl;

    // cout << "Recovered primers: ";
    // for (auto it = recovered_primers.begin(); it != recovered_primers.end(); it++) {
    //     cout << (*it) << " ";
    // }
    // cout << endl;

    PrintStatistics();
}

void VariableLength::PrintStatistics() {
    cout << "n_strand: " << n_strand << endl << flush;
    cout << "n_primer: " << n_primer << endl << flush;
    int n_recovered = recovered_primers.size();
    int n_discarded = discarded_primers.size();
    int n_free = n_primer - n_discarded - n_recovered;
    cout << "n[free, recovered, discarded] = " << n_free << ' ' << n_recovered << ' ' << n_discarded << endl << flush;
    cout << "Available primer ratio before VarLen = " << (double)n_free/n_primer << endl << flush;
    cout << "Available primer ratio after VarLen = " << (double)(n_free+n_recovered)/n_primer << endl << flush;
}