#include<fstream>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <climits>
#include <sstream>
#include <algorithm>
#include "global.h"
#include "variable_length.h"
#include "tool.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>
using namespace std;


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
    cout << "VarLen: Processing " << path << endl;
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

        if (primer_collision_num_.find(primer_id)==primer_collision_num_.end()){
            primer_collision_num_.emplace(primer_id,1);
        } else{
            primer_collision_num_.find(primer_id)->second++;
        }

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
    }

    sort(collision_linear_order.begin(), collision_linear_order.end(), CollisionPositionCMP);

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

            assert(cut_last <= cut);
            if (IsBlindSpot(cut-cut_last)) {
                if (primer_id_last == primer_id) {
                    discarded_primers.insert(primer_id);
                } else {
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

bool sortbyfirst_asending(const pair<int,PrimerID> &a, const pair<int,PrimerID> &b){
    return a.first<b.first;
}

bool sortbyfirst_descending(const pair<int,PrimerID> &a, const pair<int,PrimerID> &b){
    return a.first>b.first;
}

void VariableLength::Cut() {
    for (int i = 0; i < all_files.size(); i++) {
        ReadCollisions(all_files[i]);
    }

    int ideal_capacity = 1.55*1000000*200/2; //devide by 2 since it's a primer not a primer pair
    for (auto it : primer_collision_num_){
        double capacity = ideal_capacity - it.second*4*100; // we are using 64GB, collision num * 4 to scale to 200+GB
        capacity/=1024; // in case overflow
        if (capacity<0){
            discarded_primers.insert(it.first);
            continue;
        }
        primer_capacity_.emplace(it.first,capacity);
    }

    // prescreen self-confict primers / capacity < 0
    for (auto it = discarded_primers.begin(); it != discarded_primers.end(); it++) {
        primer_confilct_list.erase(*it);
        primer_collision_num_.erase(*it);
        primer_capacity_.erase(*it);
    }

    if (g_var_len_algorithm==1){
        for (auto it : primer_collision_num_){
            primer_process_order.push_back(make_pair(it.second, it.first));
        }
        sort(primer_process_order.begin(), primer_process_order.end(),sortbyfirst_asending);
    } else if (g_var_len_algorithm==2){
        for (auto it = primer_confilct_list.begin(); it != primer_confilct_list.end(); it++) {
            primer_process_order.push_back(make_pair(it->second.size(), it->first));
        }
        sort(primer_process_order.begin(), primer_process_order.end(),sortbyfirst_asending);
    } else if (g_var_len_algorithm==3){
        // calculate capacity diff, push to primer_process_order
        for (auto it : primer_capacity_){
            double capacity_diff=it.second;
            double tmp = capacity_diff;
            for (auto conflict_primer : primer_confilct_list[it.first]){
                capacity_diff -= primer_capacity_[conflict_primer];
            }
            if (capacity_diff>tmp) cout<<"large capacity_diff?? overflow"<<endl;
            primer_process_order.push_back(make_pair(capacity_diff, it.first));
        }
        sort(primer_process_order.begin(), primer_process_order.end(),sortbyfirst_descending);
    }

    //initialize look ahead window
    unordered_set<PrimerID> look_ahead_window;

    // used to mark the primer that is swapped with its conflict, those primer has been processed already but still stay in primer_process_order
    unordered_set<PrimerID> swapped;

    for (int i = 0; i < primer_process_order.size(); i++) {
        PrimerID current_primer_id = primer_process_order[i].second;

        if(discarded_primers.find(current_primer_id) != discarded_primers.end()) {
            continue;
        }

        if(swapped.find(current_primer_id) != swapped.end()) {
            continue;
        }

        look_ahead_window.clear();

        if (g_var_len_algorithm==3){

            int j=i+1;
            int count=0;
            while (count<=200) {
                if (j>=primer_process_order.size()) break;

                if(discarded_primers.find(primer_process_order[j].second) != discarded_primers.end()) {
                    j++;
                    continue;
                }

                if(swapped.find(primer_process_order[j].second) != swapped.end()) {
                    j++;
                    continue;
                }

                look_ahead_window.insert(primer_process_order[j].second);
                j++;
                count++;
            }
            //get current primers' conflicts in the window
            auto conflicts = primer_confilct_list.find(current_primer_id)->second;
            unordered_set<PrimerID> window_conflicts;
            //cout<<"before "<<conflicts.size()<<endl;
            for (auto it:look_ahead_window){
                if (conflicts.find(it)!=conflicts.end()){
                    window_conflicts.insert(it);
                }
            }
            //cout<<"after "<<window_conflicts.size()<<endl;

            // calculate window capacity diff for current primer
            int window_capacity_diff_of_current_primer = primer_capacity_.find(current_primer_id)->second;
            //cout<<"primer "<<current_primer_id<<" window capacity diff(1):"<<window_capacity_diff_of_current_primer<<endl;
            for (auto it:window_conflicts){
                window_capacity_diff_of_current_primer -= primer_capacity_.find(it)->second;
            }
            //cout<<"primer "<<current_primer_id<<" window capacity diff(2):"<<window_capacity_diff_of_current_primer<<endl;

            //only if capacity diff < 0 then swap
            if (window_capacity_diff_of_current_primer<0){
                vector<pair<int,PrimerID>> window_capacity_diff_of_window_conflicts;

                for (auto it:window_conflicts){
                    //get conflict primers' conflicts in the window
                    auto conflicts_of_window_conflicts = primer_confilct_list.find(it)->second;
                    unordered_set<PrimerID> window_conflicts_of_window_conflicts;
                    for (auto n:look_ahead_window){
                        if (conflicts_of_window_conflicts.find(n)!=conflicts_of_window_conflicts.end()){
                            window_conflicts_of_window_conflicts.insert(n);
                        }
                    }
                    //calculate window capacity diff for conflict primer
                    int window_capacity_diff_of_conflit_primer = primer_capacity_.find(it)->second;
                    for (auto n:window_conflicts_of_window_conflicts){
                        window_capacity_diff_of_conflit_primer -= primer_capacity_.find(n)->second;
                    }
                    window_capacity_diff_of_window_conflicts.push_back(make_pair(window_capacity_diff_of_conflit_primer,it));
                }

                sort(window_capacity_diff_of_window_conflicts.begin(), window_capacity_diff_of_window_conflicts.end(),sortbyfirst_descending);

                if(window_capacity_diff_of_window_conflicts.begin()->first > window_capacity_diff_of_current_primer){
                    //cout<<"and swap with primer "<<window_capacity_diff_of_window_conflicts.begin()->second<<"  capacity_diff:"<<window_capacity_diff_of_window_conflicts.begin()->first<<endl;
                    //swap currrent primer and the max conflict primer
                    current_primer_id=window_capacity_diff_of_window_conflicts.begin()->second;
                    swapped.insert(window_capacity_diff_of_window_conflicts.begin()->second);
                }
            }
        }

        // recover current primer = push its conflict to discarded
        unordered_set<PrimerID> &conflicts = primer_confilct_list[current_primer_id];

        for (auto it = conflicts.begin(); it != conflicts.end(); it++){
            discarded_primers.insert(*it);
        }
    }

    PrintStatistics();
}

void VariableLength::PrintStatistics() {
    cout << "n_strand: " << n_strand << endl << flush;
    cout << "n_primer: " << n_primer << endl << flush;
    int n_discarded = discarded_primers.size();
    int n_recovered = primer_collision_num_.size() - n_discarded;
    int n_free = n_primer - n_discarded - n_recovered;
    cout << "n[free, recovered, discarded] = " << n_free << ' ' << n_recovered << ' ' << n_discarded << endl << flush;
    cout << "Available primer ratio before VarLen = " << (double)n_free/n_primer << endl << flush;
    cout << "Available primer ratio after VarLen = " << (double)(n_free+n_recovered)/n_primer << endl << flush;

    double recovered_capacity=0;
    for (auto i:primer_capacity_) {
        if (discarded_primers.find(i.first)==discarded_primers.end()){
            recovered_capacity+=i.second;
        }
    }
    recovered_capacity = recovered_capacity*1.57/8/1024/1024; // has devided by 1024*1024 already in case of over flow
    cout<<"recovered_capacity = " << recovered_capacity<<endl;


}
