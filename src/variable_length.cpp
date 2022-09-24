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

bool Homopolymers(string primer) {
    char last_one = primer[0];
    int max_length_of_homo=1;
    int cur_length_of_homo=1;

    for (int i = 1; i < 4; ++i) {
        if(primer[i] == last_one){
            cur_length_of_homo++;
            if(cur_length_of_homo>max_length_of_homo)
            {
                max_length_of_homo = cur_length_of_homo;
            }
        }
        else  cur_length_of_homo = 1;
        last_one = primer[i];
    }
    return max_length_of_homo;
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

    /*// read all primer library in
    ifstream primer_library;
    primer_library.open("/home/eason/CLionProjects/Primer-payload-collision/primer28K");
    unordered_map<PrimerID, string> primerlibrary;

    while (getline(primer_library, line)){
        PrimerID primer_id;
        primer_id = stoul(line.substr(7));

        getline(primer_library, line);
        primerlibrary.emplace(primer_id,line);
    }*/

    bool if_checked = false;
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
                if_checked= true;
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
    int total_collided_primer = primer_collision_num_.size();
    

    long long total_collision_num=0;

    ofstream primer_collision_num;
    primer_collision_num.open ("primer_collision_num.csv",ios::out | ios::trunc);
    vector<int> primer_distribution(200010,0);
    for(auto n:primer_collision_num_){
        int collision_num = n.second;
        total_collision_num+=collision_num;
        primer_distribution[collision_num]++;
    }

    for(int i=0; i < primer_distribution.size(); i++){
        primer_collision_num<<i<<","<<primer_distribution[i]<<endl;
    }
    primer_collision_num.close();
    /*cout<<"total collision: "<<total_collision_num<<"  collided primer: "<<total_collided_primer<<endl;
    cout<<"avg collision per primer: "<<total_collision_num/(1.0*total_collided_primer)<<endl;*/


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

    for (int i = 0; i < primer_process_order.size(); i++) {
        PrimerID current_primer_id = primer_process_order[i].second;

        if(discarded_primers.find(current_primer_id) != discarded_primers.end()) {
            continue;
        }

        // recover current primer = push its conflict to discarded
        unordered_set<PrimerID> &conflicts = primer_confilct_list[current_primer_id];

        for (auto it = conflicts.begin(); it != conflicts.end(); it++){
            discarded_primers.insert(*it);
        }
    }

    PrintStatistics(total_collided_primer);
}

void VariableLength::PrintStatistics(int total_collided_primer) {
    cout << "n_strand: " << n_strand << endl << flush;
    cout << "n_primer: " << n_primer << endl << flush;
    int n_discarded = discarded_primers.size();
    int n_recovered = total_collided_primer - n_discarded;
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
