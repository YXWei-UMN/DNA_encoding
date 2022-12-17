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

    vector<int> combinations_strand200_4len = {150,160,190,200,300,310,320,340,350,360,380,390,400,450};
    vector<int> combinations_strand200_4len_150_160_180_200 = {150,160,180,200,300,310,320,330,340,350,360,380,400,450,460,470,480,490,500,510,520,530,540,550,560,580,600};
    vector<int> combinations_strand200_4len_150_170_180_200 = {150,170,180,200,300,320,330,340,350,360,370,380,400,450,470,480,490,500,510,520,530,540,550,560,570,580,600,620};
    vector<int> combinations_strand200_4len_150_170_190_200 = {150,170,190,200,300,320,340,350,360,370,380,390,400,450,470,490,500,510,520,530,540,550,560,570,580,590,600,620,640};
    vector<int> combinations_strand200_4len_150_180_190_200 = {150,180,190,200,300,330,340,350,360,370,380,390,400,450,480,490,500,510,520,530,540,550,560,570,580,590,600,630};

    vector<int> combinations_strand200_4len_160_170_180_200 = {160,170,180,200,320,330,340,350,360,370,380,400,480,490,500,510,520,530,540,550,560,570,580,600,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,800};
    vector<int> combinations_strand200_4len_160_170_190_200 = {160,170,190,200,320,330,340,350,360,370,380,390,400,480,490,500,510,520,530,540,550,560,570,580,590,600,640};
    vector<int> combinations_strand200_4len_160_180_190_200 = {160,180,190,200,320,340,350,360,370,380,390,400,480,500,510,520,530,540,550,560,570,580,590,600,640,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,820};
    vector<int> combinations_strand200_4len_170_180_190_200 = {170,180,190,200,340,350,360,370,380,390,400,510,520,530,540,550,560,570,580,590,600,680};

    vector<int> combinations_strand200_4len_160_170_180_190_200 = {160,170,180,190,200,320,330,340,350,360,370,380,390,400,480,490,500,510,520,530,540,550,560,570,580,590,600,640};
    vector<int> combinations_strand200_4len_150_160_170_180_190_200 = {150,160,170,180,190,200,300,310,320,330,340,350,360,370,380,390,400,450};



    vector<int> combinations_strand400_4len = {300,310,360,400,600,610,620,660,670,700,710,720,760,800,900,910,920,930,960,970,980,1000,1010,1020,1030,1060,1070,1080,1100,1110,1120,1160,1200,1210,1220,1230,1240,1260,1270,1280,1290,1300,1310,1320,1330,1340,1360,1370,1380,1390,1400,1410,1420,1430,1440,1460,1470,1480,1500};
    vector<int> combinations_strand400_8len = {300,310,340,360,370,380,390,400,600,610,620,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,900};

    vector<int> combinations_strand800_4len = {610,620,690,800,1220,1230,1240,1300,1310,1380,1410,1420,1490,1600,1830,1840,1850,1860,1910,1920,1930,1990,2000,2020,2030,2040,2070,2100,2110,2180,2210,2220,2290,2400,2440,2450,2460,2470,2480,2520,2530,2540,2550,2600,2610,2620,2630,2640,2650,2660,2680,2690,2710,2720,2730,2760,2790,2800,2820,2830,2840,2870,2900,2910,2980,3010,3020,3050,3060,3070,3080,3090,3100,3130,3140,3150,3160,3170,3200,3210,3220,3230,3240,3250,3260,3270,3280,3290,3300,3310,3320,3330,3340,3350,3370,3380,3400,3410,3420,3430,3440,3450,3460,3480,3490,3510,3520,3530,3560,3590,3600,3620,3630,3640,3660,3670,3680,3690,3700,3710,3720,3740,3750,3760,3770,3780,3790,3810,3820,3830,3840,3850,3860,3870,3880,3890,3900,3910,3920,3930,3940,3950,3960,3970,3980,3990,4000,4010,4020,4030,4040,4050,4060,4070,4080,4090,4100,4110,4120,4130,4140,4150,4170,4180,4200};
    vector<int> combinations_strand800_8len = {600,610,640,660,750,760,790,800,1200,1210,1220,1240,1250,1260,1270,1280,1300,1320,1350,1360,1370,1390,1400,1410,1420,1430,1440,1450,1460,1500,1510,1520,1540,1550,1560,1580,1590,1600,1800};

    vector<int> *coms = &combinations_strand400_4len;;
    if (g_len_group=="400_4") coms = &combinations_strand400_4len;
    else if (g_len_group=="400_8") coms = &combinations_strand400_8len;
    else if (g_len_group=="800_4") coms = &combinations_strand800_4len;
    else if (g_len_group=="800_8") coms = & combinations_strand800_8len;
    else if (g_len_group=="200_4") coms = & combinations_strand200_4len;
    else if (g_len_group=="150_160_180_200") coms = & combinations_strand200_4len_150_160_180_200;
    else if (g_len_group=="150_170_180_200") coms = & combinations_strand200_4len_150_170_180_200;
    else if (g_len_group=="150_170_190_200") coms = & combinations_strand200_4len_150_170_190_200;
    else if (g_len_group=="150_180_190_200") coms = & combinations_strand200_4len_150_180_190_200;

    else if (g_len_group=="160_170_180_200") coms = & combinations_strand200_4len_160_170_180_200;
    else if (g_len_group=="160_170_190_200") coms = & combinations_strand200_4len_160_170_190_200;
    else if (g_len_group=="160_180_190_200") coms = & combinations_strand200_4len_160_180_190_200;
    else if (g_len_group=="170_180_190_200") coms = & combinations_strand200_4len_170_180_190_200;
    else if (g_len_group=="160_170_180_190_200") coms = & combinations_strand200_4len_160_170_180_190_200;
    else if (g_len_group=="150_160_170_180_190_200") coms = & combinations_strand200_4len_150_160_170_180_190_200;
    else cerr << "what len group?" << endl << flush;

    blind_spot = vector<bool>(coms->back()/10, true);
    for (int n: (*coms)) {
        blind_spot[n/10]=false;
    }


    // if (g_len_group=="400_4"){
    //     // blind_spot = vector<bool>(combinations_strand400_4len.back()/10, true);
    //     for (int i = 0;i<=150;i++){
    //         blind_spot.push_back(true);
    //     }
    //     for (auto n:combinations_strand400_4len){
    //         blind_spot[n/10]=false;
    //     }
    // }else if (g_len_group=="400_8"){
    //     for (int i = 0;i<=90;i++){
    //         blind_spot.push_back(true);
    //     }
    //     for (auto n:combinations_strand400_8len){
    //         blind_spot[n/10]=false;
    //     }
    // } else if (g_len_group=="800_4"){
    //     for (int i = 0;i<=420;i++){
    //         blind_spot.push_back(true);
    //     }
    //     for (auto n:combinations_strand800_4len){
    //         blind_spot[n/10]=false;
    //     }
    // }else if (g_len_group=="800_8"){
    //     for (int i = 0;i<=180;i++){
    //         blind_spot.push_back(true);
    //     }
    //     for (auto n:combinations_strand800_8len){
    //         blind_spot[n/10]=false;
    //     }
    // } else cout<<"what len group?"<<endl;

    //Initialize blind spot


    /*blind_spot[30]=false;
    blind_spot[31]=false;
    blind_spot[36]=false;
    blind_spot[40]=false;
    blind_spot[60]=false;
    blind_spot[61]=false;
    blind_spot[62]=false;
    blind_spot[66]=false;
    blind_spot[67]=false;
    blind_spot[70]=false;
    blind_spot[71]=false;
    blind_spot[72]=false;
    blind_spot[76]=false;
    blind_spot[80]=false;
    blind_spot[90]=false;
    blind_spot[91]=false;

    blind_spot[92]=false;
    blind_spot[93]=false;
    blind_spot[96]=false;
    blind_spot[97]=false;
    blind_spot[98]=false;
    blind_spot[100]=false;
    blind_spot[101]=false;
    blind_spot[102]=false;
    blind_spot[103]=false;
    blind_spot[106]=false;
    blind_spot[107]=false;
    blind_spot[108]=false;
    blind_spot[110]=false;
    blind_spot[111]=false;
    blind_spot[112]=false;
    blind_spot[116]=false;

    blind_spot[120]=false;
    blind_spot[121]=false;
    blind_spot[122]=false;
    blind_spot[123]=false;
    blind_spot[124]=false;
    blind_spot[126]=false;
    blind_spot[127]=false;
    blind_spot[128]=false;
    blind_spot[129]=false;
    blind_spot[130]=false;
    blind_spot[131]=false;
    blind_spot[132]=false;
    blind_spot[133]=false;
    blind_spot[134]=false;
    blind_spot[136]=false;
    blind_spot[137]=false;

    blind_spot[138]=false;
    blind_spot[139]=false;
    blind_spot[140]=false;
    blind_spot[141]=false;
    blind_spot[142]=false;
    blind_spot[143]=false;
    blind_spot[144]=false;
    blind_spot[146]=false;
    blind_spot[147]=false;
    blind_spot[148]=false;
    blind_spot[150]=false;*/

}

bool VariableLength::IsBlindSpot(unsigned long long distance) {
    // blind spot start from 0, the size should minus 1
    //cout<<distance<<" ";
    if (distance >= blind_spot.size()-1) return false;
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
        total_collision_num++;
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
       /* for (int i = 0; i < 6; i++) {
            iss >> current_field;
        }
        iss >> current_field;
        unsigned int strand_start = stoul(current_field);
        iss >> current_field;
        unsigned int strand_end = stoul(current_field);
        if (strand_start > strand_end) swap(strand_start, strand_end);
        collision_linear_order.push_back(make_tuple(strand_id, strand_start, strand_end, primer_id));*/
    }
/*
    sort(collision_linear_order.begin(), collision_linear_order.end(), CollisionPositionCMP);

    for (int i = 1; i < collision_linear_order.size(); i++) {
        auto &cur_collision = collision_linear_order[i];
        StrandID strand_id = get<0>(cur_collision);
        unsigned long long start = get<1>(cur_collision);
        unsigned long long end = get<2>(cur_collision);
        
	start += strand_id*198;
       	end += strand_id*198;
	//cout<<start<<" "<<end<<" |";	
	PrimerID primer_id = get<3>(cur_collision);

        unsigned long long cut = start + (end - start) / 2;

        cut = (unsigned long long)(cut * 0.1 + 0.5);
        
        int j = i-1;
        while (j >= 0) {
            auto &last_collision = collision_linear_order[j];
            // if ()
            StrandID strand_id_last = get<0>(last_collision);
            //if (strand_id_last == strand_id) {
                unsigned long long start_last = get<1>(last_collision);
                unsigned long long end_last = get<2>(last_collision);
                start_last += strand_id_last*198;
         	end_last += strand_id_last*198;
		PrimerID primer_id_last = get<3>(last_collision);

                unsigned long long cut_last = start_last + (end_last - start_last)/2;
                cut_last = (unsigned long long)(cut_last * 0.1 + 0.5);
		//cout<<cut <<" "<<cut_last<<"\t";
                //assert(cut_last <= cut);
                if (cut - cut_last > blind_spot.size()-1) break;
		if (cut<=cut_last || IsBlindSpot(cut-cut_last)) {
		     if (primer_id_last == primer_id) {
                        discarded_primers.insert(primer_id);
                    } else {
                        primer_confilct_list[primer_id].insert(primer_id_last);
                        primer_confilct_list[primer_id_last].insert(primer_id);
                    }
                }

            //} else break;
            j--;
        }
	//cout<<endl;
    }
*/
    assert(strand_id2name.size() == strand_name2id.size());
    assert(primer_id2name.size() == primer_name2id.size());
    n_strand += strand_id2name.size();
    n_primer = all_primers.size();

}

void VariableLength::out_intermidum_result(string out_file_path) {
    /*fstream out;
    out.open(out_file_path,ios::out);

    out<<"collision_begin"<<endl;
    for (auto n:primer_collision_num_){
         out<<n.first<<" "<<n.second<<endl;
    }
    out<<"collision_end"<<endl;*/

    //new added
    int interval=1;
    ofstream primer_collision_num;
    primer_collision_num.open ("primer_collision_num.csv",ios::out | ios::trunc);
    vector<int> collision_num(15000,0);

    for (auto n:primer_collision_num_) {
        collision_num[n.second/interval]++;
    }

    for(int i=0; i < collision_num.size(); i++){
            // write into file
        primer_collision_num<<i*interval<<","<<collision_num[i]<<endl;
    }
    primer_collision_num.close();
    //end new added


 /*   out<<"conflict_begin"<<endl;
    for (auto m:primer_confilct_list){
        out<<m.first<<" ";
        for (auto mm:m.second){
            out<<mm<<" ";
        }
        out<<endl;
    }
    out<<"conflict_end"<<endl;
*/
}

void VariableLength::in_intermidum_result(string out_file_path) {

    vector<string> all_intermidum_results = listFiles(out_file_path, true);

    for (int i = 0; i < all_intermidum_results.size(); i++) {
        cout<<"read result:"<<all_intermidum_results[i]<<endl;
        fstream in;
        string line;
        in.open(all_intermidum_results[i],ios::in);
        bool collision_mode=false;
        bool conflict_mode=false;


        while(getline(in, line)) {
            if (line=="collision_begin") {collision_mode= true; continue;}
            if (line=="collision_end") {collision_mode= false; continue;}
            if (line=="conflict_begin") {conflict_mode=true; continue;}
            if (line=="conflict_end") {conflict_mode=false; continue;}

            if (collision_mode){
                istringstream iss(line);
                string current_field;
                iss >> current_field;
                int primer_id = stoi(current_field);
                iss >> current_field;
                int collision_num = stoi(current_field);

                if (primer_collision_num_.find(primer_id)==primer_collision_num_.end()){
                    primer_collision_num_.emplace(primer_id,collision_num);
                } else{
                    primer_collision_num_[primer_id]+=collision_num;
                }

            } else if (conflict_mode){
                istringstream iss(line);
                string current_field;
                iss >> current_field;
                int primer1 = stoi(current_field);

                unordered_set<PrimerID> conflicts;

                if (primer_confilct_list.find(primer1)==primer_confilct_list.end()){
                    primer_confilct_list.emplace(primer1,conflicts);
                }

                auto it = primer_confilct_list.find(primer1);

                string primer2;
                while (iss>>primer2){
                    it->second.insert(stoi(primer2));
                }

            } else cout<<"error in read intermidum result"<<endl;
        }
    }

}

bool sortbyfirst_asending(const pair<int,PrimerID> &a, const pair<int,PrimerID> &b){
    return a.first<b.first;
}

bool sortbyfirst_descending(const pair<int,PrimerID> &a, const pair<int,PrimerID> &b){
    return a.first>b.first;
}

void VariableLength::Cut() {
    int total_loss=0;

    if (g_out_intermedium_results){
        for (int i = 0; i < all_files.size(); i++) {
            ReadCollisions(all_files[i]);
        }
        int total_collided_primer = primer_collision_num_.size();
        out_intermidum_result(g_out_varlen_intermedium_result_path);
    }

    if (g_in_intermedium_results){
        in_intermidum_result(g_in_varlen_intermedium_result_path);

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

        //prescreen primers with too many collisions

        for (int i = 0; i < primer_process_order.size(); i++) {
            PrimerID current_primer_id = primer_process_order[i].second;

            if(discarded_primers.find(current_primer_id) != discarded_primers.end()) {
                continue;
            }
            if (primer_collision_num_.find(current_primer_id)==primer_collision_num_.end()){
                cout<<"weird"<<endl; //todo: why cannot find the primer 0

            }
            total_loss += primer_collision_num_.find(current_primer_id)->second;
            // recover current primer = push its conflict to discarded
            unordered_set<PrimerID> &conflicts = primer_confilct_list[current_primer_id];

            for (auto it = conflicts.begin(); it != conflicts.end(); it++){
                discarded_primers.insert(*it);
            }
        }

        cout<<"total cut collisions = "<<total_loss<<endl;
        PrintStatistics(primer_collision_num_.size());
    }

}

void VariableLength::PrintStatistics(int total_collided_primer) {
    cout << "n_strand: " << n_strand << endl << flush;
    cout << "n_primer: " << n_primer << endl << flush;
    int n_discarded = discarded_primers.size();
    int n_recovered = total_collided_primer - n_discarded;
    int n_free = 28000 - total_collided_primer;
    cout << "n[free, recovered, discarded] = " << n_free << ' ' << n_recovered << ' ' << n_discarded << endl << flush;
    cout << "Available primer ratio before VarLen = " << (double)n_free/n_primer << endl << flush;
    cout << "Available primer ratio after VarLen = " << (double)(n_free+n_recovered)/n_primer << endl << flush;
}
