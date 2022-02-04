#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <typeinfo>
#include <cstring>
#include <vector>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <climits>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <bitset>
#include <sstream>
#include <cstdlib>
#include <cassert>

using namespace std;

typedef vector<int> com;
typedef vector<com> comList;

class dense_area{
public:
    dense_area(){
        mapping = -1;
    }
    dense_area(unsigned int s, unsigned int e){
        start = s;
        end = e;
        mapping = -1; // means the mapping has not been decided
    };

    // dense_area(const dense_area& area){
    //     start = area.start;
    //     end = area.end;
    //     mapping = area.mapping;

    //     for(int i=0;i<4;i++){
    //         vector<pair<unsigned int,string>> init1 = area.mappings[i];
    //         vector<unsigned int> init2 = area.cut_point_lists[i];
    //         bool init3 = area.legal[i];
    //         mappings.push_back(init1);
    //         cut_point_lists.push_back(init2);
    //         legal.push_back(init3);
    //     }
    // }

    // dense_area operator= (dense_area const &area){
    //     cout<<"=operator"<<endl<<flush;
    //     dense_area res;
    //     res.start = area.start;
    //     res.end = area.end;
    //     res.mapping = area.mapping;

    //     for(int i=0;i<4;i++){
    //         vector<pair<unsigned int,string>> init1 = area.mappings[i];
    //         vector<unsigned int> init2 = area.cut_point_lists[i];
    //         bool init3 = area.legal[i];
    //         res.mappings.push_back(init1);
    //         res.cut_point_lists.push_back(init2);
    //         res.legal.push_back(init3);
    //     }
    //     cout<<res.cut_point_lists.size()<<endl<<flush;
    //     return res;
    // }


    unsigned int start;
    unsigned int end;
    int mapping;
    // -1: not decided; 0-3: mapping 0-3

    // pair: collision position and corresponding primer
    // vector: 4 mappings, each mapping has a list of pairs
    vector<vector<pair<unsigned int,string>>> mappings;
    
    // 4 mappings, each mapping has a ordered list of cut points, excluding the start and end
    vector<vector<unsigned int>> cut_point_lists;

    vector<bool> legal;
    // used in mapping selection
    // indicating each mapping is banned or not

    vector<vector<bool>> legal_combinations;
    comList avaiable_comList;

};
vector<dense_area> dense_area_list;
vector<unsigned int> area_starting_pos_list;
vector<bool> blind_spot;
//map<pair<unsigned int, unsigned int>,string> collision_list;// vector
vector<pair<pair<unsigned int, unsigned int>,string>> collision_list;
vector<pair<unsigned int, unsigned int>> hybrid_area;
vector<pair<unsigned int, unsigned int>> prev_hybrid_area;
unsigned int total_len = 0;
vector<string> all_files[4];
int strand_length_list[4] = {150,160,190,200}; 
unordered_map<string,int> collided_primers;
vector<unordered_map<string,int>> collided_primers_mapping;

string output_file_name = "collided_primers_hybrid_sol.txt";
string output_file_mapping[4] = {"collided_primers_with_mapping_0.txt","collided_primers_with_mapping_1.txt","collided_primers_with_mapping_2.txt","collided_primers_with_mapping_3.txt"};



unordered_map<int, comList> len_to_comList;

comList cur_comList;
com cur_com;

void construct_length_cut(int left_len){
    if (left_len==0){
        if (cur_com.size()!=0){
            cur_comList.push_back(cur_com);
        }
        return;
    }
    if (left_len<0){
        return;
    }

    for (int i = 0;i<4;i++){
        cur_com.push_back(strand_length_list[i]);
        construct_length_cut(left_len-strand_length_list[i]);
        cur_com.pop_back();
    }
}

bool isDir(string dir)
{
    struct stat fileInfo;
    stat(dir.c_str(), &fileInfo);
    if (S_ISDIR(fileInfo.st_mode)) {
        return true;
    } else {
        return false;
    }
}

vector<string> listFiles(string baseDir, bool recursive)
{
    vector<string> all_files;
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(baseDir.c_str())) == NULL) {
        cout << "[ERROR: " << errno << " ] Couldn't open " << baseDir << "." << endl;
        return all_files;
    } else {
        while ((dirp = readdir(dp)) != NULL) {
            if (dirp->d_name[0] != '.') {
                if (isDir(baseDir + dirp->d_name) == true && recursive == true) {
                    //all_files_.push_back(baseDir + dirp->d_name);
                    vector<string> tmp =listFiles(baseDir + dirp->d_name + "/", true);
                    for(auto n:tmp){
                        all_files.push_back(n);
                    }
                } else {
                    all_files.push_back(baseDir + dirp->d_name);
                }
            }
        }
        closedir(dp);
    }
    return all_files;
}

bool collision_sort(pair<pair<unsigned int, unsigned int>,string> i, pair<pair<unsigned int, unsigned int>,string> j){
    auto a = i.first.first;
    auto b = j.first.first;
    if (a<b){return true;}
    else if (a>b){return false;}
    else{
        a = i.first.second;
        b = j.first.second;
        return (a<b);
    }
}

bool dense_area_sort(const dense_area& i, const dense_area& j){
    auto len1 = i.end-i.start;
    auto len2 = j.end-j.start;
    if(len1 != len2) {
        return (len1<len2);
    }

    if (i.start!=j.start){
        return (i.start<j.start);
    }
    return (i.end<j.end);
    
}

bool collided_list_sort(pair<unsigned int,string> i ,pair<unsigned int,string> j){
    if(i.first!=j.first){
        return (i.first<j.first);
    }

    return (i.second.compare(j.second)<0);
}

bool primer_sort(pair<string,vector<pair<pair<int,int>, unsigned int>>> i, pair<string,vector<pair<pair<int,int>, unsigned int>>> j){
    auto isize = i.second.size();
    auto jsize = j.second.size();
    return (isize<jsize);
    
}

void init_blind_spot(){
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
};

// retrun true if the distance is within the blind spot list
// return false otherwise
bool check_blind_spot(int dis){
    if (dis>=450){return false;}
    return blind_spot[(dis/10)];

    // for(auto i=0;i<blind_spot.size();i++){
    //     if (dis==blind_spot[i])return true;
    // }
    // return false;
};

// cut==0 means there is no appropriate cut position
// Update: we now ignore the instances of cut==0, so the trivial cut pos would be return instead
unsigned int find_cut_pos(unsigned int start, unsigned int end){
    unsigned int cut = start/10 + 1;  
    cut = cut*10;

    while (cut<end){
        int left = cut - start + 1;
        int right = end - cut;
        if (left>=12 || right >=12){
            cut+=10;
            continue;
        }
        return cut;
    }
    return 10*(start/10+1);
    //return 0;

};

// find a trivial cut_pos for those intrinsic collision
unsigned int find_trivial_cut_pos(unsigned int start, unsigned int end){
    
    return 10*(start/10+1);

};

// starting from the "end" position,
// find a position which is legal from "start" and is also legal to "area_end"
// i.e. start->pos->area_end both the distance should be legal
unsigned int check_until_not_blind_spot(unsigned int start, unsigned int end, unsigned int area_end){
    int dis = end-start;
    for (auto i=0;i<16;++i){
        dis += 10;
        if (!check_blind_spot(dis) && dis+start<=area_end && !check_blind_spot(area_end-dis-start)){
            return (start+dis);
        } 
    }
    return area_end;
}

void hybrid_area_append(unsigned int start, unsigned int end){
    hybrid_area.push_back(make_pair(start,end));
}

// using the cut position to find the dense_area the cut is in
// then return the end pos of this area
unsigned int find_end_pos_of_area(unsigned int cut){
    auto low = lower_bound(area_starting_pos_list.begin(),area_starting_pos_list.end(),cut);
    unsigned int index;
    if(*low==cut){
        index = (low-area_starting_pos_list.begin());
        
    }
    else {
        index = (low-area_starting_pos_list.begin())-1;
    }
    if(index<prev_hybrid_area.size()){
        auto end = prev_hybrid_area[index].second;
        return end;
    }
    else{
        return total_len;
    }
}

// starting from the "cut" position, iterate backwards
// find a position which is legal from "cur_pos" and is also legal to "end_pos_of_area"
unsigned int find_prev_cut_pos(unsigned int cur_pos,unsigned int cut,unsigned int end_pos_of_area){
    unsigned int res=cur_pos;
    for (auto i = cut;i>=cur_pos;i-=10){
        if(!check_blind_spot(i-cur_pos) && !check_blind_spot(end_pos_of_area-i)){
            res = i;
            break;
        }
    }
    return res;
}

// in-place update the hybrid area
void seperate_hybrid_area(){
    prev_hybrid_area = hybrid_area;
    hybrid_area.clear();
    unsigned int cur_pos = 0;

    // illegal == 0 means the current collision is safe to be cut
    int illegal = 0;

    unsigned int collision_count_in_hybrid_area = 0;
    unsigned int collision_count_in_sparse_area = 0;

    // have to make sure that cur_pos is always legal to the end of nearest hybrid area
    for(auto i=0;i<collision_list.size();++i){
        auto start = collision_list[i].first.first;
        auto end = collision_list[i].first.second;
        auto primer = collision_list[i].second;

        auto cut = find_cut_pos(start,end);
        //cout<<cut<<endl;
        int prev_illegal = illegal;
        illegal = 0;
        if (cut==0){
            cut = find_trivial_cut_pos(start,end);
            illegal = 1;
            // means the next cut area should be hybrid
        }

        unsigned int end_pos_of_area = find_end_pos_of_area(cut);
        // if (!check_blind_spot(end_pos_of_area-cut)){
        //     cur_pos = cut;
        // }

        if (!check_blind_spot(cut-cur_pos) && end_pos_of_area>= cut && !check_blind_spot(end_pos_of_area-cut)){
            if(illegal==0){
                collision_count_in_sparse_area++;
            }
            else{
                collision_count_in_hybrid_area++;
            }

            if (prev_illegal==1){
                hybrid_area_append(cur_pos,cut);
            }

            cur_pos = cut;          
        }
        else{ // distance is illegal (either cut_pos to this collision or this collision to end of hybrid area is illegal)
                // then check every point between cut and next_cut
                // such that cut to this point and this point to next_cut are both legal
                // if there is one, cut it there and set it as the right boundary of hybrid area
                // if not, check the next collision,
                // until we can find one

            // maybe find a previous cut pos
            auto prev_cut_pos = find_prev_cut_pos(cur_pos,cut,end_pos_of_area);
            if(prev_cut_pos!=cur_pos){
                if(illegal==0){
                    collision_count_in_sparse_area++;
                }
                else{
                    collision_count_in_hybrid_area++;
                }

                if(prev_illegal==1){
                    auto right = check_until_not_blind_spot(cur_pos,cur_pos, prev_cut_pos);
                    hybrid_area_append(cur_pos,right);
                }
                cur_pos = prev_cut_pos;
            }



            //collision_count_in_hybrid_area++;
            auto left = cur_pos;
            bool flag = false;
            while(!flag){
                if (i==collision_list.size()-1){
                    auto start = collision_list[i].first.first;
                    auto end = collision_list[i].first.second;
                    cut = find_trivial_cut_pos(start,end);

                    auto right = check_until_not_blind_spot(cur_pos,cut,end_pos_of_area);
                    //cout<<cur_pos<<" "<<cut<<" "<<end_pos_of_area<<endl;
                    hybrid_area_append(left,right);
                    flag = true;
                    illegal = 0;
                }
                else{
                    auto prev_cut = cut;
                    auto next_start = collision_list[i+1].first.first;
                    auto next_end = collision_list[i+1].first.second;
                    cut = find_cut_pos(next_start,next_end);
                    
                    illegal = 0;
                    if (cut==0){
                        cut = find_trivial_cut_pos(next_start,next_end);
                        illegal = 1;
                    }

                    unsigned int next_end_pos_of_area = find_end_pos_of_area(cut);
                    if(next_end_pos_of_area != end_pos_of_area){ // which means the next collision is in different hybrid area
                        auto right = check_until_not_blind_spot(cur_pos,prev_cut,end_pos_of_area);
                        //auto right = next_end_pos_of_area;
                        hybrid_area_append(left,right);
                        cur_pos = right;
                        flag = true;
                        illegal = 0;
                        break;
                    }

                    for (auto it = prev_cut+10;it<=cut;it+=10){
                        if(!check_blind_spot(it-cur_pos)  && !check_blind_spot(cut-it)  && end_pos_of_area>=cut &&!check_blind_spot(end_pos_of_area-cut)){
                            hybrid_area_append(cur_pos,it);
                            cur_pos = cut;
                            if(illegal){
                                collision_count_in_hybrid_area++;
                            }
                            else{
                                collision_count_in_sparse_area++;
                            }
                            flag = true;
                            break;
                        }
                    }
                    i++;
                }
                collision_count_in_hybrid_area++;
            }
        }
        if (illegal!=0 && i<collision_list.size()-1){
            auto next_start = collision_list[i+1].first.first;
            auto next_end = collision_list[i+1].first.second;
            cut = find_cut_pos(next_start,next_end);  

            if (cut==0){
                cut = find_trivial_cut_pos(next_start,next_end);
                // legal = 2;
            }

            for (auto it = 0;it<4;it++){
                auto new_pos = cur_pos+strand_length_list[it];
                if (new_pos<=cut && !check_blind_spot(cut-new_pos) && new_pos<=end_pos_of_area &&!check_blind_spot(end_pos_of_area-new_pos)){
                    hybrid_area_append(cur_pos,new_pos);
                    cur_pos = new_pos;
                    illegal = 0;
                    break;
                }
            }
        }
        if(illegal!=0 && i==collision_list.size()-1){
            auto left = cut;
            auto right = check_until_not_blind_spot(cur_pos,cut,end_pos_of_area);
            hybrid_area_append(left,right);
        }

    }
    // cout<<hybrid_area.size()<<endl;
    // cout<<"start combine"<<endl;
    // combine the adjacent hybrid area
    area_starting_pos_list.clear();
    vector<pair<unsigned int, unsigned int>> new_hybrid_area;
    new_hybrid_area.clear();
    for (auto i=0;i<hybrid_area.size();i++){
        auto left = hybrid_area[i].first;
        auto right = hybrid_area[i].second;
        while((i+1)!=hybrid_area.size() && right==hybrid_area[i+1].first){
            right = hybrid_area[i+1].second;
            i++;
        }
        new_hybrid_area.push_back(make_pair(left,right));
        area_starting_pos_list.push_back(left);
    }
    hybrid_area.clear();
    hybrid_area = new_hybrid_area;
    // for (auto it = hybrid_area.begin();it<hybrid_area.end();it++){
    //     if ((it+1)!=hybrid_area.end()){
    //         if (it->second==(it+1)->first){
    //             it->second = (it+1)->second;
    //             hybrid_area.erase(it+1);
    //         }
    //     }
    // }


    // print 
    //unsigned int total_len = collision_list[collision_list.size()-1].first.second;
    unsigned int hybrid_len = 0; 
    //cout<<"Existing hybrid area"<<endl;
    //int not_150_count = 0; 
    for (auto it = hybrid_area.begin();it<hybrid_area.end();it++){
        //cout<<it->first<<" "<<it->second<<endl;
        //if (it->second-it->first != 150){not_150_count++;}
        hybrid_len += (it->second - it->first);
    }

    // calculate hybrid area before expansion
    unsigned int dense_len_before_expan = 0;
    for(auto i=0;i<collision_list.size();i++){
        auto start = collision_list[i].first.first;
        auto end = collision_list[i].first.second;
        auto cut = find_cut_pos(start,end);

        if (cut==0){
            dense_len_before_expan += 10;
            cut = find_trivial_cut_pos(start,end);
        }
        if (i!=collision_list.size()-1){
            auto next_start = collision_list[i+1].first.first;
            auto next_end = collision_list[i+1].first.second;
            auto next_cut = find_cut_pos(next_start,next_end);
            if (next_cut==0){
                next_cut = find_trivial_cut_pos(next_start,next_end);
            }
            if (check_blind_spot(next_cut-cut)){
                dense_len_before_expan += (next_cut-cut);
            }
        }
    }

    unsigned int num_of_dummy_collision = 0;
    for(auto i = 0;i<collision_list.size();i++){
        string primer = collision_list[i].second;
        string dummy = "d";
        if(primer.compare(dummy)==0){
            num_of_dummy_collision+=1;
        }
    }

    int count_smaller_than_200 = 0;
    for (auto i = 0;i<hybrid_area.size();i++){
        if (hybrid_area[i].second-hybrid_area[i].first <= 200){
            count_smaller_than_200++;
        }
    }
    

    cout<<">>>>>>>>>>>information output>>>>>>>>>>>"<<endl;\
    cout<<"num of hybrid area whose size is smaller than 200: "<<count_smaller_than_200<<endl;
    cout<<"num of hybrid area: "<<hybrid_area.size()<<endl;
    // for(auto i = 0;i<hybrid_area.size();i++){
    //     cout<<hybrid_area[i].first<<" "<<hybrid_area[i].second<<endl;
    // }
    
    //check dense area before expansion
    //cout<<"dense area before expansion: "<<dense_len_before_expan<<endl;
    //cout<<"dense area after expansion: "<<hybrid_len<<endl;
    //cout<<">>>>>>>>>>>"<<endl;
    cout<<"hybrid length: "<<hybrid_len<<endl;
    cout<<"total length: "<<total_len<<endl;
    cout<<"percentage of hybrid area of total len: "<<(double)hybrid_len/total_len<<endl;
    //cout<<">>>>>>>>>>>"<<endl;
    //cout<<"collision_count_in_hybrid_area: "<<collision_count_in_hybrid_area<<endl;
    //cout<<"collision_count_in_hybrid_area: "<<collision_list.size()-collision_count_in_sparse_area<<endl;
    //cout<<"collision_count_in_sparse_area: "<<collision_count_in_sparse_area-num_of_dummy_collision<<endl;
    cout<<"total collision count: "<<collision_list.size()-num_of_dummy_collision<<endl;
    
    //cout<<"num of dummy collisions: "<<num_of_dummy_collision<<endl;
    cout<<">>>>>>>>>>>"<<endl;
    cout<<endl;


    // insert dummy collisions
    // for(auto i=0;i<hybrid_area.size();i++){
    //     if(hybrid_area[i].first==0){
    //         continue;
    //     }
    //     unsigned int start = hybrid_area[i].first - 5;
    //     unsigned int end = hybrid_area[i].first + 5;
    //     string dummy_primer = "0";
    //     collision_list.push_back(make_pair(make_pair(start,end),dummy_primer));
    // }
}

// return the index id of hybrid area which the cut is in
int find_hybrid_area(unsigned int cut, unsigned int trivial_cut){
    auto low = lower_bound(area_starting_pos_list.begin(),area_starting_pos_list.end(),trivial_cut);
    if(*low==trivial_cut){return (low-area_starting_pos_list.begin());}
    return (low-area_starting_pos_list.begin())-1;
}

// check whether the cut is in any hybrid area
bool is_in_hybrid_area(unsigned int cut, unsigned int trivial_cut){

    //if(cut!=0){trivial_cut = cut;}
    //cout<<trivial_cut<<" "<<cut<<endl;
    //if(cut!=0 && cut!=trivial_cut){cout<<"!!!"<<endl;}

    if(hybrid_area.size()==0){return true;}

    int area_index = 0;
    // binary search in hybrid_area
    int start = 0;
    int end = hybrid_area.size()-1;

    while(start<=end){
        int mid = (start+end)/2;

        if(hybrid_area[mid].first==trivial_cut){
            if(trivial_cut==cut){
                return false;
            }
            else{
                return true;
            }
        }
        else if (hybrid_area[mid].first<trivial_cut){
            start = mid + 1;
        }
        else{
            end = mid - 1;
        }
    }
    area_index = end + 1;

    if(area_index==0){
        return false;
    }
    else{
        auto area_start = hybrid_area[area_index-1].first;
        auto area_end = hybrid_area[area_index-1].second;
        if(trivial_cut>area_start && trivial_cut<area_end){
            return true;
        }
        else{
            return false;
        }
    }
}

// read the blast result file and collect the collisions inside the hybrid areas
void construct_collision_list(int file_index, string file_name){
    if(file_index < 1 || file_index > 4){
        cerr << "wrong file index" << endl;
        return;
    }

    collision_list.clear();

    ifstream myfile;
    myfile.open (file_name);   
    string line;

    int count =0;
    while (getline(myfile,line)){
        istringstream iss(line);
        string a;
        iss>>a;
        
        if (!strcmp(a.c_str(),"#")){
            continue;
        }
        if(a.size()<=6 || a.substr(0,6)!="primer"){
            continue;
        }
        string primer = a.substr(6);
        for (auto i=0;i<8;i++){
            iss>>a;
        }
        unsigned int start = stoul(a);
        iss>>a;
        unsigned int end = stoul(a);
        if (start>end){
            unsigned int temp = end;
            end = start;
            start = temp;
        }

        // used for debug
        // count++;
        // if(count>=10000){break;}
        
        auto cut = find_cut_pos(start,end);
        auto trivial_cut = find_trivial_cut_pos(start,end);
        if(is_in_hybrid_area(cut,trivial_cut)){
            collision_list.push_back(make_pair(make_pair(start,end),primer));
        }
        // count++;
        // if(count>=100){break;}
        
        //collision_list.push_back(make_pair(make_pair(start,end),primer));

        
    }
    // insert dummy collisions
    // dummy collisions are inserted to solved the beginning boundary problem of each dense area
    for(auto i=0;i<hybrid_area.size();i++){
        if(hybrid_area[i].first==0){
            continue;
        }
        unsigned int start = hybrid_area[i].first - 5;
        unsigned int end = hybrid_area[i].first + 5;
        string dummy_primer = "d";
        collision_list.push_back(make_pair(make_pair(start,end),dummy_primer));
    }

    sort(collision_list.begin(),collision_list.end(),collision_sort);
    unsigned int len = collision_list[collision_list.size()-1].first.second;
    len = 10*(len/10+45);
    if(file_index==1){
        
        hybrid_area.push_back(make_pair(0,len));
        area_starting_pos_list.push_back(0);
        total_len = len;
        
    }

    //if(len>total_len){total_len = len;}
    // if (file_index==1){
    //     total_len = collision_list[collision_list.size()-1].first.second;
    // }

    // cout<<"start\tend\tprimer\tuncuttable_or_not\tcut_pos"<<endl;
    // for (auto it=collision_list.begin();it!=collision_list.end();it++){
    //     cout<<it->first.first<<" "<<it->first.second<<" "<<it->second<<" ";
    //     auto cut = find_cut_pos(it->first.first,it->first.second);
    //     if (cut == 0 ){cout<<" y ";cut = find_trivial_cut_pos(it->first.first,it->first.second);}
    //     else{cout<<" n ";}
    //     cout<<cut<<endl;
    // }
    // cout<<">>>>>>>>>>>>>>>"<<endl;

}

// Update: this function is not used anymore
void eliminate_collision_in_cut_area(){
    vector<pair<pair<unsigned int, unsigned int>,string>> new_collision_list;
    for (auto i=0;i<collision_list.size();i++){
        auto start = collision_list[i].first.first;
        auto end = collision_list[i].first.second;
        auto cut = find_cut_pos(start,end);
        auto trivial_cut = find_trivial_cut_pos(start,end);

        if(is_in_hybrid_area(cut,trivial_cut)){
        //if(1){
            new_collision_list.push_back(collision_list[i]);
        }
    }
    collision_list.clear();
    collision_list = new_collision_list;
}

// initial the data structure named "dense_area_list"
void initial_dense_area(){
    dense_area_list.clear();
    for(auto i = 0;i<hybrid_area.size();i++){
        dense_area area(hybrid_area[i].first,hybrid_area[i].second);
        
        vector<pair<unsigned int,string>> temp;
        vector<unsigned int> temp2;
        for(int i=0;i<4;i++){
            area.mappings.push_back(temp);
            area.cut_point_lists.push_back(temp2);
            area.legal.push_back(true);
        }

        
        //cout<<dense_area_list[i].mappings.size()<<endl;    

        area.legal_combinations.clear();
        auto len = hybrid_area[i].second - hybrid_area[i].first;
        //cout<<len<<" "<<len_to_comList[len].size()<<endl;
        auto cur_comList = len_to_comList[len];
        for (int j=0;j<cur_comList.size();j++){
            int num_of_pieces = cur_comList[j].size();
            int num = pow(4,num_of_pieces);
            vector<bool> temp (num,true);
            area.legal_combinations.push_back(temp);    
        }
        area.avaiable_comList = len_to_comList[len];
        //cout<<len_to_comList[len].size()<<endl;
        // cout<<"init: "<<area.avaiable_comList.size()<<endl;

        dense_area_list.push_back(area);
    }
    // cout<<"first area: "<<dense_area_list[0].avaiable_comList.size()<<endl;
    // cout<<"fiftheenth area: "<<dense_area_list[15].avaiable_comList.size()<<endl;
}

// read the blast result file again and collect the collisions inside the hybrid areas
void construct_dense_area(int file_index, string file_name){
    if(file_index < 1 || file_index > 4){
        cerr << "wrong file index" << endl;
        return;
    }
    
    ifstream myfile;
    myfile.open (file_name);   
    string line;

    int count = 0;
    while (getline(myfile,line)){
        

        istringstream iss(line);
        string a;
        iss>>a;
        
        if (!strcmp(a.c_str(),"#")){
            continue;
        }
        if(a.size()<=6 || a.substr(0,6)!="primer"){
            continue;
        }
        string primer = a.substr(6);
        for (auto i=0;i<8;i++){
            iss>>a;
        }
        unsigned int start = stoul(a);
        iss>>a;
        unsigned int end = stoul(a);
        if (start>end){
            unsigned int temp = end;
            end = start;
            start = temp;
        }

        // used for debug
        // count++;
        // if(count>=10000){break;}
        
        auto cut = find_cut_pos(start,end);
        auto trivial_cut = find_trivial_cut_pos(start,end);
        if(is_in_hybrid_area(cut,trivial_cut)){
            int index = find_hybrid_area(cut,trivial_cut);
            //cout<<trivial_cut<<" "<<hybrid_area[index].first<<endl;
            dense_area_list[index].mappings[file_index-1].push_back(make_pair(trivial_cut,primer));
        }       
    }

}

void init_len_to_comList(){
    len_to_comList.clear();

    for (int m = 0; m<=1500; m+=10){
        cur_comList.clear();
        cur_com.clear();
        construct_length_cut(m);
        // for (int i = 0; i<cur_comList.size();i++){
        //     for (int j=0;j<cur_comList[i].size();j++){
        //         cout<<cur_comList[i][j]<<" ";
        //     }
        //     cout<<endl;
        // }

        len_to_comList[m] = cur_comList;
        //cout<<"len = "<<m<<"; has number of combinations of: "<<cur_comList.size()<<endl;
    }
}

int main(int argc, char** argv) {

    if (argc != 5) {
        cerr << "argc must be 5" << endl;
        return -1;
    }

    for (int i = 0; i < 4; i++) {
        all_files[i] = listFiles(argv[i+1], true);
        sort(all_files[i].begin(), all_files[i].end());
    }

    init_len_to_comList();
    //cout<<"len: "<<len_to_comList[1000].size()<<endl;

    // read collided primers from a txt file
    string line;
    ifstream myfile (output_file_name);
    if (myfile.is_open()){
        while ( getline (myfile,line) ){
            if (line!=""){
                collided_primers[line] = 1;
            }
            
        }
        myfile.close();
    }

    for (int i=0;i<4;i++){
        unordered_map<string,int> init;
    
        string line;
        ifstream myfile (output_file_mapping[i]);
        if (myfile.is_open()){
            while ( getline (myfile,line) ){
                if (line!=""){
                    init[line] = 1;
                }
            }
            myfile.close();
        }
        collided_primers_mapping.push_back(init);
    }

    // output original collided primers
    cout<<"Num of original collided primers for hybrid solution: "<<collided_primers.size()<<endl;
    cout<<"Num of original collided primers for each mapping: "<<collided_primers_mapping[0].size()<<" "<<collided_primers_mapping[1].size()<<" "<<collided_primers_mapping[2].size()<<" "<<collided_primers_mapping[3].size()<<endl;
    cout<<endl;

    for (int n = 0; n<all_files[0].size();n++){
        cout<<all_files[0][n]<<endl;
        init_blind_spot();
        collision_list.clear();
        hybrid_area.clear();
        area_starting_pos_list.clear();
        total_len = 0;

        for (int i = 0;i<4;i++){

            construct_collision_list(i+1,all_files[i][n]);


            seperate_hybrid_area();

        }
        //reconstruct collision list with 4 mappings
        initial_dense_area();

        for (int i = 0;i<4;i++){
            construct_dense_area(i+1,all_files[i][n]);
        }
        // sort dense area base on its length
        if(dense_area_list.size()!=0){
            sort(dense_area_list.begin(),dense_area_list.end(),dense_area_sort);
        }


        // dense_area information output 

        cout<<"----------dense area information---------------"<<endl;
        int num_of_area_whose_len_smaller_than_1000 = 0;
        int num_of_area_whose_len_smaller_than_1500 = 0;
        int num_of_area_whose_len_smaller_than_2000 = 0;
        int num_of_area_whose_len_smaller_than_3000 = 0;
        int num_of_area_whose_len_smaller_than_5000 = 0;
        int others = 0;
        for (int j = 0;j<dense_area_list.size();j++){
            auto len = dense_area_list[j].end - dense_area_list[j].start;
            if (len<=1000){
                num_of_area_whose_len_smaller_than_1000++;
            }
            else if (len<=1500){
                num_of_area_whose_len_smaller_than_1500++;
            }
            else if (len<=2000){
                num_of_area_whose_len_smaller_than_2000++;
            }
            else if (len<=3000){
                num_of_area_whose_len_smaller_than_3000++;
            }
            else if (len<=5000){
                num_of_area_whose_len_smaller_than_5000++;
            }
            else{
                others++;
            }
        }
        cout<<"num_of_area_whose_len_smaller_than_1000: "<<num_of_area_whose_len_smaller_than_1000<<endl;
        cout<<"num_of_area_whose_len_smaller_than_1500: "<<num_of_area_whose_len_smaller_than_1500<<endl;
        cout<<"num_of_area_whose_len_smaller_than_2000: "<<num_of_area_whose_len_smaller_than_2000<<endl;
        cout<<"num_of_area_whose_len_smaller_than_3000: "<<num_of_area_whose_len_smaller_than_3000<<endl;
        cout<<"num_of_area_whose_len_smaller_than_5000: "<<num_of_area_whose_len_smaller_than_5000<<endl;
        cout<<"num_of_area_whose_len_larger_than_5000: "<<others<<endl;



        // compare with four baselines which only selection one particular mapping

        // check the collided primer before mapping selection
        for (auto i = 0;i<4;i++){
            //unordered_map<string,int> sec_temp_collided_primers=collided_primers_mapping[i];
            for (auto j = 0;j<dense_area_list.size();j++){
                auto collided_list = dense_area_list[j].mappings[i];
                for(int k = 0;k<collided_list.size();k++){
                    auto primer = collided_list[k].second;
                    collided_primers_mapping[i][primer] = 1;
                }
            }
            cout<<"for mapping "<<i<<", original num of collided primers is: "<<collided_primers_mapping[i].size()<<endl;
        }
        cout<<endl;




        
        // new mapping selection
        cout<<"Running new mapping algorithm..."<<endl;
        cout<<endl;

        vector<pair<string,vector<pair<pair<int,int>, unsigned int>>>> temp;
        unordered_map<string,int> primer_to_temp_index;
        // temp-> [string,[((dense_area_id,mapping_id), pos),(did,mid,pos),(did,mid,pos)]]


        // collect the information to initialize the data structure of "temp"


        for (auto i = 0;i<dense_area_list.size();i++){
            bool no_collision = false;
            for (int j = 0;j<4;j++){
                if(dense_area_list[i].mappings[j].size()==0){
                    dense_area_list[i].mapping = j;
                    no_collision = true;
                    break;
                }
            }
            if(no_collision){continue;}
            for (int j = 0;j<4;j++){
                auto collided_list = dense_area_list[i].mappings[j];
                for (auto k=0;k<collided_list.size();k++){
                    auto primer = collided_list[k].second;
                    auto pos = collided_list[k].first;

                    if(primer_to_temp_index.find(primer)!=primer_to_temp_index.end()){// find the primer
                        auto temp_index = primer_to_temp_index[primer];
                        temp[temp_index].second.push_back(make_pair(make_pair(i,j),pos));
                    }
                    else{// cannot find the primer
                        vector<pair<pair<int,int>, unsigned int>> init;
                        init.push_back(make_pair(make_pair(i,j),pos));
                        temp.push_back(make_pair(primer,init));
                        primer_to_temp_index[primer] = temp.size()-1;
                    }

                }
            }
        }
        
        // based on the size of second part of temp
        stable_sort(temp.begin(),temp.end(),primer_sort);

        // check if the size is ascending 
        cout<<"temp.size = "<<temp.size()<<endl;
        for (auto i=0;i<temp.size()-1;i++){
            //cout<<i<<" ";
            assert(temp[i].second.size()<=temp[i+1].second.size());
        }
        //cout<<"check"<<endl;


        // new updates 10-31
        // for each primer
        //      for each dense area
        //          list all cut+mapping combinations
        //

        unordered_map<string,int> delta_collided_primers;
        for (auto i=0;i<temp.size();i++){
            bool successful_recover = true;

            unordered_map<int,dense_area> delta_dense_area_list; // did -> modified_dense_area

            auto primer = temp[i].first;
            auto collided_did_mid_pos_list = temp[i].second;
            if (collided_primers.find(primer)!=collided_primers.end() || delta_collided_primers.find(primer)!=delta_collided_primers.end()){
                // dont need to handle the primer that has been given up
                cout<<"skip for primer "<<i<<" "<<endl;
                continue;
            }

            // "for dense area"
            for (auto j=0;j<collided_did_mid_pos_list.size();j++){
                auto did = collided_did_mid_pos_list[j].first.first;
                auto mid = collided_did_mid_pos_list[j].first.second;
                auto pos = collided_did_mid_pos_list[j].second;
                
                dense_area temp_dense_area;
                if(delta_dense_area_list.find(did)==delta_dense_area_list.end()){ // not in delta
                    temp_dense_area = dense_area_list[did];
                }
                else { // already in delta
                    temp_dense_area = delta_dense_area_list[did];
                }

                auto start = temp_dense_area.start;
                auto end = temp_dense_area.end;
                auto len = end-start;
                auto legal_combinations = temp_dense_area.legal_combinations;

                // tempororily added here
                if (len>1500){
                    break;
                }


                // find the pos of collision and ban all the conflicting combinations 
                
                auto cur_comList = temp_dense_area.avaiable_comList;
                
                for (int k=0;k<cur_comList.size();k++){
                    auto com = cur_comList[k]; // e.g. com = [200,200,200]
                    auto num = com.size(); // e.g. n = 3
                    auto bool_array = legal_combinations[k]; // e.g. bool_array = [False]*(4^3)
                    // bool array index -> mapping combinations:
                    // 0: 0 0 0 -> False
                    // 1: 0 0 1 -> False
                    // 2: 0 0 2
                    // 3: 0 0 3
                    // 4: 0 1 0
                    // ...
                    // 15: 0 3 3
                    // 16: 1 0 0
                    // ...
                    // 32: 2 0 0
                    // 33: 2 0 1
                    // 34: 2 0 2
                    // ...
                    // 63: 3 3 3 
                    
                    // TODO: what if cur_len is right at the collision pos
                    int cur_len = start;
                    int cur_index = 0;
                    bool collision_cut = false;
                    while (true){
                        cur_len += com[cur_index];
                        if (pos<cur_len){
                            break;
                        }
                        else if (pos==cur_len){ // collision can be removed by this com
                            collision_cut = true;
                            break;
                        }
                        cur_index++;
                    }
                    if (collision_cut){
                        continue;
                    }


                    // the collision would be in the (cur_index)th cutted part (starting from 0)
                    // e.g. pos=150 will return cur_index=0
                    // e.g. pos=350 will return cur_index=1
                    // e.g. pos=550 will return cur_index=2

                    // ban all the combinations where this part has the mapping equals to mid

                    for (int l = 0; l<bool_array.size();l++){
                        int divisor = pow(4,(num - cur_index - 1)); 
                        // e.g. 4^(3-0-1) = 4^2 = 16
                        // e.g. 4^(3-1-1) = 4^1 = 4
                        // e.g. 4^(3-2-1) = 4^0 = 1
                        int mapping = (l/divisor)%4; 
                        // e.g. 33/16 = 2; 2 % 4 = 2 means the first mapping is mapping2
                        // e.g. 33/4 = 8; 8 % 4 = 0 means the second mapping is mapping0
                        // e.g. 33/1 = 33; 33 % 4 = 1 means the third mapping is mapping1
                        if (mapping==mid){ // this mapping is the same as the collided one
                            bool_array[l]=false;
                        }
                    }

                    legal_combinations[k] = bool_array;

                    // check if there are still available combinations
                    int available_count = 0;
                    for (int l = 0; l<bool_array.size();l++){
                        if(bool_array[l]){
                            available_count++;
                        }
                    }
                    if (available_count==0){// remove this com from comlist
                        cur_comList.erase(cur_comList.begin()+k);
                        legal_combinations.erase(legal_combinations.begin()+k);
                        k--;
                    }
                }
                
                
                
                // check if thre are still available combinations
                if (cur_comList.size()==0){
                    successful_recover = false;
                    break;
                }

                // update the comlist and legal_combinations
                temp_dense_area.legal_combinations = legal_combinations;
                temp_dense_area.avaiable_comList = cur_comList;

                delta_dense_area_list[did] = temp_dense_area;

            }   


            if (successful_recover){
                cout<<"successful recover primer "<<i<<" "<<endl;
                for (auto it = delta_dense_area_list.begin();it!=delta_dense_area_list.end();it++){
                    auto did = it->first;
                    dense_area_list[did] = it->second;
                }
            }
            else{
                cout<<"fail to recover primer "<<i<<" "<<endl;
                collided_primers[primer] = 1;
            }         


        }

        // unordered_map<string,int> delta_collided_primers;
        // // gradually recover the primers/decide the mapping
        // for (auto i=0;i<temp.size();i++){
        //     bool successful_recover = true;

        //     unordered_map<int,dense_area> delta_dense_area_list; // did -> modified_dense_area
        //     auto primer = temp[i].first;
        //     auto collided_did_mid_pos_list = temp[i].second;
        //     if (collided_primers.find(primer)!=collided_primers.end() || delta_collided_primers.find(primer)!=delta_collided_primers.end()){
        //         // dont need to handle the primer that has been given up
        //         //cout<<"skip "<<i<<" ";
        //         continue;
        //     }

        //     for (auto j=0;j<collided_did_mid_pos_list.size();j++){
        //         // for each collision
        //         //      cut the collision
        //         //      if no -> ban the mapping
        //         //      if no -> give up the primer


        //         auto did = collided_did_mid_pos_list[j].first.first;
        //         auto mid = collided_did_mid_pos_list[j].first.second;
        //         auto pos = collided_did_mid_pos_list[j].second;


        //         dense_area temp_dense_area;
        //         if(delta_dense_area_list.find(did)==delta_dense_area_list.end()){ // not in delta
        //             temp_dense_area = dense_area_list[did];
        //         }
        //         else { // already in delta
        //             temp_dense_area = delta_dense_area_list[did];
        //         }

        //         auto start = temp_dense_area.start;
        //         auto end = temp_dense_area.end;
        //         auto cut_point_list = temp_dense_area.cut_point_lists[mid];


        //         // try cut method
        //         bool successful_cut = false;

        //         // // find the nearest cut points to collision position
        //         if (cut_point_list.size()!=0){
        //             auto low = lower_bound(cut_point_list.begin(),cut_point_list.end(),pos);
        //             int low_index = low - cut_point_list.begin();

        //             if (low_index==cut_point_list.size()){
        //                 // collision is behind all existing cut points

        //                 auto left = cut_point_list[cut_point_list.size()-1];
        //                 auto right = end;

        //                 if (!check_blind_spot(pos-left) && !check_blind_spot(right-pos)){
        //                     // insert pos to cut_point_list
        //                     cut_point_list.insert(low,pos);
        //                     successful_cut = true;
        //                 }

        //                 // else it can not be recovered by cut

        //             }
        //             else if (cut_point_list[low_index]==pos){
        //                 // if the collision is already in an existing cut point

        //                 // do nothing and move to next collision
        //                 continue;
        //             }
        //             else if (low_index == 0){
        //                 // collision is before all existing cut points
        //                 auto left = start;
        //                 auto right = cut_point_list[0];

        //                 if (!check_blind_spot(pos-left) && !check_blind_spot(right-pos)){
        //                     // insert pos to cut_point_list
        //                     cut_point_list.insert(low,pos);
        //                     successful_cut = true;
        //                 }

        //                 // else it can not be recovered by cut

        //             }
        //             else{
        //                 auto left = cut_point_list[low_index-1];
        //                 auto right = cut_point_list[low_index];
        //                 if (!check_blind_spot(pos-left) && !check_blind_spot(right-pos)){
        //                     // insert pos to cut_point_list
        //                     cut_point_list.insert(low,pos);
        //                     successful_cut = true;
        //                 }
        //             }
        //         }
        //         else{
        //             auto left = start;
        //             auto right = end;
        //             if (pos==left || pos==right){
        //                 continue;
        //             }
        //             else if (!check_blind_spot(pos-left) && !check_blind_spot(right-pos)){
        //                 cut_point_list.push_back(pos);
        //                 successful_cut = true;
        //             }
        //         }

        //         if (successful_cut){
        //             // update the cut point list
        //             temp_dense_area.cut_point_lists[mid] = cut_point_list;
        //             delta_dense_area_list[did] = temp_dense_area;
        //             continue;
        //         }

        //         // else try mapping selection method


        //         if (temp_dense_area.mapping==-1){// haven't decided yet
        //             // ban the corresponding mapping
        //             temp_dense_area.legal[mid] = false;

        //             int count_illegal = 0;
        //             int legal_id = 0;
        //             // if three of the mappings are illegal
        //             // then choose the left one
        //             for (int k=0;k<4;k++){
        //                 if (temp_dense_area.legal[k]==false){
        //                     count_illegal+=1;
        //                 }
        //                 else{
        //                     legal_id = k;
        //                 }
        //             }
        //             if (count_illegal >= 3){
        //                 temp_dense_area.mapping = legal_id;
        //             }

        //         }
        //         else if (temp_dense_area.mapping==mid){
        //             // mapping has been decided and is equal to the mapping id
        //             // then the recover fails
                    
        //             delta_collided_primers[primer] = 1;
        //             //temp_collided_primers[primer] = 1;

        //             // move to the next primer candidate and do not save the dense area list
        //             successful_recover = false;
        //             break;
        //         }

        //         // update the delta dense_area
        //         delta_dense_area_list[did] = temp_dense_area; 
        //     }


        //     if (successful_recover){
        //         for (auto it = delta_dense_area_list.begin();it!=delta_dense_area_list.end();it++){
        //             auto did = it->first;
        //             dense_area_list[did] = it->second;
        //         }
        //     }
        //     else{
        //         collided_primers[primer] = 1;
        //     }
        // }





        cout<<"num of collided primers: "<<collided_primers.size()<<endl;
        cout<<endl;

    }


    // write collided primers to a txt file
    ofstream output;
    output.open(output_file_name);
    if(output.is_open()){
        for (auto it = collided_primers.begin();it!=collided_primers.end();it++){
            output<<it->first<<endl;
        }
        output.close();
    }
    else{
        cout<<"unable to write to file"<<endl;
    }


    for (int i=0;i<4;i++){
        ofstream output;
        output.open(output_file_mapping[i]);
        if(output.is_open()){
            for (auto it = collided_primers_mapping[i].begin();it!=collided_primers_mapping[i].end();it++){
                output<<it->first<<endl;
            }
            output.close();
        }
        else{
            cout<<"unable to write to file"<<endl;
        }
    }

    return 0;
}

