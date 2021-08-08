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

using namespace std;


vector<bool> blind_spot;
//map<pair<unsigned int, unsigned int>,string> collision_list;// vector
vector<pair<pair<unsigned int, unsigned int>,string>> collision_list;
vector<pair<unsigned int, unsigned int>> hybrid_area;
unsigned int total_len;
int strand_length_list[4] = {150,160,190,200};

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
    return 0;

};

// find a trivial cut_pos for those intrinsic collision
unsigned int find_trivial_cut_pos(unsigned int start, unsigned int end){
    
    return 10*(start/10+1);

};

unsigned int check_until_not_blind_spot(unsigned int start, unsigned int end){
    int dis = end-start;
    for (auto i=0;i<16;++i){
        dis += 10;
        if (!check_blind_spot(dis)){
            return (start+dis);
        } 
    }
    return start+450;
}

void hybrid_area_append(unsigned int start, unsigned int end){
    hybrid_area.push_back(make_pair(start,end));
}

// in-place update the hybrid area
void seperate_hybrid_area(){
    hybrid_area.clear();
    unsigned int cur_pos = 0;

    // illegal == 0 means the current collision is safe to be cut
    int illegal = 0;

    unsigned int collision_count_in_hybrid_area = 0;
    unsigned int collision_count_in_sparse_area = 0;

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

        if (!check_blind_spot(cut-cur_pos)){
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
        else{ // distance is illegal
                // then check every point between cut and next_cut
                // such that cut to this point and this point to next_cut are both legal
                // if there is one, cut it there and set it as the right boundary of hybrid area
                // if not, check the next collision,
                // until we can find one


            //collision_count_in_hybrid_area++;
            auto left = cur_pos;
            bool flag = false;
            while(!flag){
                if (i==collision_list.size()-1){
                    auto right = check_until_not_blind_spot(cur_pos,cut);
                    hybrid_area_append(left,right);
                    flag = true;
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

                    for (auto it = prev_cut+10;it<=cut;it+=10){
                        if(!check_blind_spot(it-cur_pos)  && !check_blind_spot(cut-it)){
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
                if (cur_pos+strand_length_list[it]<=cut && !check_blind_spot(cut-cur_pos-strand_length_list[it])){
                    hybrid_area_append(cur_pos,cur_pos+strand_length_list[it]);
                    cur_pos = cur_pos+strand_length_list[it];
                    illegal = 0;
                    break;
                }
            }
        }
        if(illegal!=0 && i==collision_list.size()-1){
            auto left = cut;
            auto right = check_until_not_blind_spot(cur_pos,cut);
            hybrid_area_append(left,right);
        }

    }
    //cout<<"start combine"<<endl;
    // combine the adjacent hybrid area
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
    


    cout<<"num of hybrid area: "<<hybrid_area.size()<<endl;
    cout<<">>>>>>>>>>>information output>>>>>>>>>>>"<<endl;
    //check dense area before expansion
    cout<<"dense area before expansion: "<<dense_len_before_expan<<endl;
    cout<<"dense area after expansion: "<<hybrid_len<<endl;
    cout<<">>>>>>>>>>>"<<endl;
    cout<<"hybrid length: "<<hybrid_len<<endl;
    cout<<"total length: "<<total_len<<endl;
    cout<<"percentage of hybrid area of total len: "<<(double)hybrid_len/total_len<<endl;
    cout<<">>>>>>>>>>>"<<endl;
    cout<<"collision_count_in_hybrid_area: "<<collision_count_in_hybrid_area<<endl;
    //cout<<"collision_count_in_hybrid_area: "<<collision_list.size()-collision_count_in_sparse_area<<endl;
    cout<<"collision_count_in_sparse_area: "<<collision_count_in_sparse_area-num_of_dummy_collision<<endl;
    cout<<"total collision count: "<<collision_list.size()-num_of_dummy_collision<<endl;
    cout<<">>>>>>>>>>>"<<endl;
    cout<<"num of dummy collisions: "<<num_of_dummy_collision<<endl;
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

void construct_collision_list(int file_index, char** argv){
    if(file_index < 1 || file_index > 4){
        cerr << "wrong file index" << endl;
        return;
    }

    collision_list.clear();

    ifstream myfile;
    myfile.open (argv[file_index]);   
    string line;

    int count =0;
    while (getline(myfile,line)){
        

        istringstream iss(line);
        string a;
        iss>>a;
        
        if (!strcmp(a.c_str(),"#")){
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

    if(len>total_len){total_len = len;}
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

int main(int argc, char** argv) {

    if (argc != 5) {
        cerr << "argc must be 5" << endl;
        return -1;
    }
    init_blind_spot();
    collision_list.clear();
    hybrid_area.clear();
    total_len = 0;


    for (int i = 0;i<4;i++){
        construct_collision_list(i+1,argv);
        // if(i!=0){
        //     eliminate_collision_in_cut_area();
        // }
        seperate_hybrid_area();
    }
    

    // check sort function
    // cout<<"start\tend\tprimer\tuncuttable_or_not\tcut_pos"<<endl;
    // for (auto it=collision_list.begin();it!=collision_list.end();it++){
    //     cout<<it->first.first<<" "<<it->first.second<<" "<<it->second<<" ";
    //     auto cut = find_cut_pos(it->first.first,it->first.second);
    //     if (cut == 0 ){cout<<" y ";cut = find_trivial_cut_pos(it->first.first,it->first.second);}
    //     else{cout<<" n ";}
    //     cout<<cut<<endl;
    // }
    // cout<<">>>>>>>>>>>>>>>"<<endl;

    // check blind spot function
    // for (auto it=blind_spot.begin();it!=blind_spot.end();it++){
    //     cout<<*it<<endl;
    // }
    //cout<<check_blind_spot(450)<<endl;
    // cout<<find_trivial_cut_pos(79,100)<<endl;
    // cout<<find_trivial_cut_pos(80,100)<<endl;
    // cout<<find_trivial_cut_pos(81,100)<<endl;

    return 0;
}

