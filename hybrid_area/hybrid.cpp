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
#include <map>
#include <climits>

#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <bitset>
#include <sstream>
#include <cstdlib>
#include <cassert>

using namespace std;

class dense_area{
public:
    dense_area(unsigned int s, unsigned int e){
        start = s;
        end = e;
    };

    unsigned int start;
    unsigned int end;
    int mapping;

    // pair: collision position and corresponding primer
    // vector: 4 mappings, each mapping has a list of pairs
    vector<vector<pair<unsigned int,string>>> mappings;
    
};
vector<dense_area> dense_area_list;
vector<unsigned int> area_starting_pos_list;
vector<bool> blind_spot;
//map<pair<unsigned int, unsigned int>,string> collision_list;// vector
vector<pair<pair<unsigned int, unsigned int>,string>> collision_list;
vector<pair<unsigned int, unsigned int>> hybrid_area;
vector<pair<unsigned int, unsigned int>> prev_hybrid_area;
unsigned int total_len;
vector<string> all_files[4];
int strand_length_list[4] = {150,160,190,200}; 
map<string,int> collided_primers;


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

bool dense_area_sort(dense_area i, dense_area j){
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
    

    cout<<">>>>>>>>>>>information output>>>>>>>>>>>"<<endl;
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

int find_hybrid_area(unsigned int cut, unsigned int trivial_cut){
    auto low = lower_bound(area_starting_pos_list.begin(),area_starting_pos_list.end(),trivial_cut);
    if(*low==trivial_cut){return (low-area_starting_pos_list.begin());}
    return (low-area_starting_pos_list.begin())-1;
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


void initial_dense_area(){
    dense_area_list.clear();
    for(auto i = 0;i<hybrid_area.size();i++){
        dense_area area(hybrid_area[i].first,hybrid_area[i].second);
        
        vector<pair<unsigned int,string>> temp;
        for(int i=0;i<4;i++){
            area.mappings.push_back(temp);
        }

        dense_area_list.push_back(area);
        //cout<<dense_area_list[i].mappings.size()<<endl;        
    }
}

void construct_dense_area(int file_index, string file_name){
    if(file_index < 1 || file_index > 4){
        cerr << "wrong file index" << endl;
        return;
    }
    
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

        
        auto cut = find_cut_pos(start,end);
        auto trivial_cut = find_trivial_cut_pos(start,end);
        if(is_in_hybrid_area(cut,trivial_cut)){
            int index = find_hybrid_area(cut,trivial_cut);
            //cout<<trivial_cut<<" "<<hybrid_area[index].first<<endl;
            dense_area_list[index].mappings[file_index-1].push_back(make_pair(trivial_cut,primer));
        }       
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


    for (int n = 0; n<all_files[0].size();n++){
        cout<<all_files[0][n]<<endl;
        init_blind_spot();
        collision_list.clear();
        hybrid_area.clear();
        area_starting_pos_list.clear();
        total_len = 0;

        for (int i = 0;i<4;i++){

            construct_collision_list(i+1,all_files[i][n]);

            // cout<<"collision list: "<<endl;
            // for(auto i = 0;i<collision_list.size();i++){
            //     if(collision_list[i].first.first<3928000 || collision_list[i].first.first>3929500){continue;}
            //     cout<<collision_list[i].first.first<<" "<<collision_list[i].first.second<<endl;
            // }

            seperate_hybrid_area();

            // cout<<"hybrid area: "<<endl;
            // for(auto i = 0;i<hybrid_area.size();i++){
            //     if(hybrid_area[i].first<3928000 || hybrid_area[i].first>3929500){continue;}
            //     cout<<hybrid_area[i].first<<" "<<hybrid_area[i].second<<endl;
            // }
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
        // print information of dense area
        // for (auto i = 0;i<dense_area_list.size();i++){
        //     cout<<"dense area: "<<i<<endl;
        //     cout<<"start and end: "<<dense_area_list[i].start<<" "<<dense_area_list[i].end<<endl;
        //     for(int j=0;j<4;j++){
        //         cout<<"mapping "<<j<<endl;
        //         for(int k=0;k<dense_area_list[i].mappings[j].size();k++){
        //             cout<<"pos and primer: "<<dense_area_list[i].mappings[j][k].first<<" "<<dense_area_list[i].mappings[j][k].second<<endl;
        //         }
        //     }
        //     cout<<endl;
        // }

        // decide the mapping
        // collided_primers.clear();

        map<string,int> temp_collided_primers=collided_primers;
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
            vector<int> list_of_increased_primer;
            int min = INT_MAX;
            int min_index = 0;

            //find the mapping with least increased primer
            for (int j = 0;j<4;j++){
                int increased_primer = 0;
                for (int k = 0;k<dense_area_list[i].mappings[j].size();k++){
                    if(temp_collided_primers.find(dense_area_list[i].mappings[j][k].second)==temp_collided_primers.end()){
                        increased_primer += 1;
                    }
                }
                list_of_increased_primer.push_back(increased_primer);
                if(increased_primer<min){
                    min = increased_primer;
                    min_index = j;
                }
                // //used for debug
                // if(i==67){
                //     cout<<"i==61"<<endl;
                //     for (auto k=0;k<list_of_increased_primer.size();k++){
                //         cout<<list_of_increased_primer[k]<<" ";
                //     }
                //     cout<<endl;
                //     cout<<j<<" "<<increased_primer<<" "<<min_index<<endl;
                // }
            }

            dense_area_list[i].mapping = min_index;
            for(int k = 0;k<dense_area_list[i].mappings[min_index].size();k++){
                string primer = dense_area_list[i].mappings[min_index][k].second;
                temp_collided_primers[primer]=1;
            }
        }

        //map<string,int> P_abandon;
        auto final_dense_area_list = dense_area_list;
        final_dense_area_list.clear();

        // variable-length cut
        for (auto i = 0;i<dense_area_list.size();i++){
            //cout<<i<<" "<<dense_area_list[i].end-dense_area_list[i].start<<" "<<dense_area_list[i].mapping<<endl;
            int area_length = dense_area_list[i].end-dense_area_list[i].start;

            
            auto start = dense_area_list[i].start;
            auto end = dense_area_list[i].end;
            int mapping = dense_area_list[i].mapping;
            auto collided_list = dense_area_list[i].mappings[mapping];
            // for (auto j = 0;j<collided_list.size();j++){
            //     cout<<collided_list[j].first<<" "<<collided_list[j].second<<endl;
            // }

            if(area_length<=200){ // the area is mixed
                final_dense_area_list.push_back(dense_area_list[i]);
                for (auto j=0;j<collided_list.size();j++){
                    string primer = collided_list[j].second;
                    collided_primers[primer]=1; // mark it as abandoned
                }
                continue;
            }

            bool add_to_final = true;
            for (auto j = 0;j<collided_list.size();j++){
                string primer = collided_list[j].second;
                auto pos = collided_list[j].first;
                if (collided_primers.find(primer)==collided_primers.end() && !check_blind_spot(pos-start) && !check_blind_spot(end-pos)){
                    dense_area one(start,pos);
                    one.mapping = mapping;
                    dense_area two(pos,end);
                    two.mapping = mapping;
                    for (int k = 0;k<4;k++){
                        vector<pair<unsigned int,string>> temp;
                        one.mappings.push_back(temp);
                        two.mappings.push_back(temp);
                    }
                    for (auto k=0;k<collided_list.size();k++){
                        if (collided_list[k].first<pos){
                            one.mappings[mapping].push_back(collided_list[k]);
                        }
                        else if (collided_list[k].first>pos){
                            two.mappings[mapping].push_back(collided_list[k]);
                        }
                    }
                    dense_area_list.push_back(one);
                    dense_area_list.push_back(two);
                    add_to_final = false;
                    break;
                }
            }

            if(add_to_final){
                final_dense_area_list.push_back(dense_area_list[i]);
                for (auto j=0;j<collided_list.size();j++){
                    string primer = collided_list[j].second;
                    collided_primers[primer]=1; // mark it as abandoned
                }
            }
        }

        cout<<"num of collided primers: "<<collided_primers.size()<<endl;
        cout<<endl;

    }




    
    

    

    


    
    // for(auto it = collided_primers.begin();it != collided_primers.end();it++){
    //     cout<<it->first<<" "<<it->second<<endl;
    // }
    

    

    


    // for (auto i = 0;i<final_dense_area_list.size();i++){
    //     cout<<"dense area: "<<i<<endl;
    //     cout<<"start and end: "<<final_dense_area_list[i].start<<" "<<final_dense_area_list[i].end<<endl;

    //     int mapping = final_dense_area_list[i].mapping;
    //     cout<<"mapping "<<mapping<<endl;
    //     auto collided_list = final_dense_area_list[i].mappings[mapping];
    //     for (auto j = 0;j<collided_list.size();j++){
    //         cout<<"pos and primer: "<<collided_list[j].first<<" "<<collided_list[j].second<<endl;
    //     }
    //     cout<<endl;
    // }

    // different length of dense area
    // vector<int> diff_length;
    // for (auto i =0;i<4;i++){
    //     diff_length.push_back(0);    
    // }
    // for (auto i = 0; i< final_dense_area_list.size();i++){
    //     if(final_dense_area_list[i].end-final_dense_area_list[i].start)==
    // }

    

    // for (auto i = 0;i<dense_area_list.size();i++){

    // }






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

