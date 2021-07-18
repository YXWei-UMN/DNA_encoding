#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <typeinfo>
#include <cstring>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

// important assumption: 
// each collision could be cut at and only at one place
// so each collision could be presented by a cut_pos

// class collision{
//     public:
//         collision(unsigned int s, unsigned int e, string p){
//             start_pos = s;
//             end_pos = e;
//             primer = p;
//         }
//     private:
//         unsigned int start_pos;
//         unsigned int end_pos;
//         string primer;
// };



vector<int> blind_spot;
//map<pair<unsigned int, unsigned int>,string> collision_list;// vector
vector<pair<pair<unsigned int, unsigned int>,string>> collision_list;
vector<unsigned int> cut_positions;
vector<pair<unsigned int, unsigned int>> hybrid_area;


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
    for (int i=10;i<=140;i+=10){
        blind_spot.push_back(i);
    }
    blind_spot.push_back(170);
    blind_spot.push_back(180);
    for (int i=210;i<=290;i+=10){
        blind_spot.push_back(i);
    }
    blind_spot.push_back(330);
    blind_spot.push_back(370);
    for (int i=410;i<=440;i+=10){
        blind_spot.push_back(i);
    }
};

// retrun true if the distance is within the blind spot list
// return false otherwise
bool check_blind_spot(int dis){
    for(auto i=0;i<blind_spot.size();i++){
        if (dis==blind_spot[i])return true;
    }
    return false;
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

int main(int argc, char** argv) {

    if (argc != 2) {
        cerr << "argc must be 2" << endl;
        return -1;
    }

    ifstream myfile;
    myfile.open (argv[1]);   
    string line;

    //int count = 0;

    // getting the collision information
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
        collision_list.push_back(make_pair(make_pair(start,end),primer));


        // count the first 10000 lines of blast result;
        // count++;
        // if (count>10000){break;}
    }
    sort(collision_list.begin(),collision_list.end(),collision_sort);

    // check sort function
    cout<<"start\tend\tprimer\tuncuttable_or_not\tcut_pos"<<endl;
    for (auto it=collision_list.begin();it!=collision_list.end();it++){
        cout<<it->first.first<<" "<<it->first.second<<" "<<it->second<<" ";
        auto cut = find_cut_pos(it->first.first,it->first.second);
        if (cut == 0 ){cout<<" y ";cut = find_trivial_cut_pos(it->first.first,it->first.second);}
        else{cout<<" n ";}
        cout<<cut<<endl;
    }
    cout<<">>>>>>>>>>>>>>>"<<endl;

    init_blind_spot();

    // check blind spot function
    // for (auto it=blind_spot.begin();it!=blind_spot.end();it++){
    //     cout<<*it<<endl;
    // }
    //cout<<check_blind_spot(450)<<endl;
    // cout<<find_trivial_cut_pos(79,100)<<endl;
    // cout<<find_trivial_cut_pos(80,100)<<endl;
    // cout<<find_trivial_cut_pos(81,100)<<endl;


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
        int prev_illegal = illegal;
        illegal = 0;
        if (cut==0){
            cut = find_trivial_cut_pos(start,end);
            illegal = 1;
            // means the next cut area should be hybrid
        }

        if (!check_blind_spot(cut-cur_pos)){
            if (prev_illegal==1){
                hybrid_area_append(cur_pos,cut);
                collision_count_in_hybrid_area++;
            }
            else{
                collision_count_in_sparse_area++;
            }
            cur_pos = cut;

            // if (illegal==2){illegal--;collision_count_in_hybrid_area++;}
            // else if (illegal==1){
            //     hybrid_area_append(cur_pos,cut);
            //     illegal = 0;
            //     collision_count_in_hybrid_area++;
            // }
            // else{collision_count_in_sparse_area++;}
            // cur_pos=cut;
            
        }
        else{ // distance is illegal
                // then check every point between cut and next_cut
                // such that cut to this point and this point to next_cut are both legal
                // if there is one, cut it there and set it as the right boundary of hybrid area
                // if not, check the next collision,
                // until we can find one


            collision_count_in_hybrid_area++;
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
                            flag = true;
                            break;
                        }
                    }
                }
                i++;
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
    }


    unsigned int total_len = collision_list[collision_list.size()-1].first.second;
    unsigned int hybrid_len = 0; 
    cout<<"Existing hybrid area"<<endl;
    //int not_150_count = 0; 
    for (auto it = hybrid_area.begin();it<hybrid_area.end();it++){
        cout<<it->first<<" "<<it->second<<endl;
        //if (it->second-it->first != 150){not_150_count++;}
        hybrid_len += (it->second - it->first);
    }
    cout<<">>>>>>>>>>>"<<endl;
    cout<<"hybrid length: "<<hybrid_len<<endl;
    cout<<"total length: "<<total_len<<endl;
    cout<<"percentage of hybrid area of total len: "<<(double)hybrid_len/total_len<<endl;
    cout<<">>>>>>>>>>>"<<endl;
    cout<<"collision_count_in_hybrid_area: "<<collision_count_in_hybrid_area<<endl;
    cout<<"collision_count_in_sparse_area: "<<collision_count_in_sparse_area<<endl;
    cout<<"total collision count: "<<collision_list.size()<<endl;

    
    //cout<<not_150_count<<" "<<hybrid_area.size()<<endl;
    //most of them are 150 nt area

    return 0;
}

/*




        if(cut==0){
            auto trivial_cut = find_trivial_cut_pos(start,end);
            if(!check_blind_spot(trivial_cut-cur_pos)){
                cur_pos = trivial_cut;
                legal = false;
                continue;
            }
            else{
                auto left = check_back_until_not_blind_spot(cur_pos,trivial_cut);
                if (left!=cur_pos){
                    cur_pos = left;
                }
                auto right = check_until_not_blind_spot(left,trivial_cut);
                hybrid_area_append(left,right);
                continue;
            }
        }
        else{
            if (!check_blind_spot(cut-cur_pos)){// distance is legal
                cur_pos = cut;
                legal = true;
                continue;
            }
            else { // distance is illegal

            }
        }



        if (cut==0){
            auto trivial_cut = find_trivial_cut_pos(start,end);
            if (!check_blind_spot(trivial_cut-cur_pos)){
                cur_pos = trivial_cut;
                // left boundary of hybrid area is trivial cut
                // explore the right boundary of hybrid area
                if (i==collision_list.size()-1){
                    auto right_boundary = check_until_not_blind_spot(trivial_cut,trivial_cut);
                    hybrid_area_append(trivial_cut,right_boundary);
                    continue;
                }
                else{
                    auto next_start = collision_list[i+1].first.first;
                    auto next_end = collision_list[i+1].first.second;
                    auto next_cut = find_cut_pos(next_start,next_end);
                    if(next_cut==0){
                        auto next_trivial_cut = find_trivial_cut_pos(next_start,next_end);
                        if(!check_blind_spot(next_trivial_cut-trivial_cut)){
                            hybrid_area_append(trivial_cut,next_trivial_cut);
                            cur_pos = next_trivial_cut;
                            continue;
                        }
                        else{
                            auto next_right = check_until_not_blind_spot(trivial_cut,next_trivial_cut);

                        }
                        
                    }
                }

            }
            else{
                auto next_cut = check_until_not_blind_spot(cur_pos,trivial_cut);
            }

        }







        if (cut==0){
            hybrid_area_append(i,i);
            continue;
        }
        else{
            int dis = cut - cur_pos;
            if (!check_blind_spot(dis)){// distance is legal
                if (i==collision_list.size()-1){
                    cur_pos = cut;
                    continue;
                }
                else {
                    auto next_start = collision_list[i+1].first.first;
                    auto next_end = collision_list[i+1].first.second;
                    auto next_cut = find_cut_pos(next_start,next_end);
                    if (next_cut==0 || !check_blind_spot(next_cut-cut)){ // next collision is intrinsic or next cut is legal 
                        cur_pos = cut;
                        continue;
                    }
                    else{
                        if (i==collision_list.size()-2){
                            cur_pos = cut;
                            continue;
                        }
                        else {
                            auto next_next_start = collision_list[i+2].first.first;
                            auto next_next_end = collision_list[i+2].first.second;
                            auto next_next_cut = find_cut_pos(next_next_start,next_next_end);
                            if (next_next_cut==0 || !check_blind_spot(next_cut-cut)){ // next collision is intrinsic or next cut is legal 
                                cur_pos = cut;
                                continue;
                            }
                        }
                    }
                }
                
            }
            else{ // distance is illegal
                cur_pos = cut;
                hybrid_area_append(i,i);
                continue;
            }
        }
    }


    // auto it=collision_list.begin();
    // while(it!=collision_list.end()){
    //     auto start = it->first.first;
    //     auto end = it->first.second;
    //     auto cut = find_cut_pos(start,end);
    //     cout<<cut<<endl;

    //     if (cut!=0){
    //         int dis = cut - cur_pos;
    //         if (!check_blind_spot(dis)){ // distance is legal
    //             cut_positions.push_back(cut);
    //             cur_pos = cut;
    //             it++;
    //             continue;
    //         }
    //     }
    //     else{ // cut == 0
    //         hybrid_area_append(start,end);
    //     }
    //     it++;


    // }





    return 0;
}
*/
