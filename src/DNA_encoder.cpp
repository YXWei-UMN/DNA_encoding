//
// Created by eason on 2/22/21.
//
#include <math.h>
#include "DNA_encoder.h"
#include "tool.h"
#include "transformation.h"
#include "rs.hpp"
#include <cstdlib>
#include <cstring>
#include <cassert>



string DNA_encoder::Church_encoding(string digital_data) {
    string result;
    string last_three="ATC";
    string at[2]={"A","T"};
    string cg[2]={"C","G"};
    for (std::size_t i = 0; i < digital_data.size(); i++)
    {
        //cover from binary to decimal
        bitset<8> bits(digital_data.c_str()[i]);
        for (int j = 0; j < 8; j++) {
            if (bits[j]==0){
                if (last_three[0]=='A' && last_three[1]=='A' && last_three[2]=='A'){
                    result+="T";
                    last_three.erase(0,1);
                    last_three+="T";
                } else if(last_three[0]=='T' && last_three[1]=='T' && last_three[2]=='T'){
                    result+="A";
                    last_three.erase(0,1);
                    last_three+="A";
                } else{
                    string str = at[rand()%2];
                    last_three.erase(0,1);
                    last_three+=str;
                    result+=str;
                }
            }
            else if (bits[j]==1){
                if (last_three[0]=='C' && last_three[1]=='C' && last_three[2]=='C'){
                    result+="G";
                    last_three.erase(0,1);
                    last_three+="G";
                } else if(last_three[0]=='G' && last_three[1]=='G' && last_three[2]=='G'){
                    result+="C";
                    last_three.erase(0,1);
                    last_three+="C";
                } else{
                    string str = cg[rand()%2];
                    last_three.erase(0,1);
                    last_three+=str;
                    result+=str;
                }

            }
        }
    }
    return result;
}

//Forward Error Correction encoding
string DNA_encoder::FEC_encoding(string digital_data) {
    string result;
    int count = 0;
    // count==0: cluster A
    // count==1: cluster B
    for (std::size_t i = 0; i < digital_data.size(); i++)
    {
        string nts;
        bitset<8> bits(digital_data.c_str()[i]);
        int first = bits[0]*2+bits[1];
        int second = bits[2]*2+bits[3];
        int third = bits[4]*2+bits[5];
        int fourth = bits[6]*2+bits[7];
        // get first, second and fourth nts using table1
        // get third and fifth nts using table2

        nts+=FEC_table1[first];
        nts+=FEC_table1[second];

        // flag==true: correct option
        bool flag = (count==0)?true:false;
        count += (count==0)?1:-1;
        for (auto j=0;j<FEC_table2[fourth].size();j++){
            if (FEC_table1[first]==FEC_table1[second] && FEC_table1[first]==FEC_table2[fourth][j][0]){
                continue;
            }
            if (FEC_table2[fourth][j][1]==FEC_table1[third]){
                continue;
            }
            if (flag){
                nts+=FEC_table2[fourth][j][0];
                nts+=FEC_table1[third];
                nts+=FEC_table2[fourth][j][1];
                break;
            }
            else {
                flag=true;
                continue;
            }
        }
        result+=nts;
        
    }
    return result;
}
string DNA_encoder::RS_2bits_rotation(string digital_data){
    // step 0: padding to a multiple of 3
    //cout << "digital_data.size() = " << digital_data.size() << endl;
    while (digital_data.size() % 30 != 0){
        digital_data = digital_data + '\0';
    }
    //cout << "digital_data.size() = " << digital_data.size() << endl;

    // step 1: apply RS encoding on each 30-bytes block
    const int MSG_LENGTH = 30;
    const int ECC_LENGTH = 6;
    const int ENCODED_LENGTH = MSG_LENGTH + ECC_LENGTH;
    RS::ReedSolomon<MSG_LENGTH, ECC_LENGTH> rs;

    int n_rs_unit = digital_data.size() / MSG_LENGTH;
    //cout << "n_rs_unit = " << n_rs_unit << endl;
    //cout << "encoded_unit_length = " << ENCODED_LENGTH * n_rs_unit << endl;
    char rs_result[ENCODED_LENGTH * n_rs_unit];
    for (int i = 0, j = 0; i < digital_data.size(); i += MSG_LENGTH, j += ENCODED_LENGTH) {
        char message[MSG_LENGTH+1];
        message[MSG_LENGTH] = '\0';
        memcpy(message, digital_data.c_str()+i, MSG_LENGTH);
        // cout << "message = " << message << endl;
        char cur_encoded[ENCODED_LENGTH+1];
        cur_encoded[ENCODED_LENGTH] = '\0';
        rs.Encode(message, cur_encoded);
        // cout << "encoded message = " << cur_encoded << endl;
        // cout << "encoded length = " << strlen(cur_encoded) << endl;
        memcpy(rs_result + j, cur_encoded, ENCODED_LENGTH);
    }

//rotation code
    string result;
    for (std::size_t i = 0; i < sizeof(rs_result); i++)
    {
        bitset<8> bits(rs_result[i]);
        for (int j = 0; j < 8; j+=2) {
            if (bits[j]==0){
                if (bits[j+1]==0) {
                    last_bit=rotating_encoding_table_[0][last_bit];
                    result+=last_bit;
                }
                else if (bits[j+1]==1) {
                    last_bit=rotating_encoding_table_[1][last_bit];
                    result+=last_bit;
                }
            }
            else if (bits[j]==1){
                if (bits[j+1]==0) {
                    last_bit=rotating_encoding_table_[2][last_bit];
                    result+=last_bit;
                }
                else if (bits[j+1]==1) {
                    last_bit=rotating_encoding_table_[3][last_bit];
                    result+=last_bit;
                }
            }
        }
    }
    return result;
}
string DNA_encoder::RS_rotation(string digital_data) {
    // step 0: padding to a multiple of 3
    //cout << "digital_data.size() = " << digital_data.size() << endl;
    while (digital_data.size() % 30 != 0){
        digital_data = digital_data + '\0';
    }
    //cout << "digital_data.size() = " << digital_data.size() << endl;

    // step 1: apply RS encoding on each 30-bytes block
    const int MSG_LENGTH = 30;
    const int ECC_LENGTH = 6;
    const int ENCODED_LENGTH = MSG_LENGTH + ECC_LENGTH;
    RS::ReedSolomon<MSG_LENGTH, ECC_LENGTH> rs;

    int n_rs_unit = digital_data.size() / MSG_LENGTH;
    //cout << "n_rs_unit = " << n_rs_unit << endl;
    //cout << "encoded_unit_length = " << ENCODED_LENGTH * n_rs_unit << endl;
    char rs_result[ENCODED_LENGTH * n_rs_unit];
    for (int i = 0, j = 0; i < digital_data.size(); i += MSG_LENGTH, j += ENCODED_LENGTH) {
        char message[MSG_LENGTH+1];
        message[MSG_LENGTH] = '\0';
        memcpy(message, digital_data.c_str()+i, MSG_LENGTH);
        // cout << "message = " << message << endl;
        char cur_encoded[ENCODED_LENGTH+1];
        cur_encoded[ENCODED_LENGTH] = '\0';
        rs.Encode(message, cur_encoded);
        // cout << "encoded message = " << cur_encoded << endl;
        // cout << "encoded length = " << strlen(cur_encoded) << endl;
        memcpy(rs_result + j, cur_encoded, ENCODED_LENGTH);
    }

//rotation code
    string result;
    for (std::size_t i = 0; i < sizeof(rs_result); i+=2)
    {
        //cover from binary to decimal
        bitset<8> bits(rs_result[i]);
        int decimal_1 = bits.to_ulong();
        bitset<8> bits_2(rs_result[i+1]);
        int decimal_2 = bits_2.to_ulong();
        int decimal = decimal_2*1024+decimal_1;

        //covert to ternary & rotate encoding
        for (int i=0; i<11; i++){
            last_bit=rotating_encoding_table_[decimal%3][last_bit];
            result+=last_bit;
            decimal/=3;
        }
    }
    return result;
}

string DNA_encoder::RS_GF47(string digital_data) {
    // step 0: padding to a multiple of 3
    //cout << "digital_data.size() = " << digital_data.size() << endl;
    while (digital_data.size() % 30 != 0){
        digital_data = digital_data + '\0';
    }
    //cout << "digital_data.size() = " << digital_data.size() << endl;

    // step 1: apply RS encoding on each 30-bytes block
    const int MSG_LENGTH = 30;
    const int ECC_LENGTH = 6;
    const int ENCODED_LENGTH = MSG_LENGTH + ECC_LENGTH;
    RS::ReedSolomon<MSG_LENGTH, ECC_LENGTH> rs;

    int n_rs_unit = digital_data.size() / MSG_LENGTH;
    //cout << "n_rs_unit = " << n_rs_unit << endl;
    //cout << "encoded_unit_length = " << ENCODED_LENGTH * n_rs_unit << endl;
    char rs_result[ENCODED_LENGTH * n_rs_unit];
    for (int i = 0, j = 0; i < digital_data.size(); i += MSG_LENGTH, j += ENCODED_LENGTH) {
        char message[MSG_LENGTH+1];
        message[MSG_LENGTH] = '\0';
        memcpy(message, digital_data.c_str()+i, MSG_LENGTH);
        // cout << "message = " << message << endl;
        char cur_encoded[ENCODED_LENGTH+1];
        cur_encoded[ENCODED_LENGTH] = '\0';
        rs.Encode(message, cur_encoded);
        // cout << "encoded message = " << cur_encoded << endl;
        // cout << "encoded length = " << strlen(cur_encoded) << endl;
        memcpy(rs_result + j, cur_encoded, ENCODED_LENGTH);
    }


    // step 2: apply mapping from encoded bytes to nt
    int total_rs_len = ENCODED_LENGTH * n_rs_unit;
    int result_len = total_rs_len / 2 * 3 * 3; // 2 -> two 256-base number to use; 3 -> three 47-base number generated; 3 -> three nt for each 47-base number
    //cout << "result_len = " << result_len << endl;
    char result_cstring[result_len+1];
    result_cstring[result_len] = '\0';

    int j = 0;
    for (int i = 0; i < total_rs_len; i += 2) {
        // convert from 256^2 to 47^3
        unsigned int middle_num = (((unsigned int)rs_result[i])&0xff) * 256 + (((unsigned int)rs_result[i+1])&0xff);
        // cout << "l r m = " << hex << (unsigned int)(char)(rs_result[i]) << '\t' << (unsigned int)(rs_result[i+1]) << '\t' << middle_num << endl << flush;
        unsigned int base47_num[3];
        base47_num[2] = middle_num % 47;
        middle_num /= 47;
        base47_num[1] = middle_num % 47;
        middle_num /= 47;
        base47_num[0] = middle_num % 47;
        middle_num /= 47;
        assert(middle_num == 0);

        // mapping to nt
        for (int k1 = 0; k1 < 3; k1++) {
            for (int k2 = 0; k2 < 3; k2++) {
                // cout << j << ' ' << base47_num[k1] << ' ' << endl;;
                result_cstring[j++] = RS_table[base47_num[k1]][k2];
                //cout << j-1 << ' ' << RS_table[base47_num[k1]][k2] << endl;
            }
        }

    }
    assert(j == result_len);

    string result(result_cstring);
    return result;
}

string DNA_encoder::RS_GF25(string digital_data) {
    // step 0: padding to a multiple of 3
    //cout << "digital_data.size() = " << digital_data.size() << endl;
    while (digital_data.size() % 30 != 0){
        digital_data = digital_data + '\0';
    }
    //cout << "digital_data.size() = " << digital_data.size() << endl;

    // step 1: apply RS encoding on each 30-bytes block
    const int MSG_LENGTH = 30;
    const int ECC_LENGTH = 6;
    const int ENCODED_LENGTH = MSG_LENGTH + ECC_LENGTH;
    RS::ReedSolomon<MSG_LENGTH, ECC_LENGTH> rs;

    int n_rs_unit = digital_data.size() / MSG_LENGTH;
    //cout << "n_rs_unit = " << n_rs_unit << endl;
    //cout << "encoded_unit_length = " << ENCODED_LENGTH * n_rs_unit << endl;
    char rs_result[ENCODED_LENGTH * n_rs_unit];
    for (int i = 0, j = 0; i < digital_data.size(); i += MSG_LENGTH, j += ENCODED_LENGTH) {
        char message[MSG_LENGTH+1];
        message[MSG_LENGTH] = '\0';
        memcpy(message, digital_data.c_str()+i, MSG_LENGTH);
        // cout << "message = " << message << endl;
        char cur_encoded[ENCODED_LENGTH+1];
        cur_encoded[ENCODED_LENGTH] = '\0';
        rs.Encode(message, cur_encoded);
        // cout << "encoded message = " << cur_encoded << endl;
        // cout << "encoded length = " << strlen(cur_encoded) << endl;
        memcpy(rs_result + j, cur_encoded, ENCODED_LENGTH);
    }


    // step 2: apply mapping from encoded bytes to nt
    // convert from 2^8 to 24^2
    int total_rs_len = ENCODED_LENGTH * n_rs_unit;
    vector<char> result_cstring;

    for (int i = 0; i < total_rs_len; i++) {

        unsigned int middle_num = (((unsigned int)rs_result[i+1])&0xff);
        // cout << "l r m = " << hex << (unsigned int)(char)(rs_result[i]) << '\t' << (unsigned int)(rs_result[i+1]) << '\t' << middle_num << endl << flush;
        unsigned int base24_num[2];
        base24_num[1] = middle_num % 25;
        middle_num /= 25;
        base24_num[0] = middle_num % 25;
        middle_num /= 25;
        assert(middle_num == 0);

        // mapping to nt
        for (int k1 = 0; k1 < 2; k1++) {
            for (int k2 = 0; k2 < 3; k2++) {
                // cout << j << ' ' << base47_num[k1] << ' ' << endl;;
                result_cstring.push_back(RS_table[base24_num[k1]][k2]);
                //cout << j-1 << ' ' << RS_table[base47_num[k1]][k2] << endl;
            }
        }
    }

    string result(result_cstring.begin(),result_cstring.end());
    return result;
}

//direct mapping 00-A 01-T 10-C 11-G
string DNA_encoder::direct_encoding(string digital_data) {
    string result;
    for (std::size_t i = 0; i < digital_data.size(); i++)
    {
        //cover from binary to decimal
        bitset<8> bits(digital_data.c_str()[i]);
        for (int j = 0; j < 8; j+=2) {
            if (bits[j]==0){
                if (bits[j+1]==0){
                    result+="A";
                } else if (bits[j+1]==1){
                    result+="T";
                } else {
                    cerr<<"bit not 1 or 0"<<endl;
                }
            }
            else if (bits[j]==1){
                if (bits[j+1]==0){
                    result+="C";
                } else if (bits[j+1]==1){
                    result+="G";
                } else {
                    cerr<<"bit not 1 or 0"<<endl;
                }
            }
        }
    }
    return result;
}
string DNA_encoder::complementary_sequence(string str) {
    string complementary_str;
    for(auto n:str){
        switch (n){
            case 'A':{
                complementary_str.insert (0, 1, 'T');
                break;
            }
            case 'T':{
                complementary_str.insert (0, 1, 'A');
                break;
            }
            case 'C':{
                complementary_str.insert (0, 1, 'G');
                break;
            }
            case 'G':{
                complementary_str.insert (0, 1, 'C');
                break;
            }
        }
    }
    return complementary_str;
}

string DNA_encoder::exhaustive_search_triplets() {
    string result;
    vector<pair<int,string>> all_triplet_candidates_and_score;


    list<string> homo1;
    //homo
    for (auto n:NT_triplets_){
        if (last_17nt_.back()!=n[0])
            homo1.push_back(n);   // since n has no homo, as long as
    }

    for (auto m:homo1){
        int matches=0;
        string complement_str1 = complementary_sequence(m);

        //str 1
        auto it = last_17nt_.begin();
        for (int i = 0; i < 15; ++i) {
            auto tmp_it = it;
            if (*tmp_it != complement_str1[0]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complement_str1[1]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complement_str1[2]) {
                it++;
                continue;
            }
            matches++;
            break;
        }


        //str 2
        string str2;
        str2.insert (0, 1, m[0]);
        it = last_17nt_.end();
        it--;
        str2.insert (0, 1, *it);
        it--;
        str2.insert (0, 1, *it);
        string complement_str2 = complementary_sequence(str2);

        it = last_17nt_.begin();
        for (int i = 0; i < 15; ++i) {
            auto tmp_it = it;
            if (*tmp_it != complement_str2[0]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complement_str2[1]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complement_str2[2]) {
                it++;
                continue;
            }
            matches++;
            break;
        }


        //str 3
        string str3;
        str3.insert (0, 1, m[1]);
        str3.insert (0, 1, m[0]);
        it = last_17nt_.end();
        it--;
        str3.insert (0, 1, *it);
        string complement_str3 = complementary_sequence(str3);

        it = last_17nt_.begin();
        for (int i = 0; i < 15; ++i) {
            auto tmp_it = it;
            if (*tmp_it != complement_str3[0]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complement_str3[1]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complement_str3[2]) {
                it++;
                continue;
            }
            matches++;
            break;
        }

        //cout<<"triplet "<<m<<" matches "<<matches<<endl;
        all_triplet_candidates_and_score.push_back(make_pair(matches,m));
    }

    sort(all_triplet_candidates_and_score.begin(), all_triplet_candidates_and_score.end());

    int rad = rand()%8;
    result=all_triplet_candidates_and_score[rad].second;
    //cout<<"select top "<<rad<<" "<<result<<"   its matches are"<<all_triplet_candidates_and_score[rad].first<<result<<endl;

    //replace last 17
    for (int i = 0; i < result.size(); ++i) {
        last_17nt_.emplace_back(result[i]);
        last_17nt_.erase(last_17nt_.begin());
    }
    return result;
}

string DNA_encoder::partial_search_triplets() {
    string result;
    vector<pair<int,string>> all_triplet_candidates_and_score;
    uint8_t binary=rand()%8;
    list<string> candidates = three_bits_NT_triplets_candidates_[binary];

    list<string> homo1;
    //homo
    for (auto n:candidates){
        if (last_17nt_.back()!=n[0])
            homo1.push_back(n);   // since n has no homo, as long as
    }

    unordered_map<string,list<string>> concatenate_triplets_of_each_candidates;
    for (auto m:homo1){
        int matches=0;
        list<string> concatenate_triplets;

        concatenate_triplets.emplace_back(m);
        //string complement_str1 = complementary_sequence(m);
        string complement_str1 = complementary_triplets_table_.find(m)->second;
        //str 1
        if (triplets_in_last_17nt_.find(complement_str1)!=triplets_in_last_17nt_.end())
            matches++;


        //str 2
        string str2;
        str2.insert (0, 1, m[0]);
        auto it = last_17nt_.end();
        it--;
        str2.insert (0, 1, *it);
        it--;
        str2.insert (0, 1, *it);
        concatenate_triplets.emplace_back(str2);

        string complement_str2 = complementary_triplets_table_.find(str2)->second;
        if (triplets_in_last_17nt_.find(complement_str2)!=triplets_in_last_17nt_.end())
            matches++;




        //str 3
        string str3;
        str3.insert (0, 1, m[1]);
        str3.insert (0, 1, m[0]);
        it = last_17nt_.end();
        it--;
        str3.insert (0, 1, *it);
        concatenate_triplets.emplace_back(str3);

        string complement_str3 = complementary_triplets_table_.find(str3)->second;
        if (triplets_in_last_17nt_.find(complement_str3)!=triplets_in_last_17nt_.end())
            matches++;


        //cout<<"triplet "<<m<<" matches "<<matches<<endl;
        all_triplet_candidates_and_score.emplace_back(make_pair(matches,m));
        concatenate_triplets_of_each_candidates.emplace(m,concatenate_triplets);
    }

    sort(all_triplet_candidates_and_score.begin(), all_triplet_candidates_and_score.end());

    /*cout<<"   -----    "<<endl;
    for (auto a:all_triplet_candidates_and_score){
        cout<<a.second<<" matches"<<a.first<<endl;
    }*/


    result=all_triplet_candidates_and_score.begin()->second;
    //cout<<"select top "<<rad<<" "<<result<<"   its matches are"<<all_triplet_candidates_and_score[rad].first<<result<<endl;

    //replace last 17
    for (int i = 0; i < result.size(); ++i) {
        last_17nt_.emplace_back(result[i]);
        last_17nt_.erase(last_17nt_.begin());
    }

    //insert new triplets to triplets_in_last_17nt_
    list<string> concatenate_triplets_of_select_triplets = concatenate_triplets_of_each_candidates.find(result)->second;
    for (auto s:concatenate_triplets_of_select_triplets) {
        if (triplets_in_last_17nt_.find(s)==triplets_in_last_17nt_.end())
            triplets_in_last_17nt_.emplace(s,1);
        else
            triplets_in_last_17nt_.find(s)->second++;
    }

    //delete old triplets from triplets_in_last_17nt_
    auto it = last_17nt_.begin();
    string delete_str1, delete_str2, delete_str3;
    delete_str1.push_back(*it);

    it++;
    delete_str1.push_back(*it);
    delete_str2.push_back(*it);

    it++;
    delete_str1.push_back(*it);
    delete_str2.push_back(*it);
    delete_str3.push_back(*it);

    it++;
    delete_str2.push_back(*it);
    delete_str3.push_back(*it);

    it++;
    delete_str3.push_back(*it);

    auto ptr1 = triplets_in_last_17nt_.find(delete_str1);
    if (ptr1==triplets_in_last_17nt_.end()){
        cout<<"error in delete old triplets:"<<delete_str1<<" from triplets_in_last_17nt_";
        cout<<"if not happen at begin, must be sth wrong"<<endl;
    } else{
        if (ptr1->second==1)
            triplets_in_last_17nt_.erase(ptr1);
        else
            ptr1->second--;
    }




    auto ptr2 = triplets_in_last_17nt_.find(delete_str2);
    if (ptr2==triplets_in_last_17nt_.end()){
        cout<<"error in delete old triplets:"<<delete_str2<<" from triplets_in_last_17nt_";
        cout<<"if not happen at begin, must be sth wrong"<<endl;
    } else{
        if (ptr2->second==1)
            triplets_in_last_17nt_.erase(ptr2);
        else
            ptr2->second--;
    }




    auto ptr3 = triplets_in_last_17nt_.find(delete_str3);
    if (ptr3==triplets_in_last_17nt_.end()){
        cout<<"error in delete old triplets:"<<delete_str3<<" from triplets_in_last_17nt_";
        cout<<"if not happen at begin, must be sth wrong"<<endl;
    } else{
        if (ptr3->second==1)
            triplets_in_last_17nt_.erase(ptr3);
        else
            ptr3->second--;
    }

    return result;
}


string DNA_encoder::heuristic_encoding(string digital_data) {
    string result;

        //grab 30*8 bits a time
        for (int j = 0; j < 30; ++j) {
            bitset<8> bits1(digital_data.c_str()[j]);
        }


        for (int j = 0; j < 240; j+=g_num_bit_per_triplet) {
            //mapping between binary bits and NT triplets
            //because in experiment, we dont have to decode, we can just encode with the top one

            //list<string>triplets = sort_triplets();
            if (g_encoding_scheme==5)
                result+=exhaustive_search_triplets();
            if (g_encoding_scheme==6)
                result+=partial_search_triplets();
        }

    return result;
}

string DNA_encoder::prefix_encoding(string digital_data) {
    string result;
    for (std::size_t i = 0; i < digital_data.size(); i++)
    {
        //cover from binary to decimal
        bitset<8> bits(digital_data.c_str()[i]);
        for (int j = 0; j < 8; j+=2) {
            if (bits[j]==1){
               result+=*last_20nt_.begin();
               last_20nt_.emplace_back(*last_20nt_.begin());
               last_20nt_.erase(last_20nt_.begin());
            }
            else if (bits[j]==0){
                if (*last_20nt_.begin()=="A"){
                    result+="T";
                    last_20nt_.emplace_back("T");
                    last_20nt_.erase(last_20nt_.begin());
                } else if (*last_20nt_.begin()=="T"){
                    result+="A";
                    last_20nt_.emplace_back("A");
                    last_20nt_.erase(last_20nt_.begin());
                }else if (*last_20nt_.begin()=="C"){
                    result+="G";
                    last_20nt_.emplace_back("G");
                    last_20nt_.erase(last_20nt_.begin());
                } else if (*last_20nt_.begin()=="G"){
                    result+="C";
                    last_20nt_.emplace_back("C");
                    last_20nt_.erase(last_20nt_.begin());
                } else {
                    cerr<<"bit not 1 or 0"<<endl;
                    EXIT_FAILURE;
                }
            }
        }
    }
    return result;
}


// strand level randomize: pseudo_strand XOR encoded_strand
void DNA_encoder::randomize_XOR(string &digital_data) {
    for(int i=0; i<digital_data.size(); i+=pseudo_random_sequence_.size()){
        string original=digital_data.substr(i,pseudo_random_sequence_.size());
        string new_str;
        for (int j = 0; j < pseudo_random_sequence_.size(); ++j) {
            new_str+=original[j]^pseudo_random_sequence_[j];
        }
        digital_data.replace(i,pseudo_random_sequence_.size(),new_str);
    }
}


// fixed 200 strand
void DNA_encoder::encoding_stranding(){
    // create payload file to store encoded strands
    fstream payload_file;
    int prefix = rand();
    string payload_path = g_payload_path+to_string(prefix)+"payload";
    int num_of_GB = 1;
    unsigned long total_len = 0;
    payload_path += to_string(num_of_GB);
    payload_file.open(payload_path,ios::out);

    long int strand_num=0;
    //create reading buffer
    uint8_t buf[30];
    //go over all files to chunking and encoding
    FILE *fp;
    string nt_sequence;




/*    // Todo remember to decommand all file initialization
    fp = fopen(g_data_path.c_str(), "r");
    if (fp==NULL) {fputs ("File open error",stderr); exit (1);}
    while ( !feof(fp) ) {

        //cout<<total_len<<endl;
        size_t len = fread(buf, 1, sizeof(buf), fp);
        // total_len += len;
        uint8_t *ptr = &buf[0];

        string digital_data ((char*)ptr,len);

        string strand=heuristic_encoding(digital_data);
        payload_file<<">payload"<<strand_num++<<endl;
        // execute transformation: mapping/swap/...
        payload_file<<strand<<endl;
        total_len += strand.length();


    }
    fclose(fp);*/


    for(auto n:all_files_){
        fp = fopen(n.c_str(), "r");
        if (fp==NULL) {fputs ("File open error",stderr); exit (1);}
        while ( !feof(fp) ) {

            // TODO remember to change file size when change num_of_bits_per_triplet
            if (total_len>=(long)130*1024*1024){
                total_len = 0;
                payload_file.close();
                num_of_GB += 1;
                payload_path = g_payload_path+to_string(prefix)+"payload"+to_string(num_of_GB);
                payload_file.open(payload_path,ios::out);
            }
            //cout<<total_len<<endl;
            size_t len = fread(buf, 1, sizeof(buf), fp);
            // total_len += len;
            uint8_t *ptr = &buf[0];

            string digital_data ((char*)ptr,len);


            // rotate code
            if(g_encoding_scheme==1) {
                if (g_if_randomization) {
                    randomize_XOR(digital_data);
                }
                string strand = RS_rotation(digital_data);
                //string strand=RS_2bits_rotation(digital_data);
                //string strand=Church_encoding(digital_data);
                payload_file << ">payload" << strand_num++ << endl;
                // execute transformation: mapping/swap/...
                payload_file << strand << endl;
                total_len += sizeof(buf);
            }
            else if(g_encoding_scheme==0){
                string strand=direct_encoding(digital_data);
                payload_file<<">payload"<<strand_num++<<endl;
                // execute transformation: mapping/swap/...
                payload_file<<strand<<endl;
                total_len += sizeof(buf);
            }
            else if(g_encoding_scheme==2){
                // 200 nts correspoinding to 320 bits
                // 2 nts have to be reserved for ECC (CCITT16: 16 bits)
                // each oligo store information with size of 320-16 = 304 bits -> 38 char
                while (digital_data.size()>0){
                    string digital_strand;
                    
                    if (g_if_ECC){
                        digital_strand = digital_data.substr(0,38);

                        // add ECC code
                        char strand_char[digital_strand.length()+1];
                        strcpy(strand_char, digital_strand.c_str());
                        uint16_t ECC = CCITT16(strand_char,strlen(strand_char));

                        // convert 16-bit int to char array
                        char ECC_char[2];
                        ECC_char[0] = (ECC>>8);
                        ECC_char[1] = (ECC&0xff);
                        digital_strand += ECC_char[0];
                        digital_strand += ECC_char[1];
                    }
                    else {
                        digital_strand = digital_data.substr(0,40);
                    }
                    
                    string strand = FEC_encoding(digital_strand);
                    payload_file<<">payload"<<strand_num++<<endl;
                    // execute transformation: mapping/swap/...
                    payload_file<<strand<<endl;
                    total_len += sizeof(buf);
                    if (g_if_ECC){
                        digital_data.erase(0, 38);
                    }
                    else {
                        digital_data.erase(0, 40);
                    }
                }
                //nt_sequence=FEC_encoding(digital_data);
            }
            else if(g_encoding_scheme==4){ //prefix encoding
                string strand=prefix_encoding(digital_data);
                payload_file<<">payload"<<strand_num++<<endl;
                // execute transformation: mapping/swap/...
                payload_file<<strand<<endl;
                total_len += sizeof(buf);
            }
            else if(g_encoding_scheme==5 || g_encoding_scheme==6){ //heuristic encoding
                string strand=heuristic_encoding(digital_data);
                payload_file<<">payload"<<strand_num++<<endl;
                // execute transformation: mapping/swap/...
                payload_file<<strand<<endl;
                total_len += sizeof(buf);
            }
            else
                cout<<"no encoding scheme"<<endl;
        }
        fclose(fp);
    }
    payload_file.close();
}

// no strand
// it's the original video_data stream feed for our algorithm (to cut & transformation)
void DNA_encoder::encoding_no_strand(){

    // create payload file to store encoded strands
    fstream payload_file;
    string payload_path = g_payload_path+"payload";
    int num_of_payload_file = 1;
    unsigned long total_len = 0;
    payload_path += to_string(num_of_payload_file);
    payload_path += ".txt";
    payload_file.open(payload_path,ios::out);
    int num_of_file = 0;
    payload_file<<">payload"<<num_of_file<<endl;
    num_of_file++;
    //create chunking buffer and related structure
    uint8_t buf[1024*1024];
    //go over all files to chunking and encoding
    FILE *fp;
    string nt_sequence;


    for(auto n:all_files_){
        fp = fopen(n.c_str(), "r");
        if (fp==NULL) {fputs ("File open error",stderr); exit (1);}

        while ( !feof(fp) ) {
            size_t len = fread(buf, 1, sizeof(buf), fp);
            // total_len += len;
            uint8_t *ptr = &buf[0];

            string digital_data ((char*)ptr,len);

            if(g_encoding_scheme==1)
                nt_sequence=RS_rotation(digital_data);
                //nt_sequence=base3_rotate_encoding(digital_data);
            else if (g_encoding_scheme==2)
                nt_sequence=FEC_encoding(digital_data);
            else if (g_encoding_scheme==3)
                //nt_sequence=RS_GF47(digital_data);
                nt_sequence=RS_GF25(digital_data);
            else
                cout<<"no encoding scheme"<<endl;
            payload_file<<nt_sequence;
            total_len += sizeof(buf);

            if (total_len>=(long)5*1024*1024*1024){
                total_len = 0;
                payload_file.close();
                num_of_payload_file += 1;
                payload_path = g_payload_path+"payload"+to_string(num_of_payload_file)+".txt";
                payload_file.open(payload_path,ios::out);
                payload_file<<">payload"<<num_of_file<<endl;
                num_of_file++;
            }

        }
	payload_file<<endl;
        fclose(fp);
    }
    payload_file.close();
}

uint16_t DNA_encoder::CCITT16(char *ptr, int length)
{
   uint16_t crc = 0xffff;
   int i;
   while (length > 0)
   {
      crc = crc ^ (unsigned char) *ptr++ << 8;
	  for (i=0;i<8;i++){
        if (crc & 0x8000)
        	crc = crc << 1 ^ 0x1021;
        else
            crc = crc << 1;
	  }
	  length--;
   }
   return (crc & 0xffff);
}

void DNA_encoder::initial_FEC_table(){
    FEC_table1 = "ACGT";
    FEC_table2 = {{"AA","CC","GG","TT"},
                {"AC","CG","GT","TA"},
                {"AG","CT","GA","TC"},
                {"AT","CA","GC","TG"}};
}

void DNA_encoder::init_heuristic_encoding() {
    // init last 17 nt: TCACTATGCTCGCTATG
    // provide all NT triplets
    last_17nt_.emplace_back('T');
    last_17nt_.emplace_back('C');
    last_17nt_.emplace_back('A');
    last_17nt_.emplace_back('C');
    last_17nt_.emplace_back('T');
    last_17nt_.emplace_back('A');
    last_17nt_.emplace_back('T');
    last_17nt_.emplace_back('G');
    last_17nt_.emplace_back('C');
    last_17nt_.emplace_back('T');
    last_17nt_.emplace_back('C');
    last_17nt_.emplace_back('G');
    last_17nt_.emplace_back('C');
    last_17nt_.emplace_back('T');
    last_17nt_.emplace_back('A');
    last_17nt_.emplace_back('T');
    last_17nt_.emplace_back('G');

    //NT_triplets_.emplace_back("AAA");
    //NT_triplets_.emplace_back("AAT");
    //NT_triplets_.emplace_back("AAC");
    //NT_triplets_.emplace_back("AAG");
            NT_triplets_.emplace_back("ATA");
    //NT_triplets_.emplace_back("ATT");
            NT_triplets_.emplace_back("ATC");
            NT_triplets_.emplace_back("ATG");
            NT_triplets_.emplace_back("ACA");
                    NT_triplets_.emplace_back("ACT");
    //NT_triplets_.emplace_back("ACC");
            NT_triplets_.emplace_back("ACG");
            NT_triplets_.emplace_back("AGA");
                    NT_triplets_.emplace_back("AGT");
            NT_triplets_.emplace_back("AGC");
    //NT_triplets_.emplace_back("AGG");


    //NT_triplets_.emplace_back("TAA");
            NT_triplets_.emplace_back("TAT");
            NT_triplets_.emplace_back("TAC");
            NT_triplets_.emplace_back("TAG");
    //NT_triplets_.emplace_back("TTA");
    //NT_triplets_.emplace_back("TTT");
    //NT_triplets_.emplace_back("TTC");
    //NT_triplets_.emplace_back("TTG");
                    NT_triplets_.emplace_back("TCA");
            NT_triplets_.emplace_back("TCT");
    //NT_triplets_.emplace_back("TCC");
            NT_triplets_.emplace_back("TCG");
                    NT_triplets_.emplace_back("TGA");
            NT_triplets_.emplace_back("TGT");
            NT_triplets_.emplace_back("TGC");
    //NT_triplets_.emplace_back("TGG");

    //NT_triplets_.emplace_back("CAA");
            NT_triplets_.emplace_back("CAT");
            NT_triplets_.emplace_back("CAC");
            NT_triplets_.emplace_back("CAG");
            NT_triplets_.emplace_back("CTA");
    //NT_triplets_.emplace_back("CTT");
    NT_triplets_.emplace_back("CTC");
            NT_triplets_.emplace_back("CTG");
    //NT_triplets_.emplace_back("CCA");
    //NT_triplets_.emplace_back("CCT");
    //NT_triplets_.emplace_back("CCC");
    //NT_triplets_.emplace_back("CCG");
            NT_triplets_.emplace_back("CGA");
            NT_triplets_.emplace_back("CGT");
            NT_triplets_.emplace_back("CGC");
    //NT_triplets_.emplace_back("CGG");

    //NT_triplets_.emplace_back("GAA");
            NT_triplets_.emplace_back("GAT");
                    NT_triplets_.emplace_back("GAC");
            NT_triplets_.emplace_back("GAG");
            NT_triplets_.emplace_back("GTA");
    //NT_triplets_.emplace_back("GTT");
                    NT_triplets_.emplace_back("GTC");
            NT_triplets_.emplace_back("GTG");
            NT_triplets_.emplace_back("GCA");
            NT_triplets_.emplace_back("GCT");
    //NT_triplets_.emplace_back("GCC");
            NT_triplets_.emplace_back("GCG");
    //NT_triplets_.emplace_back("GGA");
    //NT_triplets_.emplace_back("GGT");
    //NT_triplets_.emplace_back("GGC");
    //NT_triplets_.emplace_back("GGG");

    list<string> candidates_row_0;
    candidates_row_0.emplace_back("TAC");
    candidates_row_0.emplace_back("GTA");
    candidates_row_0.emplace_back("CTG");
    candidates_row_0.emplace_back("CAG");
    candidates_row_0.emplace_back("ATA");
    three_bits_NT_triplets_candidates_.emplace(0,candidates_row_0);

    list<string> candidates_row_1;
    candidates_row_1.emplace_back("ATG");
    candidates_row_1.emplace_back("CAT");
    candidates_row_1.emplace_back("GAC");
    candidates_row_1.emplace_back("GTC");
    candidates_row_1.emplace_back("TAT");
    three_bits_NT_triplets_candidates_.emplace(1,candidates_row_1);

    list<string> candidates_row_2;
    candidates_row_2.emplace_back("ACT");
    candidates_row_2.emplace_back("AGT");
    candidates_row_2.emplace_back("TGC");
    candidates_row_2.emplace_back("GCA");
    candidates_row_2.emplace_back("CGC");
    three_bits_NT_triplets_candidates_.emplace(2,candidates_row_2);

    list<string> candidates_row_3;
    candidates_row_3.emplace_back("TCA");
    candidates_row_3.emplace_back("TGA");
    candidates_row_3.emplace_back("ACG");
    candidates_row_3.emplace_back("CGT");
    candidates_row_3.emplace_back("GCG");
    three_bits_NT_triplets_candidates_.emplace(3,candidates_row_3);

    list<string> candidates_row_4;
    candidates_row_4.emplace_back("ATC");
    candidates_row_4.emplace_back("GAT");
    candidates_row_4.emplace_back("TCG");
    candidates_row_4.emplace_back("CGA");
    three_bits_NT_triplets_candidates_.emplace(4,candidates_row_4);

    list<string> candidates_row_5;
    candidates_row_5.emplace_back("ACA");
    candidates_row_5.emplace_back("TGT");
    candidates_row_5.emplace_back("CAC");
    candidates_row_5.emplace_back("GTG");
    three_bits_NT_triplets_candidates_.emplace(5,candidates_row_5);

    list<string> candidates_row_6;
    candidates_row_6.emplace_back("TAG");
    candidates_row_6.emplace_back("CTA");
    candidates_row_6.emplace_back("AGC");
    candidates_row_6.emplace_back("GCT");
    three_bits_NT_triplets_candidates_.emplace(6,candidates_row_6);

    list<string> candidates_row_7;
    candidates_row_7.emplace_back("AGA");
    candidates_row_7.emplace_back("TCT");
    candidates_row_7.emplace_back("CTC");
    candidates_row_7.emplace_back("GAG");
    three_bits_NT_triplets_candidates_.emplace(7,candidates_row_7);

    complementary_triplets_table_.emplace("TAC","GTA");
    complementary_triplets_table_.emplace("GTA","TAC");
    complementary_triplets_table_.emplace("ATG","CAT");
    complementary_triplets_table_.emplace("CAT","ATG");
    complementary_triplets_table_.emplace("ATC","GAT");
    complementary_triplets_table_.emplace("GAT","ATC");
    complementary_triplets_table_.emplace("ACA","TGT");
    complementary_triplets_table_.emplace("TGT","ACA");
    complementary_triplets_table_.emplace("ACT","AGT");
    complementary_triplets_table_.emplace("AGT","ACT");
    complementary_triplets_table_.emplace("TCA","TGA");
    complementary_triplets_table_.emplace("TGA","TCA");
    complementary_triplets_table_.emplace("TAG","CTA");
    complementary_triplets_table_.emplace("CTA","TAG");


    complementary_triplets_table_.emplace("AGA","TCT");
    complementary_triplets_table_.emplace("TCT","AGA");
    complementary_triplets_table_.emplace("CTG","CAG");
    complementary_triplets_table_.emplace("CAG","CTG");
    complementary_triplets_table_.emplace("TCG","CGA");
    complementary_triplets_table_.emplace("CGA","TCG");
    complementary_triplets_table_.emplace("TGC","GCA");
    complementary_triplets_table_.emplace("GCA","TGC");
    complementary_triplets_table_.emplace("AGC","GCT");
    complementary_triplets_table_.emplace("GCT","AGC");
    complementary_triplets_table_.emplace("CAC","GTG");
    complementary_triplets_table_.emplace("GTG","CAC");
    complementary_triplets_table_.emplace("GAC","GTC");
    complementary_triplets_table_.emplace("GTC","GAC");
    complementary_triplets_table_.emplace("ACG","CGT");
    complementary_triplets_table_.emplace("CGT","ACG");
    complementary_triplets_table_.emplace("CTC","GAG");
    complementary_triplets_table_.emplace("GAG","CTC");
    complementary_triplets_table_.emplace("ATA","TAT");
    complementary_triplets_table_.emplace("TAT","ATA");
    complementary_triplets_table_.emplace("CGC","GCG");
    complementary_triplets_table_.emplace("GCG","CGC");

}

void DNA_encoder::initial_rotating_encoding_table() {
    /*
     *  [ 0
     *      [A: C
     *      [T: A
     *      [C: G
     *      [G: T
     *  [ 1
     *      [A: G
     *      [T: C
     *      [C: T
     *      [G: A
     *  [ 2
     *      [A: T
     *      [T: G
     *      [C: A
     *      [G: C
     * */
    unordered_map<string,string> bit_value_0;
    unordered_map<string,string> bit_value_1;
    unordered_map<string,string> bit_value_2;
    string A="A";
    string T="T";
    string C="C";
    string G="G";
    bit_value_0.emplace(A,C);
    bit_value_0.emplace(T,A);
    bit_value_0.emplace(C,G);
    bit_value_0.emplace(G,T);
    rotating_encoding_table_.push_back(bit_value_0);
    bit_value_1.emplace(A,G);
    bit_value_1.emplace(T,C);
    bit_value_1.emplace(C,T);
    bit_value_1.emplace(G,A);
    rotating_encoding_table_.push_back(bit_value_1);
    bit_value_2.emplace(A,T);
    bit_value_2.emplace(T,G);
    bit_value_2.emplace(C,A);
    bit_value_2.emplace(G,C);
    rotating_encoding_table_.push_back(bit_value_2);
}

void DNA_encoder::initial_twobits_rotating_encoding_table() {
    /*
     *  [ 0
     *      [A: C
     *      [T: A
     *      [C: G
     *      [G: T
     *  [ 1
     *      [A: G
     *      [T: C
     *      [C: T
     *      [G: A
     *  [ 2
     *      [A: T
     *      [T: G
     *      [C: A
     *      [G: C
     * */
    unordered_map<string,string> bit_value_00;
    unordered_map<string,string> bit_value_01;
    unordered_map<string,string> bit_value_10;
    unordered_map<string,string> bit_value_11;
    string A="A";
    string T="T";
    string C="C";
    string G="G";
    bit_value_00.emplace(A,T);
    bit_value_00.emplace(T,C);
    bit_value_00.emplace(C,G);
    bit_value_00.emplace(G,A);
    rotating_encoding_table_.push_back(bit_value_00);
    bit_value_01.emplace(A,C);
    bit_value_01.emplace(T,G);
    bit_value_01.emplace(C,A);
    bit_value_01.emplace(G,T);
    rotating_encoding_table_.push_back(bit_value_01);
    bit_value_10.emplace(A,G);
    bit_value_10.emplace(T,A);
    bit_value_10.emplace(C,T);
    bit_value_10.emplace(G,C);
    rotating_encoding_table_.push_back(bit_value_10);
    bit_value_11.emplace(A,A);
    bit_value_11.emplace(T,T);
    bit_value_11.emplace(C,C);
    bit_value_11.emplace(G,G);
    rotating_encoding_table_.push_back(bit_value_10);
}

void DNA_encoder::init_GF47_table() {
    RS_table[0] = "ACA"; RS_table[1] = "CCA"; RS_table[2] = "GCA";
    RS_table[3] = "TCA"; RS_table[4] = "AGA"; RS_table[5] = "CGA";
    RS_table[6] = "GGA"; RS_table[7] = "TGA"; RS_table[8] = "ATA";
    RS_table[9] = "CTA"; RS_table[10] = "GTA"; RS_table[11] = "TTA";
    RS_table[12] = "AAC"; RS_table[13] = "CAC"; RS_table[14] = "GAC";
    RS_table[15] = "TAC"; RS_table[16] = "AGC"; RS_table[17] = "CGC";
    RS_table[18] = "GGC"; RS_table[19] = "TGC"; RS_table[20] = "ATC";
    RS_table[21] = "CTC"; RS_table[22] = "GTC"; RS_table[23] = "TTC";
    RS_table[24] = "AAG"; RS_table[25] = "CAG"; RS_table[26] = "GAG";
    RS_table[27] = "TAG"; RS_table[28] = "ACG"; RS_table[29] = "CCG";
    RS_table[30] = "GCG"; RS_table[31] = "TCG"; RS_table[32] = "ATG";
    RS_table[33] = "CTG"; RS_table[34] = "GCG"; RS_table[35] = "TTG";
    RS_table[36] = "AAT"; RS_table[37] = "CAT"; RS_table[38] = "GAT";
    RS_table[39] = "TAT"; RS_table[40] = "ACT"; RS_table[41] = "CCT";
    RS_table[42] = "GCT"; RS_table[43] = "TCT"; RS_table[44] = "AGT";
    RS_table[45] = "CGT"; RS_table[46] = "GGT";
}

void DNA_encoder::init_GF25_table() {
    RS_table[0] = "GCA";
    RS_table[1] = "TCA"; RS_table[2] = "CGA";
     RS_table[3] = "TGA";
    RS_table[4] = "CTA"; RS_table[5] = "GTA";
      RS_table[6] = "GAC";
    RS_table[7] = "TAC"; RS_table[8] = "AGC";
     RS_table[9] = "TGC"; RS_table[10] = "ATC";
     RS_table[11] = "GTC";
     RS_table[12] = "CAG";
    RS_table[13] = "TAG"; RS_table[14] = "ACG";
     RS_table[15] = "TCG"; RS_table[16] = "ATG";
    RS_table[17] = "CTG";
     RS_table[18] = "CAT"; RS_table[19] = "GAT";
    RS_table[20] = "TAT"; RS_table[21] = "ACT";
    RS_table[22] = "GCT";  RS_table[23] = "AGT";
    RS_table[24] = "CGT";
}

void DNA_encoder::init_last_20nt(){
    /*for(int i=0; i<5; i++){
        last_20nt_.push_back("A");
        last_20nt_.push_back("G");
        last_20nt_.push_back("T");
        last_20nt_.push_back("C");
    }*/

    last_20nt_.emplace_back("G");
    last_20nt_.emplace_back("A");
    last_20nt_.emplace_back("C");
    last_20nt_.emplace_back("A");
    last_20nt_.emplace_back("T");
    last_20nt_.emplace_back("C");


    /*last_20nt_.emplace_back("C");
    last_20nt_.emplace_back("C");

    last_20nt_.emplace_back("A");
    last_20nt_.emplace_back("A");
    last_20nt_.emplace_back("A");
    last_20nt_.emplace_back("G");
    last_20nt_.emplace_back("G");
    last_20nt_.emplace_back("G");
    last_20nt_.emplace_back("T");
    last_20nt_.emplace_back("T");
    last_20nt_.emplace_back("T");
    last_20nt_.emplace_back("C");
    last_20nt_.emplace_back("C");
    last_20nt_.emplace_back("C");*/
}

DNA_encoder::DNA_encoder() {
    // generate 267-length pseudo random. 200 ternary code is 266.6 binary code (6:8)
    srand(std::time(nullptr)); // use current time as seed for random generator
    for (int n=0; n < 267; ++n) {
        pseudo_random_sequence_+=std::rand()%2;
    }


    // record all file's path, use to read file & encode it
    all_files_ = listFiles(g_data_path, true);

    for(auto n:all_files_){
        cout<<n<<endl;
    }
    //initial rotating_encoding_table
    initial_rotating_encoding_table();
    initial_twobits_rotating_encoding_table();
    initial_FEC_table();
    init_last_20nt();
    init_heuristic_encoding();
    //initilize RS table
    //init_GF47_table();
    init_GF25_table();

    if (g_if_pre_stranding){// output 200 nts strands
        encoding_stranding();
    }
    else{
        encoding_no_strand();
    }
    //encode with strand;
    //encoding_stranding();
    //encode without strand
    //encoding_no_strand();
    //encoding_file();
}
