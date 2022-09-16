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
#include <cmath>
#include <ctime>
#include <cstdlib>



string DNA_encoder::Church_encoding(string digital_data) {
    while (digital_data.size() % 30 != 0){
        digital_data = digital_data + '\0';
    }
    //cout << "digital_data.size() = " << digital_data.size() << endl;

    // step 1: apply 2 parity bytes on each 30-bytes block
    const int MSG_LENGTH = 30;
    const int ECC_LENGTH = 2;
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


    string result;
    string last_three="ATC";
    int overall_gc=0;
    // 0 - A/C  1- T/G
    for (std::size_t i = 0; i < sizeof(rs_result); i++)
    {
        //cover from binary to decimal
        bitset<8> bits(rs_result[i]);
        for (int j = 0; j < 8; j++) {
            if (bits[j]==0){
                if (last_three[0]=='A' && last_three[1]=='A' && last_three[2]=='A'){
                    result+="C";
                    last_three.erase(0,1);
                    last_three+="C";
                    overall_gc++;
                } else if(last_three[0]=='C' && last_three[1]=='C' && last_three[2]=='C'){
                    result+="A";
                    last_three.erase(0,1);
                    last_three+="A";
                    overall_gc--;
                } else{
                    string str;
                    if (overall_gc>=0){
                        str='A';
                        overall_gc--;
                    } else{
                        str='C';
                        overall_gc++;
                    }
                    last_three.erase(0,1);
                    last_three+=str;
                    result+=str;
                }
            }
            else if (bits[j]==1){
                if (last_three[0]=='T' && last_three[1]=='T' && last_three[2]=='T'){
                    result+="G";
                    last_three.erase(0,1);
                    last_three+="G";
                    overall_gc++;
                } else if(last_three[0]=='G' && last_three[1]=='G' && last_three[2]=='G'){
                    result+="T";
                    last_three.erase(0,1);
                    last_three+="T";
                    overall_gc--;
                } else{
                    string str;
                    if (overall_gc>=0){
                        str='T';
                        overall_gc--;
                    } else{
                        str='G';
                        overall_gc++;
                    }
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


string DNA_encoder::RS_rotation(string digital_data) {
    while (digital_data.size() % 30 != 0){
        digital_data = digital_data + '\0';
    }
    //cout << "digital_data.size() = " << digital_data.size() << endl;

    // step 1: apply 2 parity bytes on each 30-bytes block
    const int MSG_LENGTH = 30;
    const int ECC_LENGTH = 2;
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
    while (digital_data.size() % 30 != 0){
        digital_data = digital_data + '\0';
    }
    //cout << "digital_data.size() = " << digital_data.size() << endl;

    // step 1: apply 2 parity bytes on each 30-bytes block
    const int MSG_LENGTH = 30;
    const int ECC_LENGTH = 2;
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

    int result_len = sizeof(rs_result) / 2 * 3 * 3; // 2 -> two 256-base number to use; 3 -> three 47-base number generated; 3 -> three nt for each 47-base number
    //cout << "result_len = " << result_len << endl;
    char result_cstring[result_len+1];
    result_cstring[result_len] = '\0';

    int j = 0;
    for (int i = 0; i < sizeof(rs_result); i += 2) {
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

        for (int k1 = 0; k1 < 3; k1++) {
            for (int k2 = 0; k2 < 3; k2++) {
                result_cstring[j++] = RS_table[base47_num[k1]][k2];
            }
        }
    }
    assert(j == result_len);
    string result(result_cstring);
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

string DNA_encoder::exhaustive_search_triplets(int bits_triplet) {
    string result;
    vector<pair<double,string>> all_triplet_candidates_and_score;


    list<string> homo1;
    //triplets that has no homo with the last nt
    for (auto n:NT_triplets_){
        if (last_17nt_.back()!=n[0])
            homo1.push_back(n);
    }

    int GC=0;
    for (auto s:last_17nt_){
        if (s=='G' || s=='C') GC++;
    }

    for (auto m:homo1){
        double matches=0;
        string complementary_sequence = complementary_triplets_table_.find(m)->second;
        auto last = last_17nt_.end();
        last--;
        if (*last=='A')  complementary_sequence.insert (0, 1, 'T');
        else if (*last=='T') complementary_sequence.insert (0, 1, 'A');
        else if (*last=='C') complementary_sequence.insert (0, 1, 'G');
        else if (*last=='G') complementary_sequence.insert (0, 1, 'C');
        last--;
        if (*last=='A')  complementary_sequence.insert (0, 1, 'T');
        else if (*last=='T') complementary_sequence.insert (0, 1, 'A');
        else if (*last=='C') complementary_sequence.insert (0, 1, 'G');
        else if (*last=='G') complementary_sequence.insert (0, 1, 'C');
        // complementary_sequence of (last two base + current triplet)


        //start from first base of complementary_sequence
        auto it = last_17nt_.begin();
        for (int i = 0; i < 14; ++i) {
            auto tmp_it = it;
            if (*tmp_it != complementary_sequence[0]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complementary_sequence[1]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complementary_sequence[2]) {
                it++;
                continue;
            }
            matches+=3;
            //extend
            tmp_it++;
            if (*tmp_it != complementary_sequence[3]) {
                it++;
                continue;
            }
            matches+=4;
            break;
        }


        //start from 2nd base of complementary_sequence
        it = last_17nt_.begin();
        for (int i = 0; i < 14; ++i) {
            auto tmp_it = it;
            if (*tmp_it != complementary_sequence[1]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complementary_sequence[2]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complementary_sequence[3]) {
                it++;
                continue;
            }
            matches+=3;
            //extend
            tmp_it++;
            if (*tmp_it != complementary_sequence[4]) {
                it++;
                continue;
            }
            matches+=4;
            break;
        }

        //start from 3nd base of complementary_sequence (the first base of the triplet, only 3 bases remain)
        it = last_17nt_.begin();
        for (int i = 0; i < 14; ++i) {
            auto tmp_it = it;
            if (*tmp_it != complementary_sequence[2]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complementary_sequence[3]) {
                it++;
                continue;
            }
            tmp_it++;
            if (*tmp_it != complementary_sequence[4]) {
                it++;
                continue;
            }
            matches+=3;
            break;
        }

        // compare GC to further tune each candidates's weight (but in finer grain)
        int GC_in_m=0;
        for (auto l:m){
            if (l=='G' || l=='C') GC_in_m++;
        }
        matches += 0.1*abs(GC_in_m+GC-10);
        all_triplet_candidates_and_score.push_back(make_pair(matches,m));
    }

    sort(all_triplet_candidates_and_score.begin(), all_triplet_candidates_and_score.end());

    result=all_triplet_candidates_and_score[bits_triplet].second;

    //replace last 17
    for (int i = 0; i < result.size(); ++i) {
        last_17nt_.emplace_back(result[i]);
        last_17nt_.erase(last_17nt_.begin());
    }
    return result;
}

static void multithread_complementary_GC_checking(list<string>& subsequences, int GC, string candidate, list<char>& last_17nt_, unordered_map<string,string>& complementary_triplets_table_, unordered_map<string,uint8_t>& look_back_window_, double& weight)
{
        bool match1= false;
        bool match2= false;
        bool match3= false;
        bool match4= false;
        bool match5= false;



        //last two base + first base of the triplet
        string str1;
        str1.insert (0, 1, candidate[0]);
        auto it = last_17nt_.end();
        it--;
        str1.insert (0, 1, *it);
        it--;
        str1.insert (0, 1, *it);
        subsequences.emplace_back(str1);

        string complement_str1 = complementary_triplets_table_.find(str1)->second;
        if (look_back_window_.find(complement_str1)!=look_back_window_.end())
            match1= true;

        //last base + first two base of the triplet
        string str2;
        str2.insert (0, 1, candidate[1]);
        str2.insert (0, 1, candidate[0]);
        it = last_17nt_.end();
        it--;
        str2.insert (0, 1, *it);
        subsequences.emplace_back(str2);

        string complement_str2 = complementary_triplets_table_.find(str2)->second;
        if (look_back_window_.find(complement_str2)!=look_back_window_.end())
            match2= true;

        //the triplet
        subsequences.emplace_back(candidate);
        string complement_str3 = complementary_triplets_table_.find(candidate)->second;
        if (look_back_window_.find(complement_str3)!=look_back_window_.end())
            match3= true;


        //last base + triplet
        string str4 = candidate;
        it = last_17nt_.end();
        it--;
        str4.insert (0, 1, *it);
        subsequences.emplace_back(str4);

        if (match2== true && match3== true){
            string complement_str4=complement_str3;
            it = last_17nt_.end();
            it--;
            if (*it=='A')  complement_str4.insert (0, 1, 'T');
            else if (*it=='T') complement_str4.insert (0, 1, 'A');
            else if (*it=='C') complement_str4.insert (0, 1, 'G');
            else if (*it=='G') complement_str4.insert (0, 1, 'C');

            if (look_back_window_.find(complement_str4)!=look_back_window_.end())
                match4= true;
        }

        //last two bases + first two bases of triplet
        string str5 = str1;
        str5.insert (3, 1, candidate[1]);
        subsequences.emplace_back(str5);

        if (match1== true && match2== true){
            string complement_str5=complement_str1;
            if (candidate[1]=='A')  complement_str5.insert (3, 1, 'T');
            else if (candidate[1]=='T') complement_str5.insert (3, 1, 'A');
            else if (candidate[1]=='C') complement_str5.insert (3, 1, 'G');
            else if (candidate[1]=='G') complement_str5.insert (3, 1, 'C');

            if (look_back_window_.find(complement_str5)!=look_back_window_.end())
                match5= true;
        }
        double matches=0;
        if (match1) matches+=3;
        if (match2) matches+=3;
        if (match3) matches+=3;
        if (match4) matches+=4;  // ensure long complementary sequence has larger weight, even if count redundantly
        if (match5) matches+=4;


        // compare GC to further tune each candidates's weight (but in finer grain)
        int GC_in_m=0;
        for (auto l:candidate){
            if (l=='G' || l=='C') GC_in_m++;
        }

        weight = matches + 0.1*abs(GC_in_m+GC-10);
}

string DNA_encoder::multithread_partial_search_triplets(int bits_triplet) {
    string result;
    vector<pair<double,string>> all_triplet_candidates_and_score;
    unordered_map<string,list<string>> subsequences_of_each_candidates;
    list<string> candidates;

    if(g_num_bit_per_triplet==2){
        candidates = two_bits_NT_triplets_candidates_[bits_triplet];
    }else if (g_num_bit_per_triplet==3){
        candidates = three_bits_NT_triplets_candidates_[bits_triplet];
    }else if (g_num_bit_per_triplet==4){
        candidates = four_bits_NT_triplets_candidates_[bits_triplet];
    }

    vector<string> homo1;
    //homo
    for (auto n:candidates){
        if (last_17nt_.back()!=n[0])
            homo1.push_back(n);
    }

    int GC=0;
    for (auto s:last_17nt_){
        if (s=='G' || s=='C') GC++;
    }

    //at most 4 candidates
    list<string> subsequences_candidate_1;
    double weight_1=10000; //default a large unrealistic weight
    list<string> subsequences_candidate_2;
    double weight_2=10000;
    list<string> subsequences_candidate_3;
    double weight_3=10000;
    list<string> subsequences_candidate_4;
    double weight_4=10000;

    //one thread one candidate
   if (homo1.size()==4){
        clock_t time_tt1 = clock();
        thread th1(multithread_complementary_GC_checking, ref(subsequences_candidate_1), GC, homo1[0], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_1));
        thread th2(multithread_complementary_GC_checking, ref(subsequences_candidate_2), GC, homo1[1], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_2));
        thread th3(multithread_complementary_GC_checking, ref(subsequences_candidate_3), GC, homo1[2], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_3));
        thread th4(multithread_complementary_GC_checking, ref(subsequences_candidate_4), GC, homo1[3], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_4));
        th1.join(); clock_t time_tt2 = clock();
        th2.join(); clock_t time_tt3 = clock();
        th3.join(); clock_t time_tt4 = clock();
        th4.join(); clock_t time_tt5 = clock();

       cout<<time_tt2-time_tt1<<endl;
       cout<<time_tt3-time_tt2<<endl;
       cout<<time_tt4-time_tt3<<endl;
       cout<<time_tt5-time_tt4<<endl;
       cout<<"total: "<<time_tt3-time_tt1<<endl;
   }
    else if (homo1.size()==3){
        thread th1(multithread_complementary_GC_checking, ref(subsequences_candidate_1), GC, homo1[0], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_1));
        thread th2(multithread_complementary_GC_checking, ref(subsequences_candidate_2), GC, homo1[1], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_2));
        thread th3(multithread_complementary_GC_checking, ref(subsequences_candidate_3), GC, homo1[2], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_3));
        th1.join();
        th2.join();
        th3.join();
    }
    else if (homo1.size()==2){
        thread th1(multithread_complementary_GC_checking, ref(subsequences_candidate_1), GC, homo1[0], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_1));
        thread th2(multithread_complementary_GC_checking, ref(subsequences_candidate_2), GC, homo1[1], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_2));
        th1.join();
        th2.join();
    } else {
        thread th1(multithread_complementary_GC_checking, ref(subsequences_candidate_1), GC, homo1[0], ref(last_17nt_), ref(complementary_triplets_table_), ref(look_back_window_), ref(weight_1));
        th1.join();
    }

    if (homo1.size()==4){
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_1,homo1[0]));
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_2,homo1[1]));
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_3,homo1[2]));
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_4,homo1[3]));

        subsequences_of_each_candidates.emplace(homo1[0],subsequences_candidate_1);
        subsequences_of_each_candidates.emplace(homo1[1],subsequences_candidate_2);
        subsequences_of_each_candidates.emplace(homo1[2],subsequences_candidate_3);
        subsequences_of_each_candidates.emplace(homo1[3],subsequences_candidate_4);
    }
    else if (homo1.size()==3){
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_1,homo1[0]));
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_2,homo1[1]));
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_3,homo1[2]));

        subsequences_of_each_candidates.emplace(homo1[0],subsequences_candidate_1);
        subsequences_of_each_candidates.emplace(homo1[1],subsequences_candidate_2);
        subsequences_of_each_candidates.emplace(homo1[2],subsequences_candidate_3);
    }
    else if (homo1.size()==2){
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_1,homo1[0]));
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_2,homo1[1]));

        subsequences_of_each_candidates.emplace(homo1[0],subsequences_candidate_1);
        subsequences_of_each_candidates.emplace(homo1[1],subsequences_candidate_2);
    } else {
        all_triplet_candidates_and_score.emplace_back(make_pair(weight_1,homo1[0]));

        subsequences_of_each_candidates.emplace(homo1[0],subsequences_candidate_1);
    }

    sort(all_triplet_candidates_and_score.begin(), all_triplet_candidates_and_score.end());

    result=all_triplet_candidates_and_score.begin()->second;

    //replace last 17
    for (int i = 0; i < result.size(); ++i) {
        last_17nt_.emplace_back(result[i]);
        last_17nt_.erase(last_17nt_.begin());
    }

    //insert new subsequences to look back window
    list<string> subsequences_of_select_triplets = subsequences_of_each_candidates.find(result)->second;
    for (auto s:subsequences_of_select_triplets) {
        if (look_back_window_.find(s)==look_back_window_.end())
            look_back_window_.emplace(s,1);
        else
            look_back_window_.find(s)->second++;
    }

    //delete old subsequence from look back window
    auto it = last_17nt_.begin();
    string delete_str1, delete_str2, delete_str3, delete_str4, delete_str5, delete_str6;
    delete_str1.push_back(*it);

    it++;
    delete_str1.push_back(*it);
    delete_str2.push_back(*it);
    delete_str4.push_back(*it);


    it++;
    delete_str1.push_back(*it);
    delete_str2.push_back(*it);
    delete_str3.push_back(*it);
    delete_str4.push_back(*it);
    delete_str5.push_back(*it);

    it++;
    delete_str2.push_back(*it);
    delete_str3.push_back(*it);
    delete_str4.push_back(*it);
    delete_str5.push_back(*it);

    it++;
    delete_str3.push_back(*it);
    delete_str4.push_back(*it);
    delete_str5.push_back(*it);

    it++;
    delete_str5.push_back(*it);

    auto ptr1 = look_back_window_.find(delete_str1);
    if (ptr1!=look_back_window_.end()){
        if (ptr1->second==1)
            look_back_window_.erase(ptr1);
        else
            ptr1->second--;
    }

    auto ptr2 = look_back_window_.find(delete_str2);
    if (ptr2!=look_back_window_.end()){
        if (ptr2->second==1)
            look_back_window_.erase(ptr2);
        else
            ptr2->second--;
    }

    auto ptr3 = look_back_window_.find(delete_str3);
    if (ptr3!=look_back_window_.end()){
        if (ptr3->second==1)
            look_back_window_.erase(ptr3);
        else
            ptr3->second--;
    }

    auto ptr4 = look_back_window_.find(delete_str4);
    if (ptr4!=look_back_window_.end()){
        if (ptr4->second==1)
            look_back_window_.erase(ptr4);
        else
            ptr4->second--;
    }

    auto ptr5 = look_back_window_.find(delete_str5);
    if (ptr5!=look_back_window_.end()){
        if (ptr5->second==1)
            look_back_window_.erase(ptr5);
        else
            ptr5->second--;
    }

    return result;
}

string DNA_encoder::partial_search_triplets(int bits_triplet) {
    string result;
    vector<pair<double,string>> all_triplet_candidates_and_score;
    unordered_map<string,list<string>> subsequences_of_each_candidates;
    list<string> candidates;

    if(g_num_bit_per_triplet==2){
        candidates = two_bits_NT_triplets_candidates_[bits_triplet];
    }else if (g_num_bit_per_triplet==3){
        candidates = three_bits_NT_triplets_candidates_[bits_triplet];
    }else if (g_num_bit_per_triplet==4){
        candidates = four_bits_NT_triplets_candidates_[bits_triplet];
    }

    vector<string> homo1;
    //homo
    for (auto n:candidates){
        if (last_17nt_.back()!=n[0])
            homo1.push_back(n);
    }

    int GC=0;
    for (auto s:last_17nt_){
        if (s=='G' || s=='C') GC++;
    }

    clock_t time_tt1, time_tt2;

    time_tt1 = clock();
    for (auto candidate:homo1){


        bool match1= false;
        bool match2= false;
        bool match3= false;
        bool match4= false;
        bool match5= false;
        list<string> subsequences;
        double weight=10000; //default a large unrealistic weight


        //last two base + first base of the triplet
        string str1;
        str1.insert (0, 1, candidate[0]);
        auto it = last_17nt_.end();
        it--;
        str1.insert (0, 1, *it);
        it--;
        str1.insert (0, 1, *it);
        subsequences.emplace_back(str1);

        string complement_str1 = complementary_triplets_table_.find(str1)->second;
        if (look_back_window_.find(complement_str1)!=look_back_window_.end())
            match1= true;

        //last base + first two base of the triplet
        string str2;
        str2.insert (0, 1, candidate[1]);
        str2.insert (0, 1, candidate[0]);
        it = last_17nt_.end();
        it--;
        str2.insert (0, 1, *it);
        subsequences.emplace_back(str2);

        string complement_str2 = complementary_triplets_table_.find(str2)->second;
        if (look_back_window_.find(complement_str2)!=look_back_window_.end())
            match2= true;

        //the triplet
        subsequences.emplace_back(candidate);
        string complement_str3 = complementary_triplets_table_.find(candidate)->second;
        if (look_back_window_.find(complement_str3)!=look_back_window_.end())
            match3= true;


        //last base + triplet
        string str4 = candidate;
        it = last_17nt_.end();
        it--;
        str4.insert (0, 1, *it);
        subsequences.emplace_back(str4);

        if (match2== true && match3== true){
            string complement_str4=complement_str3;
            it = last_17nt_.end();
            it--;
            if (*it=='A')  complement_str4.insert (0, 1, 'T');
            else if (*it=='T') complement_str4.insert (0, 1, 'A');
            else if (*it=='C') complement_str4.insert (0, 1, 'G');
            else if (*it=='G') complement_str4.insert (0, 1, 'C');

            if (look_back_window_.find(complement_str4)!=look_back_window_.end())
                match4= true;
        }

        //last two bases + first two bases of triplet
        string str5 = str1;
        str5.insert (3, 1, candidate[1]);
        subsequences.emplace_back(str5);

        if (match1== true && match2== true){
            string complement_str5=complement_str1;
            if (candidate[1]=='A')  complement_str5.insert (3, 1, 'T');
            else if (candidate[1]=='T') complement_str5.insert (3, 1, 'A');
            else if (candidate[1]=='C') complement_str5.insert (3, 1, 'G');
            else if (candidate[1]=='G') complement_str5.insert (3, 1, 'C');

            if (look_back_window_.find(complement_str5)!=look_back_window_.end())
                match5= true;
        }
        double matches=0;
        if (match1) matches+=3;
        if (match2) matches+=3;
        if (match3) matches+=3;
        if (match4) matches+=4;  // ensure long complementary sequence has larger weight, even if count redundantly
        if (match5) matches+=4;


        // compare GC to further tune each candidates's weight (but in finer grain)
        int GC_in_m=0;
        for (auto l:candidate){
            if (l=='G' || l=='C') GC_in_m++;
        }

        weight = matches + 0.1*abs(GC_in_m+GC-10);;

        all_triplet_candidates_and_score.emplace_back(make_pair(weight,candidate));
        subsequences_of_each_candidates.emplace(candidate,subsequences);
    }
        time_tt2 = clock();
        cout<<"total "<<homo1.size()<<": "<<time_tt2-time_tt1<<endl;

    sort(all_triplet_candidates_and_score.begin(), all_triplet_candidates_and_score.end());

    result=all_triplet_candidates_and_score.begin()->second;


    //replace last 17
    for (int i = 0; i < result.size(); ++i) {
        last_17nt_.emplace_back(result[i]);
        last_17nt_.erase(last_17nt_.begin());
    }

    //insert new subsequences to look back window
    list<string> subsequences_of_select_triplets = subsequences_of_each_candidates.find(result)->second;
    for (auto s:subsequences_of_select_triplets) {
        if (look_back_window_.find(s)==look_back_window_.end())
            look_back_window_.emplace(s,1);
        else
            look_back_window_.find(s)->second++;
    }

    //delete old subsequence from look back window
    auto it = last_17nt_.begin();
    string delete_str1, delete_str2, delete_str3, delete_str4, delete_str5, delete_str6;
    delete_str1.push_back(*it);

    it++;
    delete_str1.push_back(*it);
    delete_str2.push_back(*it);
    delete_str4.push_back(*it);


    it++;
    delete_str1.push_back(*it);
    delete_str2.push_back(*it);
    delete_str3.push_back(*it);
    delete_str4.push_back(*it);
    delete_str5.push_back(*it);

    it++;
    delete_str2.push_back(*it);
    delete_str3.push_back(*it);
    delete_str4.push_back(*it);
    delete_str5.push_back(*it);

    it++;
    delete_str3.push_back(*it);
    delete_str4.push_back(*it);
    delete_str5.push_back(*it);

    it++;
    delete_str5.push_back(*it);

    auto ptr1 = look_back_window_.find(delete_str1);
    if (ptr1!=look_back_window_.end()){
        if (ptr1->second==1)
            look_back_window_.erase(ptr1);
        else
            ptr1->second--;
    }

    auto ptr2 = look_back_window_.find(delete_str2);
    if (ptr2!=look_back_window_.end()){
        if (ptr2->second==1)
            look_back_window_.erase(ptr2);
        else
            ptr2->second--;
    }

    auto ptr3 = look_back_window_.find(delete_str3);
    if (ptr3!=look_back_window_.end()){
        if (ptr3->second==1)
            look_back_window_.erase(ptr3);
        else
            ptr3->second--;
    }

    auto ptr4 = look_back_window_.find(delete_str4);
    if (ptr4!=look_back_window_.end()){
        if (ptr4->second==1)
            look_back_window_.erase(ptr4);
        else
            ptr4->second--;
    }

    auto ptr5 = look_back_window_.find(delete_str5);
    if (ptr5!=look_back_window_.end()){
        if (ptr5->second==1)
            look_back_window_.erase(ptr5);
        else
            ptr5->second--;
    }

    return result;
}


string DNA_encoder::heuristic_encoding(string digital_data) {
    while (digital_data.size() % 30 != 0){
        digital_data = digital_data + '\0';
    }
    //cout << "digital_data.size() = " << digital_data.size() << endl;

    // step 1: apply 2 parity bytes on each 30-bytes block
    const int MSG_LENGTH = 30;
    const int ECC_LENGTH = 2;
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

    string result;
    int size = digital_data.size();
        for (int j = 0; j < sizeof(rs_result); j+=3) {// 3 bytes (24 bits) a time,
            // convert from 256^3 to 8^8
            unsigned int middle_num = ((((unsigned int)rs_result[j])&0xff) * 256 + (((unsigned int)rs_result[j+1])&0xff))* 256 + (((unsigned int)rs_result[j+2])&0xff);
            // cout << "l r m = " << hex << (unsigned int)(char)(rs_result[i]) << '\t' << (unsigned int)(rs_result[i+1]) << '\t' << middle_num << endl << flush;
            unsigned int base8_num[8];
            base8_num[7] = middle_num % 8;
            middle_num /= 8;
            base8_num[6] = middle_num % 8;
            middle_num /= 8;
            base8_num[5] = middle_num % 8;
            middle_num /= 8;
            base8_num[4] = middle_num % 8;
            middle_num /= 8;
            base8_num[3] = middle_num % 8;
            middle_num /= 8;
            base8_num[2] = middle_num % 8;
            middle_num /= 8;
            base8_num[1] = middle_num % 8;
            middle_num /= 8;
            base8_num[0] = middle_num % 8;
            middle_num /= 8;
            assert(middle_num == 0);

            // need to change if using 4-bits or 2-bits
            for (int i = 0; i < 8; i++) {
                if (g_encoding_scheme==5)
                    result+=exhaustive_search_triplets(base8_num[i]);
                if (g_encoding_scheme==6)
                    result+=multithread_partial_search_triplets(base8_num[i]);
                    //result+=partial_search_triplets(base8_num[i]);
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



void DNA_encoder::encoding_stranding(){
    fstream payload_file;
    int prefix = rand();
    string payload_path = g_payload_path+to_string(prefix)+"payload";
    int file_number = 1;
    unsigned long total_len = 0;
    payload_path += to_string(file_number);
    payload_file.open(payload_path,ios::out);

    long int strand_num=0;
    uint8_t buf[10*1024];
    FILE *fp;
    string nt_sequence;

    for(auto n:all_files_){
        fp = fopen(n.c_str(), "r");
        if (fp==NULL) {fputs ("File open error",stderr); exit (1);}
        while ( !feof(fp) ) {
            if (total_len>=(long)10*1024*1024){
                total_len = 0;
                payload_file.close();
                file_number += 1;
                payload_path = g_payload_path+to_string(prefix)+"payload"+to_string(file_number);
                payload_file.open(payload_path,ios::out);
            }
            size_t len = fread(buf, 1, sizeof(buf), fp);
            uint8_t *ptr = &buf[0];
            string digital_data ((char*)ptr,len);
            // rotate code

            string strand;
            if(g_encoding_scheme==1) {
                strand = RS_rotation(digital_data);
            }
            else if(g_encoding_scheme==0){
                strand=Church_encoding(digital_data);
            }
            else if(g_encoding_scheme==2){
                // 200 nts correspoinding to 320 bits
                // 2 nts have to be reserved for ECC (CCITT16: 16 bits)
                // each oligo store information with size of 320-16 = 304 bits -> 38 char
                while (digital_data.size()>0){
                    string digital_strand;
                    
                    if (g_if_ECC){//FEC has its own ecc except RS code
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
                    
                    strand = FEC_encoding(digital_strand);
                    if (g_if_ECC){
                        digital_data.erase(0, 38);
                    }
                    else {
                        digital_data.erase(0, 40);
                    }
                }
            } else if (g_encoding_scheme==3){
                strand=RS_GF47(digital_data);
            }
            else if(g_encoding_scheme==5 || g_encoding_scheme==6 ){
                strand=heuristic_encoding(digital_data);
            }
            else
                cout<<"no encoding scheme"<<endl;

            total_len += sizeof(buf);
            for (int l=0; l<strand.size()/198; l++){
                payload_file << ">payload" << strand_num++ << endl;
                payload_file << strand.substr(l*198,198) << endl;
            }
            if (strand.size()%198>0){
                payload_file << ">payload" << strand_num++ << endl;
                payload_file << strand.substr(strand.size()/198*198,strand.size()%198) << endl;
            }

            // execute transformation: mapping/swap/...


        }
        fclose(fp);
    }
    payload_file.close();
}

void DNA_encoder::decode(){
    // note: we leave this decode as an simple example. However, more decoder are temporarily not open sourced because we are furthur improving them for future work
    fstream decode_file;
    int prefix = rand();
    string payload_path = g_payload_path+to_string(prefix)+"decode";
    int num_of_GB = 1;
    unsigned long total_len = 0;
    payload_path += to_string(num_of_GB);
    decode_file.open(payload_path,ios::out);

    fstream input_file;

    unordered_map<string,string> decode_dirctionary;
    decode_dirctionary.emplace("TAC","000");
    decode_dirctionary.emplace("GTA","000");
    decode_dirctionary.emplace("CTG","000");
    decode_dirctionary.emplace("CAG","000");
    decode_dirctionary.emplace("ATA","000");

    decode_dirctionary.emplace("ATG","001");
    decode_dirctionary.emplace("CAT","001");
    decode_dirctionary.emplace("GAC","001");
    decode_dirctionary.emplace("GTC","001");
    decode_dirctionary.emplace("TAT","001");

    decode_dirctionary.emplace("ACT","010");
    decode_dirctionary.emplace("AGT","010");
    decode_dirctionary.emplace("TGC","010");
    decode_dirctionary.emplace("GCA","010");
    decode_dirctionary.emplace("CGC","010");

    decode_dirctionary.emplace("TCA","011");
    decode_dirctionary.emplace("TGA","011");
    decode_dirctionary.emplace("ACG","011");
    decode_dirctionary.emplace("CGT","011");
    decode_dirctionary.emplace("GCG","011");

    decode_dirctionary.emplace("ATC","100");
    decode_dirctionary.emplace("GAT","100");
    decode_dirctionary.emplace("TCG","100");
    decode_dirctionary.emplace("CGA","100");

    decode_dirctionary.emplace("ACA","101");
    decode_dirctionary.emplace("TGT","101");
    decode_dirctionary.emplace("CAC","101");
    decode_dirctionary.emplace("GTG","101");

    decode_dirctionary.emplace("TAG","110");
    decode_dirctionary.emplace("CTA","110");
    decode_dirctionary.emplace("AGC","110");
    decode_dirctionary.emplace("GCT","110");

    decode_dirctionary.emplace("AGA","111");
    decode_dirctionary.emplace("TCT","111");
    decode_dirctionary.emplace("CTC","111");
    decode_dirctionary.emplace("GAG","111");

    for(auto n:all_files_){
        input_file.open(n.c_str(),ios::in);
        string line;
        while (getline(input_file,line)) {
            if (line[0]=='>') continue;

            for (int i = 0; i < line.size()-3; i+=3) {
                string triplet = line.substr(i,3);
                decode_file<<decode_dirctionary.find(triplet)->second;
            }
        }
        input_file.close();
    }
    decode_file.close();
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

    list<string> candidates_row_0_4bit;
    candidates_row_0_4bit.emplace_back("TAC");
    candidates_row_0_4bit.emplace_back("GTA");
    four_bits_NT_triplets_candidates_.emplace(0,candidates_row_0_4bit);
    list<string> candidates_row_1_4bit;
    candidates_row_1_4bit.emplace_back("CTG");
    candidates_row_1_4bit.emplace_back("CAG");
    candidates_row_1_4bit.emplace_back("ATA");
    four_bits_NT_triplets_candidates_.emplace(1,candidates_row_1_4bit);

    list<string> candidates_row_1;
    candidates_row_1.emplace_back("ATG");
    candidates_row_1.emplace_back("CAT");
    candidates_row_1.emplace_back("GAC");
    candidates_row_1.emplace_back("GTC");
    candidates_row_1.emplace_back("TAT");
    three_bits_NT_triplets_candidates_.emplace(1,candidates_row_1);

    list<string> candidates_row_2_4bit;
    candidates_row_2_4bit.emplace_back("ATG");
    candidates_row_2_4bit.emplace_back("CAT");
    four_bits_NT_triplets_candidates_.emplace(2,candidates_row_2_4bit);
    list<string> candidates_row_3_4bit;
    candidates_row_3_4bit.emplace_back("GAC");
    candidates_row_3_4bit.emplace_back("GTC");
    candidates_row_3_4bit.emplace_back("TAT");
    four_bits_NT_triplets_candidates_.emplace(3,candidates_row_3_4bit);

    list<string> candidates_row_0_2bit;
    for (auto n:candidates_row_0) candidates_row_0_2bit.push_back(n);
    for (auto n:candidates_row_1) candidates_row_0_2bit.push_back(n);
    two_bits_NT_triplets_candidates_.emplace(0,candidates_row_0_2bit);

    list<string> candidates_row_2;
    candidates_row_2.emplace_back("ACT");
    candidates_row_2.emplace_back("AGT");
    candidates_row_2.emplace_back("TGC");
    candidates_row_2.emplace_back("GCA");
    candidates_row_2.emplace_back("CGC");
    three_bits_NT_triplets_candidates_.emplace(2,candidates_row_2);

    list<string> candidates_row_4_4bit;
    candidates_row_4_4bit.emplace_back("ACT");
    candidates_row_4_4bit.emplace_back("AGT");
    candidates_row_4_4bit.emplace_back("CGC");
    four_bits_NT_triplets_candidates_.emplace(4,candidates_row_4_4bit);
    list<string> candidates_row_5_4bit;
    candidates_row_5_4bit.emplace_back("TGC");
    candidates_row_5_4bit.emplace_back("GCA");
    four_bits_NT_triplets_candidates_.emplace(5,candidates_row_5_4bit);

    list<string> candidates_row_3;
    candidates_row_3.emplace_back("TCA");
    candidates_row_3.emplace_back("TGA");
    candidates_row_3.emplace_back("ACG");
    candidates_row_3.emplace_back("CGT");
    candidates_row_3.emplace_back("GCG");
    three_bits_NT_triplets_candidates_.emplace(3,candidates_row_3);

    list<string> candidates_row_6_4bit;
    candidates_row_6_4bit.emplace_back("TCA");
    candidates_row_6_4bit.emplace_back("TGA");
    candidates_row_6_4bit.emplace_back("GCG");
    four_bits_NT_triplets_candidates_.emplace(6,candidates_row_6_4bit);
    list<string> candidates_row_7_4bit;
    candidates_row_7_4bit.emplace_back("ACG");
    candidates_row_7_4bit.emplace_back("CGT");
    four_bits_NT_triplets_candidates_.emplace(7,candidates_row_7_4bit);

    list<string> candidates_row_1_2bit;
    for (auto n:candidates_row_2) candidates_row_1_2bit.push_back(n);
    for (auto n:candidates_row_3) candidates_row_1_2bit.push_back(n);
    two_bits_NT_triplets_candidates_.emplace(1,candidates_row_1_2bit);

    list<string> candidates_row_4;
    candidates_row_4.emplace_back("ATC");
    candidates_row_4.emplace_back("GAT");
    candidates_row_4.emplace_back("TCG");
    candidates_row_4.emplace_back("CGA");
    three_bits_NT_triplets_candidates_.emplace(4,candidates_row_4);

    list<string> candidates_row_8_4bit;
    candidates_row_8_4bit.emplace_back("ATC");
    candidates_row_8_4bit.emplace_back("GAT");
    four_bits_NT_triplets_candidates_.emplace(8,candidates_row_8_4bit);
    list<string> candidates_row_9_4bit;
    candidates_row_9_4bit.emplace_back("TCG");
    candidates_row_9_4bit.emplace_back("CGA");
    four_bits_NT_triplets_candidates_.emplace(9,candidates_row_9_4bit);

    list<string> candidates_row_5;
    candidates_row_5.emplace_back("ACA");
    candidates_row_5.emplace_back("TGT");
    candidates_row_5.emplace_back("CAC");
    candidates_row_5.emplace_back("GTG");
    three_bits_NT_triplets_candidates_.emplace(5,candidates_row_5);

    list<string> candidates_row_10_4bit;
    candidates_row_10_4bit.emplace_back("ACA");
    candidates_row_10_4bit.emplace_back("TGT");
    four_bits_NT_triplets_candidates_.emplace(10,candidates_row_10_4bit);
    list<string> candidates_row_11_4bit;
    candidates_row_11_4bit.emplace_back("CAC");
    candidates_row_11_4bit.emplace_back("GTG");
    four_bits_NT_triplets_candidates_.emplace(11,candidates_row_11_4bit);

    list<string> candidates_row_2_2bit;
    for (auto n:candidates_row_4) candidates_row_2_2bit.push_back(n);
    for (auto n:candidates_row_5) candidates_row_2_2bit.push_back(n);
    two_bits_NT_triplets_candidates_.emplace(2,candidates_row_2_2bit);

    list<string> candidates_row_6;
    candidates_row_6.emplace_back("TAG");
    candidates_row_6.emplace_back("CTA");
    candidates_row_6.emplace_back("AGC");
    candidates_row_6.emplace_back("GCT");
    three_bits_NT_triplets_candidates_.emplace(6,candidates_row_6);

    list<string> candidates_row_12_4bit;
    candidates_row_12_4bit.emplace_back("TAG");
    candidates_row_12_4bit.emplace_back("CTA");
    four_bits_NT_triplets_candidates_.emplace(12,candidates_row_12_4bit);
    list<string> candidates_row_13_4bit;
    candidates_row_13_4bit.emplace_back("AGC");
    candidates_row_13_4bit.emplace_back("GCT");
    four_bits_NT_triplets_candidates_.emplace(13,candidates_row_13_4bit);

    list<string> candidates_row_7;
    candidates_row_7.emplace_back("AGA");
    candidates_row_7.emplace_back("TCT");
    candidates_row_7.emplace_back("CTC");
    candidates_row_7.emplace_back("GAG");
    three_bits_NT_triplets_candidates_.emplace(7,candidates_row_7);

    list<string> candidates_row_14_4bit;
    candidates_row_14_4bit.emplace_back("AGA");
    candidates_row_14_4bit.emplace_back("TCT");
    four_bits_NT_triplets_candidates_.emplace(14,candidates_row_14_4bit);
    list<string> candidates_row_15_4bit;
    candidates_row_15_4bit.emplace_back("CTC");
    candidates_row_15_4bit.emplace_back("GAG");
    four_bits_NT_triplets_candidates_.emplace(15,candidates_row_15_4bit);

    list<string> candidates_row_3_2bit;
    for (auto n:candidates_row_6) candidates_row_3_2bit.push_back(n);
    for (auto n:candidates_row_7) candidates_row_3_2bit.push_back(n);
    two_bits_NT_triplets_candidates_.emplace(3,candidates_row_3_2bit);

    //based on normal 3bit candidate row, further add candidates
    // candidates_row_0.emplace_back("ATA"); already has
    candidates_row_0.emplace_back("TAT");
    candidates_row_0.emplace_back("CGC");
    candidates_row_0.emplace_back("GCG");
    complementary_first_triplets_candidates_.emplace(0,candidates_row_0);

    candidates_row_1.pop_back(); // delete TAT
    candidates_row_1.emplace_back("AAA");
    candidates_row_1.emplace_back("TTT");
    candidates_row_1.emplace_back("GGG");
    candidates_row_1.emplace_back("CCC");
    complementary_first_triplets_candidates_.emplace(1,candidates_row_1);

    candidates_row_2.pop_back(); // delete CGC
    candidates_row_2.emplace_back("AAT");
    candidates_row_2.emplace_back("ATT");
    candidates_row_2.emplace_back("CCG");
    candidates_row_2.emplace_back("CGG");
    complementary_first_triplets_candidates_.emplace(2,candidates_row_2);

    candidates_row_3.pop_back(); // delete GCG
    candidates_row_3.emplace_back("AAG");
    candidates_row_3.emplace_back("CTT");
    candidates_row_3.emplace_back("CCT");
    candidates_row_3.emplace_back("AGG");
    complementary_first_triplets_candidates_.emplace(3,candidates_row_3);

    candidates_row_4.emplace_back("AAC");
    candidates_row_4.emplace_back("GTT");
    candidates_row_4.emplace_back("CCA");
    candidates_row_4.emplace_back("TGG");
    complementary_first_triplets_candidates_.emplace(4,candidates_row_4);

    candidates_row_5.emplace_back("ACC");
    candidates_row_5.emplace_back("GGT");
    candidates_row_5.emplace_back("GAA");
    candidates_row_5.emplace_back("TTC");
    complementary_first_triplets_candidates_.emplace(5,candidates_row_5);

    candidates_row_6.emplace_back("TTA");
    candidates_row_6.emplace_back("TAA");
    candidates_row_6.emplace_back("GGC");
    candidates_row_6.emplace_back("GCC");
    complementary_first_triplets_candidates_.emplace(6,candidates_row_6);

    candidates_row_7.emplace_back("TTG");
    candidates_row_7.emplace_back("CAA");
    candidates_row_7.emplace_back("TCC");
    candidates_row_7.emplace_back("GGA");
    complementary_first_triplets_candidates_.emplace(7,candidates_row_7);

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


    //following are triplets that contain homo_2/3
    //ATA TAT CGC GCG has been insert already
    complementary_triplets_table_.emplace("AAA","TTT");
    complementary_triplets_table_.emplace("TTT","AAA");
    complementary_triplets_table_.emplace("GGG","CCC");
    complementary_triplets_table_.emplace("CCC","GGG");
    complementary_triplets_table_.emplace("AAT","ATT");
    complementary_triplets_table_.emplace("ATT","AAT");
    complementary_triplets_table_.emplace("CCG","CGG");
    complementary_triplets_table_.emplace("CGG","CCG");
    complementary_triplets_table_.emplace("AAG","CTT");
    complementary_triplets_table_.emplace("CTT","AAG");
    complementary_triplets_table_.emplace("CCT","AGG");
    complementary_triplets_table_.emplace("AGG","CCT");
    complementary_triplets_table_.emplace("AAC","GTT");
    complementary_triplets_table_.emplace("GTT","AAC");
    complementary_triplets_table_.emplace("CCA","TGG");
    complementary_triplets_table_.emplace("TGG","CCA");
    complementary_triplets_table_.emplace("ACC","GGT");
    complementary_triplets_table_.emplace("GGT","ACC");
    complementary_triplets_table_.emplace("GAA","TTC");
    complementary_triplets_table_.emplace("TTC","GAA");
    complementary_triplets_table_.emplace("TTA","TAA");
    complementary_triplets_table_.emplace("TAA","TTA");
    complementary_triplets_table_.emplace("GGC","GCC");
    complementary_triplets_table_.emplace("GCC","GGC");
    complementary_triplets_table_.emplace("TTG","CAA");
    complementary_triplets_table_.emplace("CAA","TTG");
    complementary_triplets_table_.emplace("TCC","GGA");
    complementary_triplets_table_.emplace("GGA","TCC");


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
    initial_FEC_table();
    init_heuristic_encoding();
    //initilize RS table
    init_GF47_table();
    //init_GF25_table();

    if (g_program==1){
        encoding_stranding();
    } else if (g_program==5){
        decode();
    }
}


// Fountain code
vector<float> DNA_encoder::ideal_distribution(int n) {
    const float EPS = 1e-5;
    assert(n > 0);

    vector<float> possibilities;
    possibilities.push_back(0.0);
    possibilities.push_back(1.0 / n);
    for (int i = 2; i <= n; i++) {
        possibilities.push_back(1.0 / (i * (i-1)));
    }
    float sum = 0.0;
    for (auto it = possibilities.begin(); it != possibilities.end(); it++) {
        sum += *it;
    }
    assert(sum <= 1.0 + EPS);
    assert(sum + EPS >= 1.0);

    return possibilities;
}

vector<float> DNA_encoder::robust_distribution(int n) {
    const float EPS = 1e-5;
    const float ROBUST_FAILURE_PROBABILITY = 0.01;

    int m = n / 2 + 1;
    float r = (float)n / m;

    vector<float> extra_proba;
    extra_proba.push_back(0);
    for (int i = 1; i < m; i++) {
        extra_proba.push_back(1.0/(i*m));
    }

    extra_proba.push_back(log(r / ROBUST_FAILURE_PROBABILITY)/m);
    for (int i = m + 1; i <= n; i++) {
        extra_proba.push_back(0.0);
    }
    
    vector<float> possibilities;
    vector<float> ideal_probabilities = ideal_distribution(n);
    float sum = 0.0;
    for (int i = 0; i <= n; i++) {
        possibilities.push_back(ideal_probabilities[i] + extra_proba[i]);
        sum += possibilities[i];
    }
    assert(abs(sum) >= EPS);
    for (int i = 0; i <= n; i++) {
        possibilities[i] /= sum;
    }
    sum = 0.0;
    for (int i = 0; i <= n; i++) {
        sum += possibilities[i];
    }
    assert(sum + EPS >= 1.0);
    assert(sum <= 1.0 + EPS);

    return possibilities;
}


vector<string> DNA_encoder::Fountain_encoding(string digital_data, float redundancy_ratio = 0.07) {
    if (digital_data.size() % 50 != 0) {
        string append(50 - digital_data.size() % 50, '\0');
        digital_data += append;
    }

    vector< bitset<400> > input;
    for (int i = 0; i < digital_data.size(); i += 50) {
        string ori = digital_data.substr(i, 50);
        string binary = "";
        for (int i = 0; i < ori.size(); i++) {
            int num = ori[i];
            string cur_binary(8, '\0');
            for(int j = 7; j >= 0; j--) {
                if (num & 1) cur_binary[j] = '1';
                else cur_binary[j] = '0';
                num >>= 1;
            }
            binary += cur_binary;
        }
        assert(binary.size() == 400);
        bitset<400> segment(binary);
        input.push_back(segment);
    }
    assert(input.size() * 50 == digital_data.size());

    int n = input.size();
    int n_out = (int)(n * (1 + redundancy_ratio)+0.5);
    vector<float> p = robust_distribution(n);
    vector<float> p_sum;
    p_sum.resize(n+1);
    p_sum[0] = 0.0;
    for(int i = 1; i <= n; i++) {
        p_sum[i] = p_sum[i-1] + p[i];
    }
    const float EPS = 1e-5;
    assert(p_sum[n] <= 1.0 + EPS);
    assert(p_sum[n] + EPS >= 1.0);

    srand (time(NULL));
    vector<string> result;
    while(result.size() < n_out) {
        float rand_num = (float)rand() / RAND_MAX;
        assert(rand_num <= 1.0 + EPS);
        int n_segment = lower_bound(p_sum.begin(), p_sum.end(), rand_num) - p_sum.begin();
        assert(n_segment >= 0);
        assert(n_segment <= n);

        bitset<400> cur_bits;
        for (int i = 0; i < n_segment; i++) {
            int choice = rand() % n;
            cur_bits ^= input[choice];
        }

        string cur_string(200, '\0');
        int n_homo = 0;
        bool flag_no_homo = true;
        int n_GC = 0;
        for (int i = 0, j = 0; i < 400 && j < 200; i += 2, j++) {
            int convert_id = cur_bits[i] * 2 + cur_bits[i+1];
            assert(convert_id >= 0);
            assert(convert_id < 4);
            switch(convert_id) {
                case 0: cur_string[j] = 'A';
                break;
                case 1: cur_string[j] = 'G'; n_GC++;
                break;
                case 2: cur_string[j] = 'C'; n_GC++;
                break;
                case 3: cur_string[j] = 'T';
            }

            if (j > 0 && cur_string[j] == cur_string[j-1]) {
                n_homo++;
                if(n_homo >= 3) {
                    flag_no_homo = false;
                    break;
                }
            } else {
                n_homo = 0;
            }
        }

        if(flag_no_homo && n_GC >= 200 * 0.45 && n_GC <= 200 * 0.55) {
            result.push_back(cur_string);
        }
    }

    return result;
}