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

//TODO other three code

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

//RS_encoding
string DNA_encoder::ReedSolomon_encoding(string digital_data) {
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
    //cout << "rs_result" << rs_result << endl;

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


//rotate encoding without huffman tree compression
string DNA_encoder::base3_rotate_encoding(string digital_data) {
    string result;
    for (std::size_t i = 0; i < digital_data.size(); i++)
    {
        //cover from binary to decimal
        bitset<8> bits(digital_data.c_str()[i]);
        int decimal = bits.to_ulong();
        //covert to ternary & rotate encoding
        for (int i=0; i<6; i++){
            last_bit=rotating_encoding_table_[decimal%3][last_bit];
            
            result+=last_bit;
            decimal/=3;
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
    string payload_path = g_payload_path+"payload";
    int num_of_GB = 1;
    unsigned long total_len = 0;
    payload_path += to_string(num_of_GB);
    payload_path += ".txt";
    payload_file.open(payload_path,ios::out);

    long int strand_num=0;
    //create reading buffer
    uint8_t buf[1024*1024];
    //go over all files to chunking and encoding
    FILE *fp;
    string nt_sequence;

    for(auto n:all_files_){
        fp = fopen(n.c_str(), "r");
        if (fp==NULL) {fputs ("File open error",stderr); exit (1);}
        while ( !feof(fp) ) {

            // assign each GB data in different files
            // 4mins for a 280MB file
            if (total_len>=(long)5*1024*1024*1024){
                total_len = 0;
                payload_file.close();
                num_of_GB += 1;
                payload_path = g_payload_path+"payload"+to_string(num_of_GB)+".txt";
                payload_file.open(payload_path,ios::out);
            }
            //cout<<total_len<<endl;
            size_t len = fread(buf, 1, sizeof(buf), fp);
            // total_len += len;
            uint8_t *ptr = &buf[0];

            string digital_data ((char*)ptr,len);


            // rotate code
            if(g_encoding_scheme==1){
                string strand=base3_rotate_encoding(digital_data);
                payload_file<<">payload"<<strand_num++<<endl;
                // execute transformation: mapping/swap/...
                payload_file<<strand<<endl;
                total_len += strand.length();

                    
            }
            else if(g_encoding_scheme==0){
                string strand=direct_encoding(digital_data);
                payload_file<<">payload"<<strand_num++<<endl;
                // execute transformation: mapping/swap/...
                payload_file<<strand<<endl;
                total_len += strand.length();

                    
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
                    total_len += strand.length();
                    if (g_if_ECC){
                        digital_data.erase(0, 38);
                    }
                    else {
                        digital_data.erase(0, 40);
                    }

                }
                //nt_sequence=FEC_encoding(digital_data);
            }
            //else if(g_xxx_code){}
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

    //create chunking buffer and related structure
    uint8_t buf[1024*1024];
    //go over all files to chunking and encoding
    FILE *fp;
    string nt_sequence;
    int num_of_file = 0;
    for(auto n:all_files_){
        if (total_len>=(long)5*1024*1024*1024){
                total_len = 0;
                payload_file.close();
                num_of_payload_file += 1;
                payload_path = g_payload_path+"payload"+to_string(num_of_payload_file)+".txt";
                payload_file.open(payload_path,ios::out);
		num_of_file = 0;
        }

        fp = fopen(n.c_str(), "r");
        if (fp==NULL) {fputs ("File open error",stderr); exit (1);}

        payload_file<<">payload"<<num_of_file<<endl;
	num_of_file++;

        while ( !feof(fp) ) {
            size_t len = fread(buf, 1, sizeof(buf), fp);
            // total_len += len;
            uint8_t *ptr = &buf[0];

            string digital_data ((char*)ptr,len);

            if(g_encoding_scheme==1)
                nt_sequence=base3_rotate_encoding(digital_data);
            else if (g_encoding_scheme==2)
                nt_sequence=FEC_encoding(digital_data);
            else if (g_encoding_scheme==3)
                nt_sequence=ReedSolomon_encoding(digital_data);
            else
                cout<<"no encoding scheme"<<endl;
            payload_file<<nt_sequence;
            total_len += nt_sequence.length();
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


void DNA_encoder::init_RS_table() {
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

    //initilize RS table
    init_RS_table();


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