//
// Created by eason on 2/22/21.
//
#include <math.h>
#include "DNA_encoder.h"
#include "tool.h"
#include "transformation.h"

//TODO other three code
//TODO their ECC code


//rotate encoding without huffman tree compression
string DNA_encoder::base3_rotate_encoding(string digital_data) {
    string result;
    for (std::size_t i = 0; i < digital_data.size(); i++)
    {
        //cover from binary to decimal
        bitset<8> bits(digital_data.c_str()[i]);
        int base = 1;
        int num = 0;
        for (int j = 7; j >=0; --j) {
            num += base*bits[j];
            base*=2;
        }

        //covert to ternary number, defined in tool.h
        string ternary_num = convert(num);

        //in case the termary number has less than 6 bits
        // 6 bits ternary cover 8 bits binary
        if (ternary_num.size()<6){
            for (int j = ternary_num.size(); j < 6; j++){
                ternary_num="0"+ternary_num;
            }
        }

        // base 3 rotate encoding
        for(int j = 0; j<ternary_num.size(); j++){
            int bit = stoi(to_string(ternary_num[j]))-48;
            // rotating_encoding_table_ is predifined when initialize DNA_encoder class
            unordered_map<string,string> table = rotating_encoding_table_[bit];
            // first nt is define is predifined in DNA_encoder.h as string last_bit
            last_bit=table[last_bit];
            result+=last_bit;
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
    payload_file.open(g_payload_path,ios::out);

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
            size_t len = fread(buf, 1, sizeof(buf), fp);
            uint8_t *ptr = &buf[0];

            string digital_data ((char*)ptr,len);

            // randomize digital video_data (XOR)
            if (g_if_randomization){
                randomize_XOR(digital_data);
            }

            // rotate code
            if(g_base3_rotate_encoding)
                nt_sequence+=base3_rotate_encoding(digital_data);
            //else if(g_fountain_code){}
            //else if(g_xxx_code){}
            //else if(g_xxx_code){}
            else
                cout<<"no encoding scheme"<<endl;

            /*write encoded strand to payload file
              >payload x
              ATCGATCG...*/
            while (nt_sequence.size()>=200){
                string strand = nt_sequence.substr(0, 200);
                //string permutated_strand = swap(strand);
                payload_file<<">payload"<<strand_num++<<endl;

                // execute transformation: mapping/swap/...
                if (g_if_mapping){
                    string permutated_strand = mapping(strand);
                    payload_file<<permutated_strand<<endl;
                } else if(g_swap_granularity>0){
                    string permutated_strand = swap(strand);;
                    payload_file<<permutated_strand<<endl;
                } else
                    payload_file<<strand<<endl;

                nt_sequence.erase(0, 200);
            }

        }
        fclose(fp);
    }
    payload_file.close();
}

void DNA_encoder::encoding_file(){
    // create payload file to store encoded strands
    fstream payload_file;
    payload_file.open(g_payload_path,ios::out);

    long int strand_num=0;
    //create reading buffer
    uint8_t buf[1024*1024];
    //go over all files to chunking and encoding
    FILE *fp;
    string strand;
    for(auto n:all_files_){
        payload_file<<">payload"<<strand_num++<<endl;
        fp = fopen(n.c_str(), "r");
        if (fp==NULL) {fputs ("File open error",stderr); exit (1);}
        while ( !feof(fp) ) {
            size_t len = fread(buf, 1, sizeof(buf), fp);
            uint8_t *ptr = &buf[0];
            string digital_data ((char*)ptr,len);
            strand=base3_rotate_encoding(digital_data);
            payload_file<<strand;
        }
        payload_file<<endl;
        fclose(fp);
    }
    payload_file.close();
}

// no strand
// it's the original video_data stream feed for our algorithm (to cut & transformation)
void DNA_encoder::encoding_no_strand(){
    // create payload file to store encoded strands
    fstream payload_file;
    payload_file.open(g_payload_path,ios::out);
    payload_file<<">payload0"<<endl;
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
            uint8_t *ptr = &buf[0];

            string digital_data ((char*)ptr,len);

            if(g_base3_rotate_encoding)
                nt_sequence=base3_rotate_encoding(digital_data);
            else
                cout<<"no encoding scheme"<<endl;
            payload_file<<nt_sequence;
        }
        fclose(fp);
    }
    payload_file.close();
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

void DNA_encoder::listFiles(string baseDir, bool recursive)
{
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(baseDir.c_str())) == NULL) {
        cout << "[ERROR: " << errno << " ] Couldn't open " << baseDir << "." << endl;
        return;
    } else {
        while ((dirp = readdir(dp)) != NULL) {
            if (dirp->d_name != string(".") && dirp->d_name != string("..")) {
                if (isDir(baseDir + dirp->d_name) == true && recursive == true) {
                    //all_files_.push_back(baseDir + dirp->d_name);
                    listFiles(baseDir + dirp->d_name + "/", true);
                } else {
                    all_files_.push_back(baseDir + dirp->d_name);
                }
            }
        }
        closedir(dp);
    }
}

DNA_encoder::DNA_encoder() {
    // generate 267-length pseudo random. 200 ternary code is 266.6 binary code (6:8)
    srand(std::time(nullptr)); // use current time as seed for random generator
    for (int n=0; n < 267; ++n) {
        pseudo_random_sequence_+=std::rand()%2;
    }


    // record all file's path, use to read file & encode it
    listFiles(g_data_path, true);

    //initial rotating_encoding_table
    initial_rotating_encoding_table();

    //encode with strand;
    //encoding_stranding();
    //encode without strand
    //encoding_no_strand();
    encoding_file();
}