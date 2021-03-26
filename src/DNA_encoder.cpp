//
// Created by eason on 2/22/21.
//
#include <math.h>
#include "DNA_encoder.h"

bool DNA_encoder::isDir(string dir)
{
    struct stat fileInfo;
    stat(dir.c_str(), &fileInfo);
    if (S_ISDIR(fileInfo.st_mode)) {
        return true;
    } else {
        return false;
    }
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

// Function to convert a decimal
// number to a ternary number
string convertToTernary(int N)
{
    // Base case
    if (N == 0)
        return "0";

    // Finding the remainder
    // when N is divided by 3
    int x = N % 3;
    N /= 3;

    // Recursive function to
    // call the function for
    // the integer division
    // of the value N/3
    string result = convertToTernary(N);

    return result+to_string(x);
}

// Function to convert the decimal to ternary
string convert(int Decimal)
{
    // If the number is greater
    // than 0, compute the
    // ternary representation
    // of the number
    if (Decimal != 0) {
        return convertToTernary(Decimal);
    }
    else
        return "000000";
}

string DNA_encoder::base3_rotate_encoding(string digital_data) {
    string result;
    for (std::size_t i = 0; i < digital_data.size(); i++)
    {
        bitset<8> bits(digital_data.c_str()[i]);
        int base = 1;
        int num = 0;
        for (int j = 7; j >=0; --j) {
            num += base*bits[j];
            base*=2;
        }

        //covert to ternary number
        string ternary_num = convert(num);
        //in case the termary number has less than 6 bits
        // 6 bits ternary cover 8 bits binary
        /*if (ternary_num.size()<6){
            for (int j = ternary_num.size(); j < 6; j++){
                ternary_num="0"+ternary_num;
            }
        }*/

        for(int j = 0; j<ternary_num.size(); j++){
            int bit = stoi(to_string(ternary_num[j]))-48;
            unordered_map<string,string> table = rotating_encoding_table_[bit];
            last_bit=table[last_bit];
            result+=last_bit;
        }
    }

    return result;
}


void DNA_encoder::chunking_encode(){
    // create payload file to store encoded strands
    fstream payload_file;
    payload_file.open(g_payload_path,ios::out);


    //create chunking buffer and related structure
    uint8_t buf[1024*1024];
    struct rabin_t *hash;
    hash = rabin_init();
    long long int chunk_num=g_payload_size*1024*1024/g_chunk_size;
    //go over all files to chunking and encoding
    FILE *fp;
    for(auto n:all_files_){
        fp = fopen(n.c_str(), "r");
        if (fp==NULL) {fputs ("File open error",stderr); exit (1);}

        while ( !feof(fp) ) {
            size_t len = fread(buf, 1, sizeof(buf), fp);
            uint8_t *ptr = &buf[0];

            while (1) {
                // remaining is the length of cut chunk
                /*
                 * int remaining = rabin_next_chunk(hash, ptr, len);
                if (remaining < 0) {
                    break;
                }*/
                if (len<64) break;
                int remaining = 64;
                total_chunks_++;
                len -= remaining;
                ptr += remaining;
                if (g_if_dedupe){
                    if(chunk_index_.find(last_chunk.cut_fingerprint)!=chunk_index_.end())
                        continue;
                }
                bytes += last_chunk.length;
                distinct_chunks_++;
                string acc = ">payload";
                acc += to_string(distinct_chunks_);
                payload_file<<acc<<endl;
                string digital_data ((char*)ptr-remaining,remaining);
                if(g_base3_rotate_encoding)
                    payload_file<<base3_rotate_encoding(digital_data)<<endl;
                //if(other encoding scheme)
                //  payload_file<<other_encoding(ptr-remaining,remaining)<<endl;
            }

        }
        if (distinct_chunks_ >= chunk_num){
            break;
        }
        /*if (rabin_finalize(hash) != NULL) {
            if (g_if_dedupe){
                if(chunk_index_.find(last_chunk.cut_fingerprint)!=chunk_index_.end())
                    continue;
            }
            chunks++;
            bytes +=last_chunk.length;
            string acc = ">payload";
            acc += to_string(chunks);
            cout<<acc<<endl;
            payload_file<<acc<<endl;
            string digital_data ((char*)(&buf[0]+last_chunk.start),last_chunk.length);
            if(g_base3_rotate_encoding)
                payload_file<<base3_rotate_encoding(digital_data)<<endl;
        }*/
        fclose(fp);
    }




    unsigned int avg = 0;
    double dedupe_rate = 0;
    if (distinct_chunks_ > 0){
        avg = bytes / distinct_chunks_;
        dedupe_rate = (total_chunks_-distinct_chunks_)/total_chunks_;
    }

    cout<<distinct_chunks_<<" distinct chunks"<<endl;
    cout<<"average chunk size "<<avg<<endl;
    cout<<"dedupe rate "<<dedupe_rate<<endl;
    payload_file.close();
}

// no stranding
void DNA_encoder::encoding(){
    // create payload file to store encoded strands
    fstream payload_file;
    payload_file.open(g_payload_path,ios::out);
    long int strand_num=0;
    //create chunking buffer and related structure
    uint8_t buf[1024*1024];
    long long int chunk_num=g_payload_size*1024*1024/g_chunk_size;
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
                nt_sequence+=base3_rotate_encoding(digital_data);
            else
                cout<<"no encoding scheme"<<endl;

            while (nt_sequence.size()>=200){
                string strand = nt_sequence.substr(0, 200);
                payload_file<<">payload"<<strand_num++<<endl;
                payload_file<<strand<<endl;
                nt_sequence.erase(0, 200);
            }

        }
        fclose(fp);
    }
    payload_file.close();
}

void DNA_encoder::encoding_stranding(){
    // create payload file to store encoded strands
    fstream payload_file;
    payload_file.open(g_payload_path,ios::out);

    //create chunking buffer and related structure
    uint8_t buf[1024*1024];
    long long int chunk_num=g_payload_size*1024*1024/g_chunk_size;
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
                nt_sequence+=base3_rotate_encoding(digital_data);
            else
                cout<<"no encoding scheme"<<endl;

            int remain_pos = nt_sequence.length();
            while (remain_pos > g_strand_length) {
                payload_file<<">payload"<<to_string(distinct_chunks_++)<<endl;
                payload_file<<nt_sequence.substr(nt_sequence.length()-remain_pos,g_strand_length)<<endl;
                remain_pos-=g_strand_length;
            }
            // keep remain nt_sequence for next iteration
            nt_sequence=nt_sequence.substr(remain_pos);

            if (distinct_chunks_ >= chunk_num){
                break;
            }
            }
        fclose(fp);
        }

    cout<<distinct_chunks_<<" number of strands"<<endl;
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
DNA_encoder::DNA_encoder() {
    // record all file's path
    listFiles(g_data_path, true);

    //initial rotating_encoding_table
    initial_rotating_encoding_table();

    //chunk and/or encode
    if (g_if_chunk){
        //chunking_encode();
        //encoding_stranding();
        encoding();
    } else {
        cout<<"object encoding without chunking, not complete yet~"<<endl;
    }
}