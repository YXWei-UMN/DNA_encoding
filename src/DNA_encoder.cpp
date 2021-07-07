//
// Created by eason on 2/22/21.
//
#include <math.h>
#include "DNA_encoder.h"
#include "tool.h"
#include "transformation.h"
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


//TODO their ECC code


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
            if (total_len>=1024*1024*1024){
                total_len = 0;
                payload_file.close();
                num_of_GB += 1;
                payload_path = g_payload_path+"payload"+to_string(num_of_GB)+".txt";
                payload_file.open(payload_path,ios::out);
            }
            //cout<<total_len<<endl;
            size_t len = fread(buf, 1, sizeof(buf), fp);
            total_len += len;
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
            else if(g_FEC_encoding){
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
                    if (g_if_mapping){
                        string permutated_strand = mapping(strand);
                        payload_file<<permutated_strand<<endl;
                    } else if(g_swap_granularity>0){
                        string permutated_strand = swap(strand);;
                        payload_file<<permutated_strand<<endl;
                    } else
                        payload_file<<strand<<endl;

                    digital_data.erase(0, 38);
                }
                //nt_sequence=FEC_encoding(digital_data);
                
            }
            //else if(g_xxx_code){}
            else
                cout<<"no encoding scheme"<<endl;

            /*write encoded strand to payload file
              >payload x
              ATCGATCG...*/
            
            // saved for future use
            // while (nt_sequence.size()>0){
            //     string strand = nt_sequence.substr(0, 200);

            //     // add ECC code
            //     // char strand_char[256];
            //     // strcpy(strand_char, strand.c_str());

            //     //string permutated_strand = swap(strand);
            //     payload_file<<">payload"<<strand_num++<<endl;

            //     // execute transformation: mapping/swap/...
            //     if (g_if_mapping){
            //         string permutated_strand = mapping(strand);
            //         payload_file<<permutated_strand<<endl;
            //     } else if(g_swap_granularity>0){
            //         string permutated_strand = swap(strand);;
            //         payload_file<<permutated_strand<<endl;
            //     } else
            //         payload_file<<strand<<endl;

            //     nt_sequence.erase(0, 200);
            // }
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
    string payload_path = g_payload_path+"payload";
    int num_of_files = 0;
    // payload_path += to_string(num_of_files);
    // payload_path += ".txt";
    //payload_file.open(payload_path,ios::out);

    
    //payload_file<<">payload0"<<endl;
    //create chunking buffer and related structure
    uint8_t buf[1024*1024];
    //go over all files to chunking and encoding
    FILE *fp;
    string nt_sequence;
    for(auto n:all_files_){
        num_of_files += 1;
        payload_path = g_payload_path+"payload"+to_string(num_of_files)+".txt";
        payload_file.open(payload_path,ios::out);
        payload_file<<">payload0"<<endl;

        fp = fopen(n.c_str(), "r");
        if (fp==NULL) {fputs ("File open error",stderr); exit (1);}
        while ( !feof(fp) ) {
            size_t len = fread(buf, 1, sizeof(buf), fp);
            uint8_t *ptr = &buf[0];

            string digital_data ((char*)ptr,len);

            if(g_base3_rotate_encoding)
                nt_sequence=base3_rotate_encoding(digital_data);
            else if (g_FEC_encoding)
                nt_sequence=FEC_encoding(digital_data);
            else
                cout<<"no encoding scheme"<<endl;
            payload_file<<nt_sequence;
            
        }
        fclose(fp);
        payload_file.close();
    }
    
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

    initial_FEC_table();


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