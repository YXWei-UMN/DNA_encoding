#include <math.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <ctime>
#include <cstring>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <list>
#include <sstream>
#include <time.h>
#include <chrono>
#include <deque>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <bitset>
#include <utility>
#include <math.h>
#include <inttypes.h>
#include <set>
#include <unordered_set>
#include <queue>
#include <cstdlib>
#include <ctime>
using namespace std;


int main(){
    string path = "/home/wang9467/DNA/DNA_encoding/encoding_files/test.txt";
    string g_payload_path = "/home/wang9467/DNA/DNA_encoding/payload.txt";

    fstream payload_file;
    payload_file.open(g_payload_path,ios::out);

    uint8_t buf[1024*1024];
    //go over all files to chunking and encoding
    FILE *fp;
    fp = fopen(path.c_str(), "r");
    while ( !feof(fp) ) {
        size_t len = fread(buf, 1, sizeof(buf), fp);
        uint8_t *ptr = &buf[0];

        string digital_data ((char*)ptr,len);

        payload_file<<digital_data<<endl;
    }
    fclose(fp);
    payload_file.close();

}