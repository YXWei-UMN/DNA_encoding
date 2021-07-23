//
// Created by eason on 5/25/21.
//

#ifndef DNA_ENCODING_TOOL_H
#define DNA_ENCODING_TOOL_H

#include <iostream>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <bitset>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>

using namespace std;

bool isDir(string dir);


vector<string> listFiles(string baseDir, bool recursive);

#endif //DNA_ENCODING_TOOL_H
