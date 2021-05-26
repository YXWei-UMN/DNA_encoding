//
// Created by eason on 5/25/21.
//

#ifndef DNA_ENCODING_TOOL_H
#define DNA_ENCODING_TOOL_H
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

#endif //DNA_ENCODING_TOOL_H
