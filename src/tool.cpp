#include "tool.h"

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

vector<string> listFiles(string baseDir, bool recursive)
{
    vector<string> all_files;
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(baseDir.c_str())) == NULL) {
        cout << "[ERROR: " << errno << " ] Couldn't open " << baseDir << "." << endl;
        return all_files;
    } else {
        while ((dirp = readdir(dp)) != NULL) {
            if (dirp->d_name != string(".") && dirp->d_name != string("..")) {
                if (isDir(baseDir + dirp->d_name) == true && recursive == true) {
                    //all_files_.push_back(baseDir + dirp->d_name);
                    vector<string> tmp =listFiles(baseDir + dirp->d_name + "/", true);
                    all_files.assign(tmp.begin(),tmp.end());
                } else {
                    all_files.push_back(baseDir + dirp->d_name);
                }
            }
        }
        closedir(dp);
    }
    return all_files;
}