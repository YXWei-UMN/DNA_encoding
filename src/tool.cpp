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
            if (dirp->d_name[0] != '.') {
                if (isDir(baseDir + dirp->d_name) == true && recursive == true) {
                    //all_files_.push_back(baseDir + dirp->d_name);
                    vector<string> tmp =listFiles(baseDir + dirp->d_name + "/", true);
                    for(auto n:tmp){
                        all_files.push_back(n);
                    }
                } else {
                    all_files.push_back(baseDir + dirp->d_name);
                }
            }
        }
        closedir(dp);
    }
    return all_files;
}

bool CollisionPositionCMP(const Collision &c1,const Collision &c2)
{
    StrandID strand_id1 = get<0>(c1);
    unsigned int start1 = get<1>(c1);
    unsigned int end1 = get<2>(c1);
    PrimerID primer_id1 = get<3>(c1);
	unsigned int cut1 = get_collisiton_cut_point(c1);

    StrandID strand_id2 = get<0>(c2);
    unsigned int start2 = get<1>(c2);
    unsigned int end2 = get<2>(c2);
    PrimerID primer_id2 = get<3>(c2);
	unsigned int cut2 = get_collisiton_cut_point(c2);
    
    if (strand_id1 != strand_id2) return strand_id1 < strand_id2;
    if (cut1 != cut2) return cut1 < cut2;
    if (start1 != start2) return start1 < start2;
    if (end1 != end2) return end1 < end2;
    return primer_id1 < primer_id2;
}

unsigned int get_collisiton_cut_point(const Collision &c) {
    unsigned int start = get<1>(c);
    unsigned int end = get<2>(c);

    unsigned int cut = start + (end - start) / 2;
    cut = (int)(cut * 0.1 + 0.5);

    return cut * 10;
}
