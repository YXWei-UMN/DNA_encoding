#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cstdlib>

using namespace std;

const int MINN = 100;
const int MAXN = 200;

int last = 450;
int min_blind_spot = 450;
int min_uncovered = 100000;


int len[4] = {0, 0, 0, 0};

int calc_blind_spot() {
  vector<int> com(len, len+4);
  unordered_set<int> s(len, len+4);

  // for (auto i = 0; i < 4; i++) cout << com[i] << " "; cout << endl;
  // for (auto it = s.begin(); it != s.end(); it++) cout << (*it) << " "; cout << endl;

  for (int i = 0; i < 6; i++) {//doing combination 5 times will get at least 100*5 == 500, which is larger than the know min_blind_spot 450
    for (int x: com) {
      for (int y: com) {
        s.insert(x+y);
      }
    }
    // for (auto it = s.begin(); it != s.end(); it++) cout << (*it) << " "; cout << endl;

    com.clear();
    for (auto it = s.begin(); it != s.end(); it++) {
      if (*it < 750) {
        com.push_back(*it);
      } 
    }
    sort(com.begin(), com.end());
  }



  for (int i = min_blind_spot; i <= min_blind_spot+ 200; i++) {
    bool found = false;
    for (int j = i - 12; j <= i + 12; j++) {
      if (s.find(j)!=s.end()){
        found = true;
        break;
      }
    }
    if (!found) return -1;
  }

  // cout << "size = " << com.size() << endl;
  // for (auto i = 0; i < com.size(); i++) cout << com[i] << " "; cout << endl;
  // cout << "min_blind_spot = " << min_blind_spot << endl;
  

  int cur_uncovered = com[0]-1;

  int contiguous_start = com[0];
  for (int i = 1; i < com.size() && i < min_blind_spot; i++) {
    assert(com[i-1] < com[i]);
    if (com[i] - com[i-1] > 24) {
      cur_uncovered += com[i] - com[i-1] -24;
      contiguous_start = com[i];
    } else {
      int contiguous_length = com[i] - contiguous_start + 1;
      if (contiguous_length >= 200) {//200??
        min_blind_spot = min(min_blind_spot, com[i]);
      }
    }
  }
  if (min_blind_spot < last) {
    cout << "min_blind_spot = " << min_blind_spot << endl;
    last = min_blind_spot;
    // cout << "size = " << com.size() << endl;
    // for (auto i = 0; i < com.size(); i++) cout << com[i] << " "; cout << endl;
  }
  min_uncovered = min(min_uncovered, cur_uncovered);

  return 1;
}

void dfs(int step = 0) {
  if (step == 4) {
    calc_blind_spot();
    // exit(0);
    return;
  }

  for (int i = len[step-1] + 1; i + 4 < MAXN; i++) {
    len[step] = i;
    dfs(step + 1);
  }
}

int main() {
  for (int i = 100; i + 4 < MAXN; i++) {
    len[0] = i;
    dfs(1);
    if (i % 10 == 0) {
      cout << "i = " << i << "[" << min_blind_spot << "," << min_uncovered << "]\n";
    }
  }
  cout << "min_blind_spot = " << min_blind_spot << endl;
  cout << "min_uncovered = " << min_uncovered << endl;
  return 0;
}
