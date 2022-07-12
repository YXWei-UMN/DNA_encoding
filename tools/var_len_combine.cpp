#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cstdlib>

using namespace std;

const int MINN = 100;
const int MAXN = 200;

int last = 200;
int min_blind_spot = 200;
int min_uncovered = 200;

int len[4] = {100, 110, 132, 144};

int calc_blind_spot() {
  vector<int> com(len, len+4);
  unordered_set<int> s(len, len+4);

  for (int i = 0; i < 5; i++) {
    for (int x: com) {
      for (int y: com) {
        s.insert(x+y);
      }
    }

    com.clear();
    for (auto it = s.begin(); it != s.end(); it++) {
      if (*it <= 450) {
        com.push_back(*it);
      } 
    }
    sort(com.begin(), com.end());
  }

  for (int i = 200; i <= 400; i++) {
    bool found = false;
    for (int j = i - 12; j <= i; j++) {
      if (s.find(j)!=s.end()){
        found = true;
        break;
      }
    }
    if (!found) return -1;

    found = false;
    for (int j = i+1; j <= i + 12; j++) {
        if (s.find(j)!=s.end()){
        found = true;
        break;
      }
    }
    if (!found) return -1;
  }

  int cur_uncovered = com[0]-1;
  int contiguous_start = com[0];
  
  for (int i = 1; i < com.size(); i++) {
    assert(com[i-1] < com[i]);

    if (com[i] - com[i-1] > 12) {
      cur_uncovered += com[i] - com[i-1] -12;
      if (cur_uncovered >= min_uncovered) return -1;
      contiguous_start = com[i];
      if (contiguous_start > 200) return -1;
    } else {
      int contiguous_length = com[i] - contiguous_start + 1;
      if (contiguous_length >= 200) {
        min_uncovered = min(min_uncovered, cur_uncovered);
        if (min_uncovered < last) {
          last = min_uncovered;
          cout << "[uncovered=" << cur_uncovered << "] ";
          // for (int i = 0; i < 4; i++) cout << com[i] << " "; cout << endl;
          for (auto i = 0; i < com.size() && i < 20; i++) cout << com[i] << " "; cout << endl;
        }
        return 1;
      }
    }
  }


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
    if (i % 20 == 0) {
      cout << "i = " << i << "[" << min_blind_spot << "," << min_uncovered << "]\n";
    }
    len[0] = i;
    dfs(1);
  }

  cout << "min_blind_spot = " << min_blind_spot << endl;
  cout << "min_uncovered = " << min_uncovered << endl;
  return 0;
}
