#include <iostream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cstdlib>

using namespace std;

const int MINN = 100;
const int MAXN = 200;
const int LOOP_WINDOW = 400;

int min_blind_spot = 450;
int min_uncovered = 500;

const int n = 6;

int len[n];
int cnt = 0;


int calc_blind_spot() {
  vector<int> com(len, len+n);
  unordered_set<int> s(len, len+n);

  for (int i = 0; i < 6; i++) {// generate combinations 6 times, because 150*6 = 900 > 450
    for (int x: com) {
      for (int y: com) {
        s.insert(x+y);
      }
    }

    com.clear();
    for (auto it = s.begin(); it != s.end(); it++) {
      if (*it <= min_blind_spot+LOOP_WINDOW + 50) {
        com.push_back(*it);
      } 
    }
    sort(com.begin(), com.end());
  }

  for (int i = min_blind_spot; i <= min_blind_spot+LOOP_WINDOW; i++) {
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
      cur_uncovered += com[i] - com[i-1];
      contiguous_start = com[i];
      if (contiguous_start > min_blind_spot) return -1;
    } else {
      int contiguous_length = com[i] - contiguous_start + 1;
      if (contiguous_length >= LOOP_WINDOW) {
        if (contiguous_start <= min_blind_spot && cur_uncovered <= 299){
          min_blind_spot = min(min_blind_spot, contiguous_start);
          min_uncovered = min(min_uncovered, cur_uncovered);
          cout << "[" <<  (++cnt) << "]\t[blind_spot, uncovered] = " << "[" << contiguous_start << "," << cur_uncovered << "] " << n << "-choices = ";
          for (int i = 0; i < n; i++) cout << com[i] << " "; cout << endl;
          // for (auto i = 0; i < com.size() && i < 20; i++) cout << com[i] << " "; cout << endl;
        }
        return 1;
      }
    }
  }


  return 1;
}

void dfs(int step = 0) {
  if (step == n-1) {
    calc_blind_spot();
    // exit(0);
    return;
  }

  for (int i = len[step-1] + 1; i + n < MAXN; i++) {
    len[step] = i;
    dfs(step + 1);
  }
}

int main() {
  len[n-1] = 200;
  for (int i = 150; i + n < MAXN; i++) {
    // if (i % 20 == 0) {
    //   cout << "i = " << i << "[" << min_blind_spot << "," << min_uncovered << "]\n";
    // }
    len[0] = i;
    
    dfs(1);
  }

  cout << "min_blind_spot = " << min_blind_spot << endl;
  cout << "min_uncovered = " << min_uncovered << endl;
  return 0;
}
