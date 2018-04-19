
#include "Hungarian.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// the functions augment_left and augment_right work together to find an
// augmented path. They call each other, which normally could lead to an
// infinite recursion, but this is avoided as eventually either an augmented
// path will be found or no more moves will be possible. If no matching found
// then print warning and return simple increasing sequence
vector<int> augment_left(int i, vector< vector<double> > &m, vector<int> &edges_right, vector<int> &blocked_left, vector<int> &blocked_right) {
    
    blocked_left[i] = 1;
    vector<int> ret(1,-1);
    
    // search all unmatched edges
    for (int j=0; j<int(m.size()); j++) {
        if (m[i][j]==0 && edges_right[j]!=i && blocked_right[j]==0) {
            
            // if edge leads to augmented path then add current node to path and return
            ret = augment_right(j, m, edges_right, blocked_left, blocked_right);
            if (ret[0]>=0) {
                ret.push_back(i);
                return ret;
            }
        }
    }
    
    // if no more moves then return -1
    return ret;
}

vector<int> augment_right(int j, vector< vector<double> > &m, vector<int> &edges_right, vector<int> &blocked_left, vector<int> &blocked_right) {
    
    blocked_right[j] = 1;
    vector<int> ret(1);
    
    // if node j is unmatched then return j as start of augmented path
    if (edges_right[j]<0) {
        ret[0] = j;
        return ret;
    }
    
    // otherwise continue chain of augmenting
    ret = augment_left(edges_right[j], m, edges_right, blocked_left, blocked_right);
    if (ret[0]>=0) {
        ret.push_back(j);
    }
    return ret;
}

//------------------------------------------------
// carry out Hungarian algorithm to find best matching given cost matrix m. If
// no best matching found then return vector with first element -1 to indicate
// error.
vector<int> hungarian(vector< vector<double> > &m, vector<int> &edges_left, vector<int> &edges_right, vector<int> &blocked_left, vector<int> &blocked_right, const int max_reps) {
  
  // initialise objects
  int n = m.size();
  vector<double> min_col(n);
  
  // search for solution until max_reps reached
  for (int rep=0; rep<max_reps; rep++) {
    
    // reset search objects
    fill(edges_left.begin(), edges_left.end(), -1);
    fill(edges_right.begin(), edges_right.end(), -1);
    int number_assigned = 0;
    
    // subtract smallest element from all rows and columns
    min_col = m[0];
    for (int i=0; i<n; i++) {
      double min_row = min(m[i]);
      for (int j=0; j<n; j++) {
        m[i][j] -= min_row;
        if (m[i][j]<min_col[j]) {
          min_col[j] = m[i][j];
        }
      }
    }
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        m[i][j] -= min_col[j];
      }
    }
    
    // generate an initial matching
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (m[i][j]==0 && edges_right[j]<0) {
          edges_left[i] = j;
          edges_right[j] = i;
          number_assigned ++;
          break;
        }
      }
    }
    
    // if this matching is perfect then we are done
    if (number_assigned==n) {
      return edges_left;
    }
    
    // continue augmenting paths until no more possible
    bool continue_augmenting = true;
    while (continue_augmenting) {
      continue_augmenting = false;
      
      // search all unmatched nodes
      for (int i=0; i<n; i++) {
        if (edges_left[i]<0) {
          
          // attempt to find augmented path
          blocked_left = vector<int>(n);
          blocked_right = vector<int>(n);
          vector<int> path = augment_left(i, m, edges_right, blocked_left, blocked_right);
          
          // if successful then augment
          if (path[0]>=0) {
            continue_augmenting = true;
            number_assigned ++;
            for (int j=0; j<int(path.size()/2); j++) {
              edges_left[path[j*2+1]] = path[j*2];
              edges_right[path[j*2]] = path[j*2+1];
            }
            
            // if best matching found then finish
            if (number_assigned==n) {
              return edges_left;
            }
          }
        }
      }
    }
    
    // find minimum value in cost matrix, looking at all elements in which neither the row or the column is part of the minimum vertex cover
    double min_val = -log(double(0));
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (blocked_left[i]==1 && blocked_right[j]==0 && m[i][j]<min_val) {
          min_val = m[i][j];
        }
      }
    }
    
    // add or subtract this value from cost matrix as required
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (blocked_left[i]==1 && blocked_right[j]==0) {
          m[i][j] -= min_val;
        }
        if (blocked_left[i]==0 && blocked_right[j]==1) {
          m[i][j] += min_val;
        }
      }
    }
    
    // at this point we have a new cost matrix and can repeat the process from the top
    
  } // rep loop
  
  // if reached this point then not managed to find best matching
  print("Warning: Hungarian algorithm unable to find best matching");
  edges_left = seq_int(0, n-1);
  
  return edges_left;
}
