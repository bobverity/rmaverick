
// TODO
// - check if can return edges_left by reference, rather than ret by value
// - augment left and right, ret object needed vector?

#pragma once

#include <random>

//------------------------------------------------
// the functions augment_left and augment_right work together to find an augmented
// path. They call each other, which normally could lead to an infinite
// recursion, however, this is avoided as eventually either an augmented path
// will be found or no more moves will be possible. Return full path, or -1 if
// no path found.
std::vector<int> augment_left(int i, std::vector< std::vector<double> > &m, std::vector<int> &edges_right, std::vector<int> &blocked_left, std::vector<int> &blocked_right);
std::vector<int> augment_right(int j, std::vector< std::vector<double> > &m, std::vector<int> &edges_right, std::vector<int> &blocked_left, std::vector<int> &blocked_right);

//------------------------------------------------
// carry out Hungarian algorithm to find best matching given cost matrix m. If
// no matching found then print warning and return simple increasing sequence
std::vector<int> hungarian(std::vector< std::vector<double> > &m, std::vector<int> &edges_left, std::vector<int> &edges_right, std::vector<int> &blocked_left, std::vector<int> &blocked_right, const int max_reps=int(1e6));

