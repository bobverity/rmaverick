
#include "misc.h"
#include "Hungarian.h"

using namespace std;

//------------------------------------------------
// define very large/small numbers for catching overflow/underflow problems
// OVERFLO
// UNDERFLO
// DEFINED IN HEADER

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types)
// sum
// DEFINED IN HEADER

//------------------------------------------------
// mean of vector (templated for different data types)
// mean
// DEFINED IN HEADER

//------------------------------------------------
// min of vector (templated for different data types)
// min
// DEFINED IN HEADER

//------------------------------------------------
// max of vector (templated for different data types)
// max
// DEFINED IN HEADER

//------------------------------------------------
// push back multiple values to vector
// push_back_multiple
// DEFINED IN HEADER

//------------------------------------------------
// test whether value can be found in vector
// is_in_vec
// DEFINED IN HEADER

//------------------------------------------------
// return unique values in a contiguous sequence of integers. v_max gives the
// maximum possible value in v, which might be larger than the maximum actual
// value in v.
// Example1: v = {1 1 2 2 0 0} and v_max = 2
// Result: unique_int(v) = {1 2 0}
// Example2: v = {1 1 2 2 0 0} and v_max = 3
// Result: unique_int(v) = {1 2 0 3}
// Example3: v = {1 1 3 3 0 0} and v_max = 3
// Result: unique_int(v) = {1 3 0 2}
vector<int> unique_int(const vector<int> &v, int v_max) {
  
  vector<int> mask(v_max);
  vector<int> ret(v_max);
  int j = 0;
  for (int i=0; i<int(v.size()); i++) {
    if (mask[v[i]]==0) {
      ret[j] = v[i];
      mask[v[i]] = 1;
      j++;
    }
  }
  for (int k=0; k<v_max; k++) {
    if (mask[k]==0) {
      ret[j] = k;
      j++;
    }
  }
  return ret;
}

//------------------------------------------------
// return order of unique values in a contiguous sequence of integers. v_max 
// gives the maximum possible value in v, which might be larger than the maximum
// actual value in v.
// Example1: v = {1 1 2 2 0 0} and v_max = 2
// Result: order_unique_int(v) = {2 0 1}
// Example2: v = {1 1 2 2 0 0} and v_max = 3
// Result: order_unique_int(v) = {2 0 1 3}
// Example3: v = {1 1 3 3 0 0} and v_max = 3
// Result: order_unique_int(v) = {2 0 3 1}
vector<int> order_unique_int(const vector<int> &v, int v_max) {
  
  vector<int> ret(v_max, -1);
  int j = 0;
  for (int i=0; i<int(v.size()); i++) {
    if (ret[v[i]]<0) {
      ret[v[i]] = j;
      j++;
    }
  }
  for (int k=0; k<v_max; k++) {
    if (ret[k]<0) {
      ret[k] = j;
      j++;
    }
  }
  return ret;
}

//------------------------------------------------
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
// [[Rcpp::export]]
double log_sum(double logA, double logB) {
    if (logA-logB > 100) {
        return logA;
    } else if (logB-logA > 100) {
        return logB;
    }
    double ret = (logA<logB) ? logB + log(1+exp(logA-logB)) : logA + log(1+exp(logB-logA));
    return ret;
}

//------------------------------------------------
// helper function for printing a single value (templated for different data types)
// print
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
// print_vector
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
// print_matrix
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
// print_array
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of an Rcpp numeric vector
void rcpp_print_vector(Rcpp::NumericVector &x) {
  for (int i=0; i<x.length(); i++) {
    Rcpp::Rcout << x[i] << " ";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}
void rcpp_print_vector(Rcpp::IntegerVector &x) {
  for (int i=0; i<x.length(); i++) {
    Rcpp::Rcout << x[i] << " ";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of an Rcpp numeric matrix
void rcpp_print_matrix(Rcpp::NumericMatrix &x) {
  for (int i=0; i<x.nrow(); i++) {
    for (int j=0; j<x.ncol(); j++) {
      Rcpp::Rcout << x(i,j) << " ";
    }
    Rcpp::Rcout << "\n";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}
void rcpp_print_matrix(Rcpp::IntegerMatrix &x) {
  for (int i=0; i<x.nrow(); i++) {
    for (int j=0; j<x.ncol(); j++) {
      Rcpp::Rcout << x(i,j) << " ";
    }
    Rcpp::Rcout << "\n";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void print_stars(string title, int n) {
  Rcpp::Rcout << title << " ";
  for (int i=0; i<n; i++) {
    Rcpp::Rcout << "*";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
}

//------------------------------------------------
// print "foo", with option number e.g. "foo2"
void foo(int n) {
  if (n==0) {
    Rcpp::Rcout << "foo\n";
  } else {
    Rcpp::Rcout << "foo" << n << "\n";
  }
  R_FlushConsole();
}

//------------------------------------------------
// print "bar", with option number e.g. "bar2"
void bar(int n) {
  if (n==0) {
    Rcpp::Rcout << "bar\n";
  } else {
    Rcpp::Rcout << "bar" << n << "\n";
  }
  R_FlushConsole();
}

//------------------------------------------------
// print "foobar", with option number e.g. "foobar2"
void foobar(int n) {
  if (n==0) {
    Rcpp::Rcout << "foobar\n";
  } else {
    Rcpp::Rcout << "foobar" << n << "\n";
  }
  R_FlushConsole();
}

//------------------------------------------------
// analogue of R function seq() for integers
vector<int> seq_int(int from, int to, int by) {
  int n = floor((to-from)/double(by)) + 1;
  vector<int> ret(n,from);
  for (int i=1; i<n; i++) {
    from += by;
    ret[i] = from;
  }
  return ret;
}

//------------------------------------------------
// converts input from Rcpp::List format to bool format.
int rcpp_to_bool(SEXP x) {
  return Rcpp::as<bool>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to int format.
int rcpp_to_int(SEXP x) {
  return Rcpp::as<int>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to double format.
double rcpp_to_double(SEXP x) {
  return Rcpp::as<double>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to string format.
string rcpp_to_string(SEXP x) {
  return Rcpp::as<string>(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<int> format.
vector<int> rcpp_to_vector_int(SEXP x) {
  return Rcpp::as<vector<int> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<double> format.
vector<double> rcpp_to_vector_double(SEXP x) {
  return Rcpp::as<vector<double> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<string> format.
vector<string> rcpp_to_vector_string(SEXP x) {
  return Rcpp::as<vector<string> >(x);
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<bool>> format.
vector< vector<bool> > rcpp_to_mat_bool(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<bool> > x_mat(nrow);
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<bool> >(x[i]);
  }
  return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<int>> format.
vector< vector<int> > rcpp_to_mat_int(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<int> > x_mat(nrow);
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<int> >(x[i]);
  }
  return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
vector< vector<double> > rcpp_to_mat_double(Rcpp::List x) {
  int nrow = int(x.size());
  vector< vector<double> > x_mat(nrow);
  for (int i=0; i<nrow; i++) {
    x_mat[i] = Rcpp::as<vector<double> >(x[i]);
  }
  return x_mat;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
vector< vector< vector<int> > > rcpp_to_array_int(Rcpp::List x) {
  int n1 = int(x.size());
  vector< vector< vector<int> > > ret(n1);
  for (int i=0; i<n1; i++) {
    Rcpp::List x_i = x[i];
    int n2 = int(x_i.size());
    ret[i] = vector< vector<int> >(n2);
    for (int j=0; j<n2; j++) {
      ret[i][j] = Rcpp::as<vector<int> >(x_i[j]);
    }
  }
  return ret;
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
vector< vector< vector<double> > > rcpp_to_array_double(Rcpp::List x) {
  int n1 = int(x.size());
  vector< vector< vector<double> > > ret(n1);
  for (int i=0; i<n1; i++) {
    Rcpp::List x_i = x[i];
    int n2 = int(x_i.size());
    ret[i] = vector< vector<double> >(n2);
    for (int j=0; j<n2; j++) {
      ret[i][j] = Rcpp::as<vector<double> >(x_i[j]);
    }
  }
  return ret;
}
