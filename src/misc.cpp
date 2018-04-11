
#include "misc.h"

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
// add two numbers together in log space. One number (but not both) is allowed to be -inf.
double log_sum(double logA, double logB) {
    if (logA-logB > 100) {
        return logA;
    } else if (logB-logA > 100) {
        return logB;
    }
    double output = (logA<logB) ? logB + log(1+exp(logA-logB)) : logA + log(1+exp(logB-logA));
    return output;
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
