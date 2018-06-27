
#pragma once

#include <random>

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1();

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(double a=0, double b=1.0);

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(double p);

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean=0, double sd=1);

//------------------------------------------------
// density of univariate normal distribution
double dnorm1(double x, double mean=0, double sd=1, bool log_on=true);

//------------------------------------------------
// draw from univariate normal distribution and reflect to interval (a,b)
double rnorm1_interval(double mean, double sd, double a, double b);

//------------------------------------------------
// sample single value from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double p_sum=1.0);

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal 
// probability. Works on positive or negative values of a or b, and works 
// irrespective of which of a or b is larger.
int sample2(int a, int b);

//------------------------------------------------
// sample a given number of values from a vector without replacement (templated
// for different data types). Note, this function re-arranges the original
// vector (passed in by reference), and the result is stored in the first n
// elements.
template<class TYPE>
void sample3(std::vector<TYPE> &x, int n) {
  int N = x.size();
  for (int i=0; i<n; i++) {
    int y = sample2(i,N-1);
    TYPE z = x[y];
    x[y] = x[i];
    x[i] = z;
  }
}

//------------------------------------------------
// draw from gamma(shape,rate) distribution
double rgamma1(double shape, double rate);

//------------------------------------------------
// draw from beta(shape1,shape2) distribution
double rbeta1(double shape1, double shape2);

//------------------------------------------------
// probability density of beta(shape1,shape2) distribution
double dbeta1(double x, double shape1, double shape2, bool log_on=true);

//------------------------------------------------
// draw from dirichlet distribution using vector of shape parameters
std::vector<double> rdirichlet1(std::vector<double> &shape_vec);

//------------------------------------------------
// draw from dirichlet distribution using bespoke inputs
void rdirichlet2(std::vector<double> &x, std::vector<int> &alpha, double beta, double gamma);

//------------------------------------------------
// draw from Poisson(rate) distribution
int rpois1(double rate);

//------------------------------------------------
// probability mass of Poisson(rate) distribution
double dpois1(int n, double rate, bool log_on=true);

//------------------------------------------------
// draw from negative binomial distribution with mean lambda and variance
// gamma*lambda (gamma must be >1)
int rnbinom1(double lambda, double gamma);

//------------------------------------------------
// probability mass of negative binomial distribution with mean lambda and
// variance gamma*lambda (gamma must be >1)
double dnbinom1(int n, double lambda, double gamma, bool log_on=true);

//------------------------------------------------
// draw number of unique groups under Chinese restaurant process
int rCRPgroups(int n, double alpha);
