// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// generate_scaffolds_cpp
Rcpp::List generate_scaffolds_cpp(Rcpp::List args);
RcppExport SEXP _rmaverick_generate_scaffolds_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(generate_scaffolds_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// run_mcmc_cpp
Rcpp::List run_mcmc_cpp(Rcpp::List args);
RcppExport SEXP _rmaverick_run_mcmc_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(run_mcmc_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// GTI_posterior_K_sim_cpp
Rcpp::List GTI_posterior_K_sim_cpp(Rcpp::List args);
RcppExport SEXP _rmaverick_GTI_posterior_K_sim_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(GTI_posterior_K_sim_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// GTI_integrated_K_sim_cpp
Rcpp::List GTI_integrated_K_sim_cpp(Rcpp::List args);
RcppExport SEXP _rmaverick_GTI_integrated_K_sim_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(GTI_integrated_K_sim_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// log_sum
double log_sum(double logA, double logB);
RcppExport SEXP _rmaverick_log_sum(SEXP logASEXP, SEXP logBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type logA(logASEXP);
    Rcpp::traits::input_parameter< double >::type logB(logBSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum(logA, logB));
    return rcpp_result_gen;
END_RCPP
}
// call_hungarian_cpp
Rcpp::List call_hungarian_cpp(Rcpp::List args);
RcppExport SEXP _rmaverick_call_hungarian_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(call_hungarian_cpp(args));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rmaverick_generate_scaffolds_cpp", (DL_FUNC) &_rmaverick_generate_scaffolds_cpp, 1},
    {"_rmaverick_run_mcmc_cpp", (DL_FUNC) &_rmaverick_run_mcmc_cpp, 1},
    {"_rmaverick_GTI_posterior_K_sim_cpp", (DL_FUNC) &_rmaverick_GTI_posterior_K_sim_cpp, 1},
    {"_rmaverick_GTI_integrated_K_sim_cpp", (DL_FUNC) &_rmaverick_GTI_integrated_K_sim_cpp, 1},
    {"_rmaverick_log_sum", (DL_FUNC) &_rmaverick_log_sum, 2},
    {"_rmaverick_call_hungarian_cpp", (DL_FUNC) &_rmaverick_call_hungarian_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rmaverick(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
