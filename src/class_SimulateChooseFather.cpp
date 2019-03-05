/**
 class_SimulateChooseFather.cpp
 Purpose: C++ class SimulateChooseFather.
 Details: C++ implementation.
  
 @author Mikkel Meyer Andersen
 */
 
#include "malan_types.h"

#include <RcppArmadillo.h> // FIXME: Avoid Rcpp here? Only in api_* files?


/*****************************************
WFRandomFather
******************************************/
WFRandomFather::WFRandomFather(size_t population_size) {
  m_population_size = (double)population_size;
}

void WFRandomFather::update_state_new_generation() {  
  //Rcpp::Rcout << "WFRandomFather: update_state_new_generation NOOP" << std::endl;
}

int WFRandomFather::get_father_i() {
  //Rcpp::Rcout << "WFRandomFather: get_father_i" << std::endl;
  return R::runif(0, 1)*m_population_size;
}


/*****************************************
GammaVarianceRandomFather
******************************************/
GammaVarianceRandomFather::GammaVarianceRandomFather(size_t population_size, double gamma_parameter_shape, double gamma_parameter_scale) {
  m_population_size = population_size;
  m_gamma_parameter_shape = gamma_parameter_shape;
  m_gamma_parameter_scale = gamma_parameter_scale;
}

// modified from 
// https://github.com/RcppCore/RcppArmadillo/blob/master/inst/include/RcppArmadilloExtensions/sample.h
// ProbSampleReplace for size = 1
void GammaVarianceRandomFather::update_state_new_generation() {    
  //Rcpp::Rcout << "GammaVarianceRandomFather: update_state_new_generation" << std::endl;
  
  Rcpp::NumericVector fathers_prob_tmpl = Rcpp::rgamma(m_population_size, m_gamma_parameter_shape, m_gamma_parameter_scale);
  //Rcpp::Rcout << "mean[fathers_prob_tmpl] = " << Rcpp::mean(fathers_prob_tmpl) << ", var[fathers_prob_tmpl] = " << Rcpp::var(fathers_prob_tmpl) << std::endl;
  fathers_prob_tmpl = fathers_prob_tmpl / sum(fathers_prob_tmpl);    

  arma::vec fathers_prob(fathers_prob_tmpl.begin(), fathers_prob_tmpl.size(), false); // false means no copy
  arma::uvec fathers_prob_perm = arma::sort_index(fathers_prob, "descend"); //descending sort of index
  fathers_prob = arma::sort(fathers_prob, "descend");  // descending sort of prob
  fathers_prob = arma::cumsum(fathers_prob);
  
  m_fathers_prob_cum = fathers_prob;
  m_fathers_prob_perm = fathers_prob_perm;
}

int GammaVarianceRandomFather::get_father_i() {
  //Rcpp::Rcout << "GammaVarianceRandomFather: get_father_i" << std::endl;
  
  double rU = R::unif_rand(); // R's internal random number generation
  //double rU = R::runif(0, 1);

  int jj;
  size_t population_size_1 = m_population_size - 1;
  
  for (jj = 0; jj < population_size_1; ++jj) {
    if (rU <= m_fathers_prob_cum[jj]) {
      break;
    }
  }
  
  return m_fathers_prob_perm[jj];
}


