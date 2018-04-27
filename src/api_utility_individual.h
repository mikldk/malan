/**
 api_utility_individual.h
 Purpose: Header file for logic related to individuals.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */

#ifndef MALAN_UTILITY_INDV_H
#define MALAN_UTILITY_INDV_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "malan_types.h"

Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid);
int get_pid(Rcpp::XPtr<Individual> individual);
void print_individual(Rcpp::XPtr<Individual> individual);
int get_generation(Rcpp::XPtr<Individual> individual);
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual);
Rcpp::IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, Rcpp::IntegerVector pids);
Rcpp::List get_family_info(Rcpp::XPtr<Individual> individual);
int count_brothers(Rcpp::XPtr<Individual> individual);
int brothers_matching(Rcpp::XPtr<Individual> individual);
bool father_matches(Rcpp::XPtr<Individual> individual);
bool grandfather_matches(Rcpp::XPtr<Individual> individual);
int count_uncles(Rcpp::XPtr<Individual> individual);
void unset_pedigree();

#endif
