/**
 api_simulate.h
 Purpose: Header file for population simulation.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>

#include "malan_types.h"

using namespace Rcpp;

void create_father_update_simulation_state(
  int father_i, 
  int* individual_id, 
  int generation, 
  int individuals_generations_return,
  std::vector<Individual*>& fathers_generation, 
  std::unordered_map<int, Individual*>* population_map, 
  IntegerVector& individual_pids_tmp_vec,
  bool verbose_result,
  int* new_founders_left,
  List& last_k_generations_individuals);
  
void create_father_update_simulation_state_varying_size(
  int father_i, 
  int* individual_id, 
  int generation, 
  int individuals_generations_return,
  std::vector<Individual*>& fathers_generation, 
  std::unordered_map<int, Individual*>* population_map, 
  int* new_founders_left,
  List& last_k_generations_individuals);
