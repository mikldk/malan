/**
 api_simulate_helper.cpp
 Purpose: Helper for functions related to simulating populations.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <progress.hpp>

#include "malan_types.h"
#include "api_simulate.h"

using namespace Rcpp;

// Create a new father and add him to population. 
// Used in sample_geneology()
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
  List& last_k_generations_individuals) {  
  
  Individual* father = new Individual(*individual_id, generation);
  (*individual_id) = (*individual_id) + 1;
  
  fathers_generation[father_i] = father;
  (*population_map)[father->get_pid()] = father;
  
  if (verbose_result) {
    individual_pids_tmp_vec[father_i] = father->get_pid();
  }
  
  (*new_founders_left) = (*new_founders_left) + 1;

  if (generation <= individuals_generations_return) {
    Rcpp::XPtr<Individual> father_xptr(father, RCPP_XPTR_2ND_ARG);
    last_k_generations_individuals.push_back(father_xptr);
  }  
}

// Create a new father and add him to population. 
// Used in sample_geneology_varying_size()
void create_father_update_simulation_state_varying_size(
  int father_i, 
  int* individual_id, 
  int generation, 
  int individuals_generations_return,
  std::vector<Individual*>& fathers_generation, 
  std::unordered_map<int, Individual*>* population_map, 
  int* new_founders_left,
  List& last_k_generations_individuals) {  
  
  Individual* father = new Individual(*individual_id, generation);
  (*individual_id) = (*individual_id) + 1;
  
  fathers_generation[father_i] = father;
  (*population_map)[father->get_pid()] = father;

  (*new_founders_left) = (*new_founders_left) + 1;

  if (generation <= individuals_generations_return) {
    //Rcpp::Rcout << "create_father_update_simulation_state_varying_size: generation = " << generation << "; individuals_generations_return = " << individuals_generations_return << std::endl;
    
    Rcpp::XPtr<Individual> father_xptr(father, RCPP_XPTR_2ND_ARG);
    last_k_generations_individuals.push_back(father_xptr);
  }  
}

