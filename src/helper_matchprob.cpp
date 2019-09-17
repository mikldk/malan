/**
 helder_matchprob.cpp
 Purpose: C++ helper functions for match probability calculations.
 Details: C++ implementation.
  
 @author Mikkel Meyer Andersen
 */
 
#include <RcppArmadillo.h>
#include "malan_types.h"

// Heavily relies on it being a TREE, hence there is only one path connecting every pair of nodes
void meiotis_dist_all_internal(Individual* indv, 
                               int current_dist, 
                               std::unordered_map<int, int>& dists) {
  
  if (indv->dijkstra_was_visited()) {
    return;
  }

  dists[indv->get_pid()] = current_dist;
  indv->dijkstra_mark_visited();
  
  if (current_dist % 10 == 0) {
    Rcpp::checkUserInterrupt();
  }
  
  current_dist += 1;
 
  Individual* father = indv->get_father();
  if (father != nullptr) {  
    meiotis_dist_all_internal(father, current_dist, dists); 
  }
  
  std::vector<Individual*>* children = indv->get_children();
  for (auto child : *children) {
    meiotis_dist_all_internal(child, current_dist, dists); 
  }
}


//' Find meiotic dist to all other individuals in pedigree
//' 
//' @param individual Individual
//' 
//' @export
// [[Rcpp::export]]
std::unordered_map<int, int> meiotis_dist_all(Rcpp::XPtr<Individual> individual) {
  // Heavily relies on it being a TREE, hence there is only one path connecting every pair of nodes
  
  Individual* indv = individual;
    
  if (!(indv->pedigree_is_set())) {
   throw std::invalid_argument("individual's pedigree not set");
  }
  
  // Reset all individuals in pedigree
  std::vector<Individual*>* inds = indv->get_pedigree()->get_all_individuals();
  for (auto member : *inds) {
    member->dijkstra_reset();
  }
  
  
  // pid -> dist
  std::unordered_map<int, int> dists;
  
  meiotis_dist_all_internal(indv, 0, dists);
  
  return dists;
}

//' Find meiotic dist to all other individuals in pedigree
//' 
//' Return an efficient look-up container
//' 
//' @param individual Individual
//' 
//' @export
// [[Rcpp::export]]
std::unordered_map<int, int> meiotis_dist_all_lookup(Rcpp::XPtr<Individual> individual) {
  std::unordered_map<int, int> dist = meiotis_dist_all(individual);
  // FIXME
  return dist;
}

