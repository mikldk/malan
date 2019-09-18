/**
 helder_matchprob.cpp
 Purpose: C++ helper functions for match probability calculations.
 Details: C++ implementation.
  
 @author Mikkel Meyer Andersen
 */
 
#include <RcppArmadillo.h>
#include "malan_types.h"

// Heavily relies on it being a TREE, hence there is only one path connecting every pair of nodes
void meiotic_dist_all_internal_visit(Individual* indv, 
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
    meiotic_dist_all_internal_visit(father, current_dist, dists); 
  }
  
  std::vector<Individual*>* children = indv->get_children();
  for (auto child : *children) {
    meiotic_dist_all_internal_visit(child, current_dist, dists); 
  }
}

// Returns map: pid -> dist
std::unordered_map<int, int> meiotic_dist_all_internal(Rcpp::XPtr<Individual> individual) {
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
  meiotic_dist_all_internal_visit(indv, 0, dists);
  
  return dists;
}

//' Find meiotic dist to all other individuals in pedigree
//' 
//' @param individual Individual to get distances from
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix meiotic_dist_all(Rcpp::XPtr<Individual> individual) {
  // pid -> dist
  std::unordered_map<int, int> dists = meiotic_dist_all_internal(individual);
  
  int n = dists.size();
  int i = 0;
  
  Rcpp::IntegerMatrix out(n, 2);
  for (auto it = dists.begin(); it != dists.end(); ++it) {
    out(i, 0) = it->first;
    out(i, 1) = it->second;
    i += 1;
  }
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("pid", "dist");
  
  return out;
}

//' Find meiotic dist to all other individuals in pedigree
//' 
//' Return an efficient look-up container
//' 
//' @param individual Individual to get distances from
//' @param pids Only return distance to individuals with these pids
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector meiotic_dist_pids(Rcpp::XPtr<Individual> individual, 
                                      Rcpp::IntegerVector pids) {
  // pid -> dist
  std::unordered_map<int, int> dists = meiotic_dist_all_internal(individual);
  
  int n = pids.size();
  Rcpp::IntegerVector out(n);
  int idx = 0;
  
  std::unordered_map<int, int>::const_iterator got;
  
  for (int i = 0; i < pids.size(); ++i) {
    got = dists.find(pids[i]);
    
    if (got == dists.end()) {
      out[idx] = Rcpp::IntegerVector::get_na();
    } else {
      out[idx] = got->second;
    }
    
    idx += 1;
  }
  
  return out;
}


//' Find meiotic dist to all other individuals in pedigree
//' 
//' Return an efficient look-up container
//' 
//' @param individual Individual
//' @param individuals Only return distance to these individuals
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector meiotic_dist_individuals(Rcpp::XPtr<Individual> individual, 
                                             Rcpp::List individuals) {
  // pid -> dist
  std::unordered_map<int, int> dists = meiotic_dist_all_internal(individual);
  
  int n = individuals.size();
  Rcpp::IntegerVector out(n);
  int idx = 0;
  
  std::unordered_map<int, int>::const_iterator got;
  
  for (int i = 0; i < individuals.size(); ++i) {
    Rcpp::XPtr<Individual> i_xptr = Rcpp::as< Rcpp::XPtr<Individual> >(individuals[i]);
    
    got = dists.find(i_xptr->get_pid());
    
    if (got == dists.end()) {
      out[idx] = Rcpp::IntegerVector::get_na();
    } else {
      out[idx] = got->second;
    }
    
    idx += 1;
  }
  
  return out;
}

