/**
 helder_Individual.cpp
 Purpose: C++ helper functions for Individual class.
 Details: C++ implementation.
  
 @author Mikkel Meyer Andersen
 */
 
#include <RcppArmadillo.h>
//#include <Rcpp.h>

#include "malan_types.h"

bool find_path_from_root_to_dest(Individual* root, std::vector<Individual*>& path, const Individual* dest) {
  if (root == NULL || root == nullptr) {
    return false;
  }
  
  int dest_pid = dest->get_pid();

  path.push_back(root);

  if (root->get_pid() == dest_pid) {
    return true;
  }

  std::vector<Individual*>* children = root->get_children();
  for (auto child : *children) {
    if (find_path_from_root_to_dest(child, path, dest)) {
      return true;
    }
  }

  path.pop_back();
  return false;
}

// Draw autosomal genetype
// 
// @param allele_dist Allele distribution (probabilities) -- gets normalised
// @param alleles Names of alleles
// @param theta Theta correction between 0 and 1 (both included)
// 
// @return Vector of length 2 with indices of alleles
std::vector<int> draw_autosomal_genotype(
    const std::vector<double>& allele_cumdist_theta,
    const int alleles_count) {
  
  std::vector<int> geno(2);
  geno[0] = -1;
  geno[1] = -1;

  double u = R::runif(0.0, 1.0);
  bool stop = false;

  int k = 0;
  for (int i = 0; i < alleles_count; ++i) {
    for (int j = 0; j <= i; ++j) {   
      if (u <= allele_cumdist_theta[k]) {
      
        // note that j <= i so get index0
        geno[0] = j;
        geno[1] = i;
        
        stop = true;
        break;
      }
      
      k++;
    }
    
    if (stop) {
      break;
    }
  }
  
  if (geno[0] == -1 || geno[1] == -1) {
    throw std::invalid_argument("geno not set!");
  }
  
  return geno;
}

