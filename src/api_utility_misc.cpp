/**
 api_utility_misc.cpp
 Purpose: Miscellaneous logic.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */
 
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "malan_types.h"

//[[Rcpp::export]]
int pop_size(Rcpp::XPtr<Population> population) {
  return population->get_population_size();
}


//' Get all individuals in population
//' 
//' @param population Population
//'
//' @export
// [[Rcpp::export]]
Rcpp::ListOf< Rcpp::XPtr<Individual> > get_individuals(Rcpp::XPtr<Population> population) {     
  std::unordered_map<int, Individual*>* pop = population->get_population();
  int n = pop->size();
  Rcpp::List individuals(n);
  int i = 0;
  
  for (auto dest : *pop) {
    Rcpp::XPtr<Individual> indv_xptr(dest.second, RCPP_XPTR_2ND_ARG);
    individuals[i] = indv_xptr;
    
    if (i >= n) {
      Rcpp::stop("i > n");
    }
    
    i+= 1;
  }  

  return individuals;
}

//' Meiotic distribution
//' 
//' Get the distribution of number of meioses from `individual` 
//' to all individuals in `individual`'s pedigree.
//' Note the `generation_upper_bound_in_result` parameter.
//' 
//' @param individual Individual to calculate all meiotic distances from
//' @param generation_upper_bound_in_result Limit on distribution; -1 means no limit. 
//' 0 is the final generation. 1 second last generation etc.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix meioses_generation_distribution(Rcpp::XPtr<Individual> individual, 
                                                    int generation_upper_bound_in_result = -1) {  
  Individual* i = individual;
  
  if (!(i->pedigree_is_set())) {
    Rcpp::stop("Pedigree not yet set");
  }
  
  Pedigree* ped = i->get_pedigree();
  std::vector<Individual*>* family = ped->get_all_individuals();
  std::map<int, std::map<int, int> > tab;
  
  for (auto dest : *family) {    
    int generation = dest->get_generation();
    
    if (generation_upper_bound_in_result != -1 && generation > generation_upper_bound_in_result) {
      continue;
    }
    
    int dist = i->meiosis_dist_tree(dest);

    (tab[generation])[dist] += 1;    
  }
  
  int row = 0;
  for (auto const& x1 : tab) {
    row += x1.second.size();
    /*
    for (auto const& x2 : x1.second) {
      ++row;
    }
    */
  }
  Rcpp::IntegerMatrix res(row, 3);
  colnames(res) = Rcpp::CharacterVector::create("generation", "meioses", "count");
  row = 0;
  for (auto const& x1 : tab) {
    for (auto const& x2 : x1.second) {
      res(row, 0) = x1.first;
      res(row, 1) = x2.first;
      res(row, 2) = x2.second;
      ++row;    
    }
  }
  
  return res;
}




//' Size of population
//' 
//' Get the size of the population.
//' Note the `generation_upper_bound_in_result` parameter.
//' 
//' @param population Population to get size of
//' @param generation_upper_bound_in_result Limit on generation to include in count; -1 means no limit. 
//' 0 only include the final generation. 1 only second last generation etc.
//' 
//' @export
// [[Rcpp::export]]
int population_size_generation(Rcpp::XPtr<Population> population, int generation_upper_bound_in_result = -1) {  
  std::unordered_map<int, Individual*>* pop = population->get_population();
  
  int size = 0;
  
  for (auto dest : *pop) {    
    int generation = dest.second->get_generation();
    
    if (generation_upper_bound_in_result != -1 && generation > generation_upper_bound_in_result) {
      continue;
    }
    
    ++size;
  }
  
  return size;
}

//' Size of pedigree
//' 
//' Get the size of the pedigree.
//' Note the `generation_upper_bound_in_result` parameter.
//' 
//' @param pedigree Pedigree to get size of
//' @param generation_upper_bound_in_result Limit on generation to include in count; -1 means no limit. 
//' 0 only include the final generation. 1 only second last generation etc.
//' 
//' @export
// [[Rcpp::export]]
int pedigree_size_generation(Rcpp::XPtr<Pedigree> pedigree, int generation_upper_bound_in_result = -1) {  
  std::vector<Individual*>* family = pedigree->get_all_individuals();
  
  int size = 0;
  
  for (auto dest : *family) {    
    int generation = dest->get_generation();
    
    if (generation_upper_bound_in_result != -1 && generation > generation_upper_bound_in_result) {
      continue;
    }
    
    ++size;
  }
  
  return size;
}





