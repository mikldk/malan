/**
 api_load_data.cpp
Purpose: Construct populations from data.

@author Mikkel Meyer Andersen
*/

#include "malan_types.h"
#include "helper_generation.h"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;

void update_generation(Individual* indv, 
                       int generation_number, 
                       int* end_generation_number, 
                       const int modifier);

//' Construct a population from data
//' 
//' Note that individuals loaded this way does not have information about generation.
//' 
//' @param pid ID of male
//' @param pid_dad ID of male's father, 0 if not known
//' @param progress Show progress.
//' @param error_on_pid_not_found Error if pid not found
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<Population> load_individuals(IntegerVector pid, 
                                        IntegerVector pid_dad, 
                                        bool progress = true,
                                        bool error_on_pid_not_found = true) {
  
  std::unordered_map<int, Individual*>* pop = new std::unordered_map<int, Individual*>(); // pid's are garanteed to be unique
  
  Population* population = new Population(pop);
  Rcpp::XPtr<Population> population_xptr(population, RCPP_XPTR_2ND_ARG_CLEANER);
  population_xptr.attr("class") = CharacterVector::create("malan_population_abort", "externalptr");
  
  int N = pid.size();
  
  if (pid_dad.size() != N) {
    stop("all vectors (pid, pid_dad) must have same length");
  }
  
  // N for building hashmap, and N for building relations (adding children)
  //Progress p(2*N, true);
  Progress p(2*N, progress);
  
  // Build hashmap
  for (int k = 0; k < N; ++k) {
    int i_pid = pid[k];
    
    Individual* i = new Individual(i_pid);
    
    (*pop)[i->get_pid()] = i;
    
    if (k % CHECK_ABORT_EVERY == 0 && Progress::check_abort() ) {
      return population_xptr;
    }
    
    if (progress) {
      p.increment();
    }
  }
  
  // Fill out parents
  for (int k = 0; k < N; ++k) {
    //if (lim_parents_children != -1 && k > lim_parents_children) {
    //  break;
    //}
    
    int i_pid = pid[k];
    int i_pid_dad = pid_dad[k];
    
    Individual* i = (*pop)[i_pid];
    
    if (i_pid_dad > 0) {
      std::unordered_map<int, Individual*>::const_iterator f = pop->find(i_pid_dad);    
      if (f == pop->end()) {
        std::ostringstream err;
        err << "NOT FOUND: pid_dad = " << i_pid_dad << " for pid = " << i_pid << " was not found as a pid itself!";
        
        if (error_on_pid_not_found) {
          stop(err.str());
        } else {
          warning(err.str());
        }
      } else {
        Individual* father = f->second;
        father->add_child(i);
      }
    }
    
    if (k % CHECK_ABORT_EVERY == 0 && Progress::check_abort() ) {      
      return population_xptr;
    }
    
    if (progress) {
      p.increment();
    }
  }
  
  // Now, population is constructed. 
  // But the generation is not yet set for the individuals.
  // This will be done now:
  // 1) Take all founders, i.e. individuals with no father, 
  //    and set generation to 0.
  // 2) Visit all children and set generation to generation+1.
  // 3) Normally, final generation is 0, so end by "revert" the 
  //    generation by offset.
  //
  // Not possible (or feasible) to take leaves as that would 
  // mean multiple visits.
  for (std::pair<int, Individual*> indv_pair : *pop) {
    Individual* indv = indv_pair.second;
    
    if (indv->get_father() != nullptr) {
      continue;
    }
    
    // Founder:
    int end_generation_number = 0;
    update_generation(indv, 0, &end_generation_number, 1);
    // Revert numbering:
    update_generation(indv, end_generation_number - 1, &end_generation_number, -1);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////
  
  population_xptr.attr("class") = CharacterVector::create("malan_population", "externalptr");
  
  return population_xptr;
}
