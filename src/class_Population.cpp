/**
 class_Population.cpp
 Purpose: C++ class Population.
 Details: C++ implementation.
  
 @author Mikkel Meyer Andersen
 */

#include <RcppArmadillo.h>
//#include <Rcpp.h>

#include "malan_types.h"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

/*
==========================================
Individual
==========================================
*/
Population::Population(std::unordered_map<int, Individual*>* population) {
  m_population = population;  
}

Population::~Population() {
  std::unordered_map<int, Individual*> pop = *m_population;
  
  /* Remember that both Individual and Pedigree has 
   * destructors of their own to clean their own members.
   */
  
  for (auto it = pop.begin(); it != pop.end(); ++it) {
    if (it->second == nullptr) {
      continue;
    }
    
    Pedigree* ped = it->second->get_pedigree();    
    if (ped != nullptr) {
      delete ped;
    }
    
    delete (it->second);
  }
  
  delete m_population;
}

std::unordered_map<int, Individual*>* Population::get_population() const {
  return m_population;
}

Individual* Population::get_individual(int pid) const {
  std::unordered_map<int, Individual*>::const_iterator got = m_population->find(pid);
  
  if (got == m_population->end()) {
    Rcpp::Rcerr << "Individual with pid = " << pid << " not found!" << std::endl;
    Rcpp::stop("Individual not found");
  }
  
  return got->second;
}

int Population::get_population_size() const {
  return m_population->size();
}
