/**
 api_igraph.cpp
 Purpose: Convert igraph to population
 Details: C++ implementation.
  
 @author Mikkel Meyer Andersen
 */
 
#include <RcppArmadillo.h>

#include "malan_types.h"
#include "helper_generation.h"

using namespace Rcpp;

//' Generate paternal brothers population
//' 
//' @param vertices vector of vertices
//' @param edges matrix with edges
//' 
//' @return An external pointer to the population.
// [[Rcpp::export]]
Rcpp::XPtr<Population> from_igraph_rcpp(Rcpp::IntegerVector vertices, 
                                        Rcpp::IntegerMatrix edges) {
  
  // pid's are garanteed to be unique
  std::unordered_map<int, Individual*>* population_map = 
    new std::unordered_map<int, Individual*>(); 
  
  Population* population = new Population(population_map);
  Rcpp::XPtr<Population> population_xptr(population, RCPP_XPTR_2ND_ARG);
  population_xptr.attr("class") = CharacterVector::create("malan_population", "externalptr");
  
  //////////
  
  for (int pid : vertices) {
    Individual* indv = new Individual(pid);
    (*population_map)[indv->get_pid()] = indv;
  }
  
  for (size_t i = 0; i < edges.nrow(); ++i) {
    int from_pid = edges(i, 0);
    int to_pid = edges(i, 1);
    
    Individual* from = (*population_map)[from_pid];
    Individual* to = (*population_map)[to_pid];
    
    from->add_child(to);
  }
  
  return population_xptr;
}




//' Infer generation numbers from pedigrees
//' 
//' @param peds Pedigrees infered by [build_pedigrees()]
//' 
//' @return Nothing
//' 
//' @export
// [[Rcpp::export]]
void infer_generations(Rcpp::XPtr< std::vector<Pedigree*> > peds) {
  std::vector<Pedigree*> peds_vec = *peds;
  
  for (Pedigree* p : peds_vec) {
    Individual* indv = p->get_root();
    
    // Founder:
    int end_generation_number = 0;
    update_generation(indv, 0, &end_generation_number, 1);
    // Revert numbering:
    update_generation(indv, end_generation_number - 1, &end_generation_number, -1);
  }
}
