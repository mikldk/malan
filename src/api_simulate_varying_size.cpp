/**
 api_simulate_varying_size.cpp
 Purpose: Logic to simulate population of varying size.
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


//' Simulate a geneology with varying population size.
//' 
//' This function simulates a geneology with varying population size specified
//' by a vector of population sizes, one for each generation. 
//' 
//' By the backwards simulating process of the Wright-Fisher model, 
//' individuals with no descendants in the end population are not simulated 
//' If for some reason additional full generations should be simulated, 
//' the number can be specified via the \code{extra_generations_full} parameter.
//' This can for example be useful if one wants to simulate the 
//' final 3 generations although some of these may not get (male) children.
//' 
//' Let \eqn{\alpha} be the parameter of a symmetric Dirichlet distribution 
//' specifying each man's probability to be the father of an arbitrary 
//' male in the next generation. When \eqn{\alpha = 5}, a man's relative probability 
//' to be the father has 95\% probability to lie between 0.32 and 2.05, compared with a 
//' constant 1 under the standard Wright-Fisher model and the standard deviation in 
//' the number of male offspring per man is 1.10 (standard Wright-Fisher = 1).
//' 
//' This symmetric Dirichlet distribution is implemented by drawing 
//' father (unscaled) probabilities from a Gamma distribution with 
//' parameters `gamma_parameter_shape` and `gamma_parameter_scale` 
//' that are then normalised to sum to 1. 
//' To obtain a symmetric Dirichlet distribution with parameter \eqn{\alpha}, 
//' the following must be used:
//' \eqn{`gamma_parameter_shape` = \alpha}
//' and 
//' \eqn{`gamma_parameter_scale` = 1/\alpha}.
//' 
//' @param population_sizes The size of the population at each generation, `g`. 
//'        `population_sizes[g]` is the population size at generation `g`.
//'        The length of population_sizes is the number of generations being simulated.
//' @param extra_generations_full Additional full generations to be simulated.
//' @param gamma_parameter_shape Parameter related to symmetric Dirichlet distribution for each man's probability to be father. Refer to details.
//' @param gamma_parameter_scale Parameter realted to symmetric Dirichlet distribution for each man's probability to be father. Refer to details.
//' @param enable_gamma_variance_extension Enable symmetric Dirichlet (and disable standard Wright-Fisher).
//' @param progress Show progress.
//' @param individuals_generations_return How many generations back to return (pointers to) individuals for.
//' 
//' @return A malan_simulation / list with the following entries:
//' \itemize{
//'   \item \code{population}. An external pointer to the population.
//'   \item \code{generations}. Generations actually simulated, mostly useful when parameter \code{generations = -1}.
//'   \item \code{founders}. Number of founders after the simulated \code{generations}.
//'   \item \code{growth_type}. Growth type model.
//'   \item \code{sdo_type}. Standard deviation in a man's number of male offspring. StandardWF or GammaVariation depending on \code{enable_gamma_variance_extension}.
//'   \item \code{end_generation_individuals}. Pointers to individuals in end generation.
//'   \item \code{individuals_generations}. Pointers to individuals in end generation in addition to the previous \code{individuals_generations_return}.
//' }
//'
//' @seealso [sample_geneology()].
//' 
//' @import Rcpp
//' @import RcppProgress
//' @import RcppArmadillo
//' @export
// [[Rcpp::export]]
List sample_geneology_varying_size(
  IntegerVector population_sizes,
  int extra_generations_full = 0,  
  double gamma_parameter_shape = 5.0, double gamma_parameter_scale = 1.0/5.0, 
  bool enable_gamma_variance_extension = false,
  bool progress = true, 
  int individuals_generations_return = 2) {
  
  // boolean chosen like this to obey NA's
  bool all_gt_1 = is_true(all(population_sizes >= 1));
  if (!all_gt_1) {
    Rcpp::stop("Please specify only population_sizes >= 1");
  }
  
  int generations = population_sizes.length();
  
  if (generations == 0) {
    Rcpp::stop("Please specify at least 1 generation (the vector population_sizes must have length >= 1)");
  }

  // FIXME: malan_individual on Xptr indvs???
  
  if (enable_gamma_variance_extension) {
    if (gamma_parameter_shape <= 0.0) {
      Rcpp::stop("gamma_parameter_shape must be > 0.0");
    }
    if (gamma_parameter_scale <= 0.0) {
      Rcpp::stop("gamma_parameter_scale must be > 0.0");
    }
  }
  
  Progress progress_bar(generations, progress);
  
  // pid's are garanteed to be unique
  std::unordered_map<int, Individual*>* population_map = new std::unordered_map<int, Individual*>(); 
  Population* population = new Population(population_map);
  Rcpp::XPtr<Population> population_xptr(population, RCPP_XPTR_2ND_ARG);
  population_xptr.attr("class") = CharacterVector::create("malan_population", "externalptr");
  
  
  int individual_id = 1;
  std::vector<Individual*> end_generation(population_sizes[generations-1]);
  List end_generation_individuals(population_sizes[generations-1]);
  List last_k_generations_individuals;

  for (size_t i = 0; i < population_sizes[generations-1]; ++i) {
    Individual* indv = new Individual(individual_id++, 0);
    end_generation[i] = indv;    
    (*population_map)[indv->get_pid()] = indv;
    
    Rcpp::XPtr<Individual> indv_xptr(indv, RCPP_XPTR_2ND_ARG);
    end_generation_individuals[i] = indv_xptr;
    
    if (individuals_generations_return >= 0) {
      last_k_generations_individuals.push_back(indv_xptr);
    }
    
    if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      stop("Aborted");
    }
  }
  
  if (progress) {
    progress_bar.increment();
  }
  
  // Next generation  
  std::vector<Individual*> children_generation(population_sizes[generations-1]);
  for (size_t i = 0; i < population_sizes[generations-1]; ++i) {
    children_generation[i] = end_generation[i];
  }
  std::vector<Individual*> fathers_generation;
  
  int founders_left = population_sizes[generations-1];
  
  // now, find out who the fathers to the children are
  for (size_t generation = 1; generation < generations; ++generation) {
    // Init ->
    int population_size = population_sizes[generations-(generation+1)];    
    int children_population_size = population_sizes[generations-generation];

    WFRandomFather wf_random_father(population_size);
    GammaVarianceRandomFather gamma_variance_father(population_size, gamma_parameter_shape, gamma_parameter_scale);  
    SimulateChooseFather* choose_father = &wf_random_father;
    if (enable_gamma_variance_extension) {
      choose_father = &gamma_variance_father;
    }
    
    fathers_generation.clear();
    fathers_generation.resize(population_size);
    // <- Init
    
    int new_founders_left = 0;

    // clear
    for (size_t i = 0; i < population_size; ++i) {
      fathers_generation[i] = nullptr;
    }
    
    choose_father->update_state_new_generation();
    
    // now, run through children to pick each child's father
    for (size_t i = 0; i < children_population_size; ++i) {
      if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
        stop("Aborted");
      }
    
      // if a child did not have children himself, forget his ancestors
      if (children_generation[i] == nullptr) {
        continue;
      }
      
      int father_i = choose_father->get_father_i();
      
      // if this is the father's first child, create the father
      if (fathers_generation[father_i] == nullptr) {
        create_father_update_simulation_state_varying_size(father_i, &individual_id, generation, 
              individuals_generations_return, fathers_generation, population_map, 
              &new_founders_left, last_k_generations_individuals);
      }
            
      fathers_generation[father_i]->add_child(children_generation[i]);
    }
    
    // create additional fathers (without children) if needed:
    if (generation <= extra_generations_full) {
      for (size_t father_i = 0; father_i < population_size; ++father_i) {
        if (fathers_generation[father_i] != nullptr) {
          continue;
        }        
        
        // create father, no children etc.
        create_father_update_simulation_state_varying_size(father_i, &individual_id, generation, 
              individuals_generations_return, fathers_generation, population_map, 
              &new_founders_left, last_k_generations_individuals);
      }      
    }

    children_generation.clear();
    children_generation.resize(population_size);
    for (size_t i = 0; i < population_size; ++i) {
      children_generation[i] = fathers_generation[i];
    }
    
    if (Progress::check_abort()) {
      stop("Aborted");
    }
    
    if (progress) {
      progress_bar.increment();
    }
    
    founders_left = new_founders_left;
  }
  
  List res;
  res["population"] = population_xptr;
  res["generations"] = generations;
  res["founders"] = founders_left;
  res["growth_type"] = "VaryingPopulationSize";
  res["sdo_type"] = (enable_gamma_variance_extension) ? "GammaVariation" : "StandardWF";
  res["end_generation_individuals"] = end_generation_individuals;
  res["individuals_generations"] = last_k_generations_individuals;

  res.attr("class") = CharacterVector::create("malan_simulation", "list");
  
  return res;
}

