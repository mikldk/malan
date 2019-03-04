/**
 test_misc.cpp
 Purpose: Code used in testing package.
 Details: C++ implementation.
  
 @author Mikkel Meyer Andersen
 */
 
#include <RcppArmadillo.h>

#include "malan_types.h"

using namespace Rcpp;

//' Generate test population
//' 
//' @return An external pointer to the population.
// [[Rcpp::export]]
Rcpp::XPtr<Population> test_create_population() {
  
  // pid's are garanteed to be unique
  std::unordered_map<int, Individual*>* population_map = 
    new std::unordered_map<int, Individual*>(); 
  
  Population* population = new Population(population_map);
  Rcpp::XPtr<Population> population_xptr(population, RCPP_XPTR_2ND_ARG);
  population_xptr.attr("class") = CharacterVector::create("malan_population", "externalptr");
  
  //////////

  std::vector<Individual*> indvs;

  Individual* i1 = new Individual(1, 0); indvs.push_back(i1);
  Individual* i2 = new Individual(2, 0); indvs.push_back(i2);
  Individual* i3 = new Individual(3, 0); indvs.push_back(i3);
  Individual* i4 = new Individual(4, 0); indvs.push_back(i4);
  Individual* i5 = new Individual(5, 0); indvs.push_back(i5);
  
  Individual* i6 = new Individual(6, 1); indvs.push_back(i6);
  Individual* i7 = new Individual(7, 1); indvs.push_back(i7);
  Individual* i8 = new Individual(8, 1); indvs.push_back(i8);
  
  Individual* i9 = new Individual(9, 2); indvs.push_back(i9);
  Individual* i10 = new Individual(10, 2); indvs.push_back(i10);
  
  Individual* i11 = new Individual(11, 3); indvs.push_back(i11);
  
  Individual* i12 = new Individual(12, 3); indvs.push_back(i12);
  
  /*
   *     
   * G3           11              12
   *           /     \
   * G2      9        10
   *        /  \       |
   * G1    6    7      8
   *       |   /\     /\
   * G0    1   2 3   4 5
   */
  i11->add_child(i9);
  i11->add_child(i10);
  i9->add_child(i6);
  i9->add_child(i7);
  i10->add_child(i8);
  
  i6->add_child(i1);
  i7->add_child(i2);
  i7->add_child(i3);
  
  i8->add_child(i4);
  i8->add_child(i5);
  
  for (auto i : indvs) {
    (*population_map)[i->get_pid()] = i;
  }
  
  return population_xptr;
}

