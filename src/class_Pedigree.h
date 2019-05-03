/**
 class_Pedigree.h
 Purpose: Header for C++ class Pedigree.
 Details: C++ header.
  
 @author Mikkel Meyer Andersen
 */

#include <RcppArmadillo.h>

#include "malan_types.h"

#include <vector>

class Pedigree {
private:
  int m_pedigree_id;
  std::vector<Individual*>* m_all_individuals = nullptr;
  std::vector< std::pair<Individual*, Individual*>* >* m_relations = nullptr; 
  Individual* m_root = nullptr;  
  
public:
  Pedigree(int id);
  ~Pedigree();
  int get_id() const;
  void add_member(Individual* i);
  void add_relation(Individual* lhs, Individual* rhs);
  std::vector<Individual*>* get_all_individuals() const;
  std::vector< std::pair<Individual*, Individual*>* >* get_relations() const;
  
  Individual* get_root();
  
  void populate_haplotypes(
      const int loci, 
      const std::vector<double>& mutation_rates,
      const Rcpp::Function& get_founder_hap,
      const double prob_two_step = 0.0,
      const double prob_genealogical_error = 0.0);
  
  void populate_haplotypes_custom_founders(
      const std::vector<double>& mutation_rates, 
      const Rcpp::Function& get_founder_hap,
      const double prob_two_step = 0.0,
      const double prob_genealogical_error = 0.0);
  
  void populate_haplotypes_ladder_bounded(
      const std::vector<double>& mutation_rates, 
      const std::vector<int>& ladder_min, 
      const std::vector<int>& ladder_max, 
      const Rcpp::Function& get_founder_hap,
      const double prob_two_step = 0.0,
      const double prob_genealogical_error = 0.0);
  
  void populate_autosomal(
    const std::vector< std::vector<double> >& allele_conditional_cumdists_theta,
    const std::vector<double>& allele_cumdist_theta,
    const int alleles_count,
    const double mutation_rate);
};

