/**
 class_Pedigree.cpp
 Purpose: C++ class Pedigree.
 Details: C++ implementation.
  
 @author Mikkel Meyer Andersen
 */

#include "malan_types.h"

#include <RcppArmadillo.h> // FIXME: Avoid Rcpp here? Only in api_* files?

/*
==========================================
Pedigree
==========================================
*/

Pedigree::Pedigree(int id) {
  //Rcpp::Rcout << "Pedigree with id = " << id << " created" << std::endl;
  
  m_pedigree_id = id;
  
  m_all_individuals = new std::vector<Individual*>();
  m_relations = new std::vector< std::pair<Individual*, Individual*>* >();
}

Pedigree::~Pedigree() {
  //Rcpp::Rcout << "   CLEANUP PEDIGREE!\n";
  
  for (auto it = m_all_individuals->begin(); it != m_all_individuals->end(); ++it) {
    /*
    Only unsets individual's knowledge about pedigree.
    This pedigree is responsible for removing this individuals from its list of individuals and
    relations.
    */
    (*it)->unset_pedigree();
  }  
  delete m_all_individuals;
  
  for (auto it = m_relations->begin(); it != m_relations->end(); ++it) {
    delete *it;
  }
  
  delete m_relations;
}

int Pedigree::get_id() const {
  return m_pedigree_id;
}

void Pedigree::add_member(Individual* i) {
  m_all_individuals->push_back(i);
}

void Pedigree::add_relation(Individual* lhs, Individual* rhs) {
  //std::pair<Individual*, Individual*>* pair = new std::pair(lhs, rhs);
  std::pair<Individual*, Individual*>* pair = new std::pair<Individual*, Individual*>(lhs, rhs);
  m_relations->push_back(pair);
}

std::vector<Individual*>* Pedigree::get_all_individuals() const {
  return m_all_individuals;
}

std::vector< std::pair<Individual*, Individual*>* >* Pedigree::get_relations() const {
  return m_relations;
}



Individual* Pedigree::get_root() {
  if (m_root != nullptr) {
    return m_root;
  }
  
  /* FIXME: Exploits tree */
  bool root_set = false;
  
  for (auto &individual : (*m_all_individuals)) {
    if (individual->get_father() == nullptr) {
      if (root_set) {
        Rcpp::stop("Only expected one root in male pedigree!");
      } else {
        m_root = individual;
        root_set = true;
      }
      
      break;
    }
  }

  if (!root_set || m_root == nullptr) {
    Rcpp::stop("Expected a root in male pedigree!");
  }

  return m_root;
}


void Pedigree::populate_haplotypes(
    int loci, 
    std::vector<double>& mutation_rates, 
    double prob_two_step) {

  if (prob_two_step < 0.0 || prob_two_step > 1.0) {
    Rcpp::stop("prob_two_step must be between 0.0 and 1.0");
  }
  
  /* FIXME: Exploits tree */
  Individual* root = this->get_root();
  
  std::vector<int> h(loci); // initialises to 0, 0, ..., 0

  root->set_haplotype(h);
  root->pass_haplotype_to_children(true, mutation_rates, prob_two_step);
}

void Pedigree::populate_haplotypes_custom_founders(
    std::vector<double>& mutation_rates, 
    Rcpp::Function get_founder_hap, 
    double prob_two_step) {
  
  if (prob_two_step < 0.0 || prob_two_step > 1.0) {
    Rcpp::stop("prob_two_step must be between 0.0 and 1.0");
  }
  
  /* FIXME: Exploits tree */
  Individual* root = this->get_root();
  
  std::vector<int> h = Rcpp::as< std::vector<int> >( get_founder_hap() );  
  
  //Rcpp::Rcout << "Unbounded: " << std::endl;
  //Rcpp::print(Rcpp::wrap(h));

  // Test that a haplotype of proper length generated  
  if (h.size() != mutation_rates.size()) {
    Rcpp::stop("get_founder_haplotype generated haplotype with number of loci different from the number of mutation rates specified");
  }
  
  //Rf_PrintValue(Rcpp::wrap(h));

  root->set_haplotype(h);
  root->pass_haplotype_to_children(true, mutation_rates, prob_two_step);
}

// Assumes ladders are okay, see checks in api_utility_haplotypes.cpp
void Pedigree::populate_haplotypes_ladder_bounded(
    std::vector<double>& mutation_rates, 
    std::vector<int>& ladder_min, 
    std::vector<int>& ladder_max, 
    Rcpp::Function get_founder_hap,
    double prob_two_step) {
  
  if (mutation_rates.size() != ladder_min.size()) {
    Rcpp::stop("mutation_rates and ladder_min must have same length");
  }
  
  if (mutation_rates.size() != ladder_min.size()) {
    Rcpp::stop("mutation_rates and ladder_max must have same length");
  }
  
  if (ladder_min.size() != ladder_max.size()) {
    Rcpp::stop("ladder_min and ladder_max must have same length");
  }

  if (prob_two_step < 0.0 || prob_two_step > 1.0) {
    Rcpp::stop("prob_two_step must be between 0.0 and 1.0");
  }

  /* Exploits tree */
  Individual* root = this->get_root();
  
  //std::vector<int> h(mutation_rates.size()); // initialises to 0, 0, ..., 0
  
  std::vector<int> h = Rcpp::as< std::vector<int> >( get_founder_hap() );  
  
  //Rcpp::Rcout << "Bounded: " << std::endl;
  //Rcpp::print(Rcpp::wrap(h));
  
  // Test that a haplotype of proper length generated  
  if (h.size() != mutation_rates.size()) {
    Rcpp::stop("get_founder_haplotype generated haplotype with number of loci different from the number of mutation rates specified");
  }
  
  //Rf_PrintValue(Rcpp::wrap(h));
  
  root->set_haplotype(h);
  root->pass_haplotype_to_children_ladder_bounded(true, mutation_rates, ladder_min, ladder_max, prob_two_step);
}


void Pedigree::populate_autosomal(
    const std::vector< std::vector<double> >& allele_conditional_cumdists_theta,
    const std::vector<double>& allele_cumdist_theta,
    const int alleles_count,
    const double mutation_rate) {
  
  /* Exploits tree */
  Individual* root = this->get_root();

  if (alleles_count <= 0) {
    Rcpp::stop("alleles_count must have at least size 1");
  }
    
  if (allele_cumdist_theta.size() <= 0) {
    Rcpp::stop("allele_cumdist_theta must have at least size 1");
  }
    
  if (allele_conditional_cumdists_theta.size() <= 0) {
    Rcpp::stop("allele_conditional_cumdists_theta must have at least size 1");
  }

  std::vector<int> h = draw_autosomal_genotype(allele_cumdist_theta, alleles_count);
  
  root->set_haplotype(h); // Not actually haplotype, but use this slot for lower memory footprint
  root->pass_autosomal_to_children(true, allele_conditional_cumdists_theta, mutation_rate);
}


