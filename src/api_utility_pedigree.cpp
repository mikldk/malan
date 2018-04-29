/**
 api_utility_pedigree.cpp
 Purpose: Logic related to pedigrees.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "malan_types.h"

//' Get pedigree id
//' 
//' @param ped Pedigree
//' 
//' @export
// [[Rcpp::export]]
int get_pedigree_id(Rcpp::XPtr<Pedigree> ped) { 
  return ped->get_id();
}

//' Get number of pedigrees
//' 
//' @param pedigrees Pedigrees
//' 
//' @export
// [[Rcpp::export]]
int pedigrees_count(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees) {
  return pedigrees->size();
}

//' Get pedigree size
//' 
//' @param ped Pedigree
//' 
//' @export
// [[Rcpp::export]]
int pedigree_size(Rcpp::XPtr<Pedigree> ped) {  
  return ped->get_all_individuals()->size();
}

//' Get distribution of pedigree sizes
//' 
//' @param pedigrees Pedigrees
//' 
//' @export
//[[Rcpp::export]]
std::unordered_map<int, int> pedigrees_table(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees) {
  std::vector<Pedigree*>* peds = pedigrees;
  std::unordered_map<int, int> tab;
  
  for (auto it = peds->begin(); it != peds->end(); ++it) {
    tab[(*it)->get_all_individuals()->size()] += 1;
  }
  
  return tab;
}

//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int index) {  
  std::vector<Pedigree*>* peds = pedigrees;
  Pedigree* p = peds->at(index);
  
  Rcpp::XPtr<Pedigree> res(p, RCPP_XPTR_2ND_ARG); // do NOT delete pedigree when not used any more, it still exists in list of pedigrees etc.!
  res.attr("class") = Rcpp::CharacterVector::create("malan_pedigree", "externalptr");
  
  return res;
}




//[[Rcpp::export]]
void print_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  std::vector< std::pair<Individual*, Individual*>* >* rels = p->get_relations();
  
  Rcpp::Rcout << "Pedigree with " << p->get_all_individuals()->size() << " individuals:" << std::endl;
  
  for (auto i : *inds) {    
    int pid_f = (i->get_father() != NULL) ? i->get_father()->get_pid() : -1;
    
    Rcpp::Rcout << "  " << i->get_pid() << " with father " << pid_f << std::endl;
  } 
}

//' Get pids in pedigree
//' 
//' @param ped Pedigree
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector get_pids_in_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  Rcpp::IntegerVector res(inds->size());
  int i = 0;
  for (auto ind : *inds) {   
    res(i) = ind->get_pid();
    ++i;
  } 
  
  return res;
}

//' Get haplotypes in pedigree
//' 
//' @param ped Pedigree
//' 
//' @return List with haplotypes
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_haplotypes_in_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
 
  size_t N = inds->size();
  Rcpp::List haps(N); 
  
  for (size_t i = 0; i < N; ++i) {
    Individual* indv = inds->at(i);
    haps(i) = indv->get_haplotype();
  }
  
  return haps;
}

//[[Rcpp::export]]
Rcpp::CharacterMatrix get_pedigree_edgelist(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector< std::pair<Individual*, Individual*>* >* rels = p->get_relations();
  
  Rcpp::CharacterMatrix edgelist(rels->size(), 2);
  int i = 0;
  
  for (auto pair: *rels) {
    edgelist(i, 0) = std::to_string(pair->first->get_pid());
    edgelist(i, 1) = std::to_string(pair->second->get_pid());
    ++i;
  }
  
  return edgelist;
}


//' Get pedigree information as graph (mainly intended for plotting)
//' 
//' @param ped Pedigree
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_pedigree_as_graph(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  Rcpp::CharacterVector nodes(inds->size());
  
  int i = 0;
  for (auto individual : *inds) {
    nodes(i) = std::to_string(individual->get_pid());   
    ++i;
  }
  
  Rcpp::List ret;
  ret["nodes"] = nodes;
  ret["edgelist"] = get_pedigree_edgelist(ped);
  
  return ret;
}




//' Get pedigrees information in tidy format
//' 
//' @param pedigrees Pedigrees
//' 
// [[Rcpp::export]]
Rcpp::List get_pedigrees_tidy(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees) {  
  std::vector<Pedigree*>* peds = pedigrees;
  
  Rcpp::List ret_ped_ids;
  Rcpp::List ret_edgelists;
  Rcpp::List ret_haplotypes;
  Rcpp::List ret_pids;
  Rcpp::List ret_generation;  
  
  for (auto it = peds->begin(); it != peds->end(); ++it) {
    Pedigree* ped = *it;

    ret_ped_ids.push_back(ped->get_id());    
    
    std::vector< std::pair<Individual*, Individual*>* >* rels = ped->get_relations();
    
    Rcpp::IntegerMatrix edgelist(rels->size(), 2);
    int i = 0;
    
    for (auto pair: *rels) {
      edgelist(i, 0) = pair->first->get_pid();
      edgelist(i, 1) = pair->second->get_pid();
      ++i;
    }
    
    ret_edgelists.push_back(edgelist);

    
    std::vector<Individual*>* inds = ped->get_all_individuals();
    
    size_t N = inds->size();
    Rcpp::List haps(N);
    Rcpp::IntegerVector pids(N);
    Rcpp::IntegerVector generation(N);
    
    for (size_t i = 0; i < N; ++i) {
      Individual* indv = inds->at(i);
      haps(i) = indv->get_haplotype();
      pids(i) = indv->get_pid();
      generation(i) = indv->get_generation();
    }
    
    ret_haplotypes.push_back(haps);
    ret_pids.push_back(pids);
    ret_generation.push_back(generation);
  }
  
  Rcpp::List ret;
  
  ret["ped_ids"] = ret_ped_ids;
  ret["edgelists"] = ret_edgelists;
  ret["haplotypes"] = ret_haplotypes;
  ret["pids"] = ret_pids;
  ret["generations_from_final"] = ret_generation;
  
  return ret;
}

