/**
 api_utility_individual.cpp
 Purpose: Logic related to individuals.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */
 
#include "malan_types.h"
#include "api_utility_individual.h"

//' Get individual by pid
//' 
//' @param population Population
//' @param pid pid
//' 
//' @return Individual
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid) {  
  Individual* ind = population->get_individual(pid);
  Rcpp::XPtr<Individual> res(ind, RCPP_XPTR_2ND_ARG); // do NOT delete individual when not used any more, it still exists in pedigree and population etc.!
  res.attr("class") = Rcpp::CharacterVector::create("malan_individual", "externalptr");
  
  return res;
}


//' Get pid from individual
//' 
//' @param individual Individual to get pid of
//' 
//' @return pid
//' 
//' @export
// [[Rcpp::export]]
int get_pid(Rcpp::XPtr<Individual> individual) {  
  return individual->get_pid();
}

//' Print individual
//' 
//' @param individual Individual
//' 
//' @export
// [[Rcpp::export]]
void print_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  int pid_f = (i->get_father() != nullptr) ? i->get_father()->get_pid() : -1;
  std::vector<Individual*>* children = i->get_children();
  
  Rcpp::Rcout << "  pid = " << i->get_pid() << " with father pid = " << pid_f << " and";
  
  if (children->size() == 0) {
    Rcpp::Rcout << " no children" << std::endl;
  } else {
    Rcpp::Rcout << " children (n = " << children->size() << "): " << std::endl;

    for (auto child : *children) {    
      std::vector<Individual*>* child_children = child->get_children();
      
      Rcpp::Rcout << "    pid = " << child->get_pid() << " with father pid = " << pid_f << " and " <<  child_children->size() << " children" << std::endl;
    }
  }
}

//' Get individual's generation number
//' 
//' Note that generation 0 is final, end generation. 
//' 1 is second last generation etc.
//' 
//' @param individual Individual
//' 
//' @return generation
//' 
//' @export
// [[Rcpp::export]]
int get_generation(Rcpp::XPtr<Individual> individual) {  
  return individual->get_generation();
}

//' Get pedigree from individual
//' 
//' @param individual Individual
//' 
//' @return pedigree
//' 
//' @export
//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;  
  Rcpp::XPtr<Pedigree> res(i->get_pedigree(), RCPP_XPTR_2ND_ARG); // do NOT delete pedigree when not used any more, it still exists in list of pedigrees etc.!
  res.attr("class") = Rcpp::CharacterVector::create("malan_pedigree", "externalptr");
  
  return res;
}

//' Get pedigree ids from pids
//'
//' @param population Population
//' @param pids Pids
//' 
//' @return Vector with pedigree ids
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, 
                                             Rcpp::IntegerVector pids) {  
  
  std::unordered_map<int, Individual*> pop = *(population->get_population());
  
  int N = pids.size();
  Rcpp::IntegerVector pedigree_ids(N);
  
  for (int i = 0; i < N; ++i) {
    Individual* ind = population->get_individual(pids[i]);
    pedigree_ids[i] = ind->get_pedigree_id();
  }
  
  return pedigree_ids;
}



//////////////////////////////////////

//' Get individual's family information
//'
//' @param individual individual
//' 
//' @return List with family information
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_family_info(Rcpp::XPtr<Individual> individual) {  
  return Rcpp::List::create(
    Rcpp::Named("num_brothers") = count_brothers(individual),
    Rcpp::Named("num_brothers_matching") = brothers_matching(individual),
    Rcpp::Named("father_matches") = father_matches(individual),
    Rcpp::Named("grandfather_matches") = grandfather_matches(individual),
    Rcpp::Named("num_uncles") = count_uncles(individual));

}


//' Get father
//' 
//' Get individual's father
//'
//' @param individual individual
//' 
//' @return Father
//' 
//' @seealso [get_brothers()], [get_uncles()], [get_children()], [get_cousins()]
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<Individual> get_father(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  if (i->get_father() == nullptr) {
    Rcpp::stop("Individual did not have a father");
  }
  
  Individual* father = i->get_father();
  Rcpp::XPtr<Individual> father_xptr(father, RCPP_XPTR_2ND_ARG); // do NOT delete individual when not used any more, it still exists in pedigree and population etc.!
  father_xptr.attr("class") = Rcpp::CharacterVector::create("malan_individual", "externalptr");

  return father_xptr;
}

//' Get children
//' 
//' Get individual's children
//'
//' @param individual individual
//' 
//' @return List with children
//' 
//' @seealso [get_father()], [get_brothers()], [get_uncles()], [get_cousins()]
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_children(Rcpp::XPtr<Individual> individual) {  
  std::vector<Individual*>* individuals_boys = individual->get_children();
  
  Rcpp::List children;
  
  for (auto child : *individuals_boys) {
    Rcpp::XPtr<Individual> child_xptr(child, RCPP_XPTR_2ND_ARG); // do NOT delete individual when not used any more, it still exists in pedigree and population etc.!
    child_xptr.attr("class") = Rcpp::CharacterVector::create("malan_individual", "externalptr");
  
    children.push_back(child_xptr);
  }
  
  return children;
}


//' Number of brothers
//' 
//' Get individual's number of brothers
//'
//' @param individual individual
//' 
//' @return Number of brothers
//' 
//' @seealso [get_brothers()]
//' 
//' @export
// [[Rcpp::export]]
int count_brothers(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  if (i->get_father() == nullptr) {
    Rcpp::stop("Individual did not have a father");
  }
  
  int fathers_boys = i->get_father()->get_children_count();
  // exclude individual
  return (fathers_boys - 1);
}


//' Get brothers
//' 
//' Get individual's brothers
//'
//' @param individual individual
//' 
//' @return List with brothers
//' 
//' @seealso [get_father()], [get_uncles()], [get_children()], [get_cousins()]
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_brothers(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  if (i->get_father() == nullptr) {
    Rcpp::stop("Individual did not have a father");
  }
  
  Individual* father = i->get_father();
  
  std::vector<Individual*>* father_boys = father->get_children();
  
  Rcpp::List brothers;
  
  for (auto brother : *father_boys) {
    if (brother->get_pid() == individual->get_pid()) {
      continue; // exclude individual itself as a brother
    }
    
    Rcpp::XPtr<Individual> brother_xptr(brother, RCPP_XPTR_2ND_ARG); // do NOT delete individual when not used any more, it still exists in pedigree and population etc.!
    brother_xptr.attr("class") = Rcpp::CharacterVector::create("malan_individual", "externalptr");
  
    brothers.push_back(brother_xptr);
  }
  
  return brothers;
}


//' Number of brothers with matching haplotype
//' 
//' Get individual's number of brothers that matches `individual`'s haplotype
//'
//' @param individual individual
//' 
//' @return Number of brothers that matches `individual`'s haplotype
//' 
//' @export
// [[Rcpp::export]]
int brothers_matching(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  if (i->get_father() == nullptr) {
    Rcpp::stop("Individual did not have a father");
  }
  
  if (!(i->is_haplotype_set())) {
    Rcpp::stop("Individual did not have a haplotype");
  }

  std::vector<int> h = i->get_haplotype();  
  int loci = h.size();  
  
  std::vector<Individual*>* brothers = i->get_father()->get_children();
  
  if (brothers->size() == 0) {
    return 0;
  }
  
  int matching = 0;

  for (auto brother : *brothers) {
    if (brother->get_pid() == i->get_pid()) {
      continue;
    }

    if (!(brother->is_haplotype_set())) {
      Rcpp::stop("Individual's brother did not have a haplotype");
    }  
  
    std::vector<int> indv_h = brother->get_haplotype();    
    if (indv_h.size() != loci) {
      Rcpp::stop("haplotype and indv_h did not have same number of loci");
    }
    
    if (indv_h == h) {
      matching += 1;
    }
  }
  
  return matching;
}

//' Father matches
//' 
//' Does the father have the same profile as `individual`?
//'
//' @param individual individual
//' 
//' @return Whether father has the same profile as `individual` or not
//' 
//' @export
// [[Rcpp::export]]
bool father_matches(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;

  if (!(i->is_haplotype_set())) {
    Rcpp::stop("Individual did not have a haplotype");
  }
    
  if (i->get_father() == nullptr) {
    Rcpp::stop("Individual did not have a father");
  }
  
  Individual* father = i->get_father();
  
  if (!(father->is_haplotype_set())) {
    Rcpp::stop("Individual's father did not have a haplotype");
  }  
  
  std::vector<int> h = i->get_haplotype();
  std::vector<int> h_father = father->get_haplotype();    
  
  return (h.size() == h_father.size() && h == h_father);
}

//' Grandfather matches
//' 
//' Does the frandfather have the same profile as `individual`?
//'
//' @param individual individual
//' 
//' @return Whether grandfather has the same profile as `individual` or not
//' 
//' @export
// [[Rcpp::export]]
bool grandfather_matches(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;

  if (!(i->is_haplotype_set())) {
    Rcpp::stop("Individual did not have a haplotype");
  }
  
  
  if (i->get_father() == nullptr) {
    Rcpp::stop("Individual did not have a father");
  }  
  Individual* father = i->get_father();
  if (!(father->is_haplotype_set())) {
    Rcpp::stop("Individual's father did not have a haplotype");
  }
  
  // It is not sufficient to calculate father_matches(father) as
  // father and grandfather may not match

  if (father->get_father() == nullptr) {
    Rcpp::stop("Individual's father did not have a father");
  }
  Individual* grandfather = father->get_father();  
  if (!(grandfather->is_haplotype_set())) {
    Rcpp::stop("Individual's grandfather did not have a haplotype");
  }  
  
  std::vector<int> h = i->get_haplotype();
  std::vector<int> h_grandfather = grandfather->get_haplotype();    
  
  return (h.size() == h_grandfather.size() && h == h_grandfather);
}

//' Number of uncles
//' 
//' Get individual's number of uncles
//'
//' @param individual individual
//' 
//' @return Number of uncles
//' 
//' @seealso [get_uncles()]
//' 
//' @export
// [[Rcpp::export]]
int count_uncles(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  if (i->get_father() == nullptr) {
    Rcpp::stop("Individual did not have a father");
  }
  
  Individual* father = i->get_father();
  
  if (father->get_father() == nullptr) {
    Rcpp::stop("Individual's father did not have a father");
  }  
  
  int fathers_fathers_boys = father->get_father()->get_children_count();

  // exclude father
  return (fathers_fathers_boys - 1);
}

//' Get uncles
//' 
//' Get individual's uncles
//'
//' @param individual individual
//' 
//' @return List with uncles
//' 
//' @seealso [get_brothers()], [get_children()], [get_cousins()]
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_uncles(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  if (i->get_father() == nullptr) {
    Rcpp::stop("Individual did not have a father");
  }
  
  Individual* father = i->get_father();
  
  if (father->get_father() == nullptr) {
    Rcpp::stop("Individual's father did not have a father");
  }  

  Individual* grandfather = father->get_father();

  std::vector<Individual*>* grandfather_boys = grandfather->get_children();
  
  Rcpp::List uncles;
  
  for (auto uncle : *grandfather_boys) {
    if (uncle->get_pid() == father->get_pid()) {
      continue; // exclude father as uncle
    }
    
    Rcpp::XPtr<Individual> uncle_xptr(uncle, RCPP_XPTR_2ND_ARG); // do NOT delete individual when not used any more, it still exists in pedigree and population etc.!
    uncle_xptr.attr("class") = Rcpp::CharacterVector::create("malan_individual", "externalptr");
  
    uncles.push_back(uncle_xptr);
  }
  
  return uncles;
}


//' Get cousins
//' 
//' Get individual's cousins
//'
//' @param individual individual
//' 
//' @return List with cousins
//' 
//' @seealso [get_brothers()], [get_uncles()], [get_children()]
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List get_cousins(Rcpp::XPtr<Individual> individual) {  
  Rcpp::List uncles = get_uncles(individual);  
  Rcpp::List cousins;
  int n = uncles.size();
  
  for (int i = 0; i < n; ++i) {
    Rcpp::List uncles_children = get_children(uncles[i]);
    
    int m = uncles_children.size();
    
    for (int j = 0; j < m; ++j) {
      cousins.push_back(uncles_children[j]);
    }
  }
  
  return cousins;
}


