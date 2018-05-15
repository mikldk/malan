/**
 api_utility_haplotypes.cpp
 Purpose: Logic related to haplotypes.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */

//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>
#include <unordered_map>
#include <unordered_set>

#include "malan_types.h"
#include "api_utility_individual.h"


double estimate_theta_subpops_unweighted_geno_engine_HWE(
    const std::vector<int> ni, 
    std::unordered_map<int, std::vector<int>> type_counts) {
  
  
  ////////////////////////////////////////////////////
  // Calculating quantities
  ////////////////////////////////////////////////////
  
  int r = ni.size(); 
  double r_dbl = (double)r; 
  
  std::vector<double> M_within(r);
  Rcpp::NumericMatrix M_between(r, r);
  
  double MW = 0.0;
  double MB = 0.0;
  
  for (int i = 0; i < r; ++i) {
    double Mi = 0.0;
    double n_i = (double)ni[i];
    double denom = (n_i * (n_i - 1.0));
    
    for (auto ele = type_counts.begin(); ele != type_counts.end(); ++ele) {
      int nu = ele->second[i];
      Mi += nu * (nu - 1.0) / denom;
    }
    
    M_within[i] = Mi;
    MW += Mi / r_dbl;
  }
  
  for (int i1 = 0; i1 < (r-1); ++i1) {
    double ni1 = (double)ni[i1];

    for (int i2 = (i1+1); i2 < r; ++i2) {
      double Mi1i2 = 0.0;

      double ni2 = (double)ni[i2];
      double denom = ni1 * ni2;
      
      for (auto ele = type_counts.begin(); ele != type_counts.end(); ++ele) {
        int nu1 = ele->second[i1];
        int nu2 = ele->second[i2];
        
        Mi1i2 += nu1 * nu2 / denom;
      }
      
      M_between(i1, i2) = Mi1i2;
      MB += Mi1i2 / (r_dbl * (r_dbl - 1.0));
    }
  }
  
  // Factor 2 because we only used upper triangular matrix;
  // now we take all.
  MB *= 2.0;
  
  /*
   Rcpp::Rcout << "M_within:" << std::endl;
   Rcpp::print(Rcpp::wrap(M_within));
   
   Rcpp::Rcout << "M_between:" << std::endl;
   Rcpp::print(M_between);
   
   Rcpp::Rcout << "MW:" << std::endl;
   Rcpp::print(Rcpp::wrap(MW));
   
   Rcpp::Rcout << "MB:" << std::endl;
   Rcpp::print(Rcpp::wrap(MB));
  */
  
  return (MW - MB) / (1.0 - MB);
}

std::pair<int, int> get_ordered_genotype(int a1, int a2) {

  if (a1 > a2) {
    int tmp = a1;
    a1 = a2;
    a2 = tmp;
  }
  
  std::pair<int, int> geno = std::make_pair(a1, a2);
  return geno;
}

void fill_count_hashmap_theta_unweighted_HWE(
    const int r, // subpops
    const int i, // subpop index
    const int a1, const int a2, // genotype
    std::unordered_map<int, std::vector<int>>& type_counts_allele) {
  
  // a1
  if (type_counts_allele.find(a1) == type_counts_allele.end()) {
    type_counts_allele[a1].resize(r); // values defaults to 0
  }
  type_counts_allele[a1][i] += 1;  

  // a2
  if (type_counts_allele.find(a2) == type_counts_allele.end()) {
    type_counts_allele[a2].resize(r); // values defaults to 0
  }
  type_counts_allele[a2][i] += 1;
}


//' Unweighted estimate of theta from subpopulations of genotypes
//' 
//' Estimates unweighted theta for a number of subpopulations given a list of subpopulations of genotypes.
//' 
//' Based on Weir and Goudet, Genetics 2017: 
//' http://www.genetics.org/content/early/2017/05/26/genetics.116.198424
//' 
//' @param subpops List of individual genotypes
//' @param assume_HWE if the alleles themselves are used instead of genotypes
//' 
//' @return Estimate of theta
//' 
//' @export
// [[Rcpp::export]]
double estimate_theta_subpops_unweighted_genotypes(Rcpp::ListOf<Rcpp::IntegerMatrix> subpops, 
                                                   bool assume_HWE) {

  if (assume_HWE == false) {
    Rcpp::stop("Not yet implemented");
  }
    
  const int r = subpops.size();

  if (r <= 0) {
    Rcpp::stop("No subpopulations given");
  }
  
  // Types:
  std::unordered_map<int, std::vector<int>> type_counts_allele;
  std::vector<int> ni(r);
    
  ////////////////////////////////////////////////////
  // Filling count container
  ////////////////////////////////////////////////////
  for (int i = 0; i < r; ++i) {
    Rcpp::IntegerMatrix subpop = subpops[i];
    const int sample_size_i = subpop.nrow();
        
    if (sample_size_i <= 0) {
      Rcpp::stop("Subpop sample of size <= 0");
    }
    
    if (subpop.ncol() != 2) {
      Rcpp::stop("Expected exactly 2 autosomal loci");
    }

    ni[i] = 2*sample_size_i; // 2 alleles per individual
    
    for (int j = 0; j < sample_size_i; ++j) {
      Rcpp::IntegerVector hap = subpop(j, Rcpp::_);
      
      if (hap.size() != 2) {
        Rcpp::stop("Expected exactly 2 autosomal loci");
      }  
      
      fill_count_hashmap_theta_unweighted_HWE(r, i, hap[0], hap[1], type_counts_allele);
    }
  }
  
  /*
  // print
  for (int i = 0; i < r; ++i) Rcpp::Rcout << ni[i] << ", "; 
  Rcpp::Rcout <<  std::endl;
  for (auto it : type_counts_allele) { 
    Rcpp::Rcout << "    allele " << it.first << ": ";     
    for (int i = 0; i < r; ++i) Rcpp::Rcout << it.second[i] << ", ";    
    Rcpp::Rcout <<  std::endl;
  } 
  */
  
  double res = estimate_theta_subpops_unweighted_geno_engine_HWE(
    ni, type_counts_allele);
  
  return res;
}



//' Unweighted estimate of theta from subpopulations of individual ids
//' 
//' Estimates unweighted theta for a number of subpopulations given a list of pids (individual ids).
//' 
//' Based on Weir and Goudet, Genetics 2017: 
//' http://www.genetics.org/content/early/2017/05/26/genetics.116.198424
//' 
//' @param population Population obtain from simulation
//' @param subpops List of individual pids
//' @param assume_HWE if the alleles themselves are used instead of genotypes
//' 
//' @return Estimate of theta
//' 
//' @export
// [[Rcpp::export]]
double estimate_theta_subpops_unweighted_pids(Rcpp::XPtr<Population> population,
                                              Rcpp::ListOf<Rcpp::IntegerVector> subpops,
                                              bool assume_HWE) {
  if (assume_HWE == false) {
    Rcpp::stop("Not yet implemented");
  }
  
  const int r = subpops.size();
  
  if (r <= 0) {
    Rcpp::stop("No subpopulations given");
  }
  
  // Types:
  std::unordered_map<int, std::vector<int>> type_counts_allele;
  std::vector<int> ni(r);
  
  ////////////////////////////////////////////////////
  // Filling count container
  ////////////////////////////////////////////////////
  for (int i = 0; i < r; ++i) {
    Rcpp::IntegerVector subpop_pids = subpops[i];

    const int sample_size_i = subpop_pids.size();

    if (sample_size_i <= 0) {
      Rcpp::stop("Subpop sample of size <= 0");

    }
    
    ni[i] = 2*sample_size_i;
    
    for (int j = 0; j < sample_size_i; ++j) {
      const int pid = subpop_pids[j];
      const Individual* individual = population->get_individual(pid);
      
      if (!(individual->is_haplotype_set())) {
        Rcpp::stop("Haplotypes not yet set");
      }      

      const std::vector<int> hap = individual->get_haplotype();
      
      if (hap.size() != 2) {
        Rcpp::stop("Expected exactly 2 autosomal loci");
      }  
      
      
      fill_count_hashmap_theta_unweighted_HWE(r, i, hap[0], hap[1], type_counts_allele);
    }
  }
  
  double res = estimate_theta_subpops_unweighted_geno_engine_HWE(
    ni, type_counts_allele);
  
  return res;
}









// All vectors must be of size r
Rcpp::IntegerMatrix convert_map_to_matrix(const int r, const std::unordered_map<int, std::vector<int>>& map) { 
  Rcpp::CharacterVector col_nms(map.size());
  Rcpp::IntegerMatrix res(r, map.size());
  int col = 0;
  
  for (auto it : map) { 
    col_nms[col] = it.first;
    
    for (int i = 0; i < r; ++i) {
      res(i, col) = it.second[i];
    }
    
    col += 1;
  } 
  
  colnames(res) = col_nms;
  
  return res;
}

//' Get allele counts from subpopulations of genotypes
//' 
//' @param subpops List of individual genotypes
//' 
//' @return Matrix with allele counts
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix get_allele_counts_genotypes(Rcpp::ListOf<Rcpp::IntegerMatrix> subpops) {
  const int r = subpops.size();

  if (r <= 0) {
    Rcpp::stop("No subpopulations given");
  }
  
  // Types:
  std::unordered_map<int, std::vector<int>> type_counts_allele;

  for (int i = 0; i < r; ++i) {
    Rcpp::IntegerMatrix subpop = subpops[i];    
    const int sample_size_i = subpop.nrow();
    
    if (sample_size_i <= 0) {
      Rcpp::stop("Subpop sample of size <= 0");
    }
    
    if (subpop.ncol() != 2) {
      Rcpp::stop("Expected exactly 2 autosomal loci");
    }
    
    
    for (int j = 0; j < sample_size_i; ++j) {
      Rcpp::IntegerVector hap = subpop(j, Rcpp::_);
      
      if (hap.size() != 2) {
        Rcpp::stop("Expected exactly 2 autosomal loci");
      }  
      
      fill_count_hashmap_theta_unweighted_HWE(r, i, hap[0], hap[1], type_counts_allele);
    }
  }
  
  Rcpp::IntegerMatrix res = convert_map_to_matrix(r, type_counts_allele); 
  
  return res;
}


//' Get allele counts from subpopulations given by pids
//' 
//' @param population Population obtain from simulation
//' @param subpops List of individual pids
//' 
//' @return Matrix with allele counts
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix get_allele_counts_pids(Rcpp::XPtr<Population> population,
                                           Rcpp::ListOf<Rcpp::IntegerVector> subpops) {
  const int r = subpops.size();

  if (r <= 0) {
    Rcpp::stop("No subpopulations given");
  }
  
  // Types:
  std::unordered_map<int, std::vector<int>> type_counts_allele;

  for (int i = 0; i < r; ++i) {
    Rcpp::IntegerVector subpop_pids = subpops[i];

    const int sample_size_i = subpop_pids.size();

    if (sample_size_i <= 0) {
      Rcpp::stop("Subpop sample of size <= 0");
    }    
    
    for (int j = 0; j < sample_size_i; ++j) {
      const int pid = subpop_pids[j];
      const Individual* individual = population->get_individual(pid);
      
      if (!(individual->is_haplotype_set())) {
        Rcpp::stop("Haplotypes not yet set");
      }      

      const std::vector<int> hap = individual->get_haplotype();
      
      if (hap.size() != 2) {
        Rcpp::stop("Expected exactly 2 autosomal loci");
      }        
      
      fill_count_hashmap_theta_unweighted_HWE(r, i, hap[0], hap[1], type_counts_allele);
    }
  }
  
  Rcpp::IntegerMatrix res = convert_map_to_matrix(r, type_counts_allele); 
  
  return res;
}

