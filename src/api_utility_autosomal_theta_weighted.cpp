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

// [[Rcpp::export]]
std::unordered_map<int, int> hash_colisions(int p) {
  std::unordered_map<int, int> tab;  
  
  for (int i = 0; i < (p-1); ++i) {
    for (int j = (i+1); j < p; ++j) {
      std::pair<int, int> pair = std::make_pair(i, j);
      int hash = std::hash< std::pair<int, int> >()(pair);
      tab[hash] += 1;
    }
  }
  
  return tab;
}




Rcpp::List estimate_theta_1subpop(const std::unordered_map<int, double>& allele_p,
                                  const std::unordered_map<std::pair<int, int>, double>& genotype_p,
                                  const std::unordered_set<std::pair<int, int>>& genotypes_unique,
                                  const bool return_estimation_info = false) {
  Rcpp::List theta;
  
  theta["estimate"] = NA_REAL;
  theta["error"] = true;
  theta["details"] = "NA";
  theta["estimation_info"] = R_NilValue;  

  // Loop over unique genotypes
  std::unordered_set<std::pair<int, int>>::const_iterator it;
  int K = genotypes_unique.size();
  int k = 0;
  
  k = 0;
  arma::mat X(K, 1, arma::fill::none);
  arma::vec y(K, arma::fill::none);
  
  for (it = genotypes_unique.begin(); it != genotypes_unique.end(); ++it) {
    std::pair<int, int> geno = *it;
    int a1 = geno.first;
    int a2 = geno.second;
    
    // homozyg
    if (a1 == a2) {
      double p_i = allele_p.at(a1);
      double p_ii = genotype_p.at(geno);
      double p_i2 = p_i*p_i;
      X(k, 0) = p_i - p_i2;
      y(k) = p_ii - p_i2;
    } else {
      // heterozyg
      double p_i = allele_p.at(a1);
      double p_j = allele_p.at(a2);
      double p_ij = genotype_p.at(geno);
      double tmp = -2.0*p_i*p_j;
      X(k, 0) = tmp;
      y(k) = p_ij + tmp;
    }
    
    ++k;
  }

  if (return_estimation_info) {
    Rcpp::List est_info;
    est_info["X"] = Rcpp::wrap(X);
    est_info["y"] = y;
    
    k = 0;
    
    Rcpp::IntegerMatrix genotypes(K, 2);
    Rcpp::NumericVector genoptype_probs(K);
    Rcpp::NumericMatrix geno_allele_probs(K, 2);
    Rcpp::IntegerVector zygosity(K);
      
    for (it = genotypes_unique.begin(); it != genotypes_unique.end(); ++it) {
      std::pair<int, int> geno = *it;
      int a1 = geno.first;
      int a2 = geno.second;
      genotypes(k, 0) = a1;
      genotypes(k, 1) = a2;
      genoptype_probs[k] = genotype_p.at(geno);
      
      // homozyg
      if (a1 == a2) {
        zygosity[k] = 1;
        
        double p_i = allele_p.at(a1);
        geno_allele_probs(k, 0) = p_i;
        geno_allele_probs(k, 1) = p_i;
      } else {
        // heterozyg
        zygosity[k] = 2;

        double p_i = allele_p.at(a1);
        double p_j = allele_p.at(a2);
        
        geno_allele_probs(k, 0) = p_i;
        geno_allele_probs(k, 1) = p_j;
      }
      
      ++k;
    }

    est_info["genotypes"] = genotypes;
    est_info["genotypes_zygosity"] = zygosity;
    est_info["genotypes_probs"] = genoptype_probs;
    est_info["genotypes_allele_probs"] = geno_allele_probs;
    
    std::vector<int> alleles_names;
    alleles_names.reserve(allele_p.size());
    std::vector<double> alleles_probs;
    alleles_probs.reserve(allele_p.size());
    
    for (auto it = allele_p.begin(); it != allele_p.end(); ++it) {
      alleles_names.push_back(it->first);
      alleles_probs.push_back(it->second);
    }
    est_info["alleles"] = alleles_names;
    est_info["alleles_probs"] = alleles_probs;

    theta["estimation_info"] = est_info;
  }

  if (K == 1) {
    theta["estimate"] = NA_REAL;
    theta["error"] = true;
    theta["details"] = "Only one genotype observed";
    return theta;
  }
  
  // minimisze (Xb - y)^2 for b
  arma::mat Q, R;
  bool status = arma::qr_econ(Q, R, X);
  
  if (!status) {
    theta["estimate"] = NA_REAL;
    theta["error"] = true;
    theta["details"] = "Could not make QR decomposition";
  } else {
    arma::vec coef = arma::solve(R, Q.t() * y, arma::solve_opts::no_approx);
    
    if (coef[0] >= 0 && coef[0] <= 1) {
      theta["estimate"] = coef[0];
      theta["error"] = false;
      theta["details"] = "OK";
    } else {
      theta["estimate"] = coef[0];
      theta["error"] = true;
      theta["details"] = "Estimate outside range of (0, 1)";
    }
  }
  
  return theta;
}


void estimate_theta_1subpop_fill_containers(int a1,
                                            int a2,
                                            const double one_over_n,
                                            const double one_over_2n,
                                            std::unordered_map<int, double>& allele_p,
                                            std::unordered_map<std::pair<int, int>, double>& genotype_p,
                                            std::unordered_set<std::pair<int, int>>& genotypes_unique) {
  
  if (a2 < a1) {
    int tmp = a1;
    a1 = a2;
    a2 = tmp;
  }
  
  std::pair<int, int> geno = std::make_pair(a1, a2);
  genotypes_unique.insert(geno);
  
  genotype_p[geno] += one_over_n;
  
  if (a1 == a2) {
    allele_p[a1] += one_over_n; // 2*one_over_2n = one_over_n
  } else {
    allele_p[a1] += one_over_2n;
    allele_p[a2] += one_over_2n;
  }
}


//' Estimate theta from genotypes
//' 
//' Estimate theta for one subpopulation given a sample of genotypes.
//' 
//' @param genotypes Matrix of genotypes: two columns (allele1 and allele2) and a row per individual
//' @param return_estimation_info Whether to return the quantities used to estimate `theta`
//' 
//' @return List:
//' * `theta`
//'     + `estimate`: Vector of length 1 containing estimate of theta or NA if it could not be estimated
//'     + `error`: true if an error happened, false otherwise
//'     + `details`: contains description if an error happened
//'     + `estimation_info`: If `return_estimation_info = true`: a list with information used to estimate `theta`. Else `NULL`.
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_theta_1subpop_genotypes(Rcpp::IntegerMatrix genotypes, bool return_estimation_info = false) {
  int n = genotypes.nrow();
  
  if (n <= 0) {
    Rcpp::stop("genotypes cannot be empty");
  }
  
  if (genotypes.ncol() != 2) {
    Rcpp::stop("genotypes must have exactly two columns");
  }
  
  // Build count tables
  std::unordered_map<int, double> allele_p;
  std::unordered_map<std::pair<int, int>, double> genotype_p;
  std::unordered_set<std::pair<int, int>> genotypes_unique;
  
  double one_over_n = 1.0 / (double)n;
  double one_over_2n = 1.0 / (2.0 * (double)n);
  
  for (int i = 0; i < n; ++i) {
    int a1 = genotypes(i, 0);
    int a2 = genotypes(i, 1);
    
    estimate_theta_1subpop_fill_containers(a1, a2, one_over_n, one_over_2n, 
                                           allele_p, genotype_p, genotypes_unique);
  }
  
  Rcpp::List theta = estimate_theta_1subpop(allele_p, genotype_p, genotypes_unique, 
                                            return_estimation_info);
    
  return theta;
}


//' Estimate theta from individuals
//' 
//' Estimate theta for one subpopulation given a list of individuals.
//' 
//' @inheritParams estimate_theta_1subpop_genotypes
//' @param individuals Individuals to get haplotypes for.
//' 
//' @inherit estimate_theta_1subpop_genotypes return
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_theta_1subpop_individuals(Rcpp::ListOf< Rcpp::XPtr<Individual> > individuals, 
                                              bool return_estimation_info = false) {
  
  int n = individuals.size();
  
  if (n <= 0) {
    Rcpp::stop("No individuals given");
  }
  
  if (!(individuals[0]->is_haplotype_set())) {
    Rcpp::stop("Haplotypes not yet set");
  }
  
  int loci = individuals[0]->get_haplotype().size();
  
  if (loci != 2) {
    Rcpp::stop("Expected exactly 2 autosomal loci");
  }
  
  // Build count tables
  std::unordered_map<int, double> allele_p;
  std::unordered_map<std::pair<int, int>, double> genotype_p;
  std::unordered_set<std::pair<int, int>> genotypes_unique;
  
  double one_over_n = 1.0 / (double)n;
  double one_over_2n = 1.0 / (2.0 * (double)n);

  for (int i = 0; i < n; ++i) {
    Individual* individual = individuals[i];
    std::vector<int> hap = individual->get_haplotype();
    
    estimate_theta_1subpop_fill_containers(hap[0], hap[1], one_over_n, one_over_2n, 
                                           allele_p, genotype_p, genotypes_unique);
  }
  
  Rcpp::List theta = estimate_theta_1subpop(allele_p, genotype_p, genotypes_unique, 
                                            return_estimation_info);
  return theta;
}

















void print_map(std::unordered_map<int, double> x) {
  for (auto j = x.begin(); j != x.end(); ++j) { 
    Rcpp::Rcout << "    allele " << j->first << ": " << j->second << std::endl; 
  } 
}

void print_container(std::string headline, std::vector< std::unordered_map<int, double> > x) {
  Rcpp::Rcout << "===========================================\n";
  Rcpp::Rcout << headline << "\n";
  Rcpp::Rcout << "===========================================\n";
  
  for (auto i = x.begin(); i != x.end(); ++i) { 
    Rcpp::Rcout << "  subpop " << (i-x.begin()) << std::endl;
    print_map(*i);
  }
}


/*
// P_AA: Homozygous probabilities:
// P_AA[i][l]: Homozygous probability of allele l in subpopulation i.
std::vector< std::unordered_map<int, double> > P_AA(r);

// p_A: Allele probability:
// p_A[i][l]: Allele probability of allele l in subpopulation i.
*/
Rcpp::List estimate_theta_subpops_weighted_engine(
    std::vector< std::unordered_map<int, double> > P_AA,
    std::vector< std::unordered_map<int, double> > p_A,
    std::vector<double> n) {
  
  int r = P_AA.size();
  
  if (r <= 0) {
    Rcpp::stop("r <= 0");
  }
  
  if (P_AA.size() != r) {
    Rcpp::stop("P_AA.size() != r");
  }
  
  if (p_A.size() != r) {
    Rcpp::stop("p_A.size() != r");
  }
  
  if (n.size() != r) {
    Rcpp::stop("n.size() != r");
  }
  
  double r_dbl = (double)r;
  
  double n_mean = 0.0;
  double n_sum = 0.0;
  double n2_sum = 0.0;
  
  for (int i = 0; i < r; ++i) {
    // Different from subpop.size()!
    double n_i = n[i];
    n_mean += n_i / r_dbl;
    n_sum += n_i;
    n2_sum += n_i * n_i;
  }

  // So have a common container with alleles to iterate over later
  std::unordered_set<int> alleles;
  
  // First, be sure that all alleles are known:
  for (int i = 0; i < r; ++i) {
    for (auto ele = p_A[i].begin(); ele != p_A[i].end(); ++ele) { 
      alleles.insert(ele->first);
    }
  }
  
  //********************************
  // GDA2, p. 168, p_A. tilde
  //********************************
  std::unordered_map<int, double> mean_pA;
  for (int i = 0; i < r; ++i) {
    for (auto ele = alleles.begin(); ele != alleles.end(); ++ele) {
      int allele = *ele;
      double p = p_A[i][allele]; // 0 if not exists
      mean_pA[allele] += (n[i] * p) / n_sum;
    } 
  }  

  //********************************
  // GDA2, p. 173, s^2
  //********************************
  std::unordered_map<int, double> s2_A;
  for (int i = 0; i < r; ++i) {
    for (auto ele = alleles.begin(); ele != alleles.end(); ++ele) {
      int allele = *ele;
      double p = p_A[i][allele]; // 0 if not exists
      double d = p - mean_pA[allele];
      s2_A[allele] += (n[i] * d * d) / ((r_dbl - 1.0) * n_mean);
    } 
  }  

  
  ////////////////////////////////////////////////////
  // Calculating helper variables
  ////////////////////////////////////////////////////
  
  //********************************
  // GDA2, p. 178, H_A. tilde
  //********************************
  std::unordered_map<int, double> mean_H_A;
  for (auto ele = alleles.begin(); ele != alleles.end(); ++ele) {
    int allele = *ele;
    
    for (int i = 0; i < r; ++i) {
      double pA = p_A[i][allele];
      double PAA = P_AA[i][allele];
      mean_H_A[allele] += (2.0 * n[i] * (pA - PAA)) / n_sum;
    }
  }
  
  ////////////////////////////////////////////////////
  // Calculating S1, S2, S3, GDA2, p. 178-179
  ////////////////////////////////////////////////////
  
  double nc = (n_sum - (n2_sum/n_sum)) / (r_dbl - 1.0);
  
  std::unordered_map<int, double> allele_S1;
  std::unordered_map<int, double> allele_S2;
  std::unordered_map<int, double> allele_S3;
  
  std::unordered_map<int, double> allele_MSG;
  std::unordered_map<int, double> allele_MSI;
  std::unordered_map<int, double> allele_MSP;
  std::unordered_map<int, double> allele_sigmasq_G;
  std::unordered_map<int, double> allele_sigmasq_I;
  std::unordered_map<int, double> allele_sigmasq_P;
  
  for (auto ele = alleles.begin(); ele != alleles.end(); ++ele) {
    int allele = *ele;
    double tmp_s2 = s2_A[allele];
    double tmp_p = mean_pA[allele];
    double tmp_HA = mean_H_A[allele];
    
    allele_S1[allele] = tmp_s2 - (1.0 / (n_mean - 1.0)) * (tmp_p * (1.0 - tmp_p) - ((r_dbl - 1.0) / r_dbl) * tmp_s2 - 0.25*tmp_HA);
    
    double tmp_S2_p1 = (r_dbl*(n_mean - nc)/n_mean) * tmp_p * (1.0 - tmp_p);
    double tmp_S2_p2 = tmp_s2 * ( (n_mean - 1.0) + (r_dbl - 1.0)*(n_mean - nc) ) / n_mean;
    double tmp_S2_p3 = tmp_HA * r_dbl * (n_mean - nc) / (4.0 * n_mean * nc);
    allele_S2[allele] = (tmp_p * (1.0 - tmp_p)) - (n_mean / (r_dbl * (n_mean - 1.0))) * (tmp_S2_p1 - tmp_S2_p2 - tmp_S2_p3);
    allele_S3[allele] = (nc / (2.0 * n_mean)) * tmp_HA;
    
    
    // GDA2 p. 177, Tab. 5.4
    double tmp_SSP = 2.0*(r_dbl - 1.0)*n_mean*tmp_s2;
    double tmp_SSI = 2.0*r_dbl*n_mean*tmp_p*(1.0-tmp_p) - 0.5*r_dbl*n_mean*tmp_HA - tmp_SSP;
    double tmp_SSG = 0.5*r_dbl*n_mean*tmp_HA;
    
    // Divide sums of squares with the d.f.:
    double tmp_MSP = tmp_SSP / (r_dbl - 1.0);
    double tmp_MSI = tmp_SSI / (n_sum - r_dbl);
    double tmp_MSG = tmp_SSG / n_sum;
    
    allele_MSG[allele] = tmp_MSG;
    allele_MSI[allele] = tmp_MSI;
    allele_MSP[allele] = tmp_MSP;
    
    double tmp_sigmasq_G = tmp_MSG;
    double tmp_sigmasq_I = 0.5*(tmp_MSI - tmp_MSG);
    double tmp_sigmasq_P = (tmp_MSP - tmp_MSI) / (2.0 * nc);
    
    allele_sigmasq_G[allele] = tmp_sigmasq_G;
    allele_sigmasq_I[allele] = tmp_sigmasq_I;
    allele_sigmasq_P[allele] = tmp_sigmasq_P;
  }  
  
  ////////////////////////////////////////////////////
  // Calculating S1, S2, S3, GDA2, p. 178-179
  ////////////////////////////////////////////////////
  double sum_S1 = 0.0;
  double sum_S2 = 0.0;
  double sum_S3 = 0.0;
  
  for (auto ele = alleles.begin(); ele != alleles.end(); ++ele) {
    int allele = *ele;
    
    sum_S1 += allele_S1[allele];
    sum_S2 += allele_S2[allele];
    sum_S3 += allele_S3[allele];
  }
  
  double F = 1 - sum_S3/sum_S2;
  double theta = sum_S1 / sum_S2;
  double f = (F - theta) / (1 - theta);
  
  /*
   F: Wright's F_{IT}: 
   Overall inbreeding coefficient; correlation of alleles within individuals over all populations
   
   theta: Coancestry, Wright's F_{ST}:
   Correlation of alleles of different individuals in the same poulation.
   
   f: Wright's F_{IS}
   Correlation of alleles within individuals within one populatoin.
   */

  Rcpp::List res_allele;
  for (auto ele = alleles.begin(); ele != alleles.end(); ++ele) {
    int allele = *ele;
    std::string allele_str = std::to_string (allele);
    
    Rcpp::List res_tmp;
    res_tmp["allele"] = allele;
    
    res_tmp["s2"] = s2_A[allele];
    res_tmp["p_A_mean"] = mean_pA[allele];
    res_tmp["H_A_mean"] = mean_H_A[allele];

    res_tmp["MSG"] = allele_MSG[allele];
    res_tmp["MSI"] = allele_MSI[allele];
    res_tmp["MSP"] = allele_MSP[allele];
    res_tmp["sigmasq_G"] = allele_sigmasq_G[allele];
    res_tmp["sigmasq_I"] = allele_sigmasq_I[allele];
    res_tmp["sigmasq_P"] = allele_sigmasq_P[allele];
    res_tmp["S1"] = allele_S1[allele];
    res_tmp["S2"] = allele_S2[allele];
    res_tmp["S3"] = allele_S3[allele];
    
    res_tmp["theta_sigmasq"] = allele_sigmasq_P[allele] / (allele_sigmasq_P[allele] + allele_sigmasq_I[allele] + allele_sigmasq_G[allele]);
    res_tmp["theta_S1S2"] = allele_S1[allele] / allele_S2[allele];
    
    res_allele[allele_str] = res_tmp;
  }

  Rcpp::List res_additional;
  res_additional["S1"] = sum_S1;
  res_additional["S2"] = sum_S2;
  res_additional["S3"] = sum_S3;
  res_additional["n_c"] = nc;
  res_additional["n_mean"] = n_mean;
  res_additional["n_sum"] = n_sum;
  res_additional["nsq_sum"] = n2_sum;
  res_additional["allele_values"] = res_allele;
  
  Rcpp::List res;
  res["F"] = F;
  res["theta"] = theta;
  res["f"] = f;
  res["extra"] = res_additional;

  return res;
}




void fill_P_AA_p_A(int a, int b, int i, double frac1, double frac2, 
                       std::vector< std::unordered_map<int, double> >& P_AA,
                       std::vector< std::unordered_map<int, double> >& p_A
) {
  
  /*
   * 
   * GDA2, p. 18:
   * 
   * p_A[i][a]: allele frequency of allele a in subpopulation i
   * P_AA[i][a]: homozygous frequency of genotype AA in subpopulation i
   * 
   * double frac1 = 1.0 / (2.0 * sample_size_i);
   * double frac2 = 1.0 / sample_size_i;
   */
  if (a == b) {
    // Homozygous
    
    // There are 2*n alleles, and two of them are now a:
    // p_A[i][a]: there are two alleles, so 
    // 2*1.0 / (2.0 * sample_size_i) = 2*frac1 = frac2 = 1.0 / sample_size_i.
    p_A[i][a] += frac2;
    
    // There are n genotypes in the population, and 1 of them is now homozygous
    P_AA[i][a] += frac2;
  } else {
    // Heterozygous
    
    // There are 2*n alleles, and 1 is a and another is b
    p_A[i][a] += frac1;
    p_A[i][b] += frac1;
  }
}


//' Estimate F, theta, and f from subpopulations of individuals
//' 
//' Estimates F, theta, and f for a number of subpopulations given a list of individuals.
//' 
//' Based on Bruce S Weir, Genetic Data Analysis 2, 1996. (GDA2).
//' 
//' @param subpops List of subpopulations, each a list of individuals
//' @param subpops_sizes Size of each subpopulation
//' 
//' @return Estimates of F, theta, and f as well as additional information
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_theta_subpops_individuals(Rcpp::List subpops, 
                                              Rcpp::IntegerVector subpops_sizes) {

  int r = subpops.size();
  
  if (r <= 0) {
    Rcpp::stop("No subpopulations given");
  }
  
  if (subpops_sizes.size() != r) {
    Rcpp::stop("length(subpops) != length(subpops_sizes)");
  }
  
  if (any(subpops_sizes <= 0).is_true()) {
    Rcpp::stop("All subpops_sizes must be positive");
  }
  
  // H_A: Heterozygous probabilities:
  // H_A[i][l]: Heterozygous probability of allele l in subpopulation i.
  std::vector< std::unordered_map<int, double> > H_A(r);

  // P_AA: Homozygous probabilities:
  // P_AA[i][l]: Homozygous probability of allele l in subpopulation i.
  std::vector< std::unordered_map<int, double> > P_AA(r);

  // p_A: Allele probability:
  // p_A[i][l]: Allele probability of allele l in subpopulation i.
  std::vector< std::unordered_map<int, double> > p_A(r);

  std::vector<double> n(r);
  
  ////////////////////////////////////////////////////
  // Filling containers with probabilities  
  ////////////////////////////////////////////////////
  for (int i = 0; i < r; ++i) {
    Rcpp::List subpop = Rcpp::as< Rcpp::List >(subpops[i]);
    
    if (subpop.size() <= 0) {
      Rcpp::stop("Subpop sample of size <= 0");
    }
    
    if (subpops_sizes[i] <= 0) {
      Rcpp::stop("Subpop size <= 0");
    }

    n[i] = subpops_sizes[i];
    
    double sample_size_i = (double)subpop.size();
    double frac1 = 1.0 / (2.0 * sample_size_i);
    double frac2 = 1.0 / sample_size_i;
    
    for (int j = 0; j < sample_size_i; ++j) {
      Rcpp::XPtr<Individual> individual = Rcpp::as< Rcpp::XPtr<Individual> >(subpop[j]);

      if (!(individual->is_haplotype_set())) {
        Rcpp::stop("Haplotypes not yet set");
      }      
      
      std::vector<int> hap = individual->get_haplotype();
      if (hap.size() != 2) {
        Rcpp::stop("Expected exactly 2 autosomal loci");
      }
      
      fill_P_AA_p_A(hap[0], hap[1], i, frac1, frac2, P_AA, p_A);    
    }
  }
  
  Rcpp::List res = estimate_theta_subpops_weighted_engine(P_AA, p_A, n);

  return res;
}

//' Estimate F, theta, and f from subpopulations of genotypes
//' 
//' Estimates F, theta, and f for a number of subpopulations given a list of genotypes.
//' 
//' Based on Bruce S Weir, Genetic Data Analysis 2, 1996. (GDA2).
//' 
//' @param subpops List of subpopulations, each a list of individuals
//' @param subpops_sizes Size of each subpopulation
//' 
//' @return Estimates of F, theta, and f as well as additional information
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_theta_subpops_genotypes(Rcpp::ListOf<Rcpp::IntegerMatrix> subpops, 
                                            Rcpp::IntegerVector subpops_sizes) {
  
  int r = subpops.size();
  
  if (r <= 0) {
    Rcpp::stop("No subpopulations given");
  }
  
  if (subpops_sizes.size() != r) {
    Rcpp::stop("length(subpops) != length(subpops_sizes)");
  }
  
  if (any(subpops_sizes <= 0).is_true()) {
    Rcpp::stop("All subpops_sizes must be positive");
  }
  
  
  // H_A: Heterozygous probabilities:
  // H_A[i][l]: Heterozygous probability of allele l in subpopulation i.
  std::vector< std::unordered_map<int, double> > H_A(r);
  
  // P_AA: Homozygous probabilities:
  // P_AA[i][l]: Homozygous probability of allele l in subpopulation i.
  std::vector< std::unordered_map<int, double> > P_AA(r);
  
  // p_A: Allele probability:
  // p_A[i][l]: Allele probability of allele l in subpopulation i.
  std::vector< std::unordered_map<int, double> > p_A(r);
  
  std::vector<double> n(r);
  
  ////////////////////////////////////////////////////
  // Filling containers with probabilities  
  ////////////////////////////////////////////////////
  for (int i = 0; i < r; ++i) {
    Rcpp::IntegerMatrix subpop = subpops[i];
    
    if (subpop.nrow() <= 0) {
      Rcpp::stop("Subpop sample of size <= 0");
    }
    
    if (subpop.ncol() != 2) {
      Rcpp::stop("Expected exactly 2 autosomal loci");
    }
    
    if (subpops_sizes[i] <= 0) {
      Rcpp::stop("Subpop size <= 0");
    }
    
    n[i] = subpops_sizes[i];
    
    double sample_size_i = (double)subpop.nrow();
    double frac1 = 1.0 / (2.0 * sample_size_i);
    double frac2 = 1.0 / sample_size_i;
    
    for (int j = 0; j < sample_size_i; ++j) {
      Rcpp::IntegerVector hap = subpop(j, Rcpp::_);
      fill_P_AA_p_A(hap[0], hap[1], i, frac1, frac2, P_AA, p_A);
    }
  }
  
  /*
  print_container("Heterozygous", H_A);
  print_container("Homozygous", P_AA);
  print_container("Allele", p_A);
  Rcpp::Rcout << "n:\n";
  Rcpp::print(Rcpp::wrap(n));
  */
  
  Rcpp::List res = estimate_theta_subpops_weighted_engine(P_AA, p_A, n);
  
  return res;
}


//' Estimate F, theta, and f from subpopulations of individual ids
//' 
//' Estimates F, theta, and f for a number of subpopulations given a list of pids (individual ids).
//' 
//' Based on Bruce S Weir, Genetic Data Analysis 2, 1996. (GDA2).
//' 
//' @param population Population obtain from simulation
//' @param subpops List of individual pids
//' @param subpops_sizes Size of each subpopulation
//' 
//' @return Estimates of F, theta, and f as well as additional information
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List estimate_theta_subpops_pids(Rcpp::XPtr<Population> population,
                                       Rcpp::ListOf<Rcpp::IntegerVector> subpops, 
                                       Rcpp::IntegerVector subpops_sizes) {
  
  int r = subpops.size();

  if (r <= 0) {
    Rcpp::stop("No subpopulations given");
  }
  
  if (subpops_sizes.size() != r) {
    Rcpp::stop("length(subpops) != length(subpops_sizes)");
  }
  
  if (any(subpops_sizes <= 0).is_true()) {
    Rcpp::stop("All subpops_sizes must be positive");
  }
  
  // P_AA: Homozygous probabilities:
  // P_AA[i][l]: Homozygous probability of allele l in subpopulation i.
  std::vector< std::unordered_map<int, double> > P_AA(r);
  
  // p_A: Allele probability:
  // p_A[i][l]: Allele probability of allele l in subpopulation i.
  std::vector< std::unordered_map<int, double> > p_A(r);
  
  std::vector<double> n(r);
  
  ////////////////////////////////////////////////////
  // Filling containers with probabilities  
  ////////////////////////////////////////////////////
  for (int i = 0; i < r; ++i) {
    Rcpp::IntegerVector subpop_pids = subpops[i];

    if (subpop_pids.size() <= 0) {
      Rcpp::stop("Subpop sample of size <= 0");
    }
    
    if (subpops_sizes[i] <= 0) {
      Rcpp::stop("Subpop size <= 0");
    }
    
    n[i] = subpops_sizes[i];
    
    double sample_size_i = (double)subpop_pids.size();
    double frac1 = 1.0 / (2.0 * sample_size_i);
    double frac2 = 1.0 / sample_size_i;
    
    for (int j = 0; j < sample_size_i; ++j) {
      int pid = subpop_pids[j];
      Individual* individual = population->get_individual(pid);
      
      if (!(individual->is_haplotype_set())) {
        Rcpp::stop("Haplotypes not yet set");
      }      
      
      std::vector<int> hap = individual->get_haplotype();
      if (hap.size() != 2) {
        Rcpp::stop("Expected exactly 2 autosomal loci");
      }
      
      fill_P_AA_p_A(hap[0], hap[1], i, frac1, frac2, P_AA, p_A);    
    }
  }
  
  Rcpp::List res = estimate_theta_subpops_weighted_engine(P_AA, p_A, n);
  
  return res;
}











