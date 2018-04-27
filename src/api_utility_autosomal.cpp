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

//' Calculate genotype probabilities with theta
//' 
//' @param allele_dist Allele distribution (probabilities) -- gets normalised
//' @param theta Theta correction between 0 and 1 (both included)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> calc_autosomal_genotype_probs(Rcpp::NumericVector allele_dist,
                                                  double theta) {
  
  if (any(allele_dist < 0).is_true() || any(allele_dist > 1).is_true()) {
    Rcpp::stop("allele_dist's elements must be between 0 and 1, both included");
  }
  
  if (theta < 0 || theta > 1) {
    Rcpp::stop("theta must be between 0 and 1, both included");
  }
  
  std::vector<double> ps = Rcpp::as< std::vector<double> >(allele_dist);
  double ps_sum = std::accumulate(ps.begin(), ps.end(), 0.0);
  const int alleles_count = ps.size();          
  
  // Normalisation
  for (int i = 0; i < alleles_count; ++i) {
    ps[i] = ps[i] / ps_sum;
  }
  
  std::vector<double> allele_dist_theta(alleles_count * (alleles_count + 1) / 2);
  int k = 0;
                                    
  for (int i = 0; i < alleles_count; ++i) {
    for (int j = 0; j <= i; ++j) {   
      if (i == j) { // homozyg
        allele_dist_theta[k] = theta*ps[i] + (1.0-theta)*ps[i]*ps[i];        
      } else { // hetegozyg
        allele_dist_theta[k] = (1.0-theta)*2.0*ps[i]*ps[j];
      }

      k++;
    }
  }
  
  return allele_dist_theta;                                   
}

//' Calculate conditional genotype cumulative probabilities with theta
//' 
//' @param allele_dist Allele distribution (probabilities) -- gets normalised
//' @param theta Theta correction between 0 and 1 (both included)
//' 
//' @return Matrix: row i: conditional cumulative distribution of alleles given allele i
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix calc_autosomal_genotype_conditional_cumdist(
    Rcpp::NumericVector allele_dist,
    double theta) {
  
  if (any(allele_dist < 0).is_true() || any(allele_dist > 1).is_true()) {
    Rcpp::stop("allele_dist's elements must be between 0 and 1, both included");
  }
  
  if (theta < 0 || theta > 1) {
    Rcpp::stop("theta must be between 0 and 1, both included");
  }
  
  std::vector<double> ps = Rcpp::as< std::vector<double> >(allele_dist);
  double ps_sum = std::accumulate(ps.begin(), ps.end(), 0.0);
  const int alleles_count = ps.size();          
  
  // Normalisation
  for (int i = 0; i < alleles_count; ++i) {
    ps[i] = ps[i] / ps_sum;
  }
  
  Rcpp::NumericMatrix dists(alleles_count, alleles_count);

  for (int i = 0; i < alleles_count; ++i) {
    for (int j = 0; j <= i; ++j) {
      if (i == j) { // homozyg
        double p = theta*ps[i] + (1.0-theta)*ps[i]*ps[i];
        dists(i, i) = p;
      } else { // hetegozyg
        double p = (1.0-theta)*ps[i]*ps[j];
        dists(i, j) = p;
        dists(j, i) = p;
      }
    }
  }
  
  // Multiple passes, but easier to follow:
  
  // Get row i to sum to 1; ps[i] = sum(dists(i, Rcpp::_))
  for (int i = 0; i < alleles_count; ++i) {
    Rcpp::NumericVector row = dists(i, Rcpp::_) / ps[i];
    
    // cumsum, for some reason Rcpp::cumsum doesn't work...
    Rcpp::NumericVector res(row.size());
    std::partial_sum(row.begin(), row.end(), res.begin(), std::plus<double>());
    dists(i, Rcpp::_) = res;
  }
  
  return dists;                                   
}
             
//' Sample genotype with theta
//' 
//' @param allele_dist Allele distribution (probabilities) -- gets normalised
//' @param theta Theta correction between 0 and 1 (both included)
//'
//' @export
// [[Rcpp::export]]
std::vector<int> sample_autosomal_genotype(Rcpp::NumericVector allele_dist,
                                           double theta) {
                                           
  const int alleles_count = allele_dist.size();
  const std::vector<double> allele_dist_theta = calc_autosomal_genotype_probs(allele_dist, theta);
  
  std::vector<double> allele_cumdist_theta(allele_dist_theta.size());
  std::partial_sum(allele_dist_theta.begin(), allele_dist_theta.end(), allele_cumdist_theta.begin(), std::plus<double>());
  
  std::vector<int> geno = draw_autosomal_genotype(allele_cumdist_theta, alleles_count);
  
  return geno;
}
                                      
//' Populate 1-locus autosomal DNA profile in pedigrees.
//' 
//' Populate 1-locus autosomal DNA profile from founder and down in all pedigrees.
//' Note, that only alleles from ladder is assigned and 
//' that all founders draw type randomly.
//' 
//' Note, that pedigrees must first have been inferred by [build_pedigrees()].
//' 
//' @param pedigrees Pedigree list in which to populate haplotypes
//' @param allele_dist Allele distribution (probabilities) -- gets normalised
//' @param theta Theta correction between 0 and 1 (both included)
//' @param mutation_rate Mutation rate between 0 and 1 (both included)
//' @param progress Show progress
//'
//' @seealso [pedigrees_all_populate_haplotypes_custom_founders()] and 
//' [pedigrees_all_populate_haplotypes_ladder_bounded()].
//' 
//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_autosomal(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, 
                                      Rcpp::NumericVector allele_dist,
                                      double theta,
                                      double mutation_rate,
                                      bool progress = true) {  
  std::vector<Pedigree*> peds = (*pedigrees);

  // For drawing founder types ->
  const int alleles_count = allele_dist.size();
  const std::vector<double> allele_dist_theta = calc_autosomal_genotype_probs(allele_dist, theta);
  std::vector<double> allele_cumdist_theta(allele_dist_theta.size());
  std::partial_sum(allele_dist_theta.begin(), allele_dist_theta.end(), allele_cumdist_theta.begin(), std::plus<double>());
  // <- founder
  
  // For children's ->
  Rcpp::NumericMatrix cumdist_mat = calc_autosomal_genotype_conditional_cumdist(allele_dist, theta);
  
  if (cumdist_mat.nrow() != alleles_count) {
    Rcpp::stop("Unexpected error");
  }
  std::vector< std::vector<double> > cumdists(alleles_count);
  
  for (int i = 0; i < alleles_count; ++i) {
    Rcpp::NumericVector row_vec = cumdist_mat(i, Rcpp::_);
    std::vector<double> row = Rcpp::as< std::vector<double> >(row_vec);
    cumdists[i] = row;
  }
  // <- 
  
  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_autosomal(cumdists, allele_cumdist_theta, alleles_count, mutation_rate);
    
    if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      Rcpp::stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}



// boost::hash_combine
// https://stackoverflow.com/a/27952689/3446913
size_t hash_combine(size_t lhs, size_t rhs) {
  lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
  return lhs;
}

// https://stackoverflow.com/a/20602159/3446913
struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return hash_combine(x.first, x.second);
  }
};

// [[Rcpp::export]]
std::unordered_map<int, int> hash_colisions(int p) {
  std::unordered_map<int, int> tab;
  
  for (int i = 0; i < (p-1); ++i) {
    for (int j = (i+1); j < p; ++j) {
      int hash = hash_combine(i, j);
      tab[hash] += 1;
    }
  }
  
  return tab;
}




Rcpp::List estimate_theta_1subpop(const std::unordered_map<int, double>& allele_p,
                                  const std::unordered_map<std::pair<int, int>, double, pairhash>& genotype_p,
                                  const std::unordered_set<std::pair<int, int>, pairhash>& genotypes_unique,
                                  const bool return_estimation_info = false) {
  Rcpp::List theta;
  
  theta["estimate"] = NA_REAL;
  theta["error"] = true;
  theta["details"] = "NA";
  theta["estimation_info"] = R_NilValue;  

  // Loop over unique genotypes
  std::unordered_set<std::pair<int, int>, pairhash>::const_iterator it;
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
                                            std::unordered_map<std::pair<int, int>, double, pairhash>& genotype_p,
                                            std::unordered_set<std::pair<int, int>, pairhash>& genotypes_unique) {
  
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
  std::unordered_map<std::pair<int, int>, double, pairhash> genotype_p;
  std::unordered_set<std::pair<int, int>, pairhash> genotypes_unique;
  
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
  std::unordered_map<std::pair<int, int>, double, pairhash> genotype_p;
  std::unordered_set<std::pair<int, int>, pairhash> genotypes_unique;
  
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
// H_A: Heterozygous probabilities:
// H_A[i][l]: Heterozygous probability of allele l in subpopulation i.
std::vector< std::unordered_map<int, double> > H_A(r);

// P_AA: Homozygous probabilities:
// P_AA[i][l]: Homozygous probability of allele l in subpopulation i.
std::vector< std::unordered_map<int, double> > P_AA(r);

// p_A: Allele probability:
// p_A[i][l]: Allele probability of allele l in subpopulation i.
*/
Rcpp::List estimate_theta_subpops_weighted_engine(
    std::vector< std::unordered_map<int, double> > H_A,
    std::vector< std::unordered_map<int, double> > P_AA,
    std::vector< std::unordered_map<int, double> > p_A,
    std::vector<double> n) {
  
  int r = H_A.size();
  
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
  
  double n_mean = 0.0;
  double n_sum = 0.0;
  double n2_sum = 0.0;
  
  for (int i = 0; i < r; ++i) {
    // Different from subpop.size()!
    double n_i = n[i];
    n_mean += n_i / (double)r;
    n_sum += n_i;
    n2_sum += n_i * n_i;
  }

  ////////////////////////////////////////////////////
  // Calculating helper variables
  ////////////////////////////////////////////////////
  
  //********************************
  // GDA2, p. 178, H_A. tilde
  //********************************
  std::unordered_map<int, double> mean_H_A;
  for (int i = 0; i < r; ++i) {
    for (auto ele = H_A[i].begin(); ele != H_A[i].end(); ++ele) {
      int allele = ele->first;      
      double HAi = ele->second;
      
      mean_H_A[allele] += (n[i] * HAi) / n_sum;
    } 
  }  
  //Rcpp::Rcout << "mean_H_A\n";
  //print_map(mean_H_A);
  
  /*
  // Gives the same as mean_H_A
  std::unordered_map<int, double> mean_H_A_type2;
  for (int i = 0; i < r; ++i) {
  for (auto ele = H_A[i].begin(); ele != H_A[i].end(); ++ele) {
  int allele = ele->first;      
  double HAi = ele->second;
  alleles.insert(allele);
  
  mean_H_A_type2[allele] += (2.0 * n[i] * (p_A[i][allele] - P_AA[i][allele]) ) / n_sum;
  } 
  }  
  Rcpp::Rcout << "mean_H_A_type2\n";
  print_map(mean_H_A_type2);
  */
  
  // So have a common container with alleles to iterate over later
  std::unordered_set<int> alleles;
  
  //********************************
  // GDA2, p. 168, p_A. tilde
  //********************************
  std::unordered_map<int, double> mean_pA;
  for (int i = 0; i < r; ++i) {
    for (auto ele = p_A[i].begin(); ele != p_A[i].end(); ++ele) { 
      alleles.insert(ele->first);
      mean_pA[ ele->first ] += (n[i] * ele->second) / n_sum;
    } 
  }  
  //Rcpp::Rcout << "mean_pA\n";
  //print_map(mean_pA);
  
  //********************************
  // GDA2, p. 173, s^2
  //********************************
  std::unordered_map<int, double> s2_A;
  for (int i = 0; i < r; ++i) {
    for (auto ele = p_A[i].begin(); ele != p_A[i].end(); ++ele) {
      alleles.insert(ele->first);
      double d = ele->second - mean_pA[ele->first];
      s2_A[ ele->first ] += (n[i] * d * d) / ((r-1.0) * n_mean);
    } 
  }  
  //Rcpp::Rcout << "s2_A\n";
  //print_map(s2_A);
  
  ////////////////////////////////////////////////////
  // Calculating S1, S2, S3, GDA2, p. 178-179
  ////////////////////////////////////////////////////
  
  double nc = (n_sum - n2_sum/n_sum) / (double)(r - 1);
  //Rcpp::Rcout << "n_c = " << nc << "\n";
  
  double r_dbl = (double)r;
  
  std::unordered_map<int, double> allele_S1;
  std::unordered_map<int, double> allele_S2;
  std::unordered_map<int, double> allele_S3;
  
  for (auto ele = alleles.begin(); ele != alleles.end(); ++ele) {
    int allele = *ele;
    double tmp_s2 = s2_A[allele];
    double tmp_p = mean_pA[allele];
    double tmp_HA = mean_H_A[allele];
    
    allele_S1[allele] = tmp_s2 - (1.0 / (n_mean - 1.0)) * (tmp_p * (1.0 - tmp_p) - ((r_dbl - 1.0) / r_dbl) * tmp_s2 - 0.25*tmp_HA);
    
    double tmp_S2_p1 = (r_dbl*(n_mean - nc)/n_mean) * tmp_p * (1.0 - tmp_p);
    double tmp_S2_p2 = tmp_s2 * ( (n_mean - 1) + (r_dbl - 1)*(n_mean - nc) ) / n_mean;
    double tmp_S2_p3 = tmp_HA * r_dbl * (n_mean - nc) / (4.0 * n_mean * nc);
    allele_S2[allele] = (tmp_p * (1.0 - tmp_p)) - (n_mean / (r_dbl * (n_mean - 1.0))) * (tmp_S2_p1 - tmp_S2_p2 - tmp_S2_p3);
    
    allele_S3[allele] = (nc / (2.0 * n_mean)) * tmp_HA;
  }  
  /*
  Rcpp::Rcout << "allele_S1\n";
  print_map(allele_S1);
  
  Rcpp::Rcout << "allele_S2\n";
  print_map(allele_S2);
  
  Rcpp::Rcout << "allele_S3\n";
  print_map(allele_S3);  
  */
  
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
  
  /*
   Rcpp::Rcout << "F = " << F << "\n";
   Rcpp::Rcout << "theta = " << theta << "\n";
   Rcpp::Rcout << "f = " << f << "\n";
   */
  
  Rcpp::List res;
  res["F"] = F;
  res["theta"] = theta;
  res["f"] = f;
  
  return res;
}




void fill_H_A_P_AA_p_A(int a, int b, int i, double frac1, double frac2, 
                       std::vector< std::unordered_map<int, double> >& H_A,
                       std::vector< std::unordered_map<int, double> >& P_AA,
                       std::vector< std::unordered_map<int, double> >& p_A
) {
  
  if (a == b) {
    // Homozygous
    p_A[i][a] += frac2;
    
    P_AA[i][a] += frac2;
  } else {
    // Heterozygous
    p_A[i][a] += frac1;
    p_A[i][b] += frac1;
    
    H_A[i][a] += frac2;
    H_A[i][b] += frac2;
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
//' @return  Estimates of F, theta, and f
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
    
    double sample_size_i = subpop.size();
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
      
      fill_H_A_P_AA_p_A(hap[0], hap[1], i, frac1, frac2, H_A, P_AA, p_A);    
    }
  }
  
  /*  
  print_container("Heterozygous", H_A);
  print_container("Homozygous", P_AA);
  print_container("Allele", p_A);  
  Rcpp::Rcout << "n_mean = " << n_mean << std::endl;
  */
  
  Rcpp::List res = estimate_theta_subpops_weighted_engine(H_A, P_AA, p_A, n);

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
//' @return  Estimates of F, theta, and f
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
    
    double sample_size_i = subpop.nrow();
    double frac1 = 1.0 / (2.0 * sample_size_i);
    double frac2 = 1.0 / sample_size_i;
    
    for (int j = 0; j < sample_size_i; ++j) {
      Rcpp::IntegerVector hap = subpop(j, Rcpp::_);
      fill_H_A_P_AA_p_A(hap[0], hap[1], i, frac1, frac2, H_A, P_AA, p_A);
    }
  }
  
  /*  
   print_container("Heterozygous", H_A);
   print_container("Homozygous", P_AA);
   print_container("Allele", p_A);  
   Rcpp::Rcout << "n_mean = " << n_mean << std::endl;
   */
  
  Rcpp::List res = estimate_theta_subpops_weighted_engine(H_A, P_AA, p_A, n);
  
  return res;
}


