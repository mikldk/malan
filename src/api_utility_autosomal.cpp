/**
 api_utility_haplotypes.cpp
 Purpose: Logic related to haplotypes.
 Details: API between R user and C++ logic.
  
 @author Mikkel Meyer Andersen
 */

//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
//#include <Rcpp.h>


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
                                      
//' Populate 1-locus autosomal DNA profile in pedigrees with single-step mutation model.
//' 
//' Populate 1-locus autosomal DNA profile from founder and down in all pedigrees.
//' Note, that only alleles from ladder is assigned and 
//' that all founders draw type randomly.
//' 
//' Note, that pedigrees must first have been inferred by [build_pedigrees()].
//' 
//' @param pedigrees Pedigree list in which to populate genotypes
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



//' Populate 1-locus autosomal DNA profile in pedigrees with infinite alleles mutation model.
//' 
//' Populate 1-locus autosomal DNA profile from founder and down in all pedigrees.
//' Note, that all founders have type 0 to begin with.
//' 
//' The maternal allele is taken by random from 
//' the `2*N[g]` alleles in the previous generation consisting of `N[g]` males
//' with descendants in the live population.
//' 
//' This is also why this is not using pedigrees but instead the population.
//' 
//' Note, that pedigrees need not be inferred.
//' 
//' @param population Population in which to populate genotypes
//' @param mutation_rate Mutation rate between 0 and 1 (both included)
//' @param progress Show progress
//'
//' @seealso [pedigrees_all_populate_haplotypes_custom_founders()] and 
//' [pedigrees_all_populate_haplotypes_ladder_bounded()].
//' 
//' @export
// [[Rcpp::export]]
void population_populate_autosomal_infinite_alleles(
    Rcpp::XPtr<Population> population, 
    double mutation_rate,
    bool progress = true) {  
  
  // FIXME: Maternal alleles: Only from individuals with descendants...
  
  // Put individuals into generations
  std::unordered_map<int, Individual*> population_map = *(population->get_population());
  
  int max_generation = -1;

  for (std::pair<int, Individual*> element : population_map) {
    Individual* ind = element.second;
    int gen = ind->get_generation();
    
    if (gen > max_generation) {
      max_generation = gen;
    }
  }
  
  if (max_generation == -1) {
    Rcpp::stop("Unexpected error");
  }
  
  std::vector< std::vector<Individual*> > generation_individuals(max_generation + 1);
  
  for (std::pair<int, Individual*> element : population_map) {
    Individual* ind = element.second;
    int gen = ind->get_generation();
    generation_individuals[gen].push_back(ind);
  }
  
  //////////////////////////////////////////////
  // Then iterate over generations
  //////////////////////////////////////////////
  
  //--------------------------------------------
  // First generation.
  //--------------------------------------------
  // Note, that last generation is 0, so max_generation is 
  // actually the first generation.
  int current_allele = 0;
  
  std::unordered_map<int, int> allele_counts;
  allele_counts[current_allele] = 1;
  
  //Rcpp::print(Rcpp::wrap(max_generation));
  for (Individual* ind : generation_individuals[max_generation]) {
    std::vector<int> h = {current_allele, current_allele};
    ind->set_haplotype(h);
  }
  
  if (max_generation == 0) {
    return;
  }
  
  //--------------------------------------------
  // Second generation and forward
  //--------------------------------------------
  for (int i = max_generation - 1; i >= 0; --i) {
    //Rcpp::print(Rcpp::wrap(i));
    
    // Previous generation's alleles
    int alleles_count = 0;
    for (std::pair<int, int> p : allele_counts) {
      alleles_count += p.second;
    }
    double alleles_count_dbl = (double)alleles_count;
    
    double p_sum = 0.0;
    std::vector<int> alleles_support;
    std::vector<double> alleles_cumdist;
    for (std::pair<int, int> p : allele_counts) {
      p_sum += p.second / alleles_count_dbl;
      alleles_support.push_back(p.first);
      alleles_cumdist.push_back(p_sum);
    }
    
    //if (i >= max_generation - 2) {
    if (!true) {
      Rcpp::Rcout << "===========" << std::endl;
      Rcpp::Rcout << i << std::endl;
      Rcpp::Rcout << "===========" << std::endl;
      Rcpp::print(Rcpp::wrap(allele_counts));
      Rcpp::print(Rcpp::wrap(alleles_cumdist));
      Rcpp::print(Rcpp::wrap(alleles_support));
    }
    
    allele_counts.clear();
    
    for (Individual* ind : generation_individuals[i]) {
      int father_allele = -1;
      int mother_allele = -1;
      
      // Mutation father_allele; if so, mutate to new kind, 
      // doesn't depend on father's actual allele!
      if (R::runif(0.0, 1.0) < mutation_rate) {
        // Mutation happened
        if ((current_allele + 1) < current_allele) {
          Rcpp::stop("Overflow for current_allele happened!");
        }
        current_allele += 1;
        father_allele = current_allele;
      } else {
        std::vector<int> geno_father = ind->get_father()->get_haplotype();
        father_allele = (R::runif(0.0, 1.0) < 0.5) ? geno_father[0] : geno_father[1];
      }

      // Mutations mother_allele; if so, mutate to new kind, 
      // doesn't depend on mothers's actual allele!
      if (R::runif(0.0, 1.0) < mutation_rate) {
        // Mutation happened
        if ((current_allele + 1) < current_allele) {
          Rcpp::stop("Overflow for current_allele happened!");
        }
        current_allele += 1;
        mother_allele = current_allele;
      } else {
        const double u_mother_allele = R::runif(0.0, 1.0); // in (0, 1), both 0/1 excluded
        mother_allele = alleles_support[0];
        
        if (u_mother_allele > alleles_cumdist[0]) {
          for (int i = 1; i < alleles_support.size(); ++i) {
            if (u_mother_allele <= alleles_cumdist[i]) {
              mother_allele = alleles_support[i];
              break;
            }
          }
        }
      }
      
      std::vector<int> geno = { father_allele, mother_allele };
      
      if (geno[1] <= geno[0]) {
        int tmp = geno[0];
        geno[0] = geno[1];
        geno[1] = tmp;
      }

      ind->set_haplotype(geno);
      
      allele_counts[father_allele] += 1;
      allele_counts[mother_allele] += 1;
    }
  }
}

