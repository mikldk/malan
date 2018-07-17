#include <unordered_set>

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

Rcpp::List get_mixture_analyse_failure(std::string fail_reason) {
  Rcpp::List terms_i;
  
  terms_i["fail"] = true;
  terms_i["fail_reason"] = fail_reason;
  terms_i["terms_Hp"] = Rcpp::IntegerVector::create(NA_INTEGER);
  terms_i["terms_Hp_count"] = Rcpp::IntegerVector::create(NA_INTEGER);
  terms_i["LR_contrib_Hp"] = Rcpp::IntegerVector::create(NA_INTEGER);
  terms_i["terms_Hd"] = Rcpp::IntegerVector::create(NA_INTEGER);
  terms_i["terms_Hd_count"] = Rcpp::IntegerVector::create(NA_INTEGER);
  terms_i["LR_contrib_Hd"] = Rcpp::IntegerVector::create(NA_INTEGER);
  terms_i["number_contributors"] = Rcpp::IntegerVector::create(NA_INTEGER);
  
  return terms_i;
}

/* 
  Given haplotypes by haps_in_mixture_indices and haps_in_mixture (with counts haps_in_mixture_counts), 
  do they together with known_contributors give the mixture?
  If so, add the counts of the haplotypes in haps_in_mixture_indices to 
  counts_contributor_sets.
*/
void analyse_set(const int num_contributors, 
                 const std::vector<int>& haps_in_mixture_indices, 
                 const std::vector< std::vector<int> >& haps_in_mixture,
                 const std::vector<int>& haps_in_mixture_counts, 
                 const std::vector<int>& poi_profile,
                 const std::vector< std::vector<int> >& known_contributors,
                 const std::vector< std::unordered_set<int> >& mixture,
                 std::vector< std::vector<int> >& counts_contributor_sets) {

  //Rcpp::Rcout << "Hello...." << std::endl;
  //Rcpp::print(Rcpp::wrap(haps_in_mixture_indices));  
  
  if ((haps_in_mixture_indices.size() + known_contributors.size()) != num_contributors) {
    Rcpp::stop("# unknown + # known != K");
  }
  
  if (haps_in_mixture.size() != haps_in_mixture_counts.size()) {
    Rcpp::stop("haps_in_mixture.size() != haps_in_mixture_counts.size()");
  }

  std::vector< std::vector<int> > suggested_contributors;
    
  for (int set_i = 0; set_i < haps_in_mixture_indices.size(); set_i++) {
    int hap_index = haps_in_mixture_indices[set_i];    
    const std::vector<int> h = haps_in_mixture[hap_index];
    //Rcpp::print(Rcpp::wrap(h));
    suggested_contributors.push_back(h);
  }

  const int loci = haps_in_mixture[0].size();
  
  bool set_gives_mixtures = true;
  
  //Rcpp::Rcout << "=========================" << std::endl;
  
  std::vector< std::unordered_set<int> > tmp(loci);
  
  for (int locus = 0; locus < loci; ++locus) {   
    std::unordered_set<int> set_alleles;
    
    // known:
    for (int set_i = 0; set_i < known_contributors.size(); ++set_i) {
      set_alleles.insert(known_contributors[set_i][locus]);
    }
        
    // unknown:
    for (int set_i = 0; set_i < suggested_contributors.size(); ++set_i) {
      set_alleles.insert(suggested_contributors[set_i][locus]);
    }
    
    tmp[locus] = set_alleles;
    
    //Rcpp::Rcout << "Test..." << std::endl;
    //Rcpp::print(Rcpp::wrap(set_alleles));
    //Rcpp::print(Rcpp::wrap(mixture[locus]));  
  
    if (mixture[locus] != set_alleles) {
      set_gives_mixtures = false;
      break;
    }
  }
  
  if (!set_gives_mixtures) {
    return;
  }
  
  // Gives mixture, fill counts:
  std::vector<int> counts;
  
  for (int level = 0; level < haps_in_mixture_indices.size(); level++) {
    int hap_index = haps_in_mixture_indices[level];    
    
    std::vector<int> rnd_profile = haps_in_mixture[hap_index];
    int n_u = haps_in_mixture_counts[hap_index];
    
    // Defence's hypothesis: 
    //   K random males different from the suspect, 
    //   hence ignore the suspect by substracting one
    if (poi_profile.size() > 0 && poi_profile == rnd_profile) {
      n_u = n_u - 1;
    }
    
    counts.push_back(n_u);
  }
  
  /*
  Rcpp::Rcout << "################################################" << std::endl;
  Rcpp::Rcout << "counts:" << std::endl;
  Rcpp::print(Rcpp::wrap(counts));  
  Rcpp::Rcout << "haps_in_mix_indices:" << std::endl;
  Rcpp::print(Rcpp::wrap(haps_in_mixture_indices));  
  for (int locus = 0; locus < loci; ++locus) {
    Rcpp::Rcout << "Mixture:" << std::endl;
    Rcpp::print(Rcpp::wrap(mixture[locus]));
    Rcpp::Rcout << "Contributors:" << std::endl;
    Rcpp::print(Rcpp::wrap(tmp[locus]));
    Rcpp::Rcout << std::endl;
  }
  Rcpp::Rcout << "################################################" << std::endl;      
  */

  counts_contributor_sets.push_back(counts);                            
}

// https://stackoverflow.com/a/19406536
void nested_loop_operation(std::vector<int> counters, 
                           std::vector<int>& length, 
                           int level,
                           const int num_contributors, 
                           const std::vector< std::vector<int> >& haps_in_mixture,
                           const std::vector<int>& haps_in_mixture_counts, 
                           const std::vector<int>& poi_profile,
                           const std::vector< std::vector<int> >& known_contributors,
                           const std::vector< std::unordered_set<int> >& mixture,
                           std::vector< std::vector<int> >& counts_contributor_sets);

void nested_loop_operation(std::vector<int> counters, 
                           std::vector<int>& length, 
                           int level,
                           const int num_contributors, 
                           const std::vector< std::vector<int> >& haps_in_mixture,
                           const std::vector<int>& haps_in_mixture_counts, 
                           const std::vector<int>& poi_profile,
                           const std::vector< std::vector<int> >& known_contributors,
                           const std::vector< std::unordered_set<int> >& mixture,
                           std::vector< std::vector<int> >& counts_contributor_sets) {
  if ((counters.size() + known_contributors.size()) != num_contributors) {
    Rcpp::stop("# unknown + # known != K");
  }
                           
  if (level == counters.size()) {
    /*
    We want 
      haps_in_mixture_indices[0] < haps_in_mixture_indices[1] < ...
      i.e.
      quit if
        haps_in_mixture_indices[0] >= haps_in_mixture_indices[1]
        or
        haps_in_mixture_indices[1] >= haps_in_mixture_indices[2]
        or 
        ...
    Can maybe be done elsewhere, but for now this is done...
    */    
    for (int set_i = 1; set_i < counters.size(); set_i++) {
      if (counters[set_i - 1] >= counters[set_i]) {
        return; // exit function, don't call
      }
    }
    
    analyse_set(num_contributors, 
                counters, haps_in_mixture, haps_in_mixture_counts, 
                poi_profile, 
                known_contributors, mixture, counts_contributor_sets);                              
  } else {
    for (counters[level] = 0; counters[level] < length[level]; counters[level]++) {
       nested_loop_operation(counters, length, level + 1,
                             num_contributors, 
                             haps_in_mixture, haps_in_mixture_counts,
                             poi_profile,
                             known_contributors, mixture, 
                             counts_contributor_sets);
    }
  }
}



//' Analyse mixture results
//' 
//' Calculate LR-like quantities by haplotype counts.
//' 
//' NOTE: Only takes up to 9 contributors!
//' 
//' @param mix_res Mixture result from [mixture_info_by_individuals_2pers()], 
//' [mixture_info_by_individuals_3pers()], [mixture_info_by_individuals_4pers()], 
//' [mixture_info_by_individuals_5pers()]
//' @param unique_haps_in_mixture Included unique haplotypes to use as elements in contributor sets.
//' @param unique_haps_in_mixture_counts Population counts of the included haplotypes
//' 
//' @return A list with numeric quantities
// [[Rcpp::export]]
Rcpp::List analyse_mixture_result(Rcpp::List& mix_res, 
                                  const Rcpp::List& unique_haps_in_mixture, 
                                  const Rcpp::List& unique_haps_in_mixture_counts) { 
  
  // donor 1: suspect, donor 2+3+...+K: unknown
  
  // No. contributors:
  std::vector< std::vector<int> > true_donors;

  if (!mix_res.containsElementNamed("donor1_profile")) {
    Rcpp::stop("No donor1 as expected");
  }
  
  std::vector<int> poi_profile = mix_res["donor1_profile"];
  
  // FIXME: Stupid, but fast
  if (mix_res.containsElementNamed("donor1_profile")) true_donors.push_back(mix_res["donor1_profile"]);
  if (mix_res.containsElementNamed("donor2_profile")) true_donors.push_back(mix_res["donor2_profile"]);
  if (mix_res.containsElementNamed("donor3_profile")) true_donors.push_back(mix_res["donor3_profile"]);
  if (mix_res.containsElementNamed("donor4_profile")) true_donors.push_back(mix_res["donor4_profile"]);
  if (mix_res.containsElementNamed("donor5_profile")) true_donors.push_back(mix_res["donor5_profile"]);
  if (mix_res.containsElementNamed("donor6_profile")) true_donors.push_back(mix_res["donor6_profile"]);
  if (mix_res.containsElementNamed("donor7_profile")) true_donors.push_back(mix_res["donor7_profile"]);
  if (mix_res.containsElementNamed("donor8_profile")) true_donors.push_back(mix_res["donor8_profile"]);
  if (mix_res.containsElementNamed("donor9_profile")) true_donors.push_back(mix_res["donor9_profile"]);
  if (mix_res.containsElementNamed("donor10_profile")) {
    Rcpp::List terms_i = get_mixture_analyse_failure("Only up to 9 contributors are supported.");      
    return terms_i;
  }

  int K = true_donors.size();
  
  if (K <= 1) {
    Rcpp::List terms_i = get_mixture_analyse_failure("Need at least 2 donors.");
    return terms_i;
  }
  
  int loci = true_donors[0].size();  
  for (int donor_i = 1; donor_i < K; ++donor_i) {
    if (loci != true_donors[donor_i].size()) {
      Rcpp::List terms_i = get_mixture_analyse_failure("Not all donors have same number of loci.");
      return terms_i;
    }
  }
  
  int LR_contrib_Hp = 0;
  int LR_contrib_Hd = 0;
  
  ////////////////////////////////////////////////////////////////////////////
  // Check if this is detectable as a K persons' mixture:
  // at least one locus must have K distinct alleles (so that all contributors are different)
  bool mixture_detectable = false;
  std::vector< std::unordered_set<int> > mixture(loci);
  
  for (int locus = 0; locus < loci; ++locus) {
    std::unordered_set<int> alleles; // counts number of unique alleles / allele masking
    
    for (int donor_i = 0; donor_i < K; ++donor_i) {
      alleles.insert(true_donors[donor_i][locus]);
    }

    // K alleles detected
    if (alleles.size() == K) {
      mixture_detectable = true;
      // must not break; mixture needs to be filled in case all is okay
    }
    
    mixture[locus] = alleles;
  }    
  
  if (!mixture_detectable) {
    Rcpp::List terms_i = get_mixture_analyse_failure("Allele masking makes the true number of donors undetectable.");      
    return terms_i;
  }
  ////////////////////////////////////////////////////////////////////////////
  
  // Now, all hs are different because one locus has at least K different alleles
  
  int m = unique_haps_in_mixture.size();

  if (unique_haps_in_mixture_counts.size() != m) {
    Rcpp::List terms_i = get_mixture_analyse_failure("unique_haps_in_mixture_counts.size() != m.");
    return terms_i;
  }

  if (m < K) {
    Rcpp::List terms_i = get_mixture_analyse_failure("m (number of unique included haplotypes) < K (number of contributors).");
    return terms_i;
  }
  
  std::vector< std::vector<int> > haps_in_mixture(m);
  std::vector<int> haps_in_mixture_counts(m);
  
  for (int hap_i = 0; hap_i < m; ++hap_i) {
    haps_in_mixture[hap_i] = Rcpp::as<std::vector<int>>(unique_haps_in_mixture[hap_i]);
    haps_in_mixture_counts[hap_i] = unique_haps_in_mixture_counts[hap_i];
  }

  /*
  ****************************************************************************
  * DENOMINATOR: Defence's hypothesis
  * CHECKING IF R+S+T = CSP.
  ****************************************************************************
  */
  /*
   Originally, an explicit loop nesting could be 
   used when number of contributors (K) are assumed known and fixed:
   
   for (int k1 = 0; k1 < (m-2); ++k1) {
     Rcpp::IntegerVector R1 = haps[k1];
     
     for (int k2 = (k1+1); k2 < (m-1); ++k2) {
       Rcpp::IntegerVector R2 = haps[k2];
        
       for (int k3 = (k2+1); k3 < m; ++k3) {
         Rcpp::IntegerVector R3 = haps[k3];
         ...
       }
     }
   }
  */
  std::vector< std::vector<int> > terms_Hd;
  std::vector<int> counters_Hd(K); // init 0
  std::vector<int> length_Hd(K);

  for (int donor_i = 0; donor_i < K; ++donor_i) {
    length_Hd[donor_i] = m - (K - (donor_i + 1)); 
    /*
    
    m: haps.size(); 
    K: number of contributors

    K = 2: (loop to m-1 and inner to m)
      donor_i = 0: m - (2 - (0 + 1)) = m - (2 - 1) = m - 1
      donor_i = 1: m - (2 - (1 + 1)) = m - (2 - 2) = m
    
    K = 3: (loop to m-2, m-1 and inner to m)
      donor_i = 0: m - (3 - (0 + 1)) = m - (3 - 1) = m - 2
      donor_i = 1: m - (3 - (1 + 1)) = m - (3 - 2) = m - 1
      donor_i = 2: m - (3 - (2 + 1)) = m - (3 - 3) = m
    */
  }       

  nested_loop_operation(counters_Hd, length_Hd, 
                        0, // level
                        K, // num_contributors
                        haps_in_mixture, haps_in_mixture_counts,
                        poi_profile, 
                        std::vector< std::vector<int> >(), // known_contributors for denom is none
                        mixture, 
                        terms_Hd); // counts_contributor_sets
  
  for (int term_i = 0; term_i < terms_Hd.size(); ++term_i) {
    int term = std::accumulate(terms_Hd[term_i].begin(), terms_Hd[term_i].end(), 1, std::multiplies<int>());
    LR_contrib_Hd += term;
  }
  
  /*
  ****************************************************************************
  * NUMERATOR: Prosecutor's hypothesis
  * CHECKING IF 
  *    donor1+U+V = CSP.
  ****************************************************************************
  */    
  // As above, but with Km1 = K - 1 instead of K
  int Km1 = K - 1;
  std::vector< std::vector<int> > terms_Hp;
  std::vector<int> counters_Hp(Km1); // init 0
  std::vector<int> length_Hp(Km1);
  for (int donor_i = 0; donor_i < Km1; ++donor_i) {
    length_Hp[donor_i] = m - (Km1 - (donor_i + 1)); 
  }
  
  const std::vector< std::vector<int> > known_contribs = { true_donors[0] };
  nested_loop_operation(counters_Hp, length_Hp, 
                        0, // level
                        K, // num_contributors
                        haps_in_mixture, haps_in_mixture_counts,
                        std::vector<int>(), // poi_profile := empty, no need to possibly decrease by one
                        known_contribs, // known_contributors for num simply assumed the first donor
                        mixture, 
                        terms_Hp); // counts_contributor_sets
  
  for (int term_i = 0; term_i < terms_Hp.size(); ++term_i) {
    int term = std::accumulate(terms_Hp[term_i].begin(), terms_Hp[term_i].end(), 1, std::multiplies<int>());
    LR_contrib_Hp += term;
  }
  
  /*
  ****************************************************************************
  * BUILDING RETURN VALUE
  ****************************************************************************
  */          
  Rcpp::List terms_i;
  terms_i["fail"] = false;
  terms_i["fail_reason"] = Rcpp::CharacterVector::create(NA_STRING); // vs. NA_STRING?
  terms_i["terms_Hp"] = terms_Hp;
  terms_i["terms_Hp_count"] = terms_Hp.size();
  terms_i["LR_contrib_Hp"] = LR_contrib_Hp;
  terms_i["terms_Hd"] = terms_Hd;
  terms_i["terms_Hd_count"] = terms_Hd.size();
  terms_i["LR_contrib_Hd"] = LR_contrib_Hd;
  terms_i["number_contributors"] = K;
  
  return terms_i;
}

//' Analyse mixture results in a vectorised fashion
//' 
//' Refer to [analyse_mixture_result()] for details. 
//' Essentially, [analyse_mixture_result()] is run on each element of `mixture_results`.
//' 
//' NOTE: Only takes up to 9 contributors!
//' 
//' @param mixture_results List of `n` mixture results from [mixture_info_by_individuals_2pers()], 
//' [mixture_info_by_individuals_3pers()], [mixture_info_by_individuals_4pers()], 
//' [mixture_info_by_individuals_5pers()]
//' @param unique_haps_in_mixture_list List of `n` included unique haplotypes, one for each element in `mix_res`
//' @param unique_haps_in_mixture_counts_list List of `n` population counts of the included unique haplotypes
//' 
//' @return A list with lists of numeric quantities
// [[Rcpp::export]]
Rcpp::List analyse_mixture_results(Rcpp::List& mixture_results, 
                                   const Rcpp::List& unique_haps_in_mixture_list, 
                                   const Rcpp::List& unique_haps_in_mixture_counts_list) { 
                                        
  int n = mixture_results.size();

  Rcpp::List CSP_terms(n);
  
  /*
  Contributors must also be included in the mixture
  */

  if (unique_haps_in_mixture_list.size() != n) {
    Rcpp::stop("unique_haps_in_mixture_list.size() != n");
  }

  if (unique_haps_in_mixture_counts_list.size() != n) {
    Rcpp::stop("unique_haps_in_mixture_counts_list.size() != n");
  }
    
  
  for (int i = 0; i < n; ++i) {
    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    Rcpp::List mix_res = mixture_results[i];
    Rcpp::List unique_haps_in_mixture = unique_haps_in_mixture_list[i];
    Rcpp::List unique_haps_in_mixture_counts = unique_haps_in_mixture_counts_list[i];
    Rcpp::List res = analyse_mixture_result(mix_res, unique_haps_in_mixture, unique_haps_in_mixture_counts);
    CSP_terms[i] = res;
  }

  //Rcpp::print(CSP_terms);
  
  return CSP_terms;
}

