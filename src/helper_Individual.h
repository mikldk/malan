/**
 helder_Individual.h
 Purpose: Header for C++ helper functions for Individual class.
 Details: C++ header.
  
 @author Mikkel Meyer Andersen
 */
 
#ifndef HELPER_INDIVIDUALS_H
#define HELPER_INDIVIDUALS_H

#include "malan_types.h"

#include <vector>

bool find_path_from_root_to_dest(Individual* root, std::vector<Individual*>& path, const Individual* dest);

std::vector<int> draw_autosomal_genotype(
    const std::vector<double>& allele_cumdist_theta,
    const int alleles_count);
  
#endif
