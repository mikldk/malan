/**
 malan_types.h
 Purpose: Header for C++ classes.
 Details: C++ header.
  
 @author Mikkel Meyer Andersen
 */
 
#ifndef MALAN_TYPES_H
#define MALAN_TYPES_H

#define CHECK_ABORT_EVERY 10000


/*
Rcpp::Xptr<p: pointer, set_delete_finalizer: bool>

Rcpp:
set_delete_finalizer: if set to true, a finalizer will be registered for the external pointer. 
The finalizer is called when the Xptr is garbage collected. 
The finalizer is merely a call to the delete operator or the pointer so you need to make sure the pointer can be "delete"'d this way (has to be a C++ object)

Here, strategy:
All except Population in simulate* functions have RCPP_XPTR_2ND_ARG (false), 
but Population has true and cleans up.
*/

#define RCPP_XPTR_2ND_ARG false // do not call finalisers, others are responsible
#define RCPP_XPTR_2ND_ARG_CLEANER true // finaliser called, responsible for cleaning up

class Individual;
class Pedigree;
class Population;

//class SimulateChooseFather;
//class WFRandomFather;
//class GammaVarianceRandomFather;

#include "helper_Individual.h"

#include "class_Individual.h"
#include "class_Pedigree.h"
#include "class_Population.h"
#include "class_SimulateChooseFather.h"




namespace{
// a little helper that should IMHO be standardized
template<typename T>
std::size_t make_hash(const T& v){
  return std::hash<T>()(v);
}

// adapted from boost::hash_combine
void hash_combine(std::size_t& h, const std::size_t& v){
  h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
}

// hash any container
template<typename T>
struct hash_container{
  size_t operator()(const T& v) const{
    size_t h = 0;
    for( const auto& e : v ) {
      hash_combine(h, make_hash(e));
    }
    return h;
  }
};
}

namespace std{
// support for vector<T> if T is hashable
// (the T... is a required trick if the vector has a non-standard allocator)
template<>
struct hash< std::vector<int> > : hash_container< std::vector<int> > {};

struct equal_to_intvec : binary_function< std::vector<int> , std::vector<int> , bool> {
  bool operator() (const std::vector<int>& x, const std::vector<int>& y) const{
    if (x.size() != y.size()){
      return false;
    }
    
    int n = x.size();
    
    for (int i = 0; i < n; ++i) {
      if (x[i] != y[i]) {
        return false;
      }
    }
    
    return true;
  }
};

// the same for map<T,U> if T and U are hashable
template<typename... T>
struct hash<map<T...>> : hash_container<map<T...>> {};

// simply add more containers as needed
}

#endif
