/**
 helder_hash.h
 Purpose: C++ helper functions for hash custom types, here vector<int>
 Details: C++ header.
  
 @author Mikkel Meyer Andersen
 */

#ifndef HELPER_HASH_H
#define HELPER_HASH_H

#include <functional>
#include <vector>

namespace std {
  // a little helper
  template<typename T>
  std::size_t make_hash(const T& v) {
    return std::hash<T>()(v);
  }

  // hash any container
  template<typename T>
  struct hash_container {
    size_t operator()(const T& v) const {
      size_t h = 0;
      
      for (const auto& e: v) {
        size_t e_hash = make_hash(e);
        // adapted from boost::hash_combine
        h ^= e_hash + 0x9e3779b9 + (h << 6) + (h >> 2);        
      }
      
      return h;
    }
  };
  
  // hash any pair
  template<typename T>
  struct hash_pair {
    size_t operator()(const T& v) const {
      size_t lhs = v.first;
      size_t rhs = v.second;
      
      // adapted from boost::hash_combine
      lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
      return lhs;
    }
  };
  
  template<>
  struct hash< std::vector<int> > : hash_container< std::vector<int> > {};
  
  template<>
  struct hash< std::pair<int, int> > : hash_pair< std::pair<int, int> > {};
  
}

#endif

