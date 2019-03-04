/**
 class_Population.h
 Purpose: Header for C++ class Population.
 Details: C++ header.
  
 @author Mikkel Meyer Andersen
 */

#include "malan_types.h"

#include <vector>
#include <unordered_map>

class Population {
private:
  std::unordered_map<int, Individual*>* m_population = NULL;

public:
  Population(std::unordered_map<int, Individual*>* population);
  ~Population();
  std::unordered_map<int, Individual*>* get_population() const;
  int get_population_size() const;
  Individual* get_individual(int pid) const;
};


