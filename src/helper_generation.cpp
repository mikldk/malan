#include "malan_types.h"

// end_generation_number pointer easier than traversing afterwards 
// to find biggest number
void update_generation(Individual* indv, 
                       int generation_number, 
                       int* end_generation_number, 
                       const int modifier) {
  
  indv->set_generation(generation_number);
  
  std::vector<Individual*> children = *(indv->get_children());
  
  int new_generation_number = modifier + generation_number;
  
  if (new_generation_number > *end_generation_number) {
    *end_generation_number = new_generation_number;
  }
  
  for (auto child : children) {
    update_generation(child, new_generation_number, end_generation_number, modifier);
  }
}

