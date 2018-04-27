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

set_delete_finalizer: if set to true, a finalizer will be registered for the external pointer. 
The finalizer is called when the Xptr is garbage collected. 
The finalizer is merely a call to the delete operator or the pointer so you need to make sure the pointer can be "delete" d this way (has to be a C++ object)
*/
//#define RCPP_XPTR_2ND_ARG true // ensures that finaliser is called
#define RCPP_XPTR_2ND_ARG false // do not call finalisers, I guess we live with some memory leaks for now...!!! FIXME

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

#endif
