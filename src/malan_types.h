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

#endif
