/**
 class_SimulateChooseFather.h
 Purpose: Header for C++ class SimulateChooseFather.
 Details: C++ header.
  
 @author Mikkel Meyer Andersen
 */

#include <RcppArmadillo.h> // FIXME: Avoid Rcpp here? Only in api_* files?
//#include <Rcpp.h>


class SimulateChooseFather {
  public:
    virtual void update_state_new_generation() = 0;
    virtual int get_father_i() = 0;
};

class WFRandomFather: public SimulateChooseFather {
  private:
    double m_population_size;
    
  public:
    WFRandomFather(size_t population_size);
    void update_state_new_generation();
    int get_father_i();
};


class GammaVarianceRandomFather: public SimulateChooseFather {
  private:
    // fixed for entire simulate
    size_t m_population_size;
    double m_gamma_parameter_shape;
    double m_gamma_parameter_scale;
    
    // new for each generation
    arma::vec m_fathers_prob_cum;
    arma::uvec m_fathers_prob_perm;
    
  public:
    GammaVarianceRandomFather(size_t population_size, double gamma_parameter_shape, double gamma_parameter_scale);
    void update_state_new_generation();
    int get_father_i();
 };
