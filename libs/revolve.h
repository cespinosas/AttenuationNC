#ifndef REVOLVE_H
#define REVOLVE_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <set>
#include <list>
#include "alea.h"
#include "basics.h"
#include "garnacha.h"

extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/dist_vars.h>
#include <ViennaRNA/treedist.h>
#include <ViennaRNA/RNAstruct.h>
#include <ViennaRNA/stringdist.h>
#include <ViennaRNA/subopt.h>
#include <ViennaRNA/inverse.h>
#include <ViennaRNA/part_func.h>
}

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;
using std::ifstream;
using std::istream;
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;



class REvolve{
public:
  REvolve();
  REvolve(Alea& jacta);
  void set_default_parameters();
  void start_rng(Alea& jacta);
  void set_paramsG(double muratep, int popsizep, string estru, int eb);
  void set_optimum(string estru);
  void set_paramsAF(double cefefe); //
  void set_paramsLR(double rt, double lt, double kac, double kd_a, double kd_b); //
  void start_pop(GaRNAcha &founder);
  void clear_pop();
  void mutation_to_pop();
  void get_pop_wwop(); //llamar estas antes de empezar la evolución
  double stddev_wwop();
  void get_pop_waf();
  double stddev_waf();
  void get_pop_wlr();
  void get_pop_wlr_se();
  double stddev_wlr();
  void one_gen_wop();
  void one_gen_waf();
  void one_gen_wlr();

  double meanw_wop();
  double meanw_af();
  double meanw_lr();

  double maxw_wop();
  double maxw_af();
  double maxw_lr();

  void Rmut_in_pop(double &mean, double &max, double &stddev);
  void Boltzmann_in_pop(double &mean, double &max, double &stddev);
  void Jaccard_in_pop(double &mean, double &max, double &stddev); //sacar cong
  void overlap_in_pop(double &mean, double &max, double &stddev);
  void Pearson_in_pop(double &mean, double &max, double &stddev);
  void mutational_access_in_pop(double &mean, double &max, double &stddev); //sacar m_rep antes
  void plasticity_access_in_pop(double &mean, double &max, double &stddev);

  void nucleotide_div(double &mean, double &max, double &stddev);
  double heterozygosity();
  int inds_with_opt_mfe();
  string dominant_seq(int &hmany);
  void freeze_population_in_file(string nombre);
  void freeze_population_in_file();
  void congruence_in_pop();
  void m_rep_in_pop();
  
  
private:
  Alea est;
  Basics basic;
  int popsize;
  int seqsize;
  GaRNAcha *pop;
  GaRNAcha *nepop;
  double *wwop, *waf, *wlr;
  double meanwwop, meanwaf, meanwlr;
  double maxwwop, maxwaf, maxwlr;
  double minwwop, minwaf, minwlr;
  double sumwwop, sumwaf, sumwlr;
  int eband;
  string optist;
  double mu;
  bool yapG, yapAF, yapLR;
  double cff;
  double rtot, ltot, ka, kda, kdb;
  bool yapop;

};

#endif

//g++ -Wall -o ejem ejem.cc $MYLIBS/alea.o $MYLIBS/basics.o garnacha.o -lgsl -lgslcblas -lm -lRNA
//g++ -Wall -c garnacha.cc
