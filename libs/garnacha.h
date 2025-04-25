#ifndef GARNACHA_H
#define GARNACHA_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_statistics.h>
#include <set>
#include <list>
#include "alea.h"
#include "basics.h"

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


struct losfen{
  int nufen;
  string *estrus;
  double *eners;
  int *demeama;
};

struct parligR { //en orden de menor a mayor en términos de Boltzmann
  //generales
  double ka; //0.001
  double Rtot;
  double Ltot;
  double kda;
  double kdb;
  //específicos
  int numec;
  double *B;
  double *dis;
  double *kd;

};

class GaRNAcha
{
  
  friend void cue_subopt(const char  *structure, float energy, void *data);
  friend void lle_subopt(const char  *structure, float energy, void *data);
  friend int lrmod_se(const gsl_vector *x, void *params, gsl_vector *f); //static equilibrium
  friend int lrmod_se_df(const gsl_vector *x, void *params, gsl_matrix *J);
  friend int lrmod_se_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);
  friend int lrmod(const gsl_vector *x, void *params, gsl_vector *f); //static equilibrium
  friend int lrmod_df(const gsl_vector *x, void *params, gsl_matrix *J);
  friend int lrmod_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);
  
public:
  GaRNAcha();
  GaRNAcha(int tam);
  void start_rng(Alea& jacta);
  GaRNAcha(Alea& jacta);
  GaRNAcha(Alea& jacta, int tam);
  void set_alphabet();
  void set_default_parameters();
  void copy_seq(GaRNAcha &daughter); 
  void set_size(int tam);
  void rand_seq();
  void rand_seq(int *ntpercent);
  void rand_seq_exactcomp(int *ntcomp);
  void get_seq_comp(int *ntcomp);
  void clear_seq();
  void clear_size();
  void clear_mfes();
  void clear_plas_phen();
  void clear_allmfes();
  void clear_rng();
  void clear_mrep();
  void clear_mef();
  void clear_plr_sp();
  void clear_complr();
  void clear_cong();
  void clear(); //all but rng
  bool already_plastic_phen();
  void print_seq(ostream& fs);
  void set_seq(string unasec);
  void set_seq(const char *unasec);
  void set_seq_allow_lc(string unasec);
  void set_seq_allow_lc(const char *unasec);
  void seq_to_upper();
  void force_inverse_fold(string target);
  double inv_fold(string target);
  double inv_fold(string target, string inseq);
  int bp_distance(string st1, string st2);
  int bp_distance(string st1);
  int bp_distance(const char *st1, const char *st2);
  int bp_distance(const char *st1);
  void get_mfe_str();
  int seq_size();
  double get_temperature();
  void set_temperature(double t);
  int hamming_d_seq(const char *s1, const char *s2);
  int hamming_d_seq(string s1, string s2);
  int hamming_d_seq(const char *s1);
  int hamming_d_seq(string s1);
  double tree_d(string stru1, string stru2);
  double tree_d(const char *stru1, const char *stru2);
  double tree_d(string stru1);
  double tree_d(const char *stru1);
  void print_mfestr(ostream& fs);
  double eval_structure(string stru);
  double eval_structure(const char *stru);
  bool is_this_mfestr(string lastr);
  void mutate_seq(int pos, char nuc);
  void mutate_seq_fd(int pos); //force diff
  void rand_mut();
  void rand_mut_fd();
  void set_e_band(int eb); //eb como en Vienna 300,200
  int e_band();
  void plastic_phenotypes();
  void plastic_phenotypes(int eb);
  int count_compatible_structures();
  string one_compatible_structure(int cual);
  void print_phenotypes(ostream& sal);
  void get_all_mfe_strs();
  int number_of_mfestrs();
  string an_mfes(int cual);
  void print_all_mfe_strs(ostream& sal);
  string nt_sequence();
  string mfe_structure();
  double mfenergy();
  double m_neighbors_with_mfe_str(string target); //mrobustness
  double m_neighbors_with_mfe_str(const char *targ);
  double mutational_access(string whstr); //fracdemutsquevanafen*3size
  void m_repertoire(); //incluye estructura wt
  bool already_m_rep();
  void copy_m_rep(GaRNAcha &froms);
  void print_mrep(ostream& sal);
  double count_in_mrep(int cual);
  string struct_in_mrep(int cual);
  int size_of_mrep();
  void partition_function();
  bool already_partition_function();
  double p_Boltzmann(int cual);
  double p_Boltzmann(string unphen);
  double plastic_robustness();
  void get_monomer_eq_freqs();
  double return_cff();
  void set_cff(double ncff);
  double w_no_plast(const char *optstr);
  double w_no_plast(string optstr);
  double w_plast_AF(const char *optstr);
  double w_plast_AF(string optstr);
  void next_sampled_seq_w(int steps, const char *estru);
  void next_sampled_seq_w(int steps, string estru);
  void sample_w(ostream& fs, int cuantas, int multip, string estru);
  void sample_w(ostream& fs, int cuantas, int multip, const char *estru);
  string phen_rand_mut();
  string phen_rand_mut_with_seq(string &laseq);
  int maccesible_strucs_from_struc(int vueltas, string esdesde, string *accest, int *contando, int numaxest);
  int maccesible_strucs_from_struc(int vueltas, const char *esdesde, string *accest, int *contando, int numaxest);
  int seqs_neighboring_struct(ostream& fs, string oldst, string newstr, int tries);

  
//  double w_of_struct(const char *estru, const char *optstr);
  //func para tamanio
  void set_def_plr(double rt, double lt, double kac, double kd_a, double kd_b); //puede cambiarlos
 // void set_plr(double rt, double lt, double kac, double kd_a, double kd_b, string opt);
  void set_spec_plr(string opt);
  bool already_spec_plr();
  double get_eqcomp(); //no incluye set_spec_plr
  double get_eqcomp_num(); //no incluye set_spec_plr
  double get_eqcomp_se(); //no incluye set_spec_plr
  bool already_get_eqcomp();
  double get_compLRT();
  double get_compLR(int i);
//  double amount_of_structure(int i); //hecha pa version anterior, creo que no es necesaria. Si sí, actualizar
  double dist_to_opt(int i);
  double get_kd(int i);
  
  void congruence_analysis();
  bool already_cong();
  void copy_cong_an(GaRNAcha &froms);
  void cong_in_one_phen(int cual, double &demu, double &depl);
  double jaccard();
  double overlap();
  double pearson_c();
  int size_phen_taccess();

private:
  char alphabet[4]; // = {'A', 'C', 'G', 'U'};
  char *seq;
  Alea est;
  Basics basic;
  bool yarng;
  bool yaseq;
  bool yasize;
  bool yamfe;
  bool yaband;
  bool yaplas;
  bool yaallmfe; 
  bool yamrep;
  bool yapafu;
  losfen allmfes;
  losfen otipos;
  double mfe;
  char *mfestr;
  int band; //kcal * 100
  int size;
  int mrep_size;
  double *mrep_counts;
  string *mrep;
  double pafu;
  
  double kb;
  double e;
  double *meqfreqs;
  bool yamef;
  
  bool yaplr;
  bool yawlr;
  bool yaplrd;
  parligR paramLR;
  
  double *compLR; // de menor a mayor en p de Boltzmann
  bool yacompLR;
  double compLRT;
  
  double cff; //constant for fitness. Fitness component: 1/(cff+d) con d, distancia a target
  double jacc,olap,pear;
  int sunion;
  double *conmu, *conpl;
  bool yacong;
  
};

#endif

//g++ -Wall -o ejem ejem.cc $MYLIBS/alea.o $MYLIBS/basics.o garnacha.o -lgsl -lgslcblas -lm -lRNA
//g++ -Wall -c garnacha.cc




