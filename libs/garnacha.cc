/*#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <set>
#include <list>
#include <string>
#include "alea.h"
#include "basics.h"
 */
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
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;


//Friends:

void cue_subopt(const char  *structure, float energy, void *data) {
  if (structure)
    (*((int *)data))++;
  return;
}

void lle_subopt(const char  *structure, float energy, void *data) {
  if (structure) {
    (*((losfen *)data)).estrus[(*((losfen *)data)).nufen] = structure;
    (*((losfen *)data)).eners[(*((losfen *)data)).nufen] = energy;
    (*((losfen *)data)).nufen++;
  }
  return;
}

int lrmod_se(const gsl_vector *x, void *params, gsl_vector *f) {
  double ka, *kd, *B, *dis, Rtot, Ltot;
  Rtot = ((struct parligR *) params)->Rtot;
  Ltot = ((struct parligR *) params)->Ltot;
  ka = ((struct parligR *) params)->ka;
  int numec = ((struct parligR *) params)->numec;
  kd = new double[numec];
  B = new double[numec];
  dis = new double[numec];
  int i;
  for (i=0; i < numec; i++) {
    kd[i]  = ((struct parligR *) params)->kd[i];
    B[i]  = ((struct parligR *) params)->B[i];
    B[i] *= Rtot;
    dis[i]  = ((struct parligR *) params)->dis[i];
  }
  
  double *comp;
  comp = new double[numec];
  for (i=0; i< numec; i++)
    comp[i] = gsl_vector_get(x, i);
  
  double sutoc = 0; //suma de todos los complejos
  for (i=0; i< numec; i++)
    sutoc += comp[i];
  
  double *y;
  y = new double[numec];
  for (i=0; i< numec; i++)
    y[i] = (ka*(B[i] - comp[i])*(Ltot - sutoc)) - (kd[i]*comp[i]);
    
  for (i=0; i< numec; i++)
    gsl_vector_set(f, i, y[i]);
  
  delete [] y;
  delete [] comp;
  delete [] kd;
  delete [] B;
  delete [] dis;
  return GSL_SUCCESS;
}

int lrmod_se_df(const gsl_vector *x, void *params, gsl_matrix *J) {
  double ka, *kd, *B, *dis, Rtot, Ltot;
  Rtot = ((struct parligR *) params)->Rtot;
  Ltot = ((struct parligR *) params)->Ltot;
  ka = ((struct parligR *) params)->ka;
  int numec = ((struct parligR *) params)->numec;
  kd = new double[numec];
  B = new double[numec];
  dis = new double[numec];
  int i,j;
  for (i=0; i < numec; i++) {
    kd[i]  = ((struct parligR *) params)->kd[i];
    B[i]  = ((struct parligR *) params)->B[i];
    B[i] *= Rtot;
    dis[i]  = ((struct parligR *) params)->dis[i];
  }
  
  double *comp;
  comp = new double[numec];
  for (i=0; i< numec; i++)
    comp[i] = gsl_vector_get(x, i);
  
  double sutoc = 0; //suma de todos los complejos
  for (i=0; i< numec; i++)
    sutoc += comp[i];
  
  double **df;
  df = new double*[numec];
  for (i=0; i < numec; i++)
    df[i] = new double[numec];
  
  for (i=0; i < numec; i++) {
    for (j=0; j < numec; j++) {
      if (i==j) {
        df[i][j] = (ka*(comp[i] + sutoc - B[i] - Ltot)) -kd[i];
      } else {
        df[i][j] = ka*(comp[i] - B[i]);
      }
    }
  }
  
  for (i=0; i < numec; i++)
    for (j=0; j < numec; j++)
      gsl_matrix_set(J, i, j, df[i][j]);
  
  for (i=0; i < numec; i++)
    delete [] df[i];
  delete [] df;
  delete [] comp;
  delete [] kd;
  delete [] B;
  delete [] dis;
  
  return GSL_SUCCESS;
}

int lrmod_se_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
  lrmod_se(x, params, f);
  lrmod_se_df(x, params, J);
  return GSL_SUCCESS;
}

//de

int lrmod(const gsl_vector *x, void *params, gsl_vector *f) {
  double ka, *kd, *B, *dis, Rtot, Ltot;
  Rtot = ((struct parligR *) params)->Rtot;
  Ltot = ((struct parligR *) params)->Ltot;
  ka = ((struct parligR *) params)->ka;
  int numec = ((struct parligR *) params)->numec;
  kd = new double[numec];
  B = new double[numec];
  dis = new double[numec];
  int i;
  for (i=0; i < numec; i++) {
    kd[i]  = ((struct parligR *) params)->kd[i];
    B[i]  = ((struct parligR *) params)->B[i];
    dis[i]  = ((struct parligR *) params)->dis[i];
  }
  
  double *comp;
  comp = new double[numec];
  for (i=0; i< numec; i++)
    comp[i] = gsl_vector_get(x, i);
  
  double sutoc = 0; //suma de todos los complejos
  for (i=0; i< numec; i++)
    sutoc += comp[i];
  
  double *y;
  y = new double[numec];
  for (i=0; i< numec; i++)
    y[i] = (ka*B[i]*(Rtot - sutoc)*(Ltot - sutoc)) - (kd[i]*comp[i]);
  
  for (i=0; i< numec; i++)
    gsl_vector_set(f, i, y[i]);
  
  delete [] y;
  delete [] comp;
  delete [] kd;
  delete [] B;
  delete [] dis;
  return GSL_SUCCESS;
}

int lrmod_df(const gsl_vector *x, void *params, gsl_matrix *J) {
  double ka, *kd, *B, *dis, Rtot, Ltot;
  Rtot = ((struct parligR *) params)->Rtot;
  Ltot = ((struct parligR *) params)->Ltot;
  ka = ((struct parligR *) params)->ka;
  int numec = ((struct parligR *) params)->numec;
  kd = new double[numec];
  B = new double[numec];
  dis = new double[numec];
  int i,j;
  for (i=0; i < numec; i++) {
    kd[i]  = ((struct parligR *) params)->kd[i];
    B[i]  = ((struct parligR *) params)->B[i];
    dis[i]  = ((struct parligR *) params)->dis[i];
  }
  
  double *comp;
  comp = new double[numec];
  for (i=0; i< numec; i++)
    comp[i] = gsl_vector_get(x, i);
  
  double sutoc = 0; //suma de todos los complejos
  for (i=0; i< numec; i++)
    sutoc += comp[i];
  
  double **df;
  df = new double*[numec];
  for (i=0; i < numec; i++)
    df[i] = new double[numec];
  
  for (i=0; i < numec; i++) {
    for (j=0; j < numec; j++) {
      if (i==j) {
        df[i][j] = (ka*B[i]*((2*sutoc) - Rtot - Ltot)) -kd[i];
      } else {
        df[i][j] = ka*B[i]*((2*sutoc) - Rtot - Ltot);
      }
    }
  }
  
  for (i=0; i < numec; i++)
    for (j=0; j < numec; j++)
      gsl_matrix_set(J, i, j, df[i][j]);
  
  for (i=0; i < numec; i++)
    delete [] df[i];
  delete [] df;
  delete [] comp;
  delete [] kd;
  delete [] B;
  delete [] dis;
  
  return GSL_SUCCESS;
}

int lrmod_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
  lrmod(x, params, f);
  lrmod_df(x, params, J);
  return GSL_SUCCESS;
}


//public


GaRNAcha::GaRNAcha() {
  set_default_parameters();
}

GaRNAcha::GaRNAcha(int tam) {
  set_default_parameters();
  set_size(tam);
}

void GaRNAcha::start_rng(Alea& jacta) //assigns rng to jacta
{
  est = jacta;
  basic.start_rng(est);
  yarng = true;
}

GaRNAcha::GaRNAcha(Alea& jacta) {
  set_default_parameters();
  start_rng(jacta);
}

GaRNAcha::GaRNAcha(Alea& jacta, int tam) {
  set_default_parameters();
  start_rng(jacta);
  set_size(tam);
}

void GaRNAcha::set_alphabet() {
  alphabet[0] = 'A';
  alphabet[1] = 'C';
  alphabet[2] = 'G';
  alphabet[3] = 'U';
}

void GaRNAcha::set_default_parameters() {
  kb = 0.001987204259;
  e = 2.71828182845904523536;
  cff = 0.01;
  yaseq = false;
  size = -1;
  yasize = false;
  mfe = 99999;
  yamfe = false;
  yarng = false;
  band = -1;
  yaband = false;
  set_alphabet();
  yaplas = false;
  yaallmfe = false;
  yamrep = false;
  mrep_size = -1;
  yapafu= false;
  yamef = false;
  yaplr = false;
  yawlr = false;
  yaplrd = false;
  yacompLR = false;
  yacong = false;
}

void GaRNAcha::copy_seq(GaRNAcha &daughter) {
  daughter.start_rng(est);
  daughter.set_default_parameters();
  daughter.set_seq(seq);
}

void GaRNAcha::set_size(int tam) {
  if (yasize) {
    cout << "[Error]: Sequence length had been declared when you called GaRNAcha::set_size.\n";
    exit(1);
  }
  size = tam;
  yasize = true;
}

void GaRNAcha::rand_seq() {
  if (yaseq) {
    cout << "[Error]: Sequence already existed when you called GaRNAcha::rand_seq.\n";
    exit(1);
  }
  if (!yasize) {
    cout << "[Error]: Sequence length had not been declared when you called GaRNAcha::rand_seq.\n";
    exit(1);
  }
//  seq = vrna_random_string(size, "ACGU"); //si declaro tiene que ser char *seq
  int i;
  seq = new char[size+1];
  for (i=0; i < size; i++)
    seq[i] = alphabet[est.randint(0,4)];
  seq[size] = '\0';
  yaseq = true;
}

void GaRNAcha::rand_seq(int *ntpercent) {
  if (yaseq) {
    cout << "[Error]: Sequence already existed when you called GaRNAcha::rand_seq.\n";
    exit(1);
  }
  if (!yasize) {
    cout << "[Error]: Sequence length had not been declared when you called GaRNAcha::rand_seq.\n";
    exit(1);
  }
  if ((ntpercent[0]+ntpercent[1]+ntpercent[2]+ntpercent[3]) != 100) {
    cout << "[Error]: Nt percentages do not add up to 100 when you called GaRNAcha::rand_seq.\n";
    exit(1);
  }
  int i, cual, j, van;
  seq = new char[size+1];
  for (i=0; i < size; i++) {
    van = 0;
    cual = est.randint(0,100);
    for (j=0; j < 4; j++) {
      van += ntpercent[i];
      if (cual < van) {
        seq[i] = alphabet[j];
        break;
      }
    }
  }
  seq[size] = '\0';
  yaseq = true;
}

void GaRNAcha::rand_seq_exactcomp(int *ntcomp) {
  if (yaseq) {
    cout << "[Error]: Sequence already existed when you called GaRNAcha::rand_seq_exactcomp.\n";
    exit(1);
  }
  if (!yasize) {
    cout << "[Error]: Sequence length had not been declared when you called GaRNAcha::rand_seq_exactcomp.\n";
    exit(1);
  }
  if ((ntcomp[0]+ntcomp[1]+ntcomp[2]+ntcomp[3]) != size) {
    cout << "[Error]: Nt composition does not add up to sequence length when you called GaRNAcha::rand_seq_exactcomp.\n";
    exit(1);
  }
  int *wcomp;
  wcomp = new int[4];
  basic.copy_vector(ntcomp, wcomp, 4);
  int i, cual;
  seq = new char[size+1];
  for (i=0; i < size; i++) {
    do {
      cual = est.randint(0, 4);
    } while (wcomp[cual] <=0 );
    seq[i] = alphabet[cual];
    wcomp[cual]--;
  }
  seq[size] = '\0';
  delete [] wcomp;
  yaseq = true;
}

void GaRNAcha::get_seq_comp(int *ntcomp) {
  if (!yaseq) {
    cout << "[Error]: Sequence did not exist when you called GaRNAcha::get_seq_comp.\n";
    exit(1);
  }
  if (!yasize) {
    cout << "[Error]: Sequence length had not been declared when you called GaRNAcha::get_seq_comp.\n";
    exit(1);
  }
  basic.fillv0(ntcomp, 4);
  int i, j;
  for (i=0; i < 4; i++)
    for (j=0; j < size; j++)
      if (seq[j] == alphabet[i])
        ntcomp[i]++;
}

void GaRNAcha::clear_seq() {
  if (!yaseq) {
    cout << "[Error]: Sequence did not exist when you called GaRNAcha::clear_seq.\n";
    exit(1);
  }
  if (yamfe)
    clear_mfes();
  yaseq = false;
  delete [] seq;
  
  if (yamfe)
    clear_mfes();
  if (yaplas)
    clear_plas_phen();
  if (yaallmfe)
    clear_allmfes();
  if (yamrep)
    clear_mrep();
  if (yamef)
    clear_mef();
  if (yaplr)
    clear_plr_sp();
  if (yacompLR)
    clear_complr();
  if (yacong)
    clear_cong();
  yapafu = false;
}

void GaRNAcha::clear_size() {
  if (!yasize) {
    cout << "[Error]: Sequence length had not been declared when you called GaRNAcha::clear_size.\n";
    exit(1);
  }
  size = -1;
  yasize = false;
}

void GaRNAcha::clear_mfes() {
  if (!yamfe) {
    cout << "[Error]: MFE had not been obtained when you called GaRNAcha::clear_mfes.\n";
    exit(1);
  }
  mfe = 99999;
  delete [] mfestr;
  yamfe = false;
}

void GaRNAcha::clear_plas_phen() {
  otipos.nufen = 0;
  delete [] otipos.estrus;
  delete [] otipos.eners;
  delete [] otipos.demeama;
  yaplas = false;
}

void GaRNAcha::clear_allmfes() {
  allmfes.nufen = 0;
  delete [] allmfes.estrus;
  delete [] allmfes.eners;
  delete [] allmfes.demeama;
  yaallmfe = false;
}

void GaRNAcha::clear_rng() {
  est.close_rng();
}

void GaRNAcha::clear_mrep() {
  mrep_size = -1;
  delete [] mrep_counts;
  delete [] mrep;
  yamrep = false;
}

void GaRNAcha::clear_mef() {
  delete [] meqfreqs;
  yamef = false;
}

void GaRNAcha::clear_plr_sp() {
  delete [] paramLR.B;
  delete [] paramLR.dis;
  delete [] paramLR.kd;
  paramLR.numec = 0;
  yaplr = false;
}

void GaRNAcha::clear_complr() {
  delete [] compLR;
  compLRT = 0;
  yacompLR = false;
}

void GaRNAcha::clear_cong() {
  delete [] conmu;
  delete [] conpl;
  yacong = false;
}

void GaRNAcha::clear() {
  if (yaseq)
    clear_seq();
  if (yasize)
    clear_size();
  if (yamfe)
    clear_mfes();
  if (yaplas)
    clear_plas_phen();
  if (yaallmfe)
    clear_allmfes();
  if (yamrep)
    clear_mrep();
  yapafu = false;
  if (yamef)
    clear_mef();
  if (yaplr)
    clear_plr_sp();
  if (yacompLR)
    clear_complr();
  if (yacong)
    clear_cong();
}

bool GaRNAcha::already_plastic_phen() {
  return yaplas;
}

void GaRNAcha::print_seq(ostream& fs) {
  if (!yaseq) {
    cout << "[Error]: Sequence did not exist when you called GaRNAcha::print_seq.\n";
    exit(1);
  }
  fs << seq;
}

void GaRNAcha::set_seq(string unasec) {
  if (yasize && (size != unasec.length())) {
    cout << "[Error]: Sequence length of new sequence did not coincide with expected size when you called GaRNAcha::set_seq.\n";
    exit(1);
  }
  if (yaseq){
    cout << "[Error]: Sequence already existed when you called GaRNAcha::set_seq.\n";
    exit(1);
  }
  if (!yasize) {
    size = unasec.length();
    yasize = true;
  }
  seq = new char[size+1];
  strcpy(seq, unasec.c_str());
  yaseq= true;
  seq_to_upper();
}

void GaRNAcha::set_seq(const char *unasec) {
  if (yasize && (size != strlen(unasec))) {
    cout << "[Error]: Sequence length of new sequence did not coincide with expected size when you called GaRNAcha::set_seq.\n";
    exit(1);
  }
  if (yaseq){
    cout << "[Error]: Sequence already existed when you called GaRNAcha::set_seq.\n";
    exit(1);
  }
  if (!yasize) {
    size = strlen(unasec);
    yasize = true;
  }
  seq = new char[size+1];
  strcpy(seq, unasec);
  yaseq= true;
  seq_to_upper();
}

void GaRNAcha::set_seq_allow_lc(string unasec) {
  if (yasize && (size != unasec.length())) {
    cout << "[Error]: Sequence length of new sequence did not coincide with expected size when you called GaRNAcha::set_seq.\n";
    exit(1);
  }
  if (yaseq){
    cout << "[Error]: Sequence already existed when you called GaRNAcha::set_seq.\n";
    exit(1);
  }
  if (!yasize) {
    size = unasec.length();
    yasize = true;
  }
  seq = new char[size+1];
  yaseq= true;
  strcpy(seq, unasec.c_str());
}

void GaRNAcha::set_seq_allow_lc(const char *unasec) {
  if (yasize && (size != strlen(unasec))) {
    cout << "[Error]: Sequence length of new sequence did not coincide with expected size when you called GaRNAcha::set_seq.\n";
    exit(1);
  }
  if (yaseq){
    cout << "[Error]: Sequence already existed when you called GaRNAcha::set_seq.\n";
    exit(1);
  }
  if (!yasize) {
    size = strlen(unasec);
    yasize = true;
  }
  seq = new char[size+1];
  yaseq= true;
  strcpy(seq, unasec);
}

void GaRNAcha::seq_to_upper() {
  if (!yaseq) {
    cout << "[Error]: Sequence had not been declared when you called GaRNAcha::seq_to_upper.\n";
    exit(1);
  }
  vrna_seq_toupper(seq);
}

void GaRNAcha::force_inverse_fold(string target) {
  if (yaseq) {
    cout << "[Error]: Sequence already defined when you called GaRNAcha::inverse_fold.\n";
    exit(1);
  }
  if (!yasize) {
    size = target.length();
    yasize = true;
  }
  bool malo = true;
  int i, cont = 0;
  double d;
  do {
    rand_seq();
    d = inv_fold(target);
    if (d==0)
      malo = false;
    else {
      get_all_mfe_strs();
      for (i=0; i < allmfes.nufen; i++)
        if (allmfes.estrus[i] == target) {
          malo = false;
          break;
        }
      clear_allmfes();
    }
    if (malo)
      clear_seq();
    cont++;
  } while (malo && (cont < 3000));
  if (malo) {
    cout << "[Error]: GaRNAcha::force_inverse_fold was unable to find a sequence that produced the target structure.\n";
    exit(1);
  }
}

double GaRNAcha::inv_fold(string target) {
  //returns distance to target in terms of energy
  if (yasize && (size != target.length())) {
    cout << "[Error]: Sequence length of new target structure did not coincide with expected size when you called GaRNAcha::inv_fold.\n";
    exit(1);
  }
  if (!yasize && !yaseq) {
    size = target.length();
    yasize = true;
  }
  if (!yaseq)
    rand_seq();
  give_up = 1; //strict solutions only
  double distot = inverse_fold(seq, target.c_str());
  return distot;
}

double GaRNAcha::inv_fold(string target, string inseq) {
  if (yasize && (size != target.length())) {
    cout << "[Error]: Sequence length of new target structure did not coincide with expected size when you called GaRNAcha::inv_fold.\n";
    exit(1);
  }
  if (target.length() != inseq.length()) {
    cout << "[Error]: Target structure length and initial sequence length were different when you called GaRNAcha::inv_fold.\n";
    exit(1);
  }
  if (yaseq) {
    cout << "[Error]: Sequence has been already declared when you called GaRNAcha::inv_fold.\n";
    exit(1);
  }
  if (!yasize && !yaseq) {
    size = target.length();
    yasize = true;
  }
  set_seq_allow_lc(inseq);
  give_up = 1; //strict solutions only
  double distot = inverse_fold(seq, target.c_str());
  return distot;
}

int GaRNAcha::bp_distance(string st1, string st2) {
  return vrna_bp_distance(st1.c_str(), st2.c_str());
}

int GaRNAcha::bp_distance(string st1) {
  if (!yamfe) {
    cout << "[Error]: The MFE structure had not been obtained when you called GaRNAcha::bp_distance.\n";
  }
  return bp_distance(st1.c_str(), mfestr);
}

int GaRNAcha::bp_distance(const char *st1, const char *st2) {
  return vrna_bp_distance(st1, st2);
}

int GaRNAcha::bp_distance(const char *st1) {
  if (!yamfe) {
    cout << "[Error]: The MFE structure had not been obtained when you called GaRNAcha::bp_distance.\n";
  }
  return bp_distance(st1, mfestr);
}

void GaRNAcha::get_mfe_str() {
  if (!yaseq) {
    cout << "[Error]: Sequence did not exist when you called GaRNAcha::get_mfe_str.\n";
    exit(1);
  }
  if (yamfe) {
    cout << "[Error]: MFE had already been obtained when you called GaRNAcha::get_mfe_str.\n";
    exit(1);
  }
  mfestr = new char[size+1];
  mfe = vrna_fold(seq, mfestr);
  yamfe = true;
  return;
}

int GaRNAcha::seq_size() {
  if (!yaseq) {
    cout << "[Error]: Sequence had not been declared when you called GaRNAcha::seq_size.\n";
    exit(1);
  }
  return size;
}

double GaRNAcha::get_temperature() {
  return vrna_md_defaults_temperature_get();
}

void GaRNAcha::set_temperature(double t) {
  vrna_md_defaults_temperature(t);
}

int GaRNAcha::hamming_d_seq(const char *s1, const char *s2)  {
  return vrna_hamming_distance(s1, s2);
}

int GaRNAcha::hamming_d_seq(string s1, string s2)  {
  return vrna_hamming_distance(s1.c_str(), s2.c_str());
}

int GaRNAcha::hamming_d_seq(const char *s1)  {
  if (!yaseq) {
    cout << "[Error]: Sequence had not been declared when you called GaRNAcha::hamming_d_seq.\n";
    exit(1);
  }
  return vrna_hamming_distance(s1, seq);
}

int GaRNAcha::hamming_d_seq(string s1)  {
  if (!yaseq) {
    cout << "[Error]: Sequence had not been declared when you called GaRNAcha::hamming_d_seq.\n";
    exit(1);
  }
  return vrna_hamming_distance(s1.c_str(), seq);
}

double GaRNAcha::tree_d(string stru1, string stru2) { //
  char *xtru;
  xtru = expand_Full(stru1.c_str());
  Tree *beard1 = make_tree(xtru);
  free(xtru);
  xtru = expand_Full(stru2.c_str());
  Tree *beard2 = make_tree(xtru);
  free(xtru);
  double tridi = tree_edit_distance(beard1, beard2);
  free_tree(beard1);
  free_tree(beard2);
  return tridi;
}

double GaRNAcha::tree_d(const char *stru1, const char *stru2) {
  char *xtru;
  xtru = expand_Full(stru1);
  Tree *beard1 = make_tree(xtru);
  free(xtru);
  xtru = expand_Full(stru2);
  Tree *beard2 = make_tree(xtru);
  free(xtru);
  double tridi = tree_edit_distance(beard1, beard2);
  free_tree(beard1);
  free_tree(beard2);
  return tridi;
}

double GaRNAcha::tree_d(string stru1) { //
  if (!yamfe) {
    cout << "[Error]: MFE had not been obtained when you called GaRNAcha::tree_d.\n";
    exit(1);
  }
  char *xtru;
  xtru = expand_Full(stru1.c_str());
  Tree *beard1 = make_tree(xtru);
  free(xtru);
  xtru = expand_Full(mfestr);
  Tree *beard2 = make_tree(xtru);
  free(xtru);
  double tridi = tree_edit_distance(beard1, beard2);
  free_tree(beard1);
  free_tree(beard2);
  return tridi;
}

double GaRNAcha::tree_d(const char *stru1) {
  if (!yamfe) {
    cout << "[Error]: MFE had not been obtained when you called GaRNAcha::tree_d.\n";
    exit(1);
  }
  char *xtru;
  xtru = expand_Full(stru1);
  Tree *beard1 = make_tree(xtru);
  free(xtru);
  xtru = expand_Full(mfestr);
  Tree *beard2 = make_tree(xtru);
  free(xtru);
  double tridi = tree_edit_distance(beard1, beard2);
  free_tree(beard1);
  free_tree(beard2);
  return tridi;
}

void GaRNAcha::print_mfestr(ostream& fs) {
  if (!yamfe) {
    cout << "[Error]: MFE had not been obtained when you called GaRNAcha::print_mfestr.\n";
    exit(1);
  }
  fs << mfestr;
}

double GaRNAcha::eval_structure(string stru) {
  if (!yaseq) {
    cout << "[Error]: Sequence had not been declared when you called GaRNAcha::eval_structure.\n";
    exit(1);
  }
  return vrna_eval_structure_simple(seq, stru.c_str());
}

double GaRNAcha::eval_structure(const char *stru) {
  if (!yaseq) {
    cout << "[Error]: Sequence had not been declared when you called GaRNAcha::eval_structure.\n";
    exit(1);
  }
  return vrna_eval_structure_simple(seq, stru);
}

bool GaRNAcha::is_this_mfestr(string lastr) {
  bool resp = false;
  int i;
  if (!yamfe)
    get_mfe_str();
  if (bp_distance(lastr)==0) {
    resp = true;
  } else {
    if(eval_structure(lastr) != mfe) {
      resp = false;
    } else {
      if (!yaallmfe)
        get_all_mfe_strs();
      resp = false;
      for (i=0; i < allmfes.nufen; i++) {
        if (bp_distance(lastr, an_mfes(i))==0) {
          resp= true;
          break;
        }
      }
    }
  }
  return resp;
}

void GaRNAcha::mutate_seq(int pos, char nuc) {
  if (!yaseq) {
    cout << "[Error]: Sequence had not been declared when you called GaRNAcha::mutate_seq.\n";
    exit(1);
  }
  if ((pos >= size) || (pos < 0)) {
    cout << "[Error]: Nucleotide was outside of sequence range when you called GaRNAcha::mutate_seq.\n";
    exit(1);
  }
  if (yamfe)
    clear_mfes();
  if (yaplas)
    clear_plas_phen();
  if (yaallmfe)
    clear_allmfes();
  if (yamrep)
    clear_mrep();
  if (yamef)
    clear_mef();
  if (yaplr)
    clear_plr_sp();
  if (yacompLR)
    clear_complr();
  if (yacong)
    clear_cong();
  
  yapafu = false;
  seq[pos] = nuc;
}

void GaRNAcha::mutate_seq_fd(int pos) { //force diff
  char oldnt = seq[pos];
  char nwnt;
  do {
    nwnt = alphabet[est.randint(0, 4)];
  } while (nwnt == oldnt);
  mutate_seq(pos, nwnt);
}

void GaRNAcha::rand_mut() {
  int pos = est.randint(0, size);
  char nwnt = alphabet[est.randint(0, 4)];
  mutate_seq(pos, nwnt);
}

void GaRNAcha::rand_mut_fd() {
  int pos = est.randint(0, size);
  mutate_seq_fd(pos);
}

void GaRNAcha::set_e_band(int eb) {
  band = eb;
  yaband = true;
  return;
}

int GaRNAcha::e_band() {
  return band;
}

void GaRNAcha::plastic_phenotypes() {
  if (!yaband) {
    cout << "[Error]: Energy band was not set when you called GaRNAcha::plastic_phenotypes().\n";
    exit(1);
  }
  vrna_fold_compound_t *foco = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
  int abaco = 0;
  vrna_subopt_cb(foco, band, &cue_subopt, (void *)&abaco);
  otipos.nufen = 0;
  otipos.estrus = new string[abaco];
  otipos.eners = new double[abaco];
  vrna_subopt_cb(foco, band, &lle_subopt, (void *)&otipos);
  if (abaco != otipos.nufen) {//quitardespues
    cout << "Ahchirrion\n";
    exit(1);
  }
  otipos.demeama = new int[abaco];
  basic.sort_vector(otipos.eners, otipos.demeama, otipos.nufen, -1000);
  yaplas = true;
  vrna_fold_compound_free(foco);
  return;
}


void GaRNAcha::plastic_phenotypes(int eb) {
  set_e_band(eb);
  plastic_phenotypes();
  return;
}

int GaRNAcha::count_compatible_structures() {
  if (!yaplas) {
    cout << "[Error]: Alternative phenotypes had not been evaluated when you called GaRNAcha::count_compatible_structures().\n";
    exit(1);
  }
  return otipos.nufen;
}

string GaRNAcha::one_compatible_structure(int cual) {
  if (!yaplas) {
    cout << "[Error]: Alternative phenotypes had not been evaluated when you called GaRNAcha::one_compatible_structure().\n";
    exit(1);
  }
  return otipos.estrus[otipos.demeama[cual]];
}


void GaRNAcha::print_phenotypes(ostream& sal) {
  if (!yaplas) {
    cout << "[Error]: Alternative phenotypes had not been evaluated when you called GaRNAcha::print_phenotypes().\n";
    exit(1);
  }
  sal << "Number of alternative phenotypes within a band of " << band/100.0 <<" kcal: " << otipos.nufen << endl;
  int i;
  for (i=0; i < otipos.nufen; i++) {
    sal << i+1 << "\t" << otipos.estrus[otipos.demeama[i]] << "\t" << otipos.eners[otipos.demeama[i]] << endl;
  }
  return;
}

void GaRNAcha::get_all_mfe_strs() {
  vrna_fold_compound_t *foco = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
  int abaco = 0;
  vrna_subopt_cb(foco, 0, &cue_subopt, (void *)&abaco);
  allmfes.nufen = 0;
  allmfes.estrus = new string[abaco];
  allmfes.eners = new double[abaco];
  vrna_subopt_cb(foco, 0, &lle_subopt, (void *)&allmfes);
  if (abaco != allmfes.nufen) {//quitardespues
    cout << "Ahchirrion\n";
    exit(1);
  }
  allmfes.demeama = new int[abaco];
  basic.sort_vector(allmfes.eners, allmfes.demeama, allmfes.nufen, -1000);
  yaallmfe = true;
  vrna_fold_compound_free(foco);
  return;
}

int GaRNAcha::number_of_mfestrs() {
  if (!yaallmfe) {
    cout << "[Error]: All mfe structures had not been evaluated when you called GaRNAcha::number_of_mfestrs().\n";
    exit(1);
  }
  return allmfes.nufen;
}

string GaRNAcha::an_mfes(int cual) {
  if (!yaallmfe) {
    cout << "[Error]: All mfe structures had not been evaluated when you called GaRNAcha::an_mfes().\n";
    exit(1);
  }
  return allmfes.estrus[cual];
}

void GaRNAcha::print_all_mfe_strs(ostream& sal) {
  if (!yaallmfe) {
    cout << "[Error]: All mfe structures had not been evaluated when you called GaRNAcha::print_all_mfe_strs().\n";
    exit(1);
  }
  sal << "Number of mfe structures: " << allmfes.nufen << endl;
  int i;
  for (i=0; i < allmfes.nufen; i++) {
    sal << i+1 << "\t" << allmfes.estrus[allmfes.demeama[i]] << "\t" << allmfes.eners[allmfes.demeama[i]] << endl;
  }
  return;
}

string GaRNAcha::nt_sequence() {
  if (!yaseq) {
    cout << "[Error]: Sequence had not been declared when you called GaRNAcha::nt_sequence.\n";
    exit(1);
  }
  string resul(seq);
  return resul;
}

string GaRNAcha::mfe_structure() {
  if (!yamfe) {
    cout << "[Error]: MFE had not been obtained when you called GaRNAcha::mfe_structure.\n";
    exit(1);
  }
  string resul(mfestr);
  return resul;
}

double GaRNAcha::mfenergy() {
  if (!yamfe) {
    cout << "[Error]: MFE had not been obtained when you called GaRNAcha::mfe_structure.\n";
    exit(1);
  }
  return mfe;
}

double GaRNAcha::m_neighbors_with_mfe_str(string target) {
  if (!yaseq) {
    cout << "[Error]: Sequence did not exist when you called GaRNAcha::m_neighbors_with_mfe_str.\n";
    exit(1);
  }
  GaRNAcha patrab;
  copy_seq(patrab);
  int i, j, abac=0;
  char viej;
  for (i=0; i < size; i++) {
    viej = seq[i];
    for (j=0; j < 4; j++) {
      if (viej != alphabet[j]) {
        patrab.mutate_seq(i, alphabet[j]);
        if (patrab.is_this_mfestr(target))
          abac++;
      }
    }
    patrab.mutate_seq(i, viej);
  }
  patrab.clear();
  return abac/(3.0*size);
}

double GaRNAcha::m_neighbors_with_mfe_str(const char *targ) {
  string target(targ);
  return m_neighbors_with_mfe_str(target);
}

double GaRNAcha::mutational_access(string whstr) {
  if (!yamrep) {
    cout << "[Error]: Mutational repertoire had not been assessed when you called GaRNAcha::mutational_access_from_mrep.\n";
    exit(1);
  }
  int i,j=-1;
  for (i=0; i < mrep_size; i++) {
    if (whstr == mrep[i]) {
      j=i;
      break;
    }
  }
  if (j==-1)
    return 0;
  else
    return mrep_counts[j];
}

void GaRNAcha::m_repertoire() {
  if (!yaseq) {
    cout << "[Error]: Sequence did not exist when you called GaRNAcha::m_repertoire.\n";
    exit(1);
  }
  if (yamrep) {
    cout << "[Error]: Mutational repertoire already existed when you called GaRNAcha::m_repertoire.\n";
    exit(1);
  }
  mrep_counts = new double[size*20];
  mrep = new string[size*20];
  mrep_size = 0;
  GaRNAcha patrab;
  copy_seq(patrab);
  string estra;
  bool eya;
  int i, j, k, l;
  char viej;
  for (i=0; i < size; i++) {
    viej = seq[i];
    for (j=0; j < 4; j++) {
      if (viej != alphabet[j]) {
        patrab.mutate_seq(i, alphabet[j]);
        patrab.get_all_mfe_strs();
        for (l=0; l < patrab.allmfes.nufen; l++) {
          estra = patrab.allmfes.estrus[l];
          eya = false;
          for (k=0; k < mrep_size; k++) {
            if (estra == mrep[k]) {
              mrep_counts[k] += 1.0/patrab.allmfes.nufen;
              eya = true;
              break;
            }
          }
          if (!eya) {
            mrep[mrep_size] = estra;
            mrep_counts[mrep_size] = 1.0/patrab.allmfes.nufen;
            mrep_size++;
            if (mrep_size >= (size*20)) {
              cout << "[Error]: Number of accessible mutant structures exceeds the space allocated for it in GaRNAcha::m_repertoire().\n";
              exit(1);
            }
          }
        }
        patrab.clear_allmfes();
      }
    }
    patrab.mutate_seq(i, viej);
  }
  
  patrab.clear();
  yamrep = true;
  return;
}

bool GaRNAcha::already_m_rep() {
  return yamrep;
}

void GaRNAcha::copy_m_rep(GaRNAcha &froms) {
  if (yamrep) {
    cout << "[Error]: Mutational repertoire had already been assessed when you called GaRNAcha::copy_m_rep.\n";
    exit(1);
  }
  mrep_size = froms.size_of_mrep();
  mrep_counts = new double[mrep_size];
  mrep = new string[mrep_size];
  int i;
  for (i=0; i<mrep_size; i++) {
    mrep[i] = froms.struct_in_mrep(i);
    mrep_counts[i] = froms.count_in_mrep(i);
  }
  yamrep = true;
}


void GaRNAcha::print_mrep(ostream& sal) {
  if (!yamrep) {
    cout << "[Error]: Mutational repertoire did not exist when you called GaRNAcha::print_mrep.\n";
    exit(1);
  }
  int i;
  for (i=0; i < mrep_size ; i++)
    sal << mrep[i] << "\t" << mrep_counts[i] << endl;
  return;
}

double GaRNAcha::count_in_mrep(int cual) {
  if (!yamrep) {
    cout << "[Error]: Mutational repertoire did not exist when you called GaRNAcha::count_in_mrep.\n";
    exit(1);
  }
  return mrep_counts[cual];
}

string GaRNAcha::struct_in_mrep(int cual) {
  if (!yamrep) {
    cout << "[Error]: Mutational repertoire did not exist when you called GaRNAcha::struct_in_mrep.\n";
    exit(1);
  }
  return mrep[cual];
}

int GaRNAcha::size_of_mrep() {
  if (!yamrep) {
    cout << "[Error]: Mutational repertoire did not exist when you called GaRNAcha::size_of_mrep.\n";
    exit(1);
  }
  return mrep_size;
}
//private

void GaRNAcha::partition_function() {
  if (!yaplas) {
    cout << "[Error]: Plastic phenotypes had not been determined when you called GaRNAcha::partition_function.\n";
    exit(1);
  }
  double K = 273.15 + get_temperature();
  double beta = -1.0/(kb*K);
  pafu = 0;
  int i;
  for (i=0; i < otipos.nufen; i++)
    pafu += pow(e, (beta*otipos.eners[otipos.demeama[i]]));
  yapafu = true;
  return;
}

bool GaRNAcha::already_partition_function() {
  return yapafu;
}

double GaRNAcha::p_Boltzmann(int cual) {
  if (!yapafu)
    partition_function();
  double K = 273.15 + get_temperature();
  double beta = -1.0/(kb*K);
  double bf = pow(e, (beta*otipos.eners[otipos.demeama[cual]]));
  return bf/pafu;
}

double GaRNAcha::p_Boltzmann(string unphen) {
  int i, j=-1;
  if (!yaplas)
    plastic_phenotypes();
  for (i=0; i < otipos.nufen; i++)
    if (otipos.estrus[otipos.demeama[i]] == unphen) {
      j= i;
      break;
    }
  if (j==(-1))
    return 0;
  else
    return p_Boltzmann(j);
}

double GaRNAcha::plastic_robustness() {
  return p_Boltzmann(0);
}

void GaRNAcha::get_monomer_eq_freqs() {
  if (yamef) {
    cout << "[Error]: Monomer equilibrium frequencies had already been determined when you called GaRNAcha::get_monomer_eq_freqs.\n";
    exit(1);
  }
  if (!yapafu)
    partition_function();
  int i;
  meqfreqs = new double[otipos.nufen];
  for (i=0; i < otipos.nufen; i++)
    meqfreqs[otipos.demeama[i]] = p_Boltzmann(i);
  yamef = true;
}

double GaRNAcha::return_cff() {
  return cff;
}

void GaRNAcha::set_cff(double ncff) {
  cff = ncff;
  return;
}

double GaRNAcha::w_no_plast(const char *optstr) {
  string oopptt(optstr);
  return w_no_plast(oopptt);
}

double GaRNAcha::w_no_plast(string optstr) {
  if (!yaallmfe)
    get_all_mfe_strs();
  int i;
  double tdmin=3000.0, ttem;
  for (i=0; i < number_of_mfestrs(); i++) {
    ttem = tree_d(allmfes.estrus[i], optstr);
    if (ttem < tdmin)
      tdmin = ttem;
  }
  double ladino = tdmin/(2.0*size);
  return 1.0/(cff+ladino);
}

double GaRNAcha::w_plast_AF(const char *optstr) {
  if (!yaplas)
    plastic_phenotypes();
  if (!yapafu)
    partition_function();
  string optstruc(optstr);
  int i;
  double law = 0;
  double ladino;
  for (i=0; i < otipos.nufen; i++) {
    ladino = tree_d(optstruc, otipos.estrus[otipos.demeama[i]])/(2.0*size);
    law += p_Boltzmann(i)*(1.0/(cff+ladino));
  }
  return law;
}


double GaRNAcha::w_plast_AF(string optstruc) {
  if (!yaplas)
    plastic_phenotypes();
  if (!yapafu)
    partition_function();
  int i;
  double law = 0;
  double ladino;
  for (i=0; i < otipos.nufen; i++) {
    ladino = tree_d(optstruc, otipos.estrus[otipos.demeama[i]])/(2.0*size);
    law += p_Boltzmann(i)*(1.0/(cff+ladino));
  }
  return law;
}


//double GaRNAcha::w_of_struct(const char *estru, const char *optstr) {
//
//}

void GaRNAcha::next_sampled_seq_w(int steps, const char *estru) {
  if (!yamfe)
    get_mfe_str();
  int i;
  string estri(estru);
  if (strcmp(estru, mfestr) != 0) {
    get_all_mfe_strs();
    if (allmfes.nufen == 1) {
      cout << "[Error]: Initial sequence in GaRNAcha::next_sampled_seq does not yield target structure as mfe structure.\n";
      exit(1);
    } else {
      for (i=0; i < allmfes.nufen; i++)
        if (allmfes.estrus[i] == estri)
          break;
      if (i >= allmfes.nufen) {
        cout << "[Error]: Initial sequence in GaRNAcha::next_sampled_seq does not yield target structure as mfe structure.\n";
        exit(1);
      }
    }
    clear_allmfes();
  }
  int cual, j;
  char viejo;
  for (j=0; j < steps; j++) {
    cual = est.randint(0, size);
    viejo = seq[cual];
    mutate_seq_fd(cual);
    get_all_mfe_strs();
    for (i=0; i < allmfes.nufen; i++)
      if (allmfes.estrus[i] == estri)
        break;
    if (i >= allmfes.nufen)
      mutate_seq(cual, viejo);
    else
      clear_allmfes();
  }
  return;
}

void GaRNAcha::next_sampled_seq_w(int steps, string estru) {
  next_sampled_seq_w(steps, estru.c_str());
}

void GaRNAcha::sample_w(ostream& fs, int cuantas, int multip, string estru) {
  int i;
  force_inverse_fold(estru);
  for (i=0; i < cuantas; i++) {
    next_sampled_seq_w(multip*seq_size(), estru);
    print_seq(fs);
    fs << endl;
  }
}

void GaRNAcha::sample_w(ostream& fs, int cuantas, int multip, const char *estru) {
  string laes(estru);
  sample_w(fs, cuantas, multip, laes);
}

string GaRNAcha::phen_rand_mut() {
  GaRNAcha copia;
  copy_seq(copia);
  copia.rand_mut_fd();
  copia.get_all_mfe_strs();
  int i = est.randint(0, copia.allmfes.nufen);
  string strumut = copia.allmfes.estrus[i];
  copia.clear();
  return strumut;
}

string GaRNAcha::phen_rand_mut_with_seq(string &laseq) {
  GaRNAcha copia;
  copy_seq(copia);
  copia.rand_mut_fd();
  copia.get_all_mfe_strs();
  int i = est.randint(0, copia.allmfes.nufen);
  string strumut = copia.allmfes.estrus[i];
  laseq = copia.nt_sequence();
  copia.clear();
  return strumut;
}

int GaRNAcha::maccesible_strucs_from_struc(int vueltas, string esdesde, string *accest, int *contando, int numaxest) { //regresa numero ocupadas
  int ocupadas = 0, i, j;
  string mut;
  bool esta = false;
  for (i=0; i < vueltas; i++) {
    force_inverse_fold(esdesde); //
//    next_sampled_seq(multip*seq_size(), esdesde);
    mut = phen_rand_mut();
    if (mut != esdesde) {
      for (j=0; j < ocupadas; j++) {
        if (mut == accest[j]) {
          esta = true;
          break;
        }
      }
      if (esta)
        contando[j]++;
      else {
        accest[ocupadas] = mut;
        contando[ocupadas] = 1;
        ocupadas++;
      }
    }
    if (ocupadas >= numaxest) {
      cout << "[Error]: The number of accessible structures exceeds the maximum number allowed when you called GaRNAcha::maccesible_strucs_from_struc.\n";
      exit(1);
    }
    clear_seq();
  }
  return ocupadas;
}

int GaRNAcha::maccesible_strucs_from_struc(int vueltas, const char *esdesde, string *accest, int *contando, int numaxest) {
  string ed(esdesde);
  return maccesible_strucs_from_struc(vueltas, ed, accest, contando, numaxest);
}

int GaRNAcha::seqs_neighboring_struct(ostream& fs, string oldst, string newstr, int tries) {
  string mut, laseq;
  int i, van=0;
  for (i=0; i < tries; i++) {
    force_inverse_fold(oldst);
    mut = phen_rand_mut_with_seq(laseq);
    if (mut == newstr) {
      print_seq(fs);
      cout << "\t" << laseq << endl;
      van++;
    }
    clear_seq();
  }
  return van;
}


void GaRNAcha::set_def_plr(double rt, double lt, double kac, double kd_a, double kd_b) {
  paramLR.ka = kac;
  paramLR.Rtot = rt;
  paramLR.Ltot = lt;
  paramLR.kda = kd_a; //0.001
  paramLR.kdb = kd_b; //0.5
  yaplrd = true;
  return;
}

//void GaRNAcha::set_plr(double rt, double lt, double kac, double kd_a, double kd_b, string opt) {//}, double *lasdis, double *pbs) {
//  set_def_plr(rt, lt, kac, kd_a, kd_b);
//  set_spec_plr(opt);
//}

void GaRNAcha::set_spec_plr(string opt) {
  if (!yaplrd) {
    cout << "[Error]: General parameters for the ligand-RNA model had not been set when you called GaRNAcha::set_def_plr.\n";
    exit(1);
  }
  if (yaplr) {
    cout << "[Error]: Parameters for the ligand-RNA model had been set when you called GaRNAcha::set_def_plr.\n";
    exit(1);
  }
  if (!yaplas) {
    cout << "[Error]: Plastic phenotypes had not been assessed when you called GaRNAcha::set_def_plr.\n";
    exit(1);
  }
  if (!yapafu) {
    cout << "[Error]: Partition function had not been assessed when you called GaRNAcha::set_def_plr.\n";
    exit(1);
  }
  paramLR.numec = otipos.nufen;
  paramLR.kd = new double[paramLR.numec];
  paramLR.B = new double[paramLR.numec];
  paramLR.dis = new double[paramLR.numec];
  int i;
  for (i=0; i < paramLR.numec; i++) {
    paramLR.dis[i] = tree_d(otipos.estrus[otipos.demeama[i]], opt);
    paramLR.kd[i] = paramLR.kda*(pow(2, (paramLR.kdb*paramLR.dis[i])));
    paramLR.B[i] = p_Boltzmann(i);
  }
  yaplr = true;
}

bool GaRNAcha::already_spec_plr() {
  return yaplr;
}

double GaRNAcha::get_eqcomp() {
  if (yacompLR) {
    cout << "[Error]: Equilibrium amounts of ligand-RNA complexes had already been determined when you called GaRNAcha::get_eqcomp.\n";
    exit(1);
  }
  if (!yaplr) {
    cout << "[Error]: Specific parameters for the ligand-RNA model had not been set when you called GaRNAcha::get_eqcomp.\n";
    exit(1);
  }
  int i;
  compLRT= 0;
  compLR = new double[paramLR.numec];
  double mbdechich, pamie=0; //-b de cuadratica
  for (i=0; i < paramLR.numec; i++)
    pamie += paramLR.B[i]/paramLR.kd[i];
  mbdechich =  paramLR.Rtot + paramLR.Ltot + 1.0/(paramLR.ka*pamie);
  double dent = (mbdechich*mbdechich) - (4*paramLR.Rtot * paramLR.Ltot);
  compLRT = (mbdechich - sqrt(dent))/2.0;
  double locom = (paramLR.Rtot * paramLR.Ltot) - (compLRT*(paramLR.Rtot + paramLR.Ltot)) + (compLRT*compLRT);
  for (i=0; i < paramLR.numec; i++)
    compLR[i] = locom*paramLR.ka*paramLR.B[i]/paramLR.kd[i];
  double elerr = 0;
  if ((compLRT < 0) || (compLRT > paramLR.Rtot) || (compLRT > paramLR.Ltot))
    elerr = 1;
  yacompLR = true;
  return elerr;
}

double GaRNAcha::get_eqcomp_num() {
  if (yacompLR) {
    cout << "[Error]: Equilibrium amounts of ligand-RNA complexes had already been determined when you called GaRNAcha::get_eqcomp.\n";
    exit(1);
  }
  if (!yaplr) {
    cout << "[Error]: Specific parameters for the ligand-RNA model had not been set when you called GaRNAcha::get_eqcomp.\n";
    exit(1);
  }
  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  int i;
  int status=0;
  size_t iter = 0;
  size_t nuecs = static_cast<size_t>(paramLR.numec);
  gsl_multiroot_function_fdf f = {&lrmod, &lrmod_df, &lrmod_fdf, nuecs, &paramLR};
  double *x0;
  x0 = new double[paramLR.numec];
  if (paramLR.numec < 3) {
    x0[0] = 1;
    for (i=1; i < paramLR.numec; i++)
      x0[i] = 0;
  } else {
    x0[0] = 0.5;
    x0[1] = 0.3;
    x0[2] = 0.2;
    for (i=3; i < paramLR.numec; i++)
      x0[i] = 0;
  }
  
  gsl_vector *x = gsl_vector_alloc (paramLR.numec);
  for (i=0; i < paramLR.numec; i++)
    gsl_vector_set(x, i, x0[i]);
  
  T = gsl_multiroot_fdfsolver_hybridsj; //gnewton;
  s = gsl_multiroot_fdfsolver_alloc (T, paramLR.numec);
  gsl_multiroot_fdfsolver_set (s, &f, x);
  do
  {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate (s);
    if (status)
      break;
    status = gsl_multiroot_test_residual (s->f, 1e-14);
  } while (status == GSL_CONTINUE && iter < 5000);
  if (status) {
    cout << "status = " << gsl_strerror (status) << endl;
    exit(1);
  }
    
  double sumf = 0;
  compLRT= 0;
  compLR = new double[paramLR.numec];
  for (i=0; i < paramLR.numec; i++) {
    compLR[i] = gsl_vector_get(s->x, i);
    compLRT += compLR[i];
    sumf += fabs(gsl_vector_get(s->f, i));
  }

  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free (x);
  delete [] x0;
  yacompLR = true;
  return sumf;
}

double GaRNAcha::get_eqcomp_se() {
  if (yacompLR) {
    cout << "[Error]: Equilibrium amounts of ligand-RNA complexes had already been determined when you called GaRNAcha::get_eqcomp_se.\n";
    exit(1);
  }
  if (!yaplr) {
    cout << "[Error]: Specific parameters for the ligand-RNA model had not been set when you called GaRNAcha::get_eqcomp_se.\n";
    exit(1);
  }
  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  int i;
  int status=0;
  size_t iter = 0;
  size_t nuecs = static_cast<size_t>(paramLR.numec);
  gsl_multiroot_function_fdf f = {&lrmod_se, &lrmod_se_df, &lrmod_se_fdf, nuecs, &paramLR};
  double *x0;
  x0 = new double[paramLR.numec];
  if (paramLR.numec < 3) {
    x0[0] = 1;
    for (i=1; i < paramLR.numec; i++)
      x0[i] = 0;
  } else {
    x0[0] = 0.5;
    x0[1] = 0.3;
    x0[2] = 0.2;
    for (i=3; i < paramLR.numec; i++)
      x0[i] = 0;
  }
  
  gsl_vector *x = gsl_vector_alloc (paramLR.numec);
  for (i=0; i < paramLR.numec; i++)
    gsl_vector_set(x, i, x0[i]);
  
  T = gsl_multiroot_fdfsolver_hybridsj; //gnewton;
  s = gsl_multiroot_fdfsolver_alloc (T, paramLR.numec);
  gsl_multiroot_fdfsolver_set (s, &f, x);
  do
  {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate (s);
    if (status)
      break;
    status = gsl_multiroot_test_residual (s->f, 1e-14);
  } while (status == GSL_CONTINUE && iter < 5000);
  if (status) {
    cout << "status = " << gsl_strerror (status) << endl;
    exit(1);
  }
    
  double sumf = 0;
  compLRT= 0;
  compLR = new double[paramLR.numec];
  for (i=0; i < paramLR.numec; i++) {
    compLR[i] = gsl_vector_get(s->x, i);
    compLRT += compLR[i];
    sumf += fabs(gsl_vector_get(s->f, i));
  }

  gsl_multiroot_fdfsolver_free (s);
  gsl_vector_free (x);
  delete [] x0;
  yacompLR = true;
  return sumf;
}

bool GaRNAcha::already_get_eqcomp() {
  return yacompLR;
}

double GaRNAcha::get_compLRT() {
  if (!yacompLR) {
    cout << "[Error]: Equilibrium amounts of ligand-RNA complexes had not been assessed when you called GaRNAcha::get_compLRT.\n";
    exit(1);
  }
  return compLRT;
}

double GaRNAcha::get_compLR(int i) {
  if (!yacompLR) {
    cout << "[Error]: Equilibrium amounts of ligand-RNA complexes had not been assessed when you called GaRNAcha::get_compLR.\n";
    exit(1);
  }
  if (i >= paramLR.numec) {
    cout << "[Error]: You asked for the equilibrium amount of a ligand-RNA complex that does not exist when you called GaRNAcha::get_compLR.\n";
    exit(1);
  }
  return compLR[i];
}

//esta estaba hecha como para la versión con equilibrio estático. ¿Es necesaria para la nueva versión? Creo que no.
//double GaRNAcha::amount_of_structure(int i) {
//  if (i >= paramLR.numec) {
//    cout << "[Error]: You asked for the total amount of an RNA structure that does not exist when you called GaRNAcha::amount_of_structure.\n";
//    exit(1);
//  }
//  if (!yaplr) {
//    cout << "[Error]: Parameters for the ligand-RNA model had not been set when you called GaRNAcha::amount_of_structure.\n";
//    exit(1);
//  }
//  return paramLR.B[i]*paramLR.Rtot;
//}

double GaRNAcha::dist_to_opt(int i) {
  if (i >= paramLR.numec) {
    cout << "[Error]: You asked for the distance to optimum of an RNA structure that does not exist when you called GaRNAcha::dist_to_opt.\n";
    exit(1);
  }
  if (!yaplr) {
    cout << "[Error]: Parameters for the ligand-RNA model had not been set when you called GaRNAcha::dist_to_opt.\n";
    exit(1);
  }
  return paramLR.dis[i];
}

double GaRNAcha::get_kd(int i) {
  if (i >= paramLR.numec) {
    cout << "[Error]: You asked for the dissociation constant of an RNA structure that does not exist when you called GaRNAcha::get_kd.\n";
    exit(1);
  }
  if (!yaplr) {
    cout << "[Error]: Parameters for the ligand-RNA model had not been set when you called GaRNAcha::get_kd.\n";
    exit(1);
  }
  return paramLR.kd[i];
}

void GaRNAcha::congruence_analysis() {
  if (!yamrep) {
    cout << "[Error]: Mutational repertoire had not been assessed when you called GaRNAcha::congruence_analysis.\n";
    exit(1);
  }
  if (!yaplas) {
    plastic_phenotypes();
//    cout << "[Error]: Plastic repertoire had not been assessed when you called GaRNAcha::congruence_analysis.\n";
//    exit(1);
  }
  if (yacong) {
    cout << "[Error]: Congruence had already been assessed when you called GaRNAcha::congruence_analysis.\n";
    exit(1);
  }
  int min, inter, i,j,l,k;
  if (size_of_mrep() < count_compatible_structures())
    min = size_of_mrep();
  else
    min = count_compatible_structures();
  inter =0;
  string fostru;
  for (i=0; i < size_of_mrep(); i++) {
    fostru = struct_in_mrep(i);
    for (j=0; j < count_compatible_structures(); j++) {
      if (fostru == one_compatible_structure(j)) {
        inter++;
        break;
      }
    }
  }
  jacc = (inter*1.0)/(size_of_mrep() + count_compatible_structures() - (inter*1.0));
  olap =(inter*1.0)/(min*1.0);
  sunion = (size_of_mrep() + count_compatible_structures() - (inter*1.0));
  conmu = new double[sunion];
  conpl = new double[sunion];
  l=0;
  for (j=0; j < size_of_mrep(); j++) {
    fostru = struct_in_mrep(j);
    conmu[l] = count_in_mrep(j);
    conpl[l] = 0;
    for (k=0; k < count_compatible_structures(); k++) {
      if (fostru== one_compatible_structure(k)) {
        conpl[l]=p_Boltzmann(k);
        break;
      }
    }
    l++;
  }
  bool yapa;
  for (k=0; k < count_compatible_structures(); k++) {
    fostru = one_compatible_structure(k);
    yapa = false;
    for (j=0; j < size_of_mrep(); j++) {
      if (fostru==struct_in_mrep(j)) {
        yapa=true;
        break;
      }
    }
    if (!yapa) {
      conmu[l]=0;
      conpl[l]=p_Boltzmann(k);
      l++;
    }
  }
  pear = gsl_stats_correlation(conmu, 1, conpl, 1, sunion);
  yacong = true;
}

bool GaRNAcha::already_cong() {
  return yacong;
}

void GaRNAcha::copy_cong_an(GaRNAcha &froms) {
  if (!froms.already_cong()) {
    cout << "[Error]: Congruence analysis had not been performed when you called GaRNAcha.copy_cong_an.\n";
    exit(1);
  }
  if (yacong) {
    cout << "[Error]: Congruence analysis had already been performed when you called GaRNAcha.copy_cong_an.\n";
    exit(1);
  }
  sunion = froms.size_phen_taccess();
  conmu = new double[sunion];
  conpl = new double[sunion];
  int i;
  for (i=0; i < sunion; i++)
    froms.cong_in_one_phen(i, conmu[i], conpl[i]);
  jacc = froms.jaccard();
  olap = froms.overlap();
  pear = froms.pearson_c();
  yacong = true;
}

void GaRNAcha::cong_in_one_phen(int cual, double &demu, double &depl) {
  if (!yacong) {
    cout << "[Error]: Congruence analysis had not been performed when you called GaRNAcha::cong_in_one_phen.\n";
    exit(1);
  }
  if (cual >= sunion) {
    cout << "[Error]: The phenotype you requested when calling GaRNAcha::cong_in_one_phen does not exist.\n";
    exit(1);
  }
  demu = conmu[cual];
  depl = conpl[cual];
}

double GaRNAcha::jaccard() {
  if (!yacong) {
    cout << "[Error]: Congruence analysis had not been performed when you called GaRNAcha::jaccard.\n";
    exit(1);
  }
  return jacc;
}

double GaRNAcha::overlap() {
  if (!yacong) {
    cout << "[Error]: Congruence analysis had not been performed when you called GaRNAcha::overlap.\n";
    exit(1);
  }
  return olap;
}

double GaRNAcha::pearson_c() {
  if (!yacong) {
    cout << "[Error]: Congruence analysis had not been performed when you called GaRNAcha::pearson_c.\n";
    exit(1);
  }
  return pear;
}

int GaRNAcha::size_phen_taccess() {
  if (!yacong) {
    cout << "[Error]: Congruence analysis had not been performed when you called GaRNAcha::size_phen_taccess.\n";
    exit(1);
  }
  return sunion;
  
}




