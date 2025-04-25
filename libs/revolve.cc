#include <cstdlib>
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
#include "garnacha.h"
#include "revolve.h"

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

REvolve::REvolve() {
  set_default_parameters();
}

REvolve::REvolve(Alea& jacta) {
  set_default_parameters();
  start_rng(jacta);
}

void REvolve::set_default_parameters() {
  yapG = false;
  yapAF = false;
  yapLR = false;
  yapop = false;
}


void REvolve::start_rng(Alea& jacta) {
  est = jacta;
  basic.start_rng(est);
}

void REvolve::set_paramsG(double muratep, int popsizep, string estru, int eb) {
  mu = muratep;
  popsize = popsizep;
  set_optimum(estru);
  eband = eb;
  yapG = true;
}

void REvolve::set_optimum(string estru) {
  optist = estru;
}

void REvolve::set_paramsAF(double cefefe) {
  cff = cefefe;
  yapAF = true;
}

void REvolve::set_paramsLR(double rt, double lt, double kac, double kd_a, double kd_b) {
  rtot = rt;
  ltot = lt;
  ka = kac;
  kda = kd_a;
  kdb = kd_b;
  yapLR = true;
}

void REvolve::start_pop(GaRNAcha &founder) {
  if (yapop) {
    cout << "[Error]: Population had already been set when you called REvolve::start_pop.\n";
    exit(1);
  }
  if (!yapG) {
    cout << "[Error]: Parameters had not been set when you called REvolve::start_pop().\n";
    exit(1);
  }
  seqsize = founder.seq_size();
  pop = new GaRNAcha[popsize];
  nepop = new GaRNAcha[popsize];
  wwop = new double[popsize];
  waf = new double[popsize];
  wlr = new double[popsize];
  basic.fillv0(wwop, popsize);
  basic.fillv0(waf, popsize);
  basic.fillv0(wlr, popsize);
  meanwwop =0;
  maxwwop = 0;
  meanwaf=0;
  meanwlr=0;
  maxwaf=0;
  maxwlr=0;
  minwwop=0;
  minwlr=0;
  minwaf=0;
  int i;
  for (i=0; i<popsize; i++) {
    pop[i].start_rng(est);
    nepop[i].start_rng(est);
    founder.copy_seq(pop[i]);
  }
  mutation_to_pop();
  yapop = true;
}

void REvolve::clear_pop() {
  delete [] pop;
  delete [] nepop;
  delete [] wwop;
  delete [] waf;
  delete [] wlr;
  yapop=false;
}

void REvolve::mutation_to_pop() {
  int i,numut;
  set<int> quemuts;
  set<int>::iterator it;
  for (i=0; i<popsize; i++) {
    numut = est.ran_binomial(mu, seqsize);
    quemuts.clear();
    while (quemuts.size() < numut)
      quemuts.insert(est.randint(0,seqsize));
    for (it = quemuts.begin(); it != quemuts.end(); ++it)
      pop[i].mutate_seq_fd(*it);
  }
  quemuts.clear();
}

void REvolve::get_pop_wwop() {
  int i, j;
  sumwwop=0;
  maxwwop=0;
  minwwop=1000;
  bool antes;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j <i ; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        wwop[i] = wwop[j];
        break;
      }
    }
    if (!antes)
      wwop[i] = pop[i].w_no_plast(optist);
    sumwwop += wwop[i];
    if (wwop[i] > maxwwop)
      maxwwop = wwop[i];
    if (wwop[i] < minwwop)
      minwwop = wwop[i];
  }
  meanwwop = sumwwop/(popsize*1.0);
}

double REvolve::stddev_wwop() {
  return basic.get_pop_stddev(wwop, popsize, meanwwop);
}

void REvolve::get_pop_waf() {
  if (!yapop) {
    cout << "[Error]: Population had not been set when you called REvolve::get_pop_waf.\n";
    exit(1);
  }
  if (!yapAF) {
    cout << "[Error]: Parameters had not been set when you called REvolve::get_pop_waf().\n";
    exit(1);
  }
  int i, j;
  sumwaf=0;
  maxwaf=0;
  minwaf=1000;
  bool antes;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j <i ; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        waf[i] = waf[j];
      }
    }
    if (!antes) {
      pop[i].set_e_band(eband);
      pop[i].plastic_phenotypes();
      pop[i].partition_function();
      pop[i].set_cff(cff);
      waf[i] = pop[i].w_plast_AF(optist);
      pop[i].clear_plas_phen();
    }
    sumwaf += waf[i];
    if (waf[i] > maxwaf)
      maxwaf = waf[i];
    if (waf[i] < minwaf)
      minwaf = waf[i];
  }
  meanwaf = sumwaf/(popsize*1.0);
}

double REvolve::stddev_waf() {
  return basic.get_pop_stddev(waf, popsize, meanwaf);
}

void REvolve::get_pop_wlr() {
  if (!yapop) {
    cout << "[Error]: Population had not been set when you called REvolve::get_pop_wlr.\n";
    exit(1);
  }
  if (!yapLR) {
    cout << "[Error]: Parameters had not been set when you called REvolve::get_pop_wlr().\n";
    exit(1);
  }
  int i, j;
  sumwlr=0;
  maxwlr=0;
  minwlr=1000;
  bool antes;
  double mistakos;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j <i ; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        wlr[i] = wlr[j];
      }
    }
    if (!antes) {
        if (!pop[i].already_plastic_phen()) {
          pop[i].set_e_band(eband);
          pop[i]. plastic_phenotypes();
        }
      if (!pop[i].already_partition_function())
        pop[i].partition_function();
      pop[i].set_def_plr(rtot, ltot, ka, kda, kdb);
      if (!pop[i].already_spec_plr())
        pop[i].set_spec_plr(optist);
      if (!pop[i].already_get_eqcomp()) {
        mistakos = pop[i].get_eqcomp();
        if (mistakos > 0.05) {
          cout << "[Error]: method did not converge when you called REvolve::get_pop_wlr.\n";
          exit(1);
        }
      }
      wlr[i] = pop[i].get_compLRT();
      pop[i].clear_plas_phen(); //new line to avoid leak
    }
    sumwlr += wlr[i];
    if (wlr[i] > maxwlr)
      maxwlr = wlr[i];
    if (wlr[i] < minwlr)
      minwlr = wlr[i];
  }
  meanwlr = sumwlr/(popsize*1.0);
}

void REvolve::get_pop_wlr_se() {
  if (!yapop) {
    cout << "[Error]: Population had not been set when you called REvolve::get_pop_wlr_se.\n";
    exit(1);
  }
  if (!yapLR) {
    cout << "[Error]: Parameters had not been set when you called REvolve::get_pop_wlr_se.\n";
    exit(1);
  }
  int i, j;
  sumwlr=0;
  maxwlr=0;
  minwlr=1000;
  bool antes;
  double mistakos;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j <i ; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        wlr[i] = wlr[j];
      }
    }
    if (!antes) {
        if (!pop[i].already_plastic_phen()) {
          pop[i].set_e_band(eband);
          pop[i]. plastic_phenotypes();
        }
      if (!pop[i].already_partition_function())
        pop[i].partition_function();
      pop[i].set_def_plr(rtot, ltot, ka, kda, kdb);
      if (!pop[i].already_spec_plr())
        pop[i].set_spec_plr(optist);
      if (!pop[i].already_get_eqcomp()) {
        mistakos = pop[i].get_eqcomp_se();
        if (mistakos > 0.05) {
          cout << "[Error]: method did not converge when you called REvolve::get_pop_wlr_se.\n";
          exit(1);
        }
      }
      wlr[i] = pop[i].get_compLRT();
      pop[i].clear_plas_phen(); //new line to avoid leak
    }
    sumwlr += wlr[i];
    if (wlr[i] > maxwlr)
      maxwlr = wlr[i];
    if (wlr[i] < minwlr)
      minwlr = wlr[i];
  }
  meanwlr = sumwlr/(popsize*1.0);
}

double REvolve::stddev_wlr() {
  return basic.get_pop_stddev(wlr, popsize, meanwlr);
}

void REvolve::one_gen_wop() {
  int i,j;
  double x, van;
  //parametros individuales
  //copiar a nueva gen
  
  for (i=0; i < popsize; i++) {
    x = est.randreal()*sumwwop;
    van=0;
    for (j=0; j<popsize; j++) {
      van += wwop[j];
      if (x <= van)
        break;
    }
    pop[j].copy_seq(nepop[i]);
  }
  for (i=0; i < popsize; i++) {
    pop[i].clear_seq(); //newline to avoid memory leak
    nepop[i].copy_seq(pop[i]);
    nepop[i].clear_seq();
  }
  mutation_to_pop();
  //sacar adecuaciones
  get_pop_wwop();
  
}

void REvolve::one_gen_waf() {
  int i,j;
  double x, van;
  //parametros individuales
  //copiar a nueva gen
  
  for (i=0; i < popsize; i++) {
    x = est.randreal()*sumwaf;
    van=0;
    for (j=0; j<popsize; j++) {
      van += waf[j];
      if (x <= van)
        break;
    }
    pop[j].copy_seq(nepop[i]);
  }
  for (i=0; i < popsize; i++) {
    pop[i].clear_seq(); //newline to avoid memory leak
    nepop[i].copy_seq(pop[i]);
    nepop[i].clear_seq();
  }
  mutation_to_pop();
  //sacar adecuaciones
  get_pop_waf();
}

void REvolve::one_gen_wlr() {
  int i,j;
  double x, van;
  //parametros individuales
  //copiar a nueva gen
  
  for (i=0; i < popsize; i++) {
    x = est.randreal()*sumwlr;
    van=0;
    for (j=0; j<popsize; j++) {
      van += wlr[j];
      if (x <= van)
        break;
    }
    pop[j].copy_seq(nepop[i]);
  }
  for (i=0; i < popsize; i++) {
    pop[i].clear_seq(); //newline to avoid memory leak
    nepop[i].copy_seq(pop[i]);
    nepop[i].clear_seq();
  }
  mutation_to_pop();
  //sacar adecuaciones
  get_pop_wlr();
}

double REvolve::meanw_wop() {
  return meanwwop;
}

double REvolve::meanw_af() {
  return meanwaf;
}

double REvolve::meanw_lr() {
  return meanwlr;
}

double REvolve::maxw_wop() {
  return maxwwop;
}

double REvolve::maxw_af() {
  return maxwaf;
}

double REvolve::maxw_lr() {
  return maxwlr;
}

void REvolve::Rmut_in_pop(double &mean, double &max, double &stddev) {
  int i,j;
  double *vectemp;
  vectemp = new double[popsize];
  max = -1;
  bool antes;
    
  for (i=0; i < popsize; i++) {
         antes = false;
     for (j=0; j < i; j++) {
       if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        vectemp[i] = vectemp[j];
        break;
      }
    }
    if (!antes) {
      vectemp[i] = pop[i].m_neighbors_with_mfe_str(optist);
      if (vectemp[i] > max)
        max = vectemp[i];
    }
  }
  mean = basic.get_mean(vectemp, popsize);
  stddev = basic.get_pop_stddev(vectemp, popsize, mean);
   
  delete [] vectemp;
}

void REvolve::Boltzmann_in_pop(double &mean, double &max, double &stddev) {
  int i,j;
  double *vectemp;
  vectemp = new double[popsize];
  max = -1;
  bool antes;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j < i; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        vectemp[i] = vectemp[j];
        break;
      }
    }
    if (!antes) {
      vectemp[i] = pop[i].p_Boltzmann(optist);
      if (vectemp[i] > max)
        max = vectemp[i];
    }
  }
  mean = basic.get_mean(vectemp, popsize);
  stddev = basic.get_pop_stddev(vectemp, popsize, mean);
  delete [] vectemp;
}

void REvolve::Jaccard_in_pop(double &mean, double &max, double &stddev) {
  int i,j;
  double *vectemp;
  vectemp = new double[popsize];
  max = -1;
  bool antes;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j < i; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        vectemp[i] = vectemp[j];
        break;
      }
    }
    if (!antes) {
      vectemp[i] = pop[i].jaccard();
      if (vectemp[i] > max)
        max = vectemp[i];
    }
  }
  mean = basic.get_mean(vectemp, popsize);
  stddev = basic.get_pop_stddev(vectemp, popsize, mean);
  delete [] vectemp;
}


void REvolve::overlap_in_pop(double &mean, double &max, double &stddev) {
  int i,j;
  double *vectemp;
  vectemp = new double[popsize];
  max = -1;
  bool antes;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j < i; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        vectemp[i] = vectemp[j];
        break;
      }
    }
    if (!antes) {
      vectemp[i] = pop[i].overlap();
      if (vectemp[i] > max)
        max = vectemp[i];
    }
  }
  mean = basic.get_mean(vectemp, popsize);
  stddev = basic.get_pop_stddev(vectemp, popsize, mean);
  delete [] vectemp;
}

void REvolve::Pearson_in_pop(double &mean, double &max, double &stddev) {
  int i,j;
  double *vectemp;
  vectemp = new double[popsize];
  max = -1;
  bool antes;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j < i; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        vectemp[i] = vectemp[j];
        break;
      }
    }
    if (!antes) {
      vectemp[i] = pop[i].pearson_c();
      if (vectemp[i] > max)
        max = vectemp[i];
    }
  }
  mean = basic.get_mean(vectemp, popsize);
  stddev = basic.get_pop_stddev(vectemp, popsize, mean);
  delete [] vectemp;
}

void REvolve::mutational_access_in_pop(double &mean, double &max, double &stddev) {
  int i,j;
  double *vectemp;
  vectemp = new double[popsize];
  max = -1;
  bool antes;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j < i; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        vectemp[i] = vectemp[j];
        break;
      }
    }
    if (!antes) {
      vectemp[i] = pop[i].size_of_mrep();
      if (vectemp[i] > max)
        max = vectemp[i];
    }
  }
  mean = basic.get_mean(vectemp, popsize);
  stddev = basic.get_pop_stddev(vectemp, popsize, mean);
  delete [] vectemp;
}

void REvolve::plasticity_access_in_pop(double &mean, double &max, double &stddev) {
  int i,j;
  double *vectemp;
  vectemp = new double[popsize];
  max = -1;
  bool antes;
  for (i=0; i < popsize; i++) {
    antes = false;
    for (j=0; j < i; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        vectemp[i] = vectemp[j];
        break;
      }
    }
    if (!antes) {
      vectemp[i] = pop[i].count_compatible_structures();
      if (vectemp[i] > max)
        max = vectemp[i];
    }
  }
  mean = basic.get_mean(vectemp, popsize);
  stddev = basic.get_pop_stddev(vectemp, popsize, mean);
  delete [] vectemp;
}

void REvolve::nucleotide_div(double &mean, double &max, double &stddev) {
  int i, j, cont=0;
  max = -1;
  double mien;
  mean=0;
  for (i=0; i < (popsize-1); i++) {
    for (j=i+1; j < popsize; j++) {
      mien = pop[i].hamming_d_seq(pop[j].nt_sequence());
      if (mien > max)
        max = mien;
      mean += mien;
      cont++;
    }
  }
  mean /= (cont*1.0);
  stddev = 0;
  for (i=0; i < (popsize-1); i++) {
    for (j=i+1; j < popsize; j++) {
      mien = pop[i].hamming_d_seq(pop[j].nt_sequence());
      stddev+= ((mien-mean)*(mien-mean));
    }
  }
  stddev /= (cont*1.0);
  stddev = sqrt(stddev);
}

double REvolve::heterozygosity() {
  string *mientras;
  int ocu = 0, i,j;
  int *conts;
  double *frecs;
  bool yasta;
  mientras = new string[popsize];
  frecs = new double[popsize];
  conts = new int[popsize];
  basic.fillv0(conts, popsize);
  for (i=0; i < popsize; i++) {
    yasta = false;
    for (j=0; j < ocu; j++)
      if (pop[i].nt_sequence() == mientras[j]) {
        conts[j]++;
        yasta = true;
        break;
      }
    if (!yasta) {
      mientras[ocu] = pop[i].nt_sequence();
      conts[ocu] = 1;
      ocu++;
    }
  }
  for (i=0; i < ocu; i++)
    frecs[i] = conts[i]/(popsize*1.0);
  double het = 1;
  for (i=0; i < ocu; i++)
    het -= (frecs[i]*frecs[i]);
  delete [] mientras;
  delete [] frecs;
  delete [] conts;
  return het;
}

int REvolve::inds_with_opt_mfe() {
  int i, cuenta=0;
  for (i=0; i < popsize; i++) {
    if (pop[i].is_this_mfestr(optist))
      cuenta++;
  }
  return cuenta;
}

string REvolve::dominant_seq(int &hmany) {
  if (!yapop) {
    cout << "[Error]: Population had not been set when you called REvolve::dominant_seq.\n";
    exit(1);
  }
  string *lyapupob;
  lyapupob = new string[popsize];
  int *con;
  con = new int[popsize];
  int i,j, ocu = 0;
  bool ya;
  for (i=0; i < popsize; i++) {
    ya = false;
    for (j=0; j < ocu; j++) {
      if (pop[i].nt_sequence() == lyapupob[j]) {
        ya = true;
        con[j]++;
        break;
      }
    }
    if (!ya) {
      lyapupob[ocu] = pop[i].nt_sequence();
      con[ocu] = 1;
      ocu++;
    }
  }
  int mass=0,imass=-1;
  for (i=0; i < ocu; i++) {
    if (con[i] > mass) {
      mass = con[i];
      imass = i;
    }
  }
  string domi = lyapupob[imass];
  hmany = mass;
  delete [] lyapupob;
  delete [] con;
  return domi;
}


void REvolve::freeze_population_in_file(string nombre) {
  if (!yapop) {
    cout << "[Error]: Population had not been set when you called REvolve::freeze_population_in_file.\n";
    exit(1);
  }
  string *lyapupob;
  lyapupob = new string[popsize];
  int *con;
  con = new int[popsize];
  int i,j, ocu = 0;
  bool ya;
  for (i=0; i < popsize; i++) {
    ya = false;
    for (j=0; j < ocu; j++) {
      if (pop[i].nt_sequence() == lyapupob[j]) {
        ya = true;
        con[j]++;
        break;
      }
    }
    if (!ya) {
      lyapupob[ocu] = pop[i].nt_sequence();
      con[ocu] = 1;
      ocu++;
    }
  }
  ofstream fs;
  fs.open(nombre.c_str());
  fs << ocu << endl << endl;
  for (i=0; i < ocu; i++)
    fs << i+1 << "\t" << lyapupob[i] << "\t" << con[i] << endl;
  fs << endl << endl;
  
  for (i=0; i < ocu; i++) {
    for (j=0; j < ocu; j++)
      fs << pop[i].hamming_d_seq(lyapupob[i], lyapupob[j]) << "\t";
    fs << endl;
  }
  fs.close();
  delete [] lyapupob;
  delete [] con;
}

void REvolve::freeze_population_in_file() {
  if (!yapop) {
    cout << "[Error]: Population had not been set when you called REvolve::freeze_population_in_file.\n";
    exit(1);
  }
  string *lyapupob;
  lyapupob = new string[popsize];
  int *con;
  con = new int[popsize];
  int i,j, ocu = 0;
  bool ya;
  for (i=0; i < popsize; i++) {
    ya = false;
    for (j=0; j < ocu; j++) {
      if (pop[i].nt_sequence() == lyapupob[j]) {
        ya = true;
        con[j]++;
        break;
      }
    }
    if (!ya) {
      lyapupob[ocu] = pop[i].nt_sequence();
      con[ocu] = 1;
      ocu++;
    }
  }
  cout << ocu << endl << endl;
  for (i=0; i < ocu; i++)
    cout << i+1 << "\t" << lyapupob[i] << "\t" << con[i] << endl;
  cout << endl << endl;
  
  for (i=0; i < ocu; i++) {
    for (j=0; j < ocu; j++)
      cout << pop[i].hamming_d_seq(lyapupob[i], lyapupob[j]) << "\t";
    cout << endl;
  }
  delete [] lyapupob;
  delete [] con;
}

void REvolve::congruence_in_pop() {
  int i,j;
  bool antes;
  for (i=0; i<popsize; i++) {
      antes = false;
    for (j=0; j < i; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        antes = true;
        pop[i].copy_cong_an(pop[j]);
        break;
      }
    }
    if (!antes)
      pop[i].congruence_analysis();
  }
  
}

void REvolve::m_rep_in_pop() {
  int i,j;
  bool ero;
  for (i=0; i<popsize; i++) {
    ero = false;
    for (j=0; j < i; j++) {
      if (pop[i].nt_sequence() == pop[j].nt_sequence()) {
        pop[i].copy_m_rep(pop[j]);
        ero = true;
        break;
      }
    }
    if (!ero)
      pop[i].m_repertoire();
  }
}
