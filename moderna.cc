#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
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
}

/*
 install.packages("Rcpp")
 install.packages("roxygen2")
 install.packages("ggcorrplot")
*/
using namespace std;

void report(REvolve &pipol, int model, double *we, double *wafp, double *wnpp, double *wrlp, double *romup, double *pBp, double *jaccp, double *olapp, double *pearp, double *Pmutp, double *Pplap, double *nucdivp, int &nuopp, double &hetp);
void ramble_on(int model, string semi, string sedom[], int contdom[]);
string troncos(int modelin, int perio, string semi, int cortes, double **we, double **wafp, double **wnpp, double **wrlp, double **romup, double **pBp, double **jaccp, double **olapp, double **pearp, double **Pmutp, double **Pplap, double **nucdivp, int *nuopp, double *hetp);
             
  //falta añadir si modelo5 para donde guardar           int rtoto, int ltoto)

int main(int argc, char **argv) {
  /*
   modelos:
   0: ancel-fontana
   1: l-r mucho ligando 20r:80l
   2: l-r igual ligando 20r:20l
   3: l-r poco ligando  20r:5l
   4: sin plasticidad
   5: l-r custom
   */
  //mu 0.001 0.0005 0.9995 M=0.0344
  //mu 0.0001 0.9999 0.993 M=0.00698
  // 0.0002 0.9998 0.9861 M = 0.0139
  
  // moderna semilla cualst model generaciones popsize mu cff ka kda kdb eband rtot ltot thmaxwlr
  // moderna semilla cualst model generaciones popsize mu cff ka kda kdb eband

  
  int peri=1000;
  if ((argc != 15) && (argc != 12)) {
    cout << "[Error]: Wrong number of arguments.\n";
    exit(1);
  }
  int sem = atoi(argv[1]);
  string cual(argv[2]);
  Alea jacta(sem);
  int model = atoi(argv[3]);

  if ((model < 0) || (model >5)) {
    cout << "[Error]: Non-defined model.\n";
    exit(1);
  }
  if (((model==5) && (argc ==12)) || ((model != 5) && (argc==15))) {
    cout << "[Error]: Wrong number of arguments.\n";
    exit(1);
  }
  int gens = atoi(argv[4]);
  int popsize = atoi(argv[5]);
  double mu = atof(argv[6]);
  double cff, rtot, ltot, ka, kda, kdb;
  double thmaxwlr, thmaxwaf, thmaxwnp;
  cff = atof(argv[7]);
  ka = atof(argv[8]);
  kda = atof(argv[9]);
  kdb = atof(argv[10]);
  thmaxwaf = 1.0/cff;
  thmaxwnp = 1.0/cff;
  
  
  int band = atoi(argv[11]);
  if (model==5) {
    rtot = atof(argv[12]);
    ltot = atof(argv[13]);
    thmaxwlr = atof(argv[14]);
  } else {
    rtot = 20;
    if ((model%2)==0) {
      ltot = 20;
      thmaxwlr = 16;
    } else {
      if (model==1) {
        ltot = 80;
        thmaxwlr = 19.67388;
      } else {//if (model == 3) //ya sale sobrando
        ltot = 5;
        thmaxwlr = 4.69338;
      }
    }
  }
  
  Basics bas(jacta);
  string semcue = bas.inttostring(sem)+".txt";
  string optphen, founderseq;
  ifstream fe;
  bas.open_ifstream(fe, "../../../Structures/"+cual+"/"+semcue);
  fe >> optphen;
  fe.close();
  bas.open_ifstream(fe, "../../../Sequences/"+cual+"/"+semcue);
  fe >> founderseq;
  fe.close();
  GaRNAcha fundador(jacta);
  fundador.set_seq(founderseq);
  REvolve villa(jacta);

  villa.set_paramsG(mu, popsize, optphen, band);
  villa.set_paramsAF(cff);
  villa.set_paramsLR(rtot, ltot, ka, kda, kdb);
  villa.start_pop(fundador);

  if (model==4)
    villa.get_pop_wwop();
  else if (model==0) {
    villa.get_pop_waf();
  } else
    villa.get_pop_wlr();
  int i,j,t=0;
  int tp = t/peri;
  int celdas = (gens/peri)+1;
  double **we, **waf, **wnp, **wlr, **romu, **pB, **jacc, **olap, **pear, **Pmut, **Ppla, *het, **nucdiv;
  int *nuop;
  we = new double*[celdas];
  waf = new double*[celdas];
  wnp = new double*[celdas];
  wlr = new double*[celdas];
  romu = new double*[celdas];
  pB = new double*[celdas];
  jacc = new double*[celdas];
  olap = new double*[celdas];
  pear = new double*[celdas];
  Pmut = new double*[celdas];
  Ppla = new double*[celdas];
  nucdiv = new double*[celdas];
  het = new double[celdas];
  nuop = new int[celdas];
  for (i=0; i < celdas; i++) {
    we[i] = new double[3];
    waf[i] = new double[3];
    wnp[i] = new double[3];
    wlr[i] = new double[3];
    romu[i] = new double[3];
    pB[i] = new double[3];
    jacc[i] = new double[3];
    olap[i] = new double[3];
    pear[i] = new double[3];
    Pmut[i] = new double[3];
    Ppla[i] = new double[3];
    nucdiv[i] = new double[3];
  }
  report(villa, model, we[tp], waf[tp], wnp[tp], wlr[tp], romu[tp], pB[tp], jacc[tp], olap[tp], pear[tp], Pmut[tp], Ppla[tp], nucdiv[tp], nuop[tp], het[tp]);
  villa.freeze_population_in_file();
  string sedom[11];
  int contdom[11];
  int ttt;
  
  for (t=1; t<=gens; t++) {
    tp = t/peri;
    if (model==4)
      villa.one_gen_wop();
    else if (model==0) {
      villa.one_gen_waf();
    } else
      villa.one_gen_wlr();
    
    if ((t%100)==0)
      cout << t << endl;
    if ((t%peri)==0) {
      report(villa, model, we[tp], waf[tp], wnp[tp], wlr[tp], romu[tp], pB[tp], jacc[tp], olap[tp], pear[tp], Pmut[tp], Ppla[tp], nucdiv[tp], nuop[tp], het[tp]);

      cout << t << "\t";
      if (model==4)
        cout << villa.meanw_wop() << "\t" << villa.maxw_wop() << endl;
      else if (model==0) {
        cout << villa.meanw_af() << "\t" << villa.maxw_af() << endl;
      } else
        cout << villa.meanw_lr() << "\t" << villa.maxw_lr() << endl;
    }
    if ((t >= gens-1000) && ((t%100)==0)) {
      ttt = (t-(gens-1000))/100;
      sedom[ttt] = villa.dominant_seq(contdom[ttt]);
    }
  }
  for (i=0; i<celdas; i++) {
    for (j=0; j<3; j++) {
      waf[i][j]/= thmaxwaf;
      wlr[i][j]/=thmaxwlr;
      wnp[i][j]/=thmaxwnp;
      if (model==4)
        we[i][j] = wnp[i][j];
      else if (model==0)
        we[i][j] = waf[i][j];
      else
        we[i][j] = wlr[i][j];
    }
  }
  ramble_on(model, semcue, sedom, contdom);
  
  string archfpop = troncos(model, peri, semcue, celdas, we, waf, wnp, wlr, romu, pB, jacc, olap, pear, Pmut, Ppla, nucdiv, nuop, het);
  villa.freeze_population_in_file(archfpop);
  
  fundador.clear();
  villa.clear_pop();
  
  for (i=0; i < celdas; i++) {
    delete [] we[i];
    delete [] waf[i];
    delete [] wnp[i];
    delete [] wlr[i];
    delete [] romu[i];
    delete [] pB[i];
    delete [] jacc[i];
    delete [] olap[i];
    delete [] pear[i];
    delete [] Pmut[i];
    delete [] Ppla[i];
    delete [] nucdiv[i];
  }
  delete [] we;
  delete [] waf;
  delete [] wnp;
  delete [] wlr;
  delete [] romu;
  delete [] pB;
  delete [] jacc;
  delete [] olap;
  delete [] pear;
  delete [] Pmut;
  delete [] Ppla;
  delete [] nucdiv;
  delete [] het;
  delete [] nuop;
  
  return 0;
}

void report(REvolve &pipol, int modelp, double *wep, double *wafp, double *wnpp, double *wlrp, double *romup, double *pBp, double *jaccp, double *olapp, double *pearp, double *Pmutp, double *Pplap, double *nucdivp, int &nuopp, double &hetp) {
  pipol.get_pop_wwop();
  wnpp[0] = pipol.meanw_wop();
  wnpp[1] = pipol.maxw_wop();
  wnpp[2] = pipol.stddev_wwop();
  
  pipol.get_pop_waf();
  wafp[0] = pipol.meanw_af();
  wafp[1] = pipol.maxw_af();
  wafp[2] = pipol.stddev_waf();
  

  pipol.get_pop_wlr();
  

  wlrp[0] = pipol.meanw_lr();
  wlrp[1] = pipol.maxw_lr();
  wlrp[2] = pipol.stddev_wlr();
  int i;
  if (modelp==0) {
    for (i=0; i<3; i++)
      wep[i] = wafp[i];
  } else if (modelp==4) {
    for (i=0; i<3; i++)
      wep[i] = wnpp[i];
  } else {
    for (i=0; i<3; i++)
      wep[i] = wlrp[i];
  }
  pipol.Rmut_in_pop(romup[0], romup[1], romup[2]);
  pipol.m_rep_in_pop();
  pipol.Boltzmann_in_pop(pBp[0], pBp[1], pBp[2]);

  pipol.congruence_in_pop();
  pipol.Jaccard_in_pop(jaccp[0], jaccp[1], jaccp[2]);
  pipol.overlap_in_pop(olapp[0], olapp[1], olapp[2]);
  pipol.Pearson_in_pop(pearp[0], pearp[1], pearp[2]);
  pipol.mutational_access_in_pop(Pmutp[0], Pmutp[1], Pmutp[2]);
  pipol.plasticity_access_in_pop(Pplap[0], Pplap[1], Pplap[2]);
  pipol.nucleotide_div(nucdivp[0], nucdivp[1], nucdivp[2]);
  nuopp= pipol.inds_with_opt_mfe();
  hetp = pipol.heterozygosity();
}

void ramble_on(int model, string semi, string sedom[], int contdom[]) {
  string dirs[5] = {"AF_0", "LRmL_1", "LRiL_2", "LRpL_3", "NP_4"};
  string eldir = dirs[model];
  string archi = eldir+"/Results/DomSeq/"+semi;
  int i;
  ofstream fs;
  fs.open(archi.c_str());
  for (i=0; i < 11 ; i++)
    fs << (i*100)-1000 << "\t" << contdom[i] << "\t" << sedom[i] << endl;
  return;
}

string troncos(int modelin, int perio, string semi, int cortes, double **we, double **wafp, double **wnpp, double **wlrp, double **romup, double **pBp, double **jaccp, double **olapp, double **pearp, double **Pmutp, double **Pplap, double **nucdivp, int *nuopp, double *hetp) {
  string dirs[5] = {"AF_0", "LRmL_1", "LRiL_2", "LRpL_3", "NP_4"};
  string eldir = dirs[modelin];
  //string clares[14] = {"W_e", "W_af", "W_np", "W_lr", "R_mu", "pB", "Jaccard", "Olap", "Pearson", "P_mu", "P_pl", "NucDiv", "Hetzyg", "NuOpt"};
  int i,j;
  ofstream fs;
  string archi;
  
  archi = eldir+"/Results/W_e/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << we[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/W_af/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << wafp[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/W_np/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << wnpp[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/W_lr/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << wlrp[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/R_mu/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << romup[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/pB/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << pBp[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/Jaccard/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << jaccp[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/Olap/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << olapp[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/Pearson/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << pearp[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/P_mu/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << Pmutp[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/P_pl/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << Pplap[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/NucDiv/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++) {
    fs << i*perio ;
    for (j=0; j<3; j++)
      fs << "\t" << nucdivp[i][j];
    fs << endl;
  }
  fs.close();
  
  archi = eldir+"/Results/Hetzyg/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++)
    fs << i*perio << "\t" << hetp[i] << endl;
  fs.close();
  
  archi = eldir+"/Results/NuOpt/"+semi;
  fs.open(archi.c_str());
  for (i=0; i < cortes; i++)
    fs << i*perio << "\t" << nuopp[i] << endl;
  fs.close();
  
  return eldir+"/Results/FinalPop/"+semi;
}
