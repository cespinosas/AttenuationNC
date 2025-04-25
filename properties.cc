#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
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
}

/*
 install.packages("Rcpp")
 install.packages("roxygen2")
 install.packages("ggcorrplot")
*/
using namespace std;

int main(int argc, char **argv) {
  
  if (argc != 4) {
    cout << "[Error]: Wrong number of arguments.\n";
    exit(1);
  }
  string cual(argv[2]);
  Alea jacta(atoi(argv[1]));
  int cuantas = atoi(argv[3]);
  Basics bas(jacta);
  ifstream fe;
  bas.open_ifstream(fe, "../Structures/"+cual+".txt");
  string tau, useq;
  fe >> tau;
  fe.close();
  int taulen = tau.length();
  
//  int parcomb;
  double rtot=20; //cambié el orden
  double ltot[3] ={5,20,80};
  double maxwrl[3] = {4.69338, 16, 19.67388};
  double maxwaf = 100;
  //string losti[3] = {"20:5 RNA:ligand ratio", "20:20 RNA:ligand ratio", "20:80 RNA:ligand ratio"};

  ofstream fs;
  bas.run_command("mkdir Results/"+cual);

  bas.run_command("mkdir Results/"+cual+"/DetailsCounts");
  double **wrl, *waf, *pb, *romut, *fenacmut, *fenacplas, *jaccardmp, *overlapcomp, *pearsonmp, *mdfsuboptoopt, *ddfsuboptoopt;//*spearmanmp;
  int i,j,k,l,ii;

  wrl = new double*[3];
  for (i=0; i<3; i++)
    wrl[i] = new double[cuantas];
  waf = new double[cuantas]; //
  pb = new double[cuantas]; //
  romut = new double[cuantas]; //
  fenacmut = new double[cuantas]; //
  fenacplas = new double[cuantas]; //
  jaccardmp  = new double[cuantas]; //
  overlapcomp = new double[cuantas]; //
  pearsonmp = new double[cuantas];
  mdfsuboptoopt = new double[cuantas];
  ddfsuboptoopt = new double[cuantas];
  
  bool yapa;
  double r, pv;
  int inter, min, onion;
  double **conmu, **conpl;
  conmu = new double*[cuantas];
  conpl = new double*[cuantas];
  
  string strumu;
  GaRNAcha *sope;//(jacta);
  sope = new GaRNAcha[cuantas];
  bas.open_ifstream(fe, "../Sample/Sequences/"+cual+".tex");
  
  int *ditau;
  double *pbditau;
  double pbsotot;
  int fti;
  
  for (i=0; i < cuantas; i++) {
    fe >> useq;
    sope[i].set_default_parameters();
    sope[i].start_rng(jacta);
    sope[i].set_size(taulen);
    sope[i].set_seq(useq);
    sope[i].set_e_band(200);
    sope[i].plastic_phenotypes();
    sope[i].partition_function();
    fenacplas[i] = sope[i].count_compatible_structures();
    pb[i] = sope[i].plastic_robustness();
    romut[i] =sope[i].m_neighbors_with_mfe_str(tau);
    sope[i].m_repertoire();
    fenacmut[i] = sope[i].size_of_mrep();
    waf[i] = sope[i].w_plast_AF(tau)/maxwaf;
    sope[i].congruence_analysis();
    jaccardmp[i] = sope[i].jaccard();
    overlapcomp[i] = sope[i].overlap();
    onion = sope[i].size_phen_taccess();
    fti = fenacplas[i]-1;
    ditau = new int[fti];
    pbditau = new double[fti];
    k=0;
    pbsotot=0;
    for (j=0; j < fenacplas[i]; j++) {
      if (sope[i].one_compatible_structure(j) != tau) {
        ditau[k] = sope[i].tree_d(tau, sope[i].one_compatible_structure(j));
        pbditau[k] = sope[i].p_Boltzmann(j);
        pbsotot += sope[i].p_Boltzmann(j);
        k++;
      }
    }
    if (k != fti) {
      cout << "[Error]: Wrong number of plasticity accessible phenotypes.\n";
      exit(1);
    }
    for (j=0; j < fti; j++)
      pbditau[j] /= pbsotot;
    mdfsuboptoopt[i] = 0;
    for (j=0; j < fti; j++)
      mdfsuboptoopt[i] += (ditau[j]*pbditau[j]);
    ddfsuboptoopt[i] = 0;
    for (j=0; j < fti; j++)
      ddfsuboptoopt[i] += (pbditau[j]*(ditau[j]- mdfsuboptoopt[i])*(ditau[j]- mdfsuboptoopt[i]));
    ddfsuboptoopt[i] = sqrt(ddfsuboptoopt[i]);
    delete [] ditau;
    delete [] pbditau;
    //mdfsuboptoopt
    conmu[i] = new double[onion];
    conpl[i] = new double[onion];
    l=0;
    for (j=0; j < fenacmut[i]; j++) {
      strumu = sope[i].struct_in_mrep(j);
      conmu[i][l] = sope[i].count_in_mrep(j);
      conpl[i][l] = 0;
      for (k=0; k< fenacplas[i];k++) {
        if (strumu==sope[i].one_compatible_structure(k)) {
          conpl[i][l] = sope[i].p_Boltzmann(k);
          break;
        }
      }
      l++;
    }
    for (k=0; k< fenacplas[i];k++) {
      strumu = sope[i].one_compatible_structure(k);
      yapa = false;
      for (j=0; j < (fenacmut[i]); j++) {
        if (strumu==sope[i].struct_in_mrep(j)) {
          yapa = true;
          break;
        }
      }
      if (!yapa) {
        conmu[i][l] = 0;
        conpl[i][l] = sope[i].p_Boltzmann(k);
        l++;
      }
    }
    bas.open_ofstream(fs, "Results/"+cual+"/DetailsCounts/"+bas.inttostring(i)+".txt");
    fs << "Plastic\tMutations\n";
    for (j=0; j<onion; j++)
      fs << conpl[i][j] << "\t" << conmu[i][j] << endl;
    fs.close();
    pearsonmp[i] = sope[i].pearson_c();
    for (ii=0; ii<3; ii++) {
      sope[i].set_def_plr(rtot, ltot[ii], 0.001, 0.001, 0.5); // pero igual con otro ltot
      sope[i].set_spec_plr(tau); //pero igual
      pv = sope[i].get_eqcomp(); //spec
      wrl[ii][i] = sope[i].get_compLRT()/maxwrl[ii];// //con maxwrl espe
      sope[i].clear_complr();
      sope[i].clear_plr_sp();
    }
    if ((i != 0) && ((i%1500)==0))
      bas.run_command("sleep 4m");   
  }
  fe.close();
  
  
  bas.open_ofstream(fs, "Results/"+cual+"/DetailsCounts/paR.sh");
  fs << "library(ggplot2)\n";
  for (i=0; i<cuantas; i++) {
    fs << "a <- read.csv(file = \"Results/"+cual+"/DetailsCounts/" << bas.inttostring(i) << ".txt\", head=TRUE, sep=\"\\t\")\n";
    fs << "pdf(file=\"Results/"+cual+"/DetailsCounts/ScatPlasMut_" << i << ".pdf\")\n";
    fs << "ggplot(a, aes(x=Plastic, y=Mutations))+ geom_point(size=6, shape=21, alpha=0.08, colour=\"black\", fill=\"gray\", stroke=2) + ggtitle(\"Pearson's r= " << pearsonmp[i] << "\") + labs(x=\"Access through plasticity\", y=\"Access through mutations\") + geom_smooth(method=lm , color=\"black\", fill=\"gray85\", se=TRUE) + theme_classic()\n";
    fs << "dev.off()\n";
  }
  fs.close();
  bas.run_command("Rscript Results/"+cual+"/DetailsCounts/paR.sh");
  
  int g;
  double *jackal, *olapal, *pearal, *conm, *conp;
  jackal = new double[cuantas];
  olapal = new double[cuantas];
  pearal = new double[cuantas];
  for (i=0; i < cuantas; i++) {
    g = jacta.randint(0, cuantas);
    if (fenacmut[i] < fenacplas[g])
      min = fenacmut[i];
    else
      min = fenacplas[g];
    inter = 0;
    for (j=0; j < fenacmut[i]; j++) {
      strumu = sope[i].struct_in_mrep(j);
      for (k=0; k < fenacplas[g]; k++) {
        if (strumu==sope[g].one_compatible_structure(k)) {
          inter++;
          break;
        }
      }
    }
    jackal[i] = (inter*1.0)/ (fenacmut[i] + fenacplas[g] - (inter*1.0));
    olapal[i] = (inter*1.0)/(min*1.0);
    onion = fenacmut[i] + fenacplas[g] -(inter*1.0);
    conm = new double[onion];
    conp = new double[onion];
    l=0;
    for (j=0; j < fenacmut[i]; j++) {
      strumu = sope[i].struct_in_mrep(j);
      conm[l] = sope[i].count_in_mrep(j);
      conp[l] = 0;
      for (k=0; k< fenacplas[g];k++) {
        if (strumu==sope[g].one_compatible_structure(k)) {
          conp[l] = sope[g].p_Boltzmann(k);
          break;
        }
      }
      l++;
    }
    for (k=0; k< fenacplas[g];k++) {
      strumu = sope[g].one_compatible_structure(k);
      yapa = false;
      for (j=0; j < (fenacmut[i]); j++) {
        if (strumu==sope[i].struct_in_mrep(j)) {
          yapa = true;
          break;
        }
      }
      if (!yapa) {
        conm[l] = 0;
        conp[l] = sope[g].p_Boltzmann(k);
        l++;
      }
    }
//    pv = bas.pearsons_r(r, conm, conp, onion, 'g', S, df);
//    pv = pv*1.0;
    r = gsl_stats_correlation(conm, 1, conp, 1, onion);
    pearal[i] = r;
    delete [] conm;
    delete [] conp;
  }
  

  for (i=0; i < cuantas; i++) {
    sope[i].clear();
  }
  

  
  string coles[16] = {"Wrl205","Wrl2020","Wrl2080","Waf","pBoltzmann","Rmut","Pmut","Ppl","Jaccard","Olapcoef","Pearson", "MDfsubopt", "DDfsubopt", "JaccAl", "OlapAl", "PearAl"};
  string xlab[16] = {"Fitness in LR model 20:5", "Fitness in LR model 20:20","Fitness in LR model 20:80","Fitness in AF model", "Boltzmann prob. of stablest structure", "Mutational robustness", "Mutationally accesible phenotypes", "Plasticity accessible phenotypes", "Jaccard index", "Overlap coefficient", "Pearson coefficient", "Mean distance from optimum to alternative structures", "Std. dev. of the distance from optimum to alternative structures",  "Randomized Jaccard index", "Randomized overlap coefficient", "Randomized Pearson coefficient"};
  
  bas.open_ofstream(fs, "Results/"+cual+"/data.txt");
  fs << "Wrl205\tWrl2020\tWrl2080\tWaf\tpBoltzmann\tRmut\tPmut\tPpl\tJaccard\tOlapcoef\tPearson\tMDfsubopt\tDDfsubopt\tJaccAl\tOlapAl\tPearAl\n";
  for (i=0; i < cuantas; i++) {
    fs << wrl[0][i] << "\t" << wrl[1][i] << "\t" << wrl[2][i] << "\t" << waf[i] << "\t" << pb[i] << "\t" << romut[i] << "\t" << fenacmut[i] << "\t" << fenacplas[i] << "\t" << jaccardmp[i] << "\t" << overlapcomp[i] << "\t" << pearsonmp[i] << "\t" << mdfsuboptoopt[i] << "\t" << ddfsuboptoopt[i] << "\t" << jackal[i] << "\t" << olapal[i] << "\t" << pearal[i] << endl;
  }
  fs.close();
  
  double *nvec;
  nvec = new double[cuantas];
  
  double mean, sd, median, q1, q3;
  bas.open_ofstream(fs, "Results/"+cual+"/stats.txt");
  
  mean = bas.get_mean(wrl[0], cuantas);
  sd = bas.get_sample_stddev(wrl[0], cuantas, mean);
  bas.sort(wrl[0], nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[0] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  
  mean = bas.get_mean(wrl[1], cuantas);
  sd = bas.get_sample_stddev(wrl[1], cuantas, mean);
  bas.sort(wrl[1], nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[1] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  
  mean = bas.get_mean(wrl[2], cuantas);
  sd = bas.get_sample_stddev(wrl[2], cuantas, mean);
  bas.sort(wrl[2], nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[2] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(waf, cuantas);
  sd = bas.get_sample_stddev(waf, cuantas, mean);
  bas.sort(waf, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[3] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(pb, cuantas);
  sd = bas.get_sample_stddev(pb, cuantas, mean);
  bas.sort(pb, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[4] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(romut, cuantas);
  sd = bas.get_sample_stddev(romut, cuantas, mean);
  bas.sort(romut, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[5] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(fenacmut, cuantas);
  sd = bas.get_sample_stddev(fenacmut, cuantas, mean);
  bas.sort(fenacmut, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[6] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(fenacplas, cuantas);
  sd = bas.get_sample_stddev(fenacplas, cuantas, mean);
  bas.sort(fenacplas, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[7] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(jaccardmp, cuantas);
  sd = bas.get_sample_stddev(jaccardmp, cuantas, mean);
  bas.sort(jaccardmp, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[8] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(overlapcomp, cuantas);
  sd = bas.get_sample_stddev(overlapcomp, cuantas, mean);
  bas.sort(overlapcomp, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[9] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(pearsonmp, cuantas);
  sd = bas.get_sample_stddev(pearsonmp, cuantas, mean);
  bas.sort(pearsonmp, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[10] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
//  mean = bas.get_mean(spearmanmp, cuantas);
//  sd = bas.get_sample_stddev(spearmanmp, cuantas, mean);
//  bas.sort(spearmanmp, nvec, cuantas);
//  median = bas.get_midpoint(nvec, cuantas);
//  q1 = bas.get_q1(nvec, cuantas);
//  q3 = bas.get_q3(nvec, cuantas);
//  fs << "*" << xlab[9] << endl;
//  fs << "Mean: " << mean << endl;
//  fs << "Std. Dev: " << sd << endl;
//  fs << "Coeff. of variation: " << sd/mean << endl;
//  fs << "Median: " << median << endl;
//  fs << "Q1: " << q1 << endl;
//  fs << "Q3: " << q3 << endl;
//  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
//  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  //mdfsuboptoopt
  ///////////////////////////
  mean = bas.get_mean(mdfsuboptoopt, cuantas);
  sd = bas.get_sample_stddev(mdfsuboptoopt, cuantas, mean);
  bas.sort(mdfsuboptoopt, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[11] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  fs.close();
  ///////////////////////////
  mean = bas.get_mean(ddfsuboptoopt, cuantas);
  sd = bas.get_sample_stddev(ddfsuboptoopt, cuantas, mean);
  bas.sort(ddfsuboptoopt, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[12] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  fs.close();
  mean = bas.get_mean(jackal, cuantas);
  sd = bas.get_sample_stddev(jackal, cuantas, mean);
  bas.sort(jackal, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[13] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(olapal, cuantas);
  sd = bas.get_sample_stddev(olapal, cuantas, mean);
  bas.sort(olapal, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[14] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
  ///////////////////////////
  mean = bas.get_mean(pearal, cuantas);
  sd = bas.get_sample_stddev(pearal, cuantas, mean);
  bas.sort(pearal, nvec, cuantas);
  median = bas.get_midpoint(nvec, cuantas);
  q1 = bas.get_q1(nvec, cuantas);
  q3 = bas.get_q3(nvec, cuantas);
  fs << "*" << xlab[15] << endl;
  fs << "Mean: " << mean << endl;
  fs << "Std. Dev: " << sd << endl;
  fs << "Coeff. of variation: " << sd/mean << endl;
  fs << "Median: " << median << endl;
  fs << "Q1: " << q1 << endl;
  fs << "Q3: " << q3 << endl;
  fs << "Quartile coefficient of dispersion (q3-q1)/(q3+q1): " << (q3-q1)/(q3+q1) << endl;
  fs << "--------------------------------------------------------------\n";
 
  
  
  bas.open_ofstream(fs, "Results/"+cual+"/paR.sh");
  fs << "library(ggplot2)\n";
  fs << "library(Rcpp)\n";
  fs << "library(ggcorrplot)\n";
  fs << "a <- read.csv(file = \"Results/" << cual << "/data.txt\", head=TRUE, sep=\"\\t\")\n";
  
  fs << "b = a[, c(1,2,3,4,5,6,7,8,9,10,11,12,13)]\n";
  fs << "corr <- cor(b)\n";
  fs << "p.mat <- cor_pmat(b)\n";
  fs << "pdf(file=\"Results/" << cual << "/CorMatAll.pdf\")\n";
  fs << "ggcorrplot(corr, method=\"circle\", hc.order=FALSE, type=\"upper\", colors = c(\"steelblue\", \"white\", \"red\"), p.mat = p.mat)\n";
  fs << "dev.off()\n";
  
  for (i=0; i<13; i++) {
    fs << "pdf(file=\"Results/" << cual << "/Histo" << "_" << coles[i] << ".pdf\")\n";
    fs << "hist(a$" << coles[i] << ", prob=FALSE, main=\"\", xlab =\"" <<xlab[i] << "\", ylab=\"Frequency\", col=\"gray\")\n";
    fs << "dev.off()\n";
  }
  
  //lonuevo
  
  
  
  //para todas las correlaciones
  for (i =0; i < 12; i++) {
    for (j=i+1; j < 13; j++) {
      fs << "pru <- cor.test(a$" << coles[i] << ", a$" << coles[j] << ", alternative =\"two.sided\", method=\"pearson\", exact = TRUE)\n";
      fs << "capture.output(\"****" << coles[i] << " vs " << coles[j] << "\", file=\"Results/" << cual << "/correls.txt\", append = TRUE)\n";
      fs << "capture.output(pru$estimate, file=\"Results/" << cual << "/correls.txt\", append = TRUE)\n";
      fs << "capture.output(pru$p.value, file=\"Results/" << cual << "/correls.txt\", append = TRUE)\n";
      fs << "capture.output(pru$statistic, file=\"Results/" << cual << "/correls.txt\", append = TRUE)\n";
      fs << "capture.output(pru$parameter, file=\"Results/" << cual << "/correls.txt\", append = TRUE)\n";
      fs << "capture.output(\"--------------------------------\", file=\"Results/" << cual << "/correls.txt\", append = TRUE)\n";
      fs << "pdf(file=\"Results/"+cual+"/Scatcl_" << coles[i] << "Vs" << coles[j] << ".pdf\")\n";
      fs << "ggplot(a, aes(x=" << coles[i] << ", y=" << coles[j] << "))+ geom_point(size=6, shape=21, alpha=0.05, color=\"black\", fill=\"gray\", stroke=1) + labs(x=\"" << xlab[i] << "\", y=\"" << xlab[j] << "\") + geom_smooth(method=lm , color=\"black\", fill=\"gray80\", se=TRUE) + theme_classic()\n";
      fs << "dev.off()\n";
      fs << "pdf(file=\"Results/"+cual+"/Scat_" << coles[i] << "Vs" << coles[j] << ".pdf\")\n";
      fs << "ggplot(a, aes(x=" << coles[i] << ", y=" << coles[j] << "))+ geom_point(size=6, shape=21, alpha=0.05, color=\"black\", fill=\"gray\", stroke=1) + labs(x=\"" << xlab[i] << "\", y=\"" << xlab[j] << "\") + theme_classic()\n";
      fs << "dev.off()\n";
    }
    
    fs << "capture.output(\"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\", file=\"Results/" << cual << "/correls.txt\", append = TRUE)\n";
  }
  
  for (i=0; i < 4; i++) {
    fs << "pdf(file=\"Results/"+cual+"/Ill_pB_" << coles[i] << ".pdf\")\n";
    fs << "ggplot(a, aes(x= " << coles[4] << ", y=" << coles[i] << "))+ geom_point(size=6, shape=21, alpha=0.08, color=\"black\", fill=\"gray\", stroke=1) + geom_hline(yintercept = " << 1 << ", linetype=\"dotted\") + scale_x_continuous(expand=c(0,0), limits=c(0,1.02)) + scale_y_continuous(expand=c(0, 0), limits=c(0, " << 1.05 <<"), breaks = c(0, "<< 0.5 << ", "<< 1 << "), labels=c(\"0\", expression(frac(max,2)), \"max\")) + labs(x=\"" << xlab[4] << "\", y=\"" << xlab[i] << "\") + theme_classic()\n";
    fs << "dev.off()\n";
  }

  for (i=8; i<11; i++) {
    j = i+5;
    fs << "pdf(file=\"Results/"+cual+"/CongRandom_" << coles[i] << "Vs" << coles[j] << ".pdf\")\n";
    fs << "ggplot(a, aes(x=" << coles[i] << ", y=" << coles[j] << "))+ geom_point(size=6, shape=21, alpha=0.05, color=\"black\", fill=\"gray\", stroke=1) + labs(x=\"" << xlab[i] << "\", y=\"" << xlab[j] << "\") +  theme_classic() + geom_abline(color=\"black\")\n";
    fs << "dev.off()\n";
    fs << "wtest <- wilcox.test(a$" << coles[i] << ", a$" << coles[j] << ", paired=TRUE, alternative=\"greater\")\n";
    fs << "letr <- \"Wilcoxon test for " << xlab[i] << " vs " << xlab[j] <<"\"\n";
    fs << "capture.output(letr, file=\"Results/"+cual+"/CongRandom.txt\", append=TRUE)\n";
    fs << "capture.output(wtest, file=\"Results/"+cual+"/CongRandom.txt\", append=TRUE)\n";
    fs << "capture.output(\"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\", file=\"Results/"+cual+"/CongRandom.txt\", append=TRUE)\n";
  }
  
  for (i=0; i<11; i++) {
    for (j=0; j<11; j++) {
      if (i!=j) {
        fs << "limo <- lm(a$" << coles[i] << " ~ a$" << coles[j] << ")\n";
        fs << "capture.output(\"Pendiente de" << coles[i] << " como f de " << coles[j] << "\", file=\"Results/"+cual+"/SumRegLin.txt\", append=TRUE)\n";
        fs << "capture.output(summary(limo), file=\"Results/"+cual+"/SumRegLin.txt\", append=TRUE)\n";
        fs << "capture.output(\"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\", file=\"Results/"+cual+"/SumRegLin.txt\", append=TRUE)\n";
        
        fs << "capture.output(\"Pendiente de" << coles[i] << " como f de " << coles[j] << "\", file=\"Results/"+cual+"/Pendientes.txt\", append=TRUE)\n";
        fs << "capture.output(coefficients(limo)[2], file=\"Results/"+cual+"/Pendientes.txt\", append=TRUE)\n";
        fs << "capture.output(\"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\", file=\"Results/"+cual+"/Pendientes.txt\", append=TRUE)\n";
      }
    }
  }
  fs.close();
  bas.run_command("Rscript Results/"+cual+"/paR.sh");
  
  string garb;

  double mare[9][9];
  bas.open_ifstream(fe, "Results/"+cual+"/Pendientes.txt");
  
  for (i=0; i<9; i++) {
    for (j=0; j<9; j++) {
      if (i == j)
        mare[i][j] = 1;
      else {
        getline(fe, garb);
        getline(fe, garb);
        fe >> mare[i][j];
        getline(fe, garb);
        getline(fe, garb);
      }
    }
  }
  fe.close();
  
  
  
  bas.open_ofstream(fs, "Results/"+cual+"/Tpend.tex");
  
  fs << "\\documentclass{standalone}\n";
  fs << "\\usepackage[table]{xcolor}\n";
  fs << "\\begin{document}\n";
  fs << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}\\hline\n";
  fs << "Dependent&\\multicolumn{8}{c|}{Independent variable}\\\\\\cline{2-9}\n";
  fs << "variable";
  for (i=0; i<9; i++)
    fs << "&" << coles[i];
  fs << "\\\\\\hline\n";
  for (i=0; i < 9; i++) {
    fs << coles[i];
    for (j=0; j < 9; j++)
      fs << "&" << mare[i][j];
    fs << "\\\\\\hline\n";
  }
  fs << "\\end{tabular}\n";
  fs << "\\end{document}\n";
  fs.close();
  bas.run_command("pdflatex -output-directory Results/"+cual+"/ Results/"+cual+"/Tpend.tex");

  for (i=0; i < 3; i++) {
    delete [] wrl[i];
  }
  delete [] wrl;
  delete [] waf;
  delete [] pb;
  delete [] romut;
  delete [] fenacmut;
  delete [] fenacplas;
  delete [] jaccardmp;
  delete [] overlapcomp;
  delete [] pearsonmp;
//  delete [] spearmanmp;
  
  for (i=0; i < cuantas; i++) {
    delete [] conmu[i];
    delete [] conpl[i];
  }
  delete [] conmu;
  delete [] conpl;
  delete [] mdfsuboptoopt;
  delete [] ddfsuboptoopt;
  return 0;
}
