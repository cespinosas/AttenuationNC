#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "alea.h"
#include "basics.h"

using namespace std;


struct lospar {
  double *ka; //0.001
  double *kd;
  int numec;
  double *B;
  double *dis;
  double Rtot;
  double Ltot;
  double kda;
  double kdb;
};

void print_state (size_t iter, gsl_multiroot_fdfsolver *s, int nue);
int lrmod(const gsl_vector *x, void *params, gsl_vector *f);
int lrmod_df(const gsl_vector *x, void *params, gsl_matrix *J);
int lrmod_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);

int main(int argc, char **argv) {
  Alea jacta(1);
  Basics bas(jacta);
  
  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  int i,j;
  int status=0;
  size_t iter = 0;
  size_t n = 2;
  
  struct lospar p;
  p.numec = n;
  p.ka = new double[p.numec];
  p.kd = new double[p.numec];
  p.B = new double[p.numec];
  p.dis = new double[p.numec];
//  cout << "R_tot: ";
//  cin >> p.Rtot;
//  cout << "L_tot: ";
//  cin >> p.Ltot;
//  p.Rtot = 10;
//  p.Ltot = 10;
  p.kda = 0.001;
  p.kdb = 0.5;
  
  p.Rtot = 20;
  p.dis[0] = 0;
  p.dis[1] = atoi(argv[1]);
  for (i=0; i < p.numec; i++) {
    p.ka[i] = 0.001;
    p.kd[i] = p.kda*(pow(2,(p.kdb*p.dis[i])));
  }

  double x_init[2];
  gsl_vector *x = gsl_vector_alloc (n);

  //quitar
//  p.Ltot = 5;
//  p.B[0] = 0.4;
//  p.B[1] = 0.6;
  //
  
  double eles[3] = {5, 20, 80};
  
  string losti[3] = {"20:5 RNA:ligand ratio", "20:20 RNA:ligand ratio", "20:80 RNA:ligand ratio"};
  
  double maxi[3] = {4.69388, 16, 19.67388};
  
  double res[3][51][2];
  ofstream fs;
  
  for (i=0; i<3; i++) {
    p.Ltot = eles[i];
    for (j=0; j<=50; j++) {
      p.B[0] = j*0.02;
      p.B[1] = 1 - p.B[0];
      gsl_multiroot_function_fdf f = {&lrmod, &lrmod_df, &lrmod_fdf, n, &p};
      x_init[0] = 2;
      x_init[1] = 2;
      gsl_vector_set (x, 0, x_init[0]);
      gsl_vector_set (x, 1, x_init[1]);
      T = gsl_multiroot_fdfsolver_hybridsj;
      s = gsl_multiroot_fdfsolver_alloc (T, n);
      gsl_multiroot_fdfsolver_set (s, &f, x);
      do
      {
        iter++;
        status = gsl_multiroot_fdfsolver_iterate (s);
//        print_state (iter, s, p.numec);
        if (status)
          break;
        status = gsl_multiroot_test_residual (s->f, 1e-14);
      } while (status == GSL_CONTINUE && iter < 5000);
      //printf ("status = %s\n", gsl_strerror (status));
      
      res[i][j][0] = gsl_vector_get(s->x, 0);
      res[i][j][1] = gsl_vector_get(s->x, 1);
//      res[i][j][2] = res[i][j][0]+res[i][j][1];
      
      
      gsl_multiroot_fdfsolver_free (s);
      
      
    }
  }
  gsl_vector_free (x);
  for (i=0; i < 3; i++) {
    res[i][0][0] = 0;
    bas.open_ofstream(fs, "d"+bas.inttostring(p.dis[1])+"_L"+bas.inttostring(int(eles[i]))+".txt");
    fs << "pB\tComp\tCuant\n";//\tC1\tCT\n";
    for (j=0; j <=50; j++) {
      fs << j*0.02 << "\tC0\t" << res[i][j][0] << endl;
      fs << j*0.02 << "\tC1\t" << res[i][j][1] << endl;
      fs << j*0.02 << "\tCT\t" << res[i][j][0]+res[i][j][1] << endl;
    }
    fs.close();
  }
  
  bas.open_ofstream(fs, "paR.sh");
  fs << "library(ggplot2)\n";
  fs << "library(Rcpp)\n";
  fs << "library(ggcorrplot)\n";
  i=0;
  fs << "ad <- read.csv(file = \"" << "d"+bas.inttostring(p.dis[1])+"_L"+bas.inttostring(int(eles[i]))+".txt" << "\", head=TRUE, sep=\"\\t\")\n";
    fs << "pdf(file=\"" << "d"+bas.inttostring(p.dis[1])+"_L"+bas.inttostring(int(eles[i]))+".pdf\")\n";
    fs << "ggplot(ad, aes(x=pB, y=Cuant, group=Comp)) + geom_line(aes(color=Comp), size = 3) + labs(title=\"" << losti[i] << "\", x=\"Boltzmann probability of \\noptimum structure\", y=\"\") + scale_color_discrete(name = \"Complex\", labels = c(\"Optimum:Ligand\", \"Sub-optimum:Ligand\", \"Total\"))+ geom_hline(size = 3, yintercept = " << maxi[i] << ", linetype=\"dashed\") + scale_x_continuous(expand=c(0,0), limits=c(0,1.02), breaks=c(0, 0.5, 1), labels =c(0, 0.5,1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, " << maxi[i]*1.05 <<"), breaks = c(0, "<< maxi[i]/2 << ", "<< maxi[i] << "), labels=c(\"\", \"\", \"\")) + theme_classic() + theme(legend.position = c(0.6, 0.3), text=element_text(size=38), axis.line=element_line(size=2), plot.title = element_text(size=36))\n";
    fs << "dev.off()\n";
  
  i=1;
  fs << "ad <- read.csv(file = \"" << "d"+bas.inttostring(p.dis[1])+"_L"+bas.inttostring(int(eles[i]))+".txt" << "\", head=TRUE, sep=\"\\t\")\n";
    fs << "pdf(file=\"" << "d"+bas.inttostring(p.dis[1])+"_L"+bas.inttostring(int(eles[i]))+".pdf\")\n";
    fs << "ggplot(ad, aes(x=pB, y=Cuant, group=Comp)) + geom_line(aes(color=Comp), size = 3, show.legend = FALSE) + labs(title=\"" << losti[i] << "\", x=\"Boltzmann probability of \\noptimum structure\", y=\"\") + scale_color_discrete(name = \"Complex\", labels = c(\"Optimum:Ligand\", \"Sub-optimum:Ligand\", \"Total\"))+ geom_hline(size = 3, yintercept = " << maxi[i] << ", linetype=\"dashed\") + scale_x_continuous(expand=c(0,0), limits=c(0,1.02), breaks=c(0, 0.5, 1), labels =c(0, 0.5,1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, " << maxi[i]*1.05 <<"), breaks = c(0, "<< maxi[i]/2 << ", "<< maxi[i] << "), labels=c(\"\", \"\", \"\")) + theme_classic() + theme(legend.position = c(0.6, 0.3), text=element_text(size=38), axis.line=element_line(size=2), plot.title = element_text(size=36))\n";
    fs << "dev.off()\n";
  
  i=2;
  fs << "ad <- read.csv(file = \"" << "d"+bas.inttostring(p.dis[1])+"_L"+bas.inttostring(int(eles[i]))+".txt" << "\", head=TRUE, sep=\"\\t\")\n";
    fs << "pdf(file=\"" << "d"+bas.inttostring(p.dis[1])+"_L"+bas.inttostring(int(eles[i]))+".pdf\")\n";
    fs << "ggplot(ad, aes(x=pB, y=Cuant, group=Comp)) + geom_line(aes(color=Comp), size = 3, show.legend = FALSE) + labs(title=\"" << losti[i] << "\", x=\"Boltzmann probability of \\noptimum structure\", y=\"Concentration\") + scale_color_discrete(name = \"Complex\", labels = c(\"Optimum:Ligand\", \"Sub-optimum:Ligand\", \"Total\"))+ geom_hline(size = 3, yintercept = " << maxi[i] << ", linetype=\"dashed\") + scale_x_continuous(expand=c(0,0), limits=c(0,1.02), breaks=c(0, 0.5, 1), labels =c(0, 0.5,1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, " << maxi[i]*1.05 <<"), breaks = c(0, "<< maxi[i]/2 << ", "<< maxi[i] << "), labels=c(\"0\", expression(frac(max,2)), \"max\")) + theme_classic() + theme(legend.position = c(0.6, 0.3), text=element_text(size=38), axis.line=element_line(size=2), plot.title = element_text(size=36))\n";
    fs << "dev.off()\n";
  
  double wsub = 1.0/(0.01 + (p.dis[1]/(2*72.0)));
  fs << "x <- c(0,1,0,1,0,1)\n";
  fs << "cual <- c(\"opt\", \"opt\", \"sub\", \"sub\", \"tot\", \"tot\")\n";
  fs << "y <- c(0, 100, " << wsub << ", 0, " << wsub << ", 100)\n";
  fs << "af <- data.frame(x,cual,y)\n";
  fs << "pdf(file=\"" << "d"+bas.inttostring(p.dis[1])+"_AF.pdf\")\n";
  fs << "ggplot(af, aes(x=x, y=y, group=cual)) + geom_line(aes(color=cual), size = 3) + labs(title=\"Ancel-Fontana model\", x=\"Boltzmann probability of \\n optimum structure\", y=\"Fitness\") + scale_color_discrete(name =\"Contributions\", labels=c(\"Optimum\", \"Sub-optimum\", \"Total\"))+ geom_hline(size = 3, yintercept = 100, linetype=\"dashed\") + scale_x_continuous(expand=c(0,0), limits=c(0,1.02), breaks=c(0, 0.5, 1), labels =c(0, 0.5,1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, 105), breaks=c(0, 50,100), labels=c(\"0\", expression(frac(max,2)), \"max\")) + theme_classic() + theme(legend.position = c(0.6, 0.3), text=element_text(size=38), axis.line=element_line(size=2), plot.title = element_text(size=36))\n"; // legend.text=element_text(size=12), legend.title=element_text(size=14))\n";
  fs << "dev.off()\n";

  fs.close();
  
  bas.run_command("Rscript paR.sh");
  
  delete [] p.ka;
  delete [] p.kd;
  delete [] p.B;
  delete [] p.dis;
  
  return 0;
}

int lrmod(const gsl_vector *x, void *params, gsl_vector *f) {
  double *ka, *kd, *B, *dis, Rtot, Ltot;
  int numec = ((struct lospar *) params)->numec;
  Rtot = ((struct lospar *) params)->Rtot;
  Ltot = ((struct lospar *) params)->Ltot;
  ka = new double[numec];
  kd = new double[numec];
  B = new double[numec];
  dis = new double[numec];
  int i;
  for (i=0; i < numec; i++) {
    ka[i]  = ((struct lospar *) params)->ka[i];
    kd[i]  = ((struct lospar *) params)->kd[i];
    B[i]  = ((struct lospar *) params)->B[i];
//    B[i] *= Rtot;
    dis[i]  = ((struct lospar *) params)->dis[i];
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
    y[i] = (ka[i]*B[i]*(Rtot - sutoc)*(Ltot - sutoc)) - (kd[i]*comp[i]);
    
  for (i=0; i< numec; i++)
    gsl_vector_set(f, i, y[i]);
  
  delete [] y;
  delete [] comp;
  delete [] ka;
  delete [] kd;
  delete [] B;
  delete [] dis;
  return GSL_SUCCESS;
}

int lrmod_df(const gsl_vector *x, void *params, gsl_matrix *J) {
  double *ka, *kd, *B, *dis, Rtot, Ltot;
  int numec = ((struct lospar *) params)->numec;
  Rtot = ((struct lospar *) params)->Rtot;
  Ltot = ((struct lospar *) params)->Ltot;
  ka = new double[numec];
  kd = new double[numec];
  B = new double[numec];
  dis = new double[numec];
  int i,j;
  for (i=0; i < numec; i++) {
    ka[i]  = ((struct lospar *) params)->ka[i];
    kd[i]  = ((struct lospar *) params)->kd[i];
    B[i]  = ((struct lospar *) params)->B[i];
//    B[i] *= Rtot;
    dis[i]  = ((struct lospar *) params)->dis[i];
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
        df[i][j] = (ka[i]*B[i]*((2*sutoc) - Rtot - Ltot)) -kd[i];
      } else {
        df[i][j] = ka[i]*B[i]*((2*sutoc) - Rtot - Ltot);
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
  delete [] ka;
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

void print_state (size_t iter, gsl_multiroot_fdfsolver *s, int nue) {
  cout << "iter = " << iter << "  x = [";
  int i;
  for (i=0; i < nue; i++) {
    cout << gsl_vector_get(s->x, i);
    if (i < nue-1)
      cout << ", ";
  }
  cout << "], f(x) = [";
  for (i=0; i < nue; i++) {
    cout << gsl_vector_get(s->f, i);
    if (i < nue-1)
      cout << ", ";
  }
  cout << "]\n";
}
