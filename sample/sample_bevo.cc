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


using namespace std;



int main(int argc, char **argv) {
  
  if (argc != 3) {
    cout << "[Error]: Wrong number of arguments.\n";
    exit(1);
  }
  string cual(argv[2]);
  Alea jacta(atoi(argv[1]));
  Basics bas(jacta);
  ifstream fe;
  bas.open_ifstream(fe, "../Structures/"+cual+".txt");
  string tau;
  fe >> tau;
  fe.close();
  
  GaRNAcha sope(jacta, tau.length());
  
  ofstream fs;
  bas.open_ofstream_to_append(fs, "SeqForEvo/"+cual+".tex");
//  sope.sample(fs, 1000, 60, tau);
  int i;
  for (i=0; i < 500; i++) {
    sope.force_inverse_fold(tau);
    fs << sope.nt_sequence() << endl;
    sope.clear_seq();
  }
  fs.close();
  
  sope.clear();

  return 0;
}
