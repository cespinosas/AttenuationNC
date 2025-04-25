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
  
  if (argc != 2) {
    cout << "[Error]: Wrong number of arguments.\n";
    exit(1);
  }
  Alea jacta(atoi(argv[1]));
  Basics bas(jacta);
  ofstream fs;
  int i, lalong=60;
  GaRNAcha sope(jacta);
  
  for (i=0; i < 10; i++) {
    bas.open_ofstream(fs, "../Structures/ran"+bas.inttostring(i)+".txt");
    sope.set_size(lalong);
    sope.rand_seq();
    sope.get_mfe_str();
    sope.print_mfestr(fs);
    fs << endl;
    sope.clear();
    fs.close();
  }
  

  return 0;
}
