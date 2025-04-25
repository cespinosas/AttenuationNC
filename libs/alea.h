#ifndef ALEA_H
#define ALEA_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;

class Alea
{
public:
  Alea(int sem);
  Alea();
  void start_rng(int sem);
  double randreal();
  int randint(int from, int to);// [from, to)
  void close_rng();
  void check_rng(ostream& sal);
  bool toss();
  double randgauss();
  double randgauss(double mean, double stdev);
  const gsl_rng_type * typeofrng;
  gsl_rng * rng;
  int ran_binomial(double p, int n);
private:
  
};

#endif



