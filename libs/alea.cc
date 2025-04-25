#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include "alea.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;

Alea::Alea(int sem) //creates instance of class Alea and sets random seed to sem
{
  gsl_rng_env_setup();
  typeofrng = gsl_rng_ranlxs2;
  rng = gsl_rng_alloc (typeofrng);
  gsl_rng_set(rng, sem);
}

Alea::Alea()
{
}

void Alea::start_rng(int sem) //set up and sets random seed to sem
{
  gsl_rng_env_setup();
  typeofrng = gsl_rng_ranlxs2;
  rng = gsl_rng_alloc (typeofrng);
  gsl_rng_set(rng, sem);
  
}

void Alea::close_rng() //frees rng
{
  gsl_rng_free (rng);
}

double Alea::randreal() //returns random double in the [0, 1) interval
{
  double res = gsl_rng_uniform (rng);
  return res;
}

int Alea::randint(int from, int to) {   //returns random int in the [from, to) interval
  int res;
  int range = to-from;
  res = from + gsl_rng_uniform_int(rng, range);
  return res;
}

void Alea::check_rng(ostream& sal) {
  sal << gsl_rng_name (rng) << endl;
}

bool Alea::toss() { //toss a coin
  bool res = true;
  double nu = randreal();
  if (nu < 0.5)
    res = false;
  return res;
}

double Alea::randgauss() {
  double x,y,r,phi;
  x=randreal();
  y=randreal();
  if(x==0)
    x=randreal();
  r=sqrt(-2.*log(x));
  phi=2.*3.141592654*y;
  return r*cos(phi);
}

double Alea::randgauss(double mean, double stdev) {
  double res = randgauss();
  res = res*stdev;
  res = res+mean;
  return res;
}

int Alea::ran_binomial(double p, int n) {
  int k = gsl_ran_binomial(rng, p, n);
  return k;
}
