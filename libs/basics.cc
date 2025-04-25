#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <set>
#include <list>
#include <string>
#include "alea.h"
#include "basics.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <sys/stat.h>
#include <dirent.h>

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ostream;
using std::string;
using std::set;
using std::list;
using std::ifstream;
using std::ios;

Basics::Basics()
{
}

Basics::Basics(Alea& jacta) //creates instance of class Basics and assigns rng to jacta
{
	start_rng(jacta);
}

void Basics::start_rng(Alea& jacta) //assigns rng to jacta
{
  est = jacta;
}

//vectors
bool Basics::all_entries_are_finite (double *vec, int tam) {
  bool res = true;
  int i;
  for (i = 0; i < tam; i++)
    if (!std::isfinite(vec[i])) {
      res = false;
      break;
    }
  return res;
}

int Basics::number_of_finite_entries (double *vec, int tam) {
  int res = 0, i;
  for (i = 0; i < tam; i++)
    if (std::isfinite(vec[i]))
      res++;
  return res;
}

int Basics::where_is_vector_norm (int *vec, int tam) {
  int i, val, p=0;
  val = 0;
  for (i = 0; i < tam; i++)
    if (val < abs(vec[i])) {
      p = i;
      val = abs(vec[i]);
    }
  return p;
}

int Basics::where_is_vector_norm (double *vec, int tam) {
  int i, p=0;
  double val = 0;
  for (i = 0; i < tam; i++)
    if (val < fabs(vec[i])) {
      p = i;
      val = fabs(vec[i]);
    }
  return p;
}

bool Basics::eqvec(int vec1[], int s1, int vec2[], int s2) //are two vectors equal
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i]) {
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(double vec1[], int s1, double vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i]) {
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(bool vec1[], int s1, bool vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])	{
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(char vec1[], int s1, char vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])	{
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::eqvec(string vec1[], int s1, string vec2[], int s2)
{
  int i;
  bool res = true;
  if (s1 != s2)
    res = false;
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])	{
				res = false;
				break;
			}
  }
  return res;
}

bool Basics::same_elements_in_vecs(int vec1[], int s1, int vec2[], int s2)
{
	bool ans = true;
	int i;
	if (s1!=s2)
		ans = false;
	else {
		for (i=0; i<s1; i++)
			if (!vector_contains(vec2, vec1[i], s2)) {
				ans = false;
				break;
			}
	}
	return ans;
}

bool Basics::same_elements_in_vecs(double vec1[], int s1, double vec2[], int s2)
{
	bool ans = true;
	int i;
	if (s1!=s2)
		ans = false;
	else {
		for (i=0; i<s1; i++)
			if (!vector_contains(vec2, vec1[i], s2)) {
				ans = false;
				break;
			}
	}
	return ans;
}

bool Basics::same_elements_in_vecs(char vec1[], int s1, char vec2[], int s2)
{
	bool ans = true;
	int i;
	if (s1!=s2)
		ans = false;
	else {
		for (i=0; i<s1; i++)
			if (!vector_contains(vec2, vec1[i], s2)) {
				ans = false;
				break;
			}
	}
	return ans;
}

bool Basics::same_elements_in_vecs(string vec1[], int s1, string vec2[], int s2)
{
	bool ans = true;
	int i;
	if (s1!=s2)
		ans = false;
	else {
		for (i=0; i<s1; i++)
			if (!vector_contains(vec2, vec1[i], s2)) {
				ans = false;
				break;
			}
	}
	return ans;
}

int Basics::difsinvecs(int vec1[], int s1, int vec2[], int s2) //number of differences between two vectors
{
  int i;
  int dif = 0;
  if (s1 != s2)
  {
    cout << "[Error]: You can not count differences between vectors with different sizes using Basics::difsinvecs.\n";
    exit(1);
  }
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])
				dif++;
  }
  return dif;
}

int Basics::difsinvecs(double vec1[], int s1, double vec2[], int s2) //number of differences between two vectors
{
  int i;
  int dif = 0;
  if (s1 != s2)
  {
    cout << "[Error]: You can not count differences between vectors with different sizes using Basics::difsinvecs.\n";
    exit(1);
  }
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])
				dif++;
  }
  return dif;
}

int Basics::difsinvecs(char vec1[], int s1, char vec2[], int s2)
{
  int i;
  int dif = 0;
  if (s1 != s2)
  {
		cout << "[Error]: You can not count differences between vectors with different sizes using Basics::difsinvecs.\n";
		exit(1);
  }
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])
				dif++;
  }
  return dif;
}

int Basics::difsinvecs(bool vec1[], int s1, bool vec2[], int s2)
{
  int i;
  int dif = 0;
  if (s1 != s2)
  {
    cout << "[Error]: You can not count differences between vectors with different sizes using Basics::difsinvecs.\n";
    exit(1);
  }
  else {
    for (i = 0; i < s1; i++)
      if (vec1[i] != vec2[i])
				dif++;
  }
  return dif;
}

double Basics::mean_dist_vecs(int *vec1, int s1, int *vec2, int s2)
{
  if (s1 != s2) {
    cout << "[Error]: Vectors are not equal-sized in Basics::mean_dist_vecs.\n";
    exit(1);
  }
  double res = 0;
  int i;
  for (i=0; i < s1; i++)
    res = res + abs(vec1[i] - vec2[i]);
  res = res/double(2*s1);
  return res;
}

double Basics::mean_dist_vecs(double *vec1, int s1, double *vec2, int s2)
{
  if (s1 != s2) {
    cout << "[Error]: Vectors are not equal-sized in Basics::mean_dist_vecs.\n";
    exit(1);
  }
  double res = 0;
  int i;
  for (i=0; i < s1; i++)
    res = res + fabs(vec1[i] - vec2[i]);
  res = res/double(2*s1);
  return res;
}

double Basics::multinvec(int vec[], int s)
{
  int i;
  double res = 1.0;
  for (i = 0; i < s; i++)
    res = res*vec[i];
  return res;
}

double Basics::multinvec(double vec[], int s)
{
  int i;
  double res = 1.0;
  for (i = 0; i < s; i++)
    res = res*vec[i];
  return res;
}

void Basics::fillv0(int vec[], int s) //fills 1d array with 0s
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 0;
}

void Basics::fillv0(double vec[], int s) //fills 1d array with 0s
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 0;
}

void Basics::fillv0(bool vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = false;
}

void Basics::fillv0(string vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = "";  
}

void Basics::fillv1(int vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 1;
}

void Basics::fillv1(double vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = 1;
}

void Basics::fillv1(bool vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = true;
}
//

void Basics::fillvm1(int vec[], int s) //fills vector with -1
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = -1;
}

void Basics::fillvm1(double vec[], int s) //fills vector with -1
{
  int i;
  for (i = 0; i < s; i++)
    vec[i] = -1;
}

int Basics::sumatoria(int vec[], int s)
{
	int res=0;
	int i;
	for (i=0; i<s; i++)
		res = res + vec[i];
	return res;
}

double Basics::sumatoria(double vec[], int s)
{
	double res=0;
	int i;
	for (i=0; i<s; i++)
		res = res + vec[i];
	return res;
}

void Basics::sort_vector(int vecor[], int orden[], int len, int lbo)
{
	int i, min = lbo,mie;
	int *vcop;
	vcop = new int[len];
	for (i=0; i<len; i++)
		vcop[i] = vecor[i];
	for (i=0; i < len; i++) {
		mie = find_minlb(vcop, len, min);
		orden[i] = mie;
		vcop[mie] = lbo-1;
	}
	delete [] vcop;
	return;
}

void Basics::sort_vector(double vecor[], int orden[], int len, double lbo)
{
	int i, mie;
	double *vcop, min = lbo;
	vcop = new double[len];
	for (i=0; i<len; i++)
		vcop[i] = vecor[i];
	for (i=0; i < len; i++) {
		mie = find_minlb(vcop, len, min);
		orden[i] = mie;
		vcop[mie] = floor(lbo)-1;
	}
	delete [] vcop;
	return;
}

int Basics::first_in_vector_with(int pat, int *vec, int siz) {
  return first_in_vector_with(pat, vec, siz, 0);
}

int Basics::first_in_vector_with(string pat, string *vec, int siz) {
  return first_in_vector_with(pat, vec, siz, 0);
}

int Basics::first_in_vector_with(char pat, char *vec, int siz) {
  return first_in_vector_with(pat, vec, siz, 0);
}

int Basics::first_in_vector_with(bool pat, bool *vec, int siz) {
  return first_in_vector_with(pat, vec, siz, 0);
}

int Basics::first_in_vector_with(int pat, int *vec, int siz, int desde) {
  int i, res = -1;
  for (i = desde; i < siz; i++)
    if (vec[i] == pat) {
      res = i;
      break;
    }
  return res;
}

int Basics::first_in_vector_with(string pat, string *vec, int siz, int desde) {
  int i, res = -1;
  for (i = desde; i < siz; i++)
    if (vec[i] == pat) {
      res = i;
      break;
    }
  return res;
}

int Basics::first_in_vector_with(char pat, char *vec, int siz, int desde) {
  int i, res = -1;
  for (i = desde; i < siz; i++)
    if (vec[i] == pat) {
      res = i;
      break;
    }
  return res;
}

int Basics::first_in_vector_with(bool pat, bool *vec, int siz, int desde) {
  int i, res = -1;
  for (i = desde; i < siz; i++)
    if (vec[i] == pat) {
      res = i;
      break;
    }
  return res;
}

bool Basics::vector_contains(int *vect, int pat, int des, int siz)
{
	bool res = false;
	int i;
	for (i=des; i<siz; i++)
		if (vect[i]==pat) {
			res = true;
			break;
		}
	return res;
}

bool Basics::vector_contains(int *vect, int pat, int siz)
{
	bool res = vector_contains(vect, pat, 0, siz);
	return res;
}

bool Basics::vector_contains(double *vect, double pat, int des, int siz)
{
	bool res = false;
	int i;
	for (i=des; i<siz; i++)
		if (vect[i]==pat) {
			res = true;
			break;
		}
	return res;
}

bool Basics::vector_contains(double *vect, double pat, int siz)
{
	bool res = vector_contains(vect, pat, 0, siz);
	return res;
}

bool Basics::vector_contains(bool *vect, bool pat, int des, int siz)
{
	bool res = false;
	int i;
	for (i=des; i<siz; i++)
		if (vect[i]==pat) {
			res = true;
			break;
		}
	return res;
}

bool Basics::vector_contains(bool *vect, bool pat, int siz)
{
	bool res = vector_contains(vect, pat, 0, siz);
	return res;
}

bool Basics::vector_contains(string *vect, string pat, int des, int siz)
{
	bool res = false;
	int i;
	for (i=des; i<siz; i++)
		if (vect[i]==pat) {
			res = true;
			break;
		}
	return res;
}

bool Basics::vector_contains(string *vect, string pat, int siz)
{
	bool res = vector_contains(vect, pat, 0, siz);
	return res;
}


bool Basics::vector_contains(char *vect, char pat, int des, int siz)
{
	bool res = false;
	int i;
	for (i=des; i<siz; i++)
		if (vect[i]==pat) {
			res = true;
			break;
		}
	return res;
}

bool Basics::vector_contains(char *vect, char pat, int siz)
{
	bool res = vector_contains(vect, pat, 0, siz);
	return res;
}

bool Basics::last_vector(bool *ve, int siz)
{
	bool res = true;
	int i;
	for (i=0; i<siz; i++)
		if (!ve[i]){
			res = false;
			break;
		}
	return res;
}

bool Basics::last_vector(int *ve, int siz, int max)
{
	bool res = true;
	int i;
	for (i=0; i<siz; i++)
		if (ve[i] != max){
			res = false;
			break;
		}
	return res;
}

void Basics::next_in_vector(bool *ve, int siz)
{
	int i;
	if (last_vector(ve, siz)) {
		cout << "[Error]: Last vector already reached in Basics::next_in_vector.\n";
		exit(1);
	}
	for (i=0; i<siz; i++)
		if (!ve[i]) {
			ve[i] = true;
			break;
		}
	fillv0(ve, i);
}

void Basics::next_in_vector(int *ve, int siz, int min, int max, int step)
{
	int i,j;
	if (last_vector(ve, siz, max)) {
		cout << "[Error]: Last vector already reached in Basics::next_in_vector.\n";
		exit(1);
  }
  for (i=0; i<siz; i++)
		if (ve[i] < max) {
			ve[i] = ve[i]+step;
			break;
		}
  for (j=0; j<i; j++)
    ve[j]=min;
}

int Basics::bintoint(int *vect, int tam, int max) {
  int i;
  int suma = 0;
  for (i=0; i<tam; i++) {
    if (vect[i] == max)
      suma += pow(2.0, i);
  }
  return suma;
}

void Basics::inttobin(int ent, int tam, int max, int step, int *vect) {
  int elnum = ent;
  fillv0(vect, tam);
  int i;
  for (i=tam-1; i>=0; i--) {
    if (elnum >= pow(2,i)) {
      vect[i] = max;
      elnum -= pow(2,i);
    } else
      vect[i] = max-step;
  }
}

int Basics::count_in_vector(bool *ve, int siz, bool pat)
{
	int res = 0, i;
	for (i=0; i<siz; i++) {
		if (ve[i] == pat)
			res++;
	}
	return res;
}

int Basics::count_in_vector(int *ve, int siz, int pat)
{
	int res = 0, i;
	for (i=0; i<siz; i++) {
		if (ve[i] == pat)
			res++;
	}
	return res;
}

int Basics::count_in_vector(double *ve, int siz, double pat)
{
	int res = 0, i;
	for (i=0; i<siz; i++) {
		if (ve[i] == pat)
			res++;
	}
	return res;
}

void Basics::copy_vector(bool *fromthis, bool *tothis, int siz)
{
	int i;
	for (i=0; i<siz; i++)
		tothis[i] = fromthis[i];
}

void Basics::copy_vector(int *fromthis, int *tothis, int siz)
{
  int i;
  for (i=0; i<siz; i++)
    tothis[i] = fromthis[i];
}

//vector formatting
void Basics::printv(ostream& sal, int vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    sal << vec[i] << "\t";
}

void Basics::printv(ostream& sal, double vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    sal << vec[i] << "\t";
}

void Basics::printv(ostream& sal, bool vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    sal << vec[i] << "\t";
}

void Basics::printv(ostream& sal, string vec[], int s)
{
  int i;
  for (i = 0; i < s; i++)
    sal << vec[i] << "\t";
}

//vectors and matrices
int Basics::vecinmat(int** mat, int rows, int cols, int vec[], int s) //vector appears as matrix row
{
  int i;
  int res = -1;
  if (cols != s) {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols))	{
				res = i;
				break;
			}
  }
  return res;
}


int Basics::vecinmat(double** mat, int rows, int cols, double vec[], int s)
{
  int i;
  int res = -1;
  if (cols != s) {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols)){
				res = i;
				break;
			}
  }
  return res;
}

int Basics::vecinmat(bool** mat, int rows, int cols, bool vec[], int s)
{
  int i;
  int res = -1;
  if (cols != s)   {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols))	{
				res = i;
				break;
			}
  }
  return res;
}

int Basics::vecinmat(string** mat, int rows, int cols, string vec[], int s)
{
  int i;
  int res = -1;
  if (cols != s)   {
		cout << "[Error]: You can not search a vector in a matrix of different width using Basics::vecinmat.\n";
		exit(1);
	}
  else {
    for (i = 0; i < rows; i++)
      if (eqvec(vec, s, mat[i], cols))	{
				res = i;
				break;
			}
  }
  return res;
}

int Basics::dotproduct(int *vec1, int *vec2, int siz)
{
  int res = 0;
  int i;
  for (i = 0; i < siz; i++)
    res = res + (vec1[i]*vec2[i]);
  return res;
}

double Basics::dotproduct(int *vec1, double *vec2, int siz)
{
  double res = 0;
  int i;
  for (i = 0; i < siz; i++)
    res = res + (vec1[i]*vec2[i]);
  return res;
}

double Basics::dotproduct(double *vec1, int *vec2, int siz)
{
  double res = 0;
  int i;
  for (i = 0; i < siz; i++)
    res = res + (vec1[i]*vec2[i]);
  return res;
}

double Basics::dotproduct(double *vec1, double *vec2, int siz)
{
  double res = 0;
  int i;
  for (i = 0; i < siz; i++)
    res = res + (vec1[i]*vec2[i]);
  return res;
}

void Basics::matxvec(int *vec, int sizv, int leng, int wide, int **ma, int *res)
{
  if (sizv != wide) {
		cout << "[Error]: This vector x matrix product can not be calculated using Basics::matxvec.\n";
		exit(1);
	}
  fillv0(res, sizv);
  int i;
  for (i=0; i < leng; i++)
    res[i] = dotproduct(vec, ma[i], sizv);
  return;
}

void Basics::matxvec(int *vec, int sizv, int leng, int wide, double **ma, double *res)
{
  if (sizv != wide) {
		cout << "[Error]: This vector x matrix product can not be calculated using Basics::matxvec.\n";
		exit(1);
	}
  fillv0(res, sizv);
  int i;
  for (i=0; i < leng; i++)
    res[i] = dotproduct(vec, ma[i], sizv);
  return;
}

void Basics::matxvec(double *vec, int sizv, int leng, int wide, int **ma, double *res)
{
  if (sizv != wide) {
		cout << "[Error]: This vector x matrix product can not be calculated using Basics::matxvec.\n";
		exit(1);
  }
  fillv0(res, sizv);
  int i;
  for (i=0; i < leng; i++)
    res[i] = dotproduct(vec, ma[i], sizv);
  return;
}

void Basics::matxvec(double *vec, int sizv, int leng, int wide, double **ma, double *res)
{
  if (sizv != wide) {
		cout << "[Error]: This vector x matrix product can not be calculated using Basics::matxvec.\n";
		exit(1);
  }
  fillv0(res, sizv);
  int i;
  for (i=0; i < leng; i++)
    res[i] = dotproduct(vec, ma[i], sizv);
  return;
}

//matrices
void Basics::fillmat0(int** mat, int rows, int cols) //fills 2d array with 0s
{
  int i;
  for (i = 0; i < rows; i++)
    fillv0(mat[i], cols);
}

void Basics::fillmat0(double** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillv0(mat[i], cols);
}

void Basics::fillmat0(bool** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillv0(mat[i], cols);
}

void Basics::fillmatm1(int** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillvm1(mat[i], cols);
}

void Basics::fillmatm1(double** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)
    fillvm1(mat[i], cols);
}

bool Basics::is_symmetric(int** mat, int rows, int cols)
{
  int i, j;
  bool res = true;
  if (rows != cols)
    res = false;
  else {
    for (i = 0; i < rows; i++)
      for (j = 0; j < i; j++)
				if (mat[i][j] != mat[j][i]) {
					res = false;
					break;
				}
  }
  return res;
}

bool Basics::is_symmetric(double** mat, int rows, int cols)
{
  int i, j;
  bool res = true;
  if (rows != cols)
    res = false;
  else {
    for (i = 0; i < rows; i++)
      for (j = 0; j < i; j++)
				if (mat[i][j] != mat[j][i])	  {
					res = false;
					break;
				}
  }
  return res;
}

bool Basics::is_symmetric(bool** mat, int rows, int cols)
{
  int i, j;
  bool res = true;
  if (rows != cols)
    res = false;
  else {
    for (i = 0; i < rows; i++)
      for (j = 0; j < i; j++)
				if (mat[i][j] != mat[j][i])	  {
					res = false;
					break;
				}
  }
  return res;
}

bool Basics::is_symmetric(string** mat, int rows, int cols)
{
  int i, j;
  bool res = true;
  if (rows != cols)
    res = false;
  else {
    for (i = 0; i < rows; i++)
      for (j = 0; j < i; j++)
				if (mat[i][j] != mat[j][i])	  {
					res = false;
					break;
				}
  }
  return res;
}

bool Basics::last_matrix(bool **ve, int rows, int cols)
{
	bool res = true;
	int i;
	for (i=0; i< rows; i++)
		if (!last_vector(ve[i], cols)) {
			res = false;
			break;
		}
	return res;
}

bool Basics::last_matrix(int **ve, int rows, int cols, int max)
{
	bool res = true;
	int i;
	for (i=0; i< rows; i++)
		if (!last_vector(ve[i], cols, max)) {
			res = false;
			break;
		}
	return res;
}

void Basics::next_in_matrix(bool **ve, int rows, int cols)
{
	int i,j;
	if (last_matrix(ve, rows, cols)) {
		cout << "[Error]: Last matrix already reached in Basics::next_in_matrix.\n";
		exit(1);
	}
	for (i=0; i<rows; i++)
		if (!last_vector(ve[i], cols)) {
			next_in_vector(ve[i], cols);
			break;
		}
	for (j=0; j<i; j++)
		fillv0(ve[j], cols);
}


void Basics::next_in_matrix(int **ve, int rows, int cols, int min, int max, int step)
{
	int i,j,k;
	if (last_matrix(ve, rows, cols, max)) {
		cout << "[Error]: Last matrix already reached in Basics::next_in_matrix.\n";
		exit(1);
	}
	for (i=0; i<rows; i++)
		if (!last_vector(ve[i], cols, max)) {
			next_in_vector(ve[i], cols, min, max, step);
			break;
		}
	for (j=0; j<i; j++)
		for (k=0; k < cols; k++)
			ve[j][k] = min;
}

int Basics::sum_column(int **mat, int rows, int cols, int col)
{
	int res = 0;
	if (cols<=col) {
		cout << "[Error]: Non-existent column in Basics::sum_column . \n";
		exit(1);
	}
	int i;
	for (i=0; i<rows; i++)
		res = res+ mat[i][col];
	return res;
}

double Basics::sum_column(double **mat, int rows, int cols, int col)
{
	double res = 0;
	if (cols<=col) {
		cout << "[Error]: Non-existent column in Basics::sum_column . \n";
		exit(1);
	}
	int i;
	for (i=0; i<rows; i++)
		res = res+ mat[i][col];
	return res;
}

bool Basics::eqmatrix(int **mat1, int r1, int c1, int **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

bool Basics::eqmatrix(double **mat1, int r1, int c1, double **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

bool Basics::eqmatrix(char **mat1, int r1, int c1, char **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

bool Basics::eqmatrix(bool **mat1, int r1, int c1, bool **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

bool Basics::eqmatrix(string **mat1, int r1, int c1, string **mat2, int r2, int c2)
{
  bool res = true;
  int i;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (!eqvec(mat1[i], c1, mat2[i], c2)) {
        res = false;
        break;
      }
  }
  return res;
}

bool Basics::cycle_in_matrix_reps(int **mat, int rows, int cols, int reps) {
  if ((rows%reps) != 0) {
    cout << "[Error]: period not possible for assessed matrix when calling Basics::cycle_in_matrix_per.\n";
    exit(1);
  }
  int per = rows/reps;
  int i,j;
  bool res = true;
  for (i = 0; i < per; i++) {
    for (j = 1; j < reps; j++) {
      if (!eqvec(mat[i], cols, mat[(per*j)+i], cols)) {
        res = false;
        break;
      }
    }
    if (!res)
      break;
  }
  
  return res;
}

bool Basics::eqvector_rot(int *v1, int s1, int *v2, int s2){ //2024
  bool res = true;
  int i, ini = -1;
  if (s1 != s2)
    res = false;
  else {
    for (i=0; i < s1; i++)
      if (v1[0] == v2[i]) {
        ini = i;
        break;
      }
    if (ini < 0)
      res = false;
    else {
      for (i=0; i < s1; i++)
        if (v1[i] != v2[(ini+i)%s2]) {
          res = false;
          break;
        }
    }
  }
  return res;
}

bool Basics::eqmatrix_rot(int **mat1, int r1, int c1, int **mat2, int r2, int c2)
{
  bool res = true;
  int i,ini=-1;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (eqvec(mat1[0], c1, mat2[i], c2)) {
        ini = i;
        break;
      }
    if (ini < 0)
      res = false;
    else {
      for (i=0; i<r1; i++)
        if (!eqvec(mat1[i], c1, mat2[(ini+i)%r2], c2)) {
          res = false;
          break;
        }
    }
  }
  return res;
}

bool Basics::eqmatrix_rot(double **mat1, int r1, int c1, double **mat2, int r2, int c2)
{
  bool res = true;
  int i,ini=-1;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (eqvec(mat1[0], c1, mat2[i], c2)) {
        ini = i;
        break;
      }
    if (ini < 0)
      res = false;
    else {
      for (i=0; i<r1; i++)
        if (!eqvec(mat1[i], c1, mat2[(ini+i)%r2], c2)) {
          res = false;
          break;
        }
    }
  }
  return res;
}

bool Basics::eqmatrix_rot(char **mat1, int r1, int c1, char **mat2, int r2, int c2)
{
  bool res = true;
  int i,ini=-1;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (eqvec(mat1[0], c1, mat2[i], c2)) {
        ini = i;
        break;
      }
    if (ini < 0)
      res = false;
    else {
      for (i=0; i<r1; i++)
        if (!eqvec(mat1[i], c1, mat2[(ini+i)%r2], c2)) {
          res = false;
          break;
        }
    }
  }
  return res;
}

bool Basics::eqmatrix_rot(bool **mat1, int r1, int c1, bool **mat2, int r2, int c2)
{
  bool res = true;
  int i,ini=-1;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (eqvec(mat1[0], c1, mat2[i], c2)) {
        ini = i;
        break;
      }
    if (ini < 0)
      res = false;
    else {
      for (i=0; i<r1; i++)
        if (!eqvec(mat1[i], c1, mat2[(ini+i)%r2], c2)) {
          res = false;
          break;
        }
    }
  }
  return res;
}

bool Basics::eqmatrix_rot(string **mat1, int r1, int c1, string **mat2, int r2, int c2)
{
  bool res = true;
  int i,ini=-1;
  if ((r1 != r2) || (c1 != c2))
    res = false;
  else {
    for (i=0; i<r1; i++)
      if (eqvec(mat1[0], c1, mat2[i], c2)) {
        ini = i;
        break;
      }
    if (ini < 0)
      res = false;
    else {
      for (i=0; i<r1; i++)
        if (!eqvec(mat1[i], c1, mat2[(ini+i)%r2], c2)) {
          res = false;
          break;
        }
    }
  }
  return res;
}

double Basics::dist_matrices_rot(int **mat1, int r1, int c1, int **mat2, int r2, int c2)
{
  if (c1 != c2) {
    cout << "[Error]: Matrices have different number of columns in Basics::dist_matrices.\n";
    exit(1);
  }
  int i,j,k;
  double dis;
  if (r1 != r2) {
    int mul1, mul2;
    if (r1 < r2) {
      if ((r2%r1)==0) {
        mul2 = 1;
        mul1 = r2/r1;
      }
      else {
        mul2 = r1;
        mul1 = r2;
      }
    }
    else {
      if ((r1%r2)==0) {
        mul1 = 1;
        mul2 = r1/r2;
      }
      else {
        mul2 = r1;
        mul1 = r2;
      }
    }
    int **num1, **num2;
    num1 = new int*[r1*mul1];
    for (i=0; i < (r1*mul1); i++)
      num1[i] = new int[c1];
    num2 = new int*[r2*mul2];
    for (i=0; i < (r2*mul2); i++)
      num2[i] = new int[c2];
    for (i=0; i < mul1; i++)
      for (j=0; j < r1; j++)
        for (k=0; k < c1; k++)
          num1[(r1*i)+j][k] = mat1[j][k];
    for (i=0; i < mul2; i++)
      for (j=0; j < r2; j++)
        for (k=0; k < c2; k++)
          num2[(r2*i)+j][k] = mat2[j][k];
    dis = dist_matrices_rot_aux(num1, (r1*mul1), c1, num2, (r2*mul2), c2);
    for (i=0; i < (r1*mul1); i++)
      delete [] num1[i];
    delete [] num1;
    for (i=0; i < (r2*mul2); i++)
      delete [] num2[i];
    delete [] num2;
  }
  else {
    dis = dist_matrices_rot_aux(mat1, r1, c1, mat2, r2, c2);
  }
  return dis;
}

double Basics::dist_matrices_rot_aux(int **mat1, int r1, int c1, int **mat2, int r2, int c2)
{
  if ((r1 != r2) || (c1 != c2)) {
    cout << "[Error]: Matrices are not equal-sized in Basics::dist_matrices_aux.\n";
    exit(1);
  }
  int i,j;
  double *dists, d;
  dists = new double[r1];
  fillv0(dists, r1);
  for (i=0; i < r1; i++) {
    dists[i] = 0;
    for (j=0; j < r1; j++)
      dists[i] = dists[i] + mean_dist_vecs(mat1[j], c1, mat2[(j+i)%r1], c2);
  }
  d = find_min(dists, r1);
  d = d/double(r1);
  delete [] dists;
  return d;
}

void Basics::matxmat(int **mat1, int r1, int c1, int **mat2, int r2, int c2, int **mat3, int r3, int c3) {
  if ((c1 != r2) || (r1 != r3) || (c2 != c3)) {
    cout << "[Error]: The dimensions of matrices are not appropriate for Basics::matxmat. \nr1\tc1\tr2\tc2\tr3\tc3:\n";
    cout << r1 << "\t" << c1 << "\t" << r2 << "\t" << c2 << "\t" << r3 << "\t" << c3 << "\n";
    exit(1);
  }
  int i, j, k;
  for (i = 0; i < r1; i++) {
    for (j = 0; j < c2; j++) {
      mat3[i][j] = 0;
      for (k = 0; k < r2; k++)
        mat3[i][j] = mat3[i][j] + mat1[i][k]*mat2[k][j];
    }
  }
  return;
}

void Basics::matxmat(double **mat1, int r1, int c1, int **mat2, int r2, int c2, double **mat3, int r3, int c3) {
  if ((c1 != r2) || (r1 != r3) || (c2 != c3)) {
    cout << "[Error]: The dimensions of matrices are not appropriate for Basics::matxmat. \nr1\tc1\tr2\tc2\tr3\tc3:\n";
    cout << r1 << "\t" << c1 << "\t" << r2 << "\t" << c2 << "\t" << r3 << "\t" << c3 << "\n";
    exit(1);
  }
  int i, j, k;
  for (i = 0; i < r1; i++) {
    for (j = 0; j < c2; j++) {
      mat3[i][j] = 0;
      for (k = 0; k < r2; k++)
        mat3[i][j] = mat3[i][j] + mat1[i][k]*mat2[k][j];
    }
  }
  return;
}

void Basics::matxmat(int **mat1, int r1, int c1, double **mat2, int r2, int c2, double **mat3, int r3, int c3) {
  if ((c1 != r2) || (r1 != r3) || (c2 != c3)) {
    cout << "[Error]: The dimensions of matrices are not appropriate for Basics::matxmat. \nr1\tc1\tr2\tc2\tr3\tc3:\n";
    cout << r1 << "\t" << c1 << "\t" << r2 << "\t" << c2 << "\t" << r3 << "\t" << c3 << "\n";
    exit(1);
  }
  int i, j, k;
  for (i = 0; i < r1; i++) {
    for (j = 0; j < c2; j++) {
      mat3[i][j] = 0;
      for (k = 0; k < r2; k++)
        mat3[i][j] = mat3[i][j] + mat1[i][k]*mat2[k][j];
    }
  }
  return;
}

void Basics::matxmat(double **mat1, int r1, int c1, double **mat2, int r2, int c2, double **mat3, int r3, int c3) {
  if ((c1 != r2) || (r1 != r3) || (c2 != c3)) {
    cout << "[Error]: The dimensions of matrices are not appropriate for Basics::matxmat. \nr1\tc1\tr2\tc2\tr3\tc3:\n";
    cout << r1 << "\t" << c1 << "\t" << r2 << "\t" << c2 << "\t" << r3 << "\t" << c3 << "\n";
    exit(1);
  }
  int i, j, k;
  for (i = 0; i < r1; i++) {
    for (j = 0; j < c2; j++) {
      mat3[i][j] = 0;
      for (k = 0; k < r2; k++)
        mat3[i][j] = mat3[i][j] + mat1[i][k]*mat2[k][j];
    }
  }
  return;
}

void Basics::add_to_matrix(int **growing, int **addendum, int r, int c) {
  int i,j;
  for (i = 0; i < r; i++)
    for (j = 0; j < c; j++)
      growing[i][j] = growing[i][j] + addendum[i][j];
}

void Basics::add_to_matrix(double **growing, double **addendum, int r, int c) {
  int i,j;
  for (i = 0; i < r; i++)
    for (j = 0; j < c; j++)
      growing[i][j] = growing[i][j] + addendum[i][j];
}

//matrix formatting
void Basics::printm(ostream& sal, int** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)    {
		printv(sal, mat[i], cols);
		sal << endl;
	}
}

void Basics::printm(ostream& sal, double** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)    {
		printv(sal, mat[i], cols);
		sal << endl;
	}
}

void Basics::printm(ostream& sal, bool** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)    {
		printv(sal, mat[i], cols);
		sal << endl;
	}
}

void Basics::printm(ostream& sal, string** mat, int rows, int cols)
{
  int i;
  for (i = 0; i < rows; i++)    {
		printv(sal, mat[i], cols);
		sal << endl;
	}
}

void Basics::printm_latex(ostream& sal, int** mat, int rows, int cols)
{
  int i,j;
  sal << "\\left[ \\begin{array}{";
  for (i= 0; i < cols; i++)
    sal << "c";
  sal << "}\n";
  for (i=0; i < rows; i++) {
    sal << mat[i][0];
    for (j=1; j < cols; j++)
      sal << "&" << mat[i][j];
    if (i < (rows-1))
      sal << "\\\\ \n";
  }
  sal << "\n \\end{array} \\right]\n";
}

void Basics::printm_latex(ostream& sal, double** mat, int rows, int cols)
{
  int i,j;
  sal << "\\left[ \\begin{array}{";
  for (i= 0; i < cols; i++)
    sal << "c";
  sal << "}\n";
  for (i=0; i < rows; i++) {
    sal << mat[i][0];
    for (j=1; j < cols; j++)
      sal << "&" << mat[i][j];
    if (i < (rows-1))
      sal << "\\\\ \n";
  }
  sal << "\n \\end{array} \\right]\n";
}

//sets
set<int> Basics::merge(set<int> &a, set<int> &b)
{
  set<int> res;
  set<int>::iterator It;
  res = a;
  for(It=b.begin(); It!=b.end(); It++)
    res.insert(*It);
  return res;
}

set<char> Basics::merge(set<char> &a, set<char> &b)
{
  set<char> res;
  res.clear();
  set<char>::iterator It;
  res = a;
  for(It=b.begin(); It!=b.end(); It++)
    res.insert(*It);
  return res;
}

set<string> Basics::merge(set<string> &a, set<string> &b)
{
  set<string> res;
  res.clear();
  set<string>::iterator It;
  res = a;
  for(It=b.begin(); It!=b.end(); It++)
    res.insert(*It);
  return res;
}

set<int> Basics::intersect(set<int> &a, set<int> &b)
{
  set<int> res;
  res.clear();
  set<int>::iterator It;
  for(It=b.begin(); It!=b.end(); It++)
		if (contains(*It, a))
			res.insert(*It);
  return res;
}

set<char> Basics::intersect(set<char> &a, set<char> &b)
{
  set<char> res;
  res.clear();
  set<char>::iterator It;
  for(It=b.begin(); It!=b.end(); It++)
		if (contains(*It, a))
			res.insert(*It);
  return res;
}

set<string> Basics::intersect(set<string> &a, set<string> &b)
{
  set<string> res;
  res.clear();
  set<string>::iterator It;
  for(It=b.begin(); It!=b.end(); It++)
		if (contains(*It, a))
			res.insert(*It);
  return res;
}

void Basics::printset(set<int> &ise, ostream& sal)
{
  set<int>::iterator it;
  sal << "(";
  for (it=ise.begin(); it!=ise.end(); it++)
    sal << " " << *it;
  sal << ")";
}

void Basics::printset(set<string> &ise, ostream& sal)
{
  set<string>::iterator it;
  sal << "(";
  for (it=ise.begin(); it!=ise.end(); it++)
    sal << " " << *it;
  sal << ")";
}

void Basics::printsetofsets(set<set<int> > &con, ostream& sal)
{
  set<set<int> >::iterator whi;
  set<int> unc;
  sal << "[";
  for (whi=con.begin(); whi!=con.end(); whi++) {
    unc = *whi;
    printset(unc, sal);
  }
  sal << "]";
}

bool Basics::contains(char c, set<char> &sc)
{
	bool res = false;
	set<char>::iterator it;
	for (it=sc.begin(); it !=sc.end(); it++) {
		if (*it == c) {
			res = true;
			break;
		}
	}
	return res;
}

bool Basics::contains(int c, set<int> &sc)
{
	bool res = false;
	set<int>::iterator it;
	for (it=sc.begin(); it !=sc.end(); it++) {
		if (*it == c) {
			res = true;
			break;
		}
	}
	return res;
}

bool Basics::contains(string c, set<string> &sc)
{
	bool res = false;
	set<string>::iterator it;
	for (it=sc.begin(); it !=sc.end(); it++) {
		if (*it == c) {
			res = true;
			break;
		}
	}
	return res;
}

void Basics::remove_nth(set<string>& conj, int n)
{
	set<string>::iterator it = conj.begin();
	int i = 0;
	for (i=0; i<n; i++)
		++it;
	conj.erase(it);
}

string Basics::return_nth(set<string> &conj, int n)
{
	set<string>::iterator it = conj.begin();
	int i = 0;
	for (i=0; i<n; i++)
		++it;
	return *it;
}

void Basics::append_sets_of_sets(const set<set<int> > &source, set<set<int> > &target) {
  set<set<int> >::iterator It;
  for (It = source.begin(); It != source.end(); It++)
    target.insert(*It);
}

//lists
bool Basics::contains(int num, list<int> &lis)
{
  bool res = false;
  list<int>::iterator it;
  for (it = lis.begin(); it!= lis.end(); it++)
    if (*it == num) {
			res = true;
			break;
		}
  return res;
}

bool Basics::contains(string num, list<string> &lis)
{
  bool res = false;
  list<string>::iterator it;
  for (it = lis.begin(); it!= lis.end(); it++)
    if (*it == num) {
			res = true;
			break;
		}
  return res;
}

void Basics::printlist(list<int> &lis, ostream& sal)
{
	list<int>::iterator it;
	sal << "(";
	for (it=lis.begin(); it!=lis.end(); it++)
		sal << *it << " ";
	sal << ")";
}

void Basics::printlist(list<string> &lis, ostream& sal)
{
	list<string>::iterator it;
	sal << "(";
	for (it=lis.begin(); it!=lis.end(); it++)
		sal << *it << " ";
	sal << ")";
}

void Basics::printlistoflists(list<list<int> > &lis, ostream& sal)
{
	list<list<int> >::iterator whi;
	for (whi=lis.begin(); whi!=lis.end(); whi++)
		printlist(*whi, sal);
}

double Basics::productory(list<int> con)
{
	list<int>::iterator it;
	double res = 1;
	for (it=con.begin(); it != con.end(); it++)
		res = res*(*it);
	return res;
}

int Basics::factorial(int num)
{
	int i,res = 1;
	for (i=1; i<=num; i++)
		res= res*i;
	return res;
}
double Basics::binom_coeff(int ene, int ka)
{
	int i,j;
	bool change;
	double nume, deno;
	list<int> num;
	list<int> den;
	num.clear();
	den.clear();
	for (i=(ene-ka+1); i <= ene; i++)
		num.push_back(i);
	for (i = 1; i <= ka; i++)
		den.push_back(i);
	list<int>::iterator itn;
	list<int>::iterator itd;
	do {
		change = false;
		for (itn =num.begin(); itn != num.end(); itn++)
			for (itd = den.begin(); itd != den.end(); itd++)
				if (*itn == *itd) {
					itn = num.erase(itn);
					itd = den.erase(itd);
					--itn;
					--itd;
					change = true;
					break;
				}
		
		for (itn =num.begin(); itn != num.end(); itn++)
			for (itd = den.begin(); itd != den.end(); itd++)
				if ((*itn%(*itd))==0) {
					j = (*itn)/(*itd);
					itn = num.erase(itn);
					itd = den.erase(itd);
					num.insert(itn, j);
					--itn;
					--itn;
					--itd;
					change = true;
					break;
				}
		
		for (itn =num.begin(); itn != num.end(); itn++)
			for (itd = den.begin(); itd != den.end(); itd++)
				if ((*itd%(*itn))==0) {
					j = (*itd)/(*itn);
					itn = num.erase(itn);
					itd = den.erase(itd);
					den.insert(itd, j);
					--itn;
					--itd;
					--itd;
					change = true;
					break;
				}
	}while(change);
	nume = productory(num);
	deno = productory(den);
	return (nume/deno);
}

//others
void Basics::run_command(string cuerda){
  system(cuerda.c_str());
  return;
}

void Basics::create_array(int** &arr, int rows, int cols) {
  int i;
  arr = new int*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new int[cols];
  return;
}

void Basics::create_array(bool** &arr, int rows, int cols) {
  int i;
  arr = new bool*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new bool[cols];
  return;
}

void Basics::create_array(double** &arr, int rows, int cols) {
  int i;
  arr = new double*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new double[cols];
  return;
}

void Basics::create_array(char** &arr, int rows, int cols) {
  int i;
  arr = new char*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new char[cols];
  return;
}

void Basics::create_array(string** &arr, int rows, int cols) {
  int i;
  arr = new string*[rows];
  for (i = 0; i < rows; i++)
    arr[i] = new string[cols];
  return;
}

void Basics::create_array(int*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new int**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new int*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new int[cols];
  }
  return;
}

void Basics::create_array(bool*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new bool**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new bool*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new bool[cols];
  }
  return;
}

void Basics::create_array(double*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new double**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new double*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new double[cols];
  }
  return;
}

void Basics::create_array(char*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new char**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new char*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new char[cols];
  }
  return;
}

void Basics::create_array(string*** &arr, int slices, int rows, int cols) {
  int i, j;
  arr = new string**[slices];
  for (i = 0; i < slices; i++) {
    arr[i] = new string*[rows];
    for (j = 0; j < rows; j++)
      arr[i][j] = new string[cols];
  }
  return;
}

void Basics::get_dot_fig(string arch)
{
	char s[250];
	snprintf(s, 250, "dot -Tpdf %s.dot -o %s.pdf", arch.c_str(), arch.c_str());
	system(s);
}

void Basics::get_dot_fig(string arch, string ext)
{
	if (ext=="pdf")
		get_dot_fig(arch);
	else {
    if (ext=="eps") {
      char s[250];
      snprintf(s, 250, "dot -Tps2 %s.dot -o %s.eps", arch.c_str(), arch.c_str());
      system(s);
    }
		else {
			cout << "[Error]: Basics::get_dot_fig does not recognize extension " << ext << ".\n";
			exit(1);
		}
	}
}

void Basics::get_spring_dot_fig(string arch)
{
	char s[250];
	snprintf(s, 250, "fdp -Tps2 %s.dot -o %s.eps", arch.c_str(), arch.c_str());
	system(s);
	snprintf(s, 250, "pstopdf %s.eps -o %s.pdf", arch.c_str(), arch.c_str());
	system(s);
	snprintf(s, 250, "rm %s.eps", arch.c_str());
	system(s);
}

void Basics::get_spring_dot_fig(string arch, string ext)
{
	if (ext=="pdf")
		get_spring_dot_fig(arch);
	else {
		if (ext=="eps") {
			char s[250];
			snprintf(s, 250, "fdp -Tps2 %s.dot -o %s.eps", arch.c_str(), arch.c_str());
			system(s);
		}
		else {
			cout << "[Error]: Basics::get_dot_fig does not recognize extension " << ext << ".\n";
			exit(1);
		}
	}
}

int Basics::round(double x)
{
	if (x<0)
		return int(x-0.5);
	else
		return int(x+0.5);
}

void Basics::polar_to_cartesian(double radians, double radius, double& x, double& y)
{
	double ton = tan(radians);
	double den =1+(ton*ton);
	den = sqrt(den);
	x = radius/den;
	double PI = 3.14159265;
	double xpos;
	if ((radians > (PI/2.0)) && (radians < (3*(PI/2.0))))
		xpos = -1;
	else
		xpos = 1;
	if ((xpos*x) < 0)
		x = x*xpos;
	y = ton*x;
}

void Basics::open_ifstream(ifstream& fe, string nomb) //opens input file stream
{
	fe.open(nomb.c_str());
	if (!fe.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Basics::open_ifstream.\n";
		exit(1);
	}
}

void Basics::open_ofstream(ofstream& fs, string nomb) //opens output file stream
{
	fs.open(nomb.c_str());
	if (!fs.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Basics::open_ofstream.\n";
		exit(1);
	}
}

void Basics::open_ofstream_to_append(ofstream& fs, string nomb)
{
	fs.open(nomb.c_str(), ios::app);
	if (!fs.is_open()) {
		cout << "[Error]: File named \'" << nomb << "\' could not be opened using Basics::open_ofstream.\n";
		exit(1);
	}
}

bool Basics::dir_exists(string dir)
{
  struct stat buf;
  int esmu = stat(dir.c_str(), & buf);
  if ((esmu >= 0) && (S_ISDIR(buf.st_mode)))
    return true;
  else
    return false;
}


int Basics::count_files_in_dir(string dir)
{
	int res = -2;
  struct stat buf;
  int esmu = stat(dir.c_str(), & buf);
  if ((esmu >= 0) && (S_ISDIR(buf.st_mode))) {
    DIR *eldi;
    struct dirent *qbarbol;
    eldi = opendir(dir.c_str());
    if (eldi == NULL) {
      cout << "[Error]: Directory " << dir << " cannot be opened by Basics::count_files_in_dir.\n";
      exit(1);
    }
    while ((qbarbol = readdir(eldi)) != NULL) {
      res++;
    }
    closedir(eldi);
  }
  else {
    cout << "[Error]: Directory " << dir << " does not exist. Basics::count_files_in_dir.\n";
    exit(1);
  }
	return res;
}

void Basics::get_files_in_dir(string dir, string* vecnom)
{
  int i =0;
  string mi;
  struct stat buf;
  int esmu = stat(dir.c_str(), & buf);
  if ((esmu >= 0) && (S_ISDIR(buf.st_mode))) {
    DIR *eldi;
    struct dirent *qbarbol;
    eldi = opendir(dir.c_str());
    if (eldi == NULL) {
      cout << "[Error]: Directory " << dir << " cannot be opened by Basics::get_files_in_dir.\n";
      exit(1);
    }
    while ((qbarbol = readdir(eldi)) != NULL) {
      mi = qbarbol->d_name;
      if ((mi != ".") && (mi != "..")) {
        vecnom[i] = mi;
        i++;
      }
    }
    closedir(eldi);
  }
  else {
    cout << "[Error]: Directory " << dir << " does not exist. Basics::get_files_in_dir.\n";
    exit(1);
  }
}

//strings
void Basics::erase_between_braces(char pator[], int len, ostream& fs, ifstream& fe)
{
	int i, cont=0, j;
	char *s, basu;
	s = new char[len];
	for (i=0; i<len; i++)
		if (pator[i] == '{') {
			cont++;
			break;
		}
	if (cont==0) {
		cout << "[Error]: Pattern does not contain left brace in Basics::erase_between_braces.\n";
		exit(1);
	}
	cont = 0;
	while (!fe.eof()) {
		j=0;
		do {
			fe.get(s[j]);
			if (s[j] != pator[j]) {
				for (i=0; i<=j; i++)
					fs << s[i];
				break;
			}
			j++;
		} while (j<len);
		if (j==len) {
			for (i=0; i<len; i++)
				if (s[i] != pator[i]) {
					cout << "[Error]: Error in erase_between_braces.\n" ;
					exit(1);
				}
			cont = 1;
			do {
				fe.get(basu);
				if (basu == '{')
					cont++;
				if (basu == '}')
					cont--;
			} while (cont > 0);
		}
	}
	return;
}

int Basics::char_in_string(char c, const string& cue, int from, int until)
{
	int j,i = -1;
	for (j=from; j<until; j++) {
		if (cue[j]==c) {
			i = j;
			break;
		}
	}
	return i;
}

int Basics::char_in_string(char c, const string& cue)
{
	int u, res;
	u = cue.length();
	res = char_in_string(c, cue, 0, u);
	return res;
}

string Basics::del_char(char c, string ori)
{
	string res = ori;
	int j = ori.length();
	int i = char_in_string(c, res, 0, j);
	while (i >= 0) {
		res.replace(i,1,"");
		i = char_in_string(c, res, 0, j);
	}
	return res;
}

int Basics::count_char_in_string(char c, const string& cue)
{
	int res=0, i,l=cue.length();
	for (i=0; i<l; i++) {
		if (cue[i]==c)
			res++;
	}
	return res;
}

string Basics::chop_word_from_line(string& linea)
{
	string res="";
	int l,i;
	
	if (!linea.empty()) {
		while ((linea[0]==' ') || (linea[0]=='\t')) {
			l=linea.length();
			linea = linea.substr(1, l-1);
		}
		if (!linea.empty()) {
			l=linea.length();
			for (i=0; i<l; i++)
				if ((linea[i]==' ') || (linea[i]=='\t'))
					break;
			if (i < l) {
        res = linea.substr(0,i);
        linea = linea.substr(i+1,(l-(i+1)));
			}
			else {
				res = linea;
				linea = "";
			}
		}
	}
	return res;
}

bool Basics::string_in_string(string pat, string big)
{
	int i,le=pat.length();
	int leb = big.length();
	bool res=false;
	for (i=0; i<=(leb-le); i++) //before 26.01.2012 it was <. corrected so that last char in bih is searchable
		if (big[i]==pat[0])
			if ((pat.substr(1, (le-1)))==(big.substr((i+1), (le-1))))
				res = true;
	return res;
}

string Basics::inttostring(int num)
{
	char buff[50];
	int j;
	j = snprintf(buff, 50, "%d", num);
	if (j <0) {
		cout << "[Error]: Transformation of string to int was not possible using Basics::inttostring.\n";
		exit(1);
	}
	string res(buff);
	return res;
}

string Basics::doubletostring_tex(double num)
{
	char buff[50];
	int j,k,l;
	string esp,sig,resf;
	j = snprintf(buff, 50, "%g", num);
	if (j <0) {
		cout << "[Error]: Transformation of string to double was not possible using Basics::doubletostring.\n";
		exit(1);
	}
	string res(buff);
	j = char_in_string('e', res);
	if (j==(-1))
		j = char_in_string('E', res);
	if (j>=0) {
		l = res.size();
		k= char_in_string('-', res, 1, l);
		esp = res.substr(j+1, l-(j+1));
		if (esp[0]=='-')
			sig = "-";
		else {
			if (esp[0]!='+') {
				cout << "[Error]: No + sign in exponent. The error was found while performing Basics::doubletostring.\n";
				exit(1);
			}
			sig = "";
		}
		esp = esp.substr(1, (l-(j+2)));
		resf = "$"+res.substr(0, j) +"\\times10^{"+sig+esp+"}$";
	}
	else {
		k=char_in_string('.', res);
		if (k==(-1))
			resf = res;
		else {
			l = res.size();
			while ((l>1) &&(res[l-1]==0)) {
				res = res.substr(0, (l-1));
				l = res.size();
			}
			resf = res;
		}
	}
	return resf;
}

bool Basics::string_in_vector(string* vec, string ele, int tam)
{
	int i;
	bool res = false;
	for (i=0; i<tam; i++) {
		if (vec[i] == ele) {
			res = true;
			break;
		}
	}
	return res;
}

int Basics::where_string_in_vector(string* vec, string ele, int tam)
{
	int i=-1;
	for (i=0; i<tam; i++)
		if (vec[i] == ele)
			break;
	if ((i<0) || (i==tam)) {
		cout << "[Error]: String not found by Basics::where_string_in_vector.\n";
		exit(1);
	}
	
	return i;
}

bool Basics::nothing_but_spaces(const string& cue)
{
	bool res = true;
	int i, tam;
	tam = cue.size();
	if (!cue.empty()) {
		for (i=0; i<tam; i++) {
			if (isgraph(cue[i])) {//(cue[i] != ' ') && (cue[i] != '\t')) {
				res = false;
				break;
			}
		}
	}
	return res;
}

void Basics::capitalize(string& word)
{
	int i,l;
	int j;
	l = word.length();
	word[0] = toupper(word[0]);
	for (i=1; i<l; i++)
		word[i] = tolower(word[i]);
	
	j = word.find("Á", 0);
	while (j != -1) {
		word.replace(j, 2, "á", 2);
		j = word.find("Á", j);
	}
	j = word.find("É", 0);
	while (j != -1) {
		word.replace(j, 2, "é", 2);
		j = word.find("É", j);
	}
	j = word.find("Í", 0);
	while (j != -1) {
		word.replace(j, 2, "í", 2);
		j = word.find("Í", j);
	}
	j = word.find("Ó", 0);
	while (j != -1) {
		word.replace(j, 2, "ó", 2);
		j = word.find("Ó", j);
	}
	j = word.find("Ú", 0);
	while (j != -1) {
		word.replace(j, 2, "ú", 2);
		j = word.find("Ú", j);
	}
}

void Basics::capitalize_sentence(string& sent)
{
	int i, l;
	int j;
	l = sent.length();
	capitalize(sent);
	for (i=0; i<(l-1); i++) {
		if ((sent[i]==' ') && (sent[i+2] !=' ') && ((i+2) <l))
			sent[i+1] = toupper(sent[i+1]);
	}
	j = sent.find("á", 0);
	while (j != -1) {
		if (sent[j-1]==' ')
			sent.replace(j, 2, "Á", 2);
		j = sent.find("á", j+1);
	}
	j = sent.find("é", 0);
	while (j != -1) {
		if (sent[j-1]==' ')
			sent.replace(j, 2, "É", 2);
		j = sent.find("é", j+1);
	}
	j = sent.find("í", 0);
	while (j != -1) {
		if (sent[j-1]==' ')
			sent.replace(j, 2, "Í", 2);
		j = sent.find("í", j+1);
	}
	j = sent.find("ó", 0);
	while (j != -1) {
		if (sent[j-1]==' ')
			sent.replace(j, 2, "Ó", 2);
		j = sent.find("ó", j+1);
	}
	j = sent.find("ú", 0);
	while (j != -1) {
		if (sent[j-1]==' ')
			sent.replace(j, 2, "Ú", 2);
		j = sent.find("ú", j+1);
	}
}

bool Basics::eq_string(string uno, string dos, bool capsmatter)
{
	bool res = true;
	if (capsmatter) {
		if (uno!=dos)
			res = false;
	}
	else {
		if (uno.length()!=dos.length())
			res = false;
		else {
			int lon = uno.length();
			int i;
			char cu, cd;
			for (i=0; i<lon; i++) {
				if (uno[i]!=dos[i]) {
					if (isalpha(uno[i]) && isalpha(dos[i])) {
						cu = toupper(uno[i]);
						cd = toupper(dos[i]);
						if (cu!=cd) {
							res=false;
							break;
						}
					}
					else {
						res = false;
						break;
					}
				}
			}
		}
	}
  
	return res;
}

int Basics::nth_app_of_char(char c, const string& cue, int n)
{
	int i,j,ta,res = -1;
	j= char_in_string(c, cue);
	ta = cue.length();
	if (j != (-1)) {
		for (i=1; i<n; i++) {
			j = char_in_string(c, cue, j+1, ta);
			if (j==(-1))
				break;
		}
	}
	res = j;
	return res;
}


//statistics
double Basics::get_mean(int* vec, int tam)
{
  double res = 0;
  int i;
  for (i = 0; i < tam; i++)
    res = res + vec[i];
  res = res/(tam*1.0);
  return res;
}

double Basics::get_mean(double* vec, int tam)
{
  double res = 0;
  int i;
  for (i = 0; i < tam; i++)
    res = res + vec[i];
  res = res/(tam*1.0);
  return res;
}

double Basics::get_mean(int* vec, int tam, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += double(vec[*it]);
  res /= double(cua.size());
  return res;
}

double Basics::get_mean(double* vec, int tam, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += vec[*it];
  res /= double(cua.size());
  return res;
}

double Basics::get_mean_by_col(int **mat, int rows, int cols, int col)
{
  if (col >= cols) {
    cout << "[Error]: Column does not exist in matrix. Basics::get_mean_by_col.\n";
    exit(1);
  }
  double res = 0;
  int i;
  for (i= 0; i< rows; i++)
    res = res + double(mat[i][col]);
  res = res/double(rows);
  return res;
}

double Basics::get_mean_by_col(double **mat, int rows, int cols, int col)
{
  if (col >= cols) {
    cout << "[Error]: Column does not exist in matrix. Basics::get_mean_by_col.\n";
    exit(1);
  }
  double res = 0;
  int i;
  for (i= 0; i< rows; i++)
    res = res + mat[i][col];
  res = res/double(rows);
  return res;
}

double Basics::get_sample_variance(int* vec, int tam)
{
  double mean = get_mean(vec, tam);
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((mean - vec[i])*(mean - vec[i]));
  res = sum/(tam-1.0);
  return res;
}

double Basics::get_sample_variance(double* vec, int tam)
{
  double mean = get_mean(vec, tam);
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((mean - vec[i])*(mean - vec[i]));
  res = sum/(tam-1.0);
  return res;
}

double Basics::get_sample_variance(int* vec, int tam, set<int> &cua) {
  double ave = get_mean(vec, tam, cua);
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()-1.0);
  return res;
}

double Basics::get_sample_variance(double* vec, int tam, set<int> &cua) {
  double ave = get_mean(vec, tam, cua);
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()-1.0);
  return res;
}

double Basics::get_sample_variance(int* vec, int tam, double ave)
{
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((ave - vec[i])*(ave - vec[i]));
  res = sum/(tam-1.0);
  return res;
}

double Basics::get_sample_variance(double* vec, int tam, double ave)
{
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((ave - vec[i])*(ave - vec[i]));
  res = sum/(tam-1.0);
  return res;
}

double Basics::get_sample_variance(int* vec, int tam, double ave, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()-1.0);
  return res;
}

double Basics::get_sample_variance(double* vec, int tam, double ave, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()-1.0);
  return res;
}

double Basics::get_pop_variance(int* vec, int tam)
{
  double mean = get_mean(vec, tam);
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((mean - vec[i])*(mean - vec[i]));
  res = sum/(tam*1.0);
  return res;
}

double Basics::get_pop_variance(double* vec, int tam)
{
  double mean = get_mean(vec, tam);
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((mean - vec[i])*(mean - vec[i]));
  res = sum/(tam*1.0);
  return res;
}

double Basics::get_pop_variance(int* vec, int tam, set<int> &cua) {
  double ave = get_mean(vec, tam, cua);
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()*1.0);
  return res;
}

double Basics::get_pop_variance(double* vec, int tam, set<int> &cua) {
  double ave = get_mean(vec, tam, cua);
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()*1.0);
  return res;
}

double Basics::get_pop_variance(int* vec, int tam, double ave)
{
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((ave - vec[i])*(ave - vec[i]));
  res = sum/(tam*1.0);
  return res;
}

double Basics::get_pop_variance(double* vec, int tam, double ave)
{
  double sum = 0;
  double res;
  int i;
  for (i = 0; i < tam; i++)
    sum = sum + ((ave - vec[i])*(ave - vec[i]));
  res = sum/(tam*1.0);
  return res;
}

double Basics::get_pop_variance(int* vec, int tam, double ave, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()*1.0);
  return res;
}

double Basics::get_pop_variance(double* vec, int tam, double ave, set<int> &cua) {
  double res = 0;
  set<int>::iterator it;
  for(it = cua.begin(); it != cua.end(); it++)
    res += ((ave -vec[*it])*(ave -vec[*it]));
  res /= (cua.size()*1.0);
  return res;
}

double Basics::get_sample_stddev(int* vec, int tam)
{
  double var = get_sample_variance(vec, tam);
	double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(double* vec, int tam)
{
  double var = get_sample_variance(vec, tam);
	double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(int* vec, int tam, set<int> &cua)
{
  double var = get_sample_variance(vec, tam, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(double* vec, int tam, set<int> &cua)
{
  double var = get_sample_variance(vec, tam, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(int* vec, int tam, double ave)
{
  double var = get_sample_variance(vec, tam, ave);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(double* vec, int tam, double ave)
{
  double var = get_sample_variance(vec, tam, ave);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(int* vec, int tam, double ave, set<int> &cua)
{
  double var = get_sample_variance(vec, tam, ave, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stddev(double* vec, int tam, double ave, set<int> &cua)
{
  double var = get_sample_variance(vec, tam, ave, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(int* vec, int tam)
{
  double var = get_pop_variance(vec, tam);
	double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(double* vec, int tam)
{
  double var = get_pop_variance(vec, tam);
	double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(int* vec, int tam, set<int> &cua)
{
  double var = get_pop_variance(vec, tam, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(double* vec, int tam, set<int> &cua)
{
  double var = get_pop_variance(vec, tam, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(int* vec, int tam, double ave)
{
  double var = get_pop_variance(vec, tam, ave);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(double* vec, int tam, double ave)
{
  double var = get_pop_variance(vec, tam, ave);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(int* vec, int tam, double ave, set<int> &cua)
{
  double var = get_pop_variance(vec, tam, ave, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_pop_stddev(double* vec, int tam, double ave, set<int> &cua) {
  double var = get_pop_variance(vec, tam, ave, cua);
  double res = sqrt(var);
  return res;
}

double Basics::get_sample_stderr(int* vec, int tam)
{
  double stdv = get_sample_stddev(vec, tam);
  double den = sqrt(1.0*tam);
  double res = stdv/den;
  return res;
}

double Basics::get_sample_stderr(double* vec, int tam)
{
  double stdv = get_sample_stddev(vec, tam);
  double den = sqrt(1.0*tam);
  double res = stdv/den;
  return res;
}

double Basics::get_sample_stderr(int* vec, int tam, set<int> &cua)
{
  double stdv = get_sample_stddev(vec, tam, cua);
  double den = sqrt(1.0*cua.size());
  double res = stdv/den;
  return res;
}

double Basics::get_sample_stderr(double* vec, int tam, set<int> &cua)
{
  double stdv = get_sample_stddev(vec, tam, cua);
  double den = sqrt(1.0*cua.size());
  double res = stdv/den;
  return res;
}

double Basics::get_sample_stderr(int* vec, int tam, double ave)
{
  double stdv = get_sample_stddev(vec, tam, ave);
  double den = sqrt(1.0*tam);
  double res = stdv/den;
  return res;
}

double Basics::get_sample_stderr(double* vec, int tam, double ave)
{
  double stdv = get_sample_stddev(vec, tam, ave);
  double den = sqrt(1.0*tam);
  double res = stdv/den;
  return res;
}

double Basics::get_sample_stderr(int* vec, int tam, double ave, set<int> &cua)
{
  double stdv = get_sample_stddev(vec, tam, ave, cua);
  double den = sqrt(1.0*cua.size());
  double res = stdv/den;
  return res;
}

double Basics::get_sample_stderr(double* vec, int tam, double ave, set<int> &cua)
{
  double stdv = get_sample_stddev(vec, tam, ave, cua);
  double den = sqrt(1.0*cua.size());
  double res = stdv/den;
  return res;
}

double Basics::get_samp_var_err_of_mean(int* vec, int tam)
{
	double res;
	double var = get_sample_variance(vec, tam);
	res = var/double(tam);
	return res;
}

double Basics::get_samp_var_err_of_mean(double* vec, int tam)
{
	double res;
	double var = get_sample_variance(vec, tam);
	res = var/double(tam);
	return res;
}

double Basics::get_samp_var_err_of_mean(int* vec, int tam, set<int> &cua)
{
  double res;
  double var = get_sample_variance(vec, tam, cua);
  res = var/double(cua.size());
  return res;
}

double Basics::get_samp_var_err_of_mean(double* vec, int tam, set<int> &cua)
{
  double res;
  double var = get_sample_variance(vec, tam, cua);
  res = var/double(cua.size());
  return res;
}

double Basics::find_min(double* vec, int tam)
{
  double res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] < res)
      res = vec[i];
  return res;
}

int Basics::find_min(int* vec, int tam)
{
  int res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] < res)
      res = vec[i];
  return res;
}

double Basics::find_min(double* vec, int tam, set<int> &cua)
{
  if (cua.size() == 0) {
    cout << "[Error]: Null set in Basics::find_min\n";
    exit(1);
  }
  set<int>::iterator it;
  it = cua.begin();
  double res = vec[*it];
  it++;
  while (it != cua.end()) {
    if (vec[*it] < res)
      res = vec[*it];
    it++;
  }
  return res;
}

int Basics::find_min(int* vec, int tam, set<int> &cua)
{
  if (cua.size() == 0) {
    cout << "[Error]: Null set in Basics::find_min\n";
    exit(1);
  }
  set<int>::iterator it;
  it = cua.begin();
  int res = vec[*it];
  it++;
  while (it != cua.end()) {
    if (vec[*it] < res)
      res = vec[*it];
    it++;
  }
  return res;
}

int Basics::find_minlb(int* vec, int tam, int lbo)
{
	int min = 50000, i, cual = -1;
	for (i=0; i<tam; i++) {
		if ((vec[i] >= lbo) && (vec[i] < min)) {
			cual = i;
			min = vec[i];
		}
	}
	if (cual <0) {
		cout << "[Error]: Vector does not contain any number larger than " << lbo << " when calling Basics::find_minlb.\n";
		exit(1);
	}
	return cual;
}

int Basics::find_minlb(double* vec, int tam, double lbo)
{
	int i, cual = -1;
	double min = 50000;
	for (i=0; i<tam; i++) {
		if ((vec[i] >= lbo) && (vec[i] < min)) {
			cual = i;
			min = vec[i];
		}
	}
	if (cual <0) {
		cout << "[Error]: Vector does not contain any number larger than " << lbo << " when calling Basics::find_minlb.\n";
		exit(1);
	}
	return cual;
}

double Basics::find_max(double* vec, int tam)
{
  double res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] > res)
      res = vec[i];
  return res;
}

int Basics::find_max(int* vec, int tam)
{
  int res = vec[0];
  int i;
  for (i = 1; i < tam; i++)
    if (vec[i] > res)
      res = vec[i];
  return res;
}

double Basics::find_max(double* vec, int tam, set<int> &cua)
{
  if (cua.size() == 0) {
    cout << "[Error]: Null set in Basics::find_max\n";
    exit(1);
  }
  set<int>::iterator it;
  it = cua.begin();
  double res = vec[*it];
  it++;
  while (it != cua.end()) {
    if (vec[*it] > res)
      res = vec[*it];
    it++;
  }
  return res;
}

int Basics::find_max(int* vec, int tam, set<int> &cua)
{
  if (cua.size() == 0) {
    cout << "[Error]: Null set in Basics::find_max\n";
    exit(1);
  }
  set<int>::iterator it;
  it = cua.begin();
  int res = vec[*it];
  it++;
  while (it != cua.end()) {
    if (vec[*it] > res)
      res = vec[*it];
    it++;
  }
  return res;
}

int Basics::find_max_index(double* vec, int tam)
{
  int res = 0,i;
  for (i=1; i< tam; i++)
    if (vec[i] > vec[res])
      res = i;
  return res;
}

int Basics::find_max_index(int* vec, int tam)
{
  int res = 0,i;
  for (i=1; i< tam; i++)
    if (vec[i] > vec[res])
      res = i;
  return res;
}

int Basics::find_min_index(double* vec, int tam)
{
  int res = 0,i;
  for (i=1; i< tam; i++)
    if (vec[i] < vec[res])
      res = i;
  return res;
}

int Basics::find_min_index(int* vec, int tam)
{
  int res = 0,i;
  for (i=1; i< tam; i++)
    if (vec[i] < vec[res])
      res = i;
  return res;
}

void Basics::sort(int* vec, int* nvec, int tam)
{
	int min, max;
	int i, j, quedan, k;
	int* wov;
	wov = new int[tam];
	int* wov2;
	for (i=0; i < tam; i++)
		wov[i] = vec[i];
	max = find_max(wov, tam);
	quedan = tam;
	i = 0;
	while (quedan > 0) {
		min = find_min(wov, quedan);
		for (j=0; j< quedan; j++) {
			if (wov[j] == min) {
				nvec[i] = wov[j];
				i++;
				break;
			}
		}
		wov2 = new int[quedan-1];
		for (k = 0; k < j; k++)
			wov2[k] = wov[k];
		for (k=j; k < (quedan-1); k++)
			wov2[k] = wov[k+1];
		delete [] wov;
		quedan--;
		wov = new int[quedan];
		for (k=0; k < quedan; k++)
			wov[k] = wov2[k];
		delete [] wov2;
	}
	if (min != max) {
		cout << "[Error]: sort does not converge in Basics::sort.\n";
		exit(1);
	}
	delete [] wov;
	return;
}

void Basics::sort(int* vec, int* nvec, int tam, set<int> &cua) {
  int tamo = cua.size();
  int *vs, j=0;
  vs = new int[tamo];
  set<int>::iterator it;
  for (it = cua.begin(); it != cua.end(); it++) {
    vs[j] = vec[*it];
    j++;
  }
  sort(vs, nvec, tamo);
}

void Basics::sort(double* vec, double* nvec, int tam, set<int> &cua) {
  int tamo = cua.size();
  double *vs;
  int j=0;
  vs = new double[tamo];
  set<int>::iterator it;
  for (it = cua.begin(); it != cua.end(); it++) {
    vs[j] = vec[*it];
    j++;
  }
  sort(vs, nvec, tamo);
}


void Basics::sort(double* vec, double* nvec, int tam)
{
	double min, max;
	int i, j, quedan, k;
	double* wov;
	wov = new double[tam];
	double* wov2;
	for (i=0; i < tam; i++)
		wov[i] = vec[i];
	max = find_max(wov, tam);
	quedan = tam;
	i = 0;
	while (quedan > 0) {
		min = find_min(wov, quedan);
		for (j=0; j< quedan; j++) {
			if (wov[j] == min) {
				nvec[i] = wov[j];
				i++;
				break;
			}
		}
		wov2 = new double[quedan-1];
		for (k = 0; k < j; k++)
			wov2[k] = wov[k];
		for (k=j; k < (quedan-1); k++)
			wov2[k] = wov[k+1];
		delete [] wov;
		quedan--;
		wov = new double[quedan];
		for (k=0; k < quedan; k++)
			wov[k] = wov2[k];
		delete [] wov2;
	}
	if (min != max) {
		cout << "[Error]: sort does not converge in Basics::sort.\n";
		exit(1);
	}
	delete [] wov;
	return;
}

double Basics::get_median(int* vec, int tam)
{
	double med;
	int *vt;
	vt = new int[tam];
	sort(vec, vt, tam);
	if ((tam%2) == 0)
		med = (vt[tam/2] + vt[(tam/2)-1])/2.0;
	else
		med = vt[(tam-1)/2];
	delete [] vt;
	return med;
}

double Basics::get_median(double* vec, int tam)
{
	double med;
	double *vt;
	vt = new double[tam];
	sort(vec, vt, tam);
	if ((tam%2) == 0)
		med = (vt[tam/2] + vt[(tam/2)-1])/2.0;
	else
		med = vt[(tam-1)/2];
	delete [] vt;
	return med;
}

double Basics::get_median(int* vec, int tam, set<int> &cua) {
  int *vs;
  vs = new int[cua.size()];
  int j = 0;
  set<int>::iterator it;
  for (it = cua.begin(); it != cua.end(); it++) {
    vs[j] = vec[*it];
    j++;
  }
  double res = get_median(vs, cua.size());
  delete [] vs;
  return res;
}

double Basics::get_median(double* vec, int tam, set<int> &cua) {
  double *vs;
  vs = new double[cua.size()];
  int j = 0;
  set<int>::iterator it;
  for (it = cua.begin(); it != cua.end(); it++) {
    vs[j] = vec[*it];
    j++;
  }
  double res = get_median(vs, cua.size());
  delete [] vs;
  return res;
}

double Basics::get_midpoint(int* vo, int tam)//vector debe estar ordenado
{
  double res;
  if ((tam %2)==0)
    res = (vo[tam/2] + vo[(tam/2)-1])/2.0;
  else
    res = double(vo[(tam-1)/2]);
  return res;
}

double Basics::get_midpoint(double* vo, int tam)//vector debe estar ordenado
{
  double res;
  if ((tam %2)==0)
    res = (vo[tam/2] + vo[(tam/2)-1])/2.0;
  else
    res = vo[(tam-1)/2];
  return res;
}

double Basics::get_q1(int* vo, int tam)
{
  double res;
  if ((tam%2)==0)
    res = get_midpoint(vo, (tam/2));
  else
    res = get_midpoint(vo, ((tam-1)/2));
  return res;
}

double Basics::get_q1(double* vo, int tam)
{
  double res;
  if ((tam%2)==0)
    res = get_midpoint(vo, (tam/2));//no incluye mediana
  else
    res = get_midpoint(vo, ((tam-1)/2));
  return res;
}

double Basics::get_q3(int* vo, int tam)
{
  double res;
  if ((tam%2)==0) {
    if ((tam%4)==0)
      res = (vo[((tam*3)/4)] + vo[((tam*3)/4)-1])/2.0;
    else
      res = vo[((3*tam)-2)/4];
  }
  else {
    if ((((tam-1)/2)%2)==0)
      res = (vo[(tam-1)*(3/4)] + vo[((tam-1)*(3/4))+1])/2.0;
    else
      res = vo[((3*tam)-1)/4];
  }
  return res;
}

double Basics::get_q3(double* vo, int tam)
{
  double res;
  if ((tam%2)==0) {
    if ((tam%4)==0)
      res = (vo[((tam*3)/4)] + vo[((tam*3)/4)-1])/2.0;
    else
      res = vo[((3*tam)-2)/4];
  }
  else {
    if ((((tam-1)/2)%2)==0)
      res = (vo[(tam-1)*(3/4)] + vo[((tam-1)*(3/4))+1])/2.0;
    else
      res = vo[((3*tam)-1)/4];
  }
  return res;
}


double Basics::histo2(int nudata, int nubins, int* data1, int* data2, double* dis1, double* dis2, double* xax1, double* xax2)
{
	int min, max, min1, min2, max1, max2;
	min1 = find_min(data1, nudata);
	max1 = find_max(data1, nudata);
	min2 = find_min(data2, nudata);
	max2 = find_max(data2, nudata);
	if(min1 < min2)
		min = min1;
	else
		min = min2;
	if (max1 > max2)
		max = max1;
	else
		max = max2;
	double binsize = (max-min)/(nubins*1.0);
	int i, num;
	fillv0(dis1, nubins);
	fillv0(dis2, nubins);
	for (i = 0; i < nudata; i++){
		num = floor((data1[i]-min)/binsize);
		if (num >= nubins)
			num = num-1;
		dis1[num] = dis1[num]+1;
		num = floor((data2[i]-min)/binsize);
		if (num >= nubins)
			num = num-1;
		dis2[num] = dis2[num]+1;
	}
	for (i=0; i < nubins; i++) {
		dis1[i] = dis1[i]/(nudata*1.0);
		dis2[i] = dis2[i]/(nudata*1.0);
		xax1[i] = (i*binsize)+min-(binsize*0.2)+(binsize/2.0);
		xax2[i] = (i*binsize)+min+(binsize*0.2)+(binsize/2.0);
	}
	return binsize*0.4;
}

double Basics::histo(int nudata, int nubins, int* data, double* dis, double* xax, int min, int max)
{
	int i, num;
	double binsize = (max-min)/(nubins*1.0);
	fillv0(dis, nubins);
	for (i=0; i<nudata; i++) {
		num = floor((data[i]-min)/binsize);
		if (num >= nubins)
			num = num-1;
    dis[num] = dis[num]+1;
	}
	for (i=0; i<nubins; i++) {
		dis[i] = dis[i]/(nudata*1.0);
		xax[i] = (i*binsize)+min+(binsize/2.0);
	}
	return binsize;
}


int Basics::histo_int(int nudata, int nubins, int* data, double* dis, double* xax, int min, int max)
{
	int i, num;
	int binsize = 1;
	fillv0(dis, nubins);
	for (i=0; i<nudata; i++) {
		num = floor((data[i]-min)/binsize);
		if (num >= nubins)
			num = num-1;
		dis[num] = dis[num]+1;
	}
	for (i=0; i<nubins; i++) {
		dis[i] = dis[i]/(nudata*1.0);
		xax[i] = (i*binsize)+min;//+(binsize/2.0);
	}
	return binsize;
}

double Basics::histo(int nudata, int nubins, double* data, double* dis, double* xax, double min, double max)
{
	int i, num;
	double binsize = (max-min)/(nubins*1.0);
	fillv0(dis, nubins);
	for (i=0; i<nudata; i++) {
		num = floor((data[i]-min)/binsize);
		if (num >= nubins)
			num = num-1;
    dis[num] = dis[num]+1;
  }
	for (i=0; i<nubins; i++) {
		dis[i] = dis[i]/(nudata*1.0);
		xax[i] = (i*binsize)+min+(binsize/2.0);
	}
	return binsize;
}

double Basics::histo(int nudata, int nubins, int* data, double* dis, double* xax)
{
	int min, max;
	min = find_min(data, nudata);
	max = find_max(data, nudata);
	double binsize = histo(nudata, nubins, data, dis, xax, min, max);
	return binsize;
}

double Basics::histo(int nudata, int nubins, double* data, double* dis, double* xax)
{
	double min, max;
	min = find_min(data, nudata);
	max = find_max(data, nudata);
	double binsize =  histo(nudata, nubins, data, dis, xax, min, max);
	return binsize;
}

int Basics::histo_int(int nudata, int nubins, int* data, double* dis, double* xax)
{
	int min, max;
	min = find_min(data, nudata);
	max = find_max(data, nudata);
	int bs = histo_int(nudata, nubins, data, dis, xax, min, max);
	return bs;
}

double Basics::histo2(int nudata, int nubins, double* data1, double* data2, double* dis1, double* dis2, double* xax1, double* xax2)
{
	double min, max, min1, min2, max1, max2;
	min1 = find_min(data1, nudata);
	max1 = find_max(data1, nudata);
	min2 = find_min(data2, nudata);
	max2 = find_max(data2, nudata);
	if(min1 < min2)
		min = min1;
	else
		min = min2;
	if (max1 > max2)
		max = max1;
	else
		max = max2;
	double binsize = (max-min)/(nubins*1.0);
	int i, num;
	fillv0(dis1, nubins);
	fillv0(dis2, nubins);
	for (i = 0; i < nudata; i++){
		num = floor((data1[i]-min)/binsize);
		if (num >= nubins)
			num = num-1;
		dis1[num] = dis1[num]+1;
		num = floor((data2[i]-min)/binsize);
		if (num >= nubins)
			num = num-1;
		dis2[num] = dis2[num]+1;
	}
	for (i=0; i < nubins; i++) {
		dis1[i] = dis1[i]/(nudata*1.0);
		dis2[i] = dis2[i]/(nudata*1.0);
    xax1[i] = (i*binsize)+min-(binsize*0.2)+(binsize/2.0);
		xax2[i] = (i*binsize)+min+(binsize*0.2)+(binsize/2.0);
	}
	return binsize*0.4;
}


double Basics::histo2(int nudata1, int nudata2, int nubins, double* data1, double* data2, double* dis1, double* dis2, double* xax1, double* xax2)
{
  double min, max, min1, min2, max1, max2;
  min1 = find_min(data1, nudata1);
  max1 = find_max(data1, nudata1);
  min2 = find_min(data2, nudata2);
  max2 = find_max(data2, nudata2);
  if(min1 < min2)
    min = min1;
  else
    min = min2;
  if (max1 > max2)
    max = max1;
  else
    max = max2;
  double binsize = (max-min)/(nubins*1.0);
  int i, num;
  fillv0(dis1, nubins);
  fillv0(dis2, nubins);
  for (i = 0; i < nudata1; i++){
    num = floor((data1[i]-min)/binsize);
    if (num >= nubins)
      num = num-1;
    dis1[num] = dis1[num]+1;
  }
  for (i = 0; i < nudata2; i++){
    num = floor((data2[i]-min)/binsize);
    if (num >= nubins)
      num = num-1;
    dis2[num] = dis2[num]+1;
  }
  for (i=0; i < nubins; i++) {
    dis1[i] = dis1[i]/(nudata1*1.0);
    dis2[i] = dis2[i]/(nudata2*1.0);
    xax1[i] = (i*binsize)+min-(binsize*0.22);//+(binsize/2.0);
    xax2[i] = (i*binsize)+min+(binsize*0.22);//+(binsize/2.0);
  }
  return binsize*0.44;
}


int Basics::fdr(double* vps, int n, int* orden, double alpha, bool* cualessi)
{
	int i,j,res = 0;
	fillv0(orden, n);
	fillv0(cualessi, n);
	sort_vector(vps, orden, n,  0);
	j=-1;
	double y;
	for (i=(n-1); i>=0; i--) {
		y = i+1.0;
		if (!(vps[orden[i]] > ((y/n)*alpha))) {
			j = i;
			break;
		}
	}
	for (i=0; i<=j; i++)
		cualessi[orden[i]] = true;
	res = count_in_vector(cualessi, n, true);
	return res;
}

int Basics::fdr_dep(double* vps, int n, int* orden, double alpha, bool* cualessi)
{
	int i,j,res = 0;
	fillv0(orden, n);
	fillv0(cualessi, n);
	j=-1;
	sort_vector(vps, orden, n,  0);
	double y;
	for (i=(n-1); i>=0; i--) {
		y = i+1.0;
		if (!(vps[orden[i]] > ((y/n)*(alpha/sum_consec_inverse(1, (i+1)))))) {
			j = i;
			break;
		}
	}
	for (i=0; i<=j; i++)
		cualessi[orden[i]] = true;
	res = count_in_vector(cualessi, n, true);
	return res;
}

double Basics::sum_consec_inverse(int from, int to)
{
	double res = 0;
	if (from > to) {
		cout << "[Error]: The first argument in Basics::sum_consec_inverse cannot be greater than the second argument.\n";
		exit(1);
	}
	int i;
	for (i=from; i<=to; i++)
		res = res+(1.0/i);
	return res;
}

double Basics::from_Z_to_p(double zscore, char gol)
{
	double res;
	if ((gol!='g') && (gol!='G') &&(gol!='l')&&(gol!='L')) {
		cout << "[Error]: Undefined Ha in Basics::from_Z_to_p . \n";
		exit(1);
	}
	if ((gol=='g')||(gol=='G'))
		res = 1 - gsl_cdf_ugaussian_P(zscore);
	else
		res = gsl_cdf_ugaussian_P(zscore);
	return res;
}

double Basics::from_Z_to_p(double zscore)
{
	double res = 2*(gsl_cdf_ugaussian_P((-1)*(fabs(zscore))));
	return res;
}

double Basics::from_t_to_p(double t, int df, char gol) //one-tailed
{
  double res;
  if ((gol!='g') && (gol!='G') &&(gol!='l')&&(gol!='L')) {
    cout << "[Error]: Undefined Ha in Basics::from_t_to_p . \n";
    exit(1);
  }
  if ((gol=='g')||(gol=='G'))
    res = gsl_cdf_tdist_Q(t, df);
  else
    res = gsl_cdf_tdist_P(t, df);
  return res;
}

double Basics::from_t_to_p(double t, int df) //two-tailed
{
  double res = 2*(gsl_cdf_tdist_Q(fabs(t), df));
  return res;
}


double Basics::variance_binomial(double p, int n)
{
	double res = (p*(1.0-p))/n;
	return res;
}

double Basics::stddev_binomial(double p, int n)
{
	double res = sqrt(variance_binomial(p, n));
	return res;
}

double Basics::zscore_binary(double obsfreq, double expfreq, int n)
{
	double z = (obsfreq - expfreq)/(stddev_binomial(expfreq, n));
	return z;
}

double Basics::pearsons_chi(double& chi, int& df, double **obs, double **exp, int rows, int cols)
{
	double pv=0;
	int i, j;
	chi = 0;
	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			chi = chi + (((obs[i][j]-exp[i][j])*(obs[i][j]-exp[i][j]))/exp[i][j]);
		}
	}
	df = (rows-1)*(cols-1);
	pv = gsl_cdf_chisq_Q(chi, df);
	
	return pv;
}

double Basics::pearsons_chi(double& chi, int& df, double **condmat, int rows, int cols)
{
	double pv;
	int i,j;
	double **exp, *sumcol, *sumrow, tsum = 0;
	sumcol = new double[cols];
	sumrow = new double[rows];
	fillv0(sumcol, cols);
	fillv0(sumrow, rows);
	for (i=0; i<cols; i++)
		sumcol[i] = sum_column(condmat, rows, cols, i);
	for (i=0; i <rows; i++)
		sumrow[i] = sumatoria(condmat[i], cols);
	tsum = sumatoria(sumcol, cols);
	if (tsum != sumatoria(sumrow, rows)) {
		cout << "[Error]: Columns and rows do not match in Basics::pearsons_chi .\n";
		exit(1);
	}
	exp = new double*[rows];
	for (i=0; i<rows; i++)
		exp[i] = new double[cols];
	for (i=0; i< rows; i++)
		for (j=0; j<cols; j++)
			exp[i][j] = (sumrow[i]*sumcol[j])/tsum;
	pv = pearsons_chi(chi, df, condmat, exp, rows, cols);
	for (i=0; i<rows; i++)
		delete [] exp[i];
	delete [] exp;
	delete [] sumcol;
	delete [] sumrow;
	return pv;
}


/////
double Basics::yates_chi(double& chi, int& df, double **obs, double **exp, int rows, int cols)
{
	double pv=0;
	int i, j;
	chi = 0;
	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			chi = chi + (((fabs(obs[i][j]-exp[i][j])-0.5)*(fabs(obs[i][j]-exp[i][j])-0.5))/exp[i][j]);
		}
	}
	df = (rows-1)*(cols-1);
	pv = gsl_cdf_chisq_Q(chi, df);
	
	return pv;
}

double Basics::yates_chi(double& chi, int& df, double **condmat, int rows, int cols)
{
	double pv;
	int i,j;
	double **exp, *sumcol, *sumrow, tsum = 0;
	sumcol = new double[cols];
	sumrow = new double[rows];
	fillv0(sumcol, cols);
	fillv0(sumrow, rows);
	for (i=0; i<cols; i++)
		sumcol[i] = sum_column(condmat, rows, cols, i);
	for (i=0; i <rows; i++)
		sumrow[i] = sumatoria(condmat[i], cols);
	tsum = sumatoria(sumcol, cols);
	if (tsum != sumatoria(sumrow, rows)) {
		cout << "[Error]: Columns and rows do not match in Basics::pearsons_chi .\n";
		exit(1);
	}
	exp = new double*[rows];
	for (i=0; i<rows; i++)
		exp[i] = new double[cols];
	for (i=0; i< rows; i++)
		for (j=0; j<cols; j++)
			exp[i][j] = (sumrow[i]*sumcol[j])/tsum;
	pv = yates_chi(chi, df, condmat, exp, rows, cols);
	for (i=0; i<rows; i++)
		delete [] exp[i];
	delete [] exp;
	delete [] sumcol;
	delete [] sumrow;
	return pv;
}

//aux for stats
int Basics::get_nubins(int *vector, int tvec)
{
	int nubins;
	if (tvec>0)
		nubins = find_max(vector, tvec) - find_min(vector, tvec) + 1;
	else
		nubins = 0;
	return nubins;
}

double Basics::pearsons_r(double& r, double* veh, double* vev, int tam, char tgol, double& S, int& df)
{
  ofstream fs;
  ifstream fe;
  int i;
  string alti;
  if ((tgol == 't') || (tgol =='T'))
    alti = "two.sided";
  else {
    if ((tgol =='G') ||(tgol=='g'))
      alti = "greater";
    else
      alti = "less";
  }
  open_ofstream(fs, "p0o9i8u7Rscript.sh");
  fs << "x <- c(" << veh[0];
  for (i = 1; i < tam; i++)
    fs << ", " << veh[i];
  fs << ")\n";
  fs << "y <- c(" << vev[0];
  for (i = 1; i < tam; i++)
    fs << ", " << vev[i];
  fs << ")\n";
  fs << "pru <- cor.test(x,y, alternative = \"" << alti << "\", method= \"pearson\", exact = TRUE)\n";
  fs << "capture.output(pru$estimate, file=\"p0o9i8u7stats.txt\")\n";
  fs << "capture.output(pru$p.value, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$statistic, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$parameter, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "q()\n";
  fs.close();
  system("Rscript p0o9i8u7Rscript.sh");
  double p;
  string gar;
  open_ifstream(fe, "p0o9i8u7stats.txt");
  fe >> gar;
  fe >> r;
  fe >> gar;
  fe >> p;
  fe >> gar;
  fe >> S;
  fe >> gar;
  fe >> df;
  fe.close();
  system("rm p0o9i8u7Rscript.sh p0o9i8u7stats.txt");
  return p;
}

double Basics::kolmogorov_smirnov(double& D, double* v1, double* v2, int tam1, int tam2) {
  ofstream fs;
  ifstream fe;
  int i;
  open_ofstream(fs, "p0o9i8u7Rscript.sh");
  fs << "x <- c(" << v1[0];
  for (i = 1; i < tam1; i++)
    fs << ", " << v1[i];
  fs << ")\n";
  fs << "y <- c(" << v2[0];
  for (i = 1; i < tam2; i++)
    fs << ", " << v2[i];
  fs << ")\n";
  fs << "pru <- ks.test(x,y)\n";
  fs << "capture.output(pru$statistic, file=\"p0o9i8u7stats.txt\")\n";
  fs << "capture.output(pru$p.value, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "q()\n";
  fs.close();
  system("Rscript p0o9i8u7Rscript.sh");
  double p;
  string gar;
  open_ifstream(fe, "p0o9i8u7stats.txt");
  fe >> gar;
  fe >> D;
  fe >> gar;
  fe >> p;
  fe.close();
  system("rm p0o9i8u7Rscript.sh p0o9i8u7stats.txt");
  return p;
}

double Basics::kolmogorov_smirnov(double& D, int* v1, int* v2, int tam1, int tam2) {
  ofstream fs;
  ifstream fe;
  int i;
  open_ofstream(fs, "p0o9i8u7Rscript.sh");
  fs << "x <- c(" << v1[0];
  for (i = 1; i < tam1; i++)
    fs << ", " << v1[i];
  fs << ")\n";
  fs << "y <- c(" << v2[0];
  for (i = 1; i < tam2; i++)
    fs << ", " << v2[i];
  fs << ")\n";
  fs << "pru <- ks.test(x,y)\n";
  fs << "capture.output(pru$statistic, file=\"p0o9i8u7stats.txt\")\n";
  fs << "capture.output(pru$p.value, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "q()\n";
  fs.close();
  system("Rscript p0o9i8u7Rscript.sh");
  double p;
  string gar;
  open_ifstream(fe, "p0o9i8u7stats.txt");
  fe >> gar;
  fe >> D;
  fe >> gar;
  fe >> p;
  fe.close();
  system("rm p0o9i8u7Rscript.sh p0o9i8u7stats.txt");
  return p;
}

double Basics::pearsons_r(double& r, double* veh, int* vev, int tam, char tgol, double& S, int& df) {
  ofstream fs;
  ifstream fe;
  int i;
  string alti;
  if ((tgol == 't') || (tgol =='T'))
    alti = "two.sided";
  else {
    if ((tgol =='G') ||(tgol=='g'))
      alti = "greater";
    else
      alti = "less";
  }
  open_ofstream(fs, "p0o9i8u7Rscript.sh");
  fs << "x <- c(" << veh[0];
  for (i = 1; i < tam; i++)
    fs << ", " << veh[i];
  fs << ")\n";
  fs << "y <- c(" << vev[0];
  for (i = 1; i < tam; i++)
    fs << ", " << vev[i];
  fs << ")\n";
  fs << "pru <- cor.test(x,y, alternative = \"" << alti << "\", method= \"pearson\", exact = TRUE)\n";
  fs << "capture.output(pru$estimate, file=\"p0o9i8u7stats.txt\")\n";
  fs << "capture.output(pru$p.value, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$statistic, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$parameter, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "q()\n";
  fs.close();
  system("Rscript p0o9i8u7Rscript.sh");
  double p;
  string gar;
  open_ifstream(fe, "p0o9i8u7stats.txt");
  fe >> gar;
  fe >> r;
  fe >> gar;
  fe >> p;
  fe >> gar;
  fe >> S;
  fe >> gar;
  fe >> df;
  fe.close();
  system("rm p0o9i8u7Rscript.sh p0o9i8u7stats.txt");
  return p;
}

double Basics::pearsons_r(double& r, int* veh, int* vev, int tam, char tgol, double& S, int& df)
{
  ofstream fs;
  ifstream fe;
  int i;
  string alti;
  if ((tgol == 't') || (tgol =='T'))
    alti = "two.sided";
  else {
    if ((tgol =='G') ||(tgol=='g'))
      alti = "greater";
    else
      alti = "less";
  }
  open_ofstream(fs, "p0o9i8u7Rscript.sh");
  fs << "x <- c(" << veh[0];
  for (i = 1; i < tam; i++)
    fs << ", " << veh[i];
  fs << ")\n";
  fs << "y <- c(" << vev[0];
  for (i = 1; i < tam; i++)
    fs << ", " << vev[i];
  fs << ")\n";
  fs << "pru <- cor.test(x,y, alternative = \"" << alti << "\", method= \"pearson\", exact = TRUE)\n";
  fs << "capture.output(pru$estimate, file=\"p0o9i8u7stats.txt\")\n";
  fs << "capture.output(pru$p.value, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$statistic, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$parameter, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "q()\n";
  fs.close();
  system("Rscript p0o9i8u7Rscript.sh");
  double p;
  string gar;
  open_ifstream(fe, "p0o9i8u7stats.txt");
  fe >> gar;
  fe >> r;
  fe >> gar;
  fe >> p;
  fe >> gar;
  fe >> S;
  fe >> gar;
  fe >> df;
  fe.close();
  system("rm p0o9i8u7Rscript.sh p0o9i8u7stats.txt");
  return p;
}

double Basics::spearmans_rho(double& r, double* veh, double* vev, int tam, char tgol, double& S, int& df)
{
  ofstream fs;
  ifstream fe;
  int i;
  string alti;
  if ((tgol == 't') || (tgol =='T'))
    alti = "two.sided";
  else {
    if ((tgol =='G') ||(tgol=='g'))
      alti = "greater";
    else
      alti = "less";
  }
  open_ofstream(fs, "p0o9i8u7Rscript.sh");
  fs << "x <- c(" << veh[0];
  for (i = 1; i < tam; i++)
    fs << ", " << veh[i];
  fs << ")\n";
  fs << "y <- c(" << vev[0];
  for (i = 1; i < tam; i++)
    fs << ", " << vev[i];
  fs << ")\n";
  fs << "pru <- cor.test(x,y, alternative = \"" << alti << "\", method= \"spearman\", exact = TRUE)\n";
  fs << "capture.output(pru$estimate, file=\"p0o9i8u7stats.txt\")\n";
  fs << "capture.output(pru$p.value, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$statistic, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$parameter, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "q()\n";
  fs.close();
  system("Rscript p0o9i8u7Rscript.sh");
  double p;
  string gar;
  open_ifstream(fe, "p0o9i8u7stats.txt");
  fe >> gar;
  fe >> r;
  fe >> gar;
  fe >> p;
  fe >> gar;
  fe >> S;
  fe >> gar;
  fe >> df;
  fe.close();
  system("rm p0o9i8u7Rscript.sh p0o9i8u7stats.txt");
  return p;
}

double Basics::spearmans_rho(double& r, double* veh, int* vev, int tam, char tgol, double& S, int& df)
{
  ofstream fs;
  ifstream fe;
  int i;
  string alti;
  if ((tgol == 't') || (tgol =='T'))
    alti = "two.sided";
  else {
    if ((tgol =='G') ||(tgol=='g'))
      alti = "greater";
    else
      alti = "less";
  }
  open_ofstream(fs, "p0o9i8u7Rscript.sh");
  fs << "x <- c(" << veh[0];
  for (i = 1; i < tam; i++)
    fs << ", " << veh[i];
  fs << ")\n";
  fs << "y <- c(" << vev[0];
  for (i = 1; i < tam; i++)
    fs << ", " << vev[i];
  fs << ")\n";
  fs << "pru <- cor.test(x,y, alternative = \"" << alti << "\", method= \"spearman\", exact = TRUE)\n";
  fs << "capture.output(pru$estimate, file=\"p0o9i8u7stats.txt\")\n";
  fs << "capture.output(pru$p.value, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$statistic, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$parameter, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "q()\n";
  fs.close();
  system("Rscript p0o9i8u7Rscript.sh");
  double p;
  string gar;
  open_ifstream(fe, "p0o9i8u7stats.txt");
  fe >> gar;
  fe >> r;
  fe >> gar;
  fe >> p;
  fe >> gar;
  fe >> S;
  fe >> gar;
  fe >> df;
  fe.close();
  system("rm p0o9i8u7Rscript.sh p0o9i8u7stats.txt");
  return p;
}

double Basics::spearmans_rho(double& r, int* veh, int* vev, int tam, char tgol, double& S, int& df)
{
  ofstream fs;
  ifstream fe;
  int i;
  string alti;
  if ((tgol == 't') || (tgol =='T'))
    alti = "two.sided";
  else {
    if ((tgol =='G') ||(tgol=='g'))
      alti = "greater";
    else
      alti = "less";
  }
  open_ofstream(fs, "p0o9i8u7Rscript.sh");
  fs << "x <- c(" << veh[0];
  for (i = 1; i < tam; i++)
    fs << ", " << veh[i];
  fs << ")\n";
  fs << "y <- c(" << vev[0];
  for (i = 1; i < tam; i++)
    fs << ", " << vev[i];
  fs << ")\n";
  fs << "pru <- cor.test(x,y, alternative = \"" << alti << "\", method= \"spearman\", exact = TRUE)\n";
  fs << "capture.output(pru$estimate, file=\"p0o9i8u7stats.txt\")\n";
  fs << "capture.output(pru$p.value, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$statistic, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "capture.output(pru$parameter, file=\"p0o9i8u7stats.txt\", append = TRUE)\n";
  fs << "q()\n";
  fs.close();
  system("Rscript p0o9i8u7Rscript.sh");
  double p;
  string gar;
  open_ifstream(fe, "p0o9i8u7stats.txt");
  fe >> gar;
  fe >> r;
  fe >> gar;
  fe >> p;
  fe >> gar;
  fe >> S;
  fe >> gar;
  fe >> df;
  fe.close();
  system("rm p0o9i8u7Rscript.sh p0o9i8u7stats.txt");
  return p;
}

double Basics::shapiro_wilk(double& W, double* vec, int tam)
{
	ofstream fs;
	ifstream fe;
	int i;
	open_ofstream(fs, "q1w2e3r4data.txt");
	for (i=0; i< tam; i++)
		fs << vec[i] << endl;
	fs.close();
	open_ofstream(fs, "q1w2e3r4Rscript.sh");
	fs << "elve <- read.csv(file=\"q1w2e3r4data.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
	if (tam <= 5000)
		fs << "pru <- shapiro.test(elve$V1)\n";
	else
		fs << "pru <- shapiro.test(sample(elve$V1, 5000, replace=TRUE))\n";
	fs << "capture.output(pru$p.value, file=\"q1w2e3r4p-val.txt\")\n";
	fs << "capture.output(pru$statistic, file=\"q1w2e3r4stat.txt\")\n";
	fs << "q()\n";
	fs.close();
	double p;
	string gar;
	system("Rscript q1w2e3r4Rscript.sh");
	open_ifstream(fe, "q1w2e3r4stat.txt");
	fe >> gar;
	fe >> W;
	fe.close();
	open_ifstream(fe, "q1w2e3r4p-val.txt");
	fe >> gar;
	fe >> p;
	fe.close();
	system("rm q1w2e3r4data.txt q1w2e3r4Rscript.sh q1w2e3r4p-val.txt q1w2e3r4stat.txt");
	return p;
}

double Basics::shapiro_wilk(double& W, int* vec, int tam)
{
	ofstream fs;
	ifstream fe;
	int i;
	open_ofstream(fs, "q1w2e3r4data.txt");
	for (i=0; i< tam; i++)
		fs << vec[i] << endl;
	fs.close();
	open_ofstream(fs, "q1w2e3r4Rscript.sh");
	fs << "elve <- read.csv(file=\"q1w2e3r4data.txt\", sep=\"\\t\", head=FALSE, colClasses=\"numeric\")\n";
	if (tam <= 5000)
		fs << "pru <- shapiro.test(elve$V1)\n";
	else
		fs << "pru <- shapiro.test(sample(elve$V1, 5000, replace=TRUE))\n";
	fs << "capture.output(pru$p.value, file=\"q1w2e3r4p-val.txt\")\n";
	fs << "capture.output(pru$statistic, file=\"q1w2e3r4stat.txt\")\n";
	fs << "q()\n";
	fs.close();
	double p;
	string gar;
	system("Rscript q1w2e3r4Rscript.sh");
	open_ifstream(fe, "q1w2e3r4stat.txt");
	fe >> gar;
	fe >> W;
	fe.close();
	open_ifstream(fe, "q1w2e3r4p-val.txt");
	fe >> gar;
	fe >> p;
	fe.close();
	system("rm q1w2e3r4data.txt q1w2e3r4Rscript.sh q1w2e3r4p-val.txt q1w2e3r4stat.txt");
	return p;
}

double Basics::brown_forsythe(double& F, int& dfn, int& dfd, double **groupsonrows, int *tams, int numbgroups)
{
	int i,j, totnumbofobs=0;
	double *medians, **zeta, *meansz, themeanz, sumz=0, numerator=0, denominator=0, p;
	medians = new double[numbgroups];
	meansz = new double[numbgroups];
	zeta = new double*[numbgroups];
	for (j=0; j<numbgroups; j++)
		zeta[j] = new double[tams[j]];
	for (j=0; j < numbgroups; j++) {
		totnumbofobs = totnumbofobs + tams[j];
		medians[j] = get_median(groupsonrows[j], tams[j]);
		for (i=0; i<tams[j]; i++) {
			zeta[j][i] = fabs(groupsonrows[j][i] - medians[j]);
			sumz = sumz + zeta[j][i];
		}
		meansz[j] = get_mean(zeta[j], tams[j]);
	}
	themeanz = sumz/totnumbofobs;
	
	for (j=0; j<numbgroups; j++) {
		numerator = numerator+ ((tams[j]*((meansz[j] - themeanz)*(meansz[j] - themeanz)))/(numbgroups-1.0));
		for (i=0; i < tams[j]; i++)
			denominator = denominator + (((zeta[j][i] - meansz[j])*(zeta[j][i] - meansz[j]))/(totnumbofobs-numbgroups));
	}
	F = numerator/denominator;
	dfn = numbgroups - 1;
	dfd = totnumbofobs - numbgroups;
	
	delete [] medians;
	delete [] meansz;
	for (j=0; j<numbgroups; j++)
		delete [] zeta[j];
	delete [] zeta;
	p = gsl_cdf_fdist_Q(F, dfn, dfd);
	return p;
}

double Basics::levene(double& F, int& dfn, int& dfd, double **groupsonrows, int *tams, int numbgroups)
{
  int i,j, totnumbofobs=0;
  double *means, **zeta, *meansz, themeanz, sumz=0, numerator=0, denominator=0, p;
  means = new double[numbgroups];
  meansz = new double[numbgroups];
  zeta = new double*[numbgroups];
  for (j=0; j<numbgroups; j++)
    zeta[j] = new double[tams[j]];
  for (j=0; j < numbgroups; j++) {
    totnumbofobs = totnumbofobs + tams[j];
    means[j] = get_mean(groupsonrows[j], tams[j]);
    for (i=0; i<tams[j]; i++) {
      zeta[j][i] = fabs(groupsonrows[j][i] - means[j]);
      sumz = sumz + zeta[j][i];
    }
    meansz[j] = get_mean(zeta[j], tams[j]);
  }
  themeanz = sumz/totnumbofobs;
  
  for (j=0; j<numbgroups; j++) {
    numerator = numerator+ ((tams[j]*((meansz[j] - themeanz)*(meansz[j] - themeanz)))/(numbgroups-1.0));
    for (i=0; i < tams[j]; i++)
      denominator = denominator + (((zeta[j][i] - meansz[j])*(zeta[j][i] - meansz[j]))/(totnumbofobs-numbgroups));
  }
  F = numerator/denominator;
  dfn = numbgroups - 1;
  dfd = totnumbofobs - numbgroups;
  
  delete [] means;
  delete [] meansz;
  for (j=0; j<numbgroups; j++)
    delete [] zeta[j];
  delete [] zeta;
  p = gsl_cdf_fdist_Q(F, dfn, dfd);
  return p;
}
double Basics::students_t(double& T, int& df, double *vec1, int tam1, double *vec2, int tam2)
{
	double var1, var2, pstd, sted, p;
	var1 = get_sample_variance(vec1, tam1);
	var2 = get_sample_variance(vec2, tam2);
	pstd = (((tam1-1)*var1)+((tam2-1)*var2))/(tam1 +tam2 - 2);
	pstd = sqrt(pstd);	//pooled std dev
	sted = pstd*(sqrt((1.0/tam1)+(1.0/tam2))); //std error of the difference between 2 means
	T = (get_mean(vec1, tam1) - get_mean(vec2, tam2))/sted;
  if (T < 0)
    T = T*(-1);
	df = tam1 +tam2 -2;
	p = gsl_cdf_tdist_Q(T, df);
	return p;
}

double Basics::students_t(double& T, int& df, double *vec1, double *vec2, int tam)
{
	double p = students_t(T, df, vec1, tam, vec2, tam);
	return p;
}

double Basics::students_t(double& T, int& df, double *vec, int tam, double mu) {
  double var, sted, p, ave;
  ave = get_mean(vec, tam);
  var = get_sample_variance(vec, tam, ave);
  sted = sqrt(var/tam);
  T = (ave - mu)/sted;
  if (T < 0)
    T = T*(-1);
  df = tam-1;
  p = gsl_cdf_tdist_Q(T, df);
  return p;
}

double Basics::paired_students_t(double& T, int& df, double *vec1, double *vec2, int tam) {
  double *dife, p;
  dife = new double[tam];
  int i;
  for (i = 0; i < tam; i++)
    dife[i] = vec1[i] -vec2[i];
  p = students_t(T, df, dife, tam, 0);
  delete [] dife;
  return p;
}


double Basics::welch_t(double& T, int& df, double *vec1, int tam1, double *vec2, int tam2)
{
	double dfd;
	double mean1, mean2, sveom1, sveom2, p, denom;
	mean1 = get_mean(vec1, tam1);
	mean2 = get_mean(vec2, tam2);
	sveom1 = get_samp_var_err_of_mean(vec1, tam1);
	sveom2 = get_samp_var_err_of_mean(vec2, tam2);
	denom = sqrt(sveom1+sveom2);
	T = (mean1-mean2)/denom;
	dfd = ((sveom1+sveom2)*(sveom1+sveom2))/(((sveom1*sveom1)/(tam1-1))+((sveom2*sveom2)/(tam2-1)));
	df = round(dfd);
	p = gsl_cdf_tdist_Q(T, df);
	return p;
}

double Basics::mann_whitney(double *vec1, int tam1, double *vec2, int tam2)
{
	int tam = tam1+tam2, i, remin, nuti, j;
	double *vecall, *vecorall, *ranks, *preranks, mean, sum, sd, u;
	bool *belto1;
	int *orden;
	vecall = new double[tam];
	vecorall = new double[tam];
	ranks = new double[tam];
	preranks = new double[tam];
	belto1 = new bool[tam];
	fillv0(vecall, tam);
	fillv0(vecorall, tam);
	fillv0(ranks, tam);
	fillv0(belto1, tam);
	orden= new int[tam];
	fillv0(orden, tam);
	for (i=0; i < tam1; i++)
		vecorall[i] = vec1[i];
	for (i=tam1; i < tam; i++)
		vecorall[i] = vec2[i-tam1];
	remin = find_min(vecorall, tam) - 100;
	sort_vector(vecorall, orden, tam, remin);
	for (i=0; i<tam; i++) {
		vecall[i] = vecorall[orden[i]];
		if (orden[i] < tam1)
			belto1[i] = true;
		else
			belto1[i] = false;
		preranks[i] = i+1;
	}
	for (i=0; i<tam; i++) {
		nuti = 0;
		do {
			nuti++;
		} while (vecall[i+nuti]==vecall[i]);
		sum = 0;
		for (j=0; j<nuti; j++)
			sum = sum + preranks[i+j];
		for (j=0; j<nuti; j++)
			ranks[i+j] = double(sum)/double(nuti);
		i = i+j-1;
	}
	mean = (tam1*(tam+1))/2.0;
	sd = sqrt((tam1*tam2*(tam+1))/12.0);
	sum = 0;
	for (i=0; i<tam; i++)
		if (belto1[i]) {
			sum = sum + ranks[i];
		}
	u = (sum - mean)/sd;
	return u;
}


double Basics::mann_whitney(int *vec1, int tam1, int *vec2, int tam2)
{
	int tam = tam1+tam2, i, remin, nuti, j;
	double *vecall, *vecorall, *ranks, *preranks, mean, sum, sd, u;
	bool *belto1;
	int *orden;
	vecall = new double[tam];
	vecorall = new double[tam];
	ranks = new double[tam];
	preranks = new double[tam];
	belto1 = new bool[tam];
	fillv0(vecall, tam);
	fillv0(vecorall, tam);
	fillv0(ranks, tam);
	fillv0(belto1, tam);
	orden= new int[tam];
	fillv0(orden, tam);
	for (i=0; i < tam1; i++)
		vecorall[i] = vec1[i];
	for (i=tam1; i < tam; i++)
		vecorall[i] = vec2[i-tam1];
	remin = find_min(vecorall, tam) - 100;
	sort_vector(vecorall, orden, tam, remin);
	for (i=0; i<tam; i++) {
		vecall[i] = vecorall[orden[i]];
		if (orden[i] < tam1)
			belto1[i] = true;
		else
			belto1[i] = false;
		preranks[i] = i+1;
	}
	for (i=0; i<tam; i++) {
		nuti = 0;
		do {
			nuti++;
		} while (vecall[i+nuti]==vecall[i]);
		sum = 0;
		for (j=0; j<nuti; j++)
			sum = sum + preranks[i+j];
		for (j=0; j<nuti; j++)
			ranks[i+j] = double(sum)/double(nuti);
		i = i+j-1;
	}
	mean = (tam1*(tam+1))/2.0;
	sd = ((sqrt(tam1))*(sqrt(tam2))*(sqrt(tam+1)))/(sqrt(12.0));
	sum = 0;
	for (i=0; i<tam; i++)
		if (belto1[i]) {
			sum = sum + ranks[i];
		}
	u = (sum - mean)/sd;
	return u;
}

double Basics::zscore(double *vec, int tam, double mu)
{
	double z;
	z = (mu - get_mean(vec, tam))/ get_pop_stddev(vec, tam);
	return z;
}

double Basics::zscore(int *vec, int tam, int mu)
{
	double z;
	z = (double(mu - get_mean(vec, tam)))/(double(get_pop_stddev(vec, tam)));
	return z;
}

void Basics::transpose(int r, int c, double **mator, double **mattrans)
{
  int i,j;
  for (i=0; i < r; i++)
    for (j=0; j < c; j++)
      mattrans[j][i] = mator[i][j];
}

void Basics::append_idmat(int r, int c, double **mator, double **matfin)
{
  int i, j;
  for (i=0; i < r; i++)
    for (j=0; j<c; j++)
      matfin[i][j] = mator[i][j];
  for (i=r; i < (r+c); i++) {
    for (j=0; j < c; j++) {
      if (j==(i-r))
        matfin[i][j] = 1;
      else
        matfin[i][j] = 0;
    }
  }
}

void Basics::order_for_echelon(int r, int c, double **mator, double **matfin)
{
  int i,j,k, van = 0;
  bool *ya;
  ya = new bool[r];
  fillv0(ya, r);
  for (i=0; i < c; i++) {
    if (count_in_vector(ya, r, false) >0) {
      for (j=0; j<r; j++) {
        if (!ya[j]) {
          if (mator[j][i] != 0) {
            for (k=0; k < c; k++)
              matfin[van][k] = mator[j][k];
            ya[j] = true;
            van++;
          }
        }
      }
    }
    else
      break;
  }
  delete [] ya;
}

void Basics::eigenvv_sym(double **mat, int msize, double **evecs, double *evals){
  int i, j;
  gsl_matrix * mage = gsl_matrix_calloc (msize, msize);
  for (i = 0; i < msize; i++)
    for (j = 0; j < msize; j++)
      gsl_matrix_set(mage, i, j, mat[i][j]);
  
  gsl_vector *eval = gsl_vector_alloc (msize);
  gsl_matrix *evec = gsl_matrix_alloc (msize, msize);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (msize);
  int estoque = 0;
  estoque = gsl_eigen_symmv(mage, eval, evec, w);
  
  estoque = gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  
  for (i = 0; i < msize; i++) {
    evals[i] = gsl_vector_get(eval, i);
    for (j = 0; j < msize; j++)
      evecs[i][j] = gsl_matrix_get(evec, i, j);
  }
  gsl_eigen_symmv_free (w);
  gsl_matrix_free(mage);
  gsl_vector_free (eval);
  gsl_matrix_free(evec);
  i = estoque;
  return;
}

void Basics::eigenvv_sym(int **mat, int msize, double **evecs, double *evals){
  double **dmat;
  create_array(dmat, msize, msize);
  int i,j;
  for(i=0;i < msize; i++)
    for (j = 0; j < msize; j++)
      dmat[i][j] = double(mat[i][j]);
  eigenvv_sym(dmat, msize, evecs, evals);
  for (i = 0; i < msize; i++)
    delete [] dmat[i];
  delete [] dmat;
  return;
}

double Basics::get_leading_evector_symmat(double **mat, int msize, double *evector) {
  double eigenval;
  int i, j;
  gsl_matrix * mage = gsl_matrix_calloc (msize, msize);
  for (i = 0; i < msize; i++)
    for (j = 0; j < msize; j++)
      gsl_matrix_set(mage, i, j, mat[i][j]);
  gsl_vector *eval = gsl_vector_alloc (msize);
  gsl_matrix *evec = gsl_matrix_alloc (msize, msize);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (msize);
  int estoque = 0;
  estoque = gsl_eigen_symmv(mage, eval, evec, w);
  estoque = gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  eigenval = gsl_vector_get(eval, 0);
  for (i = 0; i < msize; i++)
    evector[i] = gsl_matrix_get(evec, i, 0);
  gsl_eigen_symmv_free (w);
  gsl_matrix_free(mage);
  gsl_vector_free (eval);
  gsl_matrix_free(evec);
  i=estoque;
  return eigenval;
}

double Basics::get_leading_evector_symmat(int **mat, int msize, double *evector) {
  double eval;
  double **dmat;
  create_array(dmat, msize, msize);
  int i,j;
  for(i=0;i < msize; i++)
    for (j = 0; j < msize; j++)
      dmat[i][j] = double(mat[i][j]);
  eval = get_leading_evector_symmat(dmat, msize, evector);
  for (i = 0; i < msize; i++)
    delete [] dmat[i];
  delete [] dmat;
  return eval;
}

double Basics::pow_meth(double **matA, int msize, double *x) {
  double evalmu = 0, err, tol = 0.000001;
  double *y, xp;
  y = new double[msize];
  double *verr;
  verr = new double[msize];
  int k = 0;
  int max = 3000;
  int p, i;
  for (i = 0; i < msize; i++)
    x[i] = est.randreal();
  p = where_is_vector_norm(x, msize);
  xp = x[p];
  for (i = 0; i < msize; i++)
    x[i] = x[i]/xp;
  while (k < max) {
    matxvec(x, msize, msize, msize, matA, y);
    evalmu = y[p];
    p = where_is_vector_norm(y, msize);
    if (y[p] == 0) {
      cout << "[Error]: Zero eigenvalue when calling Basics::pow_meth.\n";
      exit(1);
    }
    for (i = 0; i < msize; i++)
      verr[i] = fabs(x[i] - y[i]/y[p]);
    err = find_max(verr, msize);
    for (i = 0; i < msize; i++)
      x[i] = y[i]/y[p];
    if (err < tol)
      break;
    else
      k++;
    }
  if (k >= max) {
    cout << "[Error]: Maximum number of iterations when calling Basics::pow_meth.\n";
    exit(1);
    
  }
  return evalmu;
}

