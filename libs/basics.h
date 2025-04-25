#ifndef BASICS_H
#define BASICS_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <set>
#include <list>
#include <string>
#include "alea.h"
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

class Basics
{
public:
	Basics();
	Basics(Alea& jacta);
	void start_rng(Alea& jacta);
	
	//vectors
  
  bool all_entries_are_finite(double *vec, int tam);
  int number_of_finite_entries(double*vec, int tam);
  int where_is_vector_norm (int *vec, int tam); //posición de max en valor absoluto
  int where_is_vector_norm (double *vec, int tam);
  
	bool eqvec(int vec1[], int s1, int vec2[], int s2);
	bool eqvec(double vec1[], int s1, double vec2[], int s2);
	bool eqvec(bool vec1[], int s1, bool vec2[], int s2);
  bool eqvec(char vec1[], int s1, char vec2[], int s2);
  bool eqvec(string vec1[], int s1, string vec2[], int s2);
  
	bool same_elements_in_vecs(int vec1[], int s1, int vec2[], int s2);
	bool same_elements_in_vecs(double vec1[], int s1, double vec2[], int s2);
	bool same_elements_in_vecs(char vec1[], int s1, char vec2[], int s2);
	bool same_elements_in_vecs(string vec1[], int s1, string vec2[], int s2);
  
	int difsinvecs(int vec1[], int s1, int vec2[], int s2);
	int difsinvecs(double vec1[], int s1, double vec2[], int s2);
	int difsinvecs(char vec1[], int s1, char vec2[], int s2);
	int difsinvecs(bool vec1[], int s1, bool vec2[], int s2);
  
  double mean_dist_vecs(int *vec1, int s1, int *vec2, int s2);
  double mean_dist_vecs(double *vec1, int s1, double *vec2, int s2);
	
	double multinvec(int vec[], int s); //analogue of productory
	double multinvec(double vec[], int s); // and also for doubles
	
	void fillv0(int vec[], int s);
	void fillv0(double vec[], int s);
	void fillv0(bool vec[], int s);
  void fillv0(string vec[], int s);
	
	void fillv1(int vec[], int s);
	void fillv1(double vec[], int s);
	void fillv1(bool vec[], int s);
  
	void fillvm1(int vec[], int s);
	void fillvm1(double vec[], int s);
	
	int sumatoria(int vec[], int s);
	double sumatoria(double vec[], int s);
	
	void sort_vector(int vecor[], int orden[], int len, int lbo); // create vector with integers to next in sequence
	void sort_vector(double vecor[], int orden[], int len, double lbo); // in original vector
	
  int first_in_vector_with(int pat, int *vec, int siz);
  int first_in_vector_with(string pat, string *vec, int siz);
  int first_in_vector_with(char pat, char *vec, int siz);
  int first_in_vector_with(bool pat, bool *vec, int siz);
  
  
  int first_in_vector_with(int pat, int *vec, int siz, int desde);
  int first_in_vector_with(string pat, string *vec, int siz, int desde);
  int first_in_vector_with(char pat, char *vec, int siz, int desde);
  int first_in_vector_with(bool pat, bool *vec, int siz, int desde);
  
	bool vector_contains(int *vect, int pat, int des, int siz);
	bool vector_contains(int *vect, int pat, int siz);
	
	bool vector_contains(double *vect, double pat, int des, int siz);
	bool vector_contains(double *vect, double pat, int siz);
	
	bool vector_contains(bool *vect, bool pat, int des, int siz);
	bool vector_contains(bool *vect, bool pat, int siz);
	
	bool vector_contains(string *vect, string pat, int des, int siz);
	bool vector_contains(string *vect, string pat, int siz);
	
	bool vector_contains(char *vect, char pat, int des, int siz);
	bool vector_contains(char *vect, char pat, int siz);
	
	bool last_vector(bool *ve, int siz); //checks if vector is full of true values
	bool last_vector(int *ve, int siz, int max);
	
	void next_in_vector(bool *ve, int siz);
	void next_in_vector(int *ve, int siz, int min, int max, int step);
	
  int bintoint(int *vect, int tam, int max); //2024
  void inttobin(int ent, int tam, int max, int step, int *vect); //2024
  
	int count_in_vector(bool *ve, int siz, bool pat);
	int count_in_vector(int *ve, int siz, int pat);
  int count_in_vector(double *ve, int siz, double pat);
  
  void copy_vector(bool *fromthis, bool *tothis, int siz);
  void copy_vector(int *fromthis, int *tothis, int siz);

	//vector formatting
	void printv(ostream& sal, int vec[], int s);
	void printv(ostream& sal, double vec[], int s);
	void printv(ostream& sal, bool vec[], int s);
	void printv(ostream& sal, string vec[], int s);
	
	//vectors and matrices
	int vecinmat(int** mat, int rows, int cols, int vec[], int s);
	int vecinmat(double** mat, int rows, int cols, double vec[], int s);
	int vecinmat(bool** mat, int rows, int cols, bool vec[], int s);
	int vecinmat(string** mat, int rows, int cols, string vec[], int s);
	
	int dotproduct(int *vec1, int *vec2, int siz);
	double dotproduct(int *vec1, double *vec2, int siz);
	double dotproduct(double *vec1, int *vec2, int siz);
	double dotproduct(double *vec1, double *vec2, int siz);
	
	void matxvec(int *vec, int sizv, int leng, int wide, int **ma, int *res);
	void matxvec(int *vec, int sizv, int leng, int wide, double **ma, double *res);
	void matxvec(double *vec, int sizv, int leng, int wide, int **ma, double *res);
	void matxvec(double *vec, int sizv, int leng, int wide, double **ma, double *res);
	
	//matrices
	void fillmat0(int** mat, int rows, int cols);
	void fillmat0(double** mat, int rows, int cols);
	void fillmat0(bool** mat, int rows, int cols);
	
	void fillmatm1(int** mat, int rows, int cols);
	void fillmatm1(double** mat, int rows, int cols);
	bool is_symmetric(int** mat, int rows, int cols);
	bool is_symmetric(double** mat, int rows, int cols);
	bool is_symmetric(bool** mat, int rows, int cols);
	bool is_symmetric(string** mat, int rows, int cols);
	
	bool last_matrix(bool **ve, int rows, int cols); //checks whether matrix is full of true values
	bool last_matrix(int **ve, int rows, int cols, int max);
	
	void next_in_matrix(bool **ve, int rows, int cols);
	void next_in_matrix(int **ve, int rows, int cols, int min, int max, int step);
	
	int sum_column(int **mat, int rows, int cols, int col);
	double sum_column(double **mat, int rows, int cols, int col);
  
  bool eqmatrix(int **mat1, int r1, int c1, int **mat2, int r2, int c2);
  bool eqmatrix(double **mat1, int r1, int c1, double **mat2, int r2, int c2);
  bool eqmatrix(char **mat1, int r1, int c1, char **mat2, int r2, int c2);
  bool eqmatrix(bool **mat1, int r1, int c1, bool **mat2, int r2, int c2);
  bool eqmatrix(string **mat1, int r1, int c1, string **mat2, int r2, int c2);
  
  bool cycle_in_matrix_reps(int **mat, int rows, int cols, int reps);
  
  bool eqvector_rot(int *v1, int s1, int *v2, int s2); //2024
  
  bool eqmatrix_rot(int **mat1, int r1, int c1, int **mat2, int r2, int c2);
  bool eqmatrix_rot(double **mat1, int r1, int c1, double **mat2, int r2, int c2);
  bool eqmatrix_rot(char **mat1, int r1, int c1, char **mat2, int r2, int c2);
  bool eqmatrix_rot(bool **mat1, int r1, int c1, bool **mat2, int r2, int c2);
  bool eqmatrix_rot(string **mat1, int r1, int c1, string **mat2, int r2, int c2);

  double dist_matrices_rot(int **mat1, int r1, int c1, int **mat2, int r2, int c2);
  double dist_matrices_rot_aux(int **mat1, int r1, int c1, int **mat2, int r2, int c2); //equal-sized matrices
	
  void matxmat(int **mat1, int r1, int c1, int **mat2, int r2, int c2, int **mat3, int r3, int c3);
  void matxmat(double **mat1, int r1, int c1, int **mat2, int r2, int c2, double **mat3, int r3, int c3);
  void matxmat(int **mat1, int r1, int c1, double **mat2, int r2, int c2, double **mat3, int r3, int c3);
  void matxmat(double **mat1, int r1, int c1, double **mat2, int r2, int c2, double **mat3, int r3, int c3);

  void add_to_matrix(int **growing, int **addendum, int r, int c);
  void add_to_matrix(double **growing, double **addendum, int r, int c);
  
	//matrix formatting
  void printm(ostream& sal, int** mat, int rows, int cols);
  void printm(ostream& sal, double** mat, int rows, int cols);
  void printm(ostream& sal, bool** mat, int rows, int cols);
  void printm(ostream& sal, string** mat, int rows, int cols);
  
  void printm_latex(ostream& sal, int** mat, int rows, int cols);
  void printm_latex(ostream& sal, double** mat, int rows, int cols);

  
	//sets
	set<int> merge(set<int> &a, set<int> &b);
	set<char> merge(set<char> &a, set<char> &b);
	set<string> merge(set<string> &a, set<string> &b);
  
	
	set<int> intersect(set<int> &a, set<int> &b);
	set<char> intersect(set<char> &a, set<char> &b);
	set<string> intersect(set<string> &a, set<string> &b);
	
	void printset(set<int> &ise, ostream& sal);
	void printset(set<string> &ise, ostream& sal);
  void printsetofsets(set<set<int> > &con, ostream& sal);
	bool contains(char c, set<char> &sc);
	bool contains(int c, set<int> &sc);
	bool contains(string c, set<string> &sc);
  void remove_nth(set<string>& conj, int n);
  string return_nth(set<string> &conj, int n);
	
  void append_sets_of_sets(const set<set<int > > &source, set<set<int > > &target);
  
	
	//lists
	bool contains(int num, list<int> &lis);
  bool contains(string num, list<string> &lis);
	void printlist(list<int> &lis, ostream& sal);
  void printlist(list<string> &lis, ostream& sal);
	
	void printlistoflists(list<list<int> > &lis, ostream& sal);
	double productory(list<int> con); //analogue of multinvec, for lists
	int factorial(int num);
	double binom_coeff(int ene, int ka);
	
	//others
  void run_command(string cuerda);
  
  void create_array(int** &arr, int rows, int cols);
  void create_array(bool** &arr, int rows, int cols);
  void create_array(double** &arr, int rows, int cols);
  void create_array(char** &arr, int rows, int cols);
  void create_array(string** &arr, int rows, int cols);
  
  void create_array(int*** &arr, int slices, int rows, int cols);
  void create_array(bool*** &arr, int slices, int rows, int cols);
  void create_array(double*** &arr, int slices, int rows, int cols);
  void create_array(char*** &arr, int slices, int rows, int cols);
  void create_array(string*** &arr, int slices, int rows, int cols);
  
	void get_dot_fig(string arch);
	void get_dot_fig(string arch, string ext);
	void get_spring_dot_fig(string arch);
	void get_spring_dot_fig(string arch, string ext);
	
	
	int round(double x);
	
	void polar_to_cartesian(double radians, double radius, double& x, double& y);
	void open_ifstream(ifstream& fe, string nomb);
	void open_ofstream(ofstream& fs, string nomb);
	void open_ofstream_to_append(ofstream& fs, string nomb);
	
  bool dir_exists(string dir);
	int count_files_in_dir(string dir);
	void get_files_in_dir(string dir, string* vecnom);
	
	//strings
	void erase_between_braces(char pator[], int len, ostream& fs, ifstream& fe);
	int char_in_string(char c, const string& cue, int from, int until);
	int char_in_string(char c, const string& cue);
	string del_char(char c, string ori);
	int count_char_in_string(char c, const string& cue);
	string chop_word_from_line(string& line);
	bool string_in_string(string pat, string big);
	string inttostring(int num);
	string doubletostring_tex(double num);
	bool string_in_vector(string* vec, string ele, int tam);
  int where_string_in_vector(string* vec, string ele, int tam);
	bool nothing_but_spaces(const string& cue);
	void capitalize(string& word);
	void capitalize_sentence(string& sentence);
	bool eq_string(string uno, string dos, bool capsmatter);
	int nth_app_of_char(char c, const string& cue, int n);
	
	
	//statistics
	double get_mean(int* vec, int tam);
	double get_mean(double* vec, int tam);
  double get_mean(int* vec, int tam, set<int> &cua);
  double get_mean(double* vec, int tam, set<int> &cua);
  double get_mean_by_col(int **mat, int rows, int cols, int col);
  double get_mean_by_col(double **mat, int rows, int cols, int col);
	double get_sample_variance(int* vec, int tam);
	double get_sample_variance(double* vec, int tam);
  double get_sample_variance(int* vec, int tam, set<int> &cua);
  double get_sample_variance(double* vec, int tam, set<int> &cua);
  double get_sample_variance(int* vec, int tam, double ave);
  double get_sample_variance(double* vec, int tam, double ave);
  double get_sample_variance(int* vec, int tam, double ave, set<int> &cua);
  double get_sample_variance(double* vec, int tam, double ave, set<int> &cua);
  double get_pop_variance(int* vec, int tam);
  double get_pop_variance(double* vec, int tam);
  double get_pop_variance(int* vec, int tam, set<int> &cua);
  double get_pop_variance(double* vec, int tam, set<int> &cua);
  double get_pop_variance(int* vec, int tam, double ave);
  double get_pop_variance(double* vec, int tam, double ave);
  double get_pop_variance(int* vec, int tam, double ave, set<int> &cua);
  double get_pop_variance(double* vec, int tam, double ave, set<int> &cua);
  double get_sample_stddev(int* vec, int tam);
  double get_sample_stddev(double* vec, int tam);
	double get_sample_stddev(int* vec, int tam, set<int> &cua);
	double get_sample_stddev(double* vec, int tam, set<int> &cua);
  double get_sample_stddev(int* vec, int tam, double ave);
  double get_sample_stddev(double* vec, int tam, double ave);
  double get_sample_stddev(int* vec, int tam, double ave, set<int> &cua);
  double get_sample_stddev(double* vec, int tam, double ave, set<int> &cua);
	double get_pop_stddev(int* vec, int tam);
	double get_pop_stddev(double* vec, int tam);
  double get_pop_stddev(int* vec, int tam, set<int> &cua);
  double get_pop_stddev(double* vec, int tam, set<int> &cua);
  double get_pop_stddev(int* vec, int tam, double ave);
  double get_pop_stddev(double* vec, int tam, double ave);
  double get_pop_stddev(int* vec, int tam, double ave, set<int> &cua);
  double get_pop_stddev(double* vec, int tam, double ave, set<int> &cua);
  double get_sample_stderr(int* vec, int tam);
  double get_sample_stderr(double* vec, int tam);
  double get_sample_stderr(int* vec, int tam, set<int> &cua);
  double get_sample_stderr(double* vec, int tam, set<int> &cua);
  double get_sample_stderr(int* vec, int tam, double ave);
  double get_sample_stderr(double* vec, int tam, double ave);
  double get_sample_stderr(int* vec, int tam, double ave, set<int> &cua);
  double get_sample_stderr(double* vec, int tam, double ave, set<int> &cua);
	double get_samp_var_err_of_mean(int* vec, int tam);
	double get_samp_var_err_of_mean(double* vec, int tam);
  double get_samp_var_err_of_mean(int* vec, int tam, set<int> &cua);
  double get_samp_var_err_of_mean(double* vec, int tam, set<int> &cua);
	
	
	double find_min(double* vec, int tam);
	int find_min(int* vec, int tam);
  
  double find_min(double* vec, int tam, set<int> &cua);
  int find_min(int* vec, int tam, set<int> &cua);
  
	int find_minlb(int* vec, int tam, int lbo); //
	int find_minlb(double* vec, int tam, double lbo); //
	double find_max(double* vec, int tam);
  int find_max(int* vec, int tam);
  
  double find_max(double* vec, int tam, set<int> &cua);
  int find_max(int* vec, int tam, set<int> &cua);
  int find_max_index(double* vec, int tam);
  int find_max_index(int* vec, int tam);
  int find_min_index(double* vec, int tam);
  int find_min_index(int* vec, int tam);
	void sort(double* vec, double* nvec, int tam);
  void sort(double* vec, double* nvec, int tam, set<int> &cua);
	void sort(int* vec, int* nvec, int tam);
  void sort(int* vec, int* nvec, int tam, set<int> &cua);
	double get_median(double* vec, int tam);
	double get_median(int* vec, int tam);
  
  double get_median(double* vec, int tam, set<int> &cua);
  double get_median(int* vec, int tam, set<int> &cua);
  
  double get_midpoint(int* vo, int tam); //mediana en vector ordenado
  double get_midpoint(double* vo, int tam); //mediana en vector ordenado
  double get_q1(int* vo, int tam);
  double get_q1(double* vo, int tam);
  double get_q3(int* vo, int tam);
  double get_q3(double* vo, int tam);
  
  
  
	double histo2(int nudata, int nubins, double* data1, double* data2, double* dis1, double* dis2, double* xax1, double* xax2);
  double histo2(int nudata1, int nudata2, int nubins, double* data1, double* data2, double* dis1, double* dis2, double* xax1, double* xax2);
	double histo2(int nudata, int nubins, int* data1, int* data2, double* dis1, double* dis2, double* xax1, double* xax2);
	double histo(int nudata, int nubins, int* data, double* dis, double* xax, int min, int max);
	int histo_int(int nudata, int nubins, int* data, double* dis, double* xax, int min, int max);
	double histo(int nudata, int nubins, double* data, double* dis, double* xax, double min, double max);
	double histo(int nudata, int nubins, int* data, double* dis, double* xax);
	double histo(int nudata, int nubins, double* data, double* dis, double* xax);
	int histo_int(int nudata, int nubins, int* data, double* dis, double* xax);
	
	int fdr(double* vps, int n, int* orden, double alpha, bool* cualessi); //Based on Benjamini and Hochberg 1995
	int fdr_dep(double* vps, int n, int* orden, double alpha, bool* cualessi); //Based on Benjamini and Yerkutieli 2001
	double sum_consec_inverse(int from, int to);
	double from_Z_to_p(double zscore, char gol); //one-tailed
  double from_Z_to_p(double zscore); //two-tailed
  double from_t_to_p(double t, int df, char gol); //one-tailed
  double from_t_to_p(double t, int df); //two-tailed
	double variance_binomial(double p, int n); //theoretical var for binomial distribution //Le2003 //CHECAR EL /N... Lomax &HV llaman a esto variance error of the proportion (p.163)
	double stddev_binomial(double p, int n); //theoretical stdev for binomial distribution // Lomax &HV llaman a esto standard error of the proportion
	double zscore_binary(double obsfreq, double expfreq, int n);//Le2003, pags 208-9//revisar stddev debe ser de pob no de sample segun el otro. Hay variacion con stderr en lugar de sd cuando es una sola media... segun el otro tambien. ¿cuando se usa std err de la muestra es t?
	double pearsons_chi(double& chi, int& df, double **obs, double **exp, int rows, int cols); //Le2003, pags 223-
	double pearsons_chi(double& chi, int& df, double **condmat, int rows, int cols); //Le2003, pags 223-
	double yates_chi(double& chi, int& df, double **obs, double **exp, int rows, int cols); //Le2003, pags 226-
	double yates_chi(double& chi, int& df, double **condmat, int rows, int cols); //Le2003, pags 226-
  //aux for stats
	int get_nubins(int *vec, int tvec); //for histo_int
  
  double kolmogorov_smirnov(double& D, double* v1, double* v2, int tam1, int tam2);
  double kolmogorov_smirnov(double& D, int* v1, int* v2, int tam1, int tam2);
  double pearsons_r(double& r, double* veh, double* voh, int tam, char tgol, double& S, int& df); //pearson´s product moment correlation coefficient
  double pearsons_r(double& r, double* veh, int* voh, int tam, char tgol, double& S, int& df); //pearson´s product moment correlation coefficient
  double pearsons_r(double& r, int* veh, int* voh, int tam, char tgol, double& S, int& df); //pearson´s product moment correlation coefficient
  double spearmans_rho(double& r, double* veh, double* voh, int tam, char tgol, double& S, int& df); //spearman rho
  double spearmans_rho(double& r, double* veh, int* voh, int tam, char tgol, double& S, int& df); //spearman rho
  double spearmans_rho(double& r, int* veh, int* voh, int tam, char tgol, double& S, int& df); //spearman rho
	double shapiro_wilk(double& W, double* vect, int tam); //Shapiro-Wilk normality test, uses R, p-values higher than alpha mean that distribution is normal
	double shapiro_wilk(double& W, int* vect, int tam);
	double brown_forsythe(double& F, int& dfn, int& dfd, double **groupsonrows, int *tams, int numbgroups);//homoscedasticity test. Lomax & Hahs-Vaughn 2012. section 9.4.2
  double levene(double& F, int& dfn, int& dfd, double **groupsonrows, int *tams, int numbgroups);//homoscedasticity test. (best for symmetric moderate-tailed

	double students_t(double& T, int& df, double *vec1, int tam1, double *vec2, int tam2); //Lomax&HV. section 7.2.1
	double students_t(double& T, int& df, double *vec1, double *vec2, int tam);//equal sample size. Para una cola, v1 > v2... Si quiero dos colas, mult p por 2. Si quiero v2>v1 1-p
  double students_t(double& T, int& df, double *vec, int tam, double mu); //para difs con respecto a mu
  double paired_students_t(double& T, int& df, double *vec1, double *vec2, int tam);
  
	double welch_t(double& T, int& df, double *vec1, int tam1, double *vec2, int tam2); //Welch t test. Lomax & HV section 7.2.1
	double mann_whitney(double *vec1, int tam1, double *vec2, int tam2); // Mann-Whitney U test. Le2003, section 7.4
	double mann_whitney(int *vec1, int tam1, int *vec2, int tam2); // Mann-Whitney U test. Le2003, section 7.4
	double zscore(double *vec, int tam, double mu); //z-score. Sec. 4.2 de Lomax & HV
	double zscore(int *vec, int tam, int mu); //
	
  void transpose(int r, int c, double **mator, double **mattrans);
  void append_idmat(int r, int c, double **mator, double **matfin);
  void order_for_echelon(int r, int c, double **mator, double **matfin);
  
  ////Eigenvalues and eigen vectors
  //Symmetric matrices!!!
  void eigenvv_sym(double **mat, int msize, double **evecs, double *evals); //matrix is symmetric. mat[i] refers to the ith row of mat, not to the ith column
  void eigenvv_sym(int **mat, int msize, double **evecs, double *evals); //matrix is symmetric. mat[i] refers to the ith row of mat, not to the ith column
  double get_leading_evector_symmat(double **mat, int msize, double *evector);
  double get_leading_evector_symmat(int **mat, int msize, double *evector);

  double pow_meth(double **mat, int msize, double *evec);
  
  
  
	//para gsl y distribucion normal: si z > 37.515, entonces p < 2.62271 \times 10^{-308}
	Alea est;
	
private:
};

#endif



