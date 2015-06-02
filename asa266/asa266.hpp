double alngam ( double xvalue, int *ifault );
double alnorm ( double x, bool upper );
double alogam ( double x, int *ifault );
double digamma ( double x );
void dirichlet_check ( int n, double a[] );
void dirichlet_estimate ( int k, int n, double x[], int ix, int init, 
  double alpha[], double &rlogl, double v[], double g[], int &niter, 
  double &s, double &eps, int &ifault );
double *dirichlet_mean ( int n, double a[] );
double *dirichlet_mix_mean ( int comp_max, int comp_num, int elem_num, 
  double a[], double comp_weight[] );
double *dirichlet_mix_sample ( int comp_max, int comp_num, int elem_num, 
  double a[], double comp_weight[], int &seed, int &comp );
double *dirichlet_sample ( int n, double a[], int &seed );
double *dirichlet_variance ( int n, double a[] );
double exponential_01_sample ( int &seed );
double exponential_cdf_inv ( double cdf, double a, double b );
double gamain ( double x, double p, int *ifault );
double gamma_sample ( double a, double b, int &seed );
double gammad ( double x, double p, int *ifault );
double gammds ( double x, double p, int *ifault );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double lngamma ( double z, int *ier );
void normp ( double z, double *p, double *q, double *pdf );
void nprob ( double z, double *p, double *q, double *pdf );
double ppchi2 ( double p, double v, double g, int *ifault );
double ppnd ( double p, int *ifault );
double ppnd16 ( double p, int *ifault );
double r8_epsilon ( );
double r8_gamma_log ( double x );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_normal_01 ( int &seed );
double r8_psi ( double xx );
double r8_uniform_01 ( int &seed );
double r8_uniform_ab ( double a, double b, int &seed );
double *r8col_mean ( int m, int n, double a[] );
double *r8col_variance ( int m, int n, double a[] );
double *r8mat_mv_new ( int m, int n, double a[], double x[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double r8poly_value ( int m, double c[], double x );
double r8vec_dot_product ( int n, double a1[], double a2[] );
void r8vec_print ( int n, double a[], string title );
double r8vec_sum ( int n, double a[] );
void timestamp ( );
double trigamma ( double x, int *ifault );
