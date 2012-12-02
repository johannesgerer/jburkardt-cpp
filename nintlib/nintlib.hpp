double box_nd ( double func ( int dim_num, double x[] ), int dim_num, 
  int order, double xtab[], double weight[], int *eval_num );
int i4_huge ( void );
int i4_power ( int i, int j );
double monte_carlo_nd ( double func ( int dim_num, double x[] ), int dim_num, 
  double a[], double b[], int eval_num, int *seed );
double p5_nd ( double func ( int dim_num, double x[] ), int dim_num, 
  double a[], double b[], int *eval_num );
double r8_abs ( double x );
double r8_epsilon ( void );
double *r8vec_uniform_01_new ( int n, int *seed );
double romberg_nd ( double func ( int dim_num, double x[] ), double a[], 
  double b[], int dim_num, int sub_num[], int it_max, double tol, int *ind, 
  int *eval_num );
void sample_nd ( double func ( int dim_num, double x[] ), int k1, int k2, 
  int dim_num, double est1[], double err1[], double dev1[], double est2[], 
  double err2[], double dev2[], int *eval_num );
double sum2_nd ( double func ( int dim_num, double x[] ), double xtab[], 
  double weight[], int order[], int dim_num, int *eval_num );
void timestamp ( void );
void tuple_next ( int m1, int m2, int n, int *rank, int x[] );
