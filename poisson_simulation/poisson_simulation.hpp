int i4_min ( int i1, int i2 );
int i4vec_max ( int n, int a[] );
double i4vec_mean ( int n, int x[] );
int i4vec_min ( int n, int a[] );
void i4vec_print ( int n, int a[], string title );
double i4vec_variance ( int n, int x[] );
void poisson_fixed_events ( double lambda, int event_num, int &seed, 
  double t[], double w[] );
int poisson_fixed_time ( double lambda, double time, int &seed );
double r8_uniform_01 ( int &seed );
void r8vec_cum ( int n, double a[], double a_cum[] );
double r8vec_max ( int n, double r8vec[] );
double r8vec_mean ( int n, double x[] );
double *r8vec_midspace_new ( int n, double a, double b );
double r8vec_min ( int n, double r8vec[] );
double *r8vec_uniform_01_new ( int n, int &seed );
double r8vec_variance ( int n, double x[] );
void timestamp ( );
