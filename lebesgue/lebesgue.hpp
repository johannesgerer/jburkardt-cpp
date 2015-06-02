double *chebyshev1 ( int n );
double *chebyshev2 ( int n );
double *chebyshev3 ( int n );
double *chebyshev4 ( int n );
double *equidistant1 ( int n );
double *equidistant2 ( int n );
double *equidistant3 ( int n );
double *fejer1 ( int n );
double *fejer2 ( int n );
double *lagrange_value ( int data_num, double t_data[], int interp_num, 
  double t_interp[] );
double lebesgue_constant ( int n, double x[], int nfun, double xfun[] );
double *lebesgue_function ( int n, double x[], int nfun, double xfun[] );
void lebesgue_plot ( int n, double x[], int nfun, double xfun[], 
  string label, string filename );
double *r8vec_linspace_new ( int n, double a, double b );
double r8vec_max ( int n, double r8vec[] );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );
