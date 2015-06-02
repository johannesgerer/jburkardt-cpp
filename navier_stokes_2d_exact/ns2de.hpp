void grid_2d ( int x_num, double x_lo, double x_hi, int y_num, double y_lo, 
  double y_hi, double x[], double y[] );
void ns2de_gnuplot ( string header, int n, double x[], double y[], double u[], 
  double v[], double s );
double r8vec_amax ( int n, double a[] );
double r8vec_amin ( int n, double a[] );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );
void r8vec_linspace ( int n, double a_lo, double a_hi, double x[] );
double r8vec_norm_l2 ( int n, double a[] );
double *r8vec_uniform_ab_new ( int n, double a, double b, int &seed );
void resid_lucas ( double nu, double rho, int n, double x[], double y[], 
  double t, double ur[], double vr[], double pr[] );
void resid_spiral ( double nu, double rho, int n, double x[], double y[], 
  double t, double ur[], double vr[], double pr[] );
void resid_taylor ( double nu, double rho, int n, double x[], double y[], 
  double t, double ur[], double vr[], double pr[] );
void rhs_lucas ( double nu, double rho, int n, double x[], double y[], double t, 
  double f[], double g[], double h[] );
void rhs_spiral ( double nu, double rho, int n, double x[], double y[], 
  double t, double f[], double g[], double h[] );
void rhs_taylor ( double nu, double rho, int n, double x[], double y[], 
  double t, double f[], double g[], double h[] );
void timestamp ( );
void uvp_lucas ( double nu, double rho, int n, double x[], double y[], 
  double t, double u[], double v[], double p[] );
void uvp_spiral ( double nu, double rho, int n, double x[], double y[], 
  double t, double u[], double v[], double p[] );
void uvp_taylor ( double nu, double rho, int n, double x[], double y[], 
  double t, double u[], double v[], double p[] );
