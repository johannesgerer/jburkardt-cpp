void grid_2d ( int x_num, double x_lo, double x_hi, int y_num, double y_lo, 
  double y_hi, double x[], double y[] );
double r8vec_amax ( int n, double a[] );
double r8vec_amin ( int n, double a[] );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );
double r8vec_norm_l2 ( int n, double a[] );
double *r8vec_uniform_ab_new ( int n, double a, double b, int &seed );
void resid_stokes1 ( int n, double x[], double y[], double ur[], double vr[], 
  double pr[] );
void resid_stokes2 ( int n, double x[], double y[], double ur[], double vr[], 
  double pr[] );
void resid_stokes3 ( int n, double x[], double y[], double ur[], double vr[], 
  double pr[] );
void rhs_stokes1 ( int n, double x[], double y[], double f[], double g[], 
  double h[] );
void rhs_stokes2 ( int n, double x[], double y[], double f[], double g[], 
  double h[] );
void rhs_stokes3 ( int n, double x[], double y[], double f[], double g[], 
  double h[] );
void stokes_gnuplot ( string header, int n, double x[], double y[], double u[], 
  double v[], double s );
void timestamp ( );
void uvp_stokes1 ( int n, double x[], double y[], double u[], double v[], 
  double p[] );
void uvp_stokes2 ( int n, double x[], double y[], double u[], double v[], 
  double p[] );
void uvp_stokes3 ( int n, double x[], double y[], double u[], double v[], 
  double p[] );
