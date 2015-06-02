void grid_2d ( int x_num, double x_lo, double x_hi, int y_num, double y_lo, 
  double y_hi, double x[], double y[] );
double r8vec_amax ( int n, double a[] );
double r8vec_amin ( int n, double a[] );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );
double *r8vec_uniform_ab_new ( int n, double a, double b, int &seed );
void resid_spiral ( int n, double x[], double y[], double c, double pr[] );
void spiral_gnuplot ( string header, int n, double x[], double y[], double u[], 
  double v[], double s );
void timestamp ( );
void uv_spiral ( int n, double x[], double y[], double c, double u[], 
  double v[] );
