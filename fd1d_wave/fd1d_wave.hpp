double fd1d_wave_alpha ( int x_num, double x1, double x2, int t_num, double t1, 
  double t2, double c );
double *fd1d_wave_start ( int x_num, double x_vec[], double t, double t_delta, 
  double alpha, double u_x1 ( double t ), double u_x2 ( double t ), 
  double *ut_t1 ( int x_num, double x_vec[] ), double u1[] );
double *fd1d_wave_step ( int x_num, double t, double alpha, 
  double u_x1 ( double t ), double u_x2 ( double t ), double u1[], double u2[] );
double *piecewise_linear ( int nd, double xd[], double yd[], int nv, double xv[] );
double r8_abs ( double x );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double *r8vec_linspace_new ( int n, double a_first, double a_last );
void timestamp ( );
