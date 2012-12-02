double *fd1d_heat_explicit ( int x_num, double x[], double t, double dt, 
  double cfl, double *rhs ( int x_num, double x[], double t ), 
  void bc ( int x_num, double x[], double t, double h[] ), double h[] );
double fd1d_heat_explicit_cfl ( double k, int t_num, double t_min, double t_max, 
  int x_num, double x_min, double x_max );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double *r8vec_linspace_new ( int n, double a_first, double a_last );
void r8vec_write ( string output_filename, int n, double x[] );
void timestamp ( );
