double *fd2d_heat_steady ( int nx, int ny, double x[], double y[], 
  double d ( double x, double y ), double f ( double x, double y ) );
void interior ( int nx, int ny, double x[], double y[], 
  double d ( double x, double y ), double f ( double x, double y ), int n, 
  double a[], double rhs[] );
double r8_abs ( double x );
void r8mat_fs ( int n, double a[], double x[] );
double *r8vec_linspace_new ( int n, double a, double b );
void r8vec_mesh_2d ( int nx, int ny, double xvec[], double yvec[], 
  double xmat[], double ymat[] );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );

