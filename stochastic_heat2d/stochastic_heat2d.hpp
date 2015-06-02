double diffusivity_2d_bnt ( double dc0, double omega[], double x, 
  double y );
void interior ( double omega[], int nx, int ny, double x[], double y[], 
  double f ( double x, double y ), int n, double a[], double rhs[] );
double r8_uniform_01 ( int *seed );
void r8mat_fs ( int n, double a[], double x[] );
double r8mat_max ( int m, int n, double a[] );
double r8mat_mean ( int m, int n, double a[] );;
double *r8vec_linspace_new ( int n, double a, double b );
void r8vec_mesh_2d ( int nx, int ny, double xvec[], double yvec[], 
  double xmat[], double ymat[] );
double *r8vec_normal_01_new ( int n, int &seed );
void r8vec_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );
double *stochastic_heat2d ( double omega[], int nx, int ny, double x[], 
  double y[], double f ( double x, double y ) );
void timestamp ( );

