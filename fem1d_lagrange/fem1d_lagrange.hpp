void fem1d_lagrange_stiffness ( int x_num, double x[], int q_num, 
  double f ( double x ), double a[], double m[], double b[] );
double *lagrange_derivative ( int nd, double xd[], int ni, double xi[] );
double *lagrange_value ( int nd, double xd[], int ni, double xi[] );
void legendre_set ( int n, double x[], double w[] );
double *r8mat_fs_new ( int n, double a[], double b[] );
double *r8vec_linspace_new ( int n, double a, double b );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );
