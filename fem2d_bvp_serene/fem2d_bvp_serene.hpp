double *basis_serene ( double xq, double yq, double xw, double ys, double xe, 
  double yn, double xx[], double yy[] );
double *basis_dx_serene ( double xq, double yq, double xw, double ys, 
  double xe, double yn, double xx[], double yy[] );
double *basis_dy_serene ( double xq, double yq, double xw, double ys, 
  double xe, double yn, double xx[], double yy[] );
double *fem2d_bvp_serene ( int nx, int ny, double a ( double x, double y ), 
  double c ( double x, double y ), double f ( double x, double y ), 
  double x[], double y[], bool show11 );
int fem2d_bvp_serene_node_num ( int nx, int ny );
double fem2d_h1s_error_serene ( int nx, int ny, double x[], double y[], 
  double u[], double exact_ux ( double x, double y ), 
  double exact_uy ( double x, double y ) );
double fem2d_l1_error_serene ( int nx, int ny, double x[], double y[], 
  double u[], double exact ( double x, double y ) );
double fem2d_l2_error_serene ( int nx, int ny, double x[], double y[], 
  double u[], double exact ( double x, double y ) );
int *i4vec_zero_new ( int n );
double not1 ( double x1, double x2, double x3 );
double not1d ( double x2, double x3 );
double not2 ( double x1, double y1, double x2, double y2, double x3, double y3, 
  double x4, double y4 );
double not2dx ( double x2, double y2, double x3, double y3, double x4, double y4 );
double not2dy ( double x2, double y2, double x3, double y3, double x4, double y4 );
double r8_uniform_01 ( int &seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8mat_solve2 ( int n, double a[], double b[], int &ierror );
double *r8mat_zero_new ( int m, int n );
double *r8vec_linspace_new ( int n, double a, double b );
double r8vec_sum ( int n, double a[] );
double *r8vec_zero_new ( int n );
void timestamp ( );
double *wathen ( int nx, int ny, int n );
