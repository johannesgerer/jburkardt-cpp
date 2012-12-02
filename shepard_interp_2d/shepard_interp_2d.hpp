double r8_abs ( double x );
double *r8vec_copy_new ( int n, double a1[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double *r8vec_linspace_new ( int n, double a, double b );
double r8vec_max_dup ( int n, double r8vec[] );
double r8vec_min_dup ( int n, double r8vec[] );
double r8vec_norm_affine ( int n, double v0[], double v1[] );
double r8vec_sum_dup ( int n, double a[] );
void r8vec3_print ( int n, double a1[], double a2[], double a3[], string title );
double *shepard_interp_2d ( int nd, double xd[], double yd[], double zd[],
  double p, int ni, double xi[], double yi[] );
