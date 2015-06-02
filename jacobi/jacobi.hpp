double *dif2 ( int m, int n );
double *jacobi1 ( int n, double a[], double b[], double x[] );
double *r8mat_mv_new ( int m, int n, double a[], double x[] );
double r8mat_residual_norm ( int m, int n, double a[], double x[], double b[] );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_diff_norm_squared ( int n, double a[], double b[] );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );
