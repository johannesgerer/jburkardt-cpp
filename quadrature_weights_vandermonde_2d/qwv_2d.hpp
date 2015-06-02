int *i4vec_zero_new ( int n );
double *qwv_2d ( int t, int n, double a, double b, double c, double d, 
  double x[], double y[] );
double r8_abs ( double x );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8mat_solve2 ( int n, double a[], double b[], int &ierror );
double *r8vec_even_new ( int n, double alo, double ahi );
void r8vec_print ( int n, double a[], string title );
void r8vec_print_16 ( int n, double a[], string title );
double *r8vec_zero_new ( int n );
void r8vec2_print ( int n, double a1[], double a2[], string title );
void timestamp ( );
