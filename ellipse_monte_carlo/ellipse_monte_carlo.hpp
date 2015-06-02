void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );
double ddot ( int n, double dx[], int incx, double dy[], int incy );
double ellipse_area1 ( double a[], double r );
double ellipse_area2 ( double a, double b, double c, double d );
double *ellipse_sample ( int n, double a[], double r, int &seed );
double *monomial_value ( int m, int n, int e[], double x[] );
double r8_uniform_01 ( int &seed );
int r8po_fa ( double a[], int lda, int n );
void r8po_sl ( double a[], int lda, int n, double b[] );
void r8vec_normal_01 ( int n, int &seed, double x[] );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int &seed );
void timestamp ( );
double *uniform_in_sphere01_map ( int dim_num, int n, int &seed );

