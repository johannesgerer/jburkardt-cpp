double dznrm2 ( int n, complex <double> x[], int incx );
double dzasum ( int n, complex <double> x[], int incx );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int izamax ( int n, complex <double> x[], int incx );
bool lsame ( char ca, char cb );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_sign ( double x );
void xerbla ( char *srname, int info );
double zabs1 ( complex <double> z );
double zabs2 ( complex <double> z );
void zaxpy ( int n, complex <double> ca, complex <double> cx[], 
  int incx, complex <double> cy[], int incy );
void zcopy ( int n, complex <double> cx[], int incx, complex <double> cy[], 
  int incy );
complex <double> zdotc ( int n, complex <double> cx[], int incx, 
  complex <double> cy[], int incy );
complex <double> zdotu ( int n, complex <double> cx[], int incx, 
  complex <double> cy[], int incy );
void zdrot ( int n, complex <double> cx[], int incx, complex <double> cy[], 
  int incy, double c, double s );
void zdscal ( int n, double sa, complex <double> cx[], int incx );
double zmach ( int job );
void zrotg ( complex <double> *ca, complex <double> cb, double *c, 
  complex <double> *s );
void zscal ( int n, complex <double> ca, complex <double> cx[], int incx );
complex <double> zsign1 ( complex <double> z1, complex <double> z2 );
complex <double> zsign2 ( complex <double> z1, complex <double> z2 );
void zswap ( int n, complex <double> cx[], int incx, complex <double> cy[], 
  int incy );
