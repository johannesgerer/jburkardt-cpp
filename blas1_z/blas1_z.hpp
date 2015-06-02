double dznrm2 ( int n, complex <double> x[], int incx );
double dzasum ( int n, complex <double> x[], int incx );
int izamax ( int n, complex <double> x[], int incx );
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
void zrotg ( complex <double> *ca, complex <double> cb, double *c, 
  complex <double> *s );
void zscal ( int n, complex <double> ca, complex <double> cx[], int incx );
void zswap ( int n, complex <double> cx[], int incx, complex <double> cy[], 
  int incy );
