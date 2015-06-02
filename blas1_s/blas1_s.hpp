int isamax ( int n, float dx[], int incx );
float sasum ( int n, float x[], int incx );
void saxpy ( int n, float da, float dx[], int incx, float dy[], int incy );
void scopy ( int n, float dx[], int incx, float dy[], int incy );
float sdot ( int n, float dx[], int incx, float dy[], int incy );
float snrm2 ( int n, float x[], int incx );
void srot ( int n, float x[], int incx, float y[], int incy, float c, 
  float s );
void srotg ( float *sa, float *sb, float *c, float *s );
void sscal ( int n, float sa, float x[], int incx );
void sswap ( int n, float x[], int incx, float y[], int incy );

