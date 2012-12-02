int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int isamax ( int n, float dx[], int incx );
bool lsame ( char ca, char cb );
float r4_abs ( float x );
float r4_max ( float x, float y );
float r4_sign ( float x );
float sasum ( int n, float x[], int incx );
void saxpy ( int n, float da, float dx[], int incx, float dy[], int incy );
void scopy ( int n, float dx[], int incx, float dy[], int incy );
float sdot ( int n, float dx[], int incx, float dy[], int incy );
float smach ( int job );
float snrm2 ( int n, float x[], int incx );
void srot ( int n, float x[], int incx, float y[], int incy, float c, 
  float s );
void srotg ( float *sa, float *sb, float *c, float *s );
void sscal ( int n, float sa, float x[], int incx );
void sswap ( int n, float x[], int incx, float y[], int incy );
void xerbla ( char *srname, int info );
