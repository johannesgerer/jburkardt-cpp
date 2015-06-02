complex <float> c4_uniform_01 ( int &seed );
void c4mat_print ( int m, int n, complex <float> a[], string title );
void c4mat_print_some ( int m, int n, complex <float> a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
complex <float> *c4mat_test ( int n );
complex <float> *c4mat_test_inverse ( int n );
complex <double> c8_uniform_01 ( int &seed );
void c8mat_print ( int m, int n, complex <double> a[], string title );
void c8mat_print_some ( int m, int n, complex <double> a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
complex <double> *c8mat_test ( int n );
complex <double> *c8mat_test_inverse ( int n );
float cabs1 ( complex <float> z );
float cabs2 ( complex <float> z );
float cmach ( int job );
complex <float> csign1 ( complex <float> z1, complex <float> z2 );
complex <float> csign2 ( complex <float> z1, complex <float> z2 );
double dmach ( int job );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
bool lsame ( char ca, char cb );
float r4_abs ( float x );
float r4_max ( float x, float y );
float r4_sign ( float x );
float r4_uniform_01 ( int &seed );
float r4_uniform_ab ( float a, float b, int &seed );
void r4mat_print ( int m, int n, float a[], string title );
void r4mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, string title );
float *r4mat_test ( char trans, int lda, int m, int n );
void r4vec_print ( int n, float a[], string title );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_sign ( double x );
double r8_uniform_01 ( int &seed );
double r8_uniform_ab ( double a, double b, int &seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8mat_test ( char trans, int lda, int m, int n );
void r8vec_print ( int n, double a[], string title );
float smach ( int job );
void timestamp ( );
void xerbla ( string srname, int info );
double zabs1 ( complex <double> z );
double zabs2 ( complex <double> z );
double zmach ( int job );
complex <double> zsign1 ( complex <double> z1, complex <double> z2 );
complex <double> zsign2 ( complex <double> z1, complex <double> z2 );
