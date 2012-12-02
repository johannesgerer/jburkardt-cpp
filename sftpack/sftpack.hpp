void c4mat_print_some ( int m, int n, complex <float> a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
complex <float> *c4mat_sftb ( int n1, int n2, complex <float> y[] );
complex <float> *c4mat_sftf ( int n1, int n2, complex <float> x[] );
complex <float> *c4mat_uniform_01_new ( int m, int n, int *seed );

void c4vec_print_part ( int n, complex <float> a[], int max_print, 
  string title );
complex <float> *c4vec_sftb ( int n, complex <float> y[] );
complex <float> *c4vec_sftf ( int n, complex <float> x[] );
complex <float> *c4vec_uniform_01_new ( int n, int *seed );

void c8mat_print_some ( int m, int n, complex <double> a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
complex <double> *c8mat_sftb ( int n1, int n2, complex <double> y[] );
complex <double> *c8mat_sftf ( int n1, int n2, complex <double> x[] );
complex <double> *c8mat_uniform_01_new ( int m, int n, int *seed );

void c8vec_print_part ( int n, complex <double> a[], int max_print, 
  string title );
complex <double> *c8vec_sftb ( int n, complex <double> y[] );
complex <double> *c8vec_sftf ( int n, complex <double> x[] );
complex <double> *c8vec_uniform_01_new ( int n, int *seed );

int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );

void r4vec_print_part ( int n, float a[], int max_print, string title );
float *r4vec_sct ( int n, float x[] );
float *r4vec_sftb ( int n, float azero, float a[], float b[] );
void r4vec_sftf ( int n, float r[], float *azero, float a[], float b[] );
float *r4vec_sht ( int n, float a[]  );
float *r4vec_sqctb ( int n, float x[] );
float *r4vec_sqctf ( int n, float x[] );
float *r4vec_sqstb ( int n, float x[] );
float *r4vec_sqstf ( int n, float x[] );
float *r4vec_sst ( int n, float x[] );
float *r4vec_uniform_new ( int n, float b, float c, int *seed );

void r8vec_print_part ( int n, double a[], int max_print, string title );
double *r8vec_sct ( int n, double x[] );
double *r8vec_sftb ( int n, double azero, double a[], double b[] );
void r8vec_sftf ( int n, double r[], double *azero, double a[], double b[] );
double *r8vec_sht ( int n, double a[]  );
double *r8vec_sqctb ( int n, double x[] );
double *r8vec_sqctf ( int n, double x[] );
double *r8vec_sqstb ( int n, double x[] );
double *r8vec_sqstf ( int n, double x[] );
double *r8vec_sst ( int n, double x[] );
double *r8vec_uniform_new ( int n, double b, double c, int *seed );

void timestamp ( );
