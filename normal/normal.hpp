complex <float> c4_normal_01 ( int &seed );

complex <double> c8_normal_01 ( int &seed );

int i4_normal_ab ( float a, float b, int &seed );

long long int i8_normal_ab ( double a, double b, long long int &seed );

float r4_normal_01 ( int &seed );
float r4_normal_ab ( float a, float b, int &seed );
float r4_uniform_01 ( int &seed );
void r4mat_print ( int m, int n, float a[], string title );
void r4mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r4vec_print ( int n, float a[], string title );
float *r4vec_uniform_01_new ( int n, int &seed );

double r8_normal_01 ( int &seed );
double r8_normal_ab ( double a, double b, int &seed );
double r8_uniform_01 ( int &seed );
void r8mat_normal_01 ( int m, int n, int &seed, double x[] );
double *r8mat_normal_01_new ( int m, int n, int &seed );
void r8mat_normal_ab ( int m, int n, double a, double b, int &seed, double x[] );
double *r8mat_normal_ab_new ( int m, int n, double a, double b, int &seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8vec_normal_01 ( int n, int &seed, double x[] );
double *r8vec_normal_01_new ( int n, int &seed );
void r8vec_normal_ab ( int n, double a, double b, int &seed, double x[] );
double *r8vec_normal_ab_new ( int n, double a, double b, int &seed );
void r8vec_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );

void timestamp ( );

