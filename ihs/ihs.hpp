void covariance ( int m, int n, int x[], double &average, double &std,
  double &covc );
int get_seed ( );
int i4_log_10 ( int i );
int i4_uniform_ab ( int b, int c, int &seed );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title );
int *ihs ( int m, int n, int d, int &seed );
double r8_uniform_01 ( int &seed );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double r8vec_average ( int n, double a[] );
double r8vec_std ( int n, double a[] );
void timestamp ( );
