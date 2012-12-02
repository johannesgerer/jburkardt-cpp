double *burgers_solution ( double nu, int vxn, double vx[], int vtn, 
  double vt[] );
void hermite_ek_compute ( int n, double x[], double w[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void imtqlx ( int n, double d[], double e[], double z[] );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_sign ( double x );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double *r8vec_even_new ( int n, double alo, double ahi );
void r8vec_print ( int n, double a[], string title );
void timestamp ( );
