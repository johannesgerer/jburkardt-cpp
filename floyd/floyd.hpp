double cpu_time ( );
int i4_huge ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4mat_floyd ( int n, int a[] );
void i4mat_print ( int m, int n, int a[], string title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double r8_huge ( );
double r8_min ( double x, double y );
void r8mat_floyd ( int n, double a[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double r8vec_diff_norm ( int n, double a[], double b[] );
void timestamp ( );
