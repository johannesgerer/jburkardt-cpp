int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double r8_uniform_01 ( int &seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8pp_print ( int n, double a[], string title );
void r8pp_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
void r8utp_print ( int n, double a[], string title );
void r8utp_print_some ( int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
void rnorm ( int &seed, double &u1, double &u2 );
void timestamp ( );
double *wshrt ( double d[], int n, int np, int &seed );

