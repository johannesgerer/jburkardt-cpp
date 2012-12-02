int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double *r83_cr_fa ( int n, double a[] );
double *r83_cr_sl ( int n, double a_cr[], double b[] );
double *r83_cr_sls ( int n, double a_cr[], int nb, double b[] );
void r83_gs_sl ( int n, double a[], double b[], double x[], int it_max, 
  int job );
double *r83_mxv_new ( int n, double a[], double x[] );
void r83_print ( int n, double a[], string title );
void r83_print_some ( int n, double a[], int ilo, int jlo, int ihi, int jhi,
  string title );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
double *r8vec_indicator_new ( int n );
void r8vec_print ( int n, double a[], string title );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title );
void timestamp ( );
