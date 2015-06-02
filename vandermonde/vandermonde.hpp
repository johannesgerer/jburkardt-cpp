double *bivand1 ( int n, double alpha[], double beta[] );
double *bivand2 ( int n, double alpha[], double beta[] );
double *dvand ( int n, double alpha[], double b[] );
void dvandprg ( int n, double alpha[], double b[], double x[], double c[], 
  double m[] );
double *pvand ( int n, double alpha[], double b[] );
void pvandprg ( int n, double alpha[], double b[], double x[], double c[], 
  double m[] );
double *r8mat_mtv_new ( int m, int n, double a[], double x[] );
double *r8mat_mv_new ( int m, int n, double a[], double x[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8vec_copy_new ( int n, double a1[] );
void r8vec_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );
void timestamp ( );
double *vand1 ( int n, double x[] );
