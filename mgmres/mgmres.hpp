void atx_cr ( int n, int nz_num, double a[], int ia[], int ja[], double x[], 
  double w[] );
void atx_st ( int n, int nz_num, double a[], int ia[], int ja[], double x[], 
  double w[] );
void ax_cr ( int n, int nz_num, double a[], int ia[], int ja[], double x[], 
  double w[] );
void ax_st ( int n, int nz_num, double a[], int ia[], int ja[], double x[], 
  double w[] );
void diagonal_pointer_cr ( int n, int nz_num, int ia[], int ja[], int ua[] );
void lus_cr ( int n, int nz_num, int ia[], int ja[], double l[], int ua[], 
  double r[], double z[] );
void mgmres_st ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
  double rhs[], int itr_max, int mr, double tol_abs, double tol_rel );
void mult_givens ( double c, double s, int k, double g[] );
void pmgmres_ilu_cr ( int n, int nz_num, int ia[], int ja[], double a[], 
  double x[], double rhs[], int itr_max, int mr, double tol_abs, 
  double tol_rel );
double r8vec_dot ( int n, double a1[], double a2[] );
double *r8vec_uniform_01 ( int n, int *seed );
void rearrange_cr ( int n, int nz_num, int ia[], int ja[], double a[] );
void timestamp ( );
