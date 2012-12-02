double *cheby_t_zero ( int n );
double *cheby_u_zero ( int n );
void data_to_dif ( int ntab, double xtab[], double ytab[], double diftab[] );
double *data_to_dif_new ( int ntab, double xtab[], double ytab[] );
void data_to_dif_display ( int ntab, double xtab[], double ytab[], double diftab[] );
void data_to_r8poly ( int ntab, double xtab[], double ytab[], double c[] );
void dif_antideriv ( int ntab, double xtab[], double diftab[], int *ntab2,
  double xtab2[], double diftab2[] );
void dif_append ( int ntab, double xtab[], double diftab[], double xval, 
  double yval, int *ntab2, double xtab2[], double diftab2[] );
void dif_basis ( int ntab, double xtab[], double *diftab );
void dif_basis_i ( int ival, int ntab, double xtab[], double diftab[] );
void dif_deriv_table ( int nd, double xd[], double yd[], int *ndp,
  double xdp[], double ydp[] );
void dif_print ( int ntab, double xtab[], double diftab[], string title );
void dif_shift_x ( int nd, double xd[], double yd[], double xv );
void dif_shift_zero ( int nd, double xd[], double yd[] );
void dif_to_r8poly ( int ntab, double xtab[], double diftab[], double c[] );
double dif_val ( int ntab, double xtab[], double diftab[], double xval );
double *dif_vals ( int nd, double xd[], double yd[], int nv, double xv[] );
void nc_rule ( int norder, double a, double b, double xtab[], double weight[] );
double *lagrange_rule ( int n, double x[] );
double lagrange_sum ( int n, double x[], double w[], double y[], double xv );
double lagrange_val ( int n, double x[], double y[], double xv );
void ncc_rule ( int norder, double xtab[], double weight[] );
void nco_rule ( int norder, double xtab[], double weight[] );
void r8_swap ( double *x, double *y );
void r8poly_ant_cof ( int n, double poly_cof[], double poly_cof2[] );
double r8poly_ant_val ( int n, double poly_cof[], double xval );
void r8poly_basis ( int ntab, double xtab[], double r8poly_cof[] );
void r8poly_basis_1 ( int ival, int ntab, double xtab[], double r8poly_cof[] );
void r8poly_der_cof ( int n, double poly_cof[], double poly_cof2[] );
double r8poly_der_val ( int n, double poly_cof[], double xval );
int r8poly_order ( int na, double a[] );
void r8poly_print ( int n, double poly_cof[], string title );
void r8poly_shift ( double scale, double shift, int n, double poly_cof[] );
double r8poly_val_horner ( int n, double poly_cof[], double xval );
bool r8vec_distinct ( int n, double x[] );
void r8vec_indicator ( int n, double a[] );
void r8vec_print ( int n, double a[], string title );
void roots_to_dif ( int nroots, double roots[], int *ntab, double xtab[], 
  double diftab[] );
void roots_to_r8poly ( int nroots, double roots[], int *nc, double c[] );
void timestamp ( );
