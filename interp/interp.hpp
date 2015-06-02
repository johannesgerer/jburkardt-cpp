double *cc_abscissas ( int n );
double *cc_abscissas_ab ( double a, double b, int n );
double *f1_abscissas ( int n );
double *f1_abscissas_ab ( double a, double b, int n );
double *f2_abscissas ( int n );
double *f2_abscissas_ab ( double a, double b, int n );
double *interp_lagrange ( int m, int data_num, double t_data[], 
  double p_data[], int interp_num, double t_interp[] );
double *interp_linear ( int m, int data_num, double t_data[], double p_data[], 
  int interp_num, double t_interp[] );
double *interp_nearest ( int m, int data_num, double t_data[], double p_data[], 
  int interp_num, double t_interp[] );
double *lagrange_value ( int data_num, double t_data[], int interp_num, 
  double t_interp[] );
double *ncc_abscissas ( int n );
double *ncc_abscissas_ab ( double a, double b, int n );
double *nco_abscissas ( int n );
double *nco_abscissas_ab ( double a, double b, int n );
double *parameterize_arc_length ( int m, int data_num, double p_data[] );
double *parameterize_index ( int m, int data_num, double p_data[] );
double r8_abs ( double x );
double *r8mat_expand_linear2 ( int m, int n, double a[], int m2, int n2 );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
bool r8vec_ascends_strictly ( int n, double x[] );
void r8vec_bracket0 ( int n, double x[], double xval, int &left,
  int &right );
double *r8vec_expand_linear ( int n, double x[], int fat );
double *r8vec_expand_linear2 ( int n, double x[], int before, int fat, 
  int after );
int r8vec_sorted_nearest0 ( int n, double a[], double value );
void timestamp ( );

