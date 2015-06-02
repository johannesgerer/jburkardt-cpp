void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m );
void bandwidth_var ( int element_order, int element_num, int element_node[],
  int node_num, int var_node[], int var_num, int var[], int *ml, int *mu, 
  int *m );
void basis_11_t3 ( double t[2*3], int i, double p[2], double *qi, double *dqidx, 
  double *dqidy );
void basis_11_t3_test ( );
void basis_11_t4 ( double t[2*3], int i, double p[2], double *qi, double *dqidx, 
  double *dqidy );
void basis_11_t4_test ( );
void basis_11_t6 ( double t[2*6], int i, double p[2], double *qi, double *dqidx, 
  double *dqidy );
void basis_11_t6_test ( );
void basis_mn_q4 ( double q[2*4], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] );
void basis_mn_q4_test ( );
void basis_mn_t3 ( double t[2*3], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] );
void basis_mn_t3_test ( );
void basis_mn_t4 ( double t[2*4], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] );
void basis_mn_t4_test ( );
void basis_mn_t6 ( double t[2*6], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] );
void basis_mn_t6_test ( );
char ch_cap ( char c );
double degrees_to_radians ( double angle );
void derivative_average_t3 ( int node_num, double node_xy[], int element_num,
  int element_node[], double c[], double dcdx[], double dcdy[] );
void div_q4 ( int m, int n, double u[], double v[], double xlo, double xhi, 
  double ylo, double yhi, double div[], double vort[] );
string element_code ( int i );
void elements_eps ( string file_name, int node_num, double node_xy[], string code, 
  int element_num, bool element_mask[], int element_node[], int node_show, 
  int element_show );
int *grid_element ( string code, int element_order, int nelemx, int nelemy );
int grid_element_num ( string code, int nelemx, int nelemy );
int grid_node_num ( string code, int nelemx, int nelemy );
double *grid_nodes_01 ( int x_num, int y_num );
void grid_print ( int n, int element_num, int element_node[] );
int *grid_q4_element ( int nelemx, int nelemy );
int grid_q4_element_num ( int nelemx, int nelemy );
int grid_q4_node_num ( int nelemx, int nelemy );
int *grid_q8_element ( int nelemx, int nelemy );
int grid_q8_element_num ( int nelemx, int nelemy );
int grid_q8_node_num ( int nelemx, int nelemy );
int *grid_q9_element ( int nelemx, int nelemy );
int grid_q9_element_num ( int nelemx, int nelemy );
int grid_q9_node_num ( int nelemx, int nelemy );
int *grid_q12_element ( int nelemx, int nelemy );
int grid_q12_element_num ( int nelemx, int nelemy );
int grid_q12_node_num ( int nelemx, int nelemy );
int *grid_q16_element ( int nelemx, int nelemy );
int grid_q16_element_num ( int nelemx, int nelemy );
int grid_q16_node_num ( int nelemx, int nelemy );
int *grid_ql_element ( int nelemx, int nelemy );
int grid_ql_element_num ( int nelemx, int nelemy );
int grid_ql_node_num ( int nelemx, int nelemy );
void grid_shape_2d ( int n, double a[], int *n1, int *n2 );
int *grid_t3_element ( int nelemx, int nelemy );
int grid_t3_element_num ( int nelemx, int nelemy );
int grid_t3_node_num ( int nelemx, int nelemy );
int *grid_t4_element ( int nelemx, int nelemy );
int grid_t4_element_num ( int nelemx, int nelemy );
int grid_t4_node_num ( int nelemx, int nelemy );
int *grid_t6_element ( int nelemx, int nelemy );
int grid_t6_element_num ( int nelemx, int nelemy );
int grid_t6_node_num ( int nelemx, int nelemy );
int *grid_t10_element ( int nelemx, int nelemy );
int grid_t10_element_num ( int nelemx, int nelemy );
int grid_t10_node_num ( int nelemx, int nelemy );
void grid_test ( string code );
int grid_width ( int n, int element_num, int element_node[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void i4mat_write ( string output_filename, int m, int n, int table[] );
void i4vec_print ( int n, int a[], string title );;
void interp ( string code, int n, double r, double s, double ubase[], double *u, 
  double *dudr, double *duds );
void interp_test ( string code );
void legendre_com ( int norder, double xtab[], double weight[] );
void legendre_set ( int norder, double xtab[], double weight[] );
double *map ( string code, int n );
void map_test ( string code );
double *mass_matrix_t3 ( int node_num, int element_num, int element_node[], 
  double node_xy[] );
double *mass_matrix_t6 ( int node_num, int element_num, int element_node[], 
  double node_xy[] );
int next_boundary_node ( int node, string code );
int next_boundary_node_q4 ( int node );
int next_boundary_node_q8 ( int node );
int next_boundary_node_q9 ( int node );
int next_boundary_node_q12 ( int node );
int next_boundary_node_q16 ( int node );
int next_boundary_node_ql ( int node );
int next_boundary_node_t3 ( int node );
int next_boundary_node_t4 ( int node );
int next_boundary_node_t6 ( int node );
int next_boundary_node_t10 ( int node );
void node_reference ( string code, double r[], double s[], double *area );
void node_reference_q4 ( double r[4], double s[4], double *area );
void node_reference_q8 ( double r[8], double s[8], double *area );
void node_reference_q9 ( double r[9], double s[9], double *area );
void node_reference_q12 ( double r[12], double s[12], double *area );
void node_reference_q16 ( double r[16], double s[16], double *area );
void node_reference_ql ( double r[6], double s[6], double *area );
void node_reference_t3 ( double r[3], double s[3], double *area );
void node_reference_t4 ( double r[4], double s[4], double *area );
void node_reference_t6 ( double r[6], double s[6], double *area );
void node_reference_t10 ( double r[10], double s[10], double *area );
int ns_t6_var_count ( int element_num, int element_node[], int node_num, 
  int var_node[] );
int *ns_t6_var_set ( int element_num, int element_node[], int node_num, 
  int var_node[], int var_num );
int order_code ( string code );
void physical_to_reference_t3 ( double t[], int n, double phy[], double ref[] );
void points_plot ( string file_name, int node_num, double node_xy[], 
  bool node_label );
void poly ( string code, int rexp[], int sexp[] );
void poly_q4 ( int rexp[], int sexp[4] );
void poly_q8 ( int rexp[48], int sexp[8] );
void poly_q9 ( int rexp[9], int sexp[9] );
void poly_q12 ( int rexp[12], int sexp[12] );
void poly_q16 ( int rexp[16], int sexp[16] );
void poly_ql ( int rexp[6], int sexp[6] );
void poly_t3 ( int rexp[3], int sexp[3] );
void poly_t6 ( int rexp[6], int sexp[6] );
void poly_t10 ( int rexp[10], int sexp[10] );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_huge ( );
int r8_nint ( double x );
double r8_power ( double r, int p );
double r8_uniform_01 ( int *seed );
int r8ge_fa ( int n, double a[], int pivot[] );
double *r8ge_inverse ( int n, double a[], int pivot[] );
void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void reference_sample ( string code, int *seed, double *r, double *s );
void reference_to_physical_q4 ( double q4[2*4], int n, double rs[], 
  double xy[] );
void reference_to_physical_t3 ( double t[], int n, double ref[], double phy[] );
void reference_to_physical_t6 ( double t[], int n, double ref[], double phy[] );
bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
double serene ( string type, double ve, double vn, double vne, double vnw, 
  double vs, double vse, double vsw, double vw );
void shape ( string code, double r, double s, double t[], double dtdr[], 
  double dtds[] );
void shape_q4 ( double r, double s, double t[4], double dtdr[4], 
  double dtds[4] );
void shape_q8 ( double r, double s, double t[8], double dtdr[8], 
  double dtds[8] );
void shape_q9 ( double r, double s, double t[9], double dtdr[9], 
  double dtds[9] );
void shape_q12 ( double r, double s, double t[12], double dtdr[12], 
  double dtds[12] );
void shape_q16 ( double r, double s, double t[16], double dtdr[16], 
  double dtds[16] );
void shape_ql ( double r, double s, double t[6], double dtdr[6], 
  double dtds[6] );
void shape_t3 ( double r, double s, double t[3], double dtdr[3], 
  double dtds[3] );
void shape_t4 ( double r, double s, double t[4], double dtdr[4], 
  double dtds[4] );
void shape_t6 ( double r, double s, double t[6], double dtdr[6], 
  double dtds[6] );
void shape_t10 ( double r, double s, double t[10], double dtdr[10], 
  double dtds[10] );
void shape_test ( string code );
int sphere_grid_element_num ( string code, int nelemx, int nelemy );
int sphere_grid_node_num ( string code, int nelemx, int nelemy );
int *sphere_grid_q4_element ( int nelemx, int nelemy );
int sphere_grid_q4_element_num ( int nelemx, int nelemy );
int sphere_grid_q4_node_num ( int nelemx, int nelemy );
double *sphere_grid_q4_node_xyz ( int nelemx, int nelemy );
int *sphere_grid_q9_element ( int nelemx, int nelemy );
int sphere_grid_q9_element_num ( int nelemx, int nelemy );
int sphere_grid_q9_node_num ( int nelemx, int nelemy );
double *sphere_grid_q9_node_xyz ( int nelemx, int nelemy );
int *sphere_grid_q16_element ( int nelemx, int nelemy );
int sphere_grid_q16_element_num ( int nelemx, int nelemy );
int sphere_grid_q16_node_num ( int nelemx, int nelemy );
double *sphere_grid_q16_node_xyz ( int nelemx, int nelemy );
int *sphere_grid_t3_element ( int nelemx, int nelemy );
int sphere_grid_t3_element_num ( int nelemx, int nelemy );
int sphere_grid_t3_node_num ( int nelemx, int nelemy );
double *sphere_grid_t3_node_xyz ( int nelemx, int nelemy );
int *sphere_grid_t6_element ( int nelemx, int nelemy );
int sphere_grid_t6_element_num ( int nelemx, int nelemy );
int sphere_grid_t6_node_num ( int nelemx, int nelemy );
double *sphere_grid_t6_node_xyz ( int nelemx, int nelemy );
void timestamp ( );
void triangle_unit_set ( int rule, double rtab[], double stab[], 
  double weight[] );
int triangle_unit_size ( int rule );

