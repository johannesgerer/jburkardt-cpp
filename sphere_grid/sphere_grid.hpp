double arc_cosine ( double c );
double arc_sine ( double s );
double atan4 ( double y, double x );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
string i4_to_string ( int i4, string format );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void i4vec_copy ( int n, int a1[], int a2[] );
void icos_shape ( int point_num, int edge_num, int face_num, 
  int face_order_max, double point_coord[], int edge_point[], int face_order[],
  int face_point[] );
void icos_num ( int *point_num, int *edge_num, int *face_num, 
  int *face_order_max );
double r8_abs ( double x );
double r8_modp ( double x, double y );
double r8_uniform_01 ( int *seed );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_diff_norm ( int n, double a[], double b[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
bool r8vec_eq ( int n, double a1[], double a2[] );
double r8vec_norm ( int n, double a[] );
void r8vec_polarize ( int n, double a[], double p[], double a_normal[], 
  double a_parallel[] );
void r8vec_print ( int n, double a[], string title );
double *sphere_cubed_ijk_to_xyz_old ( int n, int i, int j, int k );
void sphere_cubed_ijk_to_xyz ( int n, int i, int j, int k, double xyz[3] );
int sphere_cubed_line_num ( int n );
double *sphere_cubed_lines ( int n, int line_num );
double *sphere_cubed_points ( int n, int ns );
void sphere_cubed_points_face ( int n, int i1, int j1, int k1, int i2, int j2, 
  int k2, int &ns, double xyz[] );
int sphere_cubed_point_num ( int n );
double sphere_distance_xyz ( double xyz1[3], double xyz2[3] );
int *sphere_grid_q4 ( int lat_num, int long_num );
int *sphere_grid_t3 ( int lat_num, int long_num );
int sphere_icos_edge_num ( int factor );
int sphere_icos_face_num ( int factor );
int sphere_icos_point_num ( int factor );
double *sphere_icos1_points ( int factor, int node_num );
double *sphere_icos2_points ( int factor, int node_num );
int sphere_line_project ( double r, double pc[3], int n, double p[], 
  int maxpnt2, double pp[], double thetamin, double thetamax );
int *sphere_ll_lines ( int nlat, int nlong, int line_num );
int sphere_ll_line_num ( int nlat, int nlong );
double *sphere_ll_points ( double r, double pc[3], int nlat, int nlong, int point_num );
int sphere_ll_point_num ( int nlat, int nlong );
int *sphere_llq_lines ( int nlat, int nlong, int line_num );
int sphere_llq_line_num ( int nlat, int nlong );
double *sphere_spiralpoints ( double r, double pc[3], int n );
double *sphere_unit_sample ( int n, int *seed );
void timestamp ( );
