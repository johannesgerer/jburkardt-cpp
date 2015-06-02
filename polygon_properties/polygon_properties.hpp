double *angle_half ( double p1[2], double p2[2], double p3[2] );
double angle_rad ( double p1[2], double p2[2], double p3[2] );
bool between ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
bool collinear ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
bool diagonal ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
bool diagonalie ( int im1, int ip1, int n, int next[], double x[], double y[] );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
bool in_cone ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
bool intersect ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
bool intersect_prop ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
bool l4_xor ( bool l1, bool l2 );
double *polygon_angles ( int n, double v[] );
double polygon_area ( int n, double v[] );
double polygon_area_2 ( int n, double v[] );
double *polygon_centroid ( int n, double v[] );
double *polygon_centroid_2 ( int n, double v[] );
bool polygon_contains_point ( int n, double v[], double p[2] );
bool polygon_contains_point_2 ( int n, double v[], double p[2] );
double polygon_diameter ( int n, double v[] );
double *polygon_expand ( int n, double v[], double h );
void polygon_inrad_data ( int n, double radin, double &area, double &radout, 
  double &side );
double polygon_integral_1 ( int n, double v[] );
double polygon_integral_x ( int n, double v[] );
double polygon_integral_xx ( int n, double v[] );
double polygon_integral_xy ( int n, double v[] );
double polygon_integral_y ( int n, double v[] );
double polygon_integral_yy ( int n, double v[] );
int polygon_is_convex ( int n, double v[] );
double polygon_lattice_area ( int i, int b );
void polygon_outrad_data ( int n, double radout, double &area, double &radin, 
  double &side );
double polygon_point_dist ( int n, double v[], double p[] );
double *polygon_point_near ( int n, double v[], double p[] );
double *polygon_sample ( int nv, double v[], int n, int &seed );
void polygon_side_data ( int n, double side, double &area, double &radin, 
  double &radout );
int *polygon_triangulate ( int n, double x[], double y[] );
double r8_degrees ( double radians );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int &seed );
int r8mat_solve ( int n, int rhs_num, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
double r8vec_norm ( int n, double a[] );
void r8vec_print ( int n, double a[], string title );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int &seed );
double segment_point_dist ( double p1[2], double p2[2], double p[2] );
double *segment_point_near ( double p1[2], double p2[2], double p[2] );
void timestamp ( );
double triangle_area ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
double *triangle_barycentric ( double t[2*3], double p[2] );
bool triangle_contains_point_1 ( double t[2*3], double p[2] );

