double arc_cosine ( double c );
double arc_sine ( double s );
double atan4 ( double y, double x );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4vec_copy ( int n, int a1[], int a2[] );
void icos_shape ( int point_num, int edge_num, int face_num, 
  int face_order_max, double point_coord[], int edge_point[], int face_order[],
  int face_point[] );
void icos_size ( int *point_num, int *edge_num, int *face_num, 
  int *face_order_max );
double r8_abs ( double x );
double r8_gamma ( double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int *seed );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_norm ( int n, double a[] );
void r8vec_polarize ( int n, double a[], double p[], double a_normal[], 
  double a_parallel[] );
double r8vec_sum ( int n, double a[] );
double sphere01_distance_xyz ( double xyz1[3], double xyz2[3] );
double sphere01_monomial_integral ( int e[3] );
double sphere01_quad_icos1c ( int factor, 
  void fun ( int n, double x[], double v[] ), int *node_num );
double sphere01_quad_icos1m ( int factor, 
  void fun ( int n, double x[], double v[] ), int *node_num );
double sphere01_quad_icos1v ( int factor, 
  void fun ( int n, double x[], double v[] ), int *node_num );
double sphere01_quad_icos2v ( int factor, 
  void fun ( int n, double x[], double v[] ), int *node_num );
double sphere01_quad_llc ( void f ( int n, double x[], double v[] ), double h, 
  int *n );
double sphere01_quad_llm ( void f ( int n, double x[], double v[] ), double h, 
  int *n );
double sphere01_quad_llv ( void f ( int n, double x[], double v[] ), double h, 
  int *n );
double sphere01_quad_mc ( void f ( int n, double x[], double v[] ), double h, 
  int *seed, int n );
int sphere01_quad_mc_size ( double h );
double *sphere01_sample ( int n, int *seed );
double sphere01_triangle_angles_to_area ( double a, double b, double c );
double *sphere01_triangle_project ( double a_xyz[3], double b_xyz[3], double c_xyz[3], 
  int f1, int f2, int f3 );
double *sphere01_triangle_project2 ( double a_xyz[3], double b_xyz[3], double c_xyz[3], 
  int f1, int f2, int f3 );
double *sphere01_triangle_sample ( int n, double v1[3], double v2[3], double v3[3], 
  int *seed );
void sphere01_triangle_sides_to_angles ( double as, double bs, double cs, 
  double *a, double *b, double *c );
void sphere01_triangle_vertices_to_angles ( double v1[3], double v2[3], 
  double v3[3], double *a, double *b, double *c );
double sphere01_triangle_vertices_to_area ( double v1[3], double v2[3], double v3[3] );
double *sphere01_triangle_vertices_to_centroid ( double v1[3], double v2[3], 
  double v3[3] );
void sphere01_triangle_vertices_to_midpoints ( double v1[3], double v2[3], double v3[3], 
  double m1[3], double m2[3], double m3[3] );
void sphere01_triangle_vertices_to_sides ( double v1[3], double v2[3], 
  double v3[3], double *as, double *bs, double *cs );
void timestamp ( );
double *tp_to_xyz ( double theta, double phi );
