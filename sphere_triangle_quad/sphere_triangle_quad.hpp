double arc_cosine ( double c );
double arc_sine ( double s );
double atan4 ( double y, double x );
int i4_power ( int i, int j );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int *seed );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_norm ( int n, double a[] );
void r8vec_polarize ( int n, double a[], double p[], double a_normal[], 
  double a_parallel[] );
double sphere01_distance_xyz ( double xyz1[3], double xyz2[3] );
double *sphere01_sample ( int n, int *seed );
double sphere01_triangle_angles_to_area ( double a, double b, double c );
double *sphere01_triangle_project ( double a_xyz[3], double b_xyz[3], 
  double c_xyz[3], int f1, int f2, int f3 );
double *sphere01_triangle_project2 ( double a_xyz[3], double b_xyz[3], 
  double c_xyz[3], int f1, int f2, int f3 );
double sphere01_triangle_quad_00 ( int n, double v1[3], double v2[3], 
  double v3[3], double f ( double x[] ), int *seed );
double sphere01_triangle_quad_01 ( double v1[3], double v2[3], double v3[3], 
  double f ( double x[] ) );
double sphere01_triangle_quad_02 ( double v1[3], double v2[3], double v3[3], 
  double f ( double x[] ) );
double sphere01_triangle_quad_03 ( double v1[3], double v2[3], double v3[3], 
  double f ( double x[] ) );
double sphere01_triangle_quad_icos1c ( double a_xyz[3], double b_xyz[3],
  double c_xyz[], int factor, double fun ( double x[] ), int *node_num );
double sphere01_triangle_quad_icos1m ( double a_xyz[3], double b_xyz[3],
  double c_xyz[], int factor, double fun ( double x[] ), int *node_num );
double sphere01_triangle_quad_icos1v ( double a_xyz[3], double b_xyz[3],
  double c_xyz[], int factor, double fun ( double x[] ), int *node_num );
double sphere01_triangle_quad_icos2v ( double a_xyz[3], double b_xyz[3],
  double c_xyz[], int factor, double fun ( double x[] ), int *node_num );
double *sphere01_triangle_sample ( int n, double v1[3], double v2[3], 
  double v3[3], int *seed );
void sphere01_triangle_sides_to_angles ( double as, double bs, double cs, 
  double *a, double *b, double *c );
double sphere01_triangle_vertices_to_area ( double v1[3], double v2[3], 
  double v3[3] );
double *sphere01_triangle_vertices_to_centroid ( double v1[3], double v2[3], 
  double v3[3] );
void sphere01_triangle_vertices_to_midpoints ( double v1[3], double v2[3], 
  double v3[3], double m1[3], double m2[3], double m3[3] );
void sphere01_triangle_vertices_to_sides ( double v1[3], double v2[3], 
  double v3[3], double *as, double *bs, double *cs );
void timestamp ( );
