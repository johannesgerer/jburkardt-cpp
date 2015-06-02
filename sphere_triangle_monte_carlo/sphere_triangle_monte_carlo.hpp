int i4_min ( int i1, int i2 );
double *monomial_value ( int m, int n, int e[], double x[] );
double r8_acos ( double c );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int &seed );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_norm ( int n, double a[] );
void r8vec_normalize ( int n, double a[] );
double r8vec_sum ( int n, double a[] );
void r8vec_transpose_print ( int n, double a[], string title );
double *sphere01_sample ( int n, int &seed );
double sphere01_triangle_angles_to_area ( double a, double b, double c );
void sphere01_triangle_sides_to_angles ( double as, double bs, double cs, 
  double &a, double &b, double &c );
double *sphere01_triangle_sample ( int n, double v1[3], double v2[3], 
  double v3[3], int &seed );
double sphere01_triangle_vertices_to_area ( double v1[3], double v2[3], 
  double v3[3] );
void sphere01_triangle_vertices_to_sides ( double v1[3], double v2[3], 
  double v3[3], double &as, double &bs, double &cs );
void timestamp ( );

