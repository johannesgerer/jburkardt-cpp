void cartesian_to_hypersphere ( int m, int n, double c[], double x[], 
  double r[], double theta[] );
double hypersphere_01_area ( int m );
void hypersphere_01_area_values ( int &n_data, int &m, double &area );
double *hypersphere_01_interior_uniform ( int m, int n, int &seed );
double *hypersphere_01_surface_uniform ( int m, int n, int &seed );
double hypersphere_01_volume ( int m );
void hypersphere_01_volume_values ( int &n_data, int &m, double &volume );
double hypersphere_area ( int m, double r );
double *hypersphere_stereograph ( int m, int n, double x[] );
double *hypersphere_stereograph_inverse ( int m, int n, double x2[] );
double *hypersphere_surface_uniform ( int m, int n, double r, double c[], 
  int &seed );
double *hypersphere_to_cartesian ( int m, int n, double c[], double r[], 
  double theta[] );
double hypersphere_volume ( int m, double r );
int i4_min ( int i1, int i2 );
double r8_uniform_01 ( int &seed );
double r8mat_norm_fro_affine ( int m, int n, double a1[], double a2[] );
double *r8mat_normal_01_new ( int m, int n, int &seed );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
double *r8vec_normal_01_new ( int n, int &seed );
void r8vec_transpose_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );
double *sphere_stereograph ( int m, int n, double p[] );
double *sphere_stereograph_inverse ( int m, int n, double q[] );
void timestamp ( );
