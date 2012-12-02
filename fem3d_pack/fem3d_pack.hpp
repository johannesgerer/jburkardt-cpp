void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m );
void bandwidth_var ( int element_order, int element_num, int element_node[],
  int node_num, int var_node[], int var_num, int var[], int *ml, int *mu, 
  int *m );
double *basis_brick8 ( int n, double p[] );
void basis_brick8_test ( );
double *basis_brick20 ( int n, double p[] );
void basis_brick20_test ( );
double *basis_brick27 ( int n, double p[] );
void basis_brick27_test ( );
void basis_mn_tet4 ( double t[3*4], int n, double p[], double phi[] );
void basis_mn_tet4_test ( );
void basis_mn_tet10 ( double t[3*4], int n, double p[], double phi[] );
void basis_mn_tet10_test ( );
int i4_max ( int i1, int i2 );
double *nodes_brick8 ( );
double *nodes_brick20 ( );
double *nodes_brick27 ( );
double *physical_to_reference_tet4 ( double t[], int n, double phy[] );
double r8_abs ( double x );
void r8ge_fss ( int n, double a[], int nb, double b[] );
double *r8mat_copy_new ( int m, int n, double a1[] );
double r8mat_det_4d ( double a[4*4] );
int r8mat_solve ( int n, int rhs_num, double a[] );
double *r8mat_mv ( int m, int n, double a[], double x[] );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int *seed );
double *reference_tet4_sample ( int n, int *seed );
double *reference_tet4_uniform ( int n, int *seed );
double *reference_tet4_uniform2 ( int n, int *seed );
double *reference_to_physical_tet4 ( double t[], int n, double ref[] );
double *tetrahedron_barycentric ( double tetra[3*4], double p[3] );
double tetrahedron_volume ( double tetra[3*4] );
void timestamp ( );
