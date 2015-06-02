int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4_swap ( int *i, int *j );
int i4_uniform_ab ( int a, int b, int *seed );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_sort_a ( int m, int n, int a[] );
void i4col_sort2_a ( int m, int n, int a[] );
int i4col_sorted_unique_count ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
void i4i4_sort_a ( int i1, int i2, int *j1, int *j2 );
void i4i4i4_sort_a ( int i1, int i2, int i3, int *j1, int *j2, int *j3 );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void i4vec_print ( int n, int a[], string title );
int i4vec_sum ( int n, int a[] );
void i4vec_zero ( int n, int a[] );
void mesh_base_one ( int node_num, int element_order, int element_num,
  int element_node[] );
void mesh_base_zero ( int node_num, int element_order, int element_num, 
  int element_node[] );
float r4_abs ( float x );
int r4_nint ( float x );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void r8_swap ( double *x, double *y );
double r8_uniform_01 ( int *seed );
double r8mat_det_4d ( double a[4*4] );
double *r8mat_mv_new ( int m, int n, double a[], double x[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
int r8mat_solve ( int n, int rhs_num, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
double *r8vec_cross_3d ( double v1[3], double v2[3] );
bool r8vec_is_nonnegative ( int n, double x[] );
bool r8vec_is_zero ( int n, double x[] );
double r8vec_length ( int dim_num, double x[] );
double r8vec_max ( int n, double dvec[] );
double r8vec_mean ( int n, double x[] );
double r8vec_min ( int n, double dvec[] );
void r8vec_print ( int n, double a[], string title );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int *seed );
double r8vec_variance ( int n, double x[] );
void r8vec_zero ( int n, double a[] );
int s_len_trim ( string s );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
int *tet_mesh_neighbor_tets ( int tetra_order, int tetra_num, 
  int tetra_node[] );
int *tet_mesh_node_order ( int tetra_order, int tetra_num, int tetra_node[], 
  int node_num );
//
//  Order 4
//
void tet_mesh_order4_adj_count ( int node_num, int tetra_num, 
  int tetra_node[], int *adj_num, int adj_row[] );
int *tet_mesh_order4_adj_set ( int node_num, int tetra_num, 
  int tetra_node[], int adj_num, int adj_row[] );
int tet_mesh_order4_boundary_face_count ( int tetra_num, int tetra_node[] );
int tet_mesh_order4_edge_count ( int tetra_num, int tetra_node[] );
void tet_mesh_order4_example_set ( int node_num, int tetra_num, 
  double node_xyz[], int tetra_node[] );
void tet_mesh_order4_example_size ( int *node_num, int *tetra_num );
void tet_mesh_order4_refine_compute ( int node_num1, int tetra_num1, 
  double node_xyz1[], int tetra_node1[], int node_num2, int tetra_num2,
   int edge_data[], double node_xyz2[], int tetra_node2[] );
void tet_mesh_order4_refine_size ( int node_num1, int tetra_num1, 
  int tetra_node1[], int *node_num2, int *tetra_num2, int edge_data[] );
void tet_mesh_order4_to_order10_compute ( int tetra_num, int tetra_node1[], 
  int node_num1, double node_xyz1[], int edge_data[], int tetra_node2[], 
  int node_num2, double node_xyz2[] );
void tet_mesh_order4_to_order10_size ( int tetra_num, int tetra_node1[], 
  int node_num1, int edge_data[], int *node_num2 );
//
//  Order 10
//
void tet_mesh_order10_adj_count ( int node_num, int tetra_num, 
  int tetra_node[], int *adj_num, int adj_row[] );
int *tet_mesh_order10_adj_set ( int node_num, int tetra_num, 
  int tetra_node[], int adj_num, int adj_row[] );
void tet_mesh_order10_example_set ( int node_num, int tetra_num, 
  double node_xyz[], int tetra_node[] );
void tet_mesh_order10_example_size ( int *node_num, int *tetra_num );
void tet_mesh_order10_to_order4_compute ( int tetra_num1, int tetra_node1[], 
  int tetra_num2, int tetra_node2[] );
void tet_mesh_order10_to_order4_size ( int node_num1, int tetra_num1,
  int *node_num2, int *tetra_num2 );
void tet_mesh_quad ( int node_num, double node_xyz[], int tet_order, 
  int tet_num, int tet_node[], 
  void f ( int n, double xyz_vec[], double fvec[] ), 
  int quad_num, double quad_xyz[], double quad_w[], double *quad_value, 
  double *region_volume );
void tet_mesh_quality1 ( int node_num, double node_xyz[], 
  int tetra_order, int tetra_num, int tetra_node[], double *value_min, 
  double *value_mean, double *value_max, double *value_var );
void tet_mesh_quality2 ( int node_num, double node_xyz[], int tetra_order, 
  int tetra_num, int tetra_node[], double *value_min, double *value_mean, 
  double *value_max, double *value_var );
void tet_mesh_quality3 ( int node_num, double node_xyz[], int tetra_order, 
  int tetra_num, int tetra_node[], double *value_min, double *value_mean, 
  double *value_max, double *value_var );
void tet_mesh_quality4 ( int node_num, double node_xyz[], int tetra_order, 
  int tetra_num, int tetra_node[], double *value_min, double *value_mean, 
  double *value_max, double *value_var );
void tet_mesh_quality5 ( int node_num, double node_xyz[], int tetra_order, 
  int tetra_num, int tetra_node[], double *value_min, double *value_mean, 
  double *value_max, double *value_var );
int tet_mesh_search_delaunay ( int node_num, double node_xyz[], int tet_order, 
  int tet_num, int tet_node[], int tet_neighbor[], double p[], int *face, 
  int *step_num );
int tet_mesh_search_naive ( int node_num, double node_xyz[],
  int tet_order, int tet_num, int tet_node[], double p[], int *step_num );
double *tetrahedron_barycentric ( double tet_xyz[], double p[] );
void tetrahedron_circumsphere_3d ( double tetra[3*4], double *r, double pc[3] );
double *tetrahedron_edge_length_3d ( double tetra[3*4] );
void tetrahedron_insphere_3d ( double tetra[3*4], double *r, double pc[3] );
void tetrahedron_order4_physical_to_reference ( double t[], int n, 
  double phy[], double ref[] );
void tetrahedron_order4_reference_to_physical ( double t[], int n, 
  double ref[], double phy[] );
double tetrahedron_quality1_3d ( double tetra[3*4] );
double tetrahedron_quality2_3d ( double tetra[3*4] );
double tetrahedron_quality3_3d ( double tetra[3*4] );
double tetrahedron_quality4_3d ( double tetra[3*4] );
void tetrahedron_reference_sample ( int n, int *seed, double p[] );
void tetrahedron_sample ( double tetra[3*4], int n, int *seed, double p[] );
double tetrahedron_volume ( double tetra[3*4] );
void timestamp ( );
