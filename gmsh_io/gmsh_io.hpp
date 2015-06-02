char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void gmsh_data_read ( string gmsh_filename, int node_dim, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] );
void gmsh_size_read ( string gmsh_filename, int &node_num, int &node_dim,
  int &element_num, int &element_order );
int *gmsh_mesh2d_element_data_example ( int element_num, int element_order );
void gmsh_mesh2d_element_size_example ( int &element_num, int &element_order );
double *gmsh_mesh2d_node_data_example ( int node_num, int node_dim );
void gmsh_mesh2d_node_size_example ( int &node_num, int &node_dim );
void gmsh_mesh1d_write ( string gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_mesh2d_write ( string gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
void gmsh_mesh3d_write ( string gmsh_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );
int *i4mat_copy_new ( int m, int n, int a1[] );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double *r8mat_copy_new ( int m, int n, double a1[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
bool s_begin ( string s1, string s2 );
int s_len_trim ( string s );
int s_to_i4 ( string s, int &last, bool &error );
double s_to_r8 ( string s, int &lchar, bool &error );
void timestamp ( );
