char ch_cap ( char ch );
int ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void ffmsh_2d_data_example ( int v_num, int e_num, int t_num, double v_xy[], 
  int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] );
void ffmsh_2d_data_print ( string title, int v_num, int e_num, int t_num, 
  double v_xy[], int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] );
void ffmsh_2d_data_read ( string ffmsh_filename, int v_num, int e_num, int t_num, 
  double v_xy[], int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] );
void ffmsh_2d_size_example ( int &v_num, int &e_num, int &t_num );
void ffmsh_2d_size_print ( string title, int v_num, int e_num, int t_num );
void ffmsh_2d_size_read ( string ffmsh_filename, int &v_num, int &e_num, 
  int &t_num );
void ffmsh_2d_write ( string ffmsh_filename, int v_num, int e_num, int t_num, 
  double v_xy[], int v_l[], int e_v[], int e_l[], int t_v[], int t_l[] );
void i4mat_copy ( int m, int n, int a1[], int a2[] );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void i4vec_copy ( int n, int a1[], int a2[] );
void i4vec_print ( int n, int a[], string title );
void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
int s_len_trim ( string s );
int s_to_i4 ( string s, int &last, bool &error );
double s_to_r8 ( string s, int &lchar, bool &error );
void timestamp ( );
