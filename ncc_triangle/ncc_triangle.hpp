void file_name_inc ( char *file_name );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
int ncc_triangle_degree ( int rule );
int ncc_triangle_order_num ( int rule );
void ncc_triangle_rule ( int rule, int order_num, double xy[], double w[] );
int ncc_triangle_rule_num ( void );
int *ncc_triangle_suborder ( int rule, int suborder_num );
int ncc_triangle_suborder_num ( int rule );
void ncc_triangle_subrule ( int rule, int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
void ncc_triangle_subrule_01 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_triangle_subrule_02 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_triangle_subrule_03 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_triangle_subrule_04 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_triangle_subrule_05 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_triangle_subrule_06 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_triangle_subrule_07 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_triangle_subrule_08 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_triangle_subrule_09 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
double r8_huge ( void );
int r8_nint ( double x );
void reference_to_physical_t3 ( double t[], int n, double ref[], double phy[] );
int s_len_trim ( char *s );
void timestamp ( void );
double triangle_area ( double t[2*3] );
void triangle_points_plot ( char *file_name, double node_xy[], int node_show, 
  int point_num, double point_xy[], int point_show );
