void file_name_inc ( char *file_name );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
double r8_huge ( );
int r8_nint ( double x );
void reference_to_physical_t3 ( double t[], int n, double ref[], double phy[] );
int s_len_trim ( char *s );
void timestamp ( );
double triangle_area ( double t[2*3] );
int triangle_nco_degree ( int rule );
int triangle_nco_order_num ( int rule );
void triangle_nco_rule ( int rule, int order_num, double xy[], double w[] );
int triangle_nco_rule_num ( );
int *triangle_nco_suborder ( int rule, int suborder_num );
int triangle_nco_suborder_num ( int rule );
void triangle_nco_subrule ( int rule, int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
void triangle_nco_subrule_01 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void triangle_nco_subrule_02 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void triangle_nco_subrule_03 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void triangle_nco_subrule_04 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void triangle_nco_subrule_05 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void triangle_nco_subrule_06 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void triangle_nco_subrule_07 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void triangle_nco_subrule_08 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void triangle_nco_subrule_09 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
