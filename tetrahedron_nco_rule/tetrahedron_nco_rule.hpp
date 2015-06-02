double r8mat_det_4d ( double a[] );
void reference_to_physical_t4 ( double t[], int n, double ref[], double phy[] );
int tetrahedron_nco_degree ( int rule );
int tetrahedron_nco_order_num ( int rule );
void tetrahedron_nco_rule ( int rule, int order_num, double xyz[], double w[] );
int tetrahedron_nco_rule_num ( );
int *tetrahedron_nco_suborder ( int rule, int suborder_num );
int tetrahedron_nco_suborder_num ( int rule );
void tetrahedron_nco_subrule ( int rule, int suborder_num, 
  double suborder_xyz[], double suborder_w[] );
void tetrahedron_nco_subrule_01 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void tetrahedron_nco_subrule_02 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void tetrahedron_nco_subrule_03 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void tetrahedron_nco_subrule_04 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void tetrahedron_nco_subrule_05 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void tetrahedron_nco_subrule_06 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void tetrahedron_nco_subrule_07 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
double tetrahedron_volume ( double t[3*4] );
void timestamp ( );

