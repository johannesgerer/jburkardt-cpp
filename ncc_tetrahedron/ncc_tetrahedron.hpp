int ncc_tetrahedron_degree ( int rule );
int ncc_tetrahedron_order_num ( int rule );
void ncc_tetrahedron_rule ( int rule, int order_num, double xyz[], double w[] );
int ncc_tetrahedron_rule_num ( );
int *ncc_tetrahedron_suborder ( int rule, int suborder_num );
int ncc_tetrahedron_suborder_num ( int rule );
void ncc_tetrahedron_subrule ( int rule, int suborder_num, 
  double suborder_xyz[], double suborder_w[] );
void ncc_tetrahedron_subrule_01 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_tetrahedron_subrule_02 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_tetrahedron_subrule_03 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_tetrahedron_subrule_04 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_tetrahedron_subrule_05 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_tetrahedron_subrule_06 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
void ncc_tetrahedron_subrule_07 ( int suborder_num, int suborder_xyz_n[],
  int *suborder_xyz_d, int suborder_w_n[], int *suborder_w_d );
double r8mat_det_4d ( double a[] );
double tetrahedron_volume ( double t[3*4] );
void reference_to_physical_t4 ( double t[], int n, double ref[], double phy[] );
void timestamp ( );

