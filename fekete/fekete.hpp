int fekete_degree ( int rule );
int fekete_order_num ( int rule );
void fekete_rule ( int rule, int order_num, double xy[], double w[] );
int fekete_rule_num ( );
int *fekete_suborder ( int rule, int suborder_num );
int fekete_suborder_num ( int rule );
void fekete_subrule ( int rule, int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
void fekete_subrule_1 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
void fekete_subrule_2 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
void fekete_subrule_3 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
void fekete_subrule_4 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
void fekete_subrule_5 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
void fekete_subrule_6 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
void fekete_subrule_7 ( int suborder_num, double suborder_xyz[], 
  double suborder_w[] );
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
char *timestring ( );
double triangle_area ( double t[2*3] );
void triangle_points_plot ( char *file_name, double node_xy[], int node_show, 
  int point_num, double point_xy[], int point_show );
