double *cc_compute_points ( int n );
void comp_next ( int n, int k, int a[], int &more, int &h,  int &t );
int i4_choose ( int n, int k );
int i4_mop ( int i );
double *lagrange_base_1d ( int nd, double xd[], int ni, double xi[] );
double *lagrange_interp_nd_grid2 ( int m, int ind[], double a[], double b[], 
  int nd );
int lagrange_interp_nd_size2 ( int m, int ind[] );
double *lagrange_interp_nd_value2 ( int m, int ind[], double a[], double b[], 
  int nd, double zd[], int ni, double xi[] );
int order_from_level_135 ( int l );
void smolyak_coefficients ( int l_max, int m, int c[], int w[] );
