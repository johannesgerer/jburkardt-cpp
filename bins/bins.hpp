void bin_search_one_2d ( int bin[2], int nset, double pset[], int nbin[2],
  int bin_start[], int bin_next[], double ptest[2], bool *found_a_neighbor,
  int *i_min, double *d_min_sq, int *compares );
void bin_to_r8_even ( int nbin, int bin, double a, double b, double *cmin,
  double *cmax );
void bin_to_r8_even2 ( int nbin, int bin, double a, double b,
  double *cmin, double *cmax );
void bin_to_r82_even ( int nbin, int bin[2], double a[2], double b[2],
  double cmin[2], double cmax[2] );
void bin_to_r82_even2 ( int nbin, int bin[2], double a[2], double b[2],
  double cmin[2], double cmax[2] );
void bin_to_r82_even3 ( int nbin[2], int bin[2], double a[2], double b[2],
  double cmin[2], double cmax[2] );
void bin_to_r83_even2 ( int nbin, int bin[3], double a[3], double b[3],
  double cmin[3], double cmax[3] );
void bin_to_r83_even3 ( int nbin[3], int bin[3], double a[3], double b[3],
  double cmin[3], double cmax[3] );
int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 );
int get_seed ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_sign ( int i );
void i4_swap ( int *i, int *j );
int i4_uniform ( int b, int c, int *seed );
int i4_wrap ( int ival, int ilo, int ihi );
void i4mat_print ( int m, int n, int a[], char *title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void i4mat_transpose_print ( int m, int n, int a[], char *title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title );
void i4vec_heap_d ( int n, int a[] );
int *i4vec_indicator ( int n );
void i4vec_print ( int n, int a[], char *title );
void i4vec_sort_heap_a ( int n, int a[] );
int i4vec_sorted_unique ( int n, int a[] );
int i4vec2_compare ( int n, int a1[], int a2[], int i, int j );
void i4vec2_sort_a ( int n, int a1[], int a2[] );
void i4vec2_sorted_unique ( int n, int a1[], int a2[], int *nuniq );
void index_box2_next_2d ( int n1, int n2, int ic, int jc, int *i, int *j,
  int *more );
void index_box2_next_3d ( int n1, int n2, int n3, int ic, int jc, int kc,
  int *i, int *j, int *k, bool *more );
int lrline ( double xu, double yu, double xv1, double yv1, double xv2,
  double yv2, double dv );
bool perm_check ( int n, int p[] );
void perm_inv ( int n, int p[] );
int points_nearest_point_naive_2d ( int nset, double pset[], double ptest[],
  double *d_min );
int points_nearest_point_naive_3d ( int nset, double pset[], double ptest[],
  double *d_min );
int points_nearest_point_naive_nd ( int ndim, int nset, double pset[],
  double ptest[], double *d_min );
int *points_nearest_points_naive_2d ( int nset, double pset[], int ntest,
  double ptest[] );
int *points_nearest_points_naive_3d ( int nset, double pset[], int ntest,
  double ptest[] );
int r4_nint ( float x );
double r8_add ( double x, double y );
double r8_epsilon ( );
double r8_huge ( );
double r8_max ( double x, double y );
int r8_to_bin_even ( int nbin, double a, double b, double c );
int r8_to_bin_even2 ( int nbin, double a, double b, double c );
double r8_uniform ( double a, double b, int *seed );
double r8_uniform_01 ( int *seed );
int *r82_to_bin_even2 ( int nbin, double a[], double b[], double c[] );
int *r82_to_bin_even3 ( int nbin[], double a[], double b[], double c[] );
void r82_uniform ( double rlo[], double rhi[], int *seed, double r[] );
void r82vec_part_quick_a ( int n, double a[], int *l, int *r );
void r82vec_permute ( int n, double a[], int p[] );
void r82vec_print ( int n, double a[], char *title );
int *r82vec_sort_heap_index_a ( int n, double a[] );
void r82vec_sort_quick_a ( int n, double a[] );
void r82vec_uniform ( int n, double alo[], double ahi[], int *seed,
  double a[] );
int *r83_to_bin_even2 ( int nbin, double a[3], double b[3], double c[3] );
int *r83_to_bin_even3 ( int nbin[3], double a[3], double b[3], double c[3],
  int bin[3] );
void r83vec_part_quick_a ( int n, double a[], int *l, int *r );
void r83vec_sort_quick_a ( int n, double a[] );
void r83vec_uniform ( int n, double alo[], double ahi[], int *seed, double a[] );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void r8mat_transpose_print ( int m, int n, double a[], char *title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title );
void r8vec_bracket ( int n, double x[], double xval, int *left,
  int *right );
bool r8vec_eq ( int n, double a1[], double a2[] );
bool r8vec_gt ( int n, double a1[], double a2[] );
bool r8vec_lt ( int n, double a1[], double a2[] );
void r8vec_part_quick_a ( int n, double a[], int *l, int *r );
void r8vec_print ( int n, double a[], char *title );
void r8vec_sort_quick_a ( int n, double a[] );
void r8vec_swap ( int n, double a1[], double a2[] );
double *r8vec_uniform ( int n, double a, double b, int *seed );
int s_len_trim ( char *s );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
int swapec ( int i, int *top, int *btri, int *bedg, int point_num,
  double point_xy[], int tri_num, int tri_vert[], int tri_nabe[],
  int stack[] );
void timestamp ( );
double triangle_area_2d ( double t[2*3] );
void triangle_sample ( double t[2*3], int n, int *seed, double p[] );
void triangulation_print ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] );
void triangulation_sample ( int node_num, double node_xy[],
  int triangle_num, int triangle_node[], int num_ran, int *seed,
  double xd[], int td[] );
void tuple_next2 ( int n, int xmin[], int xmax[], int x[], int *rank );
void vbedg ( double x, double y, int point_num, double point_xy[], int tri_num,
  int tri_vert[], int tri_nabe[], int *ltri, int *ledg, int *rtri, int *redg );
