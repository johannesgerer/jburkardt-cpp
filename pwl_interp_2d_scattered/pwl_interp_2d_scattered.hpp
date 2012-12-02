int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void i4vec_heap_d ( int n, int a[] );
int i4vec_min ( int n, int a[] );
void i4vec_sort_heap_a ( int n, int a[] );
int i4vec_sorted_unique ( int n, int a[] );
int lrline ( double xu, double yu, double xv1, double yv1, double xv2,
  double yv2, double dv );
int perm_check2 ( int n, int p[], int base );
void perm_inverse ( int n, int p[] );
double *pwl_interp_2d_scattered_value ( int nd, double xyd[], double zd[], 
  int t_num, int t[], int t_neighbor[], int ni, double xyi[] );
int r8tris2 ( int node_num, double node_xy[], int &triangle_num,
  int triangle_node[], int triangle_neighbor[] );
int swapec ( int i, int &top, int &btri, int &bedg, int node_num,
  double node_xy[], int triangle_num, int triangle_node[],
  int triangle_neighbor[], int stack[] );
void triangulation_order3_print ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] );
void triangulation_search_delaunay ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[],
  int triangle_neighbor[], double p[2], int &triangle_index, 
  double &alpha, double &beta, double &gamma, int &edge,
  int &step_num );
void vbedg ( double x, double y, int node_num, double node_xy[],
  int triangle_num, int triangle_node[], int triangle_neighbor[],
  int &ltri, int &ledg, int &rtri, int &redg );
