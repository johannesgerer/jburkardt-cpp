void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void legendre_com ( int order, double xtab[], double weight[] );
double *local_basis_1d ( int order, double node_x[], double x );
double *local_basis_prime_1d ( int order, double node_x[], double x );
double *local_fem_1d ( int order, double node_x[], double node_v[], int sample_num, 
  double sample_x[] );
double r8_uniform ( double a, double b, int *seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_sum ( int n, double a[] );
void timestamp ( );
