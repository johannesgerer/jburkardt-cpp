double e_01_2d ( int dim_num, double a[], double b[] );
double f_01_2d ( int dim_num, double x[] );
double f2 ( double x );
double f20_s ( int dim_num, double x[] );
int fibonacci ( int k );
double fibonacci_lattice_b ( int k, double f ( int dim_num, double x[] ) );
double fibonacci_lattice_q ( int k, double f ( int dim_num, double x[] ) );
double *fibonacci_lattice_q_nodes ( int k );
double fibonacci_lattice_q1 ( int k, double f ( int dim_num, double x[] ) );
double fibonacci_lattice_q2 ( int k, double f ( int dim_num, double x[] ) );
double fibonacci_lattice_q3 ( int k, double f ( int dim_num, double x[] ) );
double fibonacci_lattice_t ( int k, double f ( int dim_num, double x[] ) );
int *find_z20 ( int dim_num, int m );
void gray_next ( int n, int *change );
int i4_huge ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
void i4vec_print ( int n, int a[], string title );
double lattice ( int dim_num, int m, int z[], 
  double f ( int dim_num, double x[] ) );
double lattice_np0 ( int dim_num, int m, int z[], 
  double f ( int dim_num, double x[] ) );
double lattice_np1 ( int dim_num, int m, int z[], 
  double f ( int dim_num, double x[] ) );
void lattice_print ( int dim_num, int m, int z[], string title );
double monte_carlo ( int dim_num, int m, double f ( int dim_num, double x[] ),
  int *seed );
int prime ( int n );
double r8_abs ( double x );
double r8_huge ( );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
double *r8vec_uniform_01 ( int n, int *seed );
void timestamp ( );
void tuple_next ( int m1, int m2, int n, int *rank, int x[] );

