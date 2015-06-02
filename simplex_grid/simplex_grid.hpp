void comp_next_grlex ( int kc, int xc[] );
int *comp_random_new ( int n, int k, int &seed );
int i4_uniform_ab ( int a, int b, int &seed );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title );
int *ksub_random_new ( int n, int k, int &seed );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
int *simplex_grid_index_all ( int m, int n, int ng );
void simplex_grid_index_next ( int m, int n, int g[] );
int *simplex_grid_index_sample ( int m, int n, int &seed );
double *simplex_grid_index_to_point ( int m, int n, int ng, int g[],
  double v[] );
int simplex_grid_size ( int m, int n );
void timestamp ( );
