int get_seed ( );
int i4_uniform_ab ( int ilo, int ihi, int &seed );
void i4vec_print ( int n, int a[], string title );
double *latin_random_new ( int dim_num, int point_num, int &seed );
int *perm_uniform_new ( int n, int &seed );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

