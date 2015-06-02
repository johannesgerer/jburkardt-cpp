int i4_choose ( int n, int k );
int i4_fall ( int x, int n );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4vec_concatenate ( int n1, int a[], int n2, int b[], int c[] );
void i4vec_permute ( int n, int p[], int a[] );
int *i4vec_sort_heap_index_a ( int n, int a[] );
int i4vec_sum ( int n, int a[] );
void mono_next_grlex ( int m, int x[] );
int mono_rank_grlex ( int m, int x[] );
void mono_total_next_grlex ( int m, int n, int x[] );
int *mono_unrank_grlex ( int m, int rank );
int mono_upto_enum ( int m, int n );
double *mono_value ( int m, int n, int f[], double x[] );
void perm_check0 ( int n, int p[] );
void polynomial_add ( int o1, double c1[], int e1[], int o2, double c2[], 
  int e2[], int &o, double c[], int e[] );
void polynomial_axpy ( double s, int o1, double c1[], int e1[], int o2, double c2[], 
  int e2[], int &o, double c[], int e[] );
void polynomial_compress ( int o1, double c1[], int e1[], int &o2, double c2[], 
  int e2[] );
void polynomial_dif ( int m, int o1, double c1[], int e1[], int dif[], 
  int &o2, double c2[], int e2[] );
void polynomial_mul ( int m, int o1, double c1[], int e1[], int o2, double c2[],
  int e2[], int &o, double c[], int e[] );
void polynomial_print ( int m, int o, double c[], int e[], string title );
void polynomial_scale ( double s, int m, int o1, double c1[], int e1[] );
void polynomial_sort ( int o, double c[], int e[] );
double *polynomial_value ( int m, int o, double c[], int e[], int nx, 
  double x[] );
void r8vec_concatenate ( int n1, double a[], int n2, double b[], double c[] );
void r8vec_permute ( int n, int p[], double a[] );
void timestamp ( );
