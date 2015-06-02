double *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  double x[] );
void cc_print ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  string title );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4vec_copy_new ( int n, int a1[] );
void i4vec_dec ( int n, int a[] );
void i4vec_inc ( int n, int a[] );
int i4vec_max ( int n, int a[] );
int i4vec_min ( int n, int a[] );
void i4vec_write ( string output_filename, int n, int table[] );
int i4vec2_compare ( int n, int a1[], int a2[], int i, int j );
void i4vec2_sort_a ( int n, int a1[], int a2[] );
int i4vec2_sorted_unique_count ( int n, int a1[], int a2[] );
void i4vec2_sorted_uniqueLly ( int n1, int a1[], int b1[], int n2, int a2[], 
  int b2[] );
double r8_uniform_01 ( int &seed );
double r8vec_diff_norm ( int n, double a[], double b[] );
void r8vec_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );
void r8vec_write ( string output_filename, int n, double x[] );
double *r8vec_zero_new ( int n );
void sort_heap_external ( int n, int &indx, int &i, int &j, int isgn );
void st_data_read ( string input_filename, int nrow, int ncol, int nnzero, 
  int row[], int col[], double a[] );
void st_header_print ( int i_min, int i_max, int j_min, int j_max, int m, 
  int n, int nst );
void st_header_read ( string input_filename, int &i_min, int &i_max, int &j_min, 
  int &j_max, int &m, int &n, int &nst );
double *st_mv ( int m, int n, int nst, int ist[], int jst[], double ast[], 
  double x[] );
void st_print ( int m, int n, int nst, int ist[], int jst[], double ast[], 
  string title );
int st_to_cc_size ( int nst, int ist[], int jst[] );
void st_to_cc_index ( int nst, int ist[], int jst[], int ncc, int n, 
  int icc[], int ccc[] );
double *st_to_cc_values ( int nst, int ist[], int jst[], double ast[], int ncc, 
  int n, int icc[], int ccc[] );
void timestamp ( );
double *wathen_st ( int nx, int ny, int nz_num, int &seed, int row[], 
  int col[] );
int wathen_st_size ( int nx, int ny );
