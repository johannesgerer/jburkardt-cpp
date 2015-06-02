char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );
double ddot ( int n, double dx[], int incx, double dy[], int incy );
double dnrm2 ( int n, double x[], int incx );
void drot ( int n, double x[], int incx, double y[], int incy, double c, 
  double s );
void drotg ( double *sa, double *sb, double *c, double *s );
void dscal ( int n, double sa, double x[], int incx );
int dsvdc ( double a[], int lda, int m, int n, double s[], double e[], 
  double u[], int ldu, double v[], int ldv, double work[], int job );
void dswap ( int n, double x[], int incx, double y[], int incy );
int file_column_count ( string input_filename );
int file_row_count ( string input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_sign ( double x );
void r8col_normalize_li ( int m, int n, double a[] );
void r8col_reverse ( int m, int n, double a[] );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8mat_svd_linpack ( int m, int n, double a[], double u[], double s[], 
  double v[] );
double *r8mat_svd_low_rank ( int m, int n, int r, double u[], double s[], double v[] );
void r8row_reverse ( int m, int n, double a[] );
double *r8vec_cum0_new ( int n, double a[] );
void r8vec_print ( int n, double a[], string title );
double r8vec_sum ( int n, double a[] );
int s_len_trim ( string s );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void timestamp ( );
