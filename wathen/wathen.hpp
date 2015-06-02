void bandwidth ( int m, int n, double a[], int &b, int &l, int &d, int &u );
void cg_gb ( int n, int ml, int mu, double a[], double b[], double x[] );
void cg_ge ( int n, double a[], double b[], double x[] );
void cg_st ( int n, int nz_num, int row[], int col[], double a[], double b[], 
  double x[] );
double cpu_time ( );
void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );
double ddot ( int n, double dx[], int incx, double dy[], int incy );
int dgbfa ( double abd[], int lda, int n, int ml, int mu, int ipvt[] );
void dgbsl ( double abd[], int lda, int n, int ml, int mu, int ipvt[], 
  double b[], int job );
int dgefa ( double a[], int lda, int n, int ipvt[] );
void dgesl ( double a[], int lda, int n, int ipvt[], double b[], int job );
void dscal ( int n, double sa, double x[], int incx );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int idamax ( int n, double dx[], int incx );
double *mv_gb ( int m, int n, int ml, int mu, double a[], double x[] );
double *mv_ge ( int m, int n, double a[], double x[] );
double *mv_st ( int m, int n, int nz_num, int row[], int col[], double a[], 
  double x[] );
int nonzeros ( int m, int n, double a[] );
double r8_max ( double x, double y );
double r8_uniform_01 ( int &seed );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
double r8vec_diff_norm_li ( int n, double a[], double b[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
void r8vec_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );
void timestamp ( );
void wathen_bandwidth ( int nx, int ny, int &l, int &d, int &u );
double *wathen_gb ( int nx, int ny, int n, int &seed );
double *wathen_ge ( int nx, int ny, int n, int &seed );
int wathen_order ( int nx, int ny );
double *wathen_st ( int nx, int ny, int nz_num, int &seed, int row[], 
  int col[] );
int wathen_st_size ( int nx, int ny );
