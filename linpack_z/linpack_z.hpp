void drotg ( double *sa, double *sb, double *c, double *s );
int zchdc ( complex <double> a[], int lda, int p, int ipvt[], int job );
int zchdd ( complex <double> r[], int ldr, int p, complex <double> x[], 
  complex <double> z[], int ldz, int nz, complex <double> y[], double rho[], 
  double c[], complex <double> s[] );
void zchex ( complex <double> r[], int ldr, int p, int k, int l, 
  complex <double> z[], int ldz, int nz, double c[], complex <double> s[], int job );
void zchud ( complex <double> r[], int ldr, int p, complex <double> x[], 
  complex <double> z[], int ldz, int nz, complex <double> y[], double rho[], 
  double c[], complex <double> s[] );
double zgbco ( complex <double> abd[], int lda, int n, int ml, int mu, int ipvt[] );
void zgbdi ( complex <double> abd[], int lda, int n, int ml, int mu, int ipvt[], 
  complex <double> det[2] );
int zgbfa ( complex <double> abd[], int lda, int n, int ml, int mu, int ipvt[] );
void zgbsl ( complex <double> abd[], int lda, int n, int ml, int mu, 
  int ipvt[], complex <double> b[], int job );
double zgeco ( complex <double> a[], int lda, int n, int ipvt[] );
void zgedi ( complex <double> a[], int lda, int n, int ipvt[], 
  complex <double> det[2], int job );
int zgefa ( complex <double> a[], int lda, int n, int ipvt[] );
void zgesl ( complex <double> a[], int lda, int n, int ipvt[], 
  complex <double> b[], int job );
int zgtsl ( int n, complex <double> c[], complex <double> d[], 
  complex <double> e[], complex <double> b[] );
double zhico ( complex <double> a[], int lda, int n, int ipvt[] );
void zhidi ( complex <double> a[], int lda, int n, int ipvt[], double det[2], 
  int inert[3], int job );
int zhifa ( complex <double> a[], int lda, int n, int ipvt[] );
void zhisl ( complex <double> a[], int lda, int n, int ipvt[], 
  complex <double> b[] );
double zhpco ( complex <double> ap[], int n, int ipvt[] );
void zhpdi ( complex <double> ap[], int n, int ipvt[], double det[2], 
  int inert[3], int job );
int zhpfa ( complex <double> ap[], int n, int ipvt[] );
void zhpsl ( complex <double> ap[], int n, int ipvt[], complex <double> b[] );
double zpbco ( complex <double> abd[], int lda, int n, int m, int *info );
void zpbdi ( complex <double> abd[], int lda, int n, int m, double det[2] );
int zpbfa ( complex <double> abd[], int lda, int n, int m );
void zpbsl ( complex <double> abd[], int lda, int n, int m, complex <double> b[] );
double zpoco ( complex <double> a[], int lda, int n, int *info );
void zpodi ( complex <double> a[], int lda, int n, double det[2], int job );
int zpofa ( complex <double> a[], int lda, int n );
void zposl ( complex <double> a[], int lda, int n, complex <double> b[] );
double zppco ( complex <double> ap[], int n, int *info );
void zppdi ( complex <double> ap[], int n, double det[2], int job );
int zppfa ( complex <double> ap[], int n );
void zppsl ( complex <double> ap[], int n, complex <double> b[] );
void zptsl ( int n, complex <double> d[], complex <double> e[], 
  complex <double> b[] );
void zqrdc ( complex <double> x[], int ldx, int n, int p, 
  complex <double> qraux[], int ipvt[], int job );
int zqrsl ( complex <double> x[], int ldx, int n, int k, complex <double> qraux[], 
  complex <double> y[], complex <double> qy[], complex <double> qty[], 
  complex <double> b[], complex <double> rsd[], complex <double> xb[], int job );
double zsico ( complex <double> a[], int lda, int n, int ipvt[] );
void zsidi ( complex <double> a[], int lda, int n, int ipvt[], 
  complex <double> det[2], int job );
int zsifa ( complex <double> a[], int lda, int n, int ipvt[] );
void zsisl ( complex <double> a[], int lda, int n, int ipvt[], complex <double> b[] );
double zspco ( complex <double> ap[], int n, int ipvt[] );
void zspdi ( complex <double> ap[], int n, int ipvt[], complex <double> det[2], 
  int job );
int zspfa ( complex <double> ap[], int n, int ipvt[] );
void zspsl ( complex <double> ap[], int n, int ipvt[], complex <double> b[] );
int zsvdc ( complex <double> x[], int ldx, int n, int p, 
  complex <double> s[], complex <double> e[], complex <double> u[], int ldu, 
  complex <double> v[], int ldv, int job );
double ztrco ( complex <double> t[], int ldt, int n, int job );
int ztrdi ( complex <double> t[], int ldt, int n, complex <double> det[2], 
  int job );
int ztrsl ( complex <double> t[], int ldt, int n, complex <double> b[], int job );

