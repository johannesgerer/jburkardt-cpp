int schdc ( float a[], int lda, int p, float work[], int ipvt[], int job );
int schdd ( float r[], int ldr, int p, float x[], float z[], int ldz, 
  int nz, float y[], float rho[], float c[], float s[] );
void schex ( float r[], int ldr, int p, int k, int l, float z[], int ldz, 
  int nz, float c[], float s[], int job );
void schud ( float r[], int ldr, int p, float x[], float z[], int ldz, int nz, 
  float y[], float rho[], float c[], float s[] );
float sgbco ( float abd[], int lda, int n, int ml, int mu, int ipvt[], 
  float z[] );
void sgbdi ( float abd[], int lda, int n, int ml, int mu, int ipvt[], 
  float det[2] );
int sgbfa ( float abd[], int lda, int n, int ml, int mu, int ipvt[] );
void sgbsl ( float abd[], int lda, int n, int ml, int mu, int ipvt[], 
  float b[], int job );
float sgeco ( float a[], int lda, int n, int ipvt[], float z[] );
void sgedi ( float a[], int lda, int n, int ipvt[], float det[], 
  float work[], int job );
int sgefa ( float a[], int lda, int n, int ipvt[] );
void sgesl ( float a[], int lda, int n, int ipvt[], float b[], int job );
int sgtsl ( int n, float c[], float d[], float e[], float b[] );
float spbco ( float abd[], int lda, int n, int m, float z[] );
void spbdi ( float abd[], int lda, int n, int m, float det[] );
int spbfa ( float abd[], int lda, int n, int m );
void spbsl ( float abd[], int lda, int n, int m, float b[] );
float spoco ( float a[], int lda, int n, float z[] );
void spodi ( float a[], int lda, int n, float det[], int job );
int spofa ( float a[], int lda, int n );
void sposl ( float a[], int lda, int n, float b[] );
float sppco ( float ap[], int n, float z[] );
void sppdi ( float ap[], int n, float det[2], int job );
int sppfa ( float ap[], int n );
void sppsl ( float ap[], int n, float b[] );
void sptsl ( int n, float d[], float e[], float b[] );
void sqrdc ( float a[], int lda, int n, int p, float qraux[], int jpvt[], 
  float work[], int job );
int sqrsl ( float a[], int lda, int n, int k, float qraux[], float y[], 
  float qy[], float qty[], float b[], float rsd[], float ab[], int job );
float ssico ( float a[], int lda, int n, int kpvt[], float z[] );
void ssidi ( float a[], int lda, int n, int kpvt[], float det[2], 
  int inert[3], float work[], int job );
int ssifa ( float a[], int lda, int n, int kpvt[] );
void ssisl ( float a[], int lda, int n, int kpvt[], float b[] );
float sspco ( float ap[], int n, int kpvt[], float z[] );
void sspdi ( float ap[], int n, int kpvt[], float det[2], int inert[3], 
  float work[], int job );
int sspfa ( float ap[], int n, int kpvt[] );
void sspsl ( float ap[], int n, int kpvt[], float b[] );
int ssvdc ( float a[], int lda, int n, int p, float s[], float e[], 
  float u[], int ldu, float v[], int ldv, float work[], int job );
float strco ( float t[], int ldt, int n, float z[], int job );
int strdi ( float t[], int ldt, int n, float det[], int job );
int strsl ( float t[], int ldt, int n, float b[], int job );
