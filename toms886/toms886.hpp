void cheb ( int deg, double pt, double tcheb[] );
void dgemm ( char transa, char transb, int m, int n, int k, 
  double alpha, double a[], int lda, double b[], int ldb, double beta, 
  double c[], int ldc );
double franke ( double x, double y );
int i4_max ( int i1, int i2 );
void padua2 ( int deg, int degmax, int npd, double wpd[], double fpd[], 
  double raux1[], double raux2[], double c0[], double &esterr );
double pd2val ( int deg, int degmax, double c0[], double tg1, double tg2 );
void pdpts ( int deg, double pd1[], double pd2[], double wpd[], int &npd );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_sign ( double x );
void timestamp ( );
