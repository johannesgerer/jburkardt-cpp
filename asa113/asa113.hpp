double crswap ( double varval[], int klass[], int clsize[], int in, int ik, 
  int iv, double *critvl, int i, int j, int l, int m, int iswitch );
double crtran ( double varval[], int klass[], int clsize[], int in, int ik, 
  int iv, double *critvl, int i, int m, int l, int iswitch );
void swap ( double varval[], int klass[], int clsize[], int in, int ik, int iv, 
  double *critvl, int *ntrans, int *ifault );
void timestamp ( );
void trnsfr ( double varval[], int klass[], int clsize[], int in, int ik, 
  int iv, double *critvl, int *ntrans, int *ifault );
