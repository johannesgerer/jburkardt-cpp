int i4_min ( int i1, int i2 );
int i4cvv_iget ( int mn, int a[], int m, int roff[], int i, int j );
void i4cvv_iinc ( int mn, int a[], int m, int roff[], int i, int j, 
  int daij );
int i4cvv_indx ( int m, int roff[], int i, int j );
void i4cvv_iset ( int mn, int a[], int m, int roff[], int i, int j, 
  int aij );
int *i4cvv_nget_new ( int mn, int a[], int m, int roff[], int nn, 
  int in[], int jn[] );
void i4cvv_ninc ( int mn, int a[], int m, int roff[], int nn, int in[], 
  int jn[], int dvn[] );
int *i4cvv_nndx ( int m, int roff[], int nn, 
  int in[], int jn[] );
void i4cvv_nset ( int mn, int a[], int m, int roff[], int nn, int in[], 
  int jn[], int vn[] );
int *i4cvv_offset ( int m, int nr[] );
void i4cvv_print ( int mn, int a[], int m, int roff[], string title );
int *i4cvv_rget_new ( int mn, int a[], int m, int roff[], int i );
void i4cvv_rinc ( int mn, int a[], int m, int roff[], int i, int dai[] );
int i4cvv_rndx ( int m, int roff[], int i );
void i4cvv_rset ( int mn, int a[], int m, int roff[], int i, int ai[] );
int i4cvv_size ( int m, int nr[] );
int i4vec_max ( int n, int a[] );
void i4vec_print ( int n, int a[], string title );
void i4vec_transpose_print ( int n, int a[], string title );
double r8cvv_iget ( int mn, double a[], int m, int roff[], int i, int j );
void r8cvv_iinc ( int mn, double a[], int m, int roff[], int i, int j, 
  double daij );
int r8cvv_indx ( int m, int roff[], int i, int j );
void r8cvv_iset ( int mn, double a[], int m, int roff[], int i, int j, 
  double aij );
double *r8cvv_nget_new ( int mn, double a[], int m, int roff[], int nn, 
  int in[], int jn[] );
void r8cvv_ninc ( int mn, double a[], int m, int roff[], int nn, int in[], 
  int jn[], double dvn[] );\
int *r8cvv_nndx ( int mn, double a[], int m, int roff[], int nn, 
  int in[], int jn[] );
void r8cvv_nset ( int mn, double a[], int m, int roff[], int nn, int in[], 
  int jn[], double vn[] );
int *r8cvv_offset ( int m, int nr[] );
void r8cvv_print ( int mn, double a[], int m, int roff[], string title );
double *r8cvv_rget_new ( int mn, double a[], int m, int roff[], int i );
void r8cvv_rinc ( int mn, double a[], int m, int roff[], int i, double dai[] );
int r8cvv_rndx ( int m, int roff[], int i );
void r8cvv_rset ( int mn, double a[], int m, int roff[], int i, double ai[] );
int r8cvv_size ( int m, int nr[] );
void r8vec_print ( int n, double a[], string title );
void r8vec_transpose_print ( int n, double a[], string title );
void timestamp ( );

