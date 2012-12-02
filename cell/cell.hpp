int i4_min ( int i1, int i2 );
int i4vec_max ( int n, int a[] );
void i4vec_print ( int n, int a[], string title );
double r8cvv_iget ( int mn, double a[], int m, int roff[], int i, int j );
void r8cvv_iinc ( int mn, double a[], int m, int roff[], int i, int j, 
  double daij );
void r8cvv_iset ( int mn, double a[], int m, int roff[], int i, int j, 
  double aij );
double *r8cvv_nget_new ( int mn, double a[], int m, int roff[], int nn, 
  int in[], int jn[] );
void r8cvv_ninc ( int mn, double a[], int m, int roff[], int nn, int in[], 
  int jn[], double dvn[] );
void r8cvv_nset ( int mn, double a[], int m, int roff[], int nn, int in[], 
  int jn[], double vn[] );
int *r8cvv_offset ( int m, int nr[] );
void r8cvv_print ( int mn, double a[], int m, int roff[], string title );
double *r8cvv_rget_new ( int mn, double a[], int m, int roff[], int i );
void r8cvv_rinc ( int mn, double a[], int m, int roff[], int i, double dai[] );
void r8cvv_rset ( int mn, double a[], int m, int roff[], int i, double ai[] );
int r8cvv_size ( int m, int nr[] );
void r8vec_print ( int n, double a[], string title );
void r8vec_transpose_print ( int n, double a[], string title );
void timestamp ( );

