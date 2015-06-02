int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4mat_print ( int m, int n, int a[], char *title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
void i4vec_print ( int n, int a[], char *title );
double r8_uniform_01 ( int *seed );
void rcont ( int nrow, int ncol, int nrowt[], int ncolt[], int nsubt[], 
  int matrix[], bool *key, int *ifault );
void timestamp ( );
