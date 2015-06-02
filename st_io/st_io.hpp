int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4vec_dec ( int n, int a[] );
void i4vec_inc ( int n, int a[] );
int i4vec_max ( int n, int a[] );
int i4vec_min ( int n, int a[] );
void sort_heap_external ( int n, int &indx, int &i, int &j, int isgn );
void st_data_read ( string input_filename, int m, int n, int nnzero, 
  int row[], int col[], double a[] );
void st_header_print ( int i_min, int i_max, int j_min, int j_max, int m, 
  int n, int nst );
void st_header_read ( string input_filename, int &i_min, int &i_max, int &j_min, 
  int &j_max, int &m, int &n, int &nst );
void st_print ( int m, int n, int nnzero, int row[], int col[], 
  double a[], string title );
void st_print_some ( int row1, int row2, int col1, int col2, int nnzero,
  int row[], int col[], double a[], string title );
void st_sort_a ( int m, int n, int nnzero, int row[], int col[], 
  double a[] );
void st_transpose ( int &m, int &n, int nnzero, int row[], int col[], 
  double a[] );
void st_write ( string output_filename, int m, int n, int nnzero, 
  int row[], int col[], double a[] );
void timestamp ( );

