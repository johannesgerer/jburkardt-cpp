void cc_data_read ( string prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] );
void cc_header_read ( string prefix, int &ncc, int &n );
void cc_print ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  string title );
void cc_print_some ( int i_min, int i_max, int j_min, int j_max, int ncc, 
  int n, int icc[], int ccc[], double acc[], string title );
void cc_write ( string prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] );
int file_row_count ( string input_filename );
void i4vec_dec ( int n, int a[] );
void i4vec_inc ( int n, int a[] );
void i4vec_data_read ( string input_filename, int n, int table[] );
void i4vec_write ( string output_filename, int n, int table[] );
void r8vec_data_read ( string input_filename, int n, double table[] );
void r8vec_write ( string output_filename, int n, double x[] );
int s_len_trim ( string s );
void timestamp ( );
