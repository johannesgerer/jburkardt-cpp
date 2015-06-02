bool between ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
bool collinear ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
bool diagonal ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
bool diagonalie ( int im1, int ip1, int n, int next[], double x[], double y[] );
int file_column_count ( string input_filename );
int file_row_count ( string input_filename );
void i4mat_print ( int m, int n, int a[], string title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void i4mat_write ( string output_filename, int m, int n, int table[] );
bool in_cone ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
bool intersect ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
bool intersect_prop ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
bool l4_xor ( bool l1, bool l2 );
void l4vec_print ( int n, bool a[], string title );
int *polygon_triangulate ( int n, double x[], double y[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int &m, int &n );
int s_len_trim ( string s );
double s_to_r8 ( string s, int &lchar, int &error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void timestamp ( );
double triangle_area ( double xa, double ya, double xb, double yb, double xc, 
  double yc );

