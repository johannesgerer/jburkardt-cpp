# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <fstream>
# include <cstring>
# include <ctime>
# include <cmath>

using namespace std;

# include "umfpack.h"

int main ( );
void cc_data_read ( string prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] );
void cc_header_read ( string prefix, int &ncc, int &n );
double *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  double x[] );
int file_row_count ( string input_filename );
void i4vec_data_read ( string input_filename, int n, int a[] );
void r8vec_data_read ( string input_filename, int n, double x[] );
double r8vec_diff_norm ( int n, double a[], double b[] );
double *r8vec_uniform_01_new ( int n, int &seed );
int s_len_trim ( string s );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for UMFPACK_WEST.
//
//  Discussion:
//
//    This program uses UMFPACK to solve a linear system A*X=B for which the
//    matrix is stored, in compressed column (CC) format, in three files.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Timothy Davis,
//    UMFPACK User Guide,
//    Version 5.6.2, 25 April 2013
//    http://suitesparse.com
//
{
  double *acc;
  double *b;
  int *ccc;
  int i;
  int *icc;
  int m;
  int n;
  int ncc;
  double *null = ( double * ) NULL;
  void *Numeric;
  char prefix[] = "west";
  double r;
  int seed;
  int status;
  void *Symbolic;
  double *x1;
  double *x2;

  timestamp ( );
  cout << "\n";
  cout << "UMFPACK_WEST:\n";
  cout << "  C++ version\n";
  cout << "  Use UMFPACK to solve the sparse linear system A*x=b.\n";
  cout << "  The matrix A is stored, in CC format, in 3 files.\n";
//
//  Get the matrix size.
//
  cc_header_read ( prefix, ncc, n );
//
//  Allocate space.
//
  acc = new double[ncc];
  ccc = new int[n+1];
  icc = new int[ncc];
//
//  Read the matrix data.
//
  cc_data_read ( prefix, ncc, n, icc, ccc, acc );
//
//  Set up the solution.
//
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, seed );
//
//  Set the right hand side.
//
  m = n;
  b = cc_mv ( m, n, ncc, icc, ccc, acc, x1 );
//
//  From the matrix data, create the symbolic factorization information.
//
  status = umfpack_di_symbolic ( n, n, ccc, icc, acc, &Symbolic, null, null );
//
//  From the symbolic factorization information, carry out the numeric factorization.
//
  status = umfpack_di_numeric ( ccc, icc, acc, Symbolic, &Numeric, null, null );
//
//  Free the symbolic factorization memory.
//
  umfpack_di_free_symbolic ( &Symbolic );
//
//  Using the numeric factorization, solve the linear system.
//
  x2 = new double[n];
  status = umfpack_di_solve ( UMFPACK_A, ccc, icc, acc, x2, b, Numeric, null, null );
//
//  Free the numeric factorization.
//
  umfpack_di_free_numeric ( &Numeric );
//
//  Compute the error:
//
  r = r8vec_diff_norm ( n, x1, x2 );
  cout << "\n";
  cout << "  Residual: ||A*x-b|| = " << r << "\n";
//
//  Free memory.
//
  delete [] acc;
  delete [] b;
  delete [] ccc;
  delete [] icc;
  delete [] x1;
  delete [] x2;
//
//  Terminate.
//
  cout << "\n";
  cout << "UMFPACK_WEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void cc_data_read ( string prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] )

//****************************************************************************80
//
//  Purpose:
//
//    CC_DATA_READ reads data about a sparse matrix in CC format.
//
//  Discussion:
//
//    Three files are presumed to exist:
//    * prefix_icc.txt contains NCC ICC values;
//    * prefix_ccc.txt contains N+1 CCC values;
//    * prefix_acc.txt contains NCC ACC values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PREFIX, a common prefix for the filenames.
//
//    Input, int NCC, the number of CC elements.
//
//    Input, int N, the number of columns in the matrix.
//
//    Output, int ICC[NCC], the CC rows.
//
//    Output, int CCC[N+1], the compressed CC columns.
//
//    Output, double ACC[NCC], the CC values.
//
{
  string filename_acc;
  string filename_ccc;
  string filename_icc;

  filename_icc = prefix + "_icc.txt";
  i4vec_data_read ( filename_icc, ncc, icc );

  filename_ccc = prefix + "_ccc.txt";
  i4vec_data_read ( filename_ccc, n + 1, ccc );

  filename_acc = prefix + "_acc.txt";
  r8vec_data_read ( filename_acc, ncc, acc );

  return;
}
//****************************************************************************80

void cc_header_read ( string prefix, int &ncc, int &n )

//****************************************************************************80
//
//  Purpose:
//
//    CC_HEADER_READ reads header information about a sparse matrix in CC format.
//
//  Discussion:
//
//    Three files are presumed to exist:
//    * prefix_icc.txt contains NCC ICC values;
//    * prefix_ccc.txt contains N+1 CCC values;
//    * prefix_acc.txt contains NCC ACC values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PREFIX, a common prefix for the filenames.
//
//    Output, int &NCC, the number of CC elements.
//
//    Output, int &N, the number of columns in the matrix.
//
{
  string filename_ccc;
  string filename_icc;

  filename_icc = prefix + "_icc.txt";
  ncc = file_row_count ( filename_icc );

  filename_ccc = prefix + "_ccc.txt";
  n = file_row_count ( filename_ccc ) - 1;

  return;
}
//****************************************************************************80

double *cc_mv ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    CC_MV multiplies a CC matrix by a vector
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Iain Duff, Roger Grimes, John Lewis,
//    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
//    October 1992
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, int NCC, the number of CC values.
//
//    Input, int RCC[NCC], the CC rows.
//
//    Input, int CCC[N+1], the compressed CC columns
//
//    Input, double ACC[NCC], the CC values.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Output, double CC_MV[M], the product A*X.
//
{
  double *b;
  int i;
  int j;
  int k;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( k = ccc[j]; k < ccc[j+1]; k++ )
    {
      i = icc[k];
      b[i] = b[i] + acc[k] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

int file_row_count ( string input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_ROW_COUNT counts the number of row records in a file.
//
//  Discussion:
//
//    It does not count lines that are blank, or that begin with a
//    comment symbol '#'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  string line;
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( line[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;

  }

  input.close ( );

  return row_num;
}
//****************************************************************************80

void i4vec_data_read ( string input_filename, int n, int table[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_DATA_READ reads data from an I4VEC file.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly one value.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, int TABLE[N], the data.
//
{
  ifstream input;
  int i;
  int j;
  int l;
  string line;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "I4VEC_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    table[j] = atoi ( line.c_str ( ) );
    j = j + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void r8vec_data_read ( string input_filename, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DATA_READ reads the data from an R8VEC file.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double TABLE[N], the data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  int lchar;
  string line;
  double x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8VEC_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    table[j] = atof ( line.c_str ( ) );
    j = j + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

double r8vec_diff_norm ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], B[N], the vectors.
//
//    Output, double R8VEC_DIFF_NORM, the L2 norm of A - B.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
