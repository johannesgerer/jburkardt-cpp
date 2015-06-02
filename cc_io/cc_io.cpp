# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "cc_io.hpp"

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
//    18 July 2014
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
//    18 July 2014
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

void cc_print ( int m, int n, int ncc, int icc[], int ccc[], double acc[], 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    CC_PRINT prints a sparse matrix in CC format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in the matrix.
//
//    Input, int N, the number of columns in the matrix.
//
//    Input, int NCC, the number of CC elements.
//
//    Input, int ICC[NCC], the CC rows.
//
//    Input, int CCC[N+1], the compressed CC columns.
//
//    Input, double ACC[NCC], the CC values.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int j;
  int k;

  cout << "\n";
  cout << title << "\n";
  cout << "     #     I     J         A\n";
  cout << "  ----  ----  ----  ----------------\n";
  cout << "\n";

  if ( ccc[0] == 0 )
  {
    j = 0;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j+1] <= k )
      {
        j = j + 1;
      }
      cout << setw(4) << k << "  "
           << setw(4) << i << "  "
           << setw(4) << j << "  "
           << setw(16) << acc[k] << "\n";
    }
  }
  else
  {
    j = 1;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j] <= k + 1 )
      {
        j = j + 1;
      }
      cout << setw(4) << k + 1 << "  "
           << setw(4) << i << "  "
           << setw(4) << j << "  "
          << setw(16) << acc[k] << "\n";
    }
  }

  return;
}
//****************************************************************************80

void cc_print_some ( int i_min, int i_max, int j_min, int j_max, int ncc, 
  int n, int icc[], int ccc[], double acc[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    CC_PRINT_SOME prints some of a sparse matrix in CC format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, IMAX, the first and last rows to print.
//
//    Input, int J_MIN, J_MAX, the first and last columns 
//    to print.
//
//    Input, int NCC, the number of CC elements.
//
//    Input, int N, the number of columns.
//
//    Input, int ICC[NCC], the CC rows.
//
//    Input, int CCC[N+1], the compressed CC columns.
//
//    Input, double ACC[NCC], the CC values.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int j;
  int k;

  cout << "\n";
  cout << title << "\n";
  cout << "     #     I     J         A\n";
  cout << "  ----  ----  ----  ----------------\n";
  cout << "\n";

  if ( ccc[0] == 0 )
  {
    j = 0;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j+2] <= k )
      {
        j = j + 1;
      }
      if ( i_min <= i && i <= i_max &&
           j_min <= j && j <= j_max )
      {
        cout << setw(4) << k << "  "
             << setw(4) << i << "  "
             << setw(4) << j << "  "
             << setw(16) << acc[k] << "\n";
      }
    }
  }
  else
  {
    j = 1;
    for ( k = 0; k < ncc; k++ )
    {
      i = icc[k];
      while ( ccc[j+1] <= k + 1 )
      {
        j = j + 1;
      }
      if ( i_min <= i && i <= i_max &&
           j_min <= j && j <= j_max )
      {
        cout << setw(4) << k + 1 << "  "
             << setw(4) << i << "  "
             << setw(4) << j << "  "
             << setw(16) << acc[k] << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void cc_write ( string prefix, int ncc, int n, int icc[], int ccc[], 
  double acc[] )

//****************************************************************************80
//
//  Purpose:
//
//    CC_WRITE writes a sparse matrix in CC format to 3 files.
//
//  Discussion:
//
//    Three files will be created:
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
//    18 July 2014
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
//    Input, int ICC[NCC], the CC rows.
//
//    Input, int CCC[N+1], the compressed CC columns.
//
//    Input, double ACC[NCC], the CC values.
//
{
  string filename_acc;
  string filename_ccc;
  string filename_icc;

  filename_icc = prefix + "_icc.txt";
  i4vec_write ( filename_icc, ncc, icc );

  filename_ccc = prefix + "_ccc.txt";
  i4vec_write ( filename_ccc, n + 1, ccc );

  filename_acc = prefix + "_acc.txt";
  r8vec_write ( filename_acc, ncc, acc );

  return;
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

void i4vec_dec ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_DEC decrements an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, int A[N], the vector to be decremented.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] - 1;
  }
  return;
}
//****************************************************************************80

void i4vec_inc ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INC increments an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, int A[N], the vector to be incremented.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] + 1;
  }
  return;
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

void i4vec_write ( string output_filename, int n, int table[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_WRITE writes an I4VEC to a file.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int N, the number of points.
//
//    Input, int TABLE[N], the data.
//
{
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "I4VEC_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    output << table[j] << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
  ifstream input;
  int i;
  int j;
  int lchar;
  string line;

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

void r8vec_write ( string output_filename, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_WRITE writes an R8VEC file.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the data.
//
{
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8VEC_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    output << "  " << setw(24) << setprecision(16) << x[j] << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
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
