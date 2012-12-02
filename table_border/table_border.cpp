# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
bool ch_is_digit ( char c );
int ch_to_digit ( char c );
double *dtable_data_border_add ( int m, int n, double table[] );
double *dtable_data_read ( char *input_filename, int m, int n );
void dtable_data_write ( ofstream &output, int m, int n, double table[] );
void dtable_header_read ( char *input_filename, int *m, int *n );
void dtable_header_write ( int m, int n, char *output_filename, ofstream &output );
void dtable_print ( int m, int n, double a[], char *title );
void dtable_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
void dtable_write ( char *output_filename, int m, int n, double table[], 
  bool comment );
int file_column_count ( char *input_filename );
bool file_exist ( char *file_name );
void file_name_inc ( char *file_name );
int file_row_count ( char *input_filename );
int i4_huge ( void );
int i4_input ( char *string, bool *error );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double r8_epsilon ( void );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, bool *error );
double s_to_r8 ( char *s, int *lchar, bool *error );
bool s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
void timestamp ( void );
char *timestring ( void );


//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TABLE_BORDER.
//
//  Discussion:
//
//    TABLE_BORDER adds the "border" from TABLE data.
//
//    The TABLE data is assumed to be scalar values associated with
//    an M by N spatial grid, but which are stored in a TABLE file
//    as an M * N vector, with one number per record of the file,
//    and with columns preserved.
//
//    Although this is confusing, the data is stored in the file as 
//    though it were a vector, but logically it's really a 2D array,
//    and is treated as such once it is read in.
//
//    (Pie in the Sky): In a future version of this program, 
//    the data at each node will be allowed to be vector-valued.
//    That's one reason we are forcing the data currently to be listed
//    with just a single value per line in the file!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2005
//
//  Author:
//
//    John Burkardt
//
{
  int column_num;
  int column_num2;
  bool comment = false;
  bool error;
  int ihi;
  char input_file_name[100];
  char input_file_name_base[100];
  int input_file_num;
  int jhi;
  int last;
  int m;
  int m2;
  int n;
  int n2;
  char output_file_name[100];
  char output_file_name_base[100];
  int row_num;
  int row_num2;
  char string[100];
  double *u;
  double *v;

  timestamp ( );

  cout << "\n";
  cout << "TABLE_BORDER:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "* Read a vector representing scalar data on\n";
  cout << "  an M by N interior grid;\n";
  cout << "\n";
  cout << "* Write a vector representing that scalar data on\n";
  cout << "  the M+2 by N+2 grid with a boundary;\n";

  if ( argc <= 1 )
  {
    cout << "\n";
    cout << "Enter the name of the first input file:\n";
    cin >> input_file_name_base;
  }
  else
  {
    strcpy ( input_file_name_base, argv[1] );
  }

  strcpy ( input_file_name, input_file_name_base );
//
//  Get the output filename.
//
  if ( argc <= 2 )
  {
    cout << "\n";
    cout << "Enter the name of the first output file:\n";
    cin >> output_file_name_base;
  }
  else
  {
    strcpy ( output_file_name_base, argv[2] );
  }

  strcpy ( output_file_name, output_file_name_base );
//
//  Get the number of grid rows.
//
  if ( argc <= 3 )
  {
    m = i4_input ( "Enter the number of grid rows M", &error );

    if ( error )
    {
      cout << "\n";
      cout << "TABLE_BORDER - Fatal error!\n";
      cout << "  Error reading user input.\n";
      return 1;
    }
  }
  else
  {
    strcpy ( string, argv[3] );

    m = s_to_i4 ( string, &last, &error );

    if ( error )
    {
      cout << "\n";
      cout << "TABLE_BORDER - Fatal error!\n";
      cout << "  Error reading user input.\n";
      return 1;
    }

  }
//
//  Get the number of grid columns.
//
  if ( argc <= 4 )
  {
    n = i4_input ( "Enter the number of grid columns N", &error );

    if ( error )
    {
      cout << "\n";
      cout << "TABLE_BORDER - Fatal error!\n";
      cout << "  Error reading user input.\n";
      return 1;
    }
  }
  else
  {
    strcpy ( string, argv[4] );

    n = s_to_i4 ( string, &last, &error );

    if ( error )
    {
      cout << "\n";
      cout << "TABLE_BORDER - Fatal error!\n";
      cout << "  Error reading user input.\n";
      return 1;
    }

  }
//
//  Count the number of rows and columns in the input file.
//  For now, we expect COLUMN_NUM = 1.
//
  dtable_header_read ( input_file_name, &column_num, &row_num );
//
//  Do we accept the shape of this data?
//
  if ( column_num != 1 )
  {
    cout << "\n";
    cout << "TABLE_BORDER - Fatal error!\n";
    cout << "  Input file must have just 1 column of values.\n";
    return 1;
  }

  if ( m * n != row_num )
  {
    cout << "\n";
    cout << "TABLE_BORDER - Fatal error!\n";
    cout << "  M * N /= ROW_NUM.\n";
    cout << "  That is, the number of lines of data in the file\n";
    cout << "  which is ROW_NUM = " << row_num << "\n";
    cout << "  does not correspond to the product of\n";
    cout << "  the number of grid rows M = " << m << "\n";
    cout << "  and grid columns N = " << n << "\n";
    return 1;
  }

  cout << "\n";
  cout << "  First input file: \"" << input_file_name << "\"\n";
  cout << "\n";
  cout << "  Input data dimensions:\n";
  cout << "\n";
  cout << "  File rows ROW_NUM =       " << row_num << "\n";
  cout << "  File columns COLUMN_NUM = " << column_num << "\n";
  cout << "  Grid rows M =             " << m << "\n";
  cout << "  Grid columns N =          " << n << "\n";

  row_num2 = ( m + 2 ) * ( n + 2 );
  column_num2 = column_num;
  m2 = m + 2;
  n2 = n + 2;
  
  cout << "\n";
  cout << "  First output file: \"" << output_file_name << "\".\n";
  cout << "\n";
  cout << "  Output data dimensions:\n";
  cout << "\n";
  cout << "  File rows ROW_NUM2 =       " << row_num2    << "\n";
  cout << "  File columns COLUMN_NUM2 = " << column_num2 << "\n";
  cout << "  Grid rows M2 =             " << m2          << "\n";
  cout << "  Grid columns N2 =          " << n2          << "\n";
//
//  Allocate space.
//
  u = new double[row_num*column_num];
  v = new double[row_num2*column_num2];

  input_file_num = 0;

  for ( ; ; )
  {
    if ( !file_exist ( input_file_name ) )
    {
      break;
    }

    input_file_num = input_file_num + 1;
//
//  Read the U vectors.
//
    u = dtable_data_read ( input_file_name, column_num, row_num ) ;
//
//  Print a bit of the input array.
//
    if ( input_file_num == 1 )
    {
      ihi = i4_min ( m, 7 );
      jhi = i4_min ( n, 7 );

      dtable_print_some ( m, n, u, 1, 1, ihi, jhi, 
        "  Upper left corner of first input matrix:" );
    }
//
//  Add the border.
//
    v = dtable_data_border_add ( m, n, u );
//
//  Print a bit of the output array.
//
    if ( input_file_num == 1 )
    {
      ihi = i4_min ( m2, 7 );
      jhi = i4_min ( n2, 7 );

      dtable_print_some ( m2, n2, v, 1, 1, ihi, jhi, 
        "  Upper left corner of first output matrix:" );
    }
//
//  Write the V vectors.
//
    dtable_write ( output_file_name, column_num2, row_num2, v, comment );
//
//  Increment the file name.
//
    file_name_inc ( input_file_name );
    file_name_inc ( output_file_name );
  }

  cout << "\n";
  cout << "  Number of files processed was " << input_file_num << "\n";
//
//  Close up shop.
//
  delete [] u;
  delete [] v;

  cout << "\n";
  cout << "TABLE_BORDER:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

char ch_cap ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= c && c <= 122 ) 
  {
    c = c - 32;
  }   

  return c;
}
//****************************************************************************80*

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80*
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C1, char C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= c1 && c1 <= 122 ) 
  {
    c1 = c1 - 32;
  } 
  if ( 97 <= c2 && c2 <= 122 ) 
  {
    c2 = c2 - 32;
  }     

  return ( c1 == c2 );
}
//****************************************************************************80

bool ch_is_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_DIGIT returns TRUE if a character is a decimal digit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to be analyzed.
//
//    Output, bool CH_IS_DIGIT, is TRUE if C is a digit.
//
{
  if ( '0' <= c && c <= '9' )
  {
    return true;
  }
  else
  {
    return false;
  }
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= c && c <= '9' )
  {
    digit = c - '0';
  }
  else if ( c == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

double *dtable_data_border_add ( int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_DATA_BORDER_ADD adds a "border" to double precision table data.
//
//  Discussion:
//
//    We suppose the input data gives values of a quantity on nodes
//    in the interior of a 2D grid, and we wish to create a new table
//    with additional positions for the nodes that would be on the
//    border of the 2D grid.
//
//                  0 0 0 0 0 0
//      * * * *     0 * * * * 0
//      * * * * --> 0 * * * * 0
//      * * * *     0 * * * * 0
//                  0 0 0 0 0 0
//
//    The illustration suggests the situation in which a 3 by 4 array
//    is input, and a 5 by 6 array is to be output.
//
//    The old data is shifted to its correct positions in the new array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
//    Output, double TABLE2[(M+2)*(N+2)], the augmented table data.
//
{
  int i;
  int j;
  double *table2;

  table2 = new double[(m+2)*(n+2)];

  for ( j = 0; j < n+2; j++ )
  {
    for ( i = 0; i < m+2; i++ )
    {
      if ( i == 0 || i == m+1 || j == 0 || j == n+1 )
      {
        table2[i+j*(m+2)] = 0.0;
      }
      else
      {
        table2[i+j*(m+2)] = table[(i-1)+(j-1)*m];
      }
    }
  }

  return table2;
}
//****************************************************************************80

double *dtable_data_read ( char *input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_DATA_READ reads the data from a real TABLE file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with the '#' character are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double DTABLE_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  char line[255];
  double *table;
  double *x;

  input.open ( input_filename );

  if ( !input )
  {
    cout << "\n";
    cout << "DTABLE_DATA_READ - Fatal error!\n";
    cout << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new double[m*n];

  x = new double[m];

  j = 0;

  while ( j < n )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_r8vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80

void dtable_data_write ( ofstream &output, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_DATA_WRITE writes data to a real TABLE file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &OUTPUT, a pointer to the output stream.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  char *s;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << setw(10) << table[i+j*m] << "  ";
    }
    output << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80
 
void dtable_header_read ( char *input_filename, int *m, int *n )
 
//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_HEADER_READ reads the header from a real TABLE file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points.
//
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    cout << "\n";
    cout << "DTABLE_HEADER_READ - Fatal error!\n";
    cout << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cout << "\n";
    cout << "DTABLE_HEADER_READ - Fatal error!\n";
    cout << "  FILE_ROW_COUNT failed.\n";
    return;
  }

  return;
}
//****************************************************************************80
 
void dtable_header_write ( int m, int n, char *output_filename, 
  ofstream &output )
 
//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_HEADER_WRITE writes the header of a real TABLE file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, char *OUTPUT_FILENAME, the output filename.
//
//    Input, ofstream &OUTPUT, the output stream.
//
{
  char *s;

  s = timestring ( );

  output << "#  " << output_filename << "\n";
  output << "#  created by TABLE_BORDER.C" << "\n";
  output << "#  at " << s << "\n";
  output << "#\n";
  output << "#  Spatial dimension M = " << m << "\n";
  output << "#  Number of points N =  " << n << "\n";
  output << "#  EPSILON (unit roundoff) = " << r8_epsilon ( ) << "\n";
  output << "#\n";

  delete [] s;

  return;
}
//****************************************************************************80

void dtable_print ( int m, int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_PRINT prints a double precision matrix, with an optional title.
//
//  Discussion:
//
//    The doubly dimensioned array A is treated as a one dimensional vector,
//    stored by COLUMNS.  Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, char *TITLE, a title to be printed.
//
{
  dtable_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void dtable_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_PRINT_SOME prints some of a double precision matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, char *TITLE, a title for the matrix.
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }

  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void dtable_write ( char *output_filename, int m, int n, double table[],
  bool comment )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_WRITE writes information to a real TABLE file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
//    Input, bool COMMENT, is true if comments should be written.
{
  ofstream output;

  output.open ( output_filename );

  if ( !output )
  {
    cout << "\n";
    cout << "DTABLE_WRITE - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    return;
  }

  if ( comment )
  {
    dtable_header_write ( m, n, output_filename, output );
  }

  dtable_data_write ( output, m, n, table );

  output.close ( );

  return;
}
//****************************************************************************80

int file_column_count ( char *input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file is presumed to consist of COLUMN_NUM words, separated
//    by spaces.  There may also be some blank lines, and some comment lines,
//    which have a "#" in column 1.
//
//    The routine tries to find the first non-comment non-blank line and
//    counts the number of words in that line.
//
//    If all lines are blanks or comments, it goes back and tries to analyze
//    a comment line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  char line[256];
//
//  Open the file.
//
  input.open ( input_filename );

  if ( !input )
  {
    column_num = -1;
    cout << "\n";
    cout << "FILE_COLUMN_COUNT - Fatal error!\n";
    cout << "  Could not open the file:\n";
    cout << "  \"" << input_filename << "\"\n";
    return column_num;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      continue;
    }

    if ( line[0] == '#' )
    {
      continue;
    }

    got_one = true;
    break;

  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( input_filename );

    for ( ; ; )
    {
      input.getline ( line, sizeof ( line ) );

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( line ) == 0 )
      {
        continue;
      }

      got_one = true;
      break;

    }

  }

  input.close ( );

  if ( !got_one )
  {
    cout << "\n";
    cout << "FILE_COLUMN_COUNT - Warning!\n";
    cout << "  The file does not seem to contain any data.\n";
    return -1;
  }

  column_num = s_word_count ( line );

  return column_num;
}
//****************************************************************************80

bool file_exist ( char *file_name )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_EXIST reports whether a file exists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, character *FILE_NAME, the name of the file.
//
//    Output, bool FILE_EXIST, is TRUE if the file exists.
//
{
  ifstream file;

  file.open ( file_name, ios::in );

  if ( !file )
  {
    return false;
  }
  else
  {
    return true;
  }
}
//****************************************************************************80

void file_name_inc ( char *file_name )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_NAME_INC generates the next file name in a series.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected, and
//    if the name contains no digits, then nothing is done.
//
//  Example:
//
//      Input          Output
//      -----          ------
//      a7to11.txt     a7to12.txt  (typical case.  Last digit incremented)
//      a7to99.txt     a8to00.txt  (last digit incremented, with carry.)
//      a9to99.txt     a0to00.txt  (wrap around)
//      cat.txt        cat.txt     (no digits to increment)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, character *FILE_NAME, (a pointer to) the character string
//    to be incremented.
//
{
  char c;
  int i;
  int lens;

  lens = strlen ( file_name );

  for ( i = lens-1; 0 <= i; i-- )
  {
    c = *(file_name+i);

    if ( '0' <= c && c <= '9' )
    {

      if ( c == '9' )
      {
        c = '0';
        *(file_name+i) = c;
      }
      else
      {
        c = c + 1;
        *(file_name+i) = c;
        return;
      }
    }
  }

  return;
}
//****************************************************************************80

int file_row_count ( char *input_filename )

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
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  char line[100];
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( input_filename );

  if ( !input )
  {
    cout << "\n";
    cout << "FILE_ROW_COUNT - Fatal error!\n";
    cout << "  Could not open the input file: \"" << input_filename << "\"\n";
    return (-1);
  }

  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

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

int i4_huge ( void )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" integer.
//
{
  return 2147483647;
}
//****************************************************************************80*

int i4_input ( char *string, bool *error )

//****************************************************************************80*
//
//  Purpose:
//
//    I4_INPUT prints a prompt string and reads an integer from the user.
//
//  Discussion:
//
//    If the input line starts with a comment character ('#') or is
//    blank, the routine ignores that line, and tries to read the next one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *STRING, the prompt string.
//
//    Output, bool *ERROR, an error flag, which is true if an error occurred.
//
//    Output, integer I4_INPUT, the value input by the user.
//
{
  int last;
  char line[80];
  int value;

  *error = false;
  value = i4_huge ( );
//
//  Write the prompt.
//
  cout << "\n";
  cout << string << "\n";

  for ( ; ; )
  {
    cin.getline ( line, sizeof ( line ) );
//
//  If the line begins with a comment character, go back and read the next line.
//
    if ( line[0] == '#' )
    {
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      continue;
    }
//
//  Extract integer information from the string.
//
    value = s_to_i4 ( line, &last, error );

    if ( *error )
    {
      value = i4_huge ( );
      return value;
    }

    break;

  }

  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two integers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two integers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

double r8_epsilon ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 round off unit.
//
//  Discussion:
//
//    R8_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}
//****************************************************************************80

int s_len_trim ( char *s )

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
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n ) 
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

int s_to_i4 ( char *s, int *last, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4 reads an I4 from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a string to be examined.
//
//    Output, int *LAST, the last character of S used to make IVAL.
//
//    Output, bool *ERROR is TRUE if an error occurred.
//
//    Output, int *S_TO_I4, the integer value read from the string.
//    If the string is blank, then IVAL will be returned 0.
//
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = false;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  while ( *s ) 
  {
    c = s[i];
    i = i + 1;
//
//  Haven't read anything.
//
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read the sign, expecting digits.
//
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read at least one digit, expecting more.
//
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
//
//  If we read all the characters in the string, see if we're OK.
//
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = true;
    *last = 0;
  }

  return ival;
}
//****************************************************************************80

double s_to_r8 ( char *s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( char *s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  bool error;
  int i;
  int lchar;
  double x;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;

  }

  return error;
}
//****************************************************************************80

int s_word_count ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_COUNT counts the number of "words" in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int i;
  int nword;

  nword = 0;
  blank = true;

  while ( *s ) 
  {
    if ( *s == ' ' )
    {
      blank = true;
    }
    else if ( blank )
    {
      nword = nword + 1;
      blank = false;
    }
    *s++;
  }

  return nword;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

char *timestring ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTRING returns the current YMDHMS date as a string.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *TIMESTRING, a string containing the current YMDHMS date.
//
{
# define TIME_SIZE 40

  const struct tm *tm;
  size_t len;
  time_t now;
  char *s;

  now = time ( NULL );
  tm = localtime ( &now );

  s = new char[TIME_SIZE];

  len = strftime ( s, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  return s;
# undef TIME_SIZE
}
