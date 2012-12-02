# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <string>

using namespace std;

int main ( int argc, char *argv[] );
int file_column_count ( string input_filename );
bool file_exist ( string filename );
int file_row_count ( string input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void i4mat_write ( string output_filename, int m, int n, int table[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
int s_word_count ( string s );
void tet_mesh_order10_to_order4_compute ( int tetra_num1, int tetra_node1[], 
  int tetra_num2, int tetra_node2[] );
void tet_mesh_order10_to_order4_size ( int node_num1, int tetra_num1,
  int *node_num2, int *tetra_num2 );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TET_MESH_Q2L.
//
//  Discussion:
//
//    TET_MESH_Q2L makes a linear tet mesh from a quadratic one.
//
//    A quadratic tet mesh is assumed to consist of 10-node
//    tetrahedrons.  This routine rearranges information so as to 
//    define a 4-node tet mesh.
//
//  Usage:
//
//    tet_mesh_q2l prefix
//
//    where prefix is the common file prefix:
//
//    * prefix_nodes.txt,    the node coordinates (not needed by this program);
//    * prefix_elements.txt,    the linear element definitions.
//    * prefix_q2l_elements.txt,    the quadratic element definitions,
//                                created by the program.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 October 2009
//
//  Author:
//
//    John Burkardt
//
{
  string input_element_filename;
  int node_num1;
  int node_num2;
  string output_element_filename;
  string prefix;
  int *tetra_node1;
  int *tetra_node2;
  int tetra_num1;
  int tetra_num2;
  int tetra_order1;
  int tetra_order2 = 4;

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "TET_MESH_Q2L\n";
  cout << "  C++ version\n";
  cout << "  Read a \"quadratic\" tet mesh and\n";
  cout << "  write out a \"linear\" one.\n";
  cout << "\n";
  cout << "  Read a tet mesh of TETRA_NUM1 tetrahedrons\n";
  cout << "  using 10 nodes.\n";
  cout << "\n";
  cout << "  Create a 4 node tet mesh by breaking\n";
  cout << "  every 10 node tetrahedron into 8 smaller ones.\n";
  cout << "  Write the new linear tet mesh to a file.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
//
//  Get the filename prefix.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "TET_MESH_Q2L:\n";
    cout << "  Please enter the filename prefix.\n";

    cin >> prefix;
  }
  else 
  {
    prefix = argv[1];
  }
//
//  Create the filenames.
//
  input_element_filename = prefix + "_elements.txt";
  output_element_filename = prefix + "_q2l_elements.txt";
//
//  Read the tet mesh data.
//
  i4mat_header_read ( input_element_filename, &tetra_order1, &tetra_num1 );

  if ( tetra_order1 != 10 )
  {
    cout << "\n";
    cout << "TET_MESH_Q2L - Fatal error!\n";
    cout << "  The tet mesh must have order 10.\n";
    return 1;
  }

  cout << "\n";
  cout << "  Read the header of \"" << input_element_filename << "\".\n";
  cout << "\n";
  cout << "  Tetrahedron order = " << tetra_order1 << "\n";
  cout << "  Number of tetras  = " << tetra_num1 << "\n";

  tetra_node1 = i4mat_data_read ( input_element_filename, tetra_order1, 
    tetra_num1 );

  cout << "\n";
  cout << "  Read the data in \"" << input_element_filename << "\".\n";

  i4mat_transpose_print_some ( tetra_order1, tetra_num1, 
    tetra_node1, 1, 1, tetra_order1, 5, "  First 5 tetrahedrons:" );
//
//  Set the number of linear tetrahedrons:
//  We didn't read in the node data, so set that count to 0.
//
  node_num1 = 0;

  tet_mesh_order10_to_order4_size ( node_num1, tetra_num1, 
    &node_num2, &tetra_num2 );

  cout << "  Number of linear tetrahedrons =    " << tetra_num2 << "\n";
//
//  Allocate space.
//
  tetra_node2 = new int[tetra_order2*tetra_num2];
//
//  Convert the data.
//
  tet_mesh_order10_to_order4_compute ( tetra_num1, tetra_node1, 
    tetra_num2, tetra_node2 );

  i4mat_transpose_print_some ( tetra_order2, tetra_num2, tetra_node2, 
    1, 1, tetra_order2, 5, "  First 5 linear tetras" );
//
//  Write out the tetrahedron data for the quadratic mesh
//
  i4mat_write ( output_element_filename, tetra_order2, tetra_num2, 
    tetra_node2 );
    
  cout << "\n";
  cout << "  Created the file \"" << output_element_filename << "\".\n";
//
//  Deallocate memory.
//
  delete [] tetra_node1;
  delete [] tetra_node2;

  cout << "\n";
  cout << "TET_MESH_Q2L:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int file_column_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file are presumed to consist of COLUMN_NUM words, 
//    separated by spaces.  There may also be some blank lines, and some 
//    comment lines, which have a "#" in column 1.
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
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  string text;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    column_num = -1;
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    return column_num;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( text ) <= 0 )
    {
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
    got_one = true;
    break;
  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( filename.c_str ( ) );

    for ( ; ; )
    {
      input >> text;

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( text ) == 0 )
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
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Warning!\n";
    cerr << "  The file does not seem to contain any data.\n";
    return -1;
  }

  column_num = s_word_count ( text );

  return column_num;
}
//****************************************************************************80

bool file_exist ( string filename )

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
//    21 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, bool FILE_EXIST, is TRUE if the file exists.
//
{
  ifstream file;
  bool value;

  file.open ( filename.c_str ( ), ios::in );

  if ( !file )
  {
    value = false;
  }
  else
  {
    value = true;
  }
  return value;
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
    return (-1);
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

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
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
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
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
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int *i4mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_DATA_READ reads data from an I4MAT file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
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
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, int I4MAT_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  string line;
  int *table;
  int *x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "I4MAT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new int[m*n];

  x = new int[m];

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

    error = s_to_i4vec ( line, m, x );

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
 
void i4mat_header_read ( string input_filename, int *m, int *n )
 
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_HEADER_READ reads the header from an I4MAT file.
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
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points
//
{
  *m = file_column_count ( input_filename );
 
  if ( *m <= 0 )
  {
    cerr << "\n";
    cerr << "I4MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }
 
  *n = file_row_count ( input_filename );
 
  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "I4MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    return;
  }
 
  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2005
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
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title for the matrix.
{
# define INCX 10

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
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j << "  ";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void i4mat_write ( string output_filename, int m, int n, int table[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_WRITE writes an I4MAT file with no header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "I4MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(10) << table[i+j*m];
    }
    output << "\n";
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

int s_to_i4 ( string s, int *last, bool *error )

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
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be examined.
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

  for ( ; ; ) 
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

bool s_to_i4vec ( string s, int n, int ivec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4VEC reads an I4VEC from a string.
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
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, int IVEC[N], the values read from the string.
//
//    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    ivec[i] = s_to_i4 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
}
//****************************************************************************80

int s_word_count ( string s )

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
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int char_count;
  int i;
  int word_count;

  word_count = 0;
  blank = true;

  char_count = s.length ( );

  for ( i = 0; i < char_count; i++ )
  {
    if ( isspace ( s[i] ) )
    {
      blank = true;
    }
    else if ( blank )
    {
      word_count = word_count + 1;
      blank = false;
    }
  }

  return word_count;
}
//****************************************************************************80

void tet_mesh_order10_to_order4_compute ( int tetra_num1, int tetra_node1[], 
  int tetra_num2, int tetra_node2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER10_TO_ORDER4_COMPUTE linearizes a quadratic tet mesh.
//
//  Discussion:
//
//    A quadratic tet mesh is assumed to consist of 10-node
//    tetrahedrons.
//
//    This routine rearranges the information so as to define a 4-node
//    tet mesh.
//
//    The same nodes are used, but there are s8 times as many
//    tetrahedrons.
//
//    The node ordering for the quadratic tetrahedron is somewhat
//    arbitrary.  In the current scheme, the vertices are listed
//    first, followed by the 6 midside nodes.  Each midside node
//    may be identified by the two vertices that bracket it.  Thus,
//    the node ordering may be suggested by:
//
//     1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Anwei Liu, Barry Joe,
//    Quality Local Refinement of Tetrahedral Meshes Based
//    on 8-Subtetrahedron Subdivision,
//    Mathematics of Computation,
//    Volume 65, Number 215, July 1996, pages 1183-1200.
//
//  Parameters:
//
//    Input, int TETRA_NUM1, the number of tetrahedrons in the quadratic
//    tet mesh.
//
//    Input, int TETRA_NODE1[10*TETRA_NUM1], the indices of the nodes
//    that made up the quadratic mesh.
//
//    Input, int TETRA_NUM2, the number of tetrahedrons in the linear
//    tet mesh.  TETRA_NUM2 = 8 * TETRA_NUM1.
//
//    Output, int TETRA_NODE2[4*TETRA_NUM2], the indices of the nodes
//    that make up the linear mesh.
//
{
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int n7;
  int n8;
  int n9;
  int nx;
  int tetra1;
  int tetra2;

  tetra2 = 0;

  for ( tetra1 = 0; tetra1 < tetra_num1; tetra1++ )
  {
    n1 = tetra_node1[0+tetra1*10];
    n2 = tetra_node1[1+tetra1*10];
    n3 = tetra_node1[2+tetra1*10];
    n4 = tetra_node1[3+tetra1*10];
    n5 = tetra_node1[4+tetra1*10];
    n6 = tetra_node1[5+tetra1*10];
    n7 = tetra_node1[6+tetra1*10];
    n8 = tetra_node1[7+tetra1*10];
    n9 = tetra_node1[8+tetra1*10];
    nx = tetra_node1[9+tetra1*10];

    tetra_node2[0+tetra2*4] = n1;
    tetra_node2[1+tetra2*4] = n5;
    tetra_node2[2+tetra2*4] = n6;
    tetra_node2[3+tetra2*4] = n7;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n2;
    tetra_node2[1+tetra2*4] = n5;
    tetra_node2[2+tetra2*4] = n8;
    tetra_node2[3+tetra2*4] = n9;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n3;
    tetra_node2[1+tetra2*4] = n6;
    tetra_node2[2+tetra2*4] = n8;
    tetra_node2[3+tetra2*4] = n9;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n4;
    tetra_node2[1+tetra2*4] = n7;
    tetra_node2[2+tetra2*4] = n9;
    tetra_node2[3+tetra2*4] = nx;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n5;
    tetra_node2[1+tetra2*4] = n6;
    tetra_node2[2+tetra2*4] = n7;
    tetra_node2[3+tetra2*4] = n9;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n5;
    tetra_node2[1+tetra2*4] = n6;
    tetra_node2[2+tetra2*4] = n8;
    tetra_node2[3+tetra2*4] = n9;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n6;
    tetra_node2[1+tetra2*4] = n7;
    tetra_node2[2+tetra2*4] = n9;
    tetra_node2[3+tetra2*4] = nx;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n6;
    tetra_node2[1+tetra2*4] = n8;
    tetra_node2[2+tetra2*4] = n9;
    tetra_node2[3+tetra2*4] = nx;
    tetra2 = tetra2 + 1;
  }

  return;
}
//****************************************************************************80

void tet_mesh_order10_to_order4_size ( int node_num1, int tetra_num1,
  int *node_num2, int *tetra_num2 )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER10_TO_ORDER4_SIZE sizes a linear tet mesh from a quadratic one.
//
//  Discussion:
//
//    A linear (4 node) tet mesh can be derived from a quadratic
//    (10 node) tet mesh using the same set of nodes, but reassigning
//    the nodes of each quadratic tet among 8 linear subtets.
//
//    This routine returns the number of nodes and tetrahedra in the
//    linear mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Anwei Liu, Barry Joe,
//    Quality Local Refinement of Tetrahedral Meshes Based
//    on 8-Subtetrahedron Subdivision,
//    Mathematics of Computation,
//    Volume 65, Number 215, July 1996, pages 1183-1200.
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes in the quadratic mesh.
//
//    Input, int TETRA_NUM1, the number of tetrahedrons in the
//    quadratic mesh.
//
//    Output, int *NODE_NUM2, the number of nodes for the linear mesh.
//
//    Output, int *TETRA_NUM2, the number of tetrahedrons in the
//    linear mesh.
//
{
  *node_num2 = node_num1;
  *tetra_num2 = 8 * tetra_num1;

  return;
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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2003
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
