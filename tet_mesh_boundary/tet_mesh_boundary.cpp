# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <string>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
int file_column_count ( string input_filename );
int file_row_count ( string input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_sort_a ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
void i4i4i4_sort_a ( int i1, int i2, int i3, int *j1, int *j2, int *j3 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void i4mat_write ( string output_filename, int m, int n, int table[] );
int *i4vec_cum ( int n, int a[] );
int i4vec_sum ( int n, int a[] );
void mesh_base_zero ( int node_num, int element_order, int element_num,
  int element_node[] );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void tet_mesh_boundary_count ( int element_order, int element_num,
  int element_node[], int node_num, int *boundary_node_num,
  int *boundary_element_num, int boundary_node_mask[] );
int *tet_mesh_boundary_set ( int element_order, int element_num,
  int element_node[], int boundary_element_order, int boundary_element_num );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TET_MESH_BOUNDARY.
//
//  Discussion:
//
//    TET_MESH_BOUNDARY reads files defining a tetrahedral mesh, and
//    creates new files that define the 3D triangular mesh formed
//    by the boundary of the tetrahedral mesh.
//
//  Usage:
//
//    tet_mesh_boundary prefix
//
//    where
//
//    * prefix_nodes.txt contains the node coordinates;
//    * prefix_elements.txt contains the element definitions;
//    * 'prefix'_boundary_node_mask.txt is 0 for interior nodes,
//       1 for boundary nodes;
//    * prefix_boundary_nodes.txt contains the boundary nodes.
//    * prefix_boundary_elements.txt contains the boundary elements (faces).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 December 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *boundary_element;
  string boundary_element_filename;
  int *boundary_element_node;
  int boundary_element_num;
  int boundary_element_order;
  string boundary_node_filename;
  int *boundary_node_index;
  int *boundary_node_mask;
  string boundary_node_mask_filename;
  int boundary_node_num;
  double *boundary_node_xyz;
  int dim_num;
  string element_filename;
  int *element_neighbor;
  int *element_node;
  int element_num;
  int element_order;
  int i;
  int j;
  int j2;
  int node_num;
  string node_filename;
  double *node_xyz;
  string prefix;

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "TET_MESH_BOUNDARY\n";
  cout << "  C++ version\n";
  cout << "  Read files defining a tet mesh dataset.\n";
  cout << "\n";
  cout << "  Determine the faces that form the boundary;\n";
  cout << "  write a boundary node mask file;\n";
  cout << "  write a boundary node file;\n";
  cout << "  write a boundary element file,\n";
  cout << "  defining the boundary as a TRI_SURFACE.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
//
//  Get the filename prefix.
//
  if ( argc <= 1 )
  {
    cout << "\n";
    cout << "TET_MESH_BOUNDARY:\n";
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
  node_filename = prefix + "_nodes.txt";
  element_filename = prefix + "_elements.txt";
  boundary_node_mask_filename = prefix + "_boundary_node_mask.txt";
  boundary_node_filename = prefix + "_boundary_nodes.txt";
  boundary_element_filename = prefix + "_boundary_elements.txt";
//
//  Read the nodes.
//
  r8mat_header_read ( node_filename, &dim_num, &node_num );

  cout << "\n";
  cout << "  Read the header of \"" << node_filename << "\".\n";

  if ( dim_num != 3 )
  {
    cout << "\n";
    cout << "TET_MESH_BOUNDARY - Fatal error!\n";
    cout << "  The spatial dimension must be 3.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  Spatial dimension = " << dim_num << "\n";
  cout << "  Number of nodes   = " << node_num << "\n";

  node_xyz = r8mat_data_read ( node_filename, dim_num, node_num );

  cout << "\n";
  cout << "  Read the data in \"" << node_filename << "\".\n";

  r8mat_transpose_print_some ( dim_num, node_num,
    node_xyz, 1, 1, dim_num, 5, "  First 5 nodes:" );
//
//  Read the elements.
//
  i4mat_header_read ( element_filename, &element_order, &element_num );

  cout << "\n";
  cout << "  Read the header of \"" << element_filename << "\".\n";

  if ( element_order != 4 && element_order != 10 )
  {
    cout << "\n";
    cout << "TET_MESH_BOUNDARY - Fatal error!\n";
    cout << "  The elements must have order 4 or 10.\n";
    return 1;
  }

  cout << "\n";
  cout << "  Element order = " << element_order << "\n";
  cout << "  Number of tetras  = " << element_num << "\n";

  element_node = i4mat_data_read ( element_filename, element_order,
    element_num );

  cout << "\n";
  cout << "  Read the data in \"" << element_filename << "\".\n";

  i4mat_transpose_print_some ( element_order, element_num,
    element_node, 1, 1, element_order, 5, "  First 5 elements:" );
//
//  If the element information is 1-based, make it 0-based.
//
  mesh_base_zero ( node_num, element_order, element_num, element_node );
//
//  Count the boundary faces and nodes.
//
  boundary_node_mask = new int[node_num];

  tet_mesh_boundary_count ( element_order, element_num, element_node,
    node_num, &boundary_node_num, &boundary_element_num, boundary_node_mask );

  cout << "\n";
  cout << "  Number of faces is " << 4 * element_num << "\n";
  cout << "  Number of boundary faces is " << boundary_element_num << "\n";
  cout << "  Number of boundary nodes is " << boundary_node_num << "\n";
//
//  Set the boundary nodes and write them out.
//
  boundary_node_xyz = new double[dim_num*boundary_node_num];

  j2 = 0;
  for ( j = 0; j < node_num; j++ )
  {
    if ( boundary_node_mask[j] == 1 )
    {
      for ( i = 0; i < 3; i++ )
      {
        boundary_node_xyz[i+j2*3] = node_xyz[i+j*3];
      }
      j2 = j2 + 1;
    }
  }

  r8mat_write ( boundary_node_filename, dim_num, boundary_node_num,
    boundary_node_xyz );

  cout << "\n";
  cout << "  Wrote the boundary node coordinates to \""
       << boundary_node_filename << "\".\n";
//
//  Write the boundary node mask data.
//
  i4mat_write ( boundary_node_mask_filename, 1, node_num,
    boundary_node_mask );

  cout << "  Wrote the boundary node mask to \""
       << boundary_node_mask_filename << "\".\n";
//
//  Compute the reduced indices for the boundary nodes.
//
  boundary_node_index = i4vec_cum ( node_num, boundary_node_mask );
//
//  Subtract 1 so boundary node indices are zero-based.
//
  for ( i = 0; i < node_num; i++ )
  {
    boundary_node_index[i] = boundary_node_index[i] - 1;
  }
//
//  Set the boundary faces, applying the reduced index labels, and write them out.
//
  if ( element_order == 4 )
  {
    boundary_element_order = 3;
  }
  else if ( element_order == 10 )
  {
    boundary_element_order = 6;
  }

  boundary_element_node = tet_mesh_boundary_set ( element_order, element_num,
    element_node, boundary_element_order, boundary_element_num );

  for ( j = 0; j < boundary_element_num; j++ )
  {
    for ( i = 0; i < boundary_element_order; i++ )
    {
      boundary_element_node[i+j*boundary_element_order] =
        boundary_node_index [ boundary_element_node[i+j*boundary_element_order] ];
    }
  }

  i4mat_write ( boundary_element_filename, boundary_element_order,
    boundary_element_num, boundary_element_node );

  cout << "  Wrote the boundary face coordinates to \""
      << boundary_element_filename << "\".\n";
//
//  Deallocate memory.
//
  delete [] boundary_element_node;
  delete [] boundary_node_index;
  delete [] boundary_node_mask;
  delete [] boundary_node_xyz;
  delete [] element_node;
  delete [] node_xyz;
//
//  Terminate.
//
  cout << "\n";
  cout << "TET_MESH_BOUNDARY:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

char ch_cap ( char ch )

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
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 )
  {
    ch = ch - 32;
  }

  return ch;
}
//****************************************************************************80

bool ch_eqi ( char ch1, char ch2 )

//****************************************************************************80
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
//    Input, char CH1, CH2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= ch1 && ch1 <= 122 )
  {
    ch1 = ch1 - 32;
  }
  if ( 97 <= ch2 && ch2 <= 122 )
  {
    ch2 = ch2 - 32;
  }

  return ( ch1 == ch2 );
}
//****************************************************************************80

int ch_to_digit ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     CH  DIGIT
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
//    Input, char CH, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the character was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
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

int i4col_compare ( int m, int n, int a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_COMPARE compares columns I and J of an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 4
//
//      A = (
//        1  2  3  4
//        5  6  7  8
//        9 10 11 12 )
//
//    Output:
//
//      I4COL_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], an array of N columns of vectors of length M.
//
//    Input, int I, J, the columns to be compared.
//    I and J must be between 1 and N.
//
//    Output, int I4COL_COMPARE, the results of the comparison:
//    -1, column I < column J,
//     0, column I = column J,
//    +1, column J < column I.
//
{
  int k;
//
//  Check.
//
  if ( i < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index I = " << i << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < i )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index I = " << i << ".\n";
    exit ( 1 );
  }

  if ( j < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index J = " << j << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < j )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index J = " << j << ".\n";
    exit ( 1 );
  }

  if ( i == j )
  {
    return 0;
  }

  k = 1;

  while ( k <= m )
  {
    if ( a[k-1+(i-1)*m] < a[k-1+(j-1)*m] )
    {
      return (-1);
    }
    else if ( a[k-1+(j-1)*m] < a[k-1+(i-1)*m] )
    {
      return 1;
    }
    k = k + 1;
  }

  return 0;
}
//****************************************************************************80

void i4col_sort_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT_A ascending sorts the columns of an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//    In lexicographic order, the statement "X < Y", applied to two
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, X is less than Y if, at the first index where they
//    differ, the X value is less than the Y value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors;
//    On output, the columns of A have been sorted in ascending
//    lexicographic order.
//
{
  int i;
  int indx;
  int isgn;
  int j;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      i4col_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

void i4col_swap ( int m, int n, int a[], int icol1, int icol2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SWAP swaps two columns of an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based//  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], an array of data.
//
//    Input, int ICOL1, ICOL2, the two columns to swap.
//    These indices should be between 1 and N.
//
{
# define OFFSET 1

  int i;
  int t;
//
//  Check.
//
  if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL1 is out of range.\n";
    exit ( 1 );
  }

  if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL2 is out of range.\n";
    exit ( 1 );
  }

  if ( icol1 == icol2 )
  {
    return;
  }
  for ( i = 0; i < m; i++ )
  {
    t                     = a[i+(icol1-OFFSET)*m];
    a[i+(icol1-OFFSET)*m] = a[i+(icol2-OFFSET)*m];
    a[i+(icol2-OFFSET)*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

void i4i4i4_sort_a ( int i1, int i2, int i3, int *j1, int *j2, int *j3 )

//****************************************************************************80
//
//  Purpose:
//
//    I4I4I4_SORT_A ascending sorts a triple of I4's.
//
//  Discussion:
//
//    The program allows the reasonable call:
//
//      i4i4i4_sort_a ( i1, i2, i3, &i1, &i2, &i3 );
//
//    and this will return the reasonable result.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, I3, the values to sort.
//
//    Output, int *J1, *J2, *J3, the sorted values.
//
{
  int k1;
  int k2;
  int k3;
//
//  Copy arguments, so that the user can make "reasonable" calls like:
//
//    i4i4i4_sort_a ( i1, i2, i3, &i1, &i2, &i3 );
//
  k1 = i1;
  k2 = i2;
  k3 = i3;

  *j1 = i4_min ( i4_min ( k1, k2 ), i4_min ( k2, k3 ) );
  *j2 = i4_min ( i4_max ( k1, k2 ),
        i4_min ( i4_max ( k2, k3 ), i4_max ( k3, k1 ) ) );
  *j3 = i4_max ( i4_max ( k1, k2 ), i4_max ( k2, k3 ) );

  return;
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

int *i4vec_cum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_CUM computes the cumulutive sum of the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = (/ 1, 2, 3, 4 /)
//
//    Output:
//
//      A_CUM = (/ 1, 3, 6, 10 /)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_CUM[N], the cumulative sum of the entries of A.
//
{
  int *a_cum;
  int i;

  a_cum = new int[n+1];

  a_cum[0] = a[0];

  for ( i = 1; i <= n - 1; i++ )
  {
    a_cum[i] = a_cum[i-1] + a[i];
  }

  return a_cum;
}
//****************************************************************************80

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

void mesh_base_zero ( int node_num, int element_order, int element_num,
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_BASE_ZERO ensures that the element definition is zero-based.
//
//  Discussion:
//
//    The ELEMENT_NODE array contains nodes indices that form elements.
//    The convention for node indexing might start at 0 or at 1.
//    Since a C++ program will naturally assume a 0-based indexing, it is
//    necessary to check a given element definition and, if it is actually
//    1-based, to convert it.
//
//    This function attempts to detect 1-based node indexing and correct it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
//    definitions.
//
{
  int element;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = node_num + 1;
  node_max = -1;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      node_min = i4_min ( node_min, node );
      node_max = i4_max ( node_max, node );
    }
  }

  if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 1-based.\n";
    cout << "  This will be converted to 0-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] - 1;
      }
    }
  }
  else if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 0-based.\n";
    cout << "  No conversion is necessary.\n";
  }
  else
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
    cout << "  NODE_MIN = " << node_min << "\n";
    cout << "  NODE_MAX = " << node_max << "\n";
    cout << "  NODE_NUM = " << node_num << "\n";
  }
  return;
}
//****************************************************************************80

double *r8mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DATA_READ reads the data from an R8MAT file.
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
//    Output, double R8MAT_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  string line;
  double *table;
  double *x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8MAT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new double[m*n];

  x = new double[m];

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

void r8mat_header_read ( string input_filename, int *m, int *n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HEADER_READ reads the header from an R8MAT file.
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
//    Output, int *N, the number of points.
//
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    return;
  }

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, an optional title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file with no header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
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
//    Input, double TABLE[M*N], the table data.
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
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//  For greater precision, try
//
//    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
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

double s_to_r8 ( string s, int *lchar, bool *error )

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
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string containing the
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
//  Not if we terminated in states 1, 2, 6 or 7//
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

bool s_to_r8vec ( string s, int n, double rvec[] )

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
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
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
    rvec[i] = s_to_r8 ( s.substr(begin,length), &lchar, &error );

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

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      if ( n1 == 1 )
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 )
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 )
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 )
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 )
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
//****************************************************************************80

void tet_mesh_boundary_count ( int element_order, int element_num,
  int element_node[], int node_num, int *boundary_node_num,
  int *boundary_element_num, int boundary_node_mask[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_BOUNDARY_COUNT counts boundary faces and nodes in a tet mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the
//    nodes that make up each element.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, int *BOUNDARY_NODE_NUM, the number of boundary nodes.
//
//    Output, int *BOUNDARY_ELEMENT_NUM, the number of boundary elements.
//
//    Output, int BOUNDARY_NODE_MASK[NODE_NUM], is 0 for
//    interior nodes, 1 for boundary nodes.
//
{
  int a;
  int b;
  int c;
  int element;
  int element1;
  int element2;
  int f;
  int face;
  int face1;
  int face2;
  int *faces;
  int i;
  int j;
  int jj;
  int k;
  int l;

  faces = new int[5*(4*element_num)];
//
//  Step 1.
//  From the list of nodes forming tetrahedron T, of the form: 
//
//    (I,J,K,L)
//
//  or
//
//    (I,J,K,L,I+J,I+K,I+L,J+K,J+L,K+L),
//    (1,2,3,4, 5,  6,  7,  8,  9, 10 ),
//
//  construct the four face relations:
//
//    F = 1: (J,K,L,F,T)
//    F = 2: (I,K,L,F,T)
//    F = 3: (I,J,L,F,T)
//    F = 4: (I,J,K,F,T)
//
//  If T is actually order 10, we can retrieve the indices of the midside
//  nodes from the values of F and T.  In that case, the 4 faces are:
//
//    F = 1: 2, 3, 4, 8, 10, 9
//    F = 2: 1, 3, 4, 6, 10, 7
//    F = 3: 1, 2, 4, 5,  9, 7
//    F = 4: 1, 2, 3, 5,  8, 6
//
//  In order to make matching easier, we reorder each triple of nodes
//  into ascending order.
//
  for ( element = 0; element < element_num; element++ )
  {
    i = element_node[0+element*element_order];
    j = element_node[1+element*element_order];
    k = element_node[2+element*element_order];
    l = element_node[3+element*element_order];

    i4i4i4_sort_a ( j, k, l, &a, &b, &c );

    jj = 4 * element;
    faces[0+jj*5] = a;
    faces[1+jj*5] = b;
    faces[2+jj*5] = c;
    faces[3+jj*5] = 0;
    faces[4+jj*5] = element;

    i4i4i4_sort_a ( i, k, l, &a, &b, &c );

    jj = 4 * element + 1;
    faces[0+jj*5] = a;
    faces[1+jj*5] = b;
    faces[2+jj*5] = c;
    faces[3+jj*5] = 1;
    faces[4+jj*5] = element;

    i4i4i4_sort_a ( i, j, l, &a, &b, &c );

    jj = 4 * element + 2;
    faces[0+jj*5] = a;
    faces[1+jj*5] = b;
    faces[2+jj*5] = c;
    faces[3+jj*5] = 2;
    faces[4+jj*5] = element;

    i4i4i4_sort_a ( i, j, k, &a, &b, &c );

    jj = 4 * element + 3;
    faces[0+jj*5] = a;
    faces[1+jj*5] = b;
    faces[2+jj*5] = c;
    faces[3+jj*5] = 3;
    faces[4+jj*5] = element;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:3; the routine we call here
//  sorts on rows 1 through 5 but that won't hurt us.
//
//  What we need is to find cases where two tetrahedrons share a face.
//  By sorting the columns of the FACES array, we will put shared faces
//  next to each other.
//
  i4col_sort_a ( 5, 4*element_num, faces );
//
//  Step 3. Neighboring faces show up as consecutive columns with
//  identical first three entries.  Count columns which don't have
//  a following column that matches the first three entries.
//
  *boundary_element_num = 0;
  face = 0;

  for ( i = 0; i < node_num; i++ )
  {
    boundary_node_mask[i] = 0;
  }

  while ( face < 4 * element_num )
  {
    if ( face < 4 * element_num - 1 )
    {
      if ( faces[0+face*5] == faces[0+(face+1)*5] &&
           faces[1+face*5] == faces[1+(face+1)*5] &&
           faces[2+face*5] == faces[2+(face+1)*5] )
      {
        face = face + 2;
        continue;
      }
    }
    *boundary_element_num = *boundary_element_num + 1;
//
//  The vertices of the triangle are boundary nodes.
//
    boundary_node_mask[faces[0+face*5]] = 1;
    boundary_node_mask[faces[1+face*5]] = 1;
    boundary_node_mask[faces[2+face*5]] = 1;
//
//  For quadratic tetrahedrons, we need to add three more side nodes.
//  We retrieve the face index by F = FACES(4,*).
//  We can determine the local indices from the value of F.
//  We can determine the global indices from ELEMENT_NODE.
//
    if ( element_order == 10 )
    {
      f = faces[3+face*5];
      element = faces[4+face*5];

      if ( f == 0 )
      {
        boundary_node_mask[element_node[7+element*element_order]] = 1;
        boundary_node_mask[element_node[9+element*element_order]] = 1;
        boundary_node_mask[element_node[8+element*element_order]] = 1;
      }
      else if ( f == 1 )
      {
        boundary_node_mask[element_node[5+element*element_order]] = 1;
        boundary_node_mask[element_node[9+element*element_order]] = 1;
        boundary_node_mask[element_node[6+element*element_order]] = 1;
      }
      else if ( f == 2 )
      {
        boundary_node_mask[element_node[4+element*element_order]] = 1;
        boundary_node_mask[element_node[8+element*element_order]] = 1;
        boundary_node_mask[element_node[6+element*element_order]] = 1;
      }
      else if ( f == 3 )
      {
        boundary_node_mask[element_node[4+element*element_order]] = 1;
        boundary_node_mask[element_node[7+element*element_order]] = 1;
        boundary_node_mask[element_node[5+element*element_order]] = 1;
      }

    }
    face = face + 1;
  }

  *boundary_node_num = i4vec_sum ( node_num, boundary_node_mask );

  delete [] faces;

  return;
}
//****************************************************************************80

int *tet_mesh_boundary_set ( int element_order, int element_num,
  int element_node[], int boundary_element_order, int boundary_element_num )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_BOUNDARY_SET sets the boundary faces in a tet mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the
//    nodes that make up each element.
//
//    Input, int BOUNDARY_ELEMENT_ORDER, the order of the boundary faces.
//
//    Input, int BOUNDARY_ELEMENT_NUM, the number of boundary faces.
//
//    Output, int 
//    TET_MESH_BOUNDARY_SET[BOUNDARY_ELEMENT_ORDER*BOUNDARY_ELEMENT_NUM],
//    the nodes that make up each boundary face.
//
{
  int a;
  int b;
  int boundary_element;
  int *boundary_element_node;
  int c;
  int element;
  int element1;
  int element2;
  int f;
  int face;
  int face1;
  int face2;
  int *faces;
  int i;
  int j;
  int jj;
  int k;
  int l;

  faces = new int[5*(4*element_num)];
//
//  Step 1.
//  From the list of nodes forming tetrahedron T, of the form: 
//
//    (I,J,K,L)
//
//  or
//
//    (I,J,K,L,I+J,I+K,I+L,J+K,J+L,K+L),
//    (1,2,3,4, 5,  6,  7,  8,  9, 10 ),
//
//  construct the four face relations:
//
//    F = 1: (J,K,L,F,T)
//    F = 2: (I,K,L,F,T)
//    F = 3: (I,J,L,F,T)
//    F = 4: (I,J,K,F,T)
//
//  If T is actually order 10, we can retrieve the indices of the midside
//  nodes from the values of F and T.  In that case, the 4 faces are:
//
//    F = 1: 2, 3, 4, 8, 10, 9
//    F = 2: 1, 3, 4, 6, 10, 7
//    F = 3: 1, 2, 4, 5,  9, 7
//    F = 4: 1, 2, 3, 5,  8, 6
//
//  In order to make matching easier, we reorder each triple of nodes
//  into ascending order.
//
  for ( element = 0; element < element_num; element++ )
  {
    i = element_node[0+element*element_order];
    j = element_node[1+element*element_order];
    k = element_node[2+element*element_order];
    l = element_node[3+element*element_order];

    i4i4i4_sort_a ( j, k, l, &a, &b, &c );

    jj = 4 * element;
    faces[0+jj*5] = a;
    faces[1+jj*5] = b;
    faces[2+jj*5] = c;
    faces[3+jj*5] = 0;
    faces[4+jj*5] = element;

    i4i4i4_sort_a ( i, k, l, &a, &b, &c );

    jj = 4 * element + 1;
    faces[0+jj*5] = a;
    faces[1+jj*5] = b;
    faces[2+jj*5] = c;
    faces[3+jj*5] = 1;
    faces[4+jj*5] = element;

    i4i4i4_sort_a ( i, j, l, &a, &b, &c );

    jj = 4 * element + 2;
    faces[0+jj*5] = a;
    faces[1+jj*5] = b;
    faces[2+jj*5] = c;
    faces[3+jj*5] = 2;
    faces[4+jj*5] = element;

    i4i4i4_sort_a ( i, j, k, &a, &b, &c );

    jj = 4 * element + 3;
    faces[0+jj*5] = a;
    faces[1+jj*5] = b;
    faces[2+jj*5] = c;
    faces[3+jj*5] = 3;
    faces[4+jj*5] = element;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:3; the routine we call here
//  sorts on rows 1 through 5 but that won't hurt us.
//
//  What we need is to find cases where two tetrahedrons share a face.
//  By sorting the columns of the FACES array, we will put shared faces
//  next to each other.
//
  i4col_sort_a ( 5, 4*element_num, faces );
//
//  Step 3. Neighboring faces show up as consecutive columns with
//  identical first three entries.  Count columns which don't have
//  a following column that matches the first three entries.
//
  boundary_element = 0;
  face = 0;

  boundary_element_node = new int[boundary_element_order*boundary_element_num];

  while ( face < 4 * element_num )
  {
    if ( face < 4 * element_num - 1 )
    {
      if ( faces[0+face*5] == faces[0+(face+1)*5] &&
           faces[1+face*5] == faces[1+(face+1)*5] &&
           faces[2+face*5] == faces[2+(face+1)*5] )
      {
        face = face + 2;
        continue;
      }
    }

    f = faces[3+face*5];
    element = faces[4+face*5];

    if ( f == 0 )
    {
      boundary_element_node[0+boundary_element*boundary_element_order] = element_node[1+element*element_order];
      boundary_element_node[1+boundary_element*boundary_element_order] = element_node[2+element*element_order];
      boundary_element_node[2+boundary_element*boundary_element_order] = element_node[3+element*element_order];
    }
    else if ( f == 1 )
    {
      boundary_element_node[0+boundary_element*boundary_element_order] = element_node[0+element*element_order];
      boundary_element_node[1+boundary_element*boundary_element_order] = element_node[2+element*element_order];
      boundary_element_node[2+boundary_element*boundary_element_order] = element_node[3+element*element_order];
    }
    else if ( f == 2 )
    {
      boundary_element_node[0+boundary_element*boundary_element_order] = element_node[0+element*element_order];
      boundary_element_node[1+boundary_element*boundary_element_order] = element_node[1+element*element_order];
      boundary_element_node[2+boundary_element*boundary_element_order] = element_node[3+element*element_order];
    }
    else if ( f == 3 )
    {
      boundary_element_node[0+boundary_element*boundary_element_order] = element_node[0+element*element_order];
      boundary_element_node[1+boundary_element*boundary_element_order] = element_node[1+element*element_order];
      boundary_element_node[2+boundary_element*boundary_element_order] = element_node[2+element*element_order];
    }
//
//  For quadratic tetrahedrons, we need to add three more side nodes.
//
    if ( element_order == 10 )
    {
      if ( f == 0 )
      {
        boundary_element_node[3+boundary_element*boundary_element_order] = element_node[7+element*element_order];
        boundary_element_node[4+boundary_element*boundary_element_order] = element_node[9+element*element_order];
        boundary_element_node[5+boundary_element*boundary_element_order] = element_node[8+element*element_order];
      }
      else if ( f == 1 )
      {
        boundary_element_node[3+boundary_element*boundary_element_order] = element_node[5+element*element_order];
        boundary_element_node[4+boundary_element*boundary_element_order] = element_node[9+element*element_order];
        boundary_element_node[5+boundary_element*boundary_element_order] = element_node[6+element*element_order];
      }
      else if ( f == 2 )
      {
        boundary_element_node[3+boundary_element*boundary_element_order] = element_node[4+element*element_order];
        boundary_element_node[4+boundary_element*boundary_element_order] = element_node[8+element*element_order];
        boundary_element_node[5+boundary_element*boundary_element_order] = element_node[6+element*element_order];
      }
      else if ( f == 3 )
      {
        boundary_element_node[3+boundary_element*boundary_element_order] = element_node[4+element*element_order];
        boundary_element_node[4+boundary_element*boundary_element_order] = element_node[7+element*element_order];
        boundary_element_node[5+boundary_element*boundary_element_order] = element_node[5+element*element_order];
      }
    }

    boundary_element = boundary_element + 1;
    face = face + 1;
  }

  delete [] faces;

  return boundary_element_node;
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

