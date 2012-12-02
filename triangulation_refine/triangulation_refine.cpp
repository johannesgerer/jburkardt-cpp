# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
int file_column_count ( string filename );
int file_row_count ( string filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_sort_a ( int m, int n, int a[] );
int i4col_sorted_unique_count ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void i4mat_write ( string output_filename, int m, int n, int table[] );
void mesh_base_zero ( int node_num, int element_order, 
  int element_num, int element_node[] );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print ( int m, int n, double a[], string title );
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
void timestamp ( );
void triangulation_order3_refine_compute ( int node_num1, int triangle_num1, 
  double node_xy1[], int triangle_node1[], int node_num2, int triangle_num2, 
  int edge_data[], double node_xy2[], int triangle_node2[] );
void triangulation_order3_refine_size ( int node_num1, int triangle_num1, 
  int triangle_node1[], int *node_num2, int *triangle_num2, int edge_data[] );
void triangulation_order6_refine_compute ( int node_num1, int triangle_num1, 
  double node_xy1[], int triangle_node1[], int node_num2, int triangle_num2, 
  int edge_data[], double node_xy2[], int triangle_node2[] );
void triangulation_order6_refine_size ( int node_num1, int triangle_num1, 
  int triangle_node1[], int *node_num2, int *triangle_num2, int edge_data[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGULATION_REFINE.
//
//  Discussion:
//
//    TRIANGULATION_REFINE refines a triangulation by doubling.
//
//  Usage:
//
//    triangulation_refine prefix
//
//    where 'prefix' is the common filename prefix:
//
//    * prefix_nodes.txt contains the node coordinates,
//    * prefix_elements.txt contains the element definitions.
//    * prefix_ref_nodes.txt will contain the refined nodes;
//    * prefix_ref_elements.txt will contain the refined elements.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2009
//
//  Author:
//
//    John Burkardt
//
{
  bool debug = true;
  int dim_num;
  int *edge_data;
  string element_filename;
  string element_ref_filename;
  int i;
  int ierror = 0;
  string node_filename;
  int node_num1;
  int node_num2;
  string node_ref_filename;
  double *node_xy1;
  double *node_xy2;
  string prefix;
  int triangle_num1;
  int triangle_num2;
  int triangle_order;
  int *triangle_node1;
  int *triangle_node2;

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "TRIANGULATION_REFINE\n";
  cout << "  C++ version:\n";
  cout << "  Read a \"linear\" or \"quadratic\" triangulation\n";
  cout << "  and write out a refined triangulation.\n";
  cout << "\n";
  cout << "  In particular:\n";
  cout << "\n";
  cout << "  Read a dataset of NODE_NUM1 points in 2 dimensions.\n";
  cout << "\n";
  cout << "  Read an associated triangulation dataset of TRIANGLE_NUM1 \n";
  cout << "  triangles which use 3 or 6 nodes per triangle.\n";
  cout << "\n";
  cout << "  Subdivide each triangle into 4 triangles,\n";
  cout << "  generate new nodes as midpoints of current nodes.\n";
  cout << "\n";
  cout << "  Write out the new node and triangulation data.\n";
  cout << "\n";
  cout << "  If the input triangulation was Delaunay, then\n";
  cout << "  the output triangulation will be Delaunay.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
//
//  Get the filename prefix.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "TRIANGULATION_REFINE:\n";
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
  node_ref_filename = prefix + "_ref_nodes.txt";
  element_ref_filename = prefix + "_ref_elements.txt";
//
//  Read the node data.
//
  r8mat_header_read (  node_filename, &dim_num, &node_num1 );

  cout << "\n";
  cout << "  Read the header of \"" << node_filename << "\".\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  Number of nodes NODE_NUM1 = " << node_num1 << "\n";

  node_xy1 = r8mat_data_read ( node_filename, dim_num, node_num1 );

  cout << "\n";
  cout << "  Read the data in \"" << node_filename << "\".\n";

  r8mat_transpose_print_some ( dim_num, node_num1, node_xy1, 1, 1, dim_num, 5, 
    "  First 5 nodes:" );
//
//  Read the element data.
//
  i4mat_header_read ( element_filename, &triangle_order, &triangle_num1 );

  if ( triangle_order != 3 && triangle_order != 6 )
  {
    cout << "\n";
    cout << "TRIANGULATION_REFINE - Fatal error!\n";
    cout << "  Data is not for a 3-node or 6-node triangulation.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  Read the header of \"" << element_filename << "\".\n";
  cout << "\n";
  cout << "  Triangle order TRIANGLE_ORDER = " << triangle_order << "\n";
  cout << "  Number of triangles TRIANGLE_NUM1  = " << triangle_num1 << "\n";

  triangle_node1 = i4mat_data_read ( element_filename, 
    triangle_order, triangle_num1 );

  cout << "\n";
  cout << "  Read the data in \"" << element_filename << "\".\n";

  i4mat_transpose_print_some ( triangle_order, triangle_num1, triangle_node1, 
    1, 1, triangle_order, 5, "  First 5 triangles:" );
//
//  Detect and correct 1-based node indexing.
//
  mesh_base_zero ( node_num1, triangle_order, triangle_num1, triangle_node1 );
//
//  Determine the size of the refined mesh.
//
  edge_data = new int[5*(3*triangle_num1)];

  if ( triangle_order == 3 )
  {
    triangulation_order3_refine_size ( node_num1, triangle_num1, 
      triangle_node1, &node_num2, &triangle_num2, edge_data );
  }
  else if ( triangle_order == 6 )
  {
    triangulation_order6_refine_size ( node_num1, triangle_num1, 
      triangle_node1, &node_num2, &triangle_num2, edge_data );
  }
  cout << "\n";
  cout << "  Number of nodes in refined mesh =      " << node_num2 << "\n";
  cout << "  Number of triangles in refined mesh =  " 
       << triangle_num2 << "\n";

  node_xy2 = new double[dim_num*node_num2];
  triangle_node2 = new int[triangle_order*triangle_num2];
//
//  Compute the refined mesh.
//
  if ( triangle_order == 3 )
  {
    triangulation_order3_refine_compute ( node_num1, triangle_num1, 
      node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, node_xy2, 
      triangle_node2 );
  }
  else if ( triangle_order == 6 )
  {
    triangulation_order6_refine_compute ( node_num1, triangle_num1, 
      node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, node_xy2, 
      triangle_node2 );
  }

  if ( debug )
  {
    r8mat_transpose_print_some ( dim_num, node_num2, node_xy2, 
      1, 1, dim_num, 5, "  First 5 output nodes:" );

    i4mat_transpose_print_some ( triangle_order, triangle_num2, triangle_node2, 
      1, 1, triangle_order, 5, "  First 5 output triangles" );
  }
//
//  Write out the node and triangle data.
//
  r8mat_write ( node_ref_filename, dim_num, node_num2, node_xy2 );

  cout << "\n";
  cout << "  Wrote the refined node data to \"" << node_ref_filename << "\".\n";

  i4mat_write ( element_ref_filename, triangle_order, triangle_num2, 
    triangle_node2 );

  cout << "  Wrote the refined element data to \"" << element_ref_filename << "\".\n";

//
//  Free up memory.
//
  delete [] edge_data;
  delete [] node_xy1;
  delete [] node_xy2;
  delete [] triangle_node1;
  delete [] triangle_node2;

  cout << "\n";
  cout << "TRIANGULATION_REFINE:\n";
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
//    Most lines of the file are presumed to consist of COLUMN_NUM words, separated
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
//    30 January 2009
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
  char text[255];
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
    input.getline ( text, sizeof ( text ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( text ) == 0 )
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
      input.getline ( text, sizeof ( text ) );

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

int file_row_count ( string filename )

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
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  int record_num;
  int row_num;
  char text[255];

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the file: \"" << filename << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    input.getline ( text, sizeof ( text ) );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( text[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( text ) == 0 )
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

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If 
//      NREM = I4_MODP ( I, J ) 
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Example:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
// 
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
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
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is 
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cout << "\n";
    cout << "I4_MODP - Fatal error!\n";
    cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80*

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80*
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
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

int i4col_sorted_unique_count ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
//
//  Discussion:
//
//    The columns of the array may be ascending or descending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], a sorted array, containing
//    N columns of data.
//
//    Output, int I4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
//
{
  int i;
  int j1;
  int j2;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  unique_num = 1;
  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        unique_num = unique_num + 1;
        j1 = j2;
        break;
      }
    }
  }

  return unique_num;
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
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based!  However, a preprocessor
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

void i4mat_transpose_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 January 2005
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
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title to be printed.
//
{
  int i;
  int j;
  int jhi;
  int jlo;

  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

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
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2005
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
//
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
    cout << "  The element indexing appears to be 1-based!\n";
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
    cout << "  The element indexing appears to be 0-based!\n";
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

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
//    Input, string TITLE, an optional title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

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
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
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
    cout << "\n";

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
//    24 September 2003
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

void triangulation_order3_refine_compute ( int node_num1, int triangle_num1, 
  double node_xy1[], int triangle_node1[], int node_num2, int triangle_num2, 
  int edge_data[], double node_xy2[], int triangle_node2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_REFINE_COMPUTE computes a refined order 3 triangulation.
//
//  Discussion:
//
//    Given a triangle defined by nodes 1, 2, 3, we need to generate
//    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
//    and T4.
//
//    The task is more complicated by the fact that we are working with
//    a mesh of triangles, so that we want to create a node only once,
//    even though it may be shared by other triangles.
//
//          3
//         / \
//        /T3 \
//      13----23
//      / \T4 / \
//     /T1 \ /T2 \
//    1----12-----2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes.
//
//    Input, int TRIANGLE_NUM1, the number of triangles.
//
//    Input, double NODE_XY1[2*NODE_NUM1], the nodes.
//
//    Input, int TRIANGLE_NODE1[3*TRIANGLE_NUM1], the nodes that make up the
//    triangles.  These should be listed in counterclockwise order.
//
//    Input, int NODE_NUM2, the number of nodes in the refined mesh.
//
//    Input, int TRIANGLE_NUM2, the number of triangles in the refined mesh.
//
//    Input, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge information computed
//    by TRIANGULATION_ORDER3_REFINE_SIZE.
//
//    Output, double NODE_XY2[2*NODE_NUM2], the refined nodes.
//
//    Output, int TRIANGLE_NODE2[3*TRIANGLE_NUM2], the nodes that make up the
//    triangles in the refined mesh.
//
{
  int edge;
  int i;
  int j;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int node;
  int triangle1;
  int v1;
  int v2;
//
//  Copy the old nodes.
//
  for ( j = 0; j < node_num1; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      node_xy2[i+j*2] = node_xy1[i+j*2];
    }
  }
  for ( j = 0; j < triangle_num2; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      triangle_node2[i+j*3] = -1;
    }
  }
//
//  We can assign the existing nodes to the new triangles.
//
  for ( triangle1 = 0; triangle1 < triangle_num1; triangle1++ )
  {
    triangle_node2[0+(triangle1*4+0)*3] = triangle_node1[0+triangle1*3];
    triangle_node2[1+(triangle1*4+1)*3] = triangle_node1[1+triangle1*3];
    triangle_node2[2+(triangle1*4+2)*3] = triangle_node1[2+triangle1*3];
  }

  node = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 3 * triangle_num1; edge++ )
  {
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];
//
//  If this edge is new, create the coordinates and index for this node.
//
    if ( n1 != n1_old || n2 != n2_old )
    {

      if ( node_num2 < node )
      {
        cout << "\n";
        cout << "TRIANGLE_MESH_ORDER3_REFINE - Fatal error!\n";
        cout << "  Node index exceeds NODE_NUM2.\n";
        exit ( 1 );
      }

      for ( i = 0; i < 2; i++ )
      {
        node_xy2[i+node*2] = ( node_xy2[i+n1*2] + node_xy2[i+n2*2] ) / 2.0;
      }

      node = node + 1;

      n1_old = n1;
      n2_old = n2;
    }
//
//  Assign the node to triangles.
//
    v1 = edge_data[2+edge*5];
    v2 = edge_data[3+edge*5];
    triangle1 = edge_data[4+edge*5];

    if ( v1 == 1 && v2 == 2 )
    {
      triangle_node2[0+(triangle1*4+1)*3] = node;
      triangle_node2[1+(triangle1*4+0)*3] = node;
      triangle_node2[2+(triangle1*4+3)*3] = node;
    }
    else if ( v1 == 1 && v2 == 3 )
    {
      triangle_node2[0+(triangle1*4+2)*3] = node;
      triangle_node2[1+(triangle1*4+3)*3] = node;
      triangle_node2[2+(triangle1*4+0)*3] = node;
    }
    else if ( v1 == 2 && v2 == 3 )
    {
      triangle_node2[0+(triangle1*4+3)*3] = node;
      triangle_node2[1+(triangle1*4+2)*3] = node;
      triangle_node2[2+(triangle1*4+1)*3] = node;
    }
  }
  return;
}
//****************************************************************************80

void triangulation_order3_refine_size ( int node_num1, int triangle_num1, 
  int triangle_node1[], int *node_num2, int *triangle_num2, int edge_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_REFINE_SIZE sizes a refined order 3 triangulation.
//
//  Discussion:
//
//    Given a triangle defined by nodes 1, 2, 3, we need to generate
//    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
//    and T4.
//
//    The task is more complicated by the fact that we are working with
//    a mesh of triangles, so that we want to create a node only once,
//    even though it may be shared by other triangles.
//
//          3
//         / \
//        /T3 \
//      13----23
//      / \T4 / \
//     /T1 \ /T2 \
//    1----12-----2
//
//    This routine simply determines the sizes of the resulting node
//    and triangle arrays.
//
//    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
//    data items, one item for every edge of every triangle.  Each
//    data item records, for a given edge, the global indices
//    of the two endpoints, the local indices of the two endpoints,
//    and the index of the triangle.
//
//    Through careful sorting, it is possible to arrange this data in
//    a way that allows the proper generation of the interpolated nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes in the original mesh.
//
//    Input, int  TRIANGLE_NUM1, the number of triangles in the
//    original mesh.
//
//    Input, int TRIANGLE_NODE1[3*TRIANGLE_NUM1], the indices of the nodes
//    that form the triangles in the input mesh.
//
//    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
//
//    Output, int *TRIANGLE_NUM2, the number of triangles in the
//    refined mesh.
//
//    Output, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge data that will
//    be needed by TRIANGULATION_ORDER3_REFINE_COMPUTE.
//
{
  int a;
  int b;
  int edge;
  int i;
  int j;
  int k;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int triangle;
//
//  Step 1.
//  From the list of nodes for triangle T, of the form: (I,J,K)
//  construct the edge relations:
//
//    (I,J,1,2,T)
//    (I,K,1,3,T)
//    (J,K,2,3,T)
//
//  In order to make matching easier, we reorder each pair of nodes
//  into ascending order.
//
  for ( triangle = 0; triangle < triangle_num1; triangle++ )
  {
    i = triangle_node1[0+triangle*3];
    j = triangle_node1[1+triangle*3];
    k = triangle_node1[2+triangle*3];

    a = i4_min ( i, j );
    b = i4_max ( i, j );

    edge_data[0+5*(3*triangle+0)] = a;
    edge_data[1+5*(3*triangle+0)] = b;
    edge_data[2+5*(3*triangle+0)] = 1;
    edge_data[3+5*(3*triangle+0)] = 2;
    edge_data[4+5*(3*triangle+0)] = triangle;

    a = i4_min ( i, k );
    b = i4_max ( i, k );

    edge_data[0+5*(3*triangle+1)] = a;
    edge_data[1+5*(3*triangle+1)] = b;
    edge_data[2+5*(3*triangle+1)] = 1;
    edge_data[3+5*(3*triangle+1)] = 3;
    edge_data[4+5*(3*triangle+1)] = triangle;

    a = i4_min ( j, k );
    b = i4_max ( j, k );

    edge_data[0+5*(3*triangle+2)] = a;
    edge_data[1+5*(3*triangle+2)] = b;
    edge_data[2+5*(3*triangle+2)] = 2;
    edge_data[3+5*(3*triangle+2)] = 3;
    edge_data[4+5*(3*triangle+2)] = triangle;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:2; the routine we call here
//  sorts on the full column but that won't hurt us.
//
//  What we need is to find all cases where triangles share an edge.
//  By sorting the columns of the EDGE_DATA array, we will put shared edges
//  next to each other.
//
  i4col_sort_a ( 5, 3*triangle_num1, edge_data );
//
//  Step 3. All the triangles which share an edge show up as consecutive
//  columns with identical first two entries.  Figure out how many new
//  nodes there are, and allocate space for their coordinates.
//
  *node_num2 = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 3 * triangle_num1; edge++ )
  {
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];
    if ( n1 != n1_old || n2 != n2_old )
    {
      *node_num2 = *node_num2 + 1;
      n1_old = n1;
      n2_old = n2;
    }
  }

  *triangle_num2 = 4 * triangle_num1;

  return;
}
//****************************************************************************80

void triangulation_order6_refine_compute ( int node_num1, int triangle_num1, 
  double node_xy1[], int triangle_node1[], int node_num2, int triangle_num2, 
  int edge_data[], double node_xy2[], int triangle_node2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_REFINE_COMPUTE computes a refined order 6 triangulation.
//
//  Discussion:
//
//    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we 
//    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
//    quadratic subtriangles T1, T2, T3 and T4.
//
//    The task is more complicated by the fact that we are working with
//    a mesh of triangles, so that we want to create a node only once,
//    even though it may be shared by other triangles.  (In fact, only
//    the new nodes on the edges can be shared, and then only by at most
//    one other triangle.)
//
//            3
//           / \
//          36 35
//         / T3  \
//        6--56---5
//       / \ T4  / \
//      16 46  45  25
//     / T1  \ / T2  \
//    1--14---4--24---2
//
//    This routine is given sorted information defining the edges, and uses
//    it to build the new node and triangle arrays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes.
//
//    Input, int TRIANGLE_NUM1, the number of triangles.
//
//    Input, double NODE_XY1[2*NODE_NUM1], the nodes.
//
//    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the nodes that make up the
//    triangles.  
//
//    Input, int NODE_NUM2, the number of nodes in the refined mesh.
//
//    Input, int TRIANGLE_NUM2, the number of triangles in the refined mesh.
//
//    Input, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge information computed
//    by TRIANGULATION_ORDER6_REFINE_SIZE.
//
//    Output, double NODE_XY2[2*NODE_NUM2], the refined nodes.
//
//    Output, int TRIANGLE_NODE2[6*TRIANGLE_NUM2], the nodes that make up the
//    triangles in the refined mesh.
//
{
  int edge;
  int i;
  int j;
  int l1;
  int l2;
  int l3;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int node;
  int t1;
  int t2;
  int t3;
  int t4;
  int triangle1;
  int v1;
  int v2;
  int v3;
  int v4;
  int v5;
  int v6;
//
//  Step 1: 
//  Copy the old nodes.
//
  for ( j = 0; j < node_num1; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      node_xy2[i+j*2] = node_xy1[i+j*2];
    }
  }
  for ( j = 0; j < triangle_num2; j++ )
  {
    for ( i = 0; i < 6; i++ )
    {
      triangle_node2[i+j*6] = -1;
    }
  }
//
//  We can assign the existing nodes to the new triangles.
//
  for ( triangle1 = 0; triangle1 < triangle_num1; triangle1++ )
  {
    t1 = triangle1 * 4 + 0;
    t2 = triangle1 * 4 + 1;
    t3 = triangle1 * 4 + 2;
    t4 = triangle1 * 4 + 3;

    triangle_node2[0+t1*6] = triangle_node1[0+triangle1*6];
    triangle_node2[1+t1*6] = triangle_node1[3+triangle1*6];
    triangle_node2[2+t1*6] = triangle_node1[5+triangle1*6];

    triangle_node2[0+t2*6] = triangle_node1[3+triangle1*6];
    triangle_node2[1+t2*6] = triangle_node1[1+triangle1*6];
    triangle_node2[2+t2*6] = triangle_node1[4+triangle1*6];

    triangle_node2[0+t3*6] = triangle_node1[5+triangle1*6];
    triangle_node2[1+t3*6] = triangle_node1[4+triangle1*6];
    triangle_node2[2+t3*6] = triangle_node1[2+triangle1*6];

    triangle_node2[0+t4*6] = triangle_node1[4+triangle1*6];
    triangle_node2[1+t4*6] = triangle_node1[5+triangle1*6];
    triangle_node2[2+t4*6] = triangle_node1[3+triangle1*6];
  }
//
//  Step 2.
//  Examine sorted edge information.  The first time an edge is encountered,
//  generate two new nodes, then assign them (usually) to the four subtriangles 
//  of the two triangles that share that edge.
//
  node = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 3 * triangle_num1; edge++ )
  {
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];

    l1 = edge_data[2+edge*5];
    l3 = edge_data[3+edge*5];

    if ( l1 == 1 && l3 == 2 )
    {
      l2 = 4;
    }
    else if ( l1 == 1 && l3 == 3 )
    {
      l2 = 6;
    }
    else if ( l1 == 2 && l3 == 3 )
    {
      l2 = 5;
    }
    triangle1 = edge_data[4+edge*5];
//
//  If this is the first time we've encountered this edge,
//  create the new nodes.
//
    if ( n1 != n1_old || n2 != n2_old )
    {
      n1_old = n1;
      n2_old = n2;

      v1 = triangle_node1[l1-1+triangle1*6];
      v2 = triangle_node1[l2-1+triangle1*6];
      v3 = triangle_node1[l3-1+triangle1*6];

      for ( i = 0; i < 2; i++ )
      {
        node_xy2[i+node*2] = ( node_xy2[i+(v1-1)*2] 
                             + node_xy2[i+(v2-1)*2] ) / 2.0;
      }
      node = node + 1;
      v4 = node;

      for ( i = 0; i < 2; i++ )
      {
        node_xy2[i+node*2] = ( node_xy2[i+(v2-1)*2] 
                             + node_xy2[i+(v3-1)*2] ) / 2.0;
      }
      node = node + 1;
      v5 = node;
    }
    t1 = triangle1 * 4 + 0;
    t2 = triangle1 * 4 + 1;
    t3 = triangle1 * 4 + 2;

    if ( l1 == 1 && l3 == 2 )
    {
      if ( triangle_node1[0+triangle1*6] == v1 + 1 )
      {
        triangle_node2[3+t1*6] = v4;
        triangle_node2[3+t2*6] = v5;
      }
      else
      {
        triangle_node2[3+t1*6] = v5;
        triangle_node2[3+t2*6] = v4;
      }
    }
    else if ( l1 == 1 && l3 == 3 )
    {
      if ( triangle_node1[0+triangle1*6] == v1 + 1 )
      {
        triangle_node2[5+t1*6] = v4;
        triangle_node2[5+t3*6] = v5;
      }
      else
      {
        triangle_node2[5+t1*6] = v5;
        triangle_node2[5+t3*6] = v4;
      }
    }
    else if ( l1 == 2 && l3 == 3 )
    {
      if ( triangle_node1[1+triangle1*6] == v1 + 1 )
      {
        triangle_node2[4+t3*6] = v4;
        triangle_node2[4+t2*6] = v5;
      }
      else
      {
        triangle_node2[4+t3*6] = v5;
        triangle_node2[4+t2*6] = v4;
      }
    }
  }
//
//  Step 3.
//  Each old triangle has a single central subtriangle, for which we now
//  need to generate three new "interior" nodes.
//
  for ( triangle1 = 0; triangle1 < triangle_num1; triangle1++ )
  {
    v4 = triangle_node1[3+triangle1*6];
    v5 = triangle_node1[4+triangle1*6];
    v6 = triangle_node1[5+triangle1*6];

    t1 = triangle1 * 4 + 0;
    t2 = triangle1 * 4 + 1;
    t3 = triangle1 * 4 + 2;
    t4 = triangle1 * 4 + 3;

    node_xy2[0+node*2] = 0.5 * ( node_xy1[0+(v5-1)*2] + node_xy1[0+(v6-1)*2] );
    node_xy2[1+node*2] = 0.5 * ( node_xy1[1+(v5-1)*2] + node_xy1[1+(v6-1)*2] );
    node = node + 1;
    triangle_node2[3+t4*6] = node;
    triangle_node2[3+t3*6] = node;

    node_xy2[0+node*2] = 0.5 * ( node_xy1[0+(v6-1)*2] + node_xy1[0+(v4-1)*2] );
    node_xy2[1+node*2] = 0.5 * ( node_xy1[1+(v6-1)*2] + node_xy1[1+(v4-1)*2] );
    node = node + 1;
    triangle_node2[4+t4*6] = node;
    triangle_node2[4+t1*6] = node;

    node_xy2[0+node*2] = 0.5 * ( node_xy1[0+(v4-1)*2] + node_xy1[0+(v5-1)*2] );
    node_xy2[1+node*2] = 0.5 * ( node_xy1[1+(v4-1)*2] + node_xy1[1+(v5-1)*2] );
    node = node + 1;
    triangle_node2[5+t4*6] = node;
    triangle_node2[5+t2*6] = node; 
  }

  return;
}
//****************************************************************************80

void triangulation_order6_refine_size ( int node_num1, int triangle_num1, 
  int triangle_node1[], int *node_num2, int *triangle_num2, int edge_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_REFINE_SIZE sizes a refined order 6 triangulation.
//
//  Discussion:
//
//    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we 
//    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
//    quadratic subtriangles T1, T2, T3 and T4.
//
//    The task is more complicated by the fact that we are working with
//    a mesh of triangles, so that we want to create a node only once,
//    even though it may be shared by other triangles.  (In fact, only
//    the new nodes on the edges can be shared, and then only by at most
//    one other triangle.)
//
//            3
//           / \
//          36 35
//         / T3  \
//        6--56---5
//       / \ T4  / \
//      16 46  45  25
//     / T1  \ / T2  \
//    1--14---4--24---2
//
//    This routine determines the sizes of the resulting node and
//    triangles, and constructs an edge array that can be used to 
//    properly number the new nodes.
//
//    The primary work occurs in sorting a list related to the edges.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes in the original mesh.
//
//    Input, int  TRIANGLE_NUM1, the number of triangles in the
//    original mesh.
//
//    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the indices of the nodes
//    that form the triangles in the input mesh.
//
//    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
//
//    Output, int *TRIANGLE_NUM2, the number of triangles in the
//    refined mesh.
//
//    Output, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge data that will
//    be needed by TRIANGULATION_ORDER6_REFINE_COMPUTE.
//
{
  int a;
  int b;
  int edge;
  int i;
  int j;
  int k;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int triangle1;
//
//  Step 1.
//  From the list of nodes for triangle T, of the form: (I,J,K)
//  construct the edge relations:
//
//    (I,J,1,2,T)
//    (I,K,1,3,T)
//    (J,K,2,3,T)
//
//  In order to make matching easier, we reorder each pair of nodes
//  into ascending order.
//
  for ( triangle1 = 0; triangle1 < triangle_num1; triangle1++ )
  {
    i = triangle_node1[0+triangle1*6];
    j = triangle_node1[1+triangle1*6];
    k = triangle_node1[2+triangle1*6];

    a = i4_min ( i, j );
    b = i4_max ( i, j );

    edge_data[0+5*(3*triangle1+0)] = a;
    edge_data[1+5*(3*triangle1+0)] = b;
    edge_data[2+5*(3*triangle1+0)] = 1;
    edge_data[3+5*(3*triangle1+0)] = 2;
    edge_data[4+5*(3*triangle1+0)] = triangle1;

    a = i4_min ( i, k );
    b = i4_max ( i, k );

    edge_data[0+5*(3*triangle1+1)] = a;
    edge_data[1+5*(3*triangle1+1)] = b;
    edge_data[2+5*(3*triangle1+1)] = 1;
    edge_data[3+5*(3*triangle1+1)] = 3;
    edge_data[4+5*(3*triangle1+1)] = triangle1;

    a = i4_min ( j, k );
    b = i4_max ( j, k );

    edge_data[0+5*(3*triangle1+2)] = a;
    edge_data[1+5*(3*triangle1+2)] = b;
    edge_data[2+5*(3*triangle1+2)] = 2;
    edge_data[3+5*(3*triangle1+2)] = 3;
    edge_data[4+5*(3*triangle1+2)] = triangle1;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:2; the routine we call here
//  sorts on the full column but that won't hurt us.
//
//  What we need is to find all cases where triangles share an edge.
//  By sorting the columns of the EDGE_DATA array, we will put shared edges
//  next to each other.
//
  i4col_sort_a ( 5, 3*triangle_num1, edge_data );
//
//  Step 3. All the triangles which share an edge show up as consecutive
//  columns with identical first two entries.  Figure out how many new
//  nodes there are, and allocate space for their coordinates.
//
  *node_num2 = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 3 * triangle_num1; edge++ )
  {
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];
    if ( n1 != n1_old || n2 != n2_old )
    {
      *node_num2 = *node_num2 + 2;
      n1_old = n1;
      n2_old = n2;
    }
  }

  *node_num2 = *node_num2 + 3 * triangle_num1;

  *triangle_num2 = 4 * triangle_num1;

  return;
}
