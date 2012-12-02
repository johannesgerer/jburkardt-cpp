# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
double arc_cosine ( double c );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
double *dtable_data_read ( string input_filename, int m, int n );
void dtable_header_read ( string input_filename, int *m, int *n );
int file_column_count ( string filename );
int file_row_count ( string input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
bool line_exp_is_degenerate_nd ( int dim_num, double p1[], double p2[] );
double *line_exp_perp_2d ( double p1[2], double p2[2], double p3[2], bool *flag );
void line_exp2imp_2d ( double p1[2], double p2[2], double *a, double *b, 
  double *c );
bool line_imp_is_degenerate_2d ( double a, double b, double c );
void lines_exp_int_2d ( double p1[2], double p2[2], double p3[2], double p4[2], 
  int *ival, double p[2] );
void lines_imp_int_2d ( double a1, double b1, double c1, double a2, double b2, 
  double c2, int *ival, double p[2] );
double r8_abs ( double x );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double *r8mat_inverse_2d ( double a[] );
int r8mat_solve ( int n, int rhs_num, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void r8vec_copy ( int n, double a1[], double a2[] );
bool r8vec_eq ( int n, double a1[], double a2[] );
double r8vec_length ( int dim_num, double x[] );
void r8vec_print ( int n, double a[], string title );
int s_len_trim ( string s );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void timestamp ( );
void triangle_angles_2d ( double t[2*3], double angle[3] );
double triangle_area_2d ( double t[2*3] );
double *triangle_centroid_2d ( double t[2*3] );
void triangle_circumcircle_2d ( double t[2*3], double *r, double pc[2] );
double *triangle_edge_length_2d ( double t[2*3] );
void triangle_incircle_2d ( double t[2*3], double pc[2], double *r );
int triangle_orientation_2d ( double t[2*3] );
void triangle_orthocenter_2d ( double t[2*3], double p[2], bool *flag );
double triangle_quality_2d ( double t[2*3] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_PROPERTIES.
//
//  Discussion:
//
//    TRIANGLE_PROPERTIES reports properties of a triangle.
//
//  Usage:
//
//    triangle_properties filename
//
//    where "filename" is a file containing the coordinates of the vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
{
  double angles[3];
  double area;
  double *centroid;
  double circum_center[2];
  double circum_radius;
  int dim_num;
  double *edge_length;
  bool flag;
  int i;
  int j;
  double in_center[2];
  double in_radius;
  string node_filename;
  int node_num;
  double *node_xy;
  int orientation;
  double ortho_center[2];
  double pi = 3.141592653589793;
  double quality;

  cout << "\n";
  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_PROPERTIES:\n";
  cout << "  C++ version:\n";
  cout << "  Determine properties of a triangle.\n";

  if ( 1 < argc )
  {
    node_filename = argv[1];
  }
  else
  {
    cout << "\n";
    cout << "  Please enter the name of the node coordinate file.\n";
    cin >> node_filename;
  }
//
//  Read the node data.
//
  dtable_header_read ( node_filename, &dim_num, &node_num );

  cout << "\n";
  cout << "  Read the header of \"" << node_filename << "\".\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  Number of points NODE_NUM = " << node_num << "\n";

  if ( dim_num != 2 )
  {
    cout << "\n";
    cout << "TRIANGLE_PROPERTIES - Fatal error!\n";
    cout << "  Dataset must have spatial dimension 2.\n";
    exit ( 1 );
  }

  if ( node_num != 3 )
  {
    cout << "\n";
    cout << "TRIANGLE_PROPERTIES - Fatal error!\n";
    cout << "  Dataset must have 3 nodes.\n";
    exit ( 1 );
  }

  node_xy = dtable_data_read ( node_filename, dim_num, node_num );

  cout << "\n";
  cout << "  Read the data in \"" << node_filename << "\".\n";

  r8mat_transpose_print ( dim_num, node_num, node_xy, "  Node coordinates:" );
//
//  ANGLES
//
  triangle_angles_2d ( node_xy, angles );

  r8vec_print ( 3, angles, "  ANGLES (radians):" );

  for ( i = 0; i < 3; i++ )
  {
    angles[i] = angles[i] * 180.0 / pi;
  }
  r8vec_print ( 3, angles, "  ANGLES (degrees):" );
//
//  AREA
//
  area = triangle_area_2d ( node_xy );

  cout << "\n";
  cout << "  AREA: " << area << "\n";
//
//  CENTROID
//
  centroid = triangle_centroid_2d ( node_xy );

  cout << "\n";
  cout << "  CENTROID: "
       << "  " << setw(14) << centroid[0]
       << "  " << setw(14) << centroid[1] << "\n";
//
//  CIRCUM_CIRCLE
//
  triangle_circumcircle_2d ( node_xy, &circum_radius, circum_center );

  cout << "\n";
  cout << "  CIRCUM_RADIUS: " << circum_radius << "\n";
  cout << "  CIRCUM_CENTER: " 
       << "  " << setw(14) << circum_center[0]
       << "  " << setw(14) << circum_center[1] << "\n";
//
//  EDGE LENGTHS
//
  edge_length = triangle_edge_length_2d ( node_xy );

  r8vec_print ( 3, edge_length, "  EDGE_LENGTHS:" );
//
//  IN_CIRCLE
//
  triangle_incircle_2d ( node_xy, &in_radius, in_center );

  cout << "\n";
  cout << "  IN_RADIUS: " << in_radius << "\n";
  cout << "  IN_CENTER: " 
       << "  " << setw(14) << in_center[0]
       << "  " << setw(14) << in_center[1] << "\n";
//
//  ORIENTATION
//
  orientation = triangle_orientation_2d ( node_xy );

  cout << "\n";
  if ( orientation == 0 )
  {
    cout << "  ORIENTATION: CounterClockwise.\n";
  }
  else if ( orientation == 1 )
  {
    cout << "  ORIENTATION: Clockwise.\n";
  }
  else if ( orientation == 2 )
  {
    cout << "  ORIENTATION: Degenerate Distinct Colinear Points.\n";
  }
  else if ( orientation == 3 )
  {
    cout << "  ORIENTATION: Degenerate, at least two points identical.\n";
  }
//
//  ORTHO_CENTER
//
  triangle_orthocenter_2d ( node_xy, ortho_center, &flag );

  if ( flag )
  {
    cout << "\n";
    cout << "  ORTHO_CENTER: Could not be computed.\n";
  }
  else
  {
    cout << "\n";
    cout << "  ORTHO_CENTER: " 
         << "  " << setw(14) << ortho_center[0]
         << "  " << setw(14) << ortho_center[1] << "\n";
  }
//
//  QUALITY
//
  quality = triangle_quality_2d ( node_xy );

  cout << "\n";
  cout << "  QUALITY: " << quality << "\n";

  delete [] centroid;
  delete [] edge_length;
  delete [] node_xy;

  cout << "\n";
  cout << "TRIANGLE_PROPERTIES:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

double arc_cosine ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    ARC_COSINE computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument, the cosine of an angle.
//
//    Output, double ARC_COSINE, an angle whose cosine is C.
//
{
  double angle;
  double pi = 3.141592653589793;

  if ( c <= -1.0 )
  {
    angle = pi;
  } 
  else if ( 1.0 <= c )
  {
    angle = 0.0;
  }
  else
  {
    angle = acos ( c );
  }
  return angle;
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
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the 
//    character was 'illegal', then DIGIT is -1.
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

double *dtable_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_DATA_READ reads the data from a DTABLE file.
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
//    Output, double DTABLE_DATA_READ[M*N], the table data.
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
    cerr << "DTABLE_DATA_READ - Fatal error!\n";
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
 
void dtable_header_read ( string input_filename, int *m, int *n )
 
//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_HEADER_READ reads the header from a DTABLE file.
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
    cerr << "DTABLE_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "DTABLE_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    return;
  }

  return;
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
//****************************************************************************80

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I  I4_WRAP
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

bool line_exp_is_degenerate_nd ( int dim_num, double p1[], double p2[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_EXP_IS_DEGENERATE_ND finds if an explicit line is degenerate in ND.
//
//  Discussion:
//
//    The explicit form of a line in ND is:
//
//      the line through the points P1 and P2.
//
//    An explicit line is degenerate if the two defining points are equal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double P1[DIM_NUM], P2[DIM_NUM], two points on the line.
//
//    Output, bool LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
//    is degenerate.
//
{
  bool value;

  value = r8vec_eq ( dim_num, p1, p2 );

  return value;
}
//****************************************************************************80

double *line_exp_perp_2d ( double p1[2], double p2[2], double p3[2], 
  bool *flag )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_EXP_PERP_2D computes a line perpendicular to a line and through a point.
//
//  Discussion:
//
//    The explicit form of a line in 2D is:
//
//      the line through P1 and P2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double P1[2], P2[2], two points on the given line.
//
//    Input, double P3[2], a point not on the given line, through which the
//    perpendicular must pass.
//
//    Output, double LINE_EXP_PERP_2D[2], a point on the given line, such that the line
//    through P3 and P4 is perpendicular to the given line.
//
//    Output, bool *FLAG, is TRUE if the point could not be computed.
//
{
# define DIM_NUM 2

  double bot;
  double *p4;
  double t;

  p4 = new double[DIM_NUM];

  bot = pow ( p2[0] - p1[0], 2 ) + pow ( p2[1] - p1[1], 2 );

  if ( bot == 0.0 )
  {
    p4[0] = r8_huge ( );
    p4[1] = r8_huge ( );
    *flag = true;
    return p4;
  }
//
//  (P3-P1) dot (P2-P1) = Norm(P3-P1) * Norm(P2-P1) * Cos(Theta).
//
//  (P3-P1) dot (P2-P1) / Norm(P3-P1)**2 = normalized coordinate T
//  of the projection of (P3-P1) onto (P2-P1).
//
  t = ( ( p1[0] - p3[0] ) * ( p1[0] - p2[0] ) 
      + ( p1[1] - p3[1] ) * ( p1[1] - p2[1] ) ) / bot;

  p4[0] = p1[0] + t * ( p2[0] - p1[0] );
  p4[1] = p1[1] + t * ( p2[1] - p1[1] );

  *flag = false;

  return p4;
# undef DIM_NUM
}
//****************************************************************************80

void line_exp2imp_2d ( double p1[2], double p2[2], double *a, double *b, 
  double *c )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
//
//  Discussion:
//
//    The explicit form of a line in 2D is:
//
//      the line through P1 and P2
//
//    The implicit form of a line in 2D is:
//
//      A * X + B * Y + C = 0
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
//    Input, double P1[2], P2[2], two distinct points on the line. 
//
//    Output, double *A, *B, *C, three coefficients which describe
//    the line that passes through P1 and P2.
//
{
//
//  Take care of degenerate cases.
//
  if ( r8vec_eq ( 2, p1, p2 ) )
  {
    cout << "\n";
    cout << "LINE_EXP2IMP_2D - Fatal error!\n";
    cout << "  P1 = P2\n";
    cout << "  P1 = " << p1[0] << " " << p1[1] << "\n";
    cout << "  P2 = " << p2[0] << " " << p2[1] << "\n";
    exit ( 1 );
  }

  *a = p2[1] - p1[1];
  *b = p1[0] - p2[0];
  *c = p2[0] * p1[1] - p1[0] * p2[1];

  return;
}
//****************************************************************************80

bool line_imp_is_degenerate_2d ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_IMP_IS_DEGENERATE_2D finds if an implicit point is degenerate in 2D.
//
//  Discussion:
//
//    The implicit form of a line in 2D is:
//
//      A * X + B * Y + C = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the implicit line parameters.
//
//    Output, bool LINE_IMP_IS_DEGENERATE_2D, is true if the
//    line is degenerate.
//
{
  bool value;

  value = ( a * a + b * b == 0.0 );

  return value;
}
//****************************************************************************80

void lines_exp_int_2d ( double p1[2], double p2[2], double p3[2], double p4[2], 
  int *ival, double p[2] )

//****************************************************************************80
//
//  Purpose:
//
//    LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
//
//  Discussion:
//
//    The explicit form of a line in 2D is:
//
//      the line through P1 and P2.
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
//    Input, double P1[2], P2[2], define the first line.
//
//    Input, double P3[2], P4[2], define the second line.
//
//    Output, int *IVAL, reports on the intersection:
//    0, no intersection, the lines may be parallel or degenerate.
//    1, one intersection point, returned in P.
//    2, infinitely many intersections, the lines are identical.
//
//    Output, double P[2], if IVAl = 1, then P contains
//    the intersection point.  Otherwise, P = 0.
//
{
# define DIM_NUM 2

  double a1 = 0.0;
  double a2 = 0.0;
  double b1 = 0.0;
  double b2 = 0.0;
  double c1 = 0.0;
  double c2 = 0.0;
  double point_1 = 0.0;
  double point_2 = 0.0;

  *ival = 0;
  p[0] = 0.0;
  p[1] = 0.0;
//
//  Check whether either line is a point.
//
  if ( r8vec_eq ( DIM_NUM, p1, p2 ) )
  {
    point_1 = true;
  }
  else
  {
    point_1 = false;
  }

  if ( r8vec_eq ( DIM_NUM, p3, p4 ) )
  {
    point_2 = true;
  }
  else
  {
    point_2 = false;
  }
//
//  Convert the lines to ABC format.
//
  if ( !point_1 )
  {
    line_exp2imp_2d ( p1, p2, &a1, &b1, &c1 );
  }

  if ( !point_2 )
  {
    line_exp2imp_2d ( p3, p4, &a2, &b2, &c2 );
  }
//
//  Search for intersection of the lines.
//
  if ( point_1 && point_2 )
  {
    if ( r8vec_eq ( DIM_NUM, p1, p3 ) )
    {
      *ival = 1;
      r8vec_copy ( DIM_NUM, p1, p );
    }
  }
  else if ( point_1 )
  {
    if ( a2 * p1[0] + b2 * p1[1] == c2 )
    {
      *ival = 1;
      r8vec_copy ( DIM_NUM, p1, p );
    }
  }
  else if ( point_2 )
  {
    if ( a1 * p3[0] + b1 * p3[1] == c1 )
    {
      *ival = 1;
      r8vec_copy ( DIM_NUM, p3, p );
    }
  }
  else
  {
    lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, p );
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void lines_imp_int_2d ( double a1, double b1, double c1, double a2, double b2, 
  double c2, int *ival, double p[2] )

//****************************************************************************80
//
//  Purpose:
//
//    LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
//
//  Discussion:
//
//    The implicit form of a line in 2D is:
//
//      A * X + B * Y + C = 0
//
//    22 May 2004: Thanks to John Asmuth for pointing out that the 
//    B array was not being deallocated on exit.
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
//    Input, double A1, B1, C1, define the first line.
//    At least one of A1 and B1 must be nonzero.
//
//    Input, double A2, B2, C2, define the second line.
//    At least one of A2 and B2 must be nonzero.
//
//    Output, int *IVAL, reports on the intersection.
//    -1, both A1 and B1 were zero.
//    -2, both A2 and B2 were zero.
//     0, no intersection, the lines are parallel.
//     1, one intersection point, returned in P.
//     2, infinitely many intersections, the lines are identical.
//
//    Output, double P[2], if IVAL = 1, then P contains
//    the intersection point.  Otherwise, P = 0.
//
{
# define DIM_NUM 2

  double a[DIM_NUM*2];
  double *b;

  p[0] = 0.0;
  p[1] = 0.0;
//
//  Refuse to handle degenerate lines.
//
  if ( a1 == 0.0 && b1 == 0.0 )
  {
    *ival = - 1;
    return;
  }
  else if ( a2 == 0.0 && b2 == 0.0 )
  {
    *ival = - 2;
    return;
  }
//
//  Set up a linear system, and compute its inverse.
//
  a[0+0*2] = a1;
  a[0+1*2] = b1;
  a[1+0*2] = a2;
  a[1+1*2] = b2;

  b = r8mat_inverse_2d ( a );
//
//  If the inverse exists, then the lines intersect.
//  Multiply the inverse times -C to get the intersection point.
//
  if ( b != NULL )
  {

    *ival = 1;
    p[0] = - b[0+0*2] * c1 - b[0+1*2] * c2;
    p[1] = - b[1+0*2] * c1 - b[1+1*2] * c2;
  }
//
//  If the inverse does not exist, then the lines are parallel
//  or coincident.  Check for parallelism by seeing if the
//  C entries are in the same ratio as the A or B entries.
//
  else
  {

    *ival = 0;

    if ( a1 == 0.0 )
    {
      if ( b2 * c1 == c2 * b1 )
      {
        *ival = 2;
      }
    }
    else
    {
      if ( a2 * c1 == c2 * a1 )
      {
        *ival = 2;
      }
    }
  }

  delete [] b;

  return;
# undef DIM_NUM
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double *r8mat_inverse_2d ( double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_INVERSE_2D inverts a 2 by 2 R8MAT using Cramer's rule.
//
//  Discussion:
//
//    The two dimensional array is stored as a one dimensional vector,
//    by COLUMNS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[2*2], the matrix to be inverted.
//
//    Output, double R8MAT_INVERSE_2D[2*2], the inverse of the matrix A.
//
{
  double *b;
  double det;
//
//  Compute the determinant of A.
//
  det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];
//
//  If the determinant is zero, bail out.
//
  if ( det == 0.0 )
  {
    return NULL;
  }
//
//  Compute the entries of the inverse matrix using an explicit formula.
//
  b = new double[2*2];

  b[0+0*2] = + a[1+1*2] / det;
  b[0+1*2] = - a[0+1*2] / det;
  b[1+0*2] = - a[1+0*2] / det;
  b[1+1*2] = + a[0+0*2] / det;

  return b;
}
//****************************************************************************80

int r8mat_solve ( int n, int rhs_num, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*N]
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
//    Input, int N, the order of the matrix.
//
//    Input, int RHS_NUM, the number of right hand sides.  RHS_NUM
//    must be at least 0.
//
//    Input/output, double A[N*(N+RHS_NUM)], contains in rows and columns 1
//    to N the coefficient matrix, and in columns N+1 through
//    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
//    area has been destroyed, while the right hand sides have
//    been overwritten with the corresponding solutions.
//
//    Output, int R8MAT_SOLVE, singularity flag.
//    0, the matrix was not singular, the solutions were computed;
//    J, factorization failed on step J, and the solutions could not
//    be computed.
//
{
  double apivot;
  double factor;
  int i;
  int ipivot;
  int j;
  int k;
  double temp;

  for ( j = 0; j < n; j++ )
  {
//
//  Choose a pivot row.
//
    ipivot = j;
    apivot = a[j+j*n];

    for ( i = j; i < n; i++ )
    {
      if ( r8_abs ( apivot ) < r8_abs ( a[i+j*n] ) )
      {
        apivot = a[i+j*n];
        ipivot = i;
      }
    }

    if ( apivot == 0.0 )
    {
      return j;
    }
//
//  Interchange.
//
    for ( i = 0; i < n + rhs_num; i++ )
    {
      temp          = a[ipivot+i*n];
      a[ipivot+i*n] = a[j+i*n];
      a[j+i*n]      = temp;
    }
//
//  A(J,J) becomes 1.
//
    a[j+j*n] = 1.0;
    for ( k = j; k < n + rhs_num; k++ )
    {
      a[j+k*n] = a[j+k*n] / apivot;
    }
//
//  A(I,J) becomes 0.
//
    for ( i = 0; i < n; i++ )
    {
      if ( i != j )
      {
        factor = a[i+j*n];
        a[i+j*n] = 0.0;
        for ( k = j; k < n + rhs_num; k++ )
        {
          a[i+k*n] = a[i+k*n] - factor * a[j+k*n];
        }
      }
    }
  }

  return 0;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Input, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

bool r8vec_eq ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EQ is true two R8VEC's are equal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], two vectors to compare.
//
//    Output, bool R8VEC_EQ.
//    R8VEC_EQ is TRUE if every pair of elements A1(I) and A2(I) are equal,
//    and FALSE otherwise.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] != a2[i] )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

double r8vec_length ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LENGTH returns the Euclidean length of an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the vector.
//
//    Output, double R8VEC_LENGTH, the Euclidean length of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    value = value + pow ( x[i], 2 );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
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
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(14) << a[i]  << "\n";
  }

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
//****************************************************************************80

void triangle_angles_2d ( double t[2*3], double angle[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ANGLES_2D computes the angles of a triangle in 2D.
//
//  Discussion:
//
//    The law of cosines is used:
//
//      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
//
//    where GAMMA is the angle opposite side C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double ANGLE[3], the angles opposite
//    sides P1-P2, P2-P3 and P3-P1, in radians.
//
{
  double a;
  double b;
  double c;
  double pi = 3.141592653589793;

  a = sqrt ( pow ( t[0+1*2] - t[0+0*2], 2 ) 
           + pow ( t[1+1*2] - t[1+0*2], 2 ) );

  b = sqrt ( pow ( t[0+2*2] - t[0+1*2], 2 ) 
           + pow ( t[1+2*2] - t[1+1*2], 2 ) );

  c = sqrt ( pow ( t[0+0*2] - t[0+2*2], 2 ) 
           + pow ( t[1+0*2] - t[1+2*2], 2 ) );
//
//  Take care of a ridiculous special case.
//
  if ( a == 0.0 && b == 0.0 && c == 0.0 )
  {
    angle[0] = 2.0 * pi / 3.0;
    angle[1] = 2.0 * pi / 3.0;
    angle[2] = 2.0 * pi / 3.0;
    return;
  }

  if ( c == 0.0 || a == 0.0 )
  {
    angle[0] = pi;
  }
  else
  {
    angle[0] = arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0 * c * a ) );
  }

  if ( a == 0.0 || b == 0.0 )
  {
    angle[1] = pi;
  }
  else
  {
    angle[1] = arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0 * a * b ) );
  }

  if ( b == 0.0 || c == 0.0 )
  {
    angle[2] = pi;
  }
  else
  {
    angle[2] = arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0 * b * c ) );
  }

  return;
}
//****************************************************************************80

double triangle_area_2d ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
//
//  Discussion:
//
//    If the triangle's vertices are given in counter clockwise order,
//    the area will be positive.  If the triangle's vertices are given
//    in clockwise order, the area will be negative!
//
//    An earlier version of this routine always returned the absolute
//    value of the computed area.  I am convinced now that that is
//    a less useful result!  For instance, by returning the signed 
//    area of a triangle, it is possible to easily compute the area 
//    of a nonconvex polygon as the sum of the (possibly negative) 
//    areas of triangles formed by node 1 and successive pairs of vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA_2D, the area of the triangle.
//
{
  double area;

  area = 0.5 * ( 
    t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) + 
    t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) + 
    t[0+2*2] * ( t[1+0*2] - t[1+1*2] ) );
 
  return area;
}
//****************************************************************************80

double *triangle_centroid_2d ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CENTROID_2D computes the centroid of a triangle in 2D.
//
//  Discussion:
//
//    The centroid of a triangle can also be considered the center
//    of gravity, assuming that the triangle is made of a thin uniform
//    sheet of massy material.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_CENTROID_2D[2], the coordinates of the centroid of the triangle.
//
{
# define DIM_NUM 2

  double *centroid;

  centroid = new double[DIM_NUM];

  centroid[0] = ( t[0+0*DIM_NUM] + t[0+1*DIM_NUM] + t[0+2*DIM_NUM] ) / 3.0;
  centroid[1] = ( t[1+0*DIM_NUM] + t[1+1*DIM_NUM] + t[1+2*DIM_NUM] ) / 3.0;
 
  return centroid;
# undef DIM_NUM
}
//****************************************************************************80

void triangle_circumcircle_2d ( double t[2*3], double *r, double pc[2] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CIRCUMCIRCLE_2D computes the circumcircle of a triangle in 2D.
//
//  Discussion:
//
//    The circumcenter of a triangle is the center of the circumcircle, the
//    circle that passes through the three vertices of the triangle.
//
//    The circumcircle contains the triangle, but it is not necessarily the
//    smallest triangle to do so.
//
//    If all angles of the triangle are no greater than 90 degrees, then
//    the center of the circumscribed circle will lie inside the triangle.
//    Otherwise, the center will lie outside the triangle.
//
//    The circumcenter is the intersection of the perpendicular bisectors
//    of the sides of the triangle.
//
//    In geometry, the circumcenter of a triangle is often symbolized by "O".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double *R, PC[2], the circumradius, and the coordinates of the 
//    circumcenter of the triangle.
//
{
# define DIM_NUM 2

  double a;
  double b;
  double bot;
  double c;
  double top1;
  double top2;
//
//  Circumradius.
//
  a = sqrt ( pow ( t[0+1*2] - t[0+0*2], 2 ) + pow ( t[1+1*2] - t[1+0*2], 2 ) );
  b = sqrt ( pow ( t[0+2*2] - t[0+1*2], 2 ) + pow ( t[1+2*2] - t[1+1*2], 2 ) );
  c = sqrt ( pow ( t[0+0*2] - t[0+2*2], 2 ) + pow ( t[1+0*2] - t[1+2*2], 2 ) );

  bot = ( a + b + c ) * ( - a + b + c ) * (   a - b + c ) * (   a + b - c );

  if ( bot <= 0.0 )
  {
    *r = -1.0;
    pc[0] = 0.0;
    pc[1] = 0.0;
    return;
  }

  *r = a * b * c / sqrt ( bot );
//
//  Circumcenter.
//
  top1 =  ( t[1+1*2] - t[1+0*2] ) * c * c - ( t[1+2*2] - t[1+0*2] ) * a * a;
  top2 =  ( t[0+1*2] - t[0+0*2] ) * c * c - ( t[0+2*2] - t[0+0*2] ) * a * a;
  bot  =  ( t[1+1*2] - t[1+0*2] ) * ( t[0+2*2] - t[0+0*2] )  
        - ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2] - t[0+0*2] );

  pc[0] = t[0+0*2] + 0.5 * top1 / bot;
  pc[1] = t[1+0*2] - 0.5 * top2 / bot;

  return;
# undef DIM_NUM
}
//****************************************************************************80

double *triangle_edge_length_2d ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_EDGE_LENGTH_2D returns edge lengths of a triangle in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double TRIANGLE_EDGE_LENGTH[3], the length of the edges.
//
{
  double *edge_length;
  int j1;
  int j2;

  edge_length = new double[3];

  for ( j1 = 0; j1 < 3; j1++ )
  {
    j2 = i4_wrap ( j1 + 1, 0, 2 );
    edge_length[j1] = sqrt ( pow ( t[0+j2*2] - t[0+j1*2], 2 ) 
                           + pow ( t[1+j2*2] - t[1+j1*2], 2 ) );
  }

  return edge_length;
}
//****************************************************************************80

void triangle_incircle_2d ( double t[2*3], double pc[2], double *r )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_INCIRCLE_2D computes the inscribed circle of a triangle in 2D.
//
//  Discussion:
//
//    The inscribed circle of a triangle is the largest circle that can
//    be drawn inside the triangle.  It is tangent to all three sides,
//    and the lines from its center to the vertices bisect the angles
//    made by each vertex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double PC[2], *R, the center of the inscribed circle, and its radius.
//
{
# define DIM_NUM 2

  double perim;
  double s12;
  double s23;
  double s31;

  s12 = sqrt ( pow ( t[0+1*2] - t[0+0*2], 2 ) 
             + pow ( t[1+1*2] - t[1+0*2], 2 ) );
  s23 = sqrt ( pow ( t[0+2*2] - t[0+1*2], 2 ) 
             + pow ( t[1+2*2] - t[1+1*2], 2 ) );
  s31 = sqrt ( pow ( t[0+0*2] - t[0+2*2], 2 ) 
             + pow ( t[1+0*2] - t[1+2*2], 2 ) );

  perim = s12 + s23 + s31;

  if ( perim == 0.0 )
  {
    *r = 0.0;
    pc[0] = t[0+0*2];
    pc[1] = t[1+0*2];
  }
  else
  {
    pc[0] = ( s23 * t[0+0*2] + s31 * t[0+1*2] + s12 * t[0+2*2] ) / perim;
    pc[1] = ( s23 * t[1+0*2] + s31 * t[1+1*2] + s12 * t[1+2*2] ) / perim;

    *r = 0.5 * sqrt (
        ( - s12 + s23 + s31 )
      * ( + s12 - s23 + s31 )
      * ( + s12 + s23 - s31 ) / perim );
  }
  return;
# undef DIM_NUM
}
//****************************************************************************80

int triangle_orientation_2d ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
//
//  Discussion:
//
//    Three distinct non-colinear points in the plane define a circle.
//    If the points are visited in the order (x1,y1), (x2,y2), and then
//    (x3,y3), this motion defines a clockwise or counter clockwise
//    rotation along the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, int TRIANGLE_ORIENTATION_2D, reports if the three points lie
//    clockwise on the circle that passes through them.  The possible
//    return values are:
//    0, the points are distinct, noncolinear, and lie counter clockwise
//    on their circle.
//    1, the points are distinct, noncolinear, and lie clockwise
//    on their circle.
//    2, the points are distinct and colinear.
//    3, at least two of the points are identical.
//
{
# define DIM_NUM 2

  double det;
  int value = 0;

  if ( r8vec_eq ( 2, t+0*2, t+1*2 ) || 
       r8vec_eq ( 2, t+1*2, t+2*2 ) || 
       r8vec_eq ( 2, t+2*2, t+0*2 ) )
  {
    value = 3;
    return value;
  }

  det = ( t[0+0*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] ) 
      - ( t[0+1*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] );

  if ( det == 0.0 )
  {
    value = 2;
  }
  else if ( det < 0.0 )
  {
    value = 1;
  }
  else if ( 0.0 < det )
  {
    value = 0;
  }
  return value;
# undef DIM_NUM
}
//****************************************************************************80

void triangle_orthocenter_2d ( double t[2*3], double p[2], bool *flag )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ORTHOCENTER_2D computes the orthocenter of a triangle in 2D.
//
//  Discussion:
//
//    The orthocenter is defined as the intersection of the three altitudes
//    of a triangle.
//
//    An altitude of a triangle is the line through a vertex of the triangle
//    and perpendicular to the opposite side.
//
//    In geometry, the orthocenter of a triangle is often symbolized by "H".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double P[2], the coordinates of the orthocenter of the triangle.
//
//    Output, bool *FLAG, is TRUE if the point could not be computed.
//
{
# define DIM_NUM 2

  int ival;
  double *p23;
  double *p31;
//
//  Determine a point P23 common to the line through P2 and P3 and
//  its perpendicular through P1.
//
  p23 = line_exp_perp_2d ( t+1*2, t+2*2, t+0*2, flag );

  if ( *flag )
  {
    p[0] = r8_huge ( );
    p[1] = r8_huge ( );
    delete [] p23;
    return;
  }
//
//  Determine a point P31 common to the line through P3 and P1 and
//  its perpendicular through P2.
//
  p31 = line_exp_perp_2d ( t+2*2, t+0*2, t+1*2, flag );
  if ( *flag )
  {
    p[0] = r8_huge ( );
    p[1] = r8_huge ( );
    delete [] p23;
    delete [] p31;
    return;
  }
//
//  Determine P, the intersection of the lines through P1 and P23, and
//  through P2 and P31.
//
  lines_exp_int_2d ( t+0*2, p23, t+1*2, p31, &ival, p );

  if ( ival != 1 )
  {
    p[0] = r8_huge ( );
    p[1] = r8_huge ( );
    *flag = true;
    delete [] p23;
    delete [] p31;
    return;
  }
  delete [] p23;
  delete [] p31;

  return;
# undef DIM_NUM
}
//****************************************************************************80

double triangle_quality_2d ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_QUALITY_2D: "quality" of a triangle in 2D.
//
//  Discussion:
//
//    The quality of a triangle is 2 times the ratio of the radius of the inscribed
//    circle divided by that of the circumscribed circle.  An equilateral
//    triangle achieves the maximum possible quality of 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double TRIANGLE_QUALITY_2D, the quality of the triangle.
//
{
# define DIM_NUM 2

  double a;
  double b;
  double c;
  int i;
  double value;
//
//  Compute the length of each side.
//
  a = 0.0;
  b = 0.0;
  c = 0.0;

  for ( i = 0; i < DIM_NUM; i++ )
  {
    a = a + pow ( t[i+0*DIM_NUM] - t[i+1*DIM_NUM], 2 );
    b = b + pow ( t[i+1*DIM_NUM] - t[i+2*DIM_NUM], 2 );
    c = c + pow ( t[i+2*DIM_NUM] - t[i+0*DIM_NUM], 2 );
  }
  a = sqrt ( a );
  b = sqrt ( b );
  c = sqrt ( c );

  if ( a * b * c == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = ( - a + b + c ) * ( a - b + c ) * ( a + b - c ) 
      / ( a * b * c );
  }
  return value;
# undef DIM_NUM
}
