# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <string>

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
double r8_abs ( double x );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void r8_swap ( double *x, double *y );
double r8mat_det_4d ( double a[] );
int r8mat_solve ( int n, int rhs_num, double a[] );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
double r8vec_angle_3d ( double u[], double v[] );
double *r8vec_cross_3d ( double v1[3], double v2[3] );
double r8vec_dot ( int n, double a1[], double a2[] );
double r8vec_length ( int dim_num, double x[] );
double r8vec_max ( int n, double r8vec[] );
void r8vec_print ( int n, double a[], string title );
void r8vec_zero ( int n, double a[] );
int s_len_trim ( string s );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
double *tetrahedron_centroid_3d ( double tetra[3*4] );
void tetrahedron_circumsphere_3d ( double tetra[3*4], double *r, double pc[3] );
double *tetrahedron_dihedral_angles_3d ( double tetra[] );
double *tetrahedron_edge_length_3d ( double tetra[3*4] );
void tetrahedron_face_angles_3d ( double tetra[], double angles[] );
void tetrahedron_face_areas_3d ( double tetra[], double areas[] );
void tetrahedron_insphere_3d ( double tetra[3*4], double *r, double pc[3] );
double tetrahedron_quality1_3d ( double tetra[3*4] );
double tetrahedron_quality2_3d ( double tetra[3*4] );
double tetrahedron_quality3_3d ( double tetra[3*4] );
double tetrahedron_quality4_3d ( double tetra[3*4] );
double *tetrahedron_solid_angles_3d ( double tetra[] );
double tetrahedron_volume_3d ( double tetra[3*4] );
void timestamp ( );
void triangle_angles_3d ( double t[3*3], double angle[3] );
double triangle_area_3d ( double t[3*3] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TETRAHEDRON_PROPERTIES.
//
//
//  Discussion:
//
//    TETRAHEDRON_PROPERTIES reports properties of a tetrahedron.
//
//  Usage:
//
//    tetrahedron_properties filename
//
//    where "filename" is a file containing the coordinates of the vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
{
  double *centroid;
  double circum_center[3];
  double circum_radius;
  double *dihedral_angles;
  int dim_num;
  double *edge_length;
  double face_angles[3*4];
  double face_areas[4];
  int i;
  int j;
  double in_center[3];
  double in_radius;
  string node_filename;
  int node_num;
  double *node_xyz;
  double pi = 3.141592653589793;
  double quality1;
  double quality2;
  double quality3;
  double quality4;
  double *solid_angles;
  double volume;

  cout << "\n";
  timestamp ( );

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

  cout << "\n";
  cout << "TETRAHEDRON_PROPERTIES:\n";
  cout << "  C++ version:\n";
  cout << "  Determine properties of a tetrahedron.\n";

  if ( dim_num != 3 )
  {
    cout << "\n";
    cout << "TETRAHEDRON_PROPERTIES - Fatal error!\n";
    cout << "  Dataset must have spatial dimension 3.\n";
    exit ( 1 );
  }

  if ( node_num != 4 )
  {
    cout << "\n";
    cout << "TETRAHEDRON_PROPERTIES - Fatal error!\n";
    cout << "  Dataset must have 4 nodes.\n";
    exit ( 1 );
  }

  node_xyz = dtable_data_read ( node_filename, dim_num, node_num );

  cout << "\n";
  cout << "  Read the data in \"" << node_filename << "\".\n";

  r8mat_transpose_print ( dim_num, node_num, node_xyz, "  Node coordinates:" );
//
//  CIRCUMSPHERE
//
  tetrahedron_circumsphere_3d ( node_xyz, &circum_radius, circum_center );

  cout << "\n";
  cout << "  CIRCUM_RADIUS = " << circum_radius << "\n";
  cout << "  CIRCUM_CENTER: " 
       << "  " << setw(14) << circum_center[0]
       << "  " << setw(14) << circum_center[1]
       << "  " << setw(14) << circum_center[2] << "\n";
//
//  CENTROID
//
  centroid = tetrahedron_centroid_3d ( node_xyz );

  cout << "\n";
  cout << "  CENTROID: "
       << "  " << setw(14) << centroid[0]
       << "  " << setw(14) << centroid[1]
       << "  " << setw(14) << centroid[2] << "\n";
//
//  DIHEDRAL ANGLES
//
  dihedral_angles = tetrahedron_dihedral_angles_3d ( node_xyz );

  r8vec_print ( 6, dihedral_angles, "  DIHEDRAL_ANGLES (radians)" );

  for ( i = 0; i < 6; i++ )
  {
    dihedral_angles[i] = dihedral_angles[i] * 180.0 / pi;
  }

  r8vec_print ( 6, dihedral_angles, "  DIHEDRAL_ANGLES (degrees)" );
//
//  EDGE LENGTHS
//
  edge_length = tetrahedron_edge_length_3d ( node_xyz );

  r8vec_print ( 6, edge_length, "  EDGE_LENGTHS" );
//
//  FACE ANGLES
//
  tetrahedron_face_angles_3d ( node_xyz, face_angles );

  r8mat_transpose_print ( 3, 4, face_angles, "  FACE_ANGLES (radians)" );

  for ( j = 0; j < 4; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      face_angles[i+j*3] = face_angles[i+j*3] * 180.0 / pi;
    }
  }
  r8mat_transpose_print ( 3, 4, face_angles, "  FACE_ANGLES (degrees)" );
//
//  FACE AREAS
//
  tetrahedron_face_areas_3d ( node_xyz, face_areas );

  r8vec_print ( 4, face_areas, "  FACE_AREAS" );
//
//  INSPHERE
//
  tetrahedron_insphere_3d ( node_xyz, &in_radius, in_center );

  cout << "\n";
  cout << "  IN_RADIUS = " << in_radius << "\n";
  cout << "  IN_CENTER: " 
       << "  " << setw(14) << in_center[0]
       << "  " << setw(14) << in_center[1]
       << "  " << setw(14) << in_center[2] << "\n";
//
//  QUALITY1
//
  quality1 = tetrahedron_quality1_3d ( node_xyz );

  cout << "\n";
  cout << "  QUALITY1 = " << quality1 << "\n";
//
//  QUALITY2
//
  quality2 = tetrahedron_quality2_3d ( node_xyz );

  cout << "\n";
  cout << "  QUALITY2 = " << quality2 << "\n";
//
//  QUALITY3
//
  quality3 = tetrahedron_quality3_3d ( node_xyz );

  cout << "\n";
  cout << "  QUALITY3 = " << quality3 << "\n";
//
//  QUALITY4
//
  quality4 = tetrahedron_quality4_3d ( node_xyz );

  cout << "\n";
  cout << "  QUALITY4 = " << quality4 << "\n";
//
//  SOLID ANGLES
//
  solid_angles = tetrahedron_solid_angles_3d ( node_xyz );

  r8vec_print ( 4, solid_angles, "  SOLID_ANGLES (steradians)" );
//
//  VOLUME
//
  volume = tetrahedron_volume_3d ( node_xyz );

  cout << "\n";
  cout << "  VOLUME = " << volume << "\n";

  delete [] centroid;
  delete [] dihedral_angles;
  delete [] edge_length;
  delete [] node_xyz;
  delete [] solid_angles;

  cout << "\n";
  cout << "TETRAHEDRON_PROPERTIES:\n";
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

void r8_swap ( double *x, double *y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP switches two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, double *X, *Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  double z;

  z = *x;
  *x = *y;
  *y = z;
 
  return;
}
//****************************************************************************80

double r8mat_det_4d ( double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
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
//    10 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[4*4], the matrix whose determinant is desired.
//
//    Output, double R8MAT_DET_4D, the determinant of the matrix.
//
{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
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

double r8vec_angle_3d ( double u[], double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ANGLE_3D computes the angle between two vectors in 3D.
//
//  Modified:
//
//    07 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double U[3], V[3], the vectors.
//
//    Output, double ANGLE, the angle between the two vectors.
//
{
  double angle;
  double angle_cos;
  double u_norm;
  double uv_dot;
  double v_norm;

  uv_dot = r8vec_dot ( 3, u, v );

  u_norm = sqrt ( r8vec_dot ( 3, u, u ) );

  v_norm = sqrt ( r8vec_dot ( 3, v, v ) );

  angle_cos = uv_dot / u_norm / v_norm;

  angle = arc_cosine ( angle_cos );

  return angle;
}
//****************************************************************************80

double *r8vec_cross_3d ( double v1[3], double v2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CROSS_3D computes the cross product of two R8VEC's in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], the coordinates of the vectors.
//
//    Output, double R8VEC_CROSS_3D[3], the cross product vector.
//
{
# define DIM_NUM 3

  double *v3;

  v3 = new double[DIM_NUM];

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;

# undef DIM_NUM
}
//****************************************************************************80

double r8vec_dot ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's in ND.
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
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }

  return value;
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

double r8vec_max ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
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
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double value;

  value = - r8_huge ( );

  if ( n <= 0 ) 
  {
    return value;
  }

  for ( i = 0; i < n; i++ ) 
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
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

void r8vec_zero ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO zeroes an R8VEC.
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
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
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

double *tetrahedron_centroid_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_CENTROID_3D computes the centroid of a tetrahedron in 3D.
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
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_CENTROID_3D[3], the coordinates of the centroid.
//
{
# define DIM_NUM 3

  double *centroid;

  centroid = new double[3];

  centroid[0] = 0.25 * ( tetra[0+0*DIM_NUM] + tetra[0+1*DIM_NUM] 
                       + tetra[0+2*DIM_NUM] + tetra[0+3*DIM_NUM] );
  centroid[1] = 0.25 * ( tetra[1+0*DIM_NUM] + tetra[1+1*DIM_NUM] 
                       + tetra[1+2*DIM_NUM] + tetra[1+3*DIM_NUM] );
  centroid[2] = 0.25 * ( tetra[2+0*DIM_NUM] + tetra[2+1*DIM_NUM] 
                       + tetra[2+2*DIM_NUM] + tetra[2+3*DIM_NUM] );

  return centroid;
# undef DIM_NUM
}
//****************************************************************************80

void tetrahedron_circumsphere_3d ( double tetra[3*4], double *r, double pc[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere of a tetrahedron in 3D.
//
//  Discussion:
//
//    The circumsphere, or circumscribed sphere, of a tetrahedron is the sphere that
//    passes through the four vertices.  The circumsphere is not necessarily
//    the smallest sphere that contains the tetrahedron.
//
//    Surprisingly, the diameter of the sphere can be found by solving
//    a 3 by 3 linear system.  This is because the vectors P2 - P1,
//    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
//    right triangle with the diameter through P1.  Hence, the dot product of
//    P2 - P1 with that diameter is equal to the square of the length
//    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
//    the diameter vector originating at P1, and hence the radius and
//    center.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 August 2005
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
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double *R, PC[3], the coordinates of the center of the
//    circumscribed sphere, and its radius.  If the linear system is
//    singular, then R = -1, PC[] = 0.
//
{
# define DIM_NUM 3
# define RHS_NUM 1

  double a[DIM_NUM*(DIM_NUM+RHS_NUM)];
  int info;
//
//  Set up the linear system.
//
  a[0+0*3] = tetra[0+1*3] - tetra[0+0*3];
  a[0+1*3] = tetra[1+1*3] - tetra[1+0*3];
  a[0+2*3] = tetra[2+1*3] - tetra[2+0*3];
  a[0+3*3] = pow ( tetra[0+1*3] - tetra[0+0*3], 2 ) 
           + pow ( tetra[1+1*3] - tetra[1+0*3], 2 ) 
           + pow ( tetra[2+1*3] - tetra[2+0*3], 2 );

  a[1+0*3] = tetra[0+2*3] - tetra[0+0*3];
  a[1+1*3] = tetra[1+2*3] - tetra[1+0*3];
  a[1+2*3] = tetra[2+2*3] - tetra[2+0*3];
  a[1+3*3] = pow ( tetra[0+2*3] - tetra[0+0*3], 2 ) 
           + pow ( tetra[1+2*3] - tetra[1+0*3], 2 ) 
           + pow ( tetra[2+2*3] - tetra[2+0*3], 2 );

  a[2+0*3] = tetra[0+3*3] - tetra[0+0*3];
  a[2+1*3] = tetra[1+3*3] - tetra[1+0*3];
  a[2+2*3] = tetra[2+3*3] - tetra[2+0*3];
  a[2+3*3] = pow ( tetra[0+3*3] - tetra[0+0*3], 2 ) 
           + pow ( tetra[1+3*3] - tetra[1+0*3], 2 ) 
           + pow ( tetra[2+3*3] - tetra[2+0*3], 2 );
//
//  Solve the linear system.
//
  info = r8mat_solve ( DIM_NUM, RHS_NUM, a );
//
//  If the system was singular, return a consolation prize.
//
  if ( info != 0 )
  {
    *r = -1.0;
    r8vec_zero ( DIM_NUM, pc );
    return;
  }
//
//  Compute the radius and center.
//
  *r = 0.5 * sqrt 
    ( a[0+3*3] * a[0+3*3] 
    + a[1+3*3] * a[1+3*3] 
    + a[2+3*3] * a[2+3*3] );

  pc[0] = tetra[0+0*3] + 0.5 * a[0+3*3];
  pc[1] = tetra[1+0*3] + 0.5 * a[1+3*3];
  pc[2] = tetra[2+0*3] + 0.5 * a[2+3*3];

  return;
# undef DIM_NUM
# undef RHS_NUM
}
//****************************************************************************80

double *tetrahedron_dihedral_angles_3d ( double tetra[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_DIHEDRAL_ANGLES_3D computes dihedral angles of a tetrahedron.
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
//    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron,
//    which can be labeled as A, B, C and D.
//
//    Output, double TETRAHEDRON_DIHEDRAL_ANGLES_3D[6], the dihedral angles 
//    along the axes AB, AC, AD, BC, BD and CD, respectively.
//
{
  double ab[3];
  double *abc_normal;
  double *abd_normal;
  double ac[3];
  double *acd_normal;
  double ad[3];
  double *angle;
  double bc[3];
  double *bcd_normal;
  double bd[3];
  int i;
  double pi = 3.141592653589793;

  for ( i = 0; i < 3; i++ )
  {
    ab[i] = tetra[i+1*3] - tetra[i+0*3];
    ac[i] = tetra[i+2*3] - tetra[i+0*3];
    ad[i] = tetra[i+3*3] - tetra[i+0*3];
    bc[i] = tetra[i+2*3] - tetra[i+1*3];
    bd[i] = tetra[i+3*3] - tetra[i+1*3];
  }

  abc_normal = r8vec_cross_3d ( ac, ab );
  abd_normal = r8vec_cross_3d ( ab, ad );
  acd_normal = r8vec_cross_3d ( ad, ac );
  bcd_normal = r8vec_cross_3d ( bc, bd );

  angle = new double[6];

  angle[0] = r8vec_angle_3d ( abc_normal, abd_normal );
  angle[1] = r8vec_angle_3d ( abc_normal, acd_normal );
  angle[2] = r8vec_angle_3d ( abd_normal, acd_normal );
  angle[3] = r8vec_angle_3d ( abc_normal, bcd_normal );
  angle[4] = r8vec_angle_3d ( abd_normal, bcd_normal );
  angle[5] = r8vec_angle_3d ( acd_normal, bcd_normal );

  for ( i = 0; i < 6; i++ )
  {
    angle[i] = pi - angle[i];
  }

  delete [] abc_normal;
  delete [] abd_normal;
  delete [] acd_normal;
  delete [] bcd_normal;

  return angle;
}

//****************************************************************************80

double *tetrahedron_edge_length_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_EDGE_LENGTH_3D returns edge lengths of a tetrahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the tetrahedron vertices.
//
//    Output, double EDGE_LENGTH[6], the length of the edges.
//
{
# define DIM_NUM 3

  double *edge_length;
  int i;
  int j1;
  int j2;
  int k;
  double v[DIM_NUM];

  edge_length = new double[6];

  k = 0;
  for ( j1 = 0; j1 < 3; j1++ )
  {
    for ( j2 = j1 + 1; j2 < 4; j2++ )
    {
      for ( i = 0; i < DIM_NUM; i++ )
      {
        v[i] = tetra[i+j2*DIM_NUM] - tetra[i+j1*DIM_NUM];
      }
      edge_length[k] = r8vec_length ( DIM_NUM, v );
      k = k + 1;
    }
  }

  return edge_length;
# undef DIM_NUM
}
//****************************************************************************80

void tetrahedron_face_angles_3d ( double tetra[], double angles[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_FACE_ANGLES_3D returns the 12 face angles of a tetrahedron 3D.
//
//  Discussion:
//
//    The tetrahedron has 4 triangular faces.  This routine computes the
//    3 planar angles associated with each face.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4] the tetrahedron vertices.
//
//    Output, double ANGLES[3*4], the face angles.
//
{
  double *tri;

  tri = new double[3*3];
//
//  Face 123
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+1*3];
  tri[1+1*3] = tetra[1+1*3];
  tri[2+1*3] = tetra[2+1*3];
  tri[0+2*3] = tetra[0+2*3];
  tri[1+2*3] = tetra[1+2*3];
  tri[2+2*3] = tetra[2+2*3];

  triangle_angles_3d ( tri, angles );
//
//  Face 124
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+1*3];
  tri[1+1*3] = tetra[1+1*3];
  tri[2+1*3] = tetra[2+1*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  triangle_angles_3d ( tri, angles+3 );
//
//  Face 134
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+2*3];
  tri[1+1*3] = tetra[1+2*3];
  tri[2+1*3] = tetra[2+2*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  triangle_angles_3d ( tri, angles+6 );
//
//  Face 234
//
  tri[0+0*3] = tetra[0+1*3];
  tri[1+0*3] = tetra[1+1*3];
  tri[2+0*3] = tetra[2+1*3];
  tri[0+1*3] = tetra[0+2*3];
  tri[1+1*3] = tetra[1+2*3];
  tri[2+1*3] = tetra[2+2*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  triangle_angles_3d ( tri, angles+9 );

  delete [] tri;

  return;
}
//****************************************************************************80

void tetrahedron_face_areas_3d ( double tetra[], double areas[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_FACE_AREAS_3D returns the 4 face areas of a tetrahedron 3D.
//
//  Discussion:
//
//    The tetrahedron has 4 triangular faces.  This routine computes the
//    areas associated with each face.
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
//    Input, double TETRA[3*4] the tetrahedron vertices.
//
//    Output, double AREAS[4], the face areas.
//
{
  double *tri;

  tri = new double[3*3];
//
//  Face 123
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+1*3];
  tri[1+1*3] = tetra[1+1*3];
  tri[2+1*3] = tetra[2+1*3];
  tri[0+2*3] = tetra[0+2*3];
  tri[1+2*3] = tetra[1+2*3];
  tri[2+2*3] = tetra[2+2*3];

  areas[0] = triangle_area_3d ( tri );
//
//  Face 124
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+1*3];
  tri[1+1*3] = tetra[1+1*3];
  tri[2+1*3] = tetra[2+1*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  areas[1] = triangle_area_3d ( tri );
//
//  Face 134
//
  tri[0+0*3] = tetra[0+0*3];
  tri[1+0*3] = tetra[1+0*3];
  tri[2+0*3] = tetra[2+0*3];
  tri[0+1*3] = tetra[0+2*3];
  tri[1+1*3] = tetra[1+2*3];
  tri[2+1*3] = tetra[2+2*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  areas[2] = triangle_area_3d ( tri );
//
//  Face 234
//
  tri[0+0*3] = tetra[0+1*3];
  tri[1+0*3] = tetra[1+1*3];
  tri[2+0*3] = tetra[2+1*3];
  tri[0+1*3] = tetra[0+2*3];
  tri[1+1*3] = tetra[1+2*3];
  tri[2+1*3] = tetra[2+2*3];
  tri[0+2*3] = tetra[0+3*3];
  tri[1+2*3] = tetra[1+3*3];
  tri[2+2*3] = tetra[2+3*3];

  areas[3] = triangle_area_3d ( tri );

  delete [] tri;

  return;
}
//****************************************************************************80

void tetrahedron_insphere_3d ( double tetra[3*4], double *r, double pc[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_INSPHERE_3D finds the insphere of a tetrahedron in 3D.
//
//  Discussion:
//
//    The insphere of a tetrahedron is the inscribed sphere, which touches
//    each face of the tetrahedron at a single point.
//
//    The points of contact are the centroids of the triangular faces
//    of the tetrahedron.  Therefore, the point of contact for a face
//    can be computed as the average of the vertices of that face.
//
//    The sphere can then be determined as the unique sphere through
//    the four given centroids.
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
//  Reference:
//
//    Philip Schneider, David Eberly,
//    Geometric Tools for Computer Graphics,
//    Elsevier, 2002,
//    ISBN: 1558605940,
//    LC: T385.G6974.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double *R, PC[3], the radius and the center
//    of the sphere.
//
{
# define DIM_NUM 3

  double b[4*4];
  double gamma;
  int i;
  int j;
  double l123;
  double l124;
  double l134;
  double l234;
  double *n123;
  double *n124;
  double *n134;
  double *n234;
  double v21[DIM_NUM];
  double v31[DIM_NUM];
  double v41[DIM_NUM];
  double v32[DIM_NUM];
  double v42[DIM_NUM];
  double v43[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    v21[i] = tetra[i+1*DIM_NUM] - tetra[i+0*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v31[i] = tetra[i+2*DIM_NUM] - tetra[i+0*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v41[i] = tetra[i+3*DIM_NUM] - tetra[i+0*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v32[i] = tetra[i+2*DIM_NUM] - tetra[i+1*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v42[i] = tetra[i+3*DIM_NUM] - tetra[i+1*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v43[i] = tetra[i+3*DIM_NUM] - tetra[i+2*DIM_NUM];
  }

  n123 = r8vec_cross_3d ( v21, v31 );
  n124 = r8vec_cross_3d ( v41, v21 );
  n134 = r8vec_cross_3d ( v31, v41 );
  n234 = r8vec_cross_3d ( v42, v32 );

  l123 = r8vec_length ( DIM_NUM, n123 );
  l124 = r8vec_length ( DIM_NUM, n124 );
  l134 = r8vec_length ( DIM_NUM, n134 );
  l234 = r8vec_length ( DIM_NUM, n234 );

  delete [] n123;
  delete [] n124;
  delete [] n134;
  delete [] n234;

  for ( i = 0; i < DIM_NUM; i++ )
  {
    pc[i] = ( l234 * tetra[i+0*DIM_NUM]
            + l134 * tetra[i+1*DIM_NUM]
            + l124 * tetra[i+2*DIM_NUM]
            + l123 * tetra[i+3*DIM_NUM] )
            / ( l234 + l134 + l124 + l123 );
  }

  for ( j = 0; j < 4; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      b[i+j*4] = tetra[i+j*DIM_NUM];
    }
    b[3+j*4] = 1.0;
  }
  
  gamma = r8_abs ( r8mat_det_4d ( b ) );

  *r = gamma / ( l234 + l134 + l124 + l123 );

  return;
# undef DIM_NUM
}
//****************************************************************************80

double tetrahedron_quality1_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_QUALITY1_3D: "quality" of a tetrahedron in 3D.
//
//  Discussion:
//
//    The quality of a tetrahedron is 3.0 times the ratio of the radius of
//    the inscribed sphere divided by that of the circumscribed sphere.
//
//    An equilateral tetrahredron achieves the maximum possible quality of 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the tetrahedron vertices.
//
//    Output, double TETRAHEDRON_QUALITY1_3D, the quality of the tetrahedron.
//
{
# define DIM_NUM 3

  double pc[DIM_NUM];
  double quality;
  double r_in;
  double r_out;

  tetrahedron_circumsphere_3d ( tetra, &r_out, pc );

  tetrahedron_insphere_3d ( tetra, &r_in, pc );

  quality = 3.0 * r_in / r_out;

  return quality;
# undef DIM_NUM
}
//****************************************************************************80

double tetrahedron_quality2_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_QUALITY2_3D: "quality" of a tetrahedron in 3D.
//
//  Discussion:
//
//    The quality measure #2 of a tetrahedron is:
//
//      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
//
//    where
//
//      RIN = radius of the inscribed sphere;
//      LMAX = length of longest side of the tetrahedron.
//
//    An equilateral tetrahredron achieves the maximum possible quality of 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Desheng Wang,
//    The Optimal Centroidal Voronoi Tesselations and the Gersho's
//    Conjecture in the Three-Dimensional Space,
//    Computers and Mathematics with Applications,
//    Volume 49, 2005, pages 1355-1373.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the tetrahedron vertices.
//
//    Output, double TETRAHEDRON_QUALITY2_3D, the quality of the tetrahedron.
//
{
# define DIM_NUM 3

  double *edge_length;
  double l_max;
  double pc[DIM_NUM];
  double quality2;
  double r_in;

  edge_length = tetrahedron_edge_length_3d ( tetra );

  l_max = r8vec_max ( 6, edge_length );

  tetrahedron_insphere_3d ( tetra, &r_in, pc );

  quality2 = 2.0 * sqrt ( 6.0 ) * r_in / l_max;

  delete [] edge_length;

  return quality2;
# undef DIM_NUM
}
//****************************************************************************80

double tetrahedron_quality3_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_QUALITY3_3D computes the mean ratio of a tetrahedron.
//
//  Discussion:
//
//    This routine computes QUALITY3, the eigenvalue or mean ratio of
//    a tetrahedron.
//
//      QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
//
//    This value may be used as a shape quality measure for the tetrahedron.
//
//    For an equilateral tetrahedron, the value of this quality measure
//    will be 1.  For any other tetrahedron, the value will be between
//    0 and 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2005
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double TETRA(3,4), the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_QUALITY3_3D, the mean ratio of the tetrahedron.
//
{
# define DIM_NUM 3

  double ab[DIM_NUM];
  double ac[DIM_NUM];
  double ad[DIM_NUM];
  double bc[DIM_NUM];
  double bd[DIM_NUM];
  double cd[DIM_NUM];
  double denom;
  int i;
  double lab;
  double lac;
  double lad;
  double lbc;
  double lbd;
  double lcd;
  double quality3;
  double volume;
//
//  Compute the vectors representing the sides of the tetrahedron.
//
  for ( i = 0; i < DIM_NUM; i++ )
  {
    ab[i] = tetra[i+1*DIM_NUM] - tetra[i+0*DIM_NUM];
    ac[i] = tetra[i+2*DIM_NUM] - tetra[i+0*DIM_NUM];
    ad[i] = tetra[i+3*DIM_NUM] - tetra[i+0*DIM_NUM];
    bc[i] = tetra[i+2*DIM_NUM] - tetra[i+1*DIM_NUM];
    bd[i] = tetra[i+3*DIM_NUM] - tetra[i+1*DIM_NUM];
    cd[i] = tetra[i+3*DIM_NUM] - tetra[i+2*DIM_NUM];
  }
//
//  Compute the squares of the lengths of the sides.
//
  lab = pow ( ab[0], 2 ) + pow ( ab[1], 2 ) + pow ( ab[2], 2 );
  lac = pow ( ac[0], 2 ) + pow ( ac[1], 2 ) + pow ( ac[2], 2 );
  lad = pow ( ad[0], 2 ) + pow ( ad[1], 2 ) + pow ( ad[2], 2 );
  lbc = pow ( bc[0], 2 ) + pow ( bc[1], 2 ) + pow ( bc[2], 2 );
  lbd = pow ( bd[0], 2 ) + pow ( bd[1], 2 ) + pow ( bd[2], 2 );
  lcd = pow ( cd[0], 2 ) + pow ( cd[1], 2 ) + pow ( cd[2], 2 );
//
//  Compute the volume.
//
  volume = r8_abs ( 
      ab[0] * ( ac[1] * ad[2] - ac[2] * ad[1] ) 
    + ab[1] * ( ac[2] * ad[0] - ac[0] * ad[2] ) 
    + ab[2] * ( ac[0] * ad[1] - ac[1] * ad[0] ) ) / 6.0;

  denom = lab + lac + lad + lbc + lbd + lcd;

  if ( denom == 0.0 )
  {
    quality3 = 0.0;
  }
  else
  {
    quality3 = 12.0 * pow ( 3.0 * volume, 2.0 / 3.0 ) / denom;
  }

  return quality3;
# undef DIM_NUM
}
//****************************************************************************80

double tetrahedron_quality4_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_QUALITY4_3D computes the minimum solid angle of a tetrahedron.
//
//  Discussion:
//
//    This routine computes a quality measure for a tetrahedron, based
//    on the sine of half the minimum of the four solid angles.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2005
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double QUALITY4, the value of the quality measure.
//
{
# define DIM_NUM 3

  double ab[DIM_NUM];
  double ac[DIM_NUM];
  double ad[DIM_NUM];
  double bc[DIM_NUM];
  double bd[DIM_NUM];
  double cd[DIM_NUM];
  double denom;
  int i;
  double l1;
  double l2;
  double l3;
  double lab;
  double lac;
  double lad;
  double lbc;
  double lbd;
  double lcd;
  double quality4;
  double volume;
//
//  Compute the vectors that represent the sides.
//
  for ( i = 0; i < DIM_NUM; i++ )
  {
    ab[i] = tetra[i+1*DIM_NUM] - tetra[i+0*DIM_NUM];
    ac[i] = tetra[i+2*DIM_NUM] - tetra[i+0*DIM_NUM];
    ad[i] = tetra[i+3*DIM_NUM] - tetra[i+0*DIM_NUM];
    bc[i] = tetra[i+2*DIM_NUM] - tetra[i+1*DIM_NUM];
    bd[i] = tetra[i+3*DIM_NUM] - tetra[i+1*DIM_NUM];
    cd[i] = tetra[i+3*DIM_NUM] - tetra[i+2*DIM_NUM];
  }
//
//  Compute the lengths of the sides.
//
  lab = r8vec_length ( DIM_NUM, ab );
  lac = r8vec_length ( DIM_NUM, ac );
  lad = r8vec_length ( DIM_NUM, ad );
  lbc = r8vec_length ( DIM_NUM, bc );
  lbd = r8vec_length ( DIM_NUM, bd );
  lcd = r8vec_length ( DIM_NUM, cd );
//
//  Compute the volume.
//
  volume = r8_abs ( 
      ab[0] * ( ac[1] * ad[2] - ac[2] * ad[1] ) 
    + ab[1] * ( ac[2] * ad[0] - ac[0] * ad[2] ) 
    + ab[2] * ( ac[0] * ad[1] - ac[1] * ad[0] ) ) / 6.0;

  quality4 = 1.0;

  l1 = lab + lac;
  l2 = lab + lad;
  l3 = lac + lad;

  denom = ( l1 + lbc ) * ( l1 - lbc ) 
        * ( l2 + lbd ) * ( l2 - lbd ) 
        * ( l3 + lcd ) * ( l3 - lcd );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = r8_min ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  l1 = lab + lbc;
  l2 = lab + lbd;
  l3 = lbc + lbd;

  denom = ( l1 + lac ) * ( l1 - lac ) 
        * ( l2 + lad ) * ( l2 - lad ) 
        * ( l3 + lcd ) * ( l3 - lcd );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = r8_min ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  l1 = lac + lbc;
  l2 = lac + lcd;
  l3 = lbc + lcd;

  denom = ( l1 + lab ) * ( l1 - lab ) 
        * ( l2 + lad ) * ( l2 - lad ) 
        * ( l3 + lbd ) * ( l3 - lbd );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = r8_min ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  l1 = lad + lbd;
  l2 = lad + lcd;
  l3 = lbd + lcd;

  denom = ( l1 + lab ) * ( l1 - lab ) 
        * ( l2 + lac ) * ( l2 - lac ) 
        * ( l3 + lbc ) * ( l3 - lbc );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = r8_min ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  quality4 = quality4 * 1.5 * sqrt ( 6.0 );

  return quality4;
# undef DIM_NUM
}
//****************************************************************************80

double *tetrahedron_solid_angles_3d ( double tetra[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_SOLID_ANGLES_3D computes solid angles of a tetrahedron.
//
//  Modified:
//
//    07 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_SOLID_ANGLES_3D[4], the solid angles.
//
{
  double *angle;
  double *dihedral_angles;
  double pi = 3.141592653589793;

  dihedral_angles = tetrahedron_dihedral_angles_3d ( tetra );

  angle = new double[4];

  angle[0] = dihedral_angles[0] + dihedral_angles[1] + dihedral_angles[2] - pi;
  angle[1] = dihedral_angles[0] + dihedral_angles[3] + dihedral_angles[4] - pi;
  angle[2] = dihedral_angles[1] + dihedral_angles[3] + dihedral_angles[5] - pi;
  angle[3] = dihedral_angles[2] + dihedral_angles[4] + dihedral_angles[5] - pi;

  delete [] dihedral_angles;

  return angle;
}
//****************************************************************************80

double tetrahedron_volume_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_VOLUME_3D computes the volume of a tetrahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_VOLUME_3D, the volume of the tetrahedron.
//
{
  double a[4*4];
  int i;
  int j;
  double volume;

  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < 4; j++ )
    { 
      a[i+j*4] = tetra[i+j*3];
    }
  }

  i = 3;
  for ( j = 0; j < 4; j++ )
  {
    a[i+j*4] = 1.0;
  }

  volume = r8_abs ( r8mat_det_4d ( a ) ) / 6.0;

  return volume;
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

void triangle_angles_3d ( double t[3*3], double angle[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
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
//    Input, double T[3*3], the triangle vertices.
//
//    Output, double ANGLE[3], the angles opposite
//    sides P1-P2, P2-P3 and P3-P1, in radians.
//
{
# define DIM_NUM 3

  double a;
  double b;
  double c;
  double pi = 3.141592653589793;

  a = sqrt ( pow ( t[0+1*DIM_NUM] - t[0+0*DIM_NUM], 2 ) 
           + pow ( t[1+1*DIM_NUM] - t[1+0*DIM_NUM], 2 )
           + pow ( t[2+1*DIM_NUM] - t[2+0*DIM_NUM], 2 ) );

  b = sqrt ( pow ( t[0+2*DIM_NUM] - t[0+1*DIM_NUM], 2 ) 
           + pow ( t[1+2*DIM_NUM] - t[1+1*DIM_NUM], 2 )
           + pow ( t[2+2*DIM_NUM] - t[2+1*DIM_NUM], 2 ) );

  c = sqrt ( pow ( t[0+0*DIM_NUM] - t[0+2*DIM_NUM], 2 ) 
           + pow ( t[1+0*DIM_NUM] - t[1+2*DIM_NUM], 2 )
           + pow ( t[2+0*DIM_NUM] - t[2+2*DIM_NUM], 2 ) );
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
# undef DIM_NUM
}
//****************************************************************************80

double triangle_area_3d ( double t[3*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA_3D computes the area of a triangle in 3D.
//
//  Discussion:
//
//    This routine uses the fact that the norm of the cross product vector
//    is the area of the parallelogram they form.  The triangle they
//    form has half that area.
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
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[3*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA_3D, the area of the triangle.
//
{
# define DIM_NUM 3

  double area;
  double *cross;
  int i;
//
//  Compute the cross product vector.
//
  cross = new double[DIM_NUM];

  cross[0] = ( t[1+1*DIM_NUM] - t[1+0*DIM_NUM] ) 
           * ( t[2+2*DIM_NUM] - t[2+0*DIM_NUM] ) 
           - ( t[2+1*DIM_NUM] - t[2+0*DIM_NUM] ) 
           * ( t[1+2*DIM_NUM] - t[1+0*DIM_NUM] );

  cross[1] = ( t[2+1*DIM_NUM] - t[2+0*DIM_NUM] ) 
           * ( t[0+2*DIM_NUM] - t[0+0*DIM_NUM] ) 
           - ( t[0+1*DIM_NUM] - t[0+0*DIM_NUM] ) 
           * ( t[2+2*DIM_NUM] - t[2+0*DIM_NUM] );

  cross[2] = ( t[0+1*DIM_NUM] - t[0+0*DIM_NUM] ) 
           * ( t[1+2*DIM_NUM] - t[1+0*DIM_NUM] ) 
           - ( t[1+1*DIM_NUM] - t[1+0*DIM_NUM] ) 
           * ( t[0+2*DIM_NUM] - t[0+0*DIM_NUM] );

  area = 0.0;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    area = area + pow ( cross[i], 2 );
  }
  
  area = 0.5 * sqrt ( area );

  delete [] cross;

  return area;
# undef DIM_NUM
}
