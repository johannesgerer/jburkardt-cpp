# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "cities.hpp"

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
//    13 September 2009
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
  bool value;

  if ( 97 <= ch1 && ch1 <= 122 ) 
  {
    ch1 = ch1 - 32;
  } 
  if ( 97 <= ch2 && ch2 <= 122 ) 
  {
    ch2 = ch2 - 32;
  }     

  value = ( ch1 == ch2 );

  return value;
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

int dist_table_check ( int n, double dist_table[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIST_TABLE_CHECK checks a distance table.
//
//  Discussion:
//
//    1) All entries must be nonnegative.
//    2) Diagonal entries must be zero.
//    3) Off-diagonal entries must be symmetric.
//    4) The triangle inequality must be observed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of cities.
//
//    Input, double DIST_TABLE[N*N], the distance table.
//
//    Output, int DIST_TABLE_CHECK, the result of the check.
//    0, the matrix passed the checks.
//    1, Not all entries are nonnegative.
//    2, Not all diagonal entries are zero.
//    3, Not all off-diagonal entries are symmetric.
//    4, Not all entries satisfy the triangle inequality.
//
{
  int check;
  int i;
  int j;
  int k;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( dist_table[i+j*n] < 0.0 )
      {
        check = 1;
        return check;
      }
    }
  }

  for ( i = 0; i < n; i++ )
  {
    if ( dist_table[i+i*n] != 0.0 )
    {
      check = 2;
      return check;
    }
  }

  for ( j = 0; j < n; j ++ )
  {
    for ( i = 0; i < j; i++ )
    {
      if ( dist_table[i+j*n] != dist_table[j+i*n] )
      {
        check = 3;
        return check;
      }
    }
  }

  for ( k = 0; k < n; k++ );
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        if ( dist_table[i+j*n] + dist_table[j+k*n] < dist_table[i+k*n] )
        {
          check = 4;
          return check;
        }
      }
    }
  }

  check = 0;

  return check;
}
//****************************************************************************80

void dms_print ( int n, int lat_dms[], int long_dms[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    DMS_PRINT prints the latitude and longitude in degrees/minutes/seconds.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data items.
//
//    Input, int LAT_DMS[4*N], LONG_DMS[4*N], the latitudes
//    and longitudes, in degrees, minutes and seconds.  The fourth
//    argument is +1/-1 for North/South latitude or East/West longitude.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  cout << "    #    Latitude            Longitude\n";
  cout << "         (Deg/Min/Sec)       (Deg/Min/Sec)\n";
  cout << "  ---  ---------------      ---------------\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(3) << i
         << "  " << setw(4) << lat_dms[0+i*4]
         << "  " << setw(4) << lat_dms[1+i*4]
         << "  " << setw(4) << lat_dms[2+i*4]
         << "  "            << lat_char ( lat_dms[3+i*4] )
         << "   "
         << "  " << setw(4) << long_dms[0+i*4]
         << "  " << setw(4) << long_dms[1+i*4]
         << "  " << setw(4) << long_dms[2+i*4]
         << "  "            << long_char ( long_dms[3+i*4] ) << "\n";
  }
  return;
}
//****************************************************************************80

void dms_read ( string file_name, int n, int lat_dms[], int long_dms[] )

//****************************************************************************80
//
//  Purpose:
//
//    DMS_READ reads DMS data from a file.
//
//  Discussion:
//
//    Blank lines and comment lines (beginning with '#') are ignored.
//
//    Individual data values should be separated by spaces or commas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to read.
//
//    Input, int N, the number of data items.
//
//    Output, int LAT_DMS[4*N], LONG_DMS[4*N], the latitude and
//    longitudes, in degrees, minutes and seconds.  The fourth
//    argument is +1/-1 for North/South latitude or East/West longitude. 
//
{
  char c1;
  char c2;
  int d1;
  int d2;
  bool done;
  int i;
  ifstream input;
  int j;
  string line;
  int line_num;
  int m1;
  int m2;
  int s1;
  int s2;
  string w;

  input.open ( file_name.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "DMS_READ - Fatal error!\n";
    cerr << "  Could not open the input file.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      lat_dms[i+j*4] = i4_huge ( );
      long_dms[i+j*4] = i4_huge ( );
    }
  }
  j = 0;
  line_num = 0;

  for ( ; ; )
  {
//
//  Have we read enough data?
//
    if ( n <= j )
    {
      break;
    }
//
//  Read the next line from the file.
//
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    line_num = line_num + 1;
//
//  Skip blank lines and comment lines.
//
    if ( s_len_trim ( line ) == 0 )
    {
    }
    else if ( line[0] == '#' )
    {
    }
    else
    {
      sscanf ( line.c_str ( ), "%d  %d  %d  %c  %d  %d  %d  %c\n", 
        &d1, &m1, &s1, &c1, &d2, &m2, &s2, &c2 );

      lat_dms[0+j*4] = d1;
      lat_dms[1+j*4] = m1;
      lat_dms[2+j*4] = s1;
      if ( c1 == 'N' || c1 == 'n' )
      {
        lat_dms[3+j*4] = + 1;
      }
      else if ( c1 == 'S' || c1 == 's' )
      {
        lat_dms[3+j*4] = - 1;
      }
      else
      {
        lat_dms[3+j*4] =   0;
      }
      long_dms[0+j*4] = d2;
      long_dms[1+j*4] = m2;
      long_dms[2+j*4] = s2;
      if ( c2 == 'E' || c2 == 'E' )
      {
        long_dms[3+j*4] = + 1;
      }
      else if ( c2 == 'W' || c2 == 'W' )
      {
        long_dms[3+j*4] = - 1;
      }
      else
      {
        long_dms[3+j*4] =   0;
      }

      j = j + 1;
    }
  }

  input.close ( );

  cout << "\n";
  cout << "DMS_READ:\n";
  cout << "  Number of lines read was " << line_num << "\n";

  return;
}
//****************************************************************************80

double *dms_to_dist ( int n, int lat_dms[], int long_dms[] )

//****************************************************************************80
//
//  Purpose:
//
//    DMS_TO_DIST creates a distance table from latitudes and longitudes.
//
//  Discussion:
//
//    A distance function is used which is appropriate for the earth.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data items.
//
//    Input, int LAT_DMS[4*N], LONG_DMS[4*N], the latitude 
//    and longitude, in degrees, minutes, and seconds, for each point.
//    The fourth argument is +1/-1 for North/South latitude or 
//    East/West longitude.
//
//    Output, double DMS_TO_DIST_TABLE[N*N], the distance matrix.  Distances 
//    are measured in miles.
//
{
  double *dist_table;
  int i;
  int j;
  double value;

  dist_table = new double[n*n];

  for ( i = 0; i < n; i++ )
  {
    dist_table[i+i*n] = 0.0;
    for ( j = i + 1; j < n; j++ )
    {
      value = dms_to_distance_earth ( lat_dms+i*4, long_dms+i*4, 
        lat_dms+j*4, long_dms+j*4 );
      dist_table[i+j*n] = value;
      dist_table[j+i*n] = value;
    }
  }

  return dist_table;
}
//****************************************************************************80

double dms_to_distance_earth ( int lat1_dms[], int long1_dms[],
  int lat2_dms[], int long2_dms[] )

//****************************************************************************80
//
//  Purpose:
//
//    DMS_TO_DISTANCE_EARTH finds the distance between two points on the earth.
//
//  Discussion:
//
//    The positions of the the points are given as longitude and
//    latitude, measured in degrees, minutes, and seconds.
//
//    The distance is measured on the surface of the earth, which
//    is approximated by a sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LAT1_DMS(4), LONG1_DMS(4), the latitude and 
//    longitude of the first point.  The fourth
//    argument is +1/-1 for North/South latitude or East/West longitude.
//
//    Input, int LAT2_DMS(4), LONG2_DMS(4), the latitude and 
//    longitude of the second point.  The fourth
//    argument is +1/-1 for North/South latitude or East/West longitude.
//
//    Output, double DMS_TO_DISTANCE_EARTH, the distance between the points, 
//    in miles.
//
{
  double dist;
  double lat1_rad;
  double lat2_rad;
  double long1_rad;
  double long2_rad;
  double radius = 3958.89;
  double theta;

  lat1_rad = dms_to_radians ( lat1_dms );
  long1_rad = dms_to_radians ( long1_dms );

  lat2_rad = dms_to_radians ( lat2_dms );
  long2_rad = dms_to_radians ( long2_dms );

  theta = acos ( sin ( lat1_rad ) * sin ( lat2_rad ) 
               + cos ( lat1_rad ) * cos ( lat2_rad ) 
               * cos ( long1_rad - long2_rad ) );

  dist = radius * theta;

  return dist;
}
//****************************************************************************80

double dms_to_radians ( int dms[] )

//****************************************************************************80
//
//  Purpose:
//
//    DMS_TO_RADIANS converts degrees, minutes, seconds to radians.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DMS[4], the measurement of an angle in
//    degrees, minutes and seconds.  The fourth
//    argument is +1/-1 for North/South latitude or East/West longitude.
//
//    Output, double DMS_TO_RADIANS, the measurement of the same
//    angle in radians.
//
{
  double d;
  double pi = 3.141592653589793;
  double r;
  int s;

  s = i4_sign ( dms[3] ) * ( dms[2] + 60 * dms[1] + 3600 * dms[0] );
  d = ( double ) ( s ) / 3600.0;
  r = pi * d / 180.0;

  return r;
}
//****************************************************************************80

double *dms_to_xy ( int n, int lat_dms[], int long_dms[] )

//****************************************************************************80
//
//  Purpose:
//
//    DMS_TO_XY: Latitude/Longitude in DMS to XY coordinates.
//
//  Discussion:
//
//    Essentially, the latitude and longitude information is treated
//    as though the Earth were a cylinder.  As long as the the
//    data is relatively close on the sphere (and far from either
//    pole) the distortion will not be too severe.  If the data
//    is closely clustered, and also near the equator, the
//    positions will be relatively accurate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data items.
//
//    Input, int LAT_DMS[4*N], LONG_DMS[4*N], the latitude and
//    longitude, in degrees, minutes, and seconds, for each point.
//    The fourth argument is +1/-1 for North/South latitude or 
//    East/West longitude.
//
//    Output, double DMS_TO_XY[2*N], the point coordinates, in miles.
//
{
  int i;
  double phi;
  double *point_xy;
  double radius = 3958.89;
  double theta;

  point_xy = new double[2*n];

  for ( i = 0; i < n; i++ )
  {
    theta = dms_to_radians ( long_dms+i*4 );
    phi = dms_to_radians ( lat_dms+i*4 );
    point_xy[0+i*2] = radius * theta;
    point_xy[1+i*2] = radius * phi;
  }

  return point_xy;
}
//****************************************************************************80

void dms_write ( string file_name, int n, int lat_dms[], int long_dms[] )

//****************************************************************************80
//
//  Purpose:
//
//    DMS_WRITE writes a DMS latitude, longitude file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to write.
//
//    Input, int N, the number of data items.
//
//    Input, int LAT_DMS[4*N], LONG_DMS[4*N], the data values.
//
{
  int i;
  ofstream output;

  output.open ( file_name.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "DMS_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }

  output << "# " << file_name << "\n";
  output << "#\n";
  output << "#  Created by DMS_WRITE.\n";
  output << "#\n";
  output << "#  Latitude, Longitude in degrees, minutes, seconds\n";
  output << "#  Number of points N is " << n << "\n";
  output << "#\n";

  for ( i = 0; i < n; i++ )
  {
    output << "  " << setw(3) << lat_dms[0+i*4]
           << "  " << setw(3) << lat_dms[1+i*4]
           << "  " << setw(3) << lat_dms[2+i*4]
           << "  " << lat_char ( lat_dms[3+i*4] )
           << "  " << setw(3) << long_dms[0+i*4]
           << "  " << setw(3) << long_dms[1+i*4]
           << "  " << setw(3) << long_dms[2+i*4]
           << "  " << lat_char ( long_dms[3+i*4] ) << "\n";
  }


  output.close ( );

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

int i4_huge ( )

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
//    Output, int I4_HUGE, a "huge" I4.
//
{
  return 2147483647;
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

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
  int value;

  if ( i < 0 ) 
  {
    value = -1;
  }
  else
  {
    value = 1;
  }
  return value;
}
//****************************************************************************80

char lat_char ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    LAT_CHAR returns a character for negative or positive latitude.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, is negative for negative latitude, 
//    and positive for positive latitude.
//
//    Output, char LAT_CHAR, is 'S' for negative latitude, and
//    'N' for positive latitude.
//
{
  char value;

  if ( i < 0 )
  {
    value = 'S';
  }
  else if ( 0 < i )
  {
    value = 'N';
  }
  else
  {
    value = '?';
  }
  return value;
}
//****************************************************************************80

double ll_rad_dist_sphere ( double lat1,double  long1, double lat2, 
  double long2, double radius )

//****************************************************************************80
//
//  Purpose:
//
//    LL_RAD_DIST_SPHERE: spherical distance, latitude and longitude in radians.
//
//  Discussion:
//
//    On a sphere of given radius, the positions of two points are given as
//    longitude and latitude, in radians.
//
//    This function determines the spherical distance or great circle distance,
//    between the two points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double LAT1, LONG1, LAT2, LONG2, the latitude and
//    longitude of the two points, in radians.
//
//    Input, double RADIUS, the radius of the sphere.
//
//    Output, double LL_RAD_DIST_SPHERE, the distance between the points.
//
{
  double dist;
  double theta;

  theta = acos ( sin ( lat1 ) * sin ( lat2 ) 
               + cos ( lat1 ) * cos ( lat2 ) * cos ( long1 - long2 ) );

  dist = radius * theta;

  return dist;
}
//****************************************************************************80

char long_char ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    LONG_CHAR returns a character for negative or positive longitude.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, is negative for negative longitude, and 
//    positive for positive longitude.
//
//    Output, char LONG_CHAR, is 'W' for negative longitude, and
//    'E' for positive longitude.
//
{
  char value;

  if ( i < 0 )
  {
    value = 'W';
  }
  else if ( 0 < i )
  {
    value = 'E';
  }
  else
  {
    value = '?';
  }
  return value;
}
//****************************************************************************80

string main_read_code ( string file_main )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN_READ_CODE reads the name of the code file from the main file.
//
//  Discussion:
//
//    FILE_CODE is the name of a file containing short codes for the 
//    cities.
//
//    There MAY be a record in the main file of the form
//
//    "code  key_code.txt"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_MAIN, the name of the file to read.
//
//    Output, string MAIN_READ_CODE, the name of the code file,
//    or ' ' if no information was found.
//
{
  bool done;
  string file_code;
  ifstream input;
  string line;
  string word;

  input.open ( file_main.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MAIN_READ_CODE - Fatal error!\n";
    cerr << "  Could not open the input file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
//
//  Read the next line from the file.
//
    getline ( input, line );

    if ( input.eof ( ) )
    {
      file_code = " ";
      break;
    }
//
//  Skip blank lines and comment lines.
//
    if ( s_len_trim ( line ) == 0 )
    {
    }
    else if ( line[0] == '#' )
    {
    }
    else
    {
      done = true;

      word = word_next_read ( line, &done );

      if ( done )
      {
        continue;
      }

      if ( !s_eqi ( word, "code" ) )
      {
        continue;
      }

      file_code = word_next_read ( line, &done );

      break;
    }
  }

  input.close ( );

  return file_code;
}
//****************************************************************************80

string main_read_dist ( string file_main )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN_READ_DIST reads the name of the distance file from the main file.
//
//  Discussion:
//
//    FILE_DIST is the name of a file containing the city-to-city
//    distance matrix.
//
//    There MAY be a record in the main file of the form
//
//    "dist  key_dist.txt"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_MAIN, the name of the file to read.
//
//    Output, string MAIN_READ_DIST, the name of the distance file,
//    or " " if no information was found.
//
{
  bool done;
  string file_dist;
  ifstream input;
  string line;
  string word;

  input.open ( file_main.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MAIN_READ_DIST - Fatal error!\n";
    cerr << "  Could not open the input file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
//
//  Read the next line from the file.
//
    getline ( input, line );

    if ( input.eof ( ) )
    {
      file_dist = " ";
      break;
    }
//
//  Skip blank lines and comment lines.
//
    if ( s_len_trim ( line ) == 0 )
    {
    }
    else if ( line[0] == '#' )
    {
    }
    else
    {
      done = true;

      word = word_next_read ( line, &done );

      if ( done )
      {
        continue;
      }

      if ( !s_eqi ( word, "dist" ) )
      {
        continue;
      }

      file_dist = word_next_read ( line, &done );

      break;
    }
  }

  input.close ( );

  return file_dist;
}
//****************************************************************************80

string main_read_dms ( string file_main )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN_READ_DMS reads the name of the DMS file from the main file.
//
//  Discussion:
//
//    FILE_DMS is the name of a file containing the longitude and latitude
//    of each city in degrees/minutes/seconds.
//
//    There MAY be a record in the main file of the form
//
//    "dms  key_dms.txt"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_MAIN, the name of the file to read.
//
//    Output, string MAIN_READ_DMS, the name of the DMS file,
//    or ' ' if no information was found.
//
{
  bool done;
  string file_dms;
  ifstream input;
  string line;
  string word;

  input.open ( file_main.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MAIN_READ_DMS - Fatal error!\n";
    cerr << "  Could not open the input file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
//
//  Read the next line from the file.
//
    getline ( input, line );

    if ( input.eof ( ) )
    {
      file_dms = " ";
      break;
    }
//
//  Skip blank lines and comment lines.
//
    if ( s_len_trim ( line ) == 0 )
    {
    }
    else if ( line[0] == '#' )
    {
    }
    else
    {
      done = true;

      word = word_next_read ( line, &done );

      if ( done )
      {
        continue;
      }

      if ( !s_eqi ( word, "dms" ) )
      {
        continue;
      }

      file_dms = word_next_read ( line, &done );

      break;
    }
  }

  input.close ( );

  return file_dms;
}
//****************************************************************************80

string main_read_geom ( string file_main )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN_READ_GEOM reads the name of the geometry from the main file.
//
//  Discussion:
//
//    GEOM is the name of the geometry of the city data.
//    Typical values include:
//    none - no special geometry
//    plane - the points lie in a plane
//    sphere - the points lie on a sphere
//
//    There MAY be a record in the main file of the form
//
//    "geom  geom_value"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_MAIN, the name of the file to read.
//
//    Output, string GEOM, the name of the geometry,
//    or ' ' if no information was found.
//
{
  bool done;
  string geom;
  ifstream input;
  string line;
  string word;

  input.open ( file_main.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MAIN_READ_GEOM - Fatal error!\n";
    cerr << "  Could not open the input file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
//
//  Read the next line from the file.
//
    getline ( input, line );

    if ( input.eof ( ) )
    {
      geom = "";
      break;
    }
//
//  Skip blank lines and comment lines.
//
    if ( s_len_trim ( line ) == 0 )
    {
    }
    else if ( line[0] == '#' )
    {
    }
    else
    {
      done = true;

      word = word_next_read ( line, &done );

      if ( done )
      {
        continue;
      }

      if ( !s_eqi ( word, "geom" ) )
      {
        continue;
      }

      geom = word_next_read ( line, &done );

      break;
    }
  }

  input.close ( );

  return geom;
}
//****************************************************************************80

string main_read_name ( string file_main )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN_READ_NAME reads the name of the name file from the main file.
//
//  Discussion:
//
//    The name file, if it exists, contains a list of the city names.
//
//    There MAY be a record in the main file of the form
//
//    "name  key_name.txt"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_MAIN, the name of the file to read.
//
//    Output, string NAME, the name of the name file,
//    or ' ' if no information was found.
//
{
  bool done;
  ifstream input;
  string line;
  string name;
  string word;

  input.open ( file_main.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MAIN_READ_NAME - Fatal error!\n";
    cerr << "  Could not open the input file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
//
//  Read the next line from the file.
//
    getline ( input, line );

    if ( input.eof ( ) )
    {
      name = "";
      break;
    }
//
//  Skip blank lines and comment lines.
//
    if ( s_len_trim ( line ) == 0 )
    {
    }
    else if ( line[0] == '#' )
    {
    }
    else
    {
      done = true;

      word = word_next_read ( line, &done );

      if ( done )
      {
        continue;
      }

      if ( !s_eqi ( word, "name" ) )
      {
        continue;
      }

      name = word_next_read ( line, &done );

      break;
    }
  }

  input.close ( );

  return name;
}
//****************************************************************************80

int main_read_size ( string file_main )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN_READ_SIZE reads the problem size N from the main file.
//
//  Discussion:
//
//    The problem size is N, the number of cities.
//
//    There should always be a record in the main file of the form
//
//    "size  7"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_MAIN, the name of the file to read.
//
//    Output, int MAIN_READ_SIZE, the problem size.
//
{
  bool done;
  bool ierror;
  ifstream input;
  int length;
  string line;
  int n;
  string word;

  input.open ( file_main.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MAIN_READ_SIZE - Fatal error!\n";
    cerr << "  Could not open the input file.\n";
    exit ( 1 );
  }

  n = 0;

  for ( ; ; )
  {
//
//  Read the next line from the file.
//
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }
//
//  Skip blank lines and comment lines.
//
    if ( s_len_trim ( line ) == 0 )
    {
    }
    else if ( line[0] == '#' )
    {
    }
    else
    {
      done = true;

      word = word_next_read ( line, &done );

      if ( done )
      {
        continue;
      }

      if ( !s_eqi ( word, "size" ) )
      {
        continue;
      }

      word = word_next_read ( line, &done );

      n = s_to_i4 ( word, &length, &ierror );

      break;
    }
  }

  input.close ( );

  return n;
}
//****************************************************************************80

string main_read_xy ( string file_main )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN_READ_XY reads the name of the XY file from the main file.
//
//  Discussion:
//
//    The XY file, if it exists, contains a list of the XY coordinates
//    of cities.
//
//    There MAY be a record in the main file of the form
//
//    "xy  key_xy.txt"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_MAIN, the name of the file to read.
//
//    Output, string XY, the name of the XY file,
//    or ' ' if no information was found.
//
{
  bool done;
  ifstream input;
  string line;
  string word;
  string xy;

  input.open ( file_main.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MAIN_READ_XY - Fatal error!\n";
    cerr << "  Could not open the input file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
//
//  Read the next line from the file.
//
    getline ( input, line );

    if ( input.eof ( ) )
    {
      xy = "";
      break;
    }
//
//  Skip blank lines and comment lines.
//
    if ( s_len_trim ( line ) == 0 )
    {
    }
    else if ( line[0] == '#' )
    {
    }
    else
    {
      done = true;

      word = word_next_read ( line, &done );

      if ( done )
      {
        continue;
      }

      if ( !s_eqi ( word, "xy" ) )
      {
        continue;
      }

      xy = word_next_read ( line, &done );

      break;
    }
  }

  input.close ( );

  return xy;
}
//****************************************************************************80

double *point_to_dist_table ( int dim_num, int point_num, double point[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_TO_DIST_TABLE creates a distance table from Cartesian coordinates.
//
//  Discussion:
//
//    The euclidean distance is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the  point coordinates.
//
//    Output, double DIST_TABLE[POINT_NUM*POINT_NUM], the 
//    distance table.
//
{
  int dim;
  double *dist_table;
  int i;
  int j;

  dist_table = new double[point_num*point_num];

  for ( i = 0; i < point_num; i++ )
  {
    dist_table[i+i*point_num] = 0.0;
    for ( j = i + 1; j < point_num; j++ )
    {
      dist_table[i+j*point_num] = 0.0;
      for ( dim = 0; dim < dim_num; dim++ )
      { 
        dist_table[i+j*point_num] = dist_table[i+j*point_num] 
          + pow ( point[dim+i*dim_num] - point[dim+j*dim_num], 2 );
      }
      dist_table[i+j*point_num] = sqrt ( dist_table[i+j*point_num] );
      dist_table[j+i*point_num] = dist_table[i+j*point_num];
    }
  }
  return dist_table;
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
    value = + x;
  } 
  else
  {
    value = - x;
  }
  return value;
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
//    An R8MAT is an array of R8's.
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
//  Discussion:
//
//    An R8MAT is an array of R8's.
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

void r8mat_nint ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NINT rounds the entries of an R8MAT.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of A.
//
//    Input/output, double A[M*N], the matrix to be NINT'ed.
//
{
  int i;
  int j;
  int s;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*m] < 0.0 )
      {
        s = -1;
      }
      else
      {
        s = 1;
      }
      a[i+j*m] = s * ( int ) ( r8_abs ( a[i+j*m] ) + 0.5 );
    }
  }

  return;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2009
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
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 August 2010
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
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
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
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
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
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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
//    An R8MAT is an array of R8's.
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
//    An R8MAT is an array of R8's.
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
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
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
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
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
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

bool s_eqi ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
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
//    Input, string S1, S2, two strings.
//
//    Output, bool S_EQI, is true if the strings are equal. 
//
{
  int i;
  int nchar;
  int s1_length;
  int s2_length;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  if ( s1_length < s2_length )
  {
    nchar = s1_length;
  }
  else
  {
    nchar = s2_length;
  }
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ ) 
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) ) 
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < s1_length ) 
  {
    for ( i = nchar; i < s1_length; i++ ) 
    {
      if ( s1[i] != ' ' ) 
      {
        return false;
      }
    } 
  }
  else if ( nchar < s2_length ) 
  {
    for ( i = nchar; i < s2_length; i++ )
    {
      if ( s2[i] != ' ' ) 
      {
        return false;
      }
    } 
  }

  return true;
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

string word_next_read ( string s, bool *done )

//****************************************************************************80
//
//  Purpose:
//
//    WORD_NEXT_READ "reads" words from a string, one at a time.
//
//  Discussion:
//
//    This routine was written to process tokens in a file.
//    A token is considered to be an alphanumeric string delimited
//    by whitespace, or any of various "brackets".
//
//    The following characters are considered to be a single word,
//    whether surrounded by spaces or not:
//
//      " ( ) { } [ ]
//
//    Also, if there is a trailing comma on the word, it is stripped off.
//    This is to facilitate the reading of lists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string, presumably containing words
//    separated by spaces.
//
//    Input/output, bool *DONE.
//    On input with a fresh string, set DONE to TRUE.
//    On output, the routine sets DONE:
//      FALSE if another word was read,
//      TRUE if no more words could be read.
//
//    Output, string WORD_NEXT_READ.
//    If DONE is FALSE, then WORD contains the "next" word read.
//    If DONE is TRUE, then WORD is NULL, because there was no more to read.
//
{
  int i;
  int ilo;
  int j;
  static int lenc = 0;
  static int next = 0;
  char TAB = 9;
  string word;
  char *word_chstar;
//
//  We "remember" LENC and NEXT from the previous call.
//
//  An input value of DONE = TRUE signals a new line of text to examine.
//
  if ( *done )
  {
    next = 0;
    *done = false;
    lenc = s.length ( );
    if ( lenc <= 0 )
    {
      *done = true;
      word = "\n";;
      return word;
    }
  }
//
//  Beginning at index NEXT, search the string for the next nonblank,
//  which signals the beginning of a word.
//
  ilo = next;
//
//  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
//
  for ( ; ; )
  {
    if ( lenc < ilo )
    {
      word = "\n";
      *done = true;
      next = lenc + 1;
      return word;
    }
//
//  If the current character is blank, skip to the next one.
//
    if ( s[ilo] != ' ' && s[ilo] != TAB )
    {
      break;
    }
    ilo = ilo + 1;
  }
//
//  ILO is the index of the next nonblank character in the string.
//
//  If this initial nonblank is a special character,
//  then that's the whole word as far as we're concerned,
//  so return immediately.
//
  if ( s[ilo] == '"' )
  {
    word = """";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '(' )
  {
    word = "(";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == ')' )
  {
    word = ")";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '{' )
  {
    word = "{";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '}' )
  {
    word = "}";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '[' )
  {
    word = "[";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == ']' )
  {
    word = "]";
    next = ilo + 1;
    return word;
  }
//
//  Now search for the last contiguous character that is not a
//  blank, TAB, or special character.
//
  next = ilo + 1;

  while ( next <= lenc )
  {
    if ( s[next] == ' ' )
    {
      break;
    }
    else if ( s[next] == TAB )
    {
      break;
    }
    else if ( s[next] == '"' )
    {
      break;
    }
    else if ( s[next] == '(' )
    {
      break;
    }
    else if ( s[next] == ')' )
    {
      break;
    }
    else if ( s[next] == '{' )
    {
      break;
    }
    else if ( s[next] == '}' )
    {
      break;
    }
    else if ( s[next] == '[' )
    {
      break;
    }
    else if ( s[next] == ']' )
    {
      break;
    }

    next = next + 1;
  }
//
//  Allocate WORD, copy characters, and return.
//
  if ( s[next-1] == ',' )
  {
    word_chstar = new char[next-ilo];
    i = 0;
    for ( j = ilo; j <= next - 2; j++ )
    {
      word_chstar[i] = s[j];
      i = i + 1;
    }
    word_chstar[i] = '\0';
    word = string ( word_chstar );
    delete [] word_chstar;
  }
  else
  {
    word_chstar = new char[next+1-ilo];
    i = 0;
    for ( j = ilo; j <= next-1; j++ )
    {
      word_chstar[i] = s[j];
      i = i + 1;
    }
    word_chstar[i] = '\0';
    word = string ( word_chstar );
    delete [] word_chstar;
  }

  return word;
}
