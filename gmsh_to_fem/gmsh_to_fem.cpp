# include <cmath>
# include <cstdlib>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <cstring>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void gmsh_data_read ( string gmsh_filename, int node_dim, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] );
void gmsh_size_read ( string gmsh_filename, int &node_num, int &node_dim,
  int &element_num, int &element_order );
void i4mat_write ( string output_filename, int m, int n, int table[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void r8mat_write ( string output_filename, int m, int n, double table[] );
bool s_begin ( string s1, string s2 );
int s_len_trim ( string s );
int s_to_i4 ( string s, int &last, bool &error );
double s_to_r8 ( string s, int &lchar, bool &error );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for GMSH_TO_FEM.
//
//  Discussion:
//
//    GMSH_TO_FEM converts mesh data from GMSH to FEM format.
//
//  Usage:
//
//    gmsh_to_fem prefix
//
//    where 'prefix' is the common filename prefix:
//
//    * 'prefix'.msh contains the GMSH mesh file.
//    * 'prefix'_nodes.txt will contain the node coordinates.
//    * 'prefix'_elements.txt will contain the element node connectivity.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *element_node;
  int element_num;
  int element_order;
  string fem_element_filename;
  string fem_node_filename;
  string gmsh_filename;
  int m;
  int node_num;
  double *node_x;
  string  prefix;

  timestamp ( );
  cout << "\n";
  cout << "GMSH_TO_FEM\n";
  cout << "  C++ version:\n";
  cout << "  Read a mesh description created by GMSH:\n";
  cout << "  * \"prefix\".msh, the GMSH mesh file.\n";
  cout << "  Write out two simple FEM files listing nodes and elements.\n";
  cout << "  * \"prefix\"_nodes.txt, node coordinates.\n";
  cout << "  * \"prefix\"_elements.txt, element connectivity.\n";
//
//  Get the filename prefix.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
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
  gmsh_filename = prefix + ".msh";
  fem_node_filename = prefix + "_nodes.txt";
  fem_element_filename = prefix + "_elements.txt";
//
//  Read GMSH sizes.
//
  gmsh_size_read ( gmsh_filename, node_num, m, element_num, element_order );
//
//  Report sizes.
//
  cout << "\n";
  cout << "  Size information from GMSH:\n";
  cout << "  Spatial dimension M = " << m << "\n";
  cout << "  Number of nodes NODE_NUM = " << node_num << "\n";
  cout << "  Number of elements ELEMENT_NUM = " << element_num << "\n";
  cout << "  Element order ELEMENT_ORDER = " << element_order << "\n";
//
//  Allocate memory.
//
  node_x = new double[m*node_num];
  element_node = new int[element_order*element_num];
//
//  Read GMSH data.
//
  gmsh_data_read ( gmsh_filename, m, node_num, node_x, element_order, 
    element_num, element_node );
//
//  Write FEM data.
//
  r8mat_write ( fem_node_filename, m, node_num, node_x );

  i4mat_write ( fem_element_filename, element_order, element_num, 
    element_node );
//
//  Free memory.
//
  delete [] element_node;
  delete [] node_x;
//
//  Terminate.
//
  cout << "\n";
  cout << "GMSH_TO_FEM:\n";
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

void gmsh_data_read ( string gmsh_filename, int node_dim, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_DATA_READ reads data from a GMSH file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, character *GMSH_FILENAME, the GMSH filename.
//
//    Input, int NODE_DIM, the spatial dimension.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_X[NODE_DIM*NODE_NUM], the node coordinates.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
//    the nodes that make up each element.
//
{
  int i;
  int i4_dummy;
  bool ierror;
  int indx;
  ifstream input;
  int j;
  int k;
  int length;
  int level;
  string text;
  double x;
  double y;
  double z;

  input.open ( gmsh_filename.c_str ( ) );

  if ( ! input )
  {
    cerr << "\n";
    cerr << "GMSH_DATA_READ - Fatal error!\n";
    cerr << "  Could not open input file \"" << gmsh_filename << "\"\n";
    exit ( 1 );
  }

  level = 0;
 
  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text, "$Nodes" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      i4_dummy = s_to_i4 ( text, length, ierror );
      level = 2;
      j = 0;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text, "$EndNodes" ) )
      {
        break;
      }
      else
      {
        i4_dummy = s_to_i4 ( text, length, ierror );
        text.erase ( 0, length );
        for ( i = 0; i < node_dim; i++ )
        {
          x = s_to_r8 ( text, length, ierror );
          text.erase ( 0, length );
          node_x[i+j*node_dim] = x;
        }
        j = j + 1;
      }
    }
  }
//
//  Now read element information.
//
  level = 0;
  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text, "$Elements" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      i4_dummy = s_to_i4 ( text, length, ierror );
      level = 2;
      j = 0;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text, "$EndElements" ) )
      {
        break;
      }
      else
      {
        for ( k = 1; k <= 5; k++ )
        {
          i4_dummy = s_to_i4 ( text, length, ierror );
          text.erase ( 0, length );
        }
        for ( i = 0; i < element_order; i++ )
        {
          k = s_to_i4 ( text, length, ierror );
          text.erase ( 0, length );
          element_node[i+j*element_order] = k;
        }
        j = j + 1;
      }
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void gmsh_size_read ( string gmsh_filename, int &node_num, int &node_dim,
  int &element_num, int &element_order )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_SIZE_READ reads sizes from a GMSH file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string GMSH_FILENAME, the GMSH filename.
//
//    Output, int &NODE_NUM, the number of nodes.
//
//    Output, int &NODE_DIM, the spatial dimension.
//
//    Output, int &ELEMENT_NUM, the number of elements.
//
//    Output, int &ELEMENT_ORDER, the order of the elements.
//
{
  bool ierror;
  int indx;
  ifstream input;
  int k;
  int length;
  int level;
  const double r8_big = 1.0E+30;
  string text;
  double x;
  double x_max;
  double x_min;
  double y;
  double y_max;
  double y_min;
  double z;
  double z_max;
  double z_min;

  node_num = 0;
  node_dim = 0;

  x_max = - r8_big;
  x_min = + r8_big;
  y_max = - r8_big;
  y_min = + r8_big;
  z_max = - r8_big;
  z_min = + r8_big;

  input.open ( gmsh_filename.c_str ( ) );

  if ( ! input )
  {
    cerr << "\n";
    cerr << "GMSH_SIZE_READ - Fatal error!\n";
    cerr << "  Could not open input file \"" << gmsh_filename << "\"\n";
    exit ( 1 );
  }

  level = 0;
 
  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text, "$Nodes" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      node_num = s_to_i4 ( text, length, ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text, "$EndNodes" ) )
      {
        break;
      }
      else
      {
        indx = s_to_i4 ( text, length, ierror );
        text.erase ( 0, length );

        x = s_to_r8 ( text, length, ierror );
        x_min = r8_min ( x_min, x );
        x_max = r8_max ( x_max, x );
        text.erase ( 0, length );

        y = s_to_r8 ( text, length, ierror );
        y_min = r8_min ( y_min, y );
        y_max = r8_max ( y_max, y );
        text.erase ( 0, length );

        z = s_to_r8 ( text, length, ierror);
        z_min = r8_min ( z_min, z );
        z_max = r8_max ( z_max, z );
        text.erase ( 0, length );
      }
    }
  }
//
//  Make a very simple guess as to the dimensionality of the data.
//
  node_dim = 3;
  if ( z_max == z_min )
  {
    node_dim = 2;
    if ( y_max == y_min )
    {
      node_dim = 1;
    }
  }
//
//  Now read element information.
//
  level = 0;

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text, "$Elements" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      element_num = s_to_i4 ( text, length, ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text, "$EndElements" ) )
      {
        break;
      }
      else
      {
        k = 0;
        for ( ; ; )
        {
          indx = s_to_i4 ( text, length, ierror );
          text.erase ( 0, length );
          if ( ierror != 0 )
          {
            break;
          }
          k = k + 1;
        }
        element_order = k - 5;
        break;
      }
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void i4mat_write ( string output_filename, int m, int n, int table[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_WRITE writes an I4MAT file with no header.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
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
//    Input, int TABLE[M*N], the data.
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
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << table[i+j*m];
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

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Discussion:
//
//    The C++ math library provides the function fmax() which is preferred.
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
//  Discussion:
//
//    The C++ math library provides the function fmin() which is preferred.
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
//    Input, double TABLE[M*N], the data.
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
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << table[i+j*m];
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

bool s_begin ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_BEGIN reports whether string 1 begins with string 2, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, string S2, two strings.
//
//    Output, bool S_BEGIN, is true if S1 is the same as S2 up to
//    the end of S2, and false otherwise.
//
{
  int i;
  int n1;
  int n2;

  n1 = s1.length ( );
  n2 = s2.length ( );

  if ( n1 < n2 )
  {
    return false;
  }

  for ( i = 0; i < n2; i++ )
  {
    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return false;
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
    if ( s[n-1] != ' ' && s[n-1] != '\n' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

int s_to_i4 ( string s, int &last, bool &error )

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
//    Output, int &LAST, the last character of S used to make IVAL.
//
//    Output, bool &ERROR is TRUE if an error occurred.
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

  error = false;
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
        error = true;
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
        error = true;
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
        ival = 10 * ival + c - '0';
      }
      else
      {
        ival = isgn * ival;
        last = i - 1;
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
    last = s_len_trim ( s );
  }
  else
  {
    error = true;
    last = 0;
  }

  return ival;
}
//****************************************************************************80

double s_to_r8 ( string s, int &lchar, bool &error )

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
//    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
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
//    Output, int &LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool &ERROR, is true if an error occurred.
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
  error = false;
  r = 0.0;
  lchar = -1;
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
    c = s[lchar+1];
    lchar = lchar + 1;
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
        lchar = lchar + 1;
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
    if ( iterm == 1 || nchar <= lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && lchar + 1 == nchar )
  {
    lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    error = true;
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
//**********************************************************************

void timestamp ( )

//**********************************************************************
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
