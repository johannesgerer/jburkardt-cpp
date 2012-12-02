# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

#include "grf_io.hpp"

//****************************************************************************80

bool ch_eqi ( char c1, char c2 )

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

void grf_data_print ( int node_num, int edge_num, int edge_pointer[], 
  int edge_data[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_DATA_PRINT prints the data of a GRF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Skiena,
//    Implementing Discrete Mathematics,
//    Combinatorics and Graph Theory with Mathematica,
//    Addison-Wesley, 1990.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int EDGE_NUM, the number of edges.
//
//    Input, int EDGE_POINTER[NODE_NUM+1], pointers to 
//    the beginning of edge data for each node.
//
//    Input, int EDGE_DATA[EDGE_NUM], the edge data.
//
//    Input, double XY[2*NODE_NUM], the node coordinates.
//
{
  int edge;
  int node;

  cout << "\n";
  cout << "  Edge pointers:\n";
  cout << "\n";
  cout << "  Node     First      Last\n";
  cout << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(4) << node
         << "  " << setw(8) << edge_pointer[node]
         << "  " << setw(8) << edge_pointer[node+1] - 1 << "\n";
  }

  cout << "\n";
  cout << "  Edge data:\n";
  cout << "\n";
  cout << "  Node     Adjacent nodes\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(4) << node;
    for ( edge = edge_pointer[node]; edge <= edge_pointer[node+1] - 1; edge++ )
    {
      cout << "  " << setw(8) <<  edge_data[edge];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Node        X          Y\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(4) << node
         << "  " << setw(10) << xy[0+node*2]
         << "  " << setw(10) << xy[1+node*2] << "\n";
  }

  return;
}
//****************************************************************************80

void grf_data_read ( string input_filename, int node_num, int edge_num, 
  int edge_pointer[], int edge_data[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_DATA_READ reads the data of a GRF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Skiena,
//    Implementing Discrete Mathematics,
//    Combinatorics and Graph Theory with Mathematica,
//    Addison-Wesley, 1990.
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the file.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int EDGE_NUM, the number of edges.
//
//    Output, int EDGE_POINTER[NODE_NUM+1], pointers to 
//    the beginning of edge data for each node.
//
//    Output, int EDGE_DATA[EDGE_NUM], the edge data.
//
//    Output, double XY[2*NODE_NUM], the node coordinates.
//
{
  int edge;
  int i;
  ifstream input_unit;
  int n;
  int node;
  int node_i;
  int node_j;
  char text[255];
  int text_begin;
  int text_end;
  double xval;
  double yval;

  for ( edge = 0; edge < edge_num; edge++ )
  {
    edge_data[edge] = -1;
  }

  for ( node = 0; node < node_num + 1; node++ )
  {
    edge_pointer[node] = -1;
  }

  input_unit.open ( input_filename.c_str() );

  if ( !input_unit )
  {
    cout << "\n";
    cout << "GRF_DATA_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }
//
//  Read a line.  If it's a blank or comment, skip it.
//  Otherwise, count the number of "words", and then reread it.
//
  edge = 0;
  node = 0;
  edge_pointer[0] = 0;

  for ( node = 0; node < node_num; node++ )
  {
    text_begin = input_unit.tellg ( );

    input_unit.getline ( text, sizeof ( text ) );

    if ( input_unit.eof ( ) )
    {
      cout << "\n";
      cout << "GRF_DATA_READ - Fatal error!\n";
      cout << "  Unexpected end of information;\n";
      exit ( 1 );
    }

    if ( s_len_trim ( text ) <= 0 )
    {
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }

    n = s_word_count ( text );

    if ( n < 3 )
    {
      cout << "\n";
      cout << "GRF_DATA_READ - Fatal error!\n";
      cout << "  Record has less than 3 items.\n";
      exit ( 1 );
    }
    text_end = input_unit.tellg ( );
//
//  Back up and reread the line.
//
    input_unit.seekg ( text_begin );

    input_unit >> node_i >> xval >> yval;

    edge_pointer[node+1] = edge_pointer[node] + n - 3;

    xy[0+node*2] = xval;
    xy[1+node*2] = yval;

    for ( i = 0; i < n - 3; i++ )
    {
      input_unit >> node_j;
      edge_data[edge] = node_j;
      edge = edge + 1;
    }

    input_unit.seekg ( text_end );
  }

  input_unit.close ( );

  return;
}

//****************************************************************************80

void grf_data_write ( ofstream &output_unit, int node_num, int edge_num, 
  int edge_pointer[], int edge_data[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_DATA_WRITE writes the data of a GRF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Skiena,
//    Implementing Discrete Mathematics,
//    Combinatorics and Graph Theory with Mathematica,
//    Addison-Wesley, 1990.
//
//  Parameters:
//
//    Input, ofstream &OUTPUT_UNIT, a pointer to the GRF file.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int EDGE_NUM, the number of edges.
//
//    Input, int EDGE_POINTER[NODE_NUM+1], pointers to 
//    the beginning of edge data for each node.
//
//    Input, int EDGE_DATA[EDGE_NUM], the edge data.
//
//    Input, double XY[2*NODE_NUM], the node coordinates.
//
{
  int edge;
  int node;

  for ( node = 0; node < node_num; node++ )
  {
    output_unit 
      << "  " << setw(4) << node+1
      << "  " << setw(10) << xy[0+node*2]
      << "  " << setw(10) << xy[1+node*2];
    for ( edge = edge_pointer[node]; edge <= edge_pointer[node+1] - 1; edge++ )
    {
      output_unit << "  " << setw(8) <<  edge_data[edge];
    }
    output_unit << "\n";
  }

  return;
}
//****************************************************************************80

void grf_example ( int node_num, int edge_num, int edge_pointer[], 
  int edge_data[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_EXAMPLE sets up a GRF example.
//
//  Discussion:
//
//    The example is known as the Coxeter graph.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Skiena,
//    Implementing Discrete Mathematics,
//    Combinatorics and Graph Theory with Mathematica,
//    Addison-Wesley, 1990.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int EDGE_NUM, the number of edges.
//
//    Output, int EDGE_POINTER[NODE_NUM+1], pointers to 
//    the beginning of edge data for each node.
//
//    Output, int EDGE_DATA[EDGE_NUM], the edge data.
//
//    Output, double XY[2*NODE_NUM], the node coordinates.
//
{
# define EDGE_NUM 84
# define NODE_NUM 28

  int EDGE_DATA[EDGE_NUM] = {
     8,   2,   3, 
    14,   1,   5, 
     9,   4,   1, 
    10,   7,   3, 
    13,   2,   6, 
    12,   5,   7, 
    11,   6,   4,
    25,  20,   1, 
    24,  21,   3, 
    23,  15,   4, 
    22,  16,   7, 
    28,  17,   6, 
    27,  18,   5, 
    26,  19,   2, 
    10,  18,  19, 
    11,  19,  20, 
    12,  21,  20, 
    13,  15,  21, 
    14,  16,  15, 
     8,  17,  16, 
     9,  18,  17, 
    11,  27,  24, 
    10,  28,  25, 
     9,  26,  22, 
     8,  23,  27, 
    14,  24,  28, 
    13,  25,  22, 
    12,  26,  23 };

  int EDGE_POINTER[NODE_NUM+1] = {
    0,  3,  6,  9, 12, 15, 18, 21, 24, 27, 
   30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 
   60, 63, 66, 69, 72, 75, 78, 81, 84 };

  double XY[2*NODE_NUM] = {
    0.412,   0.984, 
    0.494,   0.984, 
    0.366,   0.926, 
    0.388,   0.862, 
    0.546,   0.926, 
    0.518,   0.860, 
    0.458,   0.818, 
    0.152,   0.684, 
    0.264,   0.682, 
    0.354,   0.680, 
    0.458,   0.670, 
    0.554,   0.672, 
    0.658,   0.668, 
    0.774,   0.692, 
    0.164,   0.450, 
    0.228,   0.448, 
    0.274,   0.390,  
    0.242,   0.330, 
    0.194,   0.278, 
    0.146,   0.328, 
    0.102,   0.390, 
    0.668,   0.472, 
    0.638,   0.416,
    0.656,   0.334, 
    0.714,   0.270, 
    0.798,   0.326, 
    0.830,   0.408, 
    0.754,   0.466 };

  i4vec_copy ( edge_num, EDGE_DATA, edge_data );
  i4vec_copy ( node_num + 1, EDGE_POINTER, edge_pointer );
  r8vec_copy ( 2 * node_num, XY, xy );

  return;
# undef EDGE_NUM
# undef NODE_NUM
}
//****************************************************************************80

void grf_example_size ( int *node_num, int *edge_num )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_EXAMPLE_SIZE sizes a GRF example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Skiena,
//    Implementing Discrete Mathematics,
//    Combinatorics and Graph Theory with Mathematica,
//    Addison-Wesley, 1990.
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *EDGE_NUM, the number of edges.
//
{
  *node_num = 28;
  *edge_num = 84;

  return;
}
//****************************************************************************80

void grf_header_print ( int node_num, int edge_num )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_HEADER_PRINT prints the header of a GRF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Skiena,
//    Implementing Discrete Mathematics,
//    Combinatorics and Graph Theory with Mathematica,
//    Addison-Wesley, 1990.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int EDGE_NUM, the number of edges.
//
{
  cout << "\n";
  cout << "  The number of nodes NODE_NUM = " << node_num << "\n";
  cout << "  The number of edges EDGE_NUM = " << edge_num << "\n";

  return;
}
//****************************************************************************80

void grf_header_read ( string input_filename, int *node_num, int *edge_num )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_HEADER_READ reads the header of a GRF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Skiena,
//    Implementing Discrete Mathematics,
//    Combinatorics and Graph Theory with Mathematica,
//    Addison-Wesley, 1990.
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the file.
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *EDGE_NUM, the number of edges.
//
{
  ifstream input_unit;
  int n;
  char text[255];

  *edge_num = 0;
  *node_num = 0;

  input_unit.open ( input_filename.c_str() );

  if ( !input_unit )
  {
    cout << "\n";
    cout << "GRF_HEADER_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    input_unit.getline ( text, sizeof ( text ) );

    if ( input_unit.eof ( ) )
    {
      break;
    }

    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }

    n = s_word_count ( text );

    if ( n < 3 )
    {
      cout << "\n";
      cout << "GRF_HEADER_READ - Fatal error!\n";
      cout << "  Illegal record has less than 3 data items\n";
      exit ( 1 );
    }
    *edge_num = *edge_num + n - 3;
    *node_num = *node_num + 1;
  }

  input_unit.close ( );

  return;
}
//****************************************************************************80

void grf_header_write ( string output_filename, ofstream &output_unit,
  int node_num, int edge_num )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_HEADER_WRITE writes the header of a GRF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the output file.
//
//    Input, ofstream &OUTPUT_UNIT, the output file unit number.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int EDGE_NUM, the number of edges.
//
{
//
//  Write the header.
//
  output_unit << "#  " << output_filename << "\n";
  output_unit << "#  created by grf_io::grf_header_write.C\n";
  output_unit << "#\n";
  output_unit << "#  Number of nodes  = " << node_num << "\n";
  output_unit << "#  Number of edges =  " << edge_num << "\n";
  output_unit << "#\n";

  return;
}
//****************************************************************************80

void grf_write ( string output_filename, int node_num, int edge_num, 
  int edge_pointer[], int edge_data[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRF_WRITE writes a GRF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the file
//    to which the data should be written.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int EDGE_NUM, the number of edges.
//
//    Input, int EDGE_POINTER[NODE_NUM+1], pointers to the
//    first edge item for each node.
//
//    Input, int EDGE_DATA[EDGE_NUM], indices of adjacent nodes.
//
//    Input, double XY[2*NODE_NUM], the node coordinates.
//
{
  ofstream output_unit;
//
//  Open the output file.
//
  output_unit.open ( output_filename.c_str() );

  if ( !output_unit )
  {
    cout << "\n";
    cout << "GRF_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << output_filename << "\".\n";
    exit ( 1 );
  }
//
//  Write the header.
//
  if ( false )
  {
    grf_header_write ( output_filename, output_unit, node_num, edge_num );
  }
//
//  Write the data.
//
  grf_data_write ( output_unit, node_num, edge_num, edge_pointer, 
    edge_data, xy );
//
//  Close the file.
//
  output_unit.close ( );

  return;
}
//****************************************************************************80

void i4vec_copy ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY copies an I4VEC.
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
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int A2[N], the copy of A1.
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

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
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
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
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
//    04 October 2003
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
