# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "medit_io.hpp"

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

void cyl248_data ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
  int hexahedron_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    CYL248_DATA defines the data for a 3D tetrahedral mesh.
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
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
//
//    Output, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
//    of each vertex.
//
//    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.
//
//    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
//
//    Output, int EDGE_LABEL[EDGES], a label for each edge.
//
//    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
//    each triangle.
//
//    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
//
//    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
//    form each quadrilateral.
//
//    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
//    each quadrilateral.
//
//    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
//    form each tetrahedron.
//
//    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
//    each tetrahedron.
//
//    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
//    each hexahedron.
//
//    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
//
{
# define TETRAHEDRONS_SAVE 248
# define TRIANGLES_SAVE 154
# define VERTICES_SAVE 92

  int i;

  int tetrahedron_vertex_save[4*TETRAHEDRONS_SAVE] = {
    23, 1, 9, 8,
    27, 9, 23, 1,
    26, 8, 23, 9,
    26, 9, 7, 8,
    2, 9, 27, 1,
    26, 9, 10, 7,
    26, 28, 7, 10,
    11, 29, 3, 2,
    7, 6, 10, 28,
    10, 6, 31, 28,
    11, 29, 30, 3,
    11, 30, 4, 3,
    11, 30, 32, 4,
    10, 6, 5, 31,
    11, 5, 4, 32,
    19, 33, 34, 20,
    39, 22, 40, 16,
    39, 17, 36, 22,
    39, 22, 16, 17,
    40, 22, 15, 16,
    12, 19, 20, 33,
    19, 20, 34, 18,
    12, 33, 20, 35,
    38, 37, 14, 21,
    36, 22, 17, 18,
    38, 14, 15, 21,
    13, 14, 37, 21,
    12, 20, 13, 35,
    80, 32, 11, 30,
    80, 28, 10, 31,
    80, 31, 59, 28,
    80, 58, 57, 26,
    80, 28, 58, 26,
    80, 59, 58, 28,
    80, 28, 26, 10,
    80, 10, 26, 9,
    80, 9, 11, 10,
    80, 9, 26, 23,
    80, 23, 26, 57,
    80, 23, 27, 9,
    80, 23, 56, 27,
    80, 30, 11, 29,
    80, 5, 10, 11,
    80, 5, 11, 32,
    80, 5, 32, 31,
    80, 31, 10, 5,
    80, 2, 11, 9,
    80, 29, 11, 2,
    80, 2, 9, 27,
    80, 27, 29, 2,
    81, 40, 39, 22,
    81, 22, 39, 36,
    81, 18, 36, 34,
    81, 34, 20, 18,
    81, 22, 36, 18,
    81, 20, 22, 18,
    81, 37, 38, 21,
    81, 20, 33, 35,
    81, 13, 21, 20,
    81, 13, 20, 35,
    81, 13, 37, 21,
    81, 35, 37, 13,
    81, 20, 21, 22,
    81, 34, 33, 20,
    81, 21, 38, 15,
    81, 38, 40, 15,
    81, 22, 21, 15,
    81, 15, 40, 22,
    82, 60, 74, 59,
    82, 74, 25, 59,
    82, 73, 72, 58,
    82, 25, 73, 58,
    82, 59, 25, 58,
    82, 58, 72, 57,
    82, 57, 80, 58,
    82, 58, 80, 59,
    83, 71, 79, 70,
    83, 70, 76, 78,
    83, 79, 76, 70,
    83, 79, 60, 76,
    83, 82, 60, 74,
    84, 54, 64, 55,
    84, 64, 65, 55,
    84, 65, 63, 55,
    84, 65, 71, 63,
    85, 29, 62, 30,
    85, 80, 29, 30,
    85, 29, 61, 62,
    85, 78, 83, 76,
    85, 78, 76, 30,
    85, 62, 78, 30,
    85, 76, 83, 60,
    85, 76, 32, 30,
    85, 32, 80, 30,
    85, 32, 76, 60,
    85, 27, 61, 29,
    85, 80, 27, 29,
    85, 83, 82, 60,
    85, 77, 78, 62,
    85, 60, 82, 59,
    85, 59, 82, 80,
    85, 32, 60, 31,
    85, 80, 32, 31,
    85, 60, 59, 31,
    85, 59, 80, 31,
    86, 51, 68, 52,
    86, 69, 68, 51,
    86, 68, 67, 52,
    86, 52, 67, 53,
    86, 67, 66, 53,
    86, 53, 66, 54,
    87, 50, 70, 49,
    87, 71, 70, 50,
    87, 63, 71, 50,
    87, 63, 84, 71,
    87, 70, 69, 49,
    87, 71, 83, 70,
    87, 49, 69, 51,
    87, 69, 86, 51,
    88, 64, 66, 73,
    88, 72, 73, 66,
    88, 72, 82, 73,
    88, 24, 72, 66,
    88, 64, 73, 25,
    88, 73, 82, 25,
    88, 66, 64, 54,
    88, 84, 54, 64,
    88, 87, 86, 84,
    88, 67, 24, 66,
    88, 66, 86, 67,
    88, 64, 25, 65,
    88, 65, 84, 64,
    88, 25, 74, 65,
    88, 25, 82, 74,
    88, 83, 87, 71,
    88, 71, 87, 84,
    88, 82, 83, 74,
    88, 74, 83, 71,
    88, 65, 74, 71,
    88, 71, 84, 65,
    89, 86, 87, 84,
    89, 39, 48, 44,
    89, 44, 49, 43,
    89, 44, 43, 36,
    89, 44, 48, 50,
    89, 48, 63, 50,
    89, 86, 84, 54,
    89, 51, 87, 86,
    89, 44, 50, 49,
    89, 50, 87, 49,
    89, 43, 49, 51,
    89, 49, 87, 51,
    89, 39, 44, 36,
    89, 36, 81, 39,
    89, 63, 48, 47,
    89, 47, 48, 40,
    89, 46, 55, 47,
    89, 38, 46, 47,
    89, 55, 63, 47,
    89, 55, 84, 63,
    89, 43, 42, 34,
    89, 43, 51, 42,
    89, 45, 53, 54,
    89, 53, 86, 54,
    89, 45, 54, 46,
    89, 42, 52, 41,
    89, 41, 52, 53,
    89, 52, 86, 53,
    89, 42, 51, 52,
    89, 51, 86, 52,
    89, 46, 54, 55,
    89, 54, 84, 55,
    90, 56, 75, 61,
    90, 24, 75, 56,
    90, 27, 56, 61,
    90, 61, 85, 27,
    90, 75, 77, 61,
    90, 80, 82, 57,
    90, 85, 82, 80,
    90, 57, 24, 56,
    90, 72, 24, 57,
    90, 57, 82, 72,
    90, 80, 56, 27,
    90, 85, 80, 27,
    91, 85, 90, 77,
    91, 86, 87, 69,
    91, 78, 77, 69,
    91, 83, 88, 82,
    91, 90, 82, 88,
    91, 67, 88, 86,
    91, 88, 87, 86,
    91, 87, 88, 83,
    91, 83, 85, 78,
    91, 78, 85, 77,
    91, 77, 75, 68,
    91, 77, 90, 75,
    91, 69, 77, 68,
    91, 68, 86, 69,
    91, 68, 75, 67,
    91, 67, 86, 68,
    91, 24, 88, 67,
    91, 90, 88, 24,
    91, 69, 87, 70,
    91, 87, 83, 70,
    91, 75, 24, 67,
    91, 75, 90, 24,
    92, 89, 46, 45,
    92, 41, 53, 45,
    92, 89, 45, 53,
    92, 89, 53, 41,
    92, 89, 41, 42,
    92, 35, 41, 45,
    92, 33, 41, 35,
    92, 35, 81, 33,
    92, 35, 45, 37,
    92, 81, 35, 37,
    92, 34, 89, 42,
    92, 81, 89, 34,
    92, 33, 42, 41,
    92, 37, 45, 46,
    92, 37, 46, 38,
    92, 81, 37, 38,
    92, 33, 34, 42,
    92, 33, 81, 34,
    83, 74, 60, 71,
    83, 60, 79, 71,
    89, 39, 40, 48,
    89, 39, 81, 40,
    89, 36, 43, 34,
    89, 34, 81, 36,
    89, 63, 87, 50,
    89, 84, 87, 63,
    54, 88, 66, 86,
    54, 88, 86, 84,
    90, 72, 88, 24,
    90, 82, 88, 72,
    38, 47, 89, 40,
    38, 89, 81, 40,
    92, 46, 89, 38,
    92, 89, 81, 38,
    80, 23, 57, 56,
    80, 57, 90, 56,
    61, 85, 62, 77,
    61, 90, 85, 77,
    82, 85, 91, 83,
    82, 90, 91, 85,
    70, 91, 78, 83,
    70, 78, 91, 69 };
  int triangle_label_save[TRIANGLES_SAVE] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3 };
  int triangle_vertex_save[3*TRIANGLES_SAVE] = {
    12, 20, 19,
    12, 13, 20,
    19, 20, 18,
    20, 22, 18,
    22, 17, 18,
    22, 16, 17,
    13, 21, 20,
    13, 14, 21,
    14, 15, 21,
    22, 15, 16,
    22, 21, 15,
    20, 21, 22,
    1, 9, 8,
    2, 9, 1,
    9, 7, 8,
    2, 11, 9,
    11, 2, 3,
    11, 3, 4,
    9, 10, 7,
    7, 10, 6,
    10, 5, 6,
    11, 4, 5,
    5, 10, 11,
    9, 11, 10,
    23, 1, 8,
    26, 23, 8,
    26, 8, 7,
    27, 1, 23,
    2, 1, 27,
    26, 7, 28,
    7, 6, 28,
    27, 29, 2,
    29, 3, 2,
    29, 30, 3,
    30, 4, 3,
    6, 31, 28,
    6, 5, 31,
    5, 32, 31,
    5, 4, 32,
    12, 19, 33,
    19, 34, 33,
    19, 18, 34,
    12, 33, 35,
    12, 35, 13,
    18, 36, 34,
    36, 18, 17,
    35, 37, 13,
    13, 37, 14,
    38, 14, 37,
    38, 15, 14,
    39, 36, 17,
    39, 17, 16,
    38, 40, 15,
    40, 16, 15,
    39, 16, 40,
    33, 41, 35,
    33, 42, 41,
    33, 34, 42,
    36, 43, 34,
    43, 42, 34,
    39, 44, 36,
    44, 43, 36,
    35, 45, 37,
    35, 41, 45,
    37, 46, 38,
    37, 45, 46,
    38, 47, 40,
    38, 46, 47,
    39, 48, 44,
    39, 40, 48,
    47, 48, 40,
    44, 49, 43,
    44, 50, 49,
    44, 48, 50,
    43, 51, 42,
    43, 49, 51,
    42, 52, 41,
    42, 51, 52,
    41, 53, 45,
    41, 52, 53,
    45, 54, 46,
    45, 53, 54,
    46, 55, 47,
    46, 54, 55,
    30, 32, 4,
    23, 56, 27,
    23, 57, 56,
    23, 26, 57,
    28, 58, 26,
    58, 57, 26,
    31, 59, 28,
    59, 58, 28,
    32, 60, 31,
    60, 59, 31,
    27, 61, 29,
    27, 56, 61,
    29, 62, 30,
    29, 61, 62,
    55, 63, 47,
    63, 48, 47,
    48, 63, 50,
    54, 64, 55,
    64, 65, 55,
    65, 63, 55,
    53, 66, 54,
    66, 64, 54,
    52, 67, 53,
    67, 66, 53,
    51, 68, 52,
    68, 67, 52,
    49, 69, 51,
    69, 68, 51,
    50, 70, 49,
    70, 69, 49,
    63, 71, 50,
    71, 70, 50,
    65, 71, 63,
    64, 25, 65,
    64, 73, 25,
    64, 66, 73,
    67, 24, 66,
    24, 72, 66,
    72, 73, 66,
    68, 75, 67,
    75, 24, 67,
    69, 77, 68,
    77, 75, 68,
    70, 78, 69,
    78, 77, 69,
    62, 78, 30,
    78, 76, 30,
    76, 32, 30,
    32, 76, 60,
    61, 77, 62,
    77, 78, 62,
    56, 75, 61,
    75, 77, 61,
    57, 24, 56,
    24, 75, 56,
    58, 72, 57,
    72, 24, 57,
    59, 25, 58,
    25, 73, 58,
    73, 72, 58,
    60, 74, 59,
    74, 25, 59,
    25, 74, 65,
    65, 74, 71,
    70, 76, 78,
    71, 79, 70,
    79, 76, 70,
    79, 60, 76,
    74, 60, 71,
    60, 79, 71 };
  int vertex_label_save[VERTICES_SAVE] = {
    3, 3, 3, 3, 3, 3, 3, 3, 2, 2,
    2, 3, 3, 3, 3, 3, 3, 3, 3, 4,
    4, 4, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0 };
  double vertex_coordinate_save[3*VERTICES_SAVE] = {
  1.0,       0.2,         0.0,
  1.0,       0.141421,    0.141421,
  1.0,       0.0,         0.2,
  1.0,      -0.141421,    0.141421,
  1.0,      -0.2,         0.0,
  1.0,      -0.141421,   -0.141421,
  1.0,       0.0,        -0.2,
  1.0,       0.141421,   -0.141421,
  1.0,       0.066163,   -0.0302872,
  1.0,      -0.0615154,  -0.0610739,
  1.0,      -0.0306985,   0.0668017,
  0.0,       0.2,         0.0,
  0.0,       0.141421,   -0.141421,
  0.0,       0.0,        -0.2,
  0.0,      -0.141421,   -0.141421,
  0.0,      -0.2,         0.0,
  0.0,      -0.141421,    0.141421,
  0.0,       0.0,         0.2,
  0.0,       0.141421,    0.141421,
  0.0,       0.0686748,   0.0255359,
  0.0,       0.0,        -0.0865993,
  0.0,      -0.0686749,   0.0255359,
  0.8816,    0.185522,   -0.0747102,
  0.642415,  0.187806,   -0.0687668,
  0.627606, -0.0696445,  -0.187482,
  0.876431,  0.0811908,  -0.182779,
  0.881613,  0.186118,    0.0732131,
  0.872048, -0.0699008,  -0.187387,
  0.878318,  0.0844232,   0.181308,
  0.845861, -0.0716063,   0.186742,
  0.866503, -0.182493,   -0.0818307,
  0.859402, -0.186751,    0.0715813,
  0.131355,  0.18477,     0.0765501,
  0.13317,   0.077694,    0.184292,
  0.130862,  0.185301,   -0.0752567,
  0.135181, -0.0749468,   0.185426,
  0.130839,  0.0781729,  -0.18409,
  0.131856, -0.0754694,  -0.185214,
  0.135683, -0.184121,    0.0780993,
  0.134207, -0.184959,   -0.0760928,
  0.261923,  0.199982,    0.00264585,
  0.263928,  0.144161,    0.138627,
  0.268645,  0.00535339,  0.199928,
  0.272346, -0.137646,    0.145098,
  0.26108,   0.144683,   -0.138082,
  0.260772,  0.00498797, -0.199938,
  0.264253, -0.139152,   -0.143655,
  0.270288, -0.199962,    0.00389323,
  0.408181, -0.0730357,   0.186187,
  0.411818, -0.184374,    0.0774991,
  0.397539,  0.080738,    0.182979,
  0.39192,   0.185619,    0.0744699,
  0.392192,  0.184438,   -0.0773479,
  0.389194,  0.0770141,  -0.184577,
  0.38786,  -0.0747817,  -0.185493,
  0.762413,  0.199986,   -0.0023425,
  0.762987,  0.151152,   -0.13097,
  0.741526,  0.0187858,  -0.199116,
  0.746899, -0.128364,   -0.153371,
  0.720076, -0.19917,    -0.0182053,
  0.7628,    0.152219,    0.129728,
  0.763882,  0.0434475,   0.195224,
  0.399903, -0.1841,     -0.0781489,
  0.506331, -0.00579066, -0.199916,
  0.514514, -0.133894,   -0.148568,
  0.526121,  0.135152,   -0.147424,
  0.517967,  0.199953,   -0.0043215,
  0.520585,  0.147847,    0.13469,
  0.533956,  0.0124181,   0.199614,
  0.558316, -0.136902,    0.145801,
  0.549126, -0.199624,   -0.0122659,
  0.657307,  0.117735,   -0.161674,
  0.611189,  0.041829,   -0.195577,
  0.631917, -0.164669,   -0.113508,
  0.641444,  0.187001,    0.0709267,
  0.720251, -0.155557,    0.125706,
  0.647345,  0.0932963,   0.176906,
  0.677484, -0.0430068,   0.195321,
  0.635293, -0.188734,    0.0661777,
  0.888023, -0.00868364, -0.00818647,
  0.112146,  0.0,        -0.0118425,
  0.676228,  0.0124197,  -0.0856487,
  0.638436, -0.0639898,   0.0525795,
  0.452586, -0.0410297,  -0.0704842,
  0.762004, -0.0188614,   0.0693717,
  0.463368,  0.0649048,   0.0262133,
  0.473921, -0.0356443,   0.0388516,
  0.557002,  0.0123705,  -0.0932599,
  0.290986, -0.0200898,   0.00857934,
  0.7038,    0.0856777,   0.0182744,
  0.576134,  0.0436218,   0.0828782,
  0.215187,  0.080855,   -0.0314946 };

  r8vec_copy ( dim * vertices, vertex_coordinate_save, vertex_coordinate );

  i4vec_copy ( vertices, vertex_label_save, vertex_label );

  i4vec_copy ( 3 * triangles, triangle_vertex_save, triangle_vertex );

  i4vec_copy ( triangles, triangle_label_save, triangle_label );

  i4vec_copy ( 4 * tetrahedrons, tetrahedron_vertex_save, tetrahedron_vertex );

  for ( i = 0; i < tetrahedrons; i++ )
  {
    tetrahedron_label[i] = 1;
  }
  return;
# undef TETRAHEDRONS_SAVE
# undef TRIANGLES_SAVE
# undef VERTICES_SAVE
}
//****************************************************************************80

void cyl248_size ( int *dim, int *vertices, int *edges, int *triangles,
  int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

//****************************************************************************80
//
//  Purpose:
//
//    CYL248_SIZE defines the sizes for a 3D tetrahedral mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Output, int *DIM, the spatial dimension, which should be 2 or 3.
//
//    Output, int *VERTICES, the number of vertices.
//
//    Output, int *EDGES, the number of edges (may be 0).
//
//    Output, int *TRIANGLES, the number of triangles (may be 0).
//
//    Output, int *QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Output, int *TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Output, int *HEXAHEDRONS, the number of hexahedrons (may be 0).
//
{
  *dim = 3;
  *vertices = 92;
  *edges = 0;
  *triangles = 154;
  *quadrilaterals = 0;
  *tetrahedrons = 248;
  *hexahedrons = 0;

  return;
}
//****************************************************************************80

void hexahexa_2x2x2_data ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
  int hexahedron_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    HEXAHEXA_2X2X2_DATA defines the data for a 3D hexahedral mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
//
//    Output, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
//    of each vertex.
//
//    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.
//
//    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
//
//    Output, int EDGE_LABEL[EDGES], a label for each edge.
//
//    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
//    each triangle.
//
//    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
//
//    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
//    form each quadrilateral.
//
//    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
//    each quadrilateral.
//
//    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
//    form each tetrahedron.
//
//    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
//    each tetrahedron.
//
//    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
//    each hexahedron.
//
//    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
//
{
# define HEXAHEDRONS_SAVE 8
# define QUADRILATERALS_SAVE 24
# define VERTICES_SAVE 27

  int hexahedron_label_save[HEXAHEDRONS_SAVE] = {
    1, 1, 1, 1, 1, 1, 1, 1 };
  int hexahedron_vertex_save[8*HEXAHEDRONS_SAVE] = {
      1,  2,  5,  4, 10, 11, 14, 13,
      2,  3,  6,  5, 11, 12, 15, 14,
      4,  5,  8,  7, 13, 14, 17, 16,
      5,  6,  9,  8, 14, 15, 18, 17,
     10, 11, 14, 13, 19, 20, 23, 22,
     11, 12, 15, 14, 20, 21, 24, 23,
     13, 14, 17, 16, 22, 23, 26, 25,
     14, 15, 18, 17, 23, 24, 27, 26 };
  int quadrilateral_label_save[QUADRILATERALS_SAVE] = {
     1, 1, 1, 1, 2, 2, 2, 2, 3, 3,
     3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
     6, 6, 6, 6 };
  int quadrilateral_vertex_save[4*QUADRILATERALS_SAVE] = {
     1,  4,  5,  2,
     2,  5,  6,  3,
     4,  7,  8,  5,
     5,  8,  9,  6,
     1,  2, 11, 10,
     2,  3, 12, 11,
    10, 11, 20, 19,
    11, 12, 21, 20,
     3,  6, 15, 12,
     6,  9, 18, 15,
    12, 15, 24, 21,
    15, 18, 27, 24,
     7, 16, 17,  8,
     8, 17, 18,  9,
    16, 25, 26, 17,
    17, 26, 27, 18,
     1, 10, 13,  4,
     4, 13, 16,  7,
    10, 19, 22, 13,
    13, 22, 25, 16,
    19, 20, 23, 22,
    20, 21, 24, 23,
    22, 23, 26, 25,
    23, 24, 27, 26 };
  int vertex_label_save[VERTICES_SAVE] = {
    5, 2, 3, 5, 1, 3, 5, 4, 4, 5,
    2, 3, 5, 0, 3, 5, 4, 4, 6, 6,
    6, 6, 6, 6, 6, 6, 6 };
  double vertex_coordinate_save[3*VERTICES_SAVE] = {
    0.0, 0.0, 0.0,
    0.5, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 0.5, 0.0,
    0.5, 0.5, 0.0,
    1.0, 0.5, 0.0,
    0.0, 1.0, 0.0,
    0.5, 1.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 0.0, 0.5,
    0.5, 0.0, 0.5,
    1.0, 0.0, 0.5,
    0.0, 0.5, 0.5,
    0.5, 0.5, 0.5,
    1.0, 0.5, 0.5,
    0.0, 1.0, 0.5,
    0.5, 1.0, 0.5,
    1.0, 1.0, 0.5,
    0.0, 0.0, 1.0,
    0.5, 0.0, 1.0,
    1.0, 0.0, 1.0,
    0.0, 0.5, 1.0,
    0.5, 0.5, 1.0,
    1.0, 0.5, 1.0,
    0.0, 1.0, 1.0,
    0.5, 1.0, 1.0,
    1.0, 1.0, 1.0 };

  r8vec_copy ( dim * vertices, vertex_coordinate_save, vertex_coordinate );
  i4vec_copy ( vertices, vertex_label_save, vertex_label );
  i4vec_copy ( 4 * quadrilaterals, quadrilateral_vertex_save, quadrilateral_vertex );
  i4vec_copy ( quadrilaterals, quadrilateral_label_save, quadrilateral_label );
  i4vec_copy ( 8 * hexahedrons, hexahedron_vertex_save, hexahedron_vertex );
  i4vec_copy ( hexahedrons, hexahedron_label_save, hexahedron_label );

  return;
# undef HEXAHEDRONS_SAVE
# undef QUADRILATERALS_SAVE
# undef VERTICES_SAVE
}
//****************************************************************************80

void hexahexa_2x2x2_size ( int *dim, int *vertices, int *edges, int *triangles,
  int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

//****************************************************************************80
//
//  Purpose:
//
//    HEXAHEXA_2X2X2_SIZE defines the sizes for a 3D hexahedral mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Output, int *DIM, the spatial dimension, which should be 2 or 3.
//
//    Output, int *VERTICES, the number of vertices.
//
//    Output, int *EDGES, the number of edges (may be 0).
//
//    Output, int *TRIANGLES, the number of triangles (may be 0).
//
//    Output, int *QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Output, int *TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Output, int *HEXAHEDRONS, the number of hexahedrons (may be 0).
//
{
  *dim = 3;
  *vertices = 27;
  *edges = 0;
  *triangles = 0;
  *quadrilaterals = 24;
  *tetrahedrons = 0;
  *hexahedrons = 8;

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

void i4vec_zero ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO zeroes an I4VEC.
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
//    Output, int A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}
//****************************************************************************80

void mesh_data_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_DATA_PRINT prints the data of a MESH dataset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
//
//    Input, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
//    of each vertex.
//
//    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.
//
//    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
//
//    Input, int EDGE_LABEL[EDGES], a label for each edge.
//
//    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
//    each triangle.
//
//    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
//
//    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
//    form each quadrilateral.
//
//    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
//    each quadrilateral.
//
//    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
//    form each tetrahedron.
//
//    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
//    each tetrahedron.
//
//    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
//    each hexahedron.
//
//    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
//
{
  int i;
  int j;

  cout << "\n";
  cout << "  Vertices:\n";
  cout << "\n";
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < dim; i++ )
    {
      cout << "  " << vertex_coordinate[i+j*dim];
    }
    cout << "  (" << vertex_label[j] << ")\n";
  }

  if ( 0 < edges )
  {
    cout << "\n";
    cout << "  Edges:\n";
    cout << "\n";
    for ( j = 0; j < edges; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        cout << "  " << edge_vertex[i+j*2];
    }
    cout << "  (" << edge_label[j] << ")\n";
    }
  }

  if ( 0 < triangles )
  {
    cout << "\n";
    cout << "  Triangles:\n";
    cout << "\n";
    for ( j = 0; j < triangles; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        cout << "  " << triangle_vertex[i+j*3];
      }
      cout << "  (" << triangle_label[j] << ")\n";
    }
  }

  if ( 0 < quadrilaterals )
  {
    cout << "\n";
    cout << "  Quadrilaterals:\n";
    cout << "\n";
    for ( j = 0; j < quadrilaterals; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        cout << "  " << quadrilateral_vertex[i+j*4];
      }
      cout << "  (" << quadrilateral_label[j] << ")\n";
    }
  }

  if ( 0 < tetrahedrons )
  {
    cout << "\n";
    cout << "  Tetrahedrons:\n";
    cout << "\n";
    for ( j = 0; j < tetrahedrons; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        cout << "  " << tetrahedron_vertex[i+j*4];
      }
      cout << "  (" << tetrahedron_label[j] << ")\n";
    }
  }

  if ( 0 < hexahedrons )
  {
    cout << "\n";
    cout << "  Hexahedrons:\n";
    cout << "\n";
    for ( j = 0; j < hexahedrons; j++ )
    {
      for ( i = 0; i < 8; i++ )
      {
        cout << "  " << hexahedron_vertex[i+j*8];
      }
      cout << "  (" << hexahedron_label[j] << ")\n";
    }
  }
  return;
}
//****************************************************************************80

void mesh_data_read ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_DATA_READ reads data from a MESH file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the MESH file.
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals
//    (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons
//    (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons
//    (may be 0).
//
//    Output, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
//    of each vertex.
//
//    Output, int VERTEX_LABEL[VERTICES], a label for each vertex.
//
//    Output, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
//
//    Output, int EDGE_LABEL[EDGES], a label for each edge.
//
//    Output, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
//    each triangle.
//
//    Output, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
//
//    Output, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
//    form each quadrilateral.
//
//    Output, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
//    each quadrilateral.
//
//    Output, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
//    form each tetrahedron.
//
//    Output, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
//    each tetrahedron.
//
//    Output, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
//    each hexahedron.
//
//    Output, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
//
{
  int dim2;
  int edge;
  int edges2;
  int hexahedron;
  int hexahedrons2;
  int i;
  int i4vec[9];
  int ierror;
  ifstream input;
  string keyword;
  int length;
  int line_num;
  int quadrilateral;
  int quadrilaterals2;
  double r8vec[9];
  int tetrahedron;
  int tetrahedrons2;
  string text;
  int triangle;
  int triangles2;
  int vertex;
  int vertices2;
//
//  Initialize everything to nothing.
//
  i4vec_zero ( edges, edge_label );
  i4vec_zero ( 2 * edges, edge_vertex );
  i4vec_zero ( hexahedrons, hexahedron_label );
  i4vec_zero ( 8 * hexahedrons, hexahedron_vertex );
  i4vec_zero ( quadrilaterals, quadrilateral_label );;
  i4vec_zero ( 4 * quadrilaterals, quadrilateral_vertex );
  i4vec_zero ( tetrahedrons, tetrahedron_label );
  i4vec_zero ( 4 * tetrahedrons, tetrahedron_vertex );
  i4vec_zero ( triangles, triangle_label );
  i4vec_zero ( 3 * triangles, triangle_vertex );
  r8vec_zero ( dim * vertices, vertex_coordinate );
  i4vec_zero ( vertices, vertex_label );
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MESH_DATA_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Read lines til you get alphanumerics and determine a "mode"
//
  line_num = 0;
  keyword = "NONE";

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    line_num = line_num + 1;

    if ( s_len_trim ( text ) == 0 )
    {
      keyword = "NONE";
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
//
//  Remove initial blanks.
//

//
//  Expecting a keyword.
//
        if ( s_eqi ( text, "CORNERS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "DIMENSION" ) )
    {
      keyword = "DIMENSION";
    }
    else if ( s_eqi ( text, "EDGES" ) )
    {
      keyword = "EDGES";
    }
    else if ( s_eqi ( text, "END" ) )
    {
      cout << "\n";
      cout << "  END statement encountered.\n";
      break;
    }
    else if ( s_eqi ( text, "HEXAHEDRA" ) ||
              s_eqi ( text, "HEXAHEDRONS" ) )
    {
      keyword = "HEXAHEDRONS";
    }
    else if ( s_begin ( text, "MESHVERSIONFORMATTED" ) )
    {
    }
    else if ( s_eqi ( text, "NORMALATQUADRILATERALVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALATTRIANGLEVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALATVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "QUADRILATERALS" ) )
    {
      keyword = "QUADRILATERALS";
    }
    else if ( s_eqi ( text, "REQUIREDEDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "REQUIREDVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "RIDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TANGENTATEDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TANGENTS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TETRAHEDRA" ) ||
              s_eqi ( text, "TETRAHEDRONS" ) )
    {
      keyword = "TETRAHEDRONS";
    }
    else if ( s_eqi ( text, "TRIANGLES" ) )
    {
      keyword = "TRIANGLES";
    }
    else if ( s_eqi ( text, "VERTICES" ) )
    {
      keyword = "VERTICES";
    }
//
//  Presumably, numeric data to be processed by keyword.
//
    else if ( s_eqi ( keyword, "DIMENSION" ) )
    {
      dim2 = atoi ( text.c_str ( ) );
      keyword = "NONE";
    }
    else if ( s_eqi ( keyword, "EDGES" ) )
    {
      edges2 = atoi ( text.c_str ( ) );
      keyword = "EDGE_VERTEX";
      edge = 0;
    }
    else if ( s_eqi ( keyword, "EDGE_VERTEX" ) )
    {
      s_to_i4vec ( text, 3, i4vec );
      for ( i = 0; i < 2; i++ )
      {
        edge_vertex[i+edge*2] = i4vec[i];
      }
      edge_label[edge] = i4vec[2];
      edge = edge + 1;
    }
    else if ( s_eqi ( keyword, "HEXAHEDRONS" ) )
    {
      hexahedrons2 = atoi ( text.c_str ( ) );
      keyword = "HEXAHEDRON_VERTEX";
      hexahedron = 0;
    }
    else if ( s_eqi ( keyword, "HEXAHEDRON_VERTEX" ) )
    {
      s_to_i4vec ( text, 9, i4vec );
      for ( i = 0; i < 8; i++ )
      {
        hexahedron_vertex[i+hexahedron*8] = i4vec[i];
      }
      hexahedron_label[hexahedron] = i4vec[8];
      hexahedron = hexahedron + 1;
    }
    else if ( s_eqi ( keyword, "QUADRILATERALS" ) )
    {
      quadrilaterals2 = atoi ( text.c_str ( ) );
      keyword = "QUADRILATERAL_VERTEX";
      quadrilateral = 0;
    }
    else if ( s_eqi ( keyword, "QUADRILATERAL_VERTEX" ) )
    {
      s_to_i4vec ( text, 5, i4vec );
      for ( i = 0; i < 4; i++ )
      {
        quadrilateral_vertex[i+quadrilateral*4] = i4vec[i];
      }
      quadrilateral_label[quadrilateral] = i4vec[4];
      quadrilateral = quadrilateral + 1;
    }
    else if ( s_eqi ( keyword, "TETRAHEDRONS" ) )
    {
      tetrahedrons2 = atoi ( text.c_str ( ) );
      keyword = "TETRAHEDRON_VERTEX";
      tetrahedron = 0;
    }
    else if ( s_eqi ( keyword, "TETRAHEDRON_VERTEX" ) )
    {
      s_to_i4vec ( text, 5, i4vec );
      for ( i = 0; i < 4; i++ )
      {
        tetrahedron_vertex[i+tetrahedron*4] = i4vec[i];
      }
      tetrahedron_label[tetrahedron] = i4vec[4];
      tetrahedron = tetrahedron + 1;
    }
    else if ( s_eqi ( keyword, "TRIANGLES" ) )
    {
      triangles2 = atoi ( text.c_str ( ) );
      keyword = "TRIANGLE_VERTEX";
      triangle = 0;
    }
    else if ( s_eqi ( keyword, "TRIANGLE_VERTEX" ) )
    {
      s_to_i4vec ( text, 4, i4vec );
      for ( i = 0; i < 3; i++ )
      {
        triangle_vertex[i+triangle*3] = i4vec[i];
      }
      triangle_label[triangle] = i4vec[3];
      triangle = triangle + 1;
    }
    else if ( s_eqi ( keyword, "VERTICES" ) )
    {
      vertices2 = atoi ( text.c_str ( ) );
      keyword = "VERTEX_COORDINATE";
      vertex = 0;
    }
    else if ( s_eqi ( keyword, "VERTEX_COORDINATE" ) )
    {
      s_to_r8vec ( text, dim + 1, r8vec );
      for ( i = 0; i < dim; i++ )
      {
        vertex_coordinate[i+vertex*dim] = r8vec[i];
      }
      vertex_label[vertex] = ( int ) r8vec[dim];
      vertex = vertex + 1;
    }
    else if ( s_eqi ( keyword, "SKIP" ) )
    {
    }
    else
    {
      cerr << "\n";
      cerr << "MESH_DATA_READ - Fatal error!\n";
      cerr << "  Could not find keyword while reading line "
           << line_num << "\n";
      cerr << "\"" << text << "\".\n";
      exit ( 1 );
    }
  }
//
//  Close the file.
//
  input.close ( );

  cout << "\n";
  cout << "  Read " << line_num << " lines from \"" << filename << "\".\n";

  return;
}
//****************************************************************************80

void mesh_size_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_SIZE_PRINT prints the sizes of an ICE dataset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
//
{
  cout << "\n";
  cout << "  Number of dimensions = " << dim << "\n";
  cout << "  Number of vertices = " << vertices << "\n";
  cout << "  Number of edges = " << edges << "\n";
  cout << "  Number of triangles = " << triangles << "\n";
  cout << "  Number of quadrilaterals = " << quadrilaterals << "\n";
  cout << "  Number of tetrahedrons = " << tetrahedrons << "\n";
  cout << "  Number of hexahedrons = " << hexahedrons << "\n";

  return;
}
//****************************************************************************80

void mesh_size_read ( string filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_SIZE_READ reads sizes from a MESH file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the MESH file.
//
//    Output, int *DIM, the spatial dimension, which should be 2 or 3.
//
//    Output, int *VERTICES, the number of vertices.
//
//    Output, int *EDGES, the number of edges (may be 0).
//
//    Output, int *TRIANGLES, the number of triangles (may be 0).
//
//    Output, int *QUADRILATERALS, the number of quadrilaterals
//    (may be 0).
//
//    Output, int *TETRAHEDRONS, the number of tetrahedrons
//    (may be 0).
//
//    Output, int *HEXAHEDRONS, the number of hexahedrons
//    (may be 0).
//
{
  int ierror;
  ifstream input;
  string keyword;
  int length;
  int line_num;
  string text;
//
//  Initialize everything to nothing.
//
  *dim = 0;
  *vertices = 0;
  *edges = 0;
  *triangles = 0;
  *quadrilaterals = 0;
  *tetrahedrons = 0;
  *hexahedrons = 0;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "MESH_SIZE_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Read lines til you get alphanumerics and determine a "mode"
//
  line_num = 0;
  keyword = "NONE";

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    line_num = line_num + 1;

    if ( s_len_trim ( text ) == 0 )
    {
      keyword = "NONE";
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
//
//  Remove initial blanks.
//

//
//  Expecting a keyword.
//
        if ( s_eqi ( text, "CORNERS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "DIMENSION" ) )
    {
      keyword = "DIMENSION";
    }
    else if ( s_eqi ( text, "EDGES" ) )
    {
      keyword = "EDGES";
    }
    else if ( s_eqi ( text, "END" ) )
    {
      cout << "\n";
      cout << "  END statement encountered.\n";
      break;
    }
    else if ( s_eqi ( text, "HEXAHEDRA" ) ||
              s_eqi ( text, "HEXAHEDRONS" ) )
    {
      keyword = "HEXAHEDRONS";
    }
    else if ( s_begin ( text, "MESHVERSIONFORMATTED" ) )
    {
    }
    else if ( s_eqi ( text, "NORMALATQUADRILATERALVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALATTRIANGLEVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALATVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "NORMALS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "QUADRILATERALS" ) )
    {
      keyword = "QUADRILATERALS";
    }
    else if ( s_eqi ( text, "REQUIREDEDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "REQUIREDVERTICES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "RIDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TANGENTATEDGES" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TANGENTS" ) )
    {
      keyword = "SKIP";
    }
    else if ( s_eqi ( text, "TETRAHEDRA" ) ||
              s_eqi ( text, "TETRAHEDRONS" ) )
    {
      keyword = "TETRAHEDRONS";
    }
    else if ( s_eqi ( text, "TRIANGLES" ) )
    {
      keyword = "TRIANGLES";
    }
    else if ( s_eqi ( text, "VERTICES" ) )
    {
      keyword = "VERTICES";
    }
//
//  Presumably, numeric data to be processed by keyword.
//
    else if ( s_eqi ( keyword, "DIMENSION" ) )
    {
      *dim = atoi ( text.c_str ( ) );
      keyword = "NONE";
    }
    else if ( s_eqi ( keyword, "EDGES" ) )
    {
      *edges = atoi ( text.c_str ( ) );
      keyword = "EDGE_VERTEX";
    }
    else if ( s_eqi ( keyword, "EDGE_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "HEXAHEDRONS" ) )
    {
      *hexahedrons = atoi ( text.c_str ( ) );
      keyword = "HEXAHEDRON_VERTEX";
    }
    else if ( s_eqi ( keyword, "HEXAHEDRON_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "QUADRILATERALS" ) )
    {
      *quadrilaterals = atoi ( text.c_str ( ) );
      keyword = "QUADRILATERAL_VERTEX";
    }
    else if ( s_eqi ( keyword, "QUADRILATERAL_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "TETRAHEDRONS" ) )
    {
      *tetrahedrons = atoi ( text.c_str ( ) );
      keyword = "TETRAHEDRON_VERTEX";
    }
    else if ( s_eqi ( keyword, "TETRAHEDRON_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "TRIANGLES" ) )
    {
      *triangles = atoi ( text.c_str ( ) );
      keyword = "TRIANGLE_VERTEX";
    }
    else if ( s_eqi ( keyword, "TRIANGLE_VERTEX" ) )
    {
    }
    else if ( s_eqi ( keyword, "VERTICES" ) )
    {
      *vertices = atoi ( text.c_str ( ) );
      keyword = "VERTEX_COORDINATE";
    }
    else if ( s_eqi ( keyword, "VERTEX_COORDINATE" ) )
    {
    }
    else if ( s_eqi ( keyword, "SKIP" ) )
    {
    }
    else
    {
      cerr << "\n";
      cerr << "MESH_SIZE_READ - Fatal error!\n";
      cerr << "  Could not find keyword while reading line "
           << line_num << "\n";
      cerr << "\"" << text << "\".\n";
      exit ( 1 );
    }
  }
//
//  Close the file.
//
  input.close ( );

  cout << "\n";
  cout << "  Read " << line_num << " lines from \"" << filename << "\".\n";

  return;
}
//****************************************************************************80

void mesh_write ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_WRITE writes mesh data to a MESH file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be created.
//    Ordinarily, the name should include the extension ".mesh".
//
//    Input, int DIM, the spatial dimension, which should be 2 or 3.
//
//    Input, int VERTICES, the number of vertices.
//
//    Input, int EDGES, the number of edges (may be 0).
//
//    Input, int TRIANGLES, the number of triangles (may be 0).
//
//    Input, int QUADRILATERALS, the number of quadrilaterals (may be 0).
//
//    Input, int TETRAHEDRONS, the number of tetrahedrons (may be 0).
//
//    Input, int HEXAHEDRONS, the number of hexahedrons (may be 0).
//
//    Input, double VERTEX_COORDINATE[DIM*VERTICES], the coordinates
//    of each vertex.
//
//    Input, int VERTEX_LABEL[VERTICES], a label for each vertex.
//
//    Input, int EDGE_VERTEX[2*EDGES], the vertices that form each edge.
//
//    Input, int EDGE_LABEL[EDGES], a label for each edge.
//
//    Input, int TRIANGLE_VERTEX[3*TRIANGLES], the vertices that form
//    each triangle.
//
//    Input, int TRIANGLE_LABEL[TRIANGLES], a label for each triangle.
//
//    Input, int QUADRILATERAL_VERTEX[4*QUADRILATERALS], the vertices that
//    form each quadrilateral.
//
//    Input, int QUADRILATERAL_LABEL[QUADRILATERALS], a label for
//    each quadrilateral.
//
//    Input, int TETRAHEDRON_VERTEX[4*TETRAHEDRONS], the vertices that
//    form each tetrahedron.
//
//    Input, int TETRAHEDRON_LABEL[TETRAHEDRONS], a label for
//    each tetrahedron.
//
//    Input, int HEXAHEDRON_VERTEX[8*HEXAHEDRONS], the vertices that form
//    each hexahedron.
//
//    Input, int HEXAHEDRON_LABEL[HEXAHEDRONS], a label for each hexahedron.
//
{
  int i;
  int j;
  ofstream output;

  output.open ( filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "MESH_WRITE - Fatal error!\n";
    cerr << "  Unable to open output file.\n";
    exit ( 1 );
  }

  output << "MeshVersionFormatted 1\n";
  output << "#  Created by mesh_write.C\n";
//
//  Vertices.
//
  output << "\n";
  output << "Vertices\n";
  output << vertices << "\n";
  for ( j = 0; j < vertices; j++ )
  {
    for ( i = 0; i < dim; i++ )
    {
      output << "  " << vertex_coordinate[i+j*dim];
    }
    output << "  " << vertex_label[j] << "\n";
  }
//
//  Edges.
//
  if ( 0 < edges )
  {
    output << "\n";
    output << "Edges\n";
    output << edges << "\n";
    for ( j = 0; j < edges; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        output << "  " << edge_vertex[i+j*2];
    }
    output << "  " << edge_label[j] << "\n";
    }
  }
//
//  Triangles.
//
  if ( 0 < triangles )
  {
    output << "\n";
    output << "Triangles\n";
    output << triangles << "\n";
    for ( j = 0; j < triangles; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        output << "  " << triangle_vertex[i+j*3];
      }
      output << "  " << triangle_label[j] << "\n";
    }
  }
//
//  Quadrilaterals.
//
  if ( 0 < quadrilaterals )
  {
    output << "\n";
    output << "Quadrilaterals\n";
    output << quadrilaterals << "\n";
    for ( j = 0; j < quadrilaterals; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        output << "  " << quadrilateral_vertex[i+j*4];
      }
      output << "  " << quadrilateral_label[j] << "\n";
    }
  }
//
//  Tetrahedra.
//
  if ( 0 < tetrahedrons )
  {
    output << "\n";
    output << "Tetrahedra\n";
    output << tetrahedrons << "\n";
    for ( j = 0; j < tetrahedrons; j++ )
    {
      for ( i = 0; i < 4; i++ )
      {
        output << "  " << tetrahedron_vertex[i+j*4];
      }
      output << "  " << tetrahedron_label[j] << "\n";
    }
  }
//
//  Hexahedra.
//
  if ( 0 < hexahedrons )
  {
    output << "\n";
    output << "Hexahedra\n";
    output << hexahedrons << "\n";
    for ( j = 0; j < hexahedrons; j++ )
    {
      for ( i = 0; i < 8; i++ )
      {
        output << "  " << hexahedron_vertex[i+j*8];
      }
      output << "  " << hexahedron_label[j] << "\n";
    }
  }
//
//  End
//
  output << "\n";
  output << "End\n";

  output.close ( );

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

bool s_begin ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_BEGIN reports whether string 1 begins with string 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2010
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
  int n;
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
//    10 October 2014
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

string s_newline_to_null ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_NEWLINE_TO_NULL replaces carriage returns or newlines by nulls.
//
//  Discussion:
//
//    The function FGETS will read a string containing a line of text read from
//    input.  However, the string will include the linefeed character '/n', or,
//    for a PC-formatted file, the carriage return and linefeed pair '/r' + '/n'.
//
//    It may be desirable that the string not contain these characters.  The
//    easiest way to deal with this is simply to replace the first instance of
//    '/r' or '/n' by a null character, which terminates the string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be modified.
//
//    Output, string S_NEWLINE_TO_NULL, the modified string.
//
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ )
  {
//
//  Handle carriage return;
//
    if ( s[i] == '\r' )
    {
      s2[i] = '\0';
      break;
    }
//
//  Handle linefeed.
//
    if ( s[i] == '\n' )
    {
      s2[i] = '\0';
      break;
    }
  }
  return s2;
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
