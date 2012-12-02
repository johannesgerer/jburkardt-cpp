# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "ice_io.h"
# include "netcdf.hpp"

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
//    Input, int DIM, the spatial dimension, which should be 3.
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
//    Output, double VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
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

  r8vec_copy ( 3 * vertices, vertex_coordinate_save, vertex_coordinate );

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
//    Output, int *DIM, the spatial dimension, which should be 3.
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

void data_print ( int dim, int vertices, int edges, int triangles,
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
//    DATA_PRINT prints the data of a MESH dataset.
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
//    Input, int DIM, the spatial dimension, which should be 3.
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
//    Input, double VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
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
    for ( i = 0; i < 3; i++ )
    {
      cout << "  " << vertex_coordinate[i+j*3];
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

void data_read ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
  int hexahedron_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    DATA_READ reads ICE data from a NETCDF file.
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
//  Reference:
//
//    Pascal Frey,
//    MEDIT: An interactive mesh visualization software,
//    Technical Report RT-0253,
//    Institut National de Recherche en Informatique et en Automatique,
//    03 December 2001.
//
//    Russ Rew,
//    The NetCDF C++ Interface Guide,
//    Unidata Program Center, August 2008.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be read.
//    Ordinarily, the name should include the extension ".nc".
//
//    Input, int DIM, the spatial dimension, which should be 3.
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
//    Output, double VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
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
//
//  Open the file in "read only" mode.
//
  NcFile dataFile ( filename.c_str ( ), NcFile::ReadOnly );

  if ( !dataFile.is_valid ( ) )
  {
    cout << "\n";
    cout << "DATA_READ: Fatal error!\n";
    cout << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Vertices.
//
  NcVar *var_vertex_coordinate = dataFile.get_var ( "Vertex_Coordinate" );
  var_vertex_coordinate->get ( &vertex_coordinate[0], dim, vertices );

  NcVar *var_vertex_label = dataFile.get_var ( "Vertex_Label" );
  var_vertex_label->get ( &vertex_label[0], vertices );
//
//  Edges.
//
  if ( 0 < edges )
  {
    NcVar *var_edge_vertex = dataFile.get_var ( "Edge_Vertex" );
    var_edge_vertex->get ( &edge_vertex[0], 2, edges );

    NcVar *var_edge_label = dataFile.get_var ( "Edge_Label" );
    var_edge_label->get ( &edge_label[0], edges );
  }
//
//  Triangles.
//
  if ( 0 < triangles )
  {
    NcVar *var_triangle_vertex = dataFile.get_var ( "Triangle_Vertex" );
    var_triangle_vertex->get ( &triangle_vertex[0], 3, triangles );

    NcVar *var_triangle_label = dataFile.get_var ( "Triangle_Label" );
    var_triangle_label->get ( &triangle_label[0], triangles );
  }
//
//  Quadrilaterals.
//
  if ( 0 < quadrilaterals )
  {
    NcVar *var_quadrilateral_vertex = dataFile.get_var ( "Quadrilateral_Vertex" );
    var_quadrilateral_vertex->get ( &quadrilateral_vertex[0], 4, quadrilaterals );

    NcVar *var_quadrilateral_label = dataFile.get_var ( "Quadrilateral_Label" );
    var_quadrilateral_label->get ( &quadrilateral_label[0], quadrilaterals );
  }
//
//  Tetrahedrons.
//
  if ( 0 < tetrahedrons )
  {
    NcVar *var_tetrahedron_vertex = dataFile.get_var ( "Tetrahedron_Vertex" );
    var_tetrahedron_vertex->get ( &tetrahedron_vertex[0], 4, tetrahedrons );

    NcVar *var_tetrahedron_label = dataFile.get_var ( "Tetrahedron_Label" );
    var_tetrahedron_label->get ( &tetrahedron_label[0], tetrahedrons );
  }
//
//  Hexahedrons.
//
  if ( 0 < hexahedrons )
  {
    NcVar *var_hexahedron_vertex = dataFile.get_var ( "Hexahedron_Vertex" );
    var_hexahedron_vertex->get ( &hexahedron_vertex[0], 8, hexahedrons );

    NcVar *var_hexahedron_label = dataFile.get_var ( "Hexahedron_Label" );
    var_hexahedron_label->get ( &hexahedron_label[0], hexahedrons );
  }
//
//  Close the file.
//
  dataFile.close ( );

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
//    Input, int DIM, the spatial dimension, which should be 3.
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
//    Output, double VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
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

  r8vec_copy ( 3 * vertices, vertex_coordinate_save, vertex_coordinate );
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
//    Output, int *DIM, the spatial dimension, which should be 3.
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

void ice_write ( std::string filename, int dim, int vertices, int edges,
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
//    ICE_WRITE writes 3D ICE sizes and data to a NETCDF file.
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
//    Ordinarily, the name should include the extension ".nc".
//
//    Input, int DIM, the spatial dimension, which should be 3.
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
//    Input, double VERTEX_COORDINATE[3*VERTICES], the XYZ coordinates
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
  NcDim *dim_dimension;
  NcDim *dim_edges;
  NcDim *dim_eight;
  NcDim *dim_four;
  NcDim *dim_hexahedrons;
  NcDim *dim_quadrilaterals;
  NcDim *dim_tetrahedrons;
  NcDim *dim_three;
  NcDim *dim_triangles;
  NcDim *dim_two;
  NcDim *dim_vertices;
  NcVar *var_edge_vertex;
  NcVar *var_edge_label;
  NcVar *var_hexahedron_vertex;
  NcVar *var_hexahedron_label;
  NcVar *var_quadrilateral_vertex;
  NcVar *var_quadrilateral_label;
  NcVar *var_tetrahedron_vertex;
  NcVar *var_tetrahedron_label;
  NcVar *var_triangle_vertex;
  NcVar *var_triangle_label;
  NcVar *var_vertex_coordinate;
  NcVar *var_vertex_label;
//
//  Create the file.
//
  NcFile dataFile ( filename.c_str ( ), NcFile::Replace );

  if ( !dataFile.is_valid ( ) )
  {
    cout << "\n";
    cout << "ICE_WRITE - Fatal error!\n";
    cout << "  Could not open the file.\n";
    exit ( 1 );
  }
//
//  Dimension information.
//
  dim_dimension = dataFile.add_dim ( "Dimension", dim );
  dim_vertices = dataFile.add_dim ( "Vertices", vertices );

  if ( 0 < edges )
  {
    dim_edges = dataFile.add_dim ( "Edges", edges );
  }

  if ( 0 < triangles )
  {
    dim_triangles = dataFile.add_dim ( "Triangles", triangles );
  }

  if ( 0 < quadrilaterals )
  {
    dim_quadrilaterals = dataFile.add_dim ( "Quadrilaterals", quadrilaterals );
  }

  if ( 0 < tetrahedrons )
  {
    dim_tetrahedrons = dataFile.add_dim ( "Tetrahedrons", tetrahedrons );
  }

  if ( 0 < hexahedrons )
  {
    dim_hexahedrons = dataFile.add_dim ( "Hexahedrons", hexahedrons );
  }

  dim_two = dataFile.add_dim ( "Two", 2 );
  dim_three = dataFile.add_dim ( "Three", 3 );
  dim_four = dataFile.add_dim ( "Four", 4 );
  dim_eight = dataFile.add_dim ( "Eight", 8 );
//
//  Define variables.
//
  var_vertex_coordinate      = dataFile.add_var ( "Vertex_Coordinate",    ncDouble, dim_three, dim_vertices );
  var_vertex_label           = dataFile.add_var ( "Vertex_Label",         ncInt,               dim_vertices );

  if ( 0 < edges )
  {
    var_edge_vertex          = dataFile.add_var ( "Edge_Vertex",          ncInt,    dim_two,   dim_edges );
    var_edge_label           = dataFile.add_var ( "Edge_Label",           ncInt,               dim_edges );
  }

  if ( 0 < triangles )
  {
    var_triangle_vertex      = dataFile.add_var ( "Triangle_Vertex",      ncInt, dim_three, dim_triangles );
    var_triangle_label       = dataFile.add_var ( "Triangle_Label",       ncInt,            dim_triangles );
  }

  if ( 0 < quadrilaterals )
  {
    var_quadrilateral_vertex = dataFile.add_var ( "Quadrilateral_Vertex", ncInt, dim_four, dim_quadrilaterals );
    var_quadrilateral_label  = dataFile.add_var ( "Quadrilateral_Label",  ncInt,           dim_quadrilaterals );
  }

  if ( 0 < tetrahedrons )
  {
    var_tetrahedron_vertex   = dataFile.add_var ( "Tetrahedron_Vertex",   ncInt, dim_four,  dim_tetrahedrons );
    var_tetrahedron_label    = dataFile.add_var ( "Tetrahedron_Label",    ncInt,            dim_tetrahedrons );
  }

  if ( 0 < hexahedrons )
  {
    var_hexahedron_vertex    = dataFile.add_var ( "Hexahedron_Vertex",    ncInt, dim_eight, dim_hexahedrons );
    var_hexahedron_label     = dataFile.add_var ( "Hexahedron_Label",     ncInt,            dim_hexahedrons );
  }
//
//  Write the data.
//
  var_vertex_coordinate->put      ( &vertex_coordinate[0],    3, vertices );
  var_vertex_label->put           ( &vertex_label[0],            vertices );
  if ( 0 < edges )
  {
    var_edge_vertex->put          ( &edge_vertex[0],          2, edges );
    var_edge_label->put           ( &edge_label[0],              edges );
  }
  if ( 0 < triangles )
  {
    var_triangle_vertex->put      ( &triangle_vertex[0],      3, triangles );
    var_triangle_label->put       ( &triangle_label[0],          triangles );
  }
  if ( 0 < quadrilaterals )
  {
    var_quadrilateral_vertex->put ( &quadrilateral_vertex[0], 4, quadrilaterals );
    var_quadrilateral_label->put  ( &quadrilateral_label[0],     quadrilaterals );
  }
  if ( 0 < tetrahedrons )
  {
    var_tetrahedron_vertex->put   ( &tetrahedron_vertex[0],   4, tetrahedrons );
    var_tetrahedron_label->put    ( &tetrahedron_label[0],       tetrahedrons );
  }
  if ( 0 < hexahedrons )
  {
    var_hexahedron_vertex->put    ( &hexahedron_vertex[0],      8, hexahedrons );
    var_hexahedron_label->put     ( &hexahedron_label[0],          hexahedrons );
  }
//
//  Close the file.
//
  dataFile.close ( );

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

void size_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons )

//****************************************************************************80
//
//  Purpose:
//
//    SIZE_PRINT prints the sizes of an ICE dataset.
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
//    Input, int DIM, the spatial dimension, which should be 3.
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
//*****************************************************************************80

void size_read ( string filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons )

//*****************************************************************************80
//
//  Purpose:
//
//    SIZE_READ reads ICE sizes from a NETCDF file.
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
//    Input, string FILENAME, the name of the file to be read.
//    Ordinarily, the name should include the extension ".nc".
//
//    Output, int *DIM, the spatial dimension, which should be 3.
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
  NcDim *dim_dimension;
  NcDim *dim_edges;
  NcDim *dim_eight;
  NcDim *dim_four;
  NcDim *dim_hexahedrons;
  NcToken dim_name;
  int dim_num;
  NcDim *dim_quadrilaterals;
  NcDim *dim_tetrahedrons;
  NcDim *dim_three;
  NcDim *dim_triangles;
  NcDim *dim_two;
  NcDim *dim_vertices;
  NcDim *dim_pointer;
  int i;
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
//  Open the file in "read only" mode.
//
  NcFile dataFile ( filename.c_str ( ), NcFile::ReadOnly );

  if ( !dataFile.is_valid ( ) )
  {
    cout << "\n";
    cout << "SIZE_READ: Fatal error!\n";
    cout << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Get the dimension information.
//
//  I would much prefer to write "0" as the size of certain dimensions, but I am not
//  allowed to, so I simply omit them from the file.
//
//  Therefore, when I open the file and try to determine dimensions, some dimensions
//  are "missing", which I would have presumed I could discover painlessly by asking
//  for pointers to them, and getting NULLs back.  But that doesn't seem to work either.
//
//  So my bonehead backup is simply to read all the dimensions by index, retrieve
//  their names, and see what I find.
//
  dim_num = dataFile.num_dims ( );

  for ( i = 0; i < dim_num; i++ )
  {
    dim_pointer = dataFile.get_dim ( i );
    dim_name = dim_pointer->name ( );

    if ( !strcmp ( dim_name, "Dimension" ) )
    {
      *dim = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Vertices" ) )
    {
      *vertices = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Edges" ) )
    {
      *edges = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Triangles" ) )
    {
      *triangles = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Quadrilaterals" ) )
    {
      *quadrilaterals = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Tetrahedrons" ) )
    {
      *tetrahedrons = dim_pointer->size ( );
    }
    else if ( !strcmp ( dim_name, "Hexahedrons" ) )
    {
      *hexahedrons = dim_pointer->size ( );
    }
    else
    {
      cout << "  Ignoring information about dimension \"" << dim_name << "\".\n";
    }
  }
//
//  Close the file.
//
  dataFile.close ( );

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
