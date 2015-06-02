# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <fstream>
# include <cstring>

using namespace std;

# include "square_hex_grid.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SQUARE_HEX_GRID_PRB.
//
//  Discussion:
//
//    SQUARE_HEX_GRID_PRB tests the SQUARE_HEX_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SQUARE_HEX_GRID_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SQUARE_HEX_GRID library.\n";
//
//  Tests for unit square.
//
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Tests for arbitrary size coordinate box.
//
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );
  test11 ( );
  test12 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SQUARE_HEX_GRID_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests HEX_GRID_01_LAYERS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 17

  int layers;
  int nodes_per_layer;
  int nodes_per_layer_test[TEST_NUM] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101, 1001, 10001 };
  int test;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For a hexagonal grid of points in the unit square,\n";
  cout << "  given NODES_PER_LAYER, the number of grid points\n";
  cout << "  along the first layer,\n";
  cout << "\n";
  cout << "  HEX_GRID_01_LAYERS computes LAYERS, the number of\n";
  cout << "  layers.\n";
  cout << "\n";
  cout << "   NODES  LAYERS\n";
  cout << "     PER\n";
  cout << "   LAYER\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];
    layers = hex_grid_01_layers ( nodes_per_layer );
    cout << "  "
         << setw(12) << nodes_per_layer << "  "
         << setw(12) << layers << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests HEX_GRID_01_H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 17

  double hx;
  double hy;
  int layers;
  int nodes_per_layer;
  int nodes_per_layer_test[TEST_NUM] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101, 1001, 10001 };
  double temp1;
  double temp2;
  int test;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For a hexagonal grid of points in the unit square,\n";
  cout << "  given NODES_PER_LAYER, the number of grid points\n";
  cout << "  along the first layer,\n";
  cout << "\n";
  cout << "  HEX_GRID_01_H computes HX and HY, the spacings\n";
  cout << "  in the row and column directions.\n";
  cout << "\n";
  cout << "    NODES    LAYERS   HX          HY\n";
  cout << "      PER\n";
  cout << "    LAYER\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];
    layers = hex_grid_01_layers ( nodes_per_layer );
    hex_grid_01_h ( nodes_per_layer, &hx, &hy );

    cout << "  "
         << setw(6)  << nodes_per_layer << "  "
         << setw(6)  << layers          << "  "
         << setw(10) << hx              << "  "
         << setw(10) << hy              << "\n";
  }

  cout << "\n";
  cout << "  LAYERS is chosen so that LAYERS-1 layers just fit\n";
  cout << "  inside the unit square, but LAYERS layers do not\n";
  cout << "\n";
  cout << "  LAYERS      HY     (LAYERS-1)*HY    LAYERS*HY\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];
    layers = hex_grid_01_layers ( nodes_per_layer );
    hex_grid_01_h ( nodes_per_layer, &hx, &hy );

    temp1 = ( double ) ( layers - 1 ) * hy;
    temp2 = ( double ) ( layers ) * hy;

    cout << "  "
         << setw(6)  << layers << "  "
         << setw(10) << hy     << "  "
         << setw(10) << temp1  << "  "
         << setw(10) << temp2  << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests HEX_GRID_01_N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 17

  int n;
  int layers;
  int nodes_per_layer;
  int nodes_per_layer_test[TEST_NUM] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101, 1001, 10001 };
  int test;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For a hexagonal grid of points in the unit square,\n";
  cout << "  given NODES_PER_LAYER, the number of grid points\n";
  cout << "  along the first layer,\n";
  cout << "\n";
  cout << "  HEX_GRID_01_N computes N, the total number of grid\n";
  cout << "  points in the unit square.\n";
  cout << "\n";
  cout << "    NODES   LAYERS           N\n";
  cout << "      PER\n";
  cout << "    LAYER\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];

    layers = hex_grid_01_layers ( nodes_per_layer );

    n = hex_grid_01_n ( nodes_per_layer );

    cout << "  "
         << setw(6)  << nodes_per_layer << "  "
         << setw(6)  << layers          << "  "
         << setw(11) << n               << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests HEX_GRID_01_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 15

  double box[2*2];
  char file_name[81];
  double hx;
  double hy;
  int n;
  int layers;
  int nodes_per_layer;
  int nodes_per_layer_test[TEST_NUM] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101 };
  double *p;
  char *string;
  int test;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For a hexagonal grid of points in the unit square,\n";
  cout << "  given NODES_PER_LAYER, the number of grid points\n";
  cout << "  along the first layer,\n";
  cout << "\n";
  cout << "  HEX_GRID_01_POINTS computes P, the coordinates\n";
  cout << "  of the points of the grid.\n";
  cout << "\n";
  cout << "  HEX_GRID_WRITE writes the data to a file.\n";

  box[0+0*2] = 0.0;
  box[1+0*2] = 0.0;
  box[0+1*2] = 1.0;
  box[1+1*2] = 1.0;

  cout << "\n";
  cout << "   NODES  LAYERS    N    Filename\n";
  cout << "     PER\n";
  cout << "   LAYER\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];
    layers = hex_grid_01_layers ( nodes_per_layer );
    hex_grid_01_h ( nodes_per_layer, &hx, &hy );
    n = hex_grid_01_n ( nodes_per_layer);

    p = hex_grid_01_points ( nodes_per_layer, layers, n );

    strcpy ( file_name, "hex_grid_01_" );
    string = i4_to_s ( nodes_per_layer );
    strcat ( file_name, string );
    delete [] string;
    strcat ( file_name, "_" );
    string = i4_to_s ( layers );
    strcat ( file_name, string );
    delete [] string;
    strcat ( file_name, "_" );
    string = i4_to_s ( n );
    strcat ( file_name, string );
    delete [] string;
    strcat ( file_name, "_data.txt" );

    cout << "  "
         << setw(6)  << nodes_per_layer << "  "
         << setw(6)  << layers          << "  "
         << setw(12) << n               << "  "
                     << file_name       << "\n";

    hex_grid_write ( n, nodes_per_layer, layers, hx, hy, box,
      p, file_name );

    delete [] p;

  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests HEX_GRID_01_APPROXIMATE_N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  int n;
  int n_goal;
  int nodes_per_layer;
  int n_goal_test[TEST_NUM] = { 100, 200, 500, 10000 };
  int test;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For a hexagonal grid of points in the unit box,\n";
  cout << "  HEX_GRID_01_APPROXIMATE_N seeks the value of\n";
  cout << "  NODES_PER_LAYER that produces a mesh of about N points.\n";
  cout << "\n";

  cout << "\n";
  cout << "  N_GOAL  NODES_PER_LAYER       N\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n_goal = n_goal_test[test];

    hex_grid_01_approximate_n ( n_goal, &nodes_per_layer, &n );

    cout << "  "
         << setw(6)  << n_goal          << "           "
         << setw(6)  << nodes_per_layer << "  "
         << setw(6)  << n               << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests HEX_GRID_01_APPROXIMATE_H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double h;
  double h_goal;
  int n;
  int nodes_per_layer;
  double h_goal_test[TEST_NUM] = { 0.10, 0.01, 0.0001 };
  int test;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  For a hexagonal grid of points in the unit box,\n";
  cout << "  HEX_GRID_01_APPROXIMATE_H seeks the value of\n";
  cout << "  NODES_PER_LAYER that produces a mesh with spacing\n";
  cout << "  that is H or less.\n";
  cout << "\n";

  cout << "\n";
  cout << "      H_GOAL      NODES_PER_LAYER      H                      N\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    h_goal = h_goal_test[test];

    hex_grid_01_approximate_h ( h_goal, &nodes_per_layer, &h );

    n = hex_grid_01_n ( nodes_per_layer );

    cout << "  "
         << setw(14) << h_goal          << "           "
         << setw(6)  << nodes_per_layer << "  "
         << setw(14) << h               << "  "
         << setw(12) << n               << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests HEX_GRID_LAYERS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 15

  double box[2*2];
  int layers;
  int nodes_per_layer;
  int nodes_per_layer_test[TEST_NUM] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101 };
  int test;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  For a hexagonal grid of points in a coordinate box,\n";
  cout << "  given NODES_PER_LAYER, the number of grid points\n";
  cout << "  along the first layer,\n";
  cout << "\n";
  cout << "  HEX_GRID_LAYERS computes LAYERS, the number of\n";
  cout << "  layers.\n";

  box[0+0*2] = 1.0;
  box[1+0*2] = 2.0;
  box[0+1*2] = 4.0;
  box[1+1*2] = 7.0;

  box_print_2d ( box );

  cout << "\n";
  cout << "   NODES  LAYERS\n";
  cout << "     PER\n";
  cout << "   LAYER\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];
    layers = hex_grid_layers ( nodes_per_layer, box );

    cout << "  "
         << setw(6) << nodes_per_layer << "  "
         << setw(6) << layers          << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests HEX_GRID_H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 15

  double box[2*2];
  double hx;
  double hy;
  int layers;
  int nodes_per_layer;
  int nodes_per_layer_test[TEST_NUM] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101 };
  double temp1;
  double temp2;
  int test;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  For a hexagonal grid of points in a coordinate box,\n";
  cout << "  given NODES_PER_LAYER, the number of grid points\n";
  cout << "  along the first layer,\n";
  cout << "\n";
  cout << "  HEX_GRID_H computes HX and HY, the spacings\n";
  cout << "  in the row and column directions.\n";

  box[0+0*2] = 1.0;
  box[1+0*2] = 2.0;
  box[0+1*2] = 4.0;
  box[1+1*2] = 7.0;

  box_print_2d ( box );

  cout << "\n";
  cout << "    NODES    LAYERS   HX          HY\n";
  cout << "      PER\n";
  cout << "    LAYER\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];
    layers = hex_grid_layers ( nodes_per_layer, box );
    hex_grid_h ( nodes_per_layer, box, &hx, &hy );

    cout << "  "
         << setw(6)  << nodes_per_layer << "  "
         << setw(6)  << layers          << "  "
         << setw(10) << hx              << "  "
         << setw(10) << hy              << "\n";
  }

  cout << "\n";
  cout << "  LAYERS is chosen so that LAYERS-1 layers just fit\n";
  cout << "  inside the unit square, but LAYERS layers do not.\n";
  cout << "\n";
  cout << "  LAYERS      HY     (LAYERS-1)*HY    LAYERS*HY\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];
    layers = hex_grid_layers ( nodes_per_layer, box );
    hex_grid_h ( nodes_per_layer, box, &hx, &hy );

    temp1 = ( double ) ( layers - 1 ) * hy;
    temp2 = ( double ) ( layers ) * hy;

    cout << "  "
         << setw(6)  << layers          << "  "
         << setw(10) << hy              << "  "
         << setw(10) << temp1           << "  "
         << setw(10) << temp2           << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests HEX_GRID_N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 15

  double box[2*2];
  int n;
  int layers;
  int nodes_per_layer;
  int nodes_per_layer_test[TEST_NUM] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 41, 81, 101 };
  int test;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  For a hexagonal grid of points in a coordinate box,\n";
  cout << "  given NODES_PER_LAYER, the number of grid points\n";
  cout << "  along the first layer,\n";
  cout << "\n";
  cout << "  HEX_GRID_N computes N, the total number of grid\n";
  cout << "  points in the coordinate box.\n";

  box[0+0*2] = 1.0;
  box[1+0*2] = 2.0;
  box[0+1*2] = 4.0;
  box[1+1*2] = 7.0;

  box_print_2d ( box );

  cout << "\n";
  cout << "    NODES   LAYERS           N\n";
  cout << "      PER\n";
  cout << "    LAYER\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];

    layers = hex_grid_layers ( nodes_per_layer, box );

    n = hex_grid_n ( nodes_per_layer, box );

    cout << "  "
         << setw(6)  << nodes_per_layer << "  "
         << setw(6)  << layers          << "  "
         << setw(12) << n               << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests HEX_GRID_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 12

  double box[2*2];
  char file_name[81];
  double hx;
  double hy;
  int n;
  int layers;
  int nodes_per_layer;
  int nodes_per_layer_test[TEST_NUM] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21 };
  double *p;
  char *string;
  int test;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  For a hexagonal grid of points in a coordinate box,\n";
  cout << "  given NODES_PER_LAYER, the number of grid points\n";
  cout << "  along the first layer,\n";
  cout << "\n";
  cout << "  HEX_GRID_POINTS computes P, the coordinates\n";
  cout << "  of the points of the grid.\n";
  cout << "\n";
  cout << "  HEX_GRID_WRITE writes the data to a file.\n";

  box[0+0*2] = 1.0;
  box[1+0*2] = 2.0;
  box[0+1*2] = 4.0;
  box[1+1*2] = 7.0;

  box_print_2d ( box );

  cout << "\n";
  cout << "   NODES  LAYERS    N    Filename\n";
  cout << "     PER\n";
  cout << "   LAYER\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    nodes_per_layer = nodes_per_layer_test[test];

    layers = hex_grid_layers ( nodes_per_layer, box );

    hex_grid_h ( nodes_per_layer, box, &hx, &hy );

    n = hex_grid_n ( nodes_per_layer, box );

    p = hex_grid_points ( nodes_per_layer, layers, n, box );

    strcpy ( file_name, "hex_grid_" );
    string = i4_to_s ( nodes_per_layer );
    strcat ( file_name, string );
    delete [] string;
    strcat ( file_name, "_" );
    string = i4_to_s ( layers );
    strcat ( file_name, string );
    delete [] string;
    strcat ( file_name, "_" );
    string = i4_to_s ( n );
    strcat ( file_name, string );
    delete [] string;
    strcat ( file_name, "_data.txt" );

    cout << "  "
         << setw(6)  << nodes_per_layer << "  "
         << setw(3)  << layers          << "  "
         << setw(6 ) << n               << "  "
                     << file_name       << "\n";

    hex_grid_write ( n, nodes_per_layer, layers, hx, hy, box,
      p, file_name );

    delete [] p;
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests HEX_GRID_APPROXIMATE_N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double box[2*2];
  int n;
  int n_goal;
  int nodes_per_layer;
  int n_goal_test[TEST_NUM] = { 100, 200, 500, 10000 };
  int test;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  For a hexagonal grid of points in a coordinate box,\n";
  cout << "  HEX_GRID_APPROXIMATE_N seeks the value of\n";
  cout << "  NODES_PER_LAYER that produces a mesh of about N points.\n";

  box[0+0*2] = 1.0;
  box[1+0*2] = 2.0;
  box[0+1*2] = 4.0;
  box[1+1*2] = 7.0;

  box_print_2d ( box );

  cout << "\n";
  cout << "  N_GOAL  NODES_PER_LAYER       N\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n_goal = n_goal_test[test];

    hex_grid_approximate_n ( box, n_goal, &nodes_per_layer, &n );

    cout << "  "
         << setw(6)  << n_goal          << "           "
         << setw(6)  << nodes_per_layer << "  "
         << setw(6)  << n               << "\n";
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests HEX_GRID_APPROXIMATE_H.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double box[2*2];
  double h;
  double h_goal;
  int n;
  int nodes_per_layer;
  double h_goal_test[TEST_NUM] = { 0.10, 0.01, 0.0001 };
  int test;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  For a hexagonal grid of points in the unit box,\n";
  cout << "  HEX_GRID_APPROXIMATE_H seeks the value of\n";
  cout << "  NODES_PER_LAYER that produces a mesh with spacing\n";
  cout << "  that is H or less.\n";
  cout << "\n";

  box[0+0*2] = 1.0;
  box[1+0*2] = 2.0;
  box[0+1*2] = 4.0;
  box[1+1*2] = 7.0;

  box_print_2d ( box );

  cout << "\n";
  cout << "      H_GOAL      NODES_PER_LAYER      H                      N\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    h_goal = h_goal_test[test];

    hex_grid_approximate_h ( box, h_goal, &nodes_per_layer, &h );

    n = hex_grid_n ( nodes_per_layer, box );

    cout << "  "
         << setw(14) << h_goal          << "           "
         << setw(6)  << nodes_per_layer << "  "
         << setw(14) << h               << "  "
         << setw(12) << n               << "\n";
  }

  return;
# undef TEST_NUM
}
