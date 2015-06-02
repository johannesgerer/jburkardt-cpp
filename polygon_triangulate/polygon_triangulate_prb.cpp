# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "polygon_triangulate.hpp"

int main ( );
void test01 ( );
void test02 ( string prefix );
void test03 ( string prefix );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POLYGON_TRIANGULATE_PRB.
//
//  Discussion:
//
//    POLYGON_TRIANGULATE_PRB tests the POLYGON_TRIANGULATE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "POLYGON_TRIANGULATE_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the POLYGON_TRIANGULATE library.\n";

  test01 ( );

  test02 ( "comb" );
  test02 ( "hand" );
  test02 ( "i18" );

  test03 ( "comb" );
  test03 ( "hand" );
  test03 ( "i18" );
//
//  Terminate.
//
  cout << "\n";
  cout << "POLYGON_TRIANGULATE_PRB\n";
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
//    TEST01 tests the "comb_10" polygon.
//
//  Discussion:
//
//    There are N-3 triangles in the triangulation.
//
//    For the first N-2 triangles, the first edge listed is always an
//    internal diagonal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 May 2014
//
{
  int n = 10;
  int *triangles;
  double x[10] = {
    8.0, 7.0, 6.0, 5.0, 4.0, 
    3.0, 2.0, 1.0, 0.0, 4.0 };
  double y[10] = {
    0.0, 10.0,  0.0, 10.0,  0.0, 
   10.0,  0.0, 10.0,  0.0, -2.0 };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Triangulate the comb_10 polygon.\n";

  triangles = polygon_triangulate ( n, x, y );

  i4mat_print ( 3, n - 2, triangles, "  Triangles" );

  delete [] triangles;

  return;
}
//****************************************************************************80

void test02 ( string prefix )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 triangulates a polygon described in a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 May 2014
//
//  Author:
//
//    John Burkardt.
//
{
  int dim_num;
  string element_filename;
  int i;
  int n;
  string node_filename;
  int triangle_num;
  int *triangles;
  double *x;
  double *xy;
  double *y;
//
//  Create filenames.
//
  node_filename = prefix + "_nodes.txt";
  element_filename = prefix + "_elements.txt";

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Read polygon coordinates in \"" << node_filename << "\"\n";
//
//  Read the node coordinates.
//
  r8mat_header_read ( node_filename, dim_num, n );

  xy = r8mat_data_read ( node_filename, 2, n );
//
//  Get the triangulation.
//
  x = new double[n];
  y = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = xy[0+i*2];
    y[i] = xy[1+i*2];
  }
  triangles = polygon_triangulate ( n, x, y );
//
//  Write the triangulation to a file.
//
  triangle_num = n - 2;
  i4mat_write ( element_filename, 3, triangle_num, triangles );

  cout << "  Write triangulation to \"" << element_filename << "\"\n";
//
//  Free memory.
//
  delete [] triangles;
  delete [] x;
  delete [] xy;
  delete [] y;

  return;
}
//****************************************************************************80

void test03 ( string prefix )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 plots a triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 May 2014
//
//  Author:
//
//    John Burkardt.
//
{
  string command_filename;
  ofstream command_unit;
  string diagonal_filename;
  ofstream diagonal_unit;
  int dim_num;
  string edge_filename;
  ofstream edge_unit;
  int i;
  int i2;
  int j;
  int j2;
  int n;
  int node;
  string node_filename;
  string plot_filename;
  int triangle_num;
  int *triangles;
  double *x;
  double *xy;
  double *y;
//
//  Create filenames.
//
  node_filename = prefix + "_nodes.txt";
  edge_filename = prefix + "_edges.txt";
  diagonal_filename = prefix + "_diagonals.txt";
  command_filename = prefix + "_commands.txt";
  plot_filename = prefix + ".png";

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Read node coordinates in \"" << node_filename << "\"\n";
//
//  Read the node coordinates.
//
  r8mat_header_read ( node_filename, dim_num, n );

  xy = r8mat_data_read ( node_filename, 2, n );
//
//  Get the triangulation.
//
  x = new double[n];
  y = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = xy[0+i*2];
    y[i] = xy[1+i*2];
  }
  triangles = polygon_triangulate ( n, x, y );
//
//  Plot the edges.
//
  edge_unit.open ( edge_filename.c_str ( ) );

  for ( j = 0; j < n + 1; j++ )
  {
    j2 = ( j % n );
    edge_unit << xy[0+j2*2] << "  "
              << xy[1+j2*2] << "\n";
  }

  edge_unit.close ( );
//
//  Plot the diagonals.
//
  diagonal_unit.open ( diagonal_filename.c_str ( ) );

  for ( j = 0; j < n - 3; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      node = triangles[i+j*3];
      diagonal_unit << xy[0+node*2] << "  "
                    << xy[1+node*2] << "\n";
    }
    diagonal_unit << "\n";
  }

  diagonal_unit.close ( );
//
//  Write the GNUPLOT command file.
//
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output \"" << plot_filename << "\"\n";
  command_unit << "set nokey\n";
  command_unit << "set size ratio 1\n";
  command_unit << "set timestamp\n";
  command_unit << "set xlabel \"<---X--->\"\n";
  command_unit << "set ylabel \"<---Y--->\n";
  command_unit << "set title \"Edges (green) and Diagonals (red)\"\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot \"" << edge_filename 
               << "\" using 1:2 lw 3 linecolor rgb \"green\",\\\n";
  command_unit << "     \"" << diagonal_filename 
               << "\" using 1:2 lw 3 linecolor rgb \"red\",\\\n";
  command_unit << "     \"" << node_filename 
               << "\" using 1:2 with points pt 7 ps 2 lc rgb \"black\"\n";

  command_unit.close ( );

  cout << "\n";
  cout << "  Write edges to \"" << edge_filename << "\"\n";
  cout << "  Write diagonals to \"" << diagonal_filename << "\"\n";
  cout << "  Write gnuplot commands to \"" << command_filename << "\"\n";
//
//  Free memory.
//
  delete [] triangles;
  delete [] x;
  delete [] xy;
  delete [] y;

  return;
}

