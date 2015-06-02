# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "wedge_grid.hpp"

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

double *wedge_grid ( int n, int ng )  

//****************************************************************************80
//
//  Purpose:
//
//    WEDGE_GRID computes grid points in the unit wedge in 3D.
//
//  Discussion:
//
//    The interior of the unit wedge in 3D is defined by the constraints:
//      0 <= X
//      0 <= Y
//           X + Y <= 1
//     -1 <= Z <= +1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//    0 <= N.
//
//    Input, int NG, the number of grid points.
//    This can be computed by WEDGE_GRID_SIZE, or else determined by
//    NG =(N+1)*((N+1)*(N+2))/2.
//
//    Output, double WEDGE+GRID[3*NG], the coordinates of the grid points.
//
{
  int i;
  double ir;
  int j;
  double jr;
  int k;
  double kr;
  double nr;
  int p;
  double *g;

  g = new double[3*ng];

  if ( n == 0 )
  {
    g[0+0*3] = 0.5;
    g[1+0*3] = 0.5;
    g[2+0*3] = 0.0;
    return g;
  }

  p = 0;
  nr = ( double ) ( n );

  for ( k = 0; k <= n; k++ )
  {
    kr = ( double ) ( 2 * k - n ) / nr;
    for ( j = 0; j <= n; j++ )
    {
      jr = ( double ) ( j ) / nr;
      for ( i = 0; i <= n - j; i++ )
      {
        ir = ( double ) ( i ) / nr;
        g[0+p*3] = ir;
        g[1+p*3] = jr;
        g[2+p*3] = kr;
        p = p + 1;
      }
    }
  }

  return g;
}
//****************************************************************************80

int wedge_grid_size ( int n )  

//****************************************************************************80
//
//  Purpose:
//
//    WEDGE_GRID_SIZE counts the points in a grid of the unit wedge in 3D.
//
//  Discussion:
//
//    The interior of the unit wedge in 3D is defined by the constraints:
//      0 <= X
//      0 <= Y
//           X + Y <= 1
//     -1 <= Z <= +1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//    0 <= N.
//
//    Output, int WEDGE_GRID_SIZE, the number of grid points.
//
{
  int ng;

  ng = ( n + 1 ) * ( ( n + 1 ) * ( n + 2 ) ) / 2;

  return ng;
}
//****************************************************************************80

void wedge_grid_plot ( int n, int ng, double g[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    WEDGE_GRID_PLOT sets up a GNUPLOT plot of a unit wedge grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, int NG, the number of nodes.
//
//    Input, double G[3*NG], the grid point coordinates.
//
//    Input, string HEADER, the header for the files.
//
{
  string command_filename;
  ofstream command_unit;
  int i;
  int j;
  int l;
  string node_filename;
  ofstream node_unit;
  string plot_filename;
  double v1[3];
  double v2[3];
  double v3[3];
  double v4[3];
  double v5[3];
  double v6[3];
  string vertex_filename;
  ofstream vertex_unit;
//
//  Create the vertex file.
//
  wedge_vertices ( v1, v2, v3, v4, v5, v6 );

  vertex_filename = header + "_vertices.txt";
  vertex_unit.open ( vertex_filename.c_str ( ) );

  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << v2[0] << "  "
              << v2[1] << "  "
              << v2[2] << "\n";
  vertex_unit << v3[0] << "  "
              << v3[1] << "  "
              << v3[2] << "\n";
  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << "\n";

  vertex_unit << v4[0] << "  "
              << v4[1] << "  "
              << v4[2] << "\n";
  vertex_unit << v5[0] << "  "
              << v5[1] << "  "
              << v5[2] << "\n";
  vertex_unit << v6[0] << "  "
              << v6[1] << "  "
              << v6[2] << "\n";
  vertex_unit << v4[0] << "  "
              << v4[1] << "  "
              << v4[2] << "\n";
  vertex_unit << "\n";

  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << v4[0] << "  "
              << v4[1] << "  "
              << v4[2] << "\n";
  vertex_unit << "\n";

  vertex_unit << v2[0] << "  "
              << v2[1] << "  "
              << v2[2] << "\n";
  vertex_unit << v5[0] << "  "
              << v5[1] << "  "
              << v5[2] << "\n";
  vertex_unit << "\n";

  vertex_unit << v3[0] << "  "
              << v3[1] << "  "
              << v3[2] << "\n";
  vertex_unit << v6[0] << "  "
              << v6[1] << "  "
              << v6[2] << "\n";
  vertex_unit << "\n";

  vertex_unit.close ( );
  cout << "\n";
  cout << "  Created vertex file '" << vertex_filename << "'\n";
//
//  Create the node file.
//
  node_filename = header + "_nodes.txt";

  node_unit.open ( node_filename.c_str ( ) );
  for ( j = 0; j < ng; j++ )
  {
    node_unit << g[0+j*3] << "  "
              << g[1+j*3] << "  "
              << g[2+j*3] << "\n";
  }
  node_unit.close ( );
  cout << "  Created node file '" << node_filename << "'\n";
//
//  Create the command file.
//
  command_filename = header + "_commands.txt";

  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";

  plot_filename = header + ".png";

  command_unit << "set output '" << plot_filename << "\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << header << "'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "#set view equal xyz\n";
  command_unit << "set view 80, 85\n";
  command_unit << "set style data lines\n";
  command_unit << "set timestamp\n";
  command_unit << "splot '" << vertex_filename << "' with lines lt 3, \\\n";
  command_unit << "      '" << node_filename << "' with points pt 7 lt 0\n";
  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  return;
}
//****************************************************************************80

void wedge_vertices ( double v1[], double v2[], double v3[], double v4[], 
  double v5[], double v6[] )

//****************************************************************************80
//
//  Purpose:
//
//    WEDGE_VERTICES returns the vertices of the unit wege.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double V1[3], V2[3], V3[3], V4[3], V5[3], V6[3],
//    the vertices.
//
{
  static double v1_save[3] = {  0.0,  0.0, -1.0 };
  static double v2_save[3] = {  1.0,  0.0, -1.0 };
  static double v3_save[3] = {  0.0,  1.0, -1.0 };
  static double v4_save[3] = {  0.0,  0.0, +1.0 };
  static double v5_save[3] = {  1.0,  0.0, +1.0 };
  static double v6_save[3] = {  0.0,  1.0, +1.0 };

  r8vec_copy ( 3, v1_save, v1 );
  r8vec_copy ( 3, v2_save, v2 );
  r8vec_copy ( 3, v3_save, v3 );
  r8vec_copy ( 3, v4_save, v4 );
  r8vec_copy ( 3, v5_save, v5 );
  r8vec_copy ( 3, v6_save, v6 );

  return;
}
