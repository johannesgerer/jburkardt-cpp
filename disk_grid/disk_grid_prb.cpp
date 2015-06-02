# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>

using namespace std;

# include "disk_grid.hpp"

int main ( );
void disk_grid_test01 ( );
void disk_grid_test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DISK_GRID_PRB.
//
//  Discussion:
//
//    DISK_GRID_PRB tests the DISK_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "DISK_GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the DISK_GRID library.\n";

  disk_grid_test01 ( );
  disk_grid_test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISK_GRID_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void disk_grid_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    DISK_GRID_TEST01 tests DISK_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  string boundary_filename = "disk_grid_test01_boundary.txt";
  ofstream boundary_unit;
  double c[2];
  double *cg;
  string command_filename = "disk_grid_test01_commands.txt";
  ofstream command_unit;
  string data_filename = "disk_grid_test01_data.txt";
  ofstream data_unit;
  string filename = "disk_grid_test01.xy";
  int i;
  int n;
  int ng;
  const double pi = 3.141592653589793;
  double r;
  double t;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  DISK_GRID can define a grid of points\n";
  cout << "  with N+1 points on a horizontal or vertical radius,\n";
  cout << "  based on any disk.\n";

  n = 20;
  r = 2.0;
  c[0] = 1.0;
  c[1] = 5.0;

  cout << "\n";
  cout << "  We use N = " << n << "\n";
  cout << "  Radius R = " << r << "\n";;
  cout << "  Center C = (" << c[0] << "," << c[1] << ")\n";

  ng = disk_grid_count ( n, r, c );

  cout << "\n";
  cout << "  Number of grid points will be " << ng << "\n";

  cg = disk_grid ( n, r, c, ng );

  r82vec_print_part ( ng, cg, 20, "  Part of the grid point array:" );
//
//  Write the coordinate data to a file.
//
  r8mat_write ( filename, 2, ng, cg );

  cout << "\n";
  cout << "  Data written to the file \"" << filename << "\"\n";
//
//  Create graphics data files.
//
  boundary_unit.open ( boundary_filename.c_str ( ) );
  for ( i = 0; i <= 50; i++ )
  {
    t = 2.0 * pi * ( double ) ( i ) / 50.0;
    boundary_unit << "  " << c[0] + r * cos ( t )
                  << "  " << c[1] + r * sin ( t ) << "\n";
  }
  boundary_unit.close ( );
  cout << "\n";
  cout << "  Created boundary file \"" << boundary_filename << "\".\n";

  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < ng; i++ )
  {
    data_unit << "  " << cg[0+i*2]
              << "  " << cg[1+i*2] << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created data file \"" << data_filename << "\"\n";
//
//  Create graphics command file.
//
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'disk_grid_test01.png'\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set title 'Disk Grid'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "set size ratio -1\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 with points lt 3 pt 3,\\\n";
  command_unit << "    '" << boundary_filename 
               << "' using 1:2 lw 3 linecolor rgb 'black'\n";
  command_unit << "quit\n";
  command_unit.close ( );

  cout << "  Created command file \"" << command_filename << "\"\n";

  delete [] cg;

  return;
}
//****************************************************************************80

void disk_grid_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    DISK_GRID_TEST02 tests DISK_GRID_FIBONACCI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  string boundary_filename = "disk_grid_test02_boundary.txt";
  ofstream boundary_unit;
  double c[2];
  string command_filename = "disk_grid_test02_commands.txt";
  ofstream command_unit;
  string data_filename = "disk_grid_test02_data.txt";
  ofstream data_unit;
  string filename = "disk_grid_test02.xy";
  double *g;
  int i;
  int n;
  const double pi = 3.141592653589793;
  double r;
  double t;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  DISK_GRID_FIBONACCI can define a grid of N points\n";
  cout << "  based on a Fibonacci spiral inside a disk.\n";

  n = 1000;
  r = 2.0;
  c[0] = 1.0;
  c[1] = 5.0;

  cout << "\n";
  cout << "  We use N = " << n << "\n";
  cout << "  Radius R = " << r << "\n";;
  cout << "  Center C = (" << c[0] << "," << c[1] << ")\n";

  g = disk_grid_fibonacci ( n, r, c );

  r82vec_print_part ( n, g, 20, "  Part of the grid point array:" );
//
//  Write the coordinate data to a file.
//
  r8mat_write ( filename, 2, n, g );

  cout << "\n";
  cout << "  Data written to the file \"" << filename << "\"\n";
//
//  Create graphics data files.
//
  boundary_unit.open ( boundary_filename.c_str ( ) );
  for ( i = 0; i <= 50; i++ )
  {
    t = 2.0 * pi * ( double ) ( i ) / 50.0;
    boundary_unit << "  " << c[0] + r * cos ( t )
                  << "  " << c[1] + r * sin ( t ) << "\n";
  }
  boundary_unit.close ( );
  cout << "\n";
  cout << "  Created boundary file \"" << boundary_filename << "\".\n";

  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    data_unit << "  " << g[0+i*2]
              << "  " << g[1+i*2] << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created data file \"" << data_filename << "\"\n";
//
//  Create graphics command file.
//
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'disk_grid_test02.png'\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set title 'Fibonacci Disk Grid'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "set size ratio -1\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 with points lt 3 pt 3,\\\n";
  command_unit << "    '" << boundary_filename 
               << "' using 1:2 lw 3 linecolor rgb 'black'\n";
  command_unit << "quit\n";
  command_unit.close ( );

  cout << "  Created command file \"" << command_filename << "\"\n";

  delete [] g;

  return;
}
