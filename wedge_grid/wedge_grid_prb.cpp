# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "wedge_grid.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WEDGE_GRID_PRB.
//
//  Discussion:
//
//    WEDGE_GRID_PRB tests the WEDGE_GRID library.
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
{
  timestamp ( );
  cout << "\n";
  cout << "WEDGE_GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the WEDGE_GRID library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "WEDGE_GRID_PRB:\n";
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
//    TEST01 tests WEDGE_GRID.
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
{
  double *g;
  int j;
  int n = 5;
  int ng;
  string output_filename;
  ofstream output_unit;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  WEDGE_GRID can define a grid of points\n";
  cout << "  with N+1 points on a side\n";
  cout << "  over the interior of the unit wedge in 3D.\n";

  cout << "\n";
  cout << "  Grid order N = " << n << "\n";

  ng = wedge_grid_size ( n );

  cout << "  Grid count NG = " << ng << "\n";

  g = wedge_grid ( n, ng );

  cout << "\n";
  cout << "     J      X                Y               Z\n";
  cout << "\n";
  for ( j = 0; j < ng; j++ )
  {
    cout << setw(6) << j << "  "
         << setw(14) << g[0+j*3] << "  "
         << setw(14) << g[1+j*3] << "  "
         << setw(14) << g[2+j*3] << "\n";
  }

  output_filename = "wedge_grid.xy";

  output_unit.open ( output_filename.c_str ( ) );
  for ( j = 0; j < ng; j++ )
  {
    output_unit << g[0+j*3] << "  "
                << g[1+j*3] << "  "
                << g[2+j*3] << "\n";
  }
  output_unit.close ( );

  cout << "\n";
  cout << "  Data written to '" << output_filename << "'\n";

  delete [] g;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests WEDGE_GRID_PLOT.
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
{
  double *g;
  string header;
  int n = 5;
  int ng;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  WEDGE_GRID_PLOT can create GNUPLOT data files\n";
  cout << "  for displaying a wedge grid.\n";

  cout << "\n";
  cout << "  Grid order N = " << n << "\n";

  ng = wedge_grid_size ( n );

  cout << "  Grid count NG = " << ng << "\n";

  g = wedge_grid ( n, ng );

  header = "wedge";

  wedge_grid_plot ( n, ng, g, header );

  delete [] g;

  return;
}

