# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "tetrahedron_grid.hpp"

int main ( );
void tetrahedron_grid_test01 ( );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TETRAHEDRON_GRID_PRB.
//
//  Discussion:
//
//    TETRAHEDRON_GRID_PRB tests the TETRAHEDRON_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TETRAHEDRON_GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the TETRAHEDRON_GRID library.\n";

  tetrahedron_grid_test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TETRAHEDRON_GRID_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void tetrahedron_grid_test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_GRID_TEST01 tests TETRAHEDRON_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "tetrahedron_grid_test01.xyz";
  int n;
  int ng;
  double t[3*4] = {
    0.0, 0.0, 0.0, 
    1.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 
    0.0, 0.0, 1.0 };
  double *tg;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  TETRAHEDRON_GRID can define a tetrahedral grid\n";
  cout << "  with N+1 points on a side, based on any tetrahedron.\n";

  n = 10;
  cout << "  N = " << n << "\n";

  ng = tetrahedron_grid_count ( n );

  r8mat_print ( 3, 4, t, "  Tetrahedron vertices:" );

  tg = tetrahedron_grid ( n, t, ng );

  r83vec_print_part ( ng, tg, 20, "  Part of the grid point array:" );

  r8mat_write ( filename, 3, ng, tg );

  cout << "\n";
  cout << "  Data written to the file \"" << filename << "\".\n";

  delete [] tg;

  return;
}
