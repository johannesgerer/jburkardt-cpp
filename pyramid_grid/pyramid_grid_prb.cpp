# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "pyramid_grid.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PYRAMID_GRID_PRB.
//
//  Discussion:
//
//    PYRAMID_GRID_PRB tests the PYRAMID_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "PYRAMID_GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PYRAMID_GRID library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PYRAMID_GRID_PRB:\n";
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
//    TEST01 tests PYRAMID_GRID_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int ng;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  PYRAMID_GRID_SIZE determines the size of a\n";
  cout << "  pyramid grid with N+1 points along each edge.\n";

  cout << "\n";
  cout << "   N    Size\n";
  cout << "\n";
  for ( n = 0; n <= 10; n++ )
  {
    ng = pyramid_grid_size ( n );
    cout << setw(4) << n << "  "
         << setw(6) << ng << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PYRAMID_UNIT_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int ng;
  double *pg;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PYRAMID_UNIT_GRID determines a unit pyramid\n";
  cout << "  grid with N+1 points along each edge.\n";

  n = 4;
  r8_print ( n, "  Grid parameter N:" );

  ng = pyramid_grid_size ( n );
  r8_print ( ng, "  Grid size NG:" );

  pg = pyramid_unit_grid ( n, ng );

  r8mat_transpose_print ( 3, ng, pg, "  Pyramid grid points:" );

  delete [] pg;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests PYRAMID_UNIT_GRID_PLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  string header;
  int n;
  int ng;
  double *pg;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  PYRAMID_UNIT_GRID_PLOT plots a unit pyramid\n";
  cout << "  grid with N+1 points along each edge.\n";

  n = 5;
  r8_print ( n, "  Grid parameter N:" );

  ng = pyramid_grid_size ( n );
  r8_print ( ng, "  Grid size NG:" );

  pg = pyramid_unit_grid ( n, ng );

  header = "pyramid_unit";
  pyramid_unit_grid_plot ( n, ng, pg, header );

  delete [] pg;

  return;
}
