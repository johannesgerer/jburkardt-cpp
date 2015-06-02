# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "triangle_grid.hpp"

int main ( );

void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_GRID_PRB.
//
//  Discussion:
//
//    TRIANGLE_GRID_PRB tests the TRIANGLE_GRID library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_GRID_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the TRIANGLE_GRID library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_GRID_PRB:\n";
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
//    TEST01 tests TRIANGLE_GRID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  int n = 10;
  int ng = ((n+1)*(n+2))/2;

  string filename;
  int j;
  ofstream output;
  double t[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.5, 0.86602540378443860 };
  double *tg;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  TRIANGLE_GRID can define a triangular grid of points\n";
  cout << "  with N+1 points on a side, based on any triangle.\n";

  cout << "\n";
  cout << "  Defining triangle:\n";
  cout << "     J      X      Y\n";
  cout << "\n";
  for ( j = 0; j < 3; j++ )
  {
    cout << "  " << setw(4) << j
         << "  " << setw(12) << t[0+j*2]
         << "  " << setw(12) << t[1+j*2] << "\n";
  }
  tg = triangle_grid ( n, t );

  cout << "\n";
  cout << "     J      X      Y\n";
  cout << "\n";
  for ( j = 0; j < ng; j++ )
  {
    cout << "  " << setw(4) << j
         << "  " << setw(12) << tg[0+j*2]
         << "  " << setw(12) << tg[1+j*2] << "\n";
  }

  filename = "triangle_grid_test01.xy";

  output.open ( filename.c_str ( ) );
  for ( j = 0; j < ng; j++ )
  {
    output << "  " << setw(12) << tg[0+j*2]
           << "  " << setw(12) << tg[1+j*2] << "\n";
  }
  output.close ( );

  cout << "\n";
  cout << "  Data written to \"" << filename << "\"\n";

  delete [] tg;

  return;
}


