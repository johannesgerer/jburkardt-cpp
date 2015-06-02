# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "brownian_motion_simulation.hpp"

int main ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BROWNIAN_MOTION_SIMULATION_PRB.
//
//  Discussion:
//
//    BROWNIAN_MOTION_SIMULATION_PRB tests the BROWNIAN_MOTION_SIMULATION library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  double d;
  double *dsq;
  string header;
  int k;
  int m;
  int n;
  int seed;
  double t;
  double *x;

  timestamp ( );
  cout << "\n";
  cout << "BROWNIAN_MOTION_SIMULATION_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BROWNIAN_MOTION_SIMULATION library.\n";
//
//  Compute the path of a particle undergoing Brownian motion.
//
  for ( m = 1; m <= 2; m++ )
  {
    n = 1001;
    d = 10.0;
    t = 1.0;
    seed = 123456789;
    x = brownian_motion_simulation ( m, n, d, t, seed );
    if ( m == 1 )
    {
      header = "motion_1d";
    }
    else if ( m == 2 )
    {
      header = "motion_2d";
    }
    brownian_motion_display ( m, n, x, header );
    delete [] x;
  }
//
//  Estimate the average displacement of the particle from the origin
//  as a function of time.
//
  for ( m = 1; m <= 3; m++ )
  {
    k = 40;
    n = 1001;
    d = 10.0;
    t = 1.0;
    seed = 123456789;

    dsq = brownian_displacement_simulation ( k, n, m, d, t, seed );
    if ( m == 1 )
    {
      header = "displacement_1d";
    }
    else if ( m == 2 )
    {
      header = "displacement_2d";
    }
    else if ( m == 3 )
    {
      header = "displacement_3d";
    }
    brownian_displacement_display ( k, n, d, t, dsq, header );
    delete [] dsq;
  }
/*
  Terminate.
*/
  cout << "\n";
  cout << "BROWNIAN_MOTION_SIMULATION_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
