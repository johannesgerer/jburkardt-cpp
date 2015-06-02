# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "curve_plot.h"

int main ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    CURVE_PLOT_PRB demonstrates how CURVE_PLOT can be used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *x;
  double *y;

  cout << "\n";
  cout << "CURVE_PLOT_PRB:\n";
  cout << "  Demonstrate how CURVE_PLOT can be used.\n";
//
//  Set up some data to plot.
//
  n = 51;
  x = new double[n];
  y = new double[n];
  for ( i = 0; i < 51; i++ )
  {
    x[i] = ( double ) ( i ) / 10.0;
    y[i] = x[i] * cos ( x[i] );
  }
//
//  Send the data to curve_plot.
//
  curve_plot ( n, x, y, "curve_plot" );
//
//  Free memory.
//
  delete [] x;
  delete [] y;

  return 0;
}
