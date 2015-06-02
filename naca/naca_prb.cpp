# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "naca.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for NACA_PRB.
//
//  Discussion:
//
//    NACA_PRB tests the NACA library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "NACA_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the NACA library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "NACA_PRB:\n";
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
//    TEST01 tests NACA4_SYMMETRIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c;
  string command_filename = "symmetric_commands.txt";
  ofstream command_unit;
  string data_filename = "symmetric_data.txt";
  int i;
  int n = 51;
  double ratio;
  double t;
  double *x;
  double x_max;
  double x_min;
  double *xy;
  double *y;
  double y_max;
  double y_min;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  NACA4_SYMMETRIC evaluates y(x) for a NACA\n";
  cout << "  symmetric airfoil defined by a 4-digit code.\n";

  c = 10.0;
  t = 0.15;
  x = r8vec_linspace_new ( n, 0.0, c );
  y = naca4_symmetric ( t, c, n, x );
//
//  Reorganize data into a single object.
//
  xy = new double[2*2*n];

  for ( i = 0; i < n; i++ )
  {
    xy[0+i*2] = x[i];
    xy[1+i*2] = -y[i];
  }
  for ( i = 0; i < n; i++ )
  {
    xy[0+(n+i)*2] = x[n-1-i];
    xy[1+(n+i)*2] = y[n-1-i];
  }
//
//  Determine size ratio.
//
  x_min = r8vec_min ( n, x );
  x_max = r8vec_max ( n, x );
  y_max = r8vec_max ( n, y );
  y_min = - y_max;
  ratio = ( y_max - y_min ) / ( x_max - x_min );
//
//  Save data to a file.
//
  r8mat_write ( data_filename, 2, 2 * n, xy );
  cout << "  Data saved in file '" << data_filename << "'\n";
//
//  Create the command file.
//
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "set term png\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "set size ratio " << ratio << "\n";
  command_unit << "set timestamp\n";
  command_unit << "unset key\n";
  command_unit << "set output 'symmetric.png'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set title 'NACA Symmetric Airfoil'\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 with lines lw 3\n";
  command_unit << "quit\n";

  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests NACA4_CAMBERED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double c;
  string command_filename = "cambered_commands.txt";
  ofstream command_unit;
  string data_filename = "cambered_data.txt";
  int i;
  double m;
  int n = 51;
  double p;
  double ratio;
  double t;
  double x_max;
  double x_min;
  double *xc;
  double *xl;
  double *xu;
  double *xy;
  double y_max;
  double y_min;
  double *yl;
  double *yu;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  NACA4_CAMBERED evaluates (xu,yu) and (xl,yl) for a NACA\n";
  cout << "  cambered airfoil defined by a 4-digit code.\n";

  m = 0.02;
  p = 0.4;
  t = 0.12;
  c = 10.0;

  xc = r8vec_linspace_new ( n, 0.0, c );

  xu = new double[n];
  xl = new double[n];
  yu = new double[n];
  yl = new double[n];

  naca4_cambered ( m, p, t, c, n, xc, xu, yu, xl, yl );
//
//  Reorganize data into a single object.
//
  xy = new double[2*2*n];

  for ( i = 0; i < n; i++ )
  {
    xy[0+i*2] = xl[i];
    xy[1+i*2] = yl[i];
  }
  for ( i = 0; i < n; i++ )
  {
    xy[0+(n+i)*2] = xu[n-1-i];
    xy[1+(n+i)*2] = yu[n-1-i];
  }
//
//  Determine size ratio.
//
  x_min = r8_min ( r8vec_min ( n, xl ), r8vec_min ( n, xu ) );
  x_max = r8_max ( r8vec_max ( n, xl ), r8vec_max ( n, xu ) );
  y_min = r8_min ( r8vec_min ( n, yl ), r8vec_min ( n, yu ) );
  y_max = r8_max ( r8vec_max ( n, yl ), r8vec_max ( n, yu ) );
  ratio = ( y_max - y_min ) / ( x_max - x_min );
//
//  Save data to a file.
//
  r8mat_write ( data_filename, 2, 2 * n, xy );
  cout << "  Data saved in file '" << data_filename << "'\n";
//
//  Create the command file.
//
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "set term png\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "set size ratio " << ratio << "\n";
  command_unit << "set timestamp\n";
  command_unit << "unset key\n";
  command_unit << "set output 'cambered.png'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set title 'NACA Cambered Airfoil'\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 with lines lw 3\n";
  command_unit << "quit\n";

  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  delete [] xc;
  delete [] xl;
  delete [] xu;
  delete [] xy;
  delete [] yl;
  delete [] yu;

  return;
}

