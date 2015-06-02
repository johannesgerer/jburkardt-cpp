# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "ode.hpp"

int main ( );
void test01 ( );
void test02 ( );
void f01 ( double t, double y[], double yp[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ODE_PRB.
//
//  Discussion:
//
//    ODE_PRB tests the ODE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ODE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the ODE library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ODE_PRB\n";
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
//    TEST01 tests ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  double abserr;
  int i;
  int iflag;
  int iwork[5];
  int neqn = 2;
  double pi = 3.141592653589793;
  double relerr;
  int step_num = 12;
  double t;
  double tout;
  double *work;
  double *y;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  ODE solves a system of ordinary differential\n";
  cout << "  equations.\n";
  cout << "\n";
  cout << "      T           Y(1)         Y(2)\n";
  cout << "\n";

  abserr = 0.00001;
  relerr = 0.00001;

  iflag = 1;

  t = 0.0;
  y = new double[neqn];
  y[0] = 1.0;
  y[1] = 0.0;

  cout << "  " << setw(8) << t
       << "  " << setw(14) << y[0]
       << "  " << setw(14) << y[1] << "\n";

  work = new double[100+21*neqn];

  for ( i = 1; i <= step_num; i++ )
  {
    tout = ( double ) ( i ) * 2.0 * pi / ( double ) ( step_num );

    ode ( f01, neqn, y, t, tout, relerr, abserr, iflag, work, iwork );

    if ( iflag != 2 )
    {
      cout << "\n";
      cout << "TEST01 - Fatal error!\n";
      cout << "  ODE returned IFLAG = " << iflag << "\n";
      break;
    }
    cout << "  " << setw(8) << t
         << "  " << setw(14) << y[0]
         << "  " << setw(14) << y[1] << "\n";
  }

  delete [] work;
  delete [] y;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests ODE by integrating in the NEGATIVE time direction.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double abserr;
  int i;
  int iflag;
  int iwork[5];
  int neqn = 2;
  double pi = 3.141592653589793;
  double relerr;
  int step_num = 12;
  double t;
  double tout;
  double *work;
  double *y;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  ODE solves a system of ordinary differential\n";
  cout << "  equations.\n";
  cout << "\n";
  cout << "  In this example, we integrate in the negative\n";
  cout << "  time direction.\n";
  cout << "\n";
  cout << "      T           Y(1)         Y(2)\n";
  cout << "\n";

  abserr = 0.00001;
  relerr = 0.00001;

  iflag = 1;

  t = 0.0;
  y = new double[neqn];
  y[0] = 1.0;
  y[1] = 0.0;

  cout << "  " << setw(8) << t
       << "  " << setw(14) << y[0]
       << "  " << setw(14) << y[1] << "\n";

  work = new double[100+21*neqn];

  for ( i = 1; i <= step_num; i++ )
  {
    tout = - ( double ) ( i ) * 2.0 * pi / ( double ) ( step_num );

    ode ( f01, neqn, y, t, tout, relerr, abserr, iflag, work, iwork );

    if ( iflag != 2 )
    {
      cout << "\n";
      cout << "TEST02 - Fatal error!\n";
      cout << "  ODE returned IFLAG = " << iflag << "\n";
      break;
    }
    cout << "  " << setw(8) << t
         << "  " << setw(14) << y[0]
         << "  " << setw(14) << y[1] << "\n";
  }

  delete [] work;
  delete [] y;

  return;
}
//****************************************************************************80

void f01 ( double t, double y[], double yp[] )

//****************************************************************************80
//
//  Purpose:
//
//    F01 supplies the right hand side of the ODE for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the time.
//
//    Input, double Y[], the dependent variable.
//
//    Output, double YP[], the value of the derivative.
//
{
  yp[0] =   y[1];
  yp[1] = - y[0];

  return;
}
