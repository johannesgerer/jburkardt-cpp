# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "zero_rc.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ZERO_RC_PRB.
//
//  Discussion:
//
//    ZERO_RC_PRB tests the ZERO_RC library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ZERO_RC_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ZERO_RC library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ZERO_RC_PRB:\n";
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
//    TEST01 tests ROOT_RC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double ferr;
  double fx;
  int i;
  int it;
  int it_max;
  double q[9];
  double x;
  double xerr;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test ROOT_RC, which searches for an\n";
  cout << "  approximate solution of F(X) = 0.\n";
  cout << "\n";
  cout << "       X              XERR            FX              FERR\n";
  cout << "\n";
//
//  Initialization.
//
  it = 0;
  it_max = 30;
  for ( i = 0; i < 9; i++ )
  {
    q[i] = 0.0;
  }
  x = - 2.1;
//
//  Each call takes one more step of improvement.
//
  for ( ; ; )
  {
    fx = cos ( x ) - x;

    if ( it == 0 )
    {
      cout << "  " << setw(14) << x
           << "  " << "              "
           << "  " << setw(14) << fx << "\n";
    }
    else
    {
      cout << "  " << setw(14) << x
           << "  " << setw(14) << xerr
           << "  " << setw(14) << fx
           << "  " << setw(14) << ferr << "\n";
    }

    x = root_rc ( x, fx, ferr, xerr, q );

    if ( ferr < 1.0E-08 )
    {
      cout << "\n";
      cout << "  Uncertainty in F(X) less than tolerance\n";
      break;
    }

    if ( xerr < 1.0E-08 )
    {
      cout << "\n";
      cout << "  Width of X interal less than tolerance.\n";
      break;
    }

    if ( it_max < it )
    {
      cout << "\n";
      cout << "  Too many iterations!'\n";
      break;
    }
    it = it + 1;     
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests ROOTS_RC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double ferr;
  double *fx;
  int i;
  int it;
  int it_max = 30;
  int j;
  int n = 4;
  double *q;
  double *x;
  double *xnew;

  fx = new double[n];
  q = new double[(2*n+2)*(n+2)];
  x = new double[n];
  xnew = new double[n];

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Test ROOTS_RC, which seeks a solution of\n";
  cout << "  the N-dimensional nonlinear system F(X) = 0.\n";
  cout << "\n";
  cout << "       FERR          X\n";
  cout << "\n";
//
//  Initialization.
//
  for ( j = 0; j < n + 2; j++ )
  {
    for ( i = 0; i < 2 * n + 2; i++ )
    {
      q[i+j*(2*n+2)] = 0.0;
    }
  }

  xnew[0] = 1.2;
  for ( i = 1; i < n; i++ )
  {
    xnew[i] = 1.0;
  }

  it = 0;

  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = xnew[i];
    }

    fx[0] = 1.0 - x[0];
    for ( i = 1; i < n; i++ )
    {
      fx[i] = 10.0 * ( x[i] - x[i-1] * x[i-1] );
    }

    if ( it == 0 )
    {
      cout << "                ";
    }
    else
    {
      cout << "  " << setw(14) << ferr;
    }
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(14) << x[i];
    }
    cout << "\n";

    roots_rc ( n, x, fx, ferr, xnew, q );

    if ( ferr < 1.0E-07 )
    {
      cout << "\n";
      cout << "  Sum of |f(x)| less than tolerance.\n";
      break;
    }

    if ( it_max < it )
    {
      cout << "\n";
      cout << "  Too many iterations!\n";
      break;
    }
    it = it + 1;
  }

  delete [] fx;
  delete [] q;
  delete [] x;
  delete [] xnew;

  return;
}
