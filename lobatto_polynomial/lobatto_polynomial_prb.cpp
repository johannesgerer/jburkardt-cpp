# include <cmath>
# include <cstdlib>
# include <cstring>
# include <iomanip>
# include <iostream>

using namespace std;

# include "lobatto_polynomial.hpp"

int main ( );
void lobatto_polynomial_value_test ( );
void lobatto_polynomial_derivative_test ( );
void lobatto_polynomial_plot_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_POLYNOMIAL_TEST tests the LOBATTO_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LOBATTO_POLYNOMIAL_TEST:\n";
  cout << "  C++ version.\n";
  cout << "  Test the LOBATTO_POLYNOMIAL library.\n";

  lobatto_polynomial_value_test ( );
  lobatto_polynomial_derivative_test ( );
  lobatto_polynomial_plot_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LOBATTO_POLYNOMIAL_TEST\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void lobatto_polynomial_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_POLYNOMIAL_VALUE_TEST tests LOBATTO_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  double e;
  double fx1;
  double fx2;
  double *l;
  int m;
  int n;
  int n_data;
  double x;
  double xvec[1];

  m = 1;

  cout << "\n";
  cout << "LOBATTO_POLYNOMIAL_VALUE_TEST:\n";
  cout << "  LOBATTO_POLYNOMIAL_VALUES stores values of\n";
  cout << "  the completed Lobatto polynomial L(n,x).\n";
  cout << "  LOBATTO_POLYNOMIAL_VALUE evaluates the completed Lobatto polynomial.\n";
  cout << "\n";
  cout << "                                       Tabulated                 Computed\n";
  cout << "     N        X                        L(N,X)                    L(N,X)        Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lobatto_polynomial_values ( n_data, n, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    xvec[0] = x;

    l = lobatto_polynomial_value ( m, n, xvec );

    fx2 = l[0+(n-1)*m];

    e = fx1 - fx2;

    cout << "  " << setw(4) << n
         << "  " << setw(12) << x
         << "  " << setw(12) << fx1
         << "  " << setw(12) << fx2
         << "  " << setw(8) << e << "\n";

    delete [] l;
  }

  return;
}
//****************************************************************************80

void lobatto_polynomial_derivative_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_POLYNOMIAL_DERIVATIVE_TEST tests LOBATTO_POLYNOMIAL_DERIVATIVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  double e;
  double fx1;
  double fx2;
  double *lp;
  int m;
  int n;
  int n_data;
  double x;
  double xvec[1];

  m = 1;

  cout << "\n";
  cout << "LOBATTO_POLYNOMIAL_DERIVATIVE_TEST:\n";
  cout << "  LOBATTO_POLYNOMIAL_DERIVATIVES stores derivatives of\n";
  cout << "  the completed Lobatto polynomial L(n,x).\n";
  cout << "  LOBATTO_POLYNOMIAL_DERIVATIVE evaluates the completed Lobatto polynomial.\n";
  cout << "\n";
  cout << "                                       Tabulated                 Computed\n";
  cout << "     N        X                        L''(N,X)                   L''(N,X)       Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lobatto_polynomial_derivatives ( n_data, n, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    xvec[0] = x;
    lp = lobatto_polynomial_derivative ( m, n, xvec );
    fx2 = lp[0+(n-1)*m];

    e = fx1 - fx2;

    cout << "  " << setw(4) << n
         << "  " << setw(12) << x
         << "  " << setw(12) << fx1
         << "  " << setw(12) << fx2
         << "  " << setw(8) << e << "\n";

    delete [] lp;
  }

  return;
}
//****************************************************************************80

void lobatto_polynomial_plot_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_POLYNOMIAL_PLOT_TEST tests LOBATTO_POLYNOMIAL_PLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  int ndx[7] = { 1, 2, 3, 4, 5, 6, 7 };
  int ndx_num = 7;
  string prefix = "test";

  cout << "\n";
  cout << "LOBATTO_POLYNOMIAL_PLOT_TEST:\n";
  cout << "  LOBATTO_POLYNOMIAL_PLOT plots Lobatto polynomials.\n";

  lobatto_polynomial_plot ( ndx_num, ndx, prefix );

  return;
}
