# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "circle_rule.hpp"

int main ( );
void test01 ( int nt );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CIRCLE_RULE_PRB.
//
//  Discussion:
//
//    CIRCLE_RULE_PRB tests the CIRCLE_RULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  int nt;

  timestamp ( );
  cout << "\n";
  cout << "CIRCLE_RULE:\n";
  cout << "  C++ version\n";
  cout << "  Test the CIRCLE_RULE library.\n";

  nt = 8;
  test01 ( nt );

  nt = 32;
  test01 ( nt );
//
//  Terminate.
//
  cout << "\n";
  cout << "CIRCLE_RULE:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CIRCLE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  int e[2];
  int e1;
  int e2;
  double exact;
  int i;
  double q;
  double r8_pi = 3.141592653589793;
  double *t;
  double *w;
  double x;
  double y;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  CIRCLE_RULE can compute a rule Q(f) for the unit circle\n";
  cout << "  using NT equally spaced angles.\n";
  cout << "  Estimate integrals I(f) where f = x^e(1) * y^e(2)\n";
  cout << "  using " << nt << " points.\n";
//
//  Compute the quadrature rule.
//
  w = new double[nt];
  t = new double[nt];

  circle_rule ( nt, w, t );
//
//  Apply it to integrands.
//
  cout << "\n";
  cout << "  E(1)  E(2)    I(f)            Q(f)\n";
  cout << "\n";
//
//  Specify a monomial.
//
  for ( e1 = 0; e1 <= 6; e1 = e1 + 2 )
  {
    e[0] = e1;

    for ( e2 = e1; e2 <= 6; e2 = e2 + 2 )
    {
      e[1] = e2;

      q = 0.0;
      for ( i = 0; i < nt; i++ )
      {
        x = cos ( t[i] );
        y = sin ( t[i] );
        q = q + w[i] * pow ( x, e[0] ) * pow ( y, e[1] );
      }

      q = 2.0 * r8_pi * q;

      exact = circle01_monomial_integral ( e );

      cout << "  " << setw(2) << e[0]
           << "  " << setw(2) << e[1]
           << "  " << setw(14) << exact
           << "  " << setw(14) << q << "\n";
    }
  }

  delete [] t;
  delete [] w;

  return;
}

