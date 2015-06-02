# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "line_nco_rule.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LINE_NCO_RULE_PRB.
//
//  Discussion:
//
//    LINE_NCO_RULE_PRB tests the LINE_NCO_RULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LINE_NCO_RULE_PRB\n";
  cout << "  C++ version:\n";
  cout << "  Test the LINE_NCO_RULE library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LINE_NCO_RULE_PRB\n";
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
//    TEST01 computes and prints NCO rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int n;
  double *w;
  double w_sum;
  double *x;

  a = -1.0;
  b = +1.0;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  LINE_NCO_RULE computes the Newton-Cotes Open (NCO) rule\n";
  cout << "  using N equally spaced points for an interval [A,B].\n";

  for ( n = 1; n <= 12; n++ )
  {
    x = new double[n];
    w = new double[n];

    line_nco_rule ( n, a, b, x, w );
    cout << "\n";
    cout << "  Newton-Cotes Open (NCO) Rule #" << n << "\n";
    cout << "   I       X(I)            W(I)\n";
    cout << "\n";
    w_sum = 0.0;
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << i
           << "  " << setw(14) << x[i]
           << "  " << setw(14) << w[i] << "\n";
      w_sum = w_sum + fabs ( w[i] );
    }
    cout << "        Sum(|W)|) =  " << setw(14) <<  w_sum << "\n";

    delete [] x;
    delete [] w;
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses NCO rules to estimate the integral of exp(x) from 0 to 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error;
  double exact;
  int i;
  int n;
  double q;
  double *w;
  double *x;

  a =  0.0;
  b = +1.0;
  exact = exp ( b ) - exp ( a );

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Use a sequence of NCO rules to compute an estimate Q\n";
  cout << "  of the integral:\n";
  cout << "    I = integral ( 0 <= x <= 1 ) exp(x) dx.\n";
  cout << "  The exact value is:\n";
  cout << "    I = " << exact << "\n";

  cout << "\n";
  cout << "   N       Q             |Q-I|\n";
  cout << "\n";

  for ( n = 1; n <= 22; n++ )
  {
    x = new double[n];
    w = new double[n];

    line_nco_rule ( n, a, b, x, w );

    q = 0.0;
    for ( i = 0; i < n; i++ )
    {
      q = q + w[i] * exp ( x[i] );
    }
    error = fabs ( exact - q );
    cout << "  " << setw(2) << n
         << "  " << setw(14) << q
         << "  " << setw(14) << error << "\n";

    delete [] x;
    delete [] w;
  }
  return;
}
