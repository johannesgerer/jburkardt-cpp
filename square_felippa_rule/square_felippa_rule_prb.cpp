# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "square_felippa_rule.hpp"

int main ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SQUARE_FELIPPA_RULE_PRB.
//
//  Discussion:
//
//    SQUARE_FELIPPA_RULE_PRB tests the SQUARE_FELIPPA_RULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  int degree_max;

  timestamp ( );
  cout << "\n";
  cout << "SQUARE_FELIPPA_RULE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SQUARE_FELIPPA_RULE library.\n";

  degree_max = 4;
  square_monomial_test ( degree_max );

  degree_max = 5;
  square_quad_test ( degree_max );
//
//  Terminate.
//
  cout << "\n";
  cout << "SQUARE_FELIPPA_RULE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
