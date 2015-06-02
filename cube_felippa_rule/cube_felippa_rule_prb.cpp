# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "cube_felippa_rule.hpp"

int main ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CUBE_FELIPPA_RULE_PRB.
//
//  Discussion:
//
//    CUBE_FELIPPA_RULE_PRB tests the CUBE_FELIPPA_RULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  int degree_max;

  timestamp ( );
  cout << "\n";
  cout << "CUBE_FELIPPA_RULE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CUBE_FELIPPA_RULE library.\n";

  degree_max = 4;
  cube_monomial_test ( degree_max );

  degree_max = 6;
  cube_quad_test ( degree_max );
//
//  Terminate.
//
  cout << "\n";
  cout << "CUBE_FELIPPA_RULE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
