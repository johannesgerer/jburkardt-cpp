# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "felippa.hpp"

int main ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FELIPPA_PRB.
//
//  Discussion:
//
//    FELIPPA_PRB calls the FELIPPA tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2009
//
//  Author:
//
//    John Burkardt
//
{
  int degree_max;

  timestamp ( );
  cout << "\n";
  cout << "FELIPPA_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FELIPPA library.\n";

  degree_max = 4;

  hexa_unit_monomial_test ( degree_max );
  line_unit_monomial_test ( degree_max );
  pyra_unit_monomial_test ( degree_max );
  quad_unit_monomial_test ( degree_max );
  tetr_unit_monomial_test ( degree_max );
  trig_unit_monomial_test ( degree_max );
  wedg_unit_monomial_test ( degree_max );

  degree_max = 6;
  hexa_unit_quad_test ( degree_max );

  degree_max = 10;
  line_unit_quad_test ( degree_max );

  degree_max = 5;
  pyra_unit_quad_test ( degree_max );

  degree_max = 10;
  quad_unit_quad_test ( degree_max );

  degree_max = 4;
  tetr_unit_quad_test ( degree_max );

  degree_max = 7;
  trig_unit_quad_test ( degree_max );

  degree_max = 8;
  wedg_unit_quad_test ( degree_max );

  wedg_unit_write_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FELIPPA_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
