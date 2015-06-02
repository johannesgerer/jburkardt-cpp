# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "toms443.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TOMS443_PRB.
//
//  Discussion:
//
//    TOMS443_PRB tests the TOMS443 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TOMS443_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TOMS443 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TOMS433_PRB\n";
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
//    TEST01 tests WEW_A
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double en;
  int n_data;
  double w1;
  double w2;
  double x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test WEW_A to evaluate\n";
  cout << "  Lambert's W function.\n";
  cout << "\n";
  cout << "          X             Exact             Computed      Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lambert_w_values ( n_data, x, w1 );

    if ( n_data <= 0 )
    {
      break;
    }

    if ( x == 0.0 )
    {
      w2 = 0.0;
    }
    else
    {
      w2 = wew_a ( x, en );
    }

    cout << setw(14) << x << "  "
         << setw(16) << w1 << "  "
         << setw(16) << w2 << "  "
         << setw(10) << fabs ( w1 - w2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests WEW_B
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double en;
  int n_data;
  double w1;
  double w2;
  double x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Test WEW_B to evaluate\n";
  cout << "  Lambert's W function.\n";
  cout << "\n";
  cout << "          X             Exact             Computed      Error\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lambert_w_values ( n_data, x, w1 );

    if ( n_data <= 0 )
    {
      break;
    }

    if ( x == 0.0 )
    {
      w2 = 0.0;
    }
    else
    {
      w2 = wew_b ( x, en );
    }

    cout << setw(14) << x << "  "
         << setw(16) << w1 << "  "
         << setw(16) << w2 << "  "
         << setw(10) << fabs ( w1 - w2 ) << "\n";
  }

  return;
}

