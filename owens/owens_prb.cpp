# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "owens.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for OWENS_PRB.
//
//  Discussion:
//
//    OWENS_PRB calls the OWENS routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "OWENS_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the OWENS library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "OWENS_PRB:\n";
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
//    TEST01 demonstrates the use of T.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double h;
  int n_data;
  double t1;
  double t2;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  T computes Owen's T function.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "             H             A      "
       << "    T                         T  \n";
  cout << "                                  "
       << "    (Tabulated)               (TFN)               DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    owen_values ( n_data, h, a, t1 );

    if ( n_data == 0 )
    {
      break;
    }

    t2 = t ( h, a );

    cout << "  " << setprecision(4) << setw(12) << h
         << "  " << setprecision(4) << setw(12) << a
         << "  " << setprecision(16) << setw(24) << t1
         << "  " << setprecision(16) << setw(24) << t2
         << "  " << setprecision(4) << setw(10) << r8_abs ( t1 - t2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 demonstrates the use of BIVNOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fxy1;
  double fxy2;
  int n_data;
  double r;
  double x;
  double y;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  BIVNOR computes the bivariate normal probability.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "          X               Y               "
       << "R           P                         P                       DIFF\n";
  cout << "                                          "
       << "           (Tabulated)               (BIVNOR)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bivariate_normal_cdf_values ( n_data, x, y, r, fxy1 );

    if ( n_data == 0 )
    {
      break;
    }

    fxy2 = bivnor ( - x, - y, r );

    cout << "  " << setprecision(4) << setw(12) << x
         << "  " << setprecision(4) << setw(12) << y
         << "  " << setprecision(4) << setw(12) << r
         << "  " << setprecision(16) << setw(24) << fxy1
         << "  " << setprecision(16) << setw(24) << fxy2
         << "  " << setprecision(4) << setw(10) << r8_abs ( fxy1 - fxy2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 demonstrates the use of ZNORM1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double fx1;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  ZNORM1 computes the normal CDF starting at 0.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "          X           P                         P                       DIFF\n";
  cout << "                     (Tabulated)               (ZNORM1)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx1 = fx1 - 0.5;

    fx2 = znorm1 ( x );

    cout << "  " << setprecision(4) << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx1
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << r8_abs ( fx1 - fx2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 demonstrates the use of ZNORM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double fx1;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  ZNORM2 computes the complementary normal CDF.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "          X           P                         P                       DIFF\n";
  cout << "                     (Tabulated)               (ZNORM2)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx1 = 1.0 - fx1;

    fx2 = znorm2 ( x );

    cout << "  " << setprecision(4) << setw(12) << x
         << "  " << setprecision(16) << setw(24) << fx1
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << r8_abs ( fx1 - fx2 ) << "\n";
  }

  return;
}
