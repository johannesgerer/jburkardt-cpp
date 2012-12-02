# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa076.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA076_PRB.
//
//  Discussion:
//
//    ASA076_PRB calls the ASA076 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA076_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA076 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA076_PRB:\n";
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
//    TEST01 demonstrates the use of TFN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2008
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
  cout << "  TFN computes Owen's T function.\n";
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

    t2 = tfn ( h, a );

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
//    TEST02 demonstrates the use of THA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double h;
  int n_data;
  double one = 1.0;
  double t1;
  double t2;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  THA computes Owen's T function.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "             H             A      "
       << "    T                         T  \n";
  cout << "                                  "
       << "    (Tabulated)               (THA)               DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    owen_values ( n_data, h, a, t1 );

    if ( n_data == 0 )
    {
      break;
    }

    t2 = tha ( h, one, a, one );

    cout << "  " << setprecision(4) << setw(12) << h
         << "  " << setprecision(4) << setw(12) << a
         << "  " << setprecision(16) << setw(24) << t1
         << "  " << setprecision(16) << setw(24) << t2
         << "  " << setprecision(4) << setw(10) << r8_abs ( t1 - t2 ) << "\n";
  }

  return;
}
