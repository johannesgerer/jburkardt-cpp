# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa066.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA066_PRB.
//
//  Discussion:
//
//    ASA066_PRB tests the ASA066 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA066_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA066 library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA066_PRB:\n";
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
//    TEST01 compares ALNORM against tabulated values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  bool upper = false;
  double x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Compare tabulated values of the normal\n";
  cout << "  Cumulative Density Function against values\n";
  cout << "  computed by ALNORM.\n";
  cout << "\n";
  cout << "         X        CDF                       CDF"
       << "                    DIFF\n";
  cout << "               (tabulated)                 (ALNORM)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = alnorm ( x, upper );

    cout << "  " << setprecision(4)  << setw(10) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4)  << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 compares NORMP against tabulated values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double p;
  double pdf;
  double q;
  double x;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Compare tabulated values of the normal\n";
  cout << "  Cumulative Density Function against values\n";
  cout << "  computed by NORMP.\n";
  cout << "\n";
  cout << "         X        CDF                       CDF"
       << "                    DIFF\n";
  cout << "               (tabulated)                 (NORMP)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    normp ( x, &p, &q, &pdf );
    fx2 = p;

    cout << "  " << setprecision(4)  << setw(10) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4)  << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 compares NPROB against tabulated values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double p;
  double pdf;
  double q;
  double x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Compare tabulated values of the normal\n";
  cout << "  Cumulative Density Function against values\n";
  cout << "  computed by NPROBP.\n";
  cout << "\n";
  cout << "         X        CDF                       CDF"
       << "                    DIFF\n";
  cout << "               (tabulated)                 (NPROB)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    nprob ( x, &p, &q, &pdf );
    fx2 = p;

    cout << "  " << setprecision(4)  << setw(10) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4)  << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
