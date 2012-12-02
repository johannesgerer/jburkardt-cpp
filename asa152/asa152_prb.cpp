# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa152.hpp"

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
//    MAIN is the main program for ASA152_PRB.
//
//  Discussion:
//
//    ASA152_PRB calls the ASA152 routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA152_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA152 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA152_PRB:\n";
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
//    TEST01 demonstrates CHYPER for cumulative probabilities.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int ifault;
  int n_data;
  bool point;
  int pop;
  int sam;
  int suc;
  int x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  CHYPER computes cumulative probablities\n";
  cout << "  of the hypergeometric PDF.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "   SAM   SUC   POP     X    "
       << "  CDF                       CDF                     DIFF\n";
  cout << "                            "
       << " (tabulated)               (CHYPER)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_cdf_values ( &n_data, &sam, &suc, &pop, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    point = false;

    fx2 = chyper ( point, sam, x, pop, suc, &ifault );

    cout << "  " << setw(4) << sam
         << "  " << setw(4) << suc
         << "  " << setw(4) << pop
         << "  " << setw(4) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << r8_abs ( fx - fx2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 demonstrates CHYPER for point probabilities.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int ifault;
  int n_data;
  bool point;
  int pop;
  int sam;
  int suc;
  int x;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  CHYPER computes point probablities\n";
  cout << "  of the hypergeometric PDF.\n";
  cout << "  Compare to tabulated values.\n";
  cout << "\n";
  cout << "   SAM   SUC   POP     X    "
       << "  PDF                       PDF                     DIFF\n";
  cout << "                            "
       << " (tabulated)               (CHYPER)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_pdf_values ( &n_data, &sam, &suc, &pop, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    point = true;

    fx2 = chyper ( point, sam, x, pop, suc, &ifault );

    cout << "  " << setw(4) << sam
         << "  " << setw(4) << suc
         << "  " << setw(4) << pop
         << "  " << setw(4) << x
         << "  " << setprecision(16) << setw(24) << fx
         << "  " << setprecision(16) << setw(24) << fx2
         << "  " << setprecision(4) << setw(10) << r8_abs ( fx - fx2 ) << "\n";
  }

  return;
}
