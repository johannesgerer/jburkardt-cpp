# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa241.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA241_PRB.
//
//  Discussion:
//
//    ASA241_PRB tests the ASA241 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA241_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA241 library.\n";

  test01 ( );
  test02 ( );

  cout << "\n";
  cout << "ASA241_PRB:\n";
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
//    TEST01 tests R4_NORMAL_01_CDF_INVERSE, NORMAL_01_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  float fx2;
  int n_data;
  double x;
  float x2;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Let FX = NormalCDF ( X ).\n";
  cout << "\n";
  cout << "  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).\n";
  cout << "\n";
  cout << "  R4_NORMAL_01_CDF_INVERSE takes the value of FX, and computes an\n";
  cout << "    estimate X2, of the corresponding input argument,\n";
  cout << "    accurate to about 7 decimal places.\n";
  cout << "\n";
  cout << "          FX                        X                        X2          DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = ( float ) ( fx );
    x2 = r4_normal_01_cdf_inverse ( fx2 );

    cout << "  " << setprecision(16) << setw(24) << fx  
         << "  " << setprecision(16) << setw(24) << x 
         << "  " << setprecision(16) << setw(24) << x2 
         << "  " << setprecision(16) << fabs ( x - x2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R8_NORMAL_01_CDF_INVERSE, NORMAL_01_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Let FX = NormalCDF ( X ).\n";
  cout << "\n";
  cout << "  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).\n";
  cout << "\n";
  cout << "  R8_NORMAL_01_CDF_INVERSE takes the value of FX, and computes an\n";
  cout << "    estimate X2, of the corresponding input argument,\n";
  cout << "    accurate to about 16 decimal places.\n";
  cout << "\n";
  cout << "          FX                        X                        X2          DIFF\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = r8_normal_01_cdf_inverse ( fx );

    cout << "  " << setprecision(16) << setw(24) << fx  
         << "  " << setprecision(16) << setw(24) << x 
         << "  " << setprecision(16) << setw(24) << x2 
         << "  " << setprecision(16) << fabs ( x - x2 ) << "\n";
  }

  return;
}
