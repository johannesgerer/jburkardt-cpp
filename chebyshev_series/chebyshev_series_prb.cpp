# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>

using namespace std;

# include "chebyshev_series.hpp"

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
//    MAIN is the main program for CHEBYSHEV_SERIES_PRB.
//
//  Discussion:
//
//    CHEBYSHEV_SERIES_PRB tests the CHEBYSHEV_SERIES library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
{
  timestamp ();
  cout << "\n";
  cout << "CHEBYSHEV_SERIES_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the CHEBYSHEV_SERIES libary.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CHEBYSHEV_SERIES_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ();

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 considers an even Chebyshev series for EXP(X).
//
//  Discussion:
//
//    Table 5 is from Clenshaw, and contains 18 terms of the Chebyshev
//    series for exp(x) over [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
{
  int i;
  double s;
  double s1;
  double s2;
  double s3;
  double table5[18] = { 
    2.53213175550401667120,
    1.13031820798497005442,
    0.27149533953407656237,
    0.04433684984866380495,
    0.00547424044209373265,
    0.00054292631191394375,
    0.00004497732295429515,
    0.00000319843646240199,
    0.00000019921248066728,
    0.00000001103677172552,
    0.00000000055058960797,
    0.00000000002497956617,
    0.00000000000103915223,
    0.00000000000003991263,
    0.00000000000000142376,
    0.00000000000000004741,
    0.00000000000000000148,
    0.00000000000000000004 };
  double x;
  double y;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  ECHEBSER3 computes a Chebyshev series approximation\n";
  cout << "  and the first three derivatives.\n";
  cout << "\n";
  cout << "  Errors of a Chebyshev series for exp(x)\n";
  cout << "\n";
  cout << "    x        err(y)       err(y')       err(y\")      err(y\"')\n";
  cout << "\n";

  for ( i = -10; i <= 10; i++ )
  {
    x = ( double ) i / 10.0;
    s = echebser3 ( x, table5, 18, s1, s2, s3 );
    y = exp ( x );
    s = s - y;
    s1 = s1 - y;
    s2 = s2 - y;
    s3 = s3 - y;

    cout << setw(5)  << x
         << setw(14) << s
         << setw(14) << s1
         << setw(14) << s2
         << setw(14) << s3 << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 considers an even Chebyshev series for COS(PI*X/2).
//
//  Discussion:
//
//    TABLE1 contains the even Chebyshev series coefficients for
//    cos(pi*x/2) over [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
{
  int i;
  double s;
  double s1;
  double s2;
  double table1[11] = {
    +0.94400243153646953490,
    -0.49940325827040708740,
    +0.02799207961754761751,
    -0.00059669519654884650,
    +0.00000670439486991684,
    -0.00000004653229589732,
    +0.00000000021934576590,
    -0.00000000000074816487,
    +0.00000000000000193230,
    -0.00000000000000000391,
    +0.00000000000000000001 };
  double x;
  double y;
  double y1;
  double y2;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  EVENCHEBSER2 computes an even Chebyshev series\n";
  cout << "  and its first two derivatives.\n";
  cout << "\n";
  cout << "  Errors of an even Chebyshev series for cos(pi*x/2):\n";
  cout << "\n";
  cout << "    x        err(y)       err(y')       err(y\")\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    x = ( double ) i / 10.0;
    s = evenchebser2 ( x, table1, 11, s1, s2 );
    sincos ( M_PI_2 * x, &y1, &y );
    y1 = - y1 * M_PI_2;
    y2 = - y * (M_PI_2 * M_PI_2);
    s = s - y;
    s1 = s1 - y1;
    s2 = s2 - y2;

    cout << setw(5)  << x
         << setw(14) << s
         << setw(14) << s1
         << setw(14) << s2 << "\n";
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 considers an odd Chebyshev series for SINH(X).
//
//  Discussion:
//
//    TABLE5ODD contains the odd Chebyshev series coefficients for
//    sinh(x) over -1 <= x <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
{
  int i;
  double s;
  double s1;
  double s2;
  double table5odd[9] = {
    1.13031820798497005442,
    0.04433684984866380495,
    0.00054292631191394375,
    0.00000319843646240199,
    0.00000001103677172552,
    0.00000000002497956617,
    0.00000000000003991263,
    0.00000000000000004741,
    0.00000000000000000004 };
  double x;
  double y;
  double y1;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  ODDCHEBSER2 computes an odd Chebyshev series approximation.\n";
  cout << "  and its first two derivatives.\n";
  cout << "\n";
  cout << "  Errors of an odd Chebyshev series y(x) approximating sinh(x):\n";
  cout << "\n";
  cout << "    x        err(y)       err(y')       err(y\")\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    x = ( double ) ( i ) / 10.0;
    s = oddchebser2 ( x, table5odd, 9, s1, s2 );
    y = sinh ( x );
    y1 = cosh ( x );
    s = s - y;
    s1 = s1 - y1;
    s2 = s2 - y;
    cout << setw(5)  << x
         << setw(14) << s
         << setw(14) << s1
         << setw(14) << s2 << "\n";
  }

  return;
}

