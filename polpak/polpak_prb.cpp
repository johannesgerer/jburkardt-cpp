# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "polpak.hpp"

int main ( );

void test001 ( );
void test002 ( );
void test003 ( );
void test0035 ( );
void test004 ( );
void test005 ( );
void test006 ( );
void test007 ( );
void test008 ( );
void test010 ( );
void test0102 ( );

void test011 ( );
void test012 ( );
void test013 ( );
void test014 ( );
void test015 ( );
void test016 ( );
void test017 ( );
void test0175 ( );
void test018 ( );
void test0185 ( );
void test019 ( );

void test020 ( );
void test021 ( );
void test0215 ( );
void test0216 ( );
void test0217 ( );
void test0218 ( );
void test024 ( );
void test0243 ( );
void test0104 ( );
void test0115 ( );
void test0265 ( );
void test028 ( );
void test0245 ( );
void test025 ( );
void test0255 ( );
void test026 ( );
void test027 ( );
void test029 ( );

void test031 ( );
void test032 ( );
void test033 ( );
void test034 ( );
void test036 ( );
void test0365 ( );
void test037 ( );
void test038 ( );
void test039 ( );

void test040 ( );
void test041 ( );
void test042 ( );
void test0425 ( );
void test023 ( );
void test043 ( );
void test044 ( );
void test045 ( );
void test046 ( );
void test047 ( );
void test048 ( );
void test049 ( );

void test050 ( );
void test0505 ( );
void test051 ( );
void test052 ( );
void test054 ( );
void test055 ( );
void test057 ( );
void test058 ( );
void test059 ( );
void test0595 ( );

void test060 ( );
void test061 ( );
void test0615 ( );
void test062 ( );
void test0623 ( );
void test0625 ( );
void test063 ( );
void test0635 ( );
void test064 ( );
void test065 ( );
void test066 ( );
void test0665 ( );
void test0667 ( );
void test067 ( );
void test0675 ( );
void test068 ( );
void test022 ( );
void test0685 ( );
void test06855 ( );
void test06856 ( );
void test069 ( );
void test0695 ( );
void test0696 ( );
void test0697 ( );

void test070 ( );
void test071 ( );
void test072 ( );
void test073 ( );
void test074 ( );
void test075 ( );
void test076 ( );
void test077 ( );
void test0773 ( );
void test0775 ( );
void test078 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POLPAK_PRB.
//
//  Discussion:
//
//    POLPAK_PRB calls the POLPAK test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "POLPAK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the POLPAK library.\n";

  test001 ( );
  test002 ( );
  test003 ( );
  test0035 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test008 ( );
  test010 ( );
  test0102 ( );
  test0104 ( );

  test011 ( );
  test0115 ( );
  test012 ( );
  test013 ( );
  test014 ( );
  test015 ( );
  test016 ( );
  test017 ( );
  test0175 ( );
  test018 ( );
  test0185 ( );
  test019 ( );

  test020 ( );
  test021 ( );
  test0215 ( );
  test0216 ( );
  test0217 ( );
  test0218 ( );
  test024 ( );
  test0243 ( );
  test0245 ( );
  test025 ( );
  test0255 ( );
  test026 ( );
  test0265 ( );
  test028 ( );
  test027 ( );
  test029 ( );

  test031 ( );
  test032 ( );
  test033 ( );
  test034 ( );
  test036 ( );
  test0365 ( );
  test037 ( );
  test038 ( );
  test039 ( );

  test040 ( );
  test041 ( );
  test042 ( );
  test0425 ( );
  test023 ( );
  test043 ( );
  test044 ( );
  test045 ( );
  test046 ( );
  test047 ( );
  test048 ( );
  test049 ( );

  test050 ( );
  test0505 ( );
  test051 ( );
  test052 ( );
  test054 ( );
  test055 ( );
  test057 ( );
  test058 ( );
  test059 ( );
  test0595 ( );

  test060 ( );
  test061 ( );
  test0615 ( );
  test062 ( );
  test0623 ( );
  test0625 ( );
  test063 ( );
  test0635 ( );
  test064 ( );
  test065 ( );
  test066 ( );
  test0665 ( );
  test0667 ( );
  test067 ( );
  test0675 ( );
  test068 ( );
  test022 ( );
  test0685 ( );
  test06855 ( );
  test06856 ( );
  test069 ( );
  test0695 ( );
  test0696 ( );
  test0697 ( );

  test070 ( );
  test071 ( );
  test072 ( );
  test073 ( );
  test074 ( );
  test075 ( );
  test076 ( );
  test077 ( );
  test0773 ( );
  test0775 ( );
  test078 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "POLPAK_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test001 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST001 tests AGM and AGM_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n_data;

  cout << "\n";
  cout << "TEST001\n";
  cout << "  AGM computes the arithmetic geometric mean.\n";
  cout << "\n";
  cout << "           A           B         "
       << "   AGM                       AGM               Diff\n";
  cout << "                             "
       << "      (Tabulated)                AGM(A,B)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    agm_values ( &n_data, &a, &b, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = agm ( a, b );

    cout << "  " << setprecision(6)  << setw(10) << a   
         << "  " << setprecision(6)  << setw(10) << b   
         << "  " << setprecision(16) << setw(24) << fx  
         << "  " << setprecision(16) << setw(24) << fx2 
         << "  " << setprecision(6)  << setw(10) << r8_abs ( fx - fx2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test002 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST002 tests AGUD and GUD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double gamma;
  int i;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST002\n";
  cout << "  AGUD computes the inverse Gudermannian;\n";
  cout << "  GUD computes the Gudermannian.\n";
  cout << "\n";
  cout << "         X     GUD(X)     AGUD(GUD(X))\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    gamma = gud ( x );
    x2 = agud ( gamma );

    cout << "  " << setw(10) << x
         << "  " << setw(10) << gamma 
         << "  " << setw(10) << x2    << "\n";
  }

  return;
}
//****************************************************************************80

void test003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST003 tests ALIGN_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M_MAX 10
# define N_MAX 10

  int i;
  int j;

  cout << "\n";
  cout << "TEST003\n";
  cout << "  ALIGN_ENUM counts the number of possible\n";
  cout << "  alignments of two biological sequences.\n";

  cout << "\n";
  cout << "  Alignment enumeration table:\n";
  cout << "\n";

  cout << "      ";
  for ( j = 0; j <= 5; j++ )
  {
    cout << setw(8) << j << "  ";
  }
  cout << "\n";
  cout << "\n";

  for ( i = 0; i <= M_MAX; i++ )
  {
    cout << "  " << setw(2) << i << "  ";
    for ( j = 0; j <= 5; j++ )
    {
      cout << setw(8) << align_enum ( i, j ) << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "      ";
  for ( j = 6; j <= N_MAX; j++ )
  {
    cout << setw(8) << j << "  ";
  }
  cout << "\n";
  cout << "\n";

  for ( i = 0; i <= M_MAX; i++ )
  {
    cout << "  " << setw(2) << i << "  ";
    for ( j = 6; j <= N_MAX; j++ )
    {
      cout << setw(8) << align_enum ( i, j ) << "  ";
    }
    cout << "\n";
  }
  return;
# undef M_MAX
# undef N_MAX
}
//****************************************************************************80

void test0035 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0035 tests ARC_COSINE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST0035\n";
  cout << "  ARC_COSINE computes the inverse cosine,\n";
  cout << "  and chops input arguments that are out of bounds.\n";
  cout << "\n";
  cout << "         X     ARC_COSINE(X)     COS(ARC_COSINE(X))\n";
  cout << "\n";

  for ( i = -5; i <= 5; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    a = arc_cosine ( x );
    x2 = cos ( a );

    cout << "  " << setw(10) << x  
         << "  " << setw(10) << a  
         << "  " << setw(10) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void test004 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST004 tests ASINH2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST004\n";
  cout << "  ASINH2 computes the inverse hyperbolic sine\n";
  cout << "  of a given value.\n";
  cout << "\n";
  cout << "         X     ASINH2(X)     SINH(ASINH2(X))\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    a = asinh2 ( x );
    x2 = sinh ( a );

    cout                   << "  "
         << setw(10) << x  << "  "
         << setw(10) << a  << "  "
         << setw(10) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests ATAN4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double x;
  double y;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  ATAN4 computes the arc-tangent given Y and X;\n";
  cout << "  ATAN2 is the system version of this routine.\n";
  cout << "\n";
  cout << "     X         Y     ATAN2(Y,X)   ATAN4(Y,X)\n";
  cout << "\n";

  x = 1.0;
  y = 0.0;
  cout << setw(10) << x << "  "
       << setw(10) << y << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << atan4 ( y, x ) << "\n";

  x = 1.0;
  y = 1.0;
  cout << setw(10) << x << "  "
       << setw(10) << y << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << atan4 ( y, x ) << "\n";

  x = 0.0;
  y = 1.0;
  cout << setw(10) << x << "  "
       << setw(10) << y << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << atan4 ( y, x ) << "\n";

  x = -1.0;
  y = 1.0;
  cout << setw(10) << x << "  "
       << setw(10) << y << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << atan4 ( y, x ) << "\n";

  x = -1.0;
  y = 0.0;
  cout << setw(10) << x << "  "
       << setw(10) << y << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << atan4 ( y, x ) << "\n";

  x = - 1.0;
  y = - 1.0;
  cout << setw(10) << x << "  "
       << setw(10) << y << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << atan4 ( y, x ) << "\n";

  x =   0.0;
  y = - 1.0;
  cout << setw(10) << x << "  "
       << setw(10) << y << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << atan4 ( y, x ) << "\n";

  x =   1.0;
  y = - 1.0;
  cout << setw(10) << x << "  "
       << setw(10) << y << "  "
       << setw(10) << atan2 ( y, x ) << "  "
       << setw(10) << atan4 ( y, x ) << "\n";

  return;
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests ATANH2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST006\n";
  cout << "  ATANH2 computes the inverse hyperbolic tangent\n";
  cout << "  of a given value.\n";
  cout << "\n";
  cout << "         X     ATANH2(X)     TANH(ATANH2(X))\n";
  cout << "\n";

  for ( i = -2; i <= 9; i++ )
  {
    x = ( ( double ) i ) / 10.0;
    a = atanh2 ( x );
    x2 = tanh ( a );

    cout                   << "  "
         << setw(10) << x  << "  "
         << setw(10) << a  << "  "
         << setw(10) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void test007 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST007 tests BELL and BELL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int *c2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST007\n";
  cout << "  BELL computes Bell numbers.\n";
  cout << "  BELL_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N  exact C(I)  computed C(I)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bell_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = new int[n+1];

    bell ( n, c2 );

    cout                     << "  "
         << setw(4) << n     << "  "
         << setw(8) << c     << "  "
         << setw(8) << c2[n] << "\n";

    delete [] c2;

  }

  return;
}
//****************************************************************************80

void test008 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST008 tests BENFORD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "TEST008\n";
  cout << "  BENFORD(I) is the Benford probability of the\n";
  cout << "  initial digit sequence I.\n";
  cout << "\n";
  cout << "     I  BENFORD(I)\n";
  cout << "\n";

  for ( i = 1; i <= 9; i++ )
  {
    cout                              << "  "
         << setw(4) << i              << "  "
         << setw(10) << benford ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void test010 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST010 tests BERNOULLI_NUMBER and BERNOULLI_NUMBER_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double c0;
  double c1[31];
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST010\n";
  cout << "  BERNOULLI_NUMBER computes Bernoulli numbers;\n";
  cout << "  BERNOULLI_NUMBER_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "   I      Exact     BERNOULLI_NUMBER\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c0 );

    if ( n_data == 0 )
    {
      break;
    }

    bernoulli_number ( n, c1 );

    cout                      << "  "
         << setw(4)  << n     << "  "
         << setw(10) << c0    << "  "
         << setw(10) << c1[n] << "\n";
  }

  return;
}
//****************************************************************************80

void test0102 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0102 tests BERNOULLI_NUMBER2 and BERNOULLI_NUMBER_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double c0;
  double c1[31];
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST0102\n";
  cout << "  BERNOULLI_NUMBER2 computes Bernoulli numbers;\n";
  cout << "  BERNOULLI_NUMBER_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "   I      Exact     BERNOULLI_NUMBER2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c0 );

    if ( n_data == 0 )
    {
      break;
    }

    bernoulli_number2 ( n, c1 );

    cout                      << "  "
         << setw(4)  << n     << "  "
         << setw(10) << c0    << "  "
         << setw(10) << c1[n] << "\n";
  }

  return;
}
//****************************************************************************80

void test0104 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0104 tests BERNOULLI_NUMBER3 and BERNOULLI_NUMBER_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double c0;
  double c1;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST0104\n";
  cout << "  BERNOULLI_NUMBER3 computes Bernoulli numbers.\n";
  cout << "  BERNOULLI_NUMBER_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "   I      Exact     BERNOULLI_NUMBER3\n";
  cout << "\n";
  
  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c0 );

    if ( n_data == 0 )
    {
      break;
    }

    c1 = bernoulli_number3 ( n );

    cout                   << "  "
         << setw(4)  << n  << "  "
         << setw(14) << c0 << "  "
         << setw(14) << c1 << "\n";

  }
 
  return;
}
//****************************************************************************80

void test011 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST011 tests BERNOULLI_POLY;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bx;
  int i;
  int n = 15;
  double x;

  x = 0.2;

  cout << "\n";
  cout << "TEST011\n";
  cout << "  BERNOULLI_POLY evaluates Bernoulli polynomials;\n";
  cout << "\n";
  cout << "  X = " << x << "\n";
  cout << "\n";
  cout << "  I          BX\n";
  cout << "\n";

  for ( i = 1; i <= n; i++ )
  {
    bx = bernoulli_poly ( i, x );

    cout                   << "  "
         << setw(6)  << i  << "  "
         << setw(10) << bx << "\n";
  }

  return;
}
//****************************************************************************80

void test0115 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0115 tests BERNOULLI_POLY2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bx;
  int i;
  int n = 15;
  double x;

  x = 0.2;
 
  cout << "\n";
  cout << "TEST0115\n";
  cout << "  BERNOULLI_POLY2 evaluates Bernoulli polynomials.\n";
  cout << "\n";
  cout << "  X = " << x << "\n";
  cout << "\n";
  cout << "  I          BX\n";
  cout << "\n";
 
  for ( i = 1; i <= n; i++ )
  {
    bx = bernoulli_poly2 ( i, x );

    cout                   << "  "
         << setw(2)  << i  << "  "
         << setw(16) << bx << "\n";
  }
 
  return;
}
//****************************************************************************80

void test012 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST012 tests BETA and BETA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fxy;
  double fxy2;
  int n_data;
  double x;
  double y;

  cout << "\n";
  cout << "TEST012:\n";
  cout << "  BETA evaluates the Beta function.\n";
  cout << "  BETA_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     X      Y        Exact F       BETA(X,Y)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_values ( &n_data, &x, &y, &fxy );

    if ( n_data == 0 )
    {
      break;
    }

    fxy2 = beta ( x, y );

    cout                     << "  "
         << setw(10) << x    << "  "
         << setw(10) << y    << "  "
         << setw(10) << fxy  << "  "
         << setw(10) << fxy2 << "\n";
  }

  return;
}
//****************************************************************************80

void test013 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST013 tests BERNSTEIN_POLY and BERNSTEIN_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double bvec[11];
  int k;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST013:\n";
  cout << "  BERNSTEIN_POLY evaluates the Bernstein polynomials.\n";
  cout << "  BERNSTEIN_POLY_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "   N   K   X   Exact   B(N,K)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bernstein_poly_values ( &n_data, &n, &k, &x, &b );

    if ( n_data == 0 )
    {
      break;
    }

    bernstein_poly ( n, x, bvec );

    cout << "  " << setw(4)  << n
         << "  " << setw(4)  << k
         << "  " << setw(7)  << x
         << "  " << setw(14) << b
         << "  " << setw(14) << bvec[k] << "\n";
  }

  return;
}
//****************************************************************************80

void test014 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST014 tests BPAB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a;
  double b;
  double bern[N+1];
  int i;
  double x;

  cout << "\n";
  cout << "TEST014\n";
  cout << "  BPAB evaluates Bernstein polynomials.\n";
  cout << "\n";

  x = 0.3;
  a = 0.0;
  b = 1.0;

  bpab ( N, x, a, b, bern );

  cout << "  The Bernstein polynomials of degree " << N << "\n";
  cout << "  based on the interval from " << a << "\n";
  cout << "  to " << b << "\n";
  cout << "  evaluated at X = " << x << "\n";
  cout << "\n";

  for ( i = 0; i <= N; i++ )
  {
    cout << "  " << setw(4)  << i       
         << "  " << setw(14) << bern[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST015 tests CARDAN and CARDAN_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 10

  double c[N_MAX+1];
  double cx1;
  double *cx2;
  int i;
  int n;
  double s;
  double x;

  s = 1.0;

  cout << "\n";
  cout << "TEST015\n";
  cout << "  CARDAN_POLY_COEF returns the coefficients of a\n";
  cout << "  Cardan polynomial.\n";
  cout << "  CARDAN evaluates a Cardan polynomial directly.\n";
  cout << "\n";
  cout << "  We use the parameter S = " << s << "\n";
  cout << "\n";
  cout << "  Table of polynomial coefficients:\n";
  cout << "\n";

  for ( n = 0; n <= N_MAX; n++ )
  {
    cardan_poly_coef ( n, s, c );
    cout << "  "
         << setw(2) << n << "  ";
    for ( i = 0; i <= n; i++ )
    {
      cout << setw(5) << c[i] << "  ";
    }
    cout << "\n";
  }

  s = 0.5;
  x = 0.25;

  cout << "\n";
  cout << "  Compare CARDAN_POLY_COEF + R8POLY_VAL_HORNER\n";
  cout << "  versus CARDAN alone.\n";
  cout << "\n";
  cout << "  Evaluate polynomials at X = " << x << "\n";
  cout << "  We use the parameter S = " << s << "\n";
  cout << "\n";
  cout << "  Order       Horner          Direct\n";
  cout << "\n";

  cx2 = cardan ( n, x, s );

  for ( n = 0; n <= N_MAX; n++ )
  {
    cardan_poly_coef ( n, s, c );

    cx1 = r8poly_value ( n + 1, c, x );

    cout << "  " << setw(2)  << n
         << "  " << setw(14) << cx1
         << "  " << setw(14) << cx2[n] << "\n";
  }
  delete [] cx2;

  return;
# undef N_MAX
}
//****************************************************************************80

void test016 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST016 tests CATALAN and CATALAN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int *c2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST016\n";
  cout << "  CATALAN computes Catalan numbers.\n";
  cout << "  CATALAN_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N  exact C(I)  computed C(I)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = new int[n+1];

    catalan ( n, c2 );

    cout << "  " << setw(4) << n
         << "  " << setw(8) << c
         << "  " << setw(8) << c2[n] << "\n";

    delete [] c2;
  }

  return;
}
//****************************************************************************80

void test017 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST017 tests CATALAN_ROW_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 10

  int c[N_MAX+1];
  int i;
  int n;
  bool next;

  cout << "\n";
  cout << "TEST017\n";
  cout << "  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.\n";
  cout << "\n";
  cout << "  First, compute row 7:\n";
  cout << "\n";

  next = false;
  n = 7;
  catalan_row_next ( next, n, c );

  cout << setw(4) << n << "  ";
  for ( i = 0; i <= n; i++ )
  {
    cout << setw(8) << c[i] << "  ";
  }
  cout << "\n";

  cout << "\n";
  cout << "  Now compute rows consecutively, one at a time:\n";
  cout << "\n";

  next = false;

  for ( n = 0; n <= N_MAX; n++ )
  {
    catalan_row_next ( next, n, c );
    next = true;

    cout << setw(4) << i << "  ";
    for ( i = 0; i <= n; i++ )
    {
      cout << setw(8) << c[i] << "  ";
    }
    cout << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test0175 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0175 tests CHARLIER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5
# define N 5

  double a;
  double a_test[TEST_NUM] = { 0.25, 0.5, 1.0, 2.0, 10.0 };
  int i;
  int j;
  int n;
  int test;
  double x;
  double value[N+1];

  cout << "\n";
  cout << "TEST0175:\n";
  cout << "  CHARLIER evaluates Charlier polynomials.\n";
  cout << "\n";
  cout << "       N      A         X        P(N,A,X)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = N;
    a = a_test[test];

    cout << "\n";

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      charlier ( n, a, x, value );

      cout << "\n";
      for ( i = 0; i <= 5; i++ )
      {

        cout << "  " << setw(6)  << i     
             << "  " << setw(8)  << a
             << "  " << setw(8)  << x
             << "  " << setw(14) << value[i] << "\n";
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test018 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST018 tests CHEBY_T_POLY and CHEBY_T_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "TEST018:\n";
  cout << "  CHEBY_T_POLY evaluates the Chebyshev T polynomial.\n";
  cout << "  CHEBY_T_POLY_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N      X        Exact F       T(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_t_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = cheby_t_poly ( 1, n, x_vec );

    cout << "  " << setw(8)  << n
         << "  " << setw(8)  << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test0185 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0185 tests CHEBY_T_POLY_ZERO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 4

  double *fx;
  int i;
  int n;
  double *z;

  cout << "\n";
  cout << "TEST0185:\n";
  cout << "  CHEBY_T_POLY_ZERO returns zeroes of T(N,X).\n";
  cout << "\n";
  cout << "       N      X        T(N,X)\n";
  cout << "\n";

  for ( n = 1; n <= N_MAX; n++ )
  {
    z = cheby_t_poly_zero ( n );
    fx = cheby_t_poly ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << n
           << "  " << setw(8) << z[i]
           << "  " << setw(14) << fx[i+n*n] << "\n";
    }
    cout << "\n";
    delete [] fx;
    delete [] z;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test019 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST019 tests CHEBY_T_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int j;
  int n = 5;

  cout << "\n";
  cout << "TEST019\n";
  cout << "  CHEBY_T_POLY_COEF determines the  polynomial coefficients\n";
  cout << "  of the Chebyshev polynomial T(n,x).\n";

  c = cheby_t_poly_coef ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    cout << "\n";
    cout << "  T(" << i << ",x)\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          cout << setw(14) << c[i+j*(n+1)] << "\n";;
        }
        else if ( j == 1 )
        {
          cout << setw(14) << c[i+j*(n+1)] << " * x\n";
        }
        else
        {
          cout << setw(14) << c[i+j*(n+1)] << " * x^" << j << "\n";
        }
      }
    }
  }
 
  delete [] c;

  return;
}
//****************************************************************************80

void test020 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST020 tests CHEBY_U_POLY and CHEBY_U_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST020:\n";
  cout << "  CHEBY_U_POLY evaluates the Chebyshev U polynomial.\n";
  cout << "  CHEBY_U_POLY_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N      X        Exact F       U(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_u_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    cheby_u_poly ( n, x, fx2 );

    cout << "  " << setw(8)  << n
         << "  " << setw(8)  << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2[n] << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test021 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST021 tests CHEBY_U_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  cout << "\n";
  cout << "TEST021\n";
  cout << "  CHEBY_U_POLY_COEF determines the polynomial coefficients\n";
  cout << "  of the Chebyshev polynomial U(n,x).\n";

  cheby_u_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    cout << "\n";
    cout << "  U(" << i << ",x)\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(N+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          cout << setw(14) << c[i+j*(N+1)] << "\n";
        }
        else if ( j == 1 )
        {
          cout << setw(14) << c[i+j*(N+1)] << " * x\n";
        }
        else
        {
          cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
        }
      }
    }
  }
 
  return;
# undef N
}
//****************************************************************************80

void test0215 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0215 tests CHEBY_U_POLY_ZERO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 4

  double fx[N_MAX+1];
  int i;
  int n;
  double *z;

  cout << "\n";
  cout << "TEST0215:\n";
  cout << "  CHEBY_U_POLY_ZERO returns zeroes of U(N)(X).\n";
  cout << "\n";
  cout << "       N      X        U(N)(X)\n";
  cout << "\n";

  for ( n = 1; n <= N_MAX; n++ )
  {
    z = cheby_u_poly_zero ( n );

    for ( i = 0; i < n; i++ )
    {
      cheby_u_poly ( n, z[i], fx );

      cout << "  " << setw(8) << n
           << "  " << setw(8) << z[i]
           << "  " << setw(14) << fx[n] << "\n";
    }
    cout << "\n";
    delete [] z;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test0216 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0216 tests CHEBYSHEV_DISCRETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5
# define N 5

  int i;
  int j;
  int m;
  int n;
  double x;
  double value[N+1];

  cout << "\n";
  cout << "TEST0216:\n";
  cout << "  CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials.\n";
  cout << "\n";
  cout << "       N      M         X        T(N,M,X)\n";

  m = 5;
  n = N;

  for ( j = 0; j <= 5; j++ )
  {
    x = ( double ) ( j ) / 2.0;

    chebyshev_discrete ( n, m, x, value );

    cout << "\n";
    for ( i = 0; i <= 5; i++ )
    {
      cout << "  " << setw(6)  << i     
           << "  " << setw(8)  << m
           << "  " << setw(8)  << x
           << "  " << setw(14) << value[i] << "\n";
    }
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test0217 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0217 tests COLLATZ_COUNT and COLLATZ_COUNT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 March 2006
//
//  Author:
//
//    John Burkardt
//
{
  int count;
  int count2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST0217:\n";
  cout << "  COLLATZ_COUNT(N) counts the length of the\n";
  cout << "  Collatz sequence beginning with N.\n";
  cout << "\n";
  cout << "       N       COUNT(N)     COUNT(N)\n";
  cout << "              (computed)    (table)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    collatz_count_values ( &n_data, &n, &count );

    if ( n_data == 0 )
    {
      break;
    }

    count2 = collatz_count ( n );

    cout << "  " << setw(8) << n
         << "  " << setw(8) << count
         << "  " << setw(8) << count2 << "\n";
  }

  return;
}
//****************************************************************************80

void test0218 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0218 tests COLLATZ_COUNT_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 April 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i_max;
  int j_max;
  int n;

  cout << "\n";
  cout << "TEST0218:\n";
  cout << "  COLLATZ_COUNT_MAX(N) returns the length of the\n";
  cout << "  longest Collatz sequence from 1 to N.\n";
  cout << "\n";
  cout << "         N     I_MAX     J_MAX\n";
  cout << "\n";

  n = 10;

  while ( n <= 100000 )
  {
    collatz_count_max ( n, &i_max, &j_max );

    cout << "  " << setw(8) << n
         << "  " << setw(8) << i_max
         << "  " << setw(8) << j_max << "\n";

    n = n * 10;
  }

  return;
}
//****************************************************************************80

void test024 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST024 tests COMB_ROW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int c[N+1];
  int i;
  int j;
  bool next;

  cout << "\n";
  cout << "TEST024\n";
  cout << "  COMB_ROW computes a row of Pascal's triangle.\n";
  cout << "\n";

  next = false;

  for ( i = 0; i <= N; i++ )
  {
    comb_row ( next, i, c );
    next = true;
    cout                 << "  "
         << setw(2) << i << "  ";
    for ( j = 0; j <= i; j++ )
    {
      cout << setw(5) << c[j];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test0243 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0243 tests COS_POWER_INT and COS_POWER_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST0243:\n";
  cout << "  COS_POWER_INT computes the integral of the N-th power\n";
  cout << "  of the cosine function.\n";
  cout << "  COS_POWER_INT_VALUES returns selected values.\n";
  cout << "\n";
  cout << "         A         B       N        Exact    Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cos_power_int_values ( n_data, a, b, n, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = cos_power_int ( a, b, n );

    cout                    << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << b   << "  "
         << setw(6)  << n   << "  "
         << setw(12) << fx  << "  "
         << setw(12) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void test0245 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0245 tests R8_FACTORIAL and R8_FACTORIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fn;
  double fn2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST0245:\n";
  cout << "  R8_FACTORIAL evaluates the factorial function.\n";
  cout << "  R8_FACTORIAL_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N       Exact F       R8_FACTORIAL(N)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    fn2 = r8_factorial ( n );

    cout                    << "  "
         << setw(4)  << n   << "  "
         << setw(14) << fn  << "  "
         << setw(14) << fn2 << "\n";

  }

  return;
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 tests ERROR_F and ERF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
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

  cout << "\n";
  cout << "TEST025:\n";
  cout << "  ERROR_F evaluates the error function.\n";
  cout << "  ERF_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     X      Exact F       ERF(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = error_f ( x );

    cout                    << "  "
         << setw(8)  << x   << "  "
         << setw(14) << fx  << "  "
         << setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test0255 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0255 tests ERROR_F_INVERSE and ERF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x1;
  double x2;

  cout << "\n";
  cout << "TEST0255:\n";
  cout << "  ERROR_F_INVERSE inverts the error function.\n";
  cout << "  ERF_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "    FX           X1           X2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x1, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = error_f_inverse ( fx );

    cout                   << "  "
         << setw(8)  << fx << "  "
         << setw(14) << x1 << "  "
         << setw(14) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void test026 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026 tests EULER_NUMBER and EULER_NUMBER_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c1;
  int c2[13];
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST026\n";
  cout << "  EULER_NUMBER computes Euler numbers.\n";
  cout << "  EULER_NUMBER_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N  exact   EULER_NUMBER\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    euler_number ( n, c2 );

    cout                      << "  "
         << setw(4)  << n     << "  "
         << setw(12) << c1    << "  "
         << setw(12) << c2[n] << "\n";

  }
 
  return;
}
//****************************************************************************80

void test0265 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0265 tests EULER_NUMBER2 and EULER_NUMBER_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c1;
  double c2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST0265\n";
  cout << "  EULER_NUMBER2 computes Euler numbers.\n";
  cout << "  EULER_NUMBER_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N  exact   EULER_NUMBER2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = euler_number2 ( n );

    cout                   << "  "
         << setw(4)  << n  << "  "
         << setw(12) << c1 << "  "
         << setw(14) << c2 << "\n";

  }
 
  return;
}
//****************************************************************************80

void test028 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST028 tests EULER_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double f;
  int i;
  int n = 15;
  double x;

  x = 0.5;
 
  cout << "\n";
  cout << "TEST028\n";
  cout << "  EULER_POLY evaluates Euler polynomials.\n";
  cout << "\n";
  cout << "  N         X              F(X)\n";
  cout << "\n";
   
  for ( i = 0; i <= n; i++ )
  {
    f = euler_poly ( i, x );

    cout                  << "  "
         << setw(2)  << i << "  "
         << setw(14) << x << "  "
         << setw(14) << f << "\n";
  }
 
  return;
}
//****************************************************************************80

void test027 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST027 tests EULERIAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int e[N*N];
  int i;
  int j;

  cout << "\n";
  cout << "TEST027\n";
  cout << "  EULERIAN evaluates Eulerian numbers.\n";
  cout << "\n";
 
  eulerian ( N, e );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << setw(6) << e[i+j*N] << "  ";
    }
    cout << "\n";
  }
 
  return;
# undef N
}
//****************************************************************************80

void test029 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST029 tests F_HOFSTADTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int f;
  int i;

  cout << "\n";
  cout << "TEST029\n";
  cout << "  F_HOFSTADTER evaluates Hofstadter's recursive\n";
  cout << "  F function.\n";
  cout << "\n";
  cout << "     N   F(N)\n";
  cout << "\n";

  for ( i = 0; i <= 30; i++ )
  {
    f = f_hofstadter ( i );

    cout                 << "  "
         << setw(6) << i << "  "
         << setw(6) << f << "\n";
  }

  return;
}
//****************************************************************************80

void test031 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST031 tests FIBONACCI_DIRECT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int f;
  int i;
  int n = 20;

  cout << "\n";
  cout << "TEST031\n";
  cout << "  FIBONACCI_DIRECT evalutes a Fibonacci number directly.\n";
  cout << "\n";
  
  for ( i = 1; i <= n; i++ )
  {
    f = fibonacci_direct ( i );

    cout                  << "  "
         << setw(6)  << i << "  "
         << setw(10) << f << "\n";
  }
 
  return;
}
//****************************************************************************80

void test032 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST032 tests FIBONACCI_FLOOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int f;
  int i;
  int n;

  cout << "\n";
  cout << "TEST032\n";
  cout << "  FIBONACCI_FLOOR computes the largest Fibonacci number\n";
  cout << "  less than or equal to a given positive integer.\n";
  cout << "\n";
  cout << "     N  Fibonacci  Index\n";
  cout << "\n";

  for ( n = 1; n <= 20; n++ )
  {
    fibonacci_floor ( n, &f, &i );

    cout                 << "  "
         << setw(6) << n << "  "
         << setw(6) << f << "  "
         << setw(6) << i << "\n";
  }
 
  return;
}
//****************************************************************************80

void test033 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST033 tests FIBONACCI_RECURSIVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int f[N];
  int i;

  cout << "\n";
  cout << "TEST033\n";
  cout << "  FIBONACCI_RECURSIVE computes the Fibonacci sequence.\n";
  cout << "\n";
 
  fibonacci_recursive ( N, f );
 
  for ( i = 1; i <= N; i++ )
  {
    cout                       << "  "
         << setw(6)  << i      << "  "
         << setw(10) << f[i-1] << "\n";
  }
 
  return;
# undef N
}
//****************************************************************************80

void test034 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST034 tests G_HOFSTADTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "TEST034\n";
  cout << "  G_HOFSTADTER evaluates Hofstadter's recursive\n";
  cout << "  G function.\n";
  cout << "\n";
  cout << "     N   G(N)\n";
  cout << "\n";

  for ( i = 0; i <= 30; i++ )
  {
    cout                                  << "  "
         << setw(6) << i                  << "  "
         << setw(6) << g_hofstadter ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void test036 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST036 tests R8_GAMMA_LOG and GAMMA_LOG_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
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

  cout << "\n";
  cout << "TEST036:\n";
  cout << "  R8_GAMMA_LOG evaluates the logarithm of the Gamma function.\n";
  cout << "  GAMMA_LOG_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     X       Exact F       R8_GAMMA_LOG(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_gamma_log ( x );

    cout                    << "  "
         << setw(8)  << x   << "  "
         << setw(10) << fx  << "  "
         << setw(10) << fx2 << "\n";

  }

  return;
}
//****************************************************************************80

void test0365 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0365 tests LGAMMA and GAMMA_LOG_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
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

  cout << "\n";
  cout << "TEST0365:\n";
  cout << "  LGAMMA is a C math library function which evaluates\n";
  cout << "  the logarithm of the Gamma function.\n";
  cout << "  GAMMA_LOG_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     X       Exact F       LGAMMA(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lgamma ( x );

    cout                    << "  "
         << setw(8)  << x   << "  "
         << setw(10) << fx  << "  "
         << setw(10) << fx2 << "\n";

  }

  return;
}
//****************************************************************************80

void test037 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST037 tests GEGENBAUER_POLY and GEGENBAUER_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *c;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST037\n";
  cout << "  GEGENBAUER_POLY evaluates the Gegenbauer polynomials.\n";
  cout << "  GEGENBAUER_POLY_VALUES returns some exact values of\n";
  cout << "  the Gegenbauer polynomials.\n";
  cout << "\n";
  cout << "        N       A       X       GPV      GEGENBAUER\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {

    gegenbauer_poly_values ( &n_data, &n, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    c = new double[n+1];

    gegenbauer_poly ( n, a, x, c );
    fx2 = c[n];

    cout                    << "  "
         << setw(6)  << n   << "  "
         << setw(10) << a   << "  "
         << setw(10) << x   << "  "
         << setw(14) << fx  << "  "
         << setw(14) << fx2 << "\n";

    delete [] c;
  }
 
  return;
}
//****************************************************************************80

void test038 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST038 tests GUD and GUD_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
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

  cout << "\n";
  cout << "TEST038:\n";
  cout << "  GUD evaluates the Gudermannian function.\n";
  cout << "  GUD_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     X      Exact F       GUD(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gud_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = gud ( x );

    cout                    << "  "
         << setw(10) << x   << "  "
         << setw(10) << fx  << "  "
         << setw(10) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test039 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST039 tests HAIL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "TEST039\n";
  cout << "  HAIL(I) computes the length of the hail sequence\n";
  cout << "  for I, also known as the 3*N+1 sequence.\n";
  cout << "\n";
  cout << "  I,  HAIL(I)\n";
  cout << "\n";

  for ( i = 1; i <= 20; i++ )
  {
    cout                          << "  "
         << setw(4) << i          << "  "
         << setw(6) << hail ( i ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void test040 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST040 tests H_HOFSTADTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "TEST040\n";
  cout << "  H_HOFSTADTER evaluates Hofstadter's recursive\n";
  cout << "  H function.\n";

  cout << "\n";
  cout << "     N   H(N)\n";
  cout << "\n";

  for ( i = 0; i <= 30; i++ )
  {
    cout                                  << "  "
         << setw(6) << i                  << "  "
         << setw(6) << h_hofstadter ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void test041 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST041 tests HERMITE_POLY and HERMITE_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST041:\n";
  cout << "  HERMITE_POLY evaluates the Hermite polynomial.\n";
  cout << "  HERMITE_POLY_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N      X        Exact F       H(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    hermite_poly ( n, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test042 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST042 tests HERMITE_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  cout << "\n";
  cout << "TEST042\n";
  cout << "  HERMITE_POLY_COEF determines Hermite polynomial coefficients.\n";

  hermite_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    cout << "\n";
    cout << "  H(" << i << ")\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        cout << setw(14) << c[i+j*(N+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
      }
    }
  }

  return;
# undef N
}
//****************************************************************************80

void test0425 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0425 tests R8_HYPER_2F1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << " TEST0425:\n";
  cout << "   R8_HYPER_2F1 evaluates the hypergeometric function 2F1.\n";
  cout << "\n";
  cout << "      A       B       C       X      ";
  cout << " 2F1                       2F1                     DIFF\n";
  cout << "                                     ";
  cout << "(tabulated)               (computed)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hyper_2f1_values ( &n_data, &a, &b, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_hyper_2f1 ( a, b, c, x );

    cout << "  " << setw(6)  << setprecision(2)  << a
         << "  " << setw(6)  << setprecision(2)  << b  
         << "  " << setw(6)  << setprecision(2)  << c  
         << "  " << setw(6)  << setprecision(2)  << x  
         << "  " << setw(24) << setprecision(16) << fx
         << "  " << setw(24) << setprecision(16) << fx2
         << "  " << setw(10) << setprecision(4)  << r8_abs ( fx - fx2 ) << "\n";
  }
  return;
}
//****************************************************************************80

void test023 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST023 tests I4_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int cnk;
  int k;
  int n;

  cout << "\n";
  cout << "TEST023\n";
  cout << "  I4_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "   N     K    CNK\n";
  cout << "\n";

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = i4_choose ( n, k );

      cout                   << "  "
           << setw(6) << n   << "  "
           << setw(6) << k   << "  "
           << setw(6) << cnk << "\n";
    }
  }

  return;
}

//****************************************************************************80

void test043 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST043 tests I4_FACTORIAL and I4_FACTORIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int fn2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST043:\n";
  cout << "  I4_FACTORIAL evaluates the factorial function.\n";
  cout << "  I4_FACTORIAL_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     X       Exact F       I4_FACTORIAL(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    fn2 = i4_factorial ( n );

    cout                    << "  "
         << setw(4)  << n   << "  "
         << setw(12) << fn  << "  "
         << setw(12) << fn2 << "\n";

  }

  return;
}
//****************************************************************************80

void test044 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST044 tests I4_FACTORIAL2 and I4_FACTORIAL2_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int fn2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST044:\n";
  cout << "  I4_FACTORIAL2 evaluates the double factorial function.\n";
  cout << "\n";
  cout << "   N   Exact  I4_FACTORIAL2(N)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial2_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    fn2 = i4_factorial2 ( n );

    cout                   << "  "
         << setw(4) << n   << "  "
         << setw(8) << fn  << "  "
         << setw(8) << fn2 << "\n";
  }

  return;
}
//****************************************************************************80

void test045 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST045 tests PARTITION_COUNT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST045:\n";
  cout << "  For the number of partitions of an integer,\n";
  cout << "  PARTITION_COUNT_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N       Exact F\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    partition_count_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                  << "  "
         << setw(10) << n << "  "
         << setw(10) << c << "\n";
  }

  return;
}
//****************************************************************************80

void test046 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST046 tests I4_PARTITION_DISTINCT_COUNT and PARTITION_DISTINCT_COUNT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int c2;
  int n;
  int n_data;
  int n_max = 20;

  cout << "\n";
  cout << "TEST046:\n";
  cout << "  For the number of partitions of an integer\n";
  cout << "  into distinct parts,\n";
  cout << "  I4_PARTITION_DISTINCT_COUNT computes any value.\n";
  cout << "  PARTITION_DISTINCT_COUNT_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N       Exact F    Q(N)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    partition_distinct_count_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    if ( n_max < n )
    {
      continue;
    }

    c2 = i4_partition_distinct_count ( n );

    cout                   << "  "
         << setw(10) << n  << "  "
         << setw(10) << c  << "  "
         << setw(10) << c2 << "\n";

  }

  return;
}
//****************************************************************************80

void test047 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST047 tests I4_POCHHAMMER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST047\n";
  cout << "  I4_POCHHAMMER evaluates the integer Pochhammer function.\n";
  cout << "\n";
  cout << "   I   J   I4_Pochhammer(I,J)\n";
  cout << "\n";

  i = 3;
  j = 3;
  k = i4_pochhammer ( i, j );

  cout                 << "  "
       << setw(4) << i << "  "
       << setw(4) << j << "  "
       << setw(4) << k << "\n";

  i = 3;
  j = 4;
  k = i4_pochhammer ( i, j );

  cout                 << "  "
       << setw(4) << i << "  "
       << setw(4) << j << "  "
       << setw(4) << k << "\n";

  i = 3;
  j = 5;
  k = i4_pochhammer ( i, j );

  cout                 << "  "
       << setw(4) << i << "  "
       << setw(4) << j << "  "
       << setw(4) << k << "\n";

  return;
}
//****************************************************************************80

void test048 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST048 tests I4_IS_TRIANGULAR, I4_TO_TRIANGLE and TRIANGLE_TO_I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i2;
  int j;
  int k;
  int k2;
  bool l;

  cout << "\n";
  cout << "TEST048\n";
  cout << "  I4_TO_TRIANGLE converts a linear index to a\n";
  cout << "  triangular one.\n";
  cout << "  TRIANGLE_TO_I4 converts a triangular index to a\n";
  cout << "  linear one.\n";
  cout << "  I4_IS_TRIANGULAR returns 0 or 1 depending on\n";
  cout << "  whether I is triangular.\n";
  cout << "\n";
  cout << "   I  =>   J   K  =>   I   0/1\n";
  cout << "\n";

  for ( i = 0; i <= 20; i++ )
  {
    i4_to_triangle ( i, &j, &k );

    i2 = triangle_to_i4 ( j, k );

    l = i4_is_triangular ( i );

    cout                  << "  "
         << setw(4) << i  << "  "
         << setw(4) << j  << "  "
         << setw(4) << k  << "  "
         << setw(4) << i2 << "  "
         << setw(1) << l  << "\n";
  }
 
  return;
}
//****************************************************************************80

void test049 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST049 tests JACOBI_POLY and JACOBI_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST049\n";
  cout << "  JACOBI_POLY evaluates the Jacobi polynomials.\n";
  cout << "  JACOBI_POLY_VALUES returns some exact values of\n";
  cout << "  the Jacobi polynomials.\n";
  cout << "\n";
  cout << "        N       A       B      X       JPV      JACOBI\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {

    jacobi_poly_values ( n_data, n, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    c = jacobi_poly ( n, a, b, x );
    fx2 = c[n];

    cout                    << "  "
         << setw(6)  << n   << "  "
         << setw(6)  << a   << "  "
         << setw(6)  << b   << "  "
         << setw(10) << x   << "  "
         << setw(14) << fx  << "  "
         << setw(14) << fx2 << "\n";

    delete [] c;
  }
 
  return;
}
//****************************************************************************80

void test050 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST050 tests JACOBI_SYMBOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_TEST 4

  int i;
  int p;
  int ptest[N_TEST] = { 3, 9, 10, 12 };
  int q;

  cout << "\n";
  cout << "TEST050\n";
  cout << "  JACOBI_SYMBOL computes the Jacobi symbol\n";
  cout << "  (Q/P), which records if Q is a quadratic\n";
  cout << "  residue modulo the number P.\n";

  for ( i = 0; i < N_TEST; i++ )
  {
    p = ptest[i];
    cout << "\n";
    cout << "Jacobi Symbols for P = " << p << "\n";
    cout << "\n";
    for ( q = 0; q <= p; q++ )
    {
      cout << "  " << setw(8) << p
           << "  " << setw(8) << q
           << "  " << setw(8) << jacobi_symbol ( q, p ) << "\n";
    }
  }

  return;
# undef N_TEST
}
//****************************************************************************80

void test0505 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0505 tests KRAWTCHOUK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 2
# define N 5

  int i;
  int j;
  int m;
  int n;
  double p;
  double p_test[TEST_NUM] = { 0.25, 0.5 };
  int test;
  double x;
  double value[N+1];

  cout << "\n";
  cout << "TEST0505:\n";
  cout << "  KRAWTCHOUK evaluates Krawtchouk polynomials.\n";
  cout << "\n";
  cout << "        N         P         X          M      K(N,P,X,M)\n";
  cout << "\n";

  m = 5;
  n = N;

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test[test];

    cout << "\n";

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      krawtchouk ( n, p, x, m, value );

      cout << "\n";
      for ( i = 0; i <= 5; i++ )
      {

        cout << "  " << setw(8)  << i     
             << "  " << setw(8)  << p
             << "  " << setw(8)  << x
             << "  " << setw(8)  << m
             << "  " << setw(14) << value[i] << "\n";
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test051 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST051 tests LAGUERRE_ASSOCIATED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6
# define N_TEST 6

  double c[N+1];
  int i;
  int j;
  int m;
  int m_test[N_TEST] = { 0, 0, 1, 2, 3, 1 };
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  cout << "\n";
  cout << "TEST051\n";
  cout << "  LAGUERRE_ASSOCIATED evaluates the associated Laguerre\n";
  cout << "  polynomials.\n";

  for ( i = 0; i < N_TEST; i++ )
  {
    m = m_test[i];
    x = x_test[i];

    cout << "\n";
    cout << "  Table of L(N,M)(X) for\n";
    cout << "\n";
    cout << "  N(max) = " << N << "\n";
    cout << "  M      = " << m << "\n";
    cout << "  X =      " << x << "\n";
    cout << "\n";
 
    laguerre_associated ( N, m, x, c );
 
    for ( j = 0; j <= N; j++ )
    {
      cout                     << "  "
           << setw(6)  << j    << "  "
           << setw(14) << c[j] << "\n";
    }
 
  }

  return;
# undef N
# undef N_TEST
}
//****************************************************************************80

void test052 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST052 tests GEN_LAGUERRE_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define N_TEST 6

  double alpha;
  double alpha_test[N_TEST] = { 0.0, 0.0, 0.1, 0.1, 0.5, 1.0 };
  double c[N+1];
  int i;
  int j;
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  cout << "\n";
  cout << "TEST052\n";
  cout << "  GEN_LAGUERRE_POLY evaluates the generalized Laguerre\n";
  cout << "  functions.\n";

  for ( i = 0; i < N_TEST; i++ )
  {

    x = x_test[i];
    alpha = alpha_test[i];

    cout << "\n";
    cout << "  Table of L(N,ALPHA)(X) for\n";
    cout << "\n";
    cout << "    N(max) = " << N << "\n";
    cout << "    ALPHA =  " << alpha << "\n";
    cout << "    X =      " << x << "\n";
    cout << "\n";
  
    gen_laguerre_poly ( N, alpha, x, c );
 
    for ( j = 0; j <= N; j++ )
    {
      cout                     << "  "
           << setw(6)  << j    << "  "
           << setw(14) << c[j] << "\n";
    }
  }
 
  return;
# undef N
# undef N_TEST
}
//****************************************************************************80

void test054 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST054 tests LAGUERRE_POLY and LAGUERRE_POLYNOMIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST054:\n";
  cout << "  LAGUERRE_POLY evaluates the Laguerre polynomial.\n";
  cout << "  LAGUERRE_POLYNOMIAL_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N      X        Exact F       L(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    laguerre_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    laguerre_poly ( n, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test055 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST055 tests LAGUERRE_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double c[(N+1)*(N+1)];
  double fact;
  int i;
  int j;

  cout << "\n";
  cout << "TEST055\n";
  cout << "  LAGUERRE_POLY_COEF determines Laguerre \n";
  cout << "  polynomial coefficients.\n";

  laguerre_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    cout << "\n";
    cout << "  L(" << i << ")\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        cout << setw(14) << c[i+j*(N+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
      }
    }
  }
 
  for ( i = 0; i <= N; i++ )
  {
    fact = r8_factorial ( i );
    cout << "\n";
    cout << "  Factorially scaled L(" << i << ")\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        cout << setw(14) << fact *c[i+j*(N+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << fact *c[i+j*(N+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << fact *c[i+j*(N+1)] << " * x^" << j << "\n";
      }
    }
  }
  return;
# undef N
}
//****************************************************************************80

void test057 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST057 tests LEGENDRE_POLY and LEGENDRE_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double fp2[N_MAX+1];
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST057:\n";
  cout << "  LEGENDRE_POLY evaluates the Legendre PN function.\n";
  cout << "  LEGENDRE_POLY_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N      X        Exact F       P(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_poly ( n, x, fx2, fp2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test058 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST058 tests LEGENDRE_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  cout << "\n";
  cout << "TEST58\n";
  cout << "  LEGENDRE_POLY_COEF determines the Legendre P \n";
  cout << "  polynomial coefficients.\n";

  legendre_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    cout << "\n";
    cout << "  P(" << i << ")\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        cout << setw(14) << c[i+j*(N+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
      }
    }
  }
 
  return;
# undef N
}
//****************************************************************************80

void test059 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST059 tests LEGENDRE_ASSOCIATED and LEGENDRE_ASSOCIATED_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  double fx2[N_MAX+1];
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST059:\n";
  cout << "  LEGENDRE_ASSOCIATED evaluates associated Legendre functions.\n";
  cout << "  LEGENDRE_ASSOCIATED_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "      N       M    X     Exact F     PNM(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_associated ( n, m, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << m      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test0595 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0595 tests LEGENDRE_ASSOCIATED_NORMALIZED and LEGENDRE_ASSOCIATED_NORMALIZED_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2010
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  double fx2[N_MAX+1];
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST0595:\n";
  cout << "  LEGENDRE_ASSOCIATED_NORMALIZED evaluates \n";
  cout << "  normalized associated Legendre functions.\n";
  cout << "  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "      N       M    X     Exact F     PNM(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_associated_normalized ( n, m, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << m      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test060 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST060 tests LEGENDRE_FUNCTION_Q and LEGENDRE_FUNCTION_Q_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST060:\n";
  cout << "  LEGENDRE_FUNCTION_Q evaluates the Legendre Q function.\n";
  cout << "  LEGENDRE_FUNCTION_Q_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N      X        Exact F       Q(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_function_q_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_function_q ( n, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test061 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST061 tests LEGENDRE_SYMBOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_TEST 4

  int i;
  int l;
  int p;
  int ptest[N_TEST] = { 7, 11, 13, 17 };
  int q;

  cout << "\n";
  cout << "TEST061\n";
  cout << "  LEGENDRE_SYMBOL computes the Legendre\n";
  cout << "  symbol (Q/P) which records whether Q is \n";
  cout << "  a quadratic residue modulo the prime P.\n";

  for ( i = 0; i < N_TEST; i++ )
  {
    p = ptest[i];
    cout << "\n";
    cout << "  Legendre Symbols for P = " << p << "\n";
    cout << "\n";
    for ( q = 0; q <= p; q++ )
    {
      cout                                        << "  "
           << setw(8) << p                        << "  "
           << setw(8) << q                        << "  "
           << setw(8) << legendre_symbol ( q, p ) << "\n";
    }
  }

  return;
# undef N_TEST
}
//****************************************************************************80

void test0615 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0615 tests LERCH and LERCH_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double fx2;
  int n_data;
  int s;
  double z;

  cout << "\n";
  cout << "TEST0615:\n";
  cout << "  LERCH evaluates the Lerch function.\n";
  cout << "  LERCH_VALUES returns some tabulated values.\n";
  cout << "\n";
  cout << "       Z       S       A         Lerch           Lerch\n";
  cout << "                             Tabulated        Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lerch_values ( &n_data, &z, &s, &a, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lerch ( z, s, a );

    cout                    << "  "
         << setw(8)  << z   << "  "
         << setw(4)  << s   << "  "
         << setw(8)  << a   << "  "
         << setw(14) << fx  << "  "
         << setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test062 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST062 tests LOCK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a[N+1];
  int i;

  cout << "\n";
  cout << "TEST062\n";
  cout << "  LOCK counts the combinations on a button lock.\n";
  cout << "\n";
  cout << "  I,  LOCK(I)\n";
  cout << "\n";

  lock ( N, a );

  for ( i = 0; i <= N; i++ )
  {
    cout                     << "  "
         << setw(4)  << i    << "  "
         << setw(10) << a[i] << "\n";
  }
 
  return;
# undef N
}
//****************************************************************************80

void test0623 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0623 tests MEIXNER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define TEST_NUM 3

  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.0 };
  double c;
  double c_test[TEST_NUM] = { 0.125, 0.25, 0.5 };
  int i;
  int j;
  int n;
  int test;
  double v[N+1];
  double x;

  cout << "\n";
  cout << " TEST0623:\n";
  cout << "  MEIXNER evaluates Meixner polynomials.\n";
  cout << "\n";
  cout << "       N      BETA         C         X        M(N,BETA,C,X)\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = N;
    beta = beta_test[test];
    c = c_test[test];

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      meixner ( n, beta, c, x, v );

      cout << "\n";

      for ( i = 0; i <= n; i++ )
      {
        cout << "  " << setw(8) << i
             << "  " << setw(8) << beta
             << "  " << setw(8) << c
             << "  " << setw(8) << x
             << "  " << setw(14) << v[i] << "\n";
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void test0625 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0625 tests MERTENS and MERTENS_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST0625\n";
  cout << "  MERTENS computes the Mertens function.\n";
  cout << "  MERTENS_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "      N   Exact   MERTENS(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
     mertens_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                              << "  "
         << setw(8)  << n             << "  "
         << setw(10) << c             << "  "
         << setw(10) << mertens ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void test063 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST063 tests MOEBIUS and MOEBIUS_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST063\n";
  cout << "  MOEBIUS computes the Moebius function.\n";
  cout << "  MOEBIUS_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "      N   Exact   MOEBIUS(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
     moebius_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                              << "  "
         << setw(8)  << n             << "  "
         << setw(10) << c             << "  "
         << setw(10) << moebius ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void test0635 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0635 tests MOTZKIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a[N+1];
  int i;

  cout << "\n";
  cout << "TEST0635\n";
  cout << "  MOTZKIN computes the Motzkin numbers A(0:N).\n";
  cout << "  A(N) counts the paths from (0,0) to (N,0).\n";
  cout << "\n";
  cout << "  I,  A(I)\n";
  cout << "\n";

  motzkin ( N, a );

  for ( i = 0; i <= N; i++ )
  {
    cout                     << "  "
         << setw(4)  << i    << "  "
         << setw(10) << a[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test064 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST064 tests OMEGA and OMEGA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST064\n";
  cout << "  OMEGA computes the OMEGA function.\n";
  cout << "  OMEGA_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "          N   Exact   OMEGA(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
    omega_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                            << "  "
         << setw(12) << n           << "  "
         << setw(10) << c           << "  "
         << setw(10) << omega ( n ) << "\n";

  }
 
  return;
}
//****************************************************************************80

void test065 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST065 tests PENTAGON_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "TEST065\n";
  cout << "  PENTAGON_NUM computes the pentagonal numbers.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout                                  << "  "
         << setw(4) << n                  << "  "
         << setw(6) << pentagon_num ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void test066 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST066 tests PHI and PHI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST066\n";
  cout << "  PHI computes the PHI function.\n";
  cout << "  PHI_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N   Exact   PHI(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
    phi_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                          << "  "
         << setw(4)  << n         << "  "
         << setw(10) << c         << "  "
         << setw(10) << phi ( n ) << "\n";

  }
 
  return;
}
//****************************************************************************80

void test0665 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0665 tests POLY_BERNOULLI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2006
//
//  Author:
//
//    John Burkardt
//
{
  int b;
  int k;
  int n;

  cout << "\n";
  cout << "TEST0665\n";
  cout << "  POLY_BERNOULLI computes the poly-Bernoulli numbers\n";
  cout << "  of negative index, B_n^(-k)\n";
  cout << "\n";
  cout << "   N   K    B_N^(-K)\n";
  cout << "\n";

  for ( k = 0; k <= 6; k++ )
  {
    cout << "\n";
    for ( n = 0; n <= 6; n++ )
    {
      b = poly_bernoulli ( n, k );

      cout << "  " << setw(2)  << n
           << "  " << setw(2)  << k
           << "  " << setw(12) << b << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test0667 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0667 tests POLY_COEF_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  int dim;
  int n;

  cout << "\n";
  cout << "TEST0667\n";
  cout << "  POLY_COEF_COUNT counts the number of coefficients\n";
  cout << "  in a polynomial of degree DEGREE and dimension DIM.\n";
  cout << "\n";
  cout << " Dimension    Degree     Count\n";

  for ( dim = 1; dim <= 10; dim = dim + 3 )
  {
    cout << "\n";
    for ( degree = 0; degree <= 5; degree++ )
    {
      cout << "  " << setw(8) << dim
           << "  " << setw(8) << degree
           << "  " << setw(8) << poly_coef_count ( dim, degree ) << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test067 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST067 tests PYRAMID_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "TEST067\n";
  cout << "  PYRAMID_NUM computes the pyramidal numbers.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout                                 << "  "
         << setw(4) << n                 << "  "
         << setw(6) << pyramid_num ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void test0675 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0675 tests R8_ACOSH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST0675\n";
  cout << "  R8_ACOSH computes the inverse hyperbolic cosine\n";
  cout << "  of a given value.\n";
  cout << "\n";
  cout << "         X  R8_ACOSH(X)      COSH(R8_ACOSH(X))\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    a = r8_acosh ( x );
    x2 = cosh ( a );

    cout                   << "  "
         << setw(10) << x  << "  "
         << setw(10) << a  << "  "
         << setw(10) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void test068 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST068 tests R8_FACTORIAL and R8_FACTORIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fn;
  int n_data;
  int n;

  cout << "\n";
  cout << "TEST068:\n";
  cout << "  R8_FACTORIAL evaluates the factorial function.\n";
  cout << "  R8_FACTORIAL_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N       Exact F       R8_FACTORIAL(N)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    cout                                  << "  "
         << setw(4)  << n                 << "  "
         << setw(14) << fn                << "  "
         << setw(14) << r8_factorial ( n ) << "\n";
  }

  return;
}
//****************************************************************************80

void test022 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST022 tests R8_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cnk;
  int k;
  int n;

  cout << "\n";
  cout << "TEST022\n";
  cout << "  R8_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "   N     K    CNK\n";
  cout << "\n";

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = r8_choose ( n, k );

      cout                   << "  "
           << setw(6) << n   << "  "
           << setw(6) << k   << "  "
           << setw(6) << cnk << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test0685 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0685 tests R8_FACTORIAL_LOG and R8_FACTORIAL_LOG_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fn;
  int n_data;
  int n;

  cout << "\n";
  cout << "TEST0685:\n";
  cout << "  R8_FACTORIAL_LOG evaluates the logarithm of the\n";
  cout << "  factorial function.\n";
  cout << "  R8_FACTORIAL_LOG_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N	   Exact F	 R8_FACTORIAL_LOG(N)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_log_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    cout                                      << "  "
         << setw(5)  << n                     << "  "
         << setw(14) << fn                    << "  "
         << setw(14) << r8_factorial_log ( n ) << "\n";

  }

  return;
}
//****************************************************************************80

void test06855 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06855 tests R8_GAMMA, GAMMA, and GAMMA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST06855:\n";
  cout << "  R8_GAMMA evaluates the Gamma function.\n";
  cout << "  GAMMA is a C MATH routine for the Gamma function.\n";
  cout << "  GAMMA_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "         X                  Gamma(X)         " 
       << "         Gamma(X)         " 
       << "         Gamma(X)\n";
  cout << "                         (Tabulated)         " 
       << "       (R8_GAMMA)         "
       << "       (GAMMA)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_gamma ( x );

    fx3 = gamma ( x );

    cout << "  " << setw(8)  << setprecision(2)  << x   
         << "  " << setw(24) << setprecision(16) << fx  
         << "  " << setw(24) << setprecision(16) << fx2 
         << "  " << setw(24) << setprecision(16) << fx3 << "\n";

  }

  return;
}
//****************************************************************************80

void test06856 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06856 tests R8_PSI and PSI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
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

  cout << "\n";
  cout << "TEST06856:\n";
  cout << "  R8_PSI evaluates the Psi function.\n";
  cout << "  PSI_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "         X                  Psi(X)           " 
       << "         Psi(X)          DIFF\n";
  cout << "                         (Tabulated)         " 
       << "       (R8_PSI)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_psi ( x );

    cout << "  " << setprecision(2) << setw(8)  << x   
         << "  " << setprecision(16) << setw(24) << fx  
         << "  " << setprecision(16) << setw(24) << fx2 
         << "  " << setprecision(4) << setw(10) << r8_abs ( fx - fx2 ) << "\n";

  }

  return;
}
//****************************************************************************80

void test069 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST069 tests SIGMA and SIGMA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST069\n";
  cout << "  SIGMA computes the SIGMA function.\n";
  cout << "  SIGMA_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N   Exact   SIGMA(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
    sigma_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                            << "  "
         << setw(4)  << n           << "  "
         << setw(10) << c           << "  "
         << setw(10) << sigma ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void test0695 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0695 tests SIN_POWER_INT and SIN_POWER_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST0695:\n";
  cout << "  SIN_POWER_INT computes the integral of the N-th power\n";
  cout << "  of the sine function.\n";
  cout << "  SIN_POWER_INT_VALUES returns selected values.\n";
  cout << "\n";
  cout << "         A         B       N        Exact    Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   sin_power_int_values ( &n_data, &a, &b, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = sin_power_int ( a, b, n );

    cout                    << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << b   << "  "
         << setw(6)  << n   << "  "
         << setw(12) << fx  << "  "
         << setw(12) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void test0696 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0696 tests SLICE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 August 2011
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_MAX 5
# define SLICE_MAX 8

  int dim_max = DIM_MAX;
  int dim_num;
  int p[DIM_MAX*SLICE_MAX];
  int piece_num;
  int slice_max = SLICE_MAX;
  int slice_num;

  cout << "\n";
  cout << "TEST0696:\n";
  cout << "  SLICE determines the maximum number of pieces created\n";
  cout << "  by SLICE_NUM slices in a DIM_NUM space.\n";

  for ( dim_num = 1; dim_num <= dim_max; dim_num++ )
  {
    for ( slice_num = 1; slice_num <= slice_max; slice_num++ )
    {
      piece_num = slice ( dim_num, slice_num );
      p[dim_num-1+(slice_num-1)*dim_max] = piece_num;
    }
  }

  i4mat_print ( dim_max, slice_max, p, "  Slice Array:" );

  return;
# undef DIM_MAX
# undef SLICE_MAX
}
//****************************************************************************80

void test0697 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0697 tests SPHERICAL_HARMONIC and SPHERICAL_HARMONIC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  double c[N_MAX+1];
  int l;
  int m;
  int n_data;
  double phi;
  double s[N_MAX+1];
  double theta;
  double yi;
  double yi2;
  double yr;
  double yr2;

  cout << "\n";
  cout << "TEST0697:\n";
  cout << "  SPHERICAL_HARMONIC evaluates spherical harmonic functions.\n";
  cout << "  SPHERICAL_HARMONIC_VALUES returns some exact values.\n";
  cout << "\n";
  cout << 
    "         N         M    THETA      PHI            YR            YI\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    spherical_harmonic_values ( &n_data, &l, &m, &theta, &phi, &yr, &yi );

    if ( n_data == 0 )
    {
      break;
    }

    spherical_harmonic ( l, m, theta, phi, c, s );

    yr2 = c[l];
    yi2 = s[l];

    cout                      << "  "
         << setw(8)  << l     << "  "
         << setw(8)  << m     << "  "
         << setw(8)  << theta << "  "
         << setw(8)  << phi   << "  "
         << setw(14) << yr    << "  "
         << setw(14) << yi    << "\n";

    cout                      << "  "
         << "        "        << "  "
         << "        "        << "  "
         << "        "        << "  "
         << "        "        << "  "
         << setw(14) << yr2   << "  "
         << setw(14) << yi2   << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test070 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST070 tests STIRLING1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M 8
# define N 8

  int i;
  int j;
  int *s1;

  cout << "\n";
  cout << "TEST070\n";
  cout << "  STIRLING1: Stirling numbers of first kind.\n";
  cout << "    Get rows 1 through " << M << "\n";
  cout << "\n";
 
  s1 = stirling1 ( M, N );
 
  for ( i = 0; i < M; i++ )
  {
    cout << setw(6) << i+1 << "  ";
    for ( j = 0; j < N; j++ )
    {
      cout << setw(6) << s1[i+j*M] << "  ";
    }
    cout << "\n";
  }

  delete [] s1;
 
  return;
# undef M
# undef N
}
//****************************************************************************80

void test071 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST071 tests STIRLING2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M 8
# define N 8

  int i;
  int j;
  int *s2;

  cout << "\n";
  cout << "TEST071\n";
  cout << "  STIRLING2: Stirling numbers of second kind.\n";
  cout << "  Get rows 1 through " << M << "\n";
  cout << "\n";
 
  s2 = stirling2 ( M, N );
 
  for ( i = 0; i < M; i++ )
  {
    cout << setw(6) << i+1 << "  ";
    for ( j = 0; j < N; j++ )
    {
      cout << setw(6) << s2[i+j*M] << "  ";
    }
    cout << "\n";
  }
 
  delete [] s2;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test072 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST072 tests TAU and TAU_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST072\n";
  cout << "  TAU computes the Tau function.\n";
  cout << "  TAU_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N  exact C(I)  computed C(I)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
    tau_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                          << "  "
         << setw(4)  << n         << "  "
         << setw(10) << c         << "  "
         << setw(10) << tau ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void test073 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST073 tests TETRAHEDRON_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "TEST073\n";
  cout << "  TETRAHEDRON_NUM computes the tetrahedron numbers.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout                                     << "  "
         << setw(4) << n                     << "  "
         << setw(6) << tetrahedron_num ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void test074 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST074 tests TRIANGLE_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "TEST074\n";
  cout << "  TRIANGLE_NUM computes the triangular numbers.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout                                  << "  "
         << setw(4) << n                  << "  "
         << setw(6) << triangle_num ( n ) << "\n";;
  }
 
  return;
}
//****************************************************************************80

void test075 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST075 tests V_HOFSTADTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int v;

  cout << "\n";
  cout << "TEST075\n";
  cout << "  V_HOFSTADTER evaluates Hofstadter's recursive\n";
  cout << "  V function.\n";
  cout << "\n";
  cout << "     N   V(N)\n";
  cout << "\n";

  for ( i = 0; i <= 30; i++ )
  {
    cout                                  << "  "
         << setw(6) << i                  << "  "
         << setw(6) << v_hofstadter ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void test076 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST076 tests VIBONACCI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 20
  int i;
  int j;
  int seed;
  int v1[N];
  int v2[N];
  int v3[N];

  cout << "\n";
  cout << "TEST076\n";
  cout << "  VIBONACCI computes a Vibonacci sequence.\n";
  cout << "\n";
  cout << "  We compute the series 3 times.\n";
  cout << "\n";
  cout << "     I      V1      V2      V3\n";
  cout << "\n";

  seed = 123456789;

  vibonacci ( N, &seed, v1 );
  vibonacci ( N, &seed, v2 );
  vibonacci ( N, &seed, v3 );

  for ( i = 0; i < N; i++ )
  {
    cout                     << "  "
         << setw(6) << i     << "  "
         << setw(6) << v1[i] << "  "
         << setw(6) << v2[i] << "  "
         << setw(6) << v3[i] << "\n";
  } 

  return;
# undef N
}
//****************************************************************************80

void test077 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST077 tests ZECKENDORF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M_MAX 20

  int i;
  int i_list[M_MAX];
  int j;
  int f_list[M_MAX];
  int f_sum;
  int m;
  int n;

  cout << "\n";
  cout << "TEST077\n";
  cout << "  ZECKENDORF computes the Zeckendorf decomposition of\n";
  cout << "  an integer N into nonconsecutive Fibonacci numbers.\n";
  cout << "\n";
  cout << "   N Sum M Parts\n";
  cout << "\n";

  for ( n = 1; n <= 100; n++ )
  {
    zeckendorf ( n, M_MAX, &m, i_list, f_list );

    cout << setw(4) << n << "  ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(4) << f_list[j] << "  ";
    }
    cout << "\n";

  }

  return;
# undef M_MAX
}
//****************************************************************************80

void test0773 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0773 tests ZERNIKE_POLY and ZERNIKE_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int m;
  int n;
  double rho;
  double z1;
  double z2;

  cout << "\n";
  cout << "TEST0773\n";
  cout << "  ZERNIKE_POLY_COEF returns the coefficients of a\n";
  cout << "  Zernike polynomial.\n";
  cout << "  ZERNIKE_POLY evaluates a Zernike polynomial directly.\n";
  cout << "\n";
  cout << "  Table of polynomial coefficients:\n";
  cout << "\n";
  cout << "   N   M\n";
  cout << "\n";

  for ( n = 0; n <= 5; n++ )
  {
    cout << "\n";
    for ( m = 0; m <= n; m++ )
    {
      c = zernike_poly_coef ( m, n );
      cout << "  " << setw(2) << n
           << "  " << setw(2) << m;
      for ( i = 0; i <= n; i++ )
      {
        cout << "  " << setw(7) << c[i];
      }
      cout << "\n";
      delete [] c;
    }
  }

  rho = 0.987654321;

  cout << "\n";
  cout << "  Z1: Compute polynomial coefficients,\n";
  cout << "  then evaluate by Horner's method;\n";
  cout << "  Z2: Evaluate directly by recursion.\n";
  cout << "\n";
  cout << "   N   M       Z1              Z2\n";
  cout << "\n";

  for ( n = 0; n <= 5; n++ )
  {
    cout << "\n";
    for ( m = 0; m <= n; m++ )
    {
      c = zernike_poly_coef ( m, n );
      z1 = r8poly_value ( n + 1, c, rho );

      z2 = zernike_poly ( m, n, rho );
      cout << "  " << setw(2)  << n
           << "  " << setw(2)  << m
           << "  " << setw(16) << z1
           << "  " << setw(16) << z2 << "\n";

      delete [] c;
    }
  }

  return;
}
//****************************************************************************80

void test0775 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0775 tests ZERNIKE_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int m;
  int n;

  cout << "\n";
  cout << "TEST0775\n";
  cout << "  ZERNIKE_POLY_COEF determines the Zernike\n";
  cout << "  polynomial coefficients.\n";

  n = 5;

  for ( m = 0; m <= n; m++ )
  {
    c = zernike_poly_coef ( m, n );
    r8poly_print ( n, c, "  Zernike polynomial" );
    delete [] c;
  }

  return;
}
//****************************************************************************80

void test078 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST078 tests ZETA and ZETA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int *c2;
  int n;
  int n_data;
  double n_real;
  double z1;
  double z2;

  cout << "\n";
  cout << "TEST078\n";
  cout << "  ZETA computes the Zeta function.\n";
  cout << "  ZETA_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "       N            exact Zeta         computed Zeta\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    zeta_values ( &n_data, &n, &z1 );

    if ( n_data == 0 )
    {
      break;
    }

    n_real = ( double ) n;

    z2 = zeta ( n_real );

    cout                   << "  "
         << setw(6)  << n  << "  "
         << setw(20) << z1 << "  "
         << setw(20) << z2 << "\n";

  }

  return;
}
