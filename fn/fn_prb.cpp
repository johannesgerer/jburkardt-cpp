# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "fn.hpp"
# include "test_values.hpp"

int main ( );
void acos_test ( );
void acosh_test ( );
void ai_test ( );
void aid_test ( );
void asin_test ( );
void asinh_test ( );
void atan_test ( );
void atan2_test ( );
void atanh_test ( );
void besi0_test ( );
void besi1_test ( );
void besj0_test ( );
void besj1_test ( );
void besk0_test ( );
void besk1_test ( );
void besy0_test ( );
void besy1_test ( );
void beta_test ( );
void betai_test ( );
void bi_test ( );
void bid_test ( );
void binom_test ( );
void cbrt_test ( );
void chi_test ( );
void chu_test ( );
void ci_test ( );
void cin_test ( );
void cinh_test ( );
void cos_test ( );
void cos_deg_test ( );
void cosh_test ( );
void cot_test ( );
void dawson_test ( );
void e1_test ( );
void ei_test ( );
void erf_test ( );
void erfc_test ( );
void exp_test ( );
void fac_test ( );
void gamma_test ( );
void gamma_inc_test ( );
void gamma_inc_tricomi_test ( );
void int_test ( );
void lbeta_test ( );
void li_test ( );
void lngam_test ( );
void log_test ( );
void log10_test ( );
void poch_test ( );
void psi_test ( );
void rand_test ( );
void shi_test ( );
void si_test ( );
void sin_test ( );
void sin_deg_test ( );
void sinh_test ( );
void spence_test ( );
void sqrt_test ( );
void tan_test ( );
void tanh_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FN_PRB.
//
//  Discussion:
//
//    FN_PRB tests the FN library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FN_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the FN library.\n";

  acos_test ( );
  acosh_test ( );
  ai_test ( );
  aid_test ( );
  asin_test ( );
  asinh_test ( );
  atan_test ( );
  atan2_test ( );
  atanh_test ( );
  besi0_test ( );
  besi1_test ( );
  besj0_test ( );
  besj1_test ( );
  besk0_test ( );
  besk1_test ( );
  besy0_test ( );
  besy1_test ( );
  beta_test ( );
  betai_test ( );
  bi_test ( );
  bid_test ( );
  binom_test ( );
  cbrt_test ( );
  chi_test ( );
  chu_test ( );
  ci_test ( );
  cin_test ( );
  cinh_test ( );
  cos_test ( );
  cos_deg_test ( );
  cosh_test ( );
  cot_test ( );
  dawson_test ( );
  e1_test ( );
  ei_test ( );
  erf_test ( );
  erfc_test ( );
  exp_test ( );
  fac_test ( );
  gamma_test ( );
  gamma_inc_test ( );
  gamma_inc_tricomi_test ( );
  int_test ( );
  lbeta_test ( );
  li_test ( );
  lngam_test ( );
  log_test ( );
  log10_test ( );
  poch_test ( );
  psi_test ( );
  rand_test ( );
  shi_test ( );
  si_test ( );
  sin_test ( );
  sin_deg_test ( );
  sinh_test ( );
  spence_test ( );
  sqrt_test ( );
  tan_test ( );
  tanh_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FN_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void acos_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ACOS_TEST tests R4_ACOS and R8_ACOS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "ACOS_TEST:\n";
  cout << "  Test ARCCOS_VALUES, R4_ACOS, R8_ACOS.\n";
  cout << "\n";
  cout << "             X      ARCCOS(X)\n";
  cout << "                   R4_ACOS(X)         Diff\n";
  cout << "                   R8_ACOS(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    arccos_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_acos ( ( float ) x );
    fx3 = r8_acos ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void acosh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ACOSH_TEST tests R4_ACOSH and R8_ACOSH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "ACOSH_TEST:\n";
  cout << "  Test ARCCOSH_VALUES, R4_ACOSH, R8_ACOSH\n";
  cout << "\n";
  cout << "             X      ARCCOSH(X)\n";
  cout << "                   R4_ACOSH(X)        Diff\n";
  cout << "                   R8_ACOSH(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    arccosh_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_acosh ( ( float ) x );
    fx3 = r8_acosh ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void ai_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AI_TEST tests R4_AI and R8_AI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "AI_TEST:\n";
  cout << "  Test AIRY_AI_VALUES, R4_AI, R8_AI.\n";
  cout << "\n";
  cout << "             X   AIRY_AI(X)\n";
  cout << "                   R4_AI(X)         Diff\n";
  cout << "                   R8_AI(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_ai ( ( float ) x );
    fx3 = r8_ai ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void aid_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AID_TEST tests R4_AID and R8_AID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "AID_TEST:\n";
  cout << "  Test AIRY_AI_PRIME_VALUES, R4_AID, R8_AID.\n";
  cout << "\n";
  cout << "             X   AIRY_AID(X)\n";
  cout << "                   R4_AID(X)         Diff\n";
  cout << "                   R8_AID(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_prime_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_aid ( ( float ) x );
    fx3 = r8_aid ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void asin_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ASIN_TEST tests R4_ASIN and R8_ASIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "ASIN_TEST:\n";
  cout << "  Test ARCSIN_VALUES, R4_ASIN, R8_ASIN.\n";
  cout << "\n";
  cout << "             X      ARCSIN(X)\n";
  cout << "                   R4_ASIN(X)         Diff\n";
  cout << "                   R8_ASIN(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    arcsin_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_asin ( ( float ) x );
    fx3 = r8_asin ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void asinh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ASINH_TEST tests R4_ASINH and R8_ASINH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "ASINH_TEST:\n";
  cout << "  Test ARCSINH_VALUES, R4_ASINH, R8_ASINH\n";
  cout << "\n";
  cout << "             X      ARCSINH(X)\n";
  cout << "                   R4_ASINH(X)        Diff\n";
  cout << "                   R8_ASINH(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    arcsinh_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_asinh ( ( float ) x );
    fx3 = r8_asinh ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void atan_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ATAN_TEST tests R4_ATAN and R8_ATAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "ATAN_TEST:\n";
  cout << "  Test ARCTAN_VALUES, R4_ATAN, R8_ATAN.\n";
  cout << "\n";
  cout << "             X      ARCTAN(X)\n";
  cout << "                   R4_ATAN(X)         Diff\n";
  cout << "                   R8_ATAN(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    arctan_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_atan ( ( float ) x );
    fx3 = r8_atan ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void atan2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ATAN2_TEST tests R4_ATAN2 and R8_ATAN2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;
  double y;

  cout << "\n";
  cout << "ATAN2_TEST:\n";
  cout << "  Test ARCTAN2_VALUES, R4_ATAN2, R8_ATAN2.\n";
  cout << "\n";
  cout << "             X             Y      ARCTAN2(Y,X)\n";
  cout << "                                 R4_ATAN2(Y,X)         Diff\n";
  cout << "                                 R8_ATAN2(Y,X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    arctan2_values ( n_data, x, y, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_atan2 ( ( float ) y, ( float ) x );
    fx3 = r8_atan2 ( y, x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << y
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void atanh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ATANH_TEST tests R4_ATANH and R8_ATANH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "ATANH_TEST:\n";
  cout << "  Test ARCTANH_VALUES, R4_ATANH, R8_ATANH\n";
  cout << "\n";
  cout << "             X      ARCTANH(X)\n";
  cout << "                   R4_ATANH(X)        Diff\n";
  cout << "                   R8_ATANH(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    arctanh_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_atanh ( ( float ) x );
    fx3 = r8_atanh ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void besi0_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESI0_TEST tests R4_BESI0 and R8_BESI0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BESI0_TEST:\n";
  cout << "  Test BESSEL_I0_VALUES, R4_BESI0, R8_BESI0\n";
  cout << "\n";
  cout << "             X      BESI0(X)\n";
  cout << "                   R4_BESI0(X)        Diff\n";
  cout << "                   R8_BESI0(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besi0 ( ( float ) x );
    fx3 = r8_besi0 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void besi1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESI1_TEST tests R4_BESI1 and R8_BESI1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BESI1_TEST:\n";
  cout << "  Test BESSEL_I1_VALUES, R4_BESI1, R8_BESI1\n";
  cout << "\n";
  cout << "             X      BESI1(X)\n";
  cout << "                   R4_BESI1(X)        Diff\n";
  cout << "                   R8_BESI1(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_i1_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besi1 ( ( float ) x );
    fx3 = r8_besi1 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void besj0_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESJ0_TEST tests R4_BESJ0 and R8_BESJ0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BESJ0_TEST:\n";
  cout << "  Test BESSEL_J0_VALUES, R4_BESJ0, R8_BESJ0\n";
  cout << "\n";
  cout << "             X      BESJ0(X)\n";
  cout << "                   R4_BESJ0(X)        Diff\n";
  cout << "                   R8_BESJ0(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besj0 ( ( float ) x );
    fx3 = r8_besj0 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void besj1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESJ1_TEST tests R4_BESJ1 and R8_BESJ1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BESJ1_TEST:\n";
  cout << "  Test BESSEL_J1_VALUES, R4_BESJ1, R8_BESJ1\n";
  cout << "\n";
  cout << "             X      BESJ1(X)\n";
  cout << "                   R4_BESJ1(X)        Diff\n";
  cout << "                   R8_BESJ1(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_j1_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besj1 ( ( float ) x );
    fx3 = r8_besj1 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void besk0_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESK0_TEST tests R4_BESK0 and R8_BESK0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BESK0_TEST:\n";
  cout << "  Test BESSEL_K0_VALUES, R4_BESK0, R8_BESK0\n";
  cout << "\n";
  cout << "             X      BESK0(X)\n";
  cout << "                   R4_BESK0(X)        Diff\n";
  cout << "                   R8_BESK0(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_k0_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besk0 ( ( float ) x );
    fx3 = r8_besk0 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void besk1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESK1_TEST tests R4_BESK1 and R8_BESK1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BESK1_TEST:\n";
  cout << "  Test BESSEL_K1_VALUES, R4_BESK1, R8_BESK1\n";
  cout << "\n";
  cout << "             X      BESK1(X)\n";
  cout << "                   R4_BESK1(X)        Diff\n";
  cout << "                   R8_BESK1(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_k1_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besk1 ( ( float ) x );
    fx3 = r8_besk1 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void besy0_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESY0_TEST tests R4_BESY0 and R8_BESY0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BESY0_TEST:\n";
  cout << "  Test BESSEL_Y0_VALUES, R4_BESY0, R8_BESY0\n";
  cout << "\n";
  cout << "             X      BESY0(X)\n";
  cout << "                   R4_BESY0(X)        Diff\n";
  cout << "                   R8_BESY0(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_y0_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besy0 ( ( float ) x );
    fx3 = r8_besy0 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void besy1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESY1_TEST tests R4_BESY1 and R8_BESY1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BESY1_TEST:\n";
  cout << "  Test BESSEL_Y1_VALUES, R4_BESY1, R8_BESY1\n";
  cout << "\n";
  cout << "             X      BESY1(X)\n";
  cout << "                   R4_BESY1(X)        Diff\n";
  cout << "                   R8_BESY1(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_y1_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_besy1 ( ( float ) x );
    fx3 = r8_besy1 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void beta_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_TEST tests R4_BETA and R8_BETA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;

  cout << "\n";
  cout << "BETA_TEST:\n";
  cout << "  Test BETA_VALUES, R4_BETA, R8_BETA.\n";
  cout << "\n";
  cout << "             X        BETA(A,B)\n";
  cout << "                   R4_BETA(A,B)       Diff\n";
  cout << "                   R8_BETA(A,B)       Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_values ( n_data, a, b, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_beta ( ( float ) a, ( float ) b );
    fx3 = r8_beta ( a, b );

    cout << "\n";
    cout << "  " << setw(14) << a
         << "  " << setw(14) << b
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void betai_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BETAI_TEST tests R4_BETAI and R8_BETAI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BETAI_TEST:\n";
  cout << "  Test BETA_INC_VALUES, R4_BETAI, R8_BETAI.\n";
  cout << "\n";
  cout << "             X        BETA(A,B,X)\n";
  cout << "                   R4_BETAI(A,B,X)       Diff\n";
  cout << "                   R8_BETAI(A,B,X)       Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( n_data, a, b, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_betai ( ( float ) x, ( float ) a, ( float ) b );
    fx3 = r8_betai ( x, a, b );

    cout << "\n";
    cout << "  " << setw(14) << a
         << "  " << setw(14) << b
         << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void bi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BI_TEST tests R4_BI and R8_BI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BI_TEST:\n";
  cout << "  Test AIRY_BI_VALUES, R4_BI, R8_BI.\n";
  cout << "\n";
  cout << "             X   AIRY_BI(X)\n";
  cout << "                   R4_BI(X)         Diff\n";
  cout << "                   R8_BI(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_bi ( ( float ) x );
    fx3 = r8_bi ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void bid_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BID_TEST tests R4_BID and R8_BID.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "BID_TEST:\n";
  cout << "  Test AIRY_BI_PRIME_VALUES, R4_BID, R8_BID.\n";
  cout << "\n";
  cout << "             X   AIRY_BID(X)\n";
  cout << "                   R4_BID(X)         Diff\n";
  cout << "                   R8_BID(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_prime_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_bid ( ( float ) x );
    fx3 = r8_bid ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void binom_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BINOM_TEST tests R4_BINOM and R8_BINOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  double diff;
  int fx1;
  float fx2;
  double fx3;
  int n_data;

  cout << "\n";
  cout << "BINOM_TEST:\n";
  cout << "  Test BINOM_VALUES, R4_BINOM, R8_BINOM.\n";
  cout << "\n";
  cout << "             X        BINOM(A,B)\n";
  cout << "                   R4_BINOM(A,B)       Diff\n";
  cout << "                   R8_BINOM(A,B)       Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    binomial_values ( n_data, a, b, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_binom ( ( float ) a, ( float ) b );
    fx3 = r8_binom ( ( double ) a, ( double ) b );

    cout << "\n";
    cout << "  " << setw(14) << a
         << "  " << setw(14) << b
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void cbrt_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CBRT_TEST tests R4_CBRT and R8_CBRT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "CBRT_TEST:\n";
  cout << "  Test CBRT_VALUES, R4_CBRT, R8_CBRT\n";
  cout << "\n";
  cout << "             X      CBRT(X)\n";
  cout << "                   R4_CBRT(X)        Diff\n";
  cout << "                   R8_CBRT(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    cbrt_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cbrt ( ( float ) x );
    fx3 = r8_cbrt ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void chi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_TEST tests R4_CHI and R8_CHI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "CHI_TEST:\n";
  cout << "  Test CHI_VALUES, R4_CHI, R8_CHI.\n";
  cout << "\n";
  cout << "             X      CHI(X)\n";
  cout << "                   R4_CHI(X)         Diff\n";
  cout << "                   R8_CHI(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_chi ( ( float ) x );
    fx3 = r8_chi ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void chu_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHU_TEST tests R4_CHU and R8_CHU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double diff;
  double fx1;
  double fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "CHU_TEST:\n";
  cout << "  Test HYPERGEOMETRIC_U_VALUES, R4_CHU, R8_CHU.\n";
  cout << "\n";
  cout << "             A               B               X     CHU(A,B,X)\n";
  cout << "                                                R4_CHU(A,B,X)";
  cout << "         Diff\n";
  cout << "                                                R8_CHU(A,B,X)";
  cout << "         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_u_values ( n_data, a, b, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_chu ( ( float ) a, ( float ) b, ( float ) x );
    fx3 = r8_chu ( a, b, x );

    cout << "\n";
    cout << "  " << setw(14) << a
         << "  " << setw(14) << b
         << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx2
         << "  " << setw(14) << r4_abs ( fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx3
         << "  " << setw(14) << r8_abs ( fx1 - fx3 ) << "\n";
  }

  return;
}
//****************************************************************************80

void ci_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CI_TEST tests R4_CI and R8_CI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "CI_TEST:\n";
  cout << "  Test CI_VALUES, R4_CI, R8_CI.\n";
  cout << "\n";
  cout << "             X      CI(X)\n";
  cout << "                   R4_CI(X)         Diff\n";
  cout << "                   R8_CI(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    ci_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_ci ( ( float ) x );
    fx3 = r8_ci ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void cin_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CIN_TEST tests R4_CIN and R8_CIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "CIN_TEST:\n";
  cout << "  Test CIN_VALUES, R4_CIN, R8_CIN.\n";
  cout << "\n";
  cout << "             X      CIN(X)\n";
  cout << "                   R4_CIN(X)         Diff\n";
  cout << "                   R8_CIN(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    cin_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cin ( ( float ) x );
    fx3 = r8_cin ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void cinh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CINH_TEST tests R4_CINH and R8_CINH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "CINH_TEST:\n";
  cout << "  Test CINH_VALUES, R4_CINH, R8_CINH.\n";
  cout << "\n";
  cout << "             X      CINH(X)\n";
  cout << "                   R4_CINH(X)         Diff\n";
  cout << "                   R8_CINH(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    cinh_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cinh ( ( float ) x );
    fx3 = r8_cinh ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void cos_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COS_TEST tests R4_COS and R8_COS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "COS_TEST:\n";
  cout << "  Test COS_VALUES, R4_COS, R8_COS.\n";
  cout << "\n";
  cout << "             X      COS(X)\n";
  cout << "                   R4_COS(X)         Diff\n";
  cout << "                   R8_COS(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    cos_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cos ( ( float ) x );
    fx3 = r8_cos ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void cos_deg_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COS_DEG_TEST tests R4_COS_DEG and R8_COS_DEG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "COS_DEG_TEST:\n";
  cout << "  Test COS_DEGREE_VALUES, R4_COS_DEG, R8_COS_DEG.\n";
  cout << "\n";
  cout << "             X      COS_DEG(X)\n";
  cout << "                   R4_COS_DEG(X)         Diff\n";
  cout << "                   R8_COS_DEG(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    cos_degree_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cos_deg ( ( float ) x );
    fx3 = r8_cos_deg ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void cosh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COSH_TEST tests R4_COSH and R8_COSH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "COSH_TEST:\n";
  cout << "  Test COSH_VALUES, R4_COSH, R8_COSH\n";
  cout << "\n";
  cout << "             X      COSH(X)\n";
  cout << "                   R4_COSH(X)        Diff\n";
  cout << "                   R8_COSH(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    cosh_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cosh ( ( float ) x );
    fx3 = r8_cosh ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void cot_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COT_TEST tests R4_COT and R8_COT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "COT_TEST:\n";
  cout << "  Test COT_VALUES, R4_COT, R8_COT.\n";
  cout << "\n";
  cout << "             X      COT(X)\n";
  cout << "                   R4_COT(X)         Diff\n";
  cout << "                   R8_COT(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    cot_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_cot ( ( float ) x );
    fx3 = r8_cot ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void dawson_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DAWSON_TEST tests R4_DAWSON and R8_DAWSON.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "DAWSON_TEST:\n";
  cout << "  Test DAWSON_VALUES, R4_DAWSON, R8_DAWSON.\n";
  cout << "\n";
  cout << "             X      DAWSON(X)\n";
  cout << "                   R4_DAWSON(X)         Diff\n";
  cout << "                   R8_DAWSON(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    dawson_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_dawson ( ( float ) x );
    fx3 = r8_dawson ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void e1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    E1_TEST tests R4_E1 and R8_E1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "E1_TEST:\n";
  cout << "  Test E1_VALUES, R4_E1, R8_E1.\n";
  cout << "\n";
  cout << "             X      E1(X)\n";
  cout << "                   R4_E1(X)         Diff\n";
  cout << "                   R8_E1(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    e1_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_e1 ( ( float ) x );
    fx3 = r8_e1 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void ei_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EI_TEST tests R4_EI and R8_EI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "EI_TEST:\n";
  cout << "  Test EI_VALUES, R4_EI, R8_EI.\n";
  cout << "\n";
  cout << "             X      EI(X)\n";
  cout << "                   R4_EI(X)         Diff\n";
  cout << "                   R8_EI(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    ei_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_ei ( ( float ) x );
    fx3 = r8_ei ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void erf_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ERF_TEST tests R4_ERF and R8_ERF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "ERF_TEST:\n";
  cout << "  Test ERF_VALUES, R4_ERF, R8_ERF.\n";
  cout << "\n";
  cout << "             X      ERF(X)\n";
  cout << "                   R4_ERF(X)         Diff\n";
  cout << "                   R8_ERF(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_erf ( ( float ) x );
    fx3 = r8_erf ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void erfc_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ERFC_TEST tests R4_ERFC and R8_ERFC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "ERFC_TEST:\n";
  cout << "  Test ERFC_VALUES, R4_ERFC, R8_ERFC.\n";
  cout << "\n";
  cout << "             X      ERFC(X)\n";
  cout << "                   R4_ERFC(X)         Diff\n";
  cout << "                   R8_ERFC(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    erfc_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_erfc ( ( float ) x );
    fx3 = r8_erfc ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void exp_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EXP_TEST tests R4_EXP and R8_EXP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "EXP_TEST:\n";
  cout << "  Test EXP_VALUES, R4_EXP, R8_EXP.\n";
  cout << "\n";
  cout << "             X      EXP(X)\n";
  cout << "                   R4_EXP(X)         Diff\n";
  cout << "                   R8_EXP(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    exp_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_exp ( ( float ) x );
    fx3 = r8_exp ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void fac_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FAC_TEST tests R4_FAC and R8_FAC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  int fx1;
  float fx2;
  double fx3;
  int n;
  int n_data;

  cout << "\n";
  cout << "FAC_TEST:\n";
  cout << "  Test FACTORIAL_VALUES, R4_FAC, R8_FAC.\n";
  cout << "\n";
  cout << "             N      FAC(N)\n";
  cout << "                   R4_FAC(N)         Diff\n";
  cout << "                   R8_FAC(N)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    factorial_values ( n_data, n, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_fac ( n );
    fx3 = r8_fac ( n );

    cout << "\n";
    cout << "  " << setw(14) << n
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_TEST tests R4_GAMMA and R8_GAMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "GAMMA_TEST:\n";
  cout << "  Test GAMMA_VALUES, R4_GAMMA, R8_GAMMA\n";
  cout << "\n";
  cout << "             X      GAMMA(X)\n";
  cout << "                   R4_GAMMA(X)        Diff\n";
  cout << "                   R8_GAMMA(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_gamma ( ( float ) x );
    fx3 = r8_gamma ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_inc_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_TEST tests R4_GAMIC and R8_GAMIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "GAMMA_INC_TEST:\n";
  cout << "  Test GAMMA_INC_VALUES, R4_GAMIC, R8_GAMIC.\n";
  cout << "\n";
  cout << "             X        GAMIC(A,X)\n";
  cout << "                   R4_GAMIC(A,X)       Diff\n";
  cout << "                   R8_GAMIC(A,X)       Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( n_data, a, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_gamic ( ( float ) a, ( float ) x );
    fx3 = r8_gamic ( a, x );

    cout << "\n";
    cout << "  " << setw(14) << a
         << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_inc_tricomi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_TRICOMI_TEST tests R4_GAMIT and R8_GAMIT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "GAMMA_INC_TRICOMI_TEST:\n";
  cout << "  Test GAMMA_INC_TRICOMI_VALUES, R4_GAMIT, R8_GAMIT.\n";
  cout << "\n";
  cout << "             X        GAMIT(A,X)\n";
  cout << "                   R4_GAMIT(A,X)       Diff\n";
  cout << "                   R8_GAMIT(A,X)       Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_tricomi_values ( n_data, a, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_gamit ( ( float ) a, ( float ) x );
    fx3 = r8_gamit ( a, x );

    cout << "\n";
    cout << "  " << setw(14) << a
         << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void int_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    INT_TEST tests R4_INT and R8_INT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "INT_TEST:\n";
  cout << "  Test INT_VALUES, R4_INT, R8_INT\n";
  cout << "\n";
  cout << "             X      INT(X)\n";
  cout << "                   R4_INT(X)         Diff\n";
  cout << "                   R8_INT(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    int_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_int ( ( float ) x );
    fx3 = r8_int ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void lbeta_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LBETA_TEST tests R4_LBETA and R8_LBETA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;

  cout << "\n";
  cout << "LBETA_TEST:\n";
  cout << "  Test BETA_LOG_VALUES, R4_LBETA, R8_LBETA.\n";
  cout << "\n";
  cout << "             X        LBETA(A,B)\n";
  cout << "                   R4_LBETA(A,B)       Diff\n";
  cout << "                   R8_LBETA(A,B)       Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_log_values ( n_data, a, b, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_lbeta ( ( float ) a, ( float ) b );
    fx3 = r8_lbeta ( a, b );

    cout << "\n";
    cout << "  " << setw(14) << a
         << "  " << setw(14) << b
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void li_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LI_TEST tests R4_LI and R8_LI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "LI_TEST:\n";
  cout << "  Test LOGARITHMIC_INTEGRAL_VALUES, R4_LI, R8_LI\n";
  cout << "\n";
  cout << "             X      LI(X)\n";
  cout << "                   R4_LI(X)        Diff\n";
  cout << "                   R8_LI(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    logarithmic_integral_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_li ( ( float ) x );
    fx3 = r8_li ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void lngam_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LNGAM_TEST tests R4_LNGAM and R8_LNGAM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "LNGAM_TEST:\n";
  cout << "  Test GAMMA_LOG_VALUES, R4_LNGAM, R8_LNGAM\n";
  cout << "\n";
  cout << "             X        LNGAM(X)\n";
  cout << "                   R4_LNGAM(X)        Diff\n";
  cout << "                   R8_LNGAM(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_lngam ( ( float ) x );
    fx3 = r8_lngam ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void log_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_TEST tests R4_LOG and R8_LOG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL logcense.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "LOG_TEST:\n";
  cout << "  Test LOG_VALUES, R4_LOG, R8_LOG\n";
  cout << "\n";
  cout << "             X      LOG(X)\n";
  cout << "                   R4_LOG(X)        Diff\n";
  cout << "                   R8_LOG(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    log_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_log ( ( float ) x );
    fx3 = r8_log ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void log10_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOG10_TEST tests R4_LOG10 and R8_LOG10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL log10cense.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "LOG10_TEST:\n";
  cout << "  Test LOG10_VALUES, R4_LOG10, R8_LOG10\n";
  cout << "\n";
  cout << "             X      LOG10(X)\n";
  cout << "                   R4_LOG10(X)        Diff\n";
  cout << "                   R8_LOG10(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    log10_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_log10 ( ( float ) x );
    fx3 = r8_log10 ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void poch_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POCH_TEST tests R4_POCH and R8_POCH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "POCH_TEST:\n";
  cout << "  Test POCHHAMMER_VALUES, R4_POCH, R8_POCH.\n";
  cout << "\n";
  cout << "             X        POCH(A,X)\n";
  cout << "                   R4_POCH(A,X)       Diff\n";
  cout << "                   R8_POCH(A,X)       Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    pochhammer_values ( n_data, a, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_poch ( ( float ) a, ( float ) x );
    fx3 = r8_poch ( a, x );

    cout << "\n";
    cout << "  " << setw(14) << a
         << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void psi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PSI_TEST tests R4_PSI and R8_PSI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL psicense.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "PSI_TEST:\n";
  cout << "  Test PSI_VALUES, R4_PSI, R8_PSI\n";
  cout << "\n";
  cout << "             X      PSI(X)\n";
  cout << "                   R4_PSI(X)        Diff\n";
  cout << "                   R8_PSI(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_psi ( ( float ) x );
    fx3 = r8_psi ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void rand_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RAND_TEST tests R4_RAND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  float average;
  int i;
  int i_value[7] = { 1, 2, 3, 4, 10, 100, 1000 };
  int k;
  float r;
  float r_value[7] = { 
    0.0004127026, 
    0.6750836372, 
    0.1614754200, 
    0.9086198807, 
    0.5527787209, 
    0.3600893021, 
    0.2176990509 };
  float variance;

  cout << "\n";
  cout << "RAND_TEST:\n";
  cout << "  Test R4_RAND.\n";
  cout << "\n";
  cout << "               I       R4_RAND        Expected\n";
  cout << "\n";

  k = 0;

  for ( i = 1; i <= 1000; i++ )
  {
    r = r4_rand ( 0.0 );

    if ( i == i_value[k] )
    {
      cout << "  " << setw(14) << i
           << "  " << setw(14) << r
           << "  " << setw(14) << r_value[k] << "\n";
      k = k + 1;
    }
  }

  average = 0.0;
  for ( i = 1; i <= 1000000; i++ )
  {
    r = r4_rand ( 0.0 );
    average = average + r;
  }
  average = average / 1000000.0;
  cout << "\n";
  cout << "     Average =  "
       << "  " << setw(14) << average
       << "  " << setw(14) << 0.5 << "\n";

  variance = 0.0;
  for ( i = 1; i <= 1000000; i++ )
  {
    r = r4_rand ( 0.0 );
    variance = variance + pow ( r - average, 2 );
  }
  variance = variance / 1000000.0;
  cout << "     Variance = "
       << "  " << setw(14) << variance
       << "  " << setw(14) << 1.0 / 12.0 << "\n";

  return;
}
//****************************************************************************80

void shi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SHI_TEST tests R4_SHI and R8_SHI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "SHI_TEST:\n";
  cout << "  Test SHI_VALUES, R4_SHI, R8_SHI.\n";
  cout << "\n";
  cout << "             X      SHI(X)\n";
  cout << "                   R4_SHI(X)         Diff\n";
  cout << "                   R8_SHI(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    shi_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_shi ( ( float ) x );
    fx3 = r8_shi ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void si_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SI_TEST tests R4_SI and R8_SI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "SI_TEST:\n";
  cout << "  Test SI_VALUES, R4_SI, R8_SI.\n";
  cout << "\n";
  cout << "             X      SI(X)\n";
  cout << "                   R4_SI(X)         Diff\n";
  cout << "                   R8_SI(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    si_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_si ( ( float ) x );
    fx3 = r8_si ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void sin_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_TEST tests R4_SIN and R8_SIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "SIN_TEST:\n";
  cout << "  Test SIN_VALUES, R4_SIN, R8_SIN.\n";
  cout << "\n";
  cout << "             X      SIN(X)\n";
  cout << "                   R4_SIN(X)         Diff\n";
  cout << "                   R8_SIN(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    sin_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_sin ( ( float ) x );
    fx3 = r8_sin ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void sin_deg_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_DEG_TEST tests R4_SIN_DEG and R8_SIN_DEG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "SIN_DEG_TEST:\n";
  cout << "  Test SIN_DEGREE_VALUES, R4_SIN_DEG, R8_SIN_DEG.\n";
  cout << "\n";
  cout << "             X      SIN_DEG(X)\n";
  cout << "                   R4_SIN_DEG(X)         Diff\n";
  cout << "                   R8_SIN_DEG(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    sin_degree_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_sin_deg ( ( float ) x );
    fx3 = r8_sin_deg ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void sinh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SINH_TEST tests R4_SINH and R8_SINH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "SINH_TEST:\n";
  cout << "  Test SINH_VALUES, R4_SINH, R8_SINH\n";
  cout << "\n";
  cout << "             X      SINH(X)\n";
  cout << "                   R4_SINH(X)        Diff\n";
  cout << "                   R8_SINH(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    sinh_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_sinh ( ( float ) x );
    fx3 = r8_sinh ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void spence_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPENCE_TEST tests R4_SPENCE and R8_SPENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "SPENCE_TEST:\n";
  cout << "  Test DILOGARITHM_VALUES, R4_SPENCE, R8_SPENCE\n";
  cout << "\n";
  cout << "             X      SPENCE(X)\n";
  cout << "                   R4_SPENCE(X)        Diff\n";
  cout << "                   R8_SPENCE(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    dilogarithm_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_spence ( ( float ) x );
    fx3 = r8_spence ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void sqrt_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SQRT_TEST tests R4_SQRT and R8_SQRT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "SQRT_TEST:\n";
  cout << "  Test SQRT_VALUES, R4_SQRT, R8_SQRT\n";
  cout << "\n";
  cout << "             X      SQRT(X)\n";
  cout << "                   R4_SQRT(X)        Diff\n";
  cout << "                   R8_SQRT(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    sqrt_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_sqrt ( ( float ) x );
    fx3 = r8_sqrt ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void tan_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TAN_TEST tests R4_TAN and R8_TAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "TAN_TEST:\n";
  cout << "  Test TAN_VALUES, R4_TAN, R8_TAN.\n";
  cout << "\n";
  cout << "             X      TAN(X)\n";
  cout << "                   R4_TAN(X)         Diff\n";
  cout << "                   R8_TAN(X)         Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    tan_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_tan ( ( float ) x );
    fx3 = r8_tan ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
//****************************************************************************80

void tanh_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_TEST tests R4_TANH and R8_TANH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  double fx1;
  float fx2;
  double fx3;
  int n_data;
  double x;

  cout << "\n";
  cout << "TANH_TEST:\n";
  cout << "  Test TANH_VALUES, R4_TANH, R8_TANH\n";
  cout << "\n";
  cout << "             X      TANH(X)\n";
  cout << "                   R4_TANH(X)        Diff\n";
  cout << "                   R8_TANH(X)        Diff\n";

  n_data = 0;

  for ( ; ; )
  {
    tanh_values ( n_data, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r4_tanh ( ( float ) x );
    fx3 = r8_tanh ( x );

    cout << "\n";
    cout << "  " << setw(14) << x
         << "  " << setw(14) << fx1 << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx2 
         << "  " << setw(14) << r4_abs ( ( float ) fx1 - fx2 ) << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << fx3 
         << "  " << setw(14) << r8_abs (           fx1 - fx3 ) << "\n";
  }
  return;
}
