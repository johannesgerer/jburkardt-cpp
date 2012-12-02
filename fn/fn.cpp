# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <complex>
# include <cmath>

using namespace std;

# include "fn.hpp"

//****************************************************************************80

complex <float> c4_cos ( complex <float> z )

//****************************************************************************80
//
//  Purpose:
//
//    C4_COS evaluates the cosine of a C4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, complex <float> Z, the argument.
//
//    Output, complex <float> C4_COS, the cosine of Z.
//
{
  float cs;
  complex <float> value;
  float x;
  float y;

  x = real ( z );
  y = imag ( z );

  cs = cos ( x );

  value = complex <float> ( cs * cosh ( y ), - sin ( x ) * sinh ( y ) );

  return value;
}
//****************************************************************************80

complex <float> c4_sin ( complex <float> z )

//****************************************************************************80
//
//  Purpose:
//
//    C4_SIN evaluates the sine of a C4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, complex <float> Z, the argument.
//
//    Output, complex <float> C4_SIN, the sine of Z.
//
{
  float sn;
  complex <float> value;
  float x;
  float y;

  x = real ( z );
  y = imag ( z );

  sn = sin ( x );

  value = complex <float> ( sn * cosh ( y ), cos ( x ) * sinh ( y ) );

  return value;
}
//****************************************************************************80

int i4_abs ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_ABS returns the absolute value of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer.
//
//    Output, int I4_ABS, the absolute value of the integer.
//
{
  int value;

  if ( 0 <= i )
  {
    value = i;
  }
  else
  {
    value = - i;
  }
  return value;
}
//****************************************************************************80

int i4_mach ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MACH returns integer machine constants.
//
//  Discussion:
//
//    Input/output unit numbers.
//
//      I4_MACH(1) = the standard input unit.
//      I4_MACH(2) = the standard output unit.
//      I4_MACH(3) = the standard punch unit.
//      I4_MACH(4) = the standard error message unit.
//
//    Words.
//
//      I4_MACH(5) = the number of bits per integer storage unit.
//      I4_MACH(6) = the number of characters per integer storage unit.
//
//    Integers.
//
//    Assume integers are represented in the S digit base A form:
//
//      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
//
//    where 0 <= X(1:S-1) < A.
//
//      I4_MACH(7) = A, the base.
//      I4_MACH(8) = S, the number of base A digits.
//      I4_MACH(9) = A^S-1, the largest integer.
//
//    Floating point numbers
//
//    Assume floating point numbers are represented in the T digit 
//    base B form:
//
//      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
//
//    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
//
//      I4_MACH(10) = B, the base.
//
//    Single precision
//
//      I4_MACH(11) = T, the number of base B digits.
//      I4_MACH(12) = EMIN, the smallest exponent E.
//      I4_MACH(13) = EMAX, the largest exponent E.
//
//    Double precision
//
//      I4_MACH(14) = T, the number of base B digits.
//      I4_MACH(15) = EMIN, the smallest exponent E.
//      I4_MACH(16) = EMAX, the largest exponent E.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, Norman Schryer,
//    Algorithm 528,
//    Framework for a Portable Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, Number 2, June 1978, page 176-188.
//
//  Parameters:
//
//    Input, int I, chooses the parameter to be returned.
//    1 <= I <= 16.
//
//    Output, int I4_MACH, the value of the chosen parameter.
//
{
  int value;

  if ( i == 1 )
  {
    value = 5;
  }
  else if ( i == 2 )
  {
    value = 6;
  }
  else if ( i == 3 )
  {
    value = 7;
  }
  else if ( i == 4 )
  {
    value = 6;
  }
  else if ( i == 5 )
  {
    value = 32;
  }
  else if ( i == 6 )
  {
    value = 4;
  }
  else if ( i == 7 )
  {
    value = 2;
  }
  else if ( i == 8 )
  {
    value = 31;
  }
  else if ( i == 9 )
  {
    value = 2147483647;
  }
  else if ( i == 10 )
  {
    value = 2;
  }
  else if ( i == 11 )
  {
    value = 24;
  }
  else if ( i == 12 )
  {
    value = -125;
  }
  else if ( i == 13 )
  {
    value = 128;
  }
  else if ( i == 14 )
  {
    value = 53;
  }
  else if ( i == 15 )
  {
    value = -1021;
  }
  else if ( i == 16 )
  {
    value = 1024;
  }
  else
  {
    cerr << "\n";
    cerr << "I4_MACH - Fatal error!\n";
    cerr << "  The input argument I is out of bounds.\n";
    cerr << "  Legal values satisfy 1 <= I <= 16.\n";
    cerr << "  I = " << i << "\n";
    value = 0;
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_pow ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POW returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POW, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POW - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POW - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

float r4_acos ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ACOS evaluates the arc-cosine of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_ACOS, the arc-cosine of X.
//
{
  static float pi2 = 1.57079632679489661923;
  float value;

  value = pi2 - r4_asin ( x );

  return value;
}
//****************************************************************************80

float r4_acosh ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ACOSH evaluates the arc-hyperbolic cosine of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_ACOSH, the arc-hyperbolic cosine of X.
//
{
  static float aln2 = 0.69314718055994530942E+00;
  float value;
  static float xmax = 0.0;

  if ( xmax == 0.0E+00 )
  {
    xmax = 1.0E+00 / sqrt ( r4_mach ( 3 ) );
  }

  if ( x < 1.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_ACOSH - Fatal error!\n";
    cerr << "  X < 1.\n";
    exit ( 1 );
  }

  if ( x < xmax )
  {
    value = log ( x + sqrt ( x * x - 1.0E+00 ) );
  }
  else
  {
    value = aln2 + log ( x );
  }

  return value;
}
//****************************************************************************80

void r4_admp ( float x, float &ampl, float &phi )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ADMP: modulus and phase of the derivative of the Airy function.
//
//  Description:
//
//    This function must only be called when X <= -1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float &AMPL, &PHI, the modulus and phase of the 
//    derivative of the Airy function.
//
{
  static float an20cs[16] = {
     0.0126732217145738027E+00,
    -0.0005212847072615621E+00,
    -0.0000052672111140370E+00,
    -0.0000001628202185026E+00,
    -0.0000000090991442687E+00,
    -0.0000000007438647126E+00,
    -0.0000000000795494752E+00,
    -0.0000000000104050944E+00,
    -0.0000000000015932426E+00,
    -0.0000000000002770648E+00,
    -0.0000000000000535343E+00,
    -0.0000000000000113062E+00,
    -0.0000000000000025772E+00,
    -0.0000000000000006278E+00,
    -0.0000000000000001621E+00,
    -0.0000000000000000441E+00 };
  static float an21cs[24] = {
     0.0198313155263169394E+00,
    -0.0029376249067087533E+00,
    -0.0001136260695958196E+00,
    -0.0000100554451087156E+00,
    -0.0000013048787116563E+00,
    -0.0000002123881993151E+00,
    -0.0000000402270833384E+00,
    -0.0000000084996745953E+00,
    -0.0000000019514839426E+00,
    -0.0000000004783865344E+00,
    -0.0000000001236733992E+00,
    -0.0000000000334137486E+00,
    -0.0000000000093702824E+00,
    -0.0000000000027130128E+00,
    -0.0000000000008075954E+00,
    -0.0000000000002463214E+00,
    -0.0000000000000767656E+00,
    -0.0000000000000243883E+00,
    -0.0000000000000078831E+00,
    -0.0000000000000025882E+00,
    -0.0000000000000008619E+00,
    -0.0000000000000002908E+00,
    -0.0000000000000000993E+00,
    -0.0000000000000000343E+00 };
  static float an22cs[33] = {
     0.0537418629629794329E+00,
    -0.0126661435859883193E+00,
    -0.0011924334106593007E+00,
    -0.0002032327627275655E+00,
    -0.0000446468963075164E+00,
    -0.0000113359036053123E+00,
    -0.0000031641352378546E+00,
    -0.0000009446708886149E+00,
    -0.0000002966562236472E+00,
    -0.0000000969118892024E+00,
    -0.0000000326822538653E+00,
    -0.0000000113144618964E+00,
    -0.0000000040042691002E+00,
    -0.0000000014440333684E+00,
    -0.0000000005292853746E+00,
    -0.0000000001967763374E+00,
    -0.0000000000740800096E+00,
    -0.0000000000282016314E+00,
    -0.0000000000108440066E+00,
    -0.0000000000042074801E+00,
    -0.0000000000016459150E+00,
    -0.0000000000006486827E+00,
    -0.0000000000002574095E+00,
    -0.0000000000001027889E+00,
    -0.0000000000000412846E+00,
    -0.0000000000000166711E+00,
    -0.0000000000000067657E+00,
    -0.0000000000000027585E+00,
    -0.0000000000000011296E+00,
    -0.0000000000000004645E+00,
    -0.0000000000000001917E+00,
    -0.0000000000000000794E+00,
    -0.0000000000000000330E+00 };
  static float aph0cs[15] = {
    -0.0855849241130933257E+00,
     0.0011214378867065261E+00,
     0.0000042721029353664E+00,
     0.0000000817607381483E+00,
     0.0000000033907645000E+00,
     0.0000000002253264423E+00,
     0.0000000000206284209E+00,
     0.0000000000023858763E+00,
     0.0000000000003301618E+00,
     0.0000000000000527010E+00,
     0.0000000000000094555E+00,
     0.0000000000000018709E+00,
     0.0000000000000004024E+00,
     0.0000000000000000930E+00,
     0.0000000000000000229E+00 };
  static float aph1cs[22] = {
    -0.1024172908077571694E+00,
     0.0071697275146591248E+00,
     0.0001209959363122329E+00,
     0.0000073361512841220E+00,
     0.0000007535382954272E+00,
     0.0000001041478171741E+00,
     0.0000000174358728519E+00,
     0.0000000033399795033E+00,
     0.0000000007073075174E+00,
     0.0000000001619187515E+00,
     0.0000000000394539982E+00,
     0.0000000000101192282E+00,
     0.0000000000027092778E+00,
     0.0000000000007523806E+00,
     0.0000000000002156369E+00,
     0.0000000000000635283E+00,
     0.0000000000000191757E+00,
     0.0000000000000059143E+00,
     0.0000000000000018597E+00,
     0.0000000000000005950E+00,
     0.0000000000000001934E+00,
     0.0000000000000000638E+00 };
  static float aph2cs[32] = {
    -0.2057088719781465107E+00,
     0.0422196961357771922E+00,
     0.0020482560511207275E+00,
     0.0002607800735165006E+00,
     0.0000474824268004729E+00,
     0.0000105102756431612E+00,
     0.0000026353534014668E+00,
     0.0000007208824863499E+00,
     0.0000002103236664473E+00,
     0.0000000644975634555E+00,
     0.0000000205802377264E+00,
     0.0000000067836273921E+00,
     0.0000000022974015284E+00,
     0.0000000007961306765E+00,
     0.0000000002813860610E+00,
     0.0000000001011749057E+00,
     0.0000000000369306738E+00,
     0.0000000000136615066E+00,
     0.0000000000051142751E+00,
     0.0000000000019351689E+00,
     0.0000000000007393607E+00,
     0.0000000000002849792E+00,
     0.0000000000001107281E+00,
     0.0000000000000433412E+00,
     0.0000000000000170801E+00,
     0.0000000000000067733E+00,
     0.0000000000000027017E+00,
     0.0000000000000010835E+00,
     0.0000000000000004367E+00,
     0.0000000000000001769E+00,
     0.0000000000000000719E+00,
     0.0000000000000000294E+00 };
  float eta;
  static int nan20 = 0;
  static int nan21 = 0;
  static int nan22 = 0;
  static int naph0 = 0;
  static int naph1 = 0;
  static int naph2 = 0;
  static float pi34 = 2.3561944901923449E+00;
  float sqrtx;
  static float xsml = 0.0;
  float z;

  if ( nan20 == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    nan20 = r4_inits ( an20cs, 16, eta );
    nan21 = r4_inits ( an21cs, 24, eta );
    nan22 = r4_inits ( an22cs, 33, eta );
    naph0 = r4_inits ( aph0cs, 15, eta );
    naph1 = r4_inits ( aph1cs, 22, eta );
    naph2 = r4_inits ( aph2cs, 32, eta );
    xsml = - pow ( 128.0E+00 / r4_mach ( 3 ), 0.3333E+00 );
  }

  if ( x <= xsml )
  {
    z = 1.0E+00;
    ampl = 0.3125E+00 + r4_csevl ( z, an20cs, nan20 );
    phi = - 0.625E+00 + r4_csevl ( z, aph0cs, naph0 );
  }
  else if ( x < - 4.0E+00 )
  {
    z = 128.0E+00 / x / x / x + 1.0E+00;
    ampl = 0.3125E+00 + r4_csevl ( z, an20cs, nan20 );
    phi = - 0.625E+00 + r4_csevl ( z, aph0cs, naph0 );
  }
  else if ( x < - 2.0E+00 )
  {
    z = ( 128.0E+00 / x / x / x + 9.0E+00 ) / 7.0E+00;
    ampl = 0.3125E+00 + r4_csevl ( z, an21cs, nan21 );
    phi = - 0.625E+00 + r4_csevl ( z, aph1cs, naph1 );
  }
  else if ( x <= - 1.0E+00 )
  {
    z = ( 16.0E+00 / x / x / x + 9.0E+00 ) / 7.0E+00;
    ampl = 0.3125E+00 + r4_csevl ( z, an22cs, nan22 );
    phi = - 0.625E+00 + r4_csevl ( z, aph2cs, naph2 );
  }
  else
  {
    cerr << "\n";
    cerr << "R4_ADMP - Fatal error!\n";
    cerr << "  - 1.0 < x.\n";
    exit ( 1 );
  }

  sqrtx = sqrt ( - x );
  ampl = sqrt ( ampl * sqrtx );
  phi = pi34 - x * sqrtx * phi;

  return;
}
//****************************************************************************80

float r4_ai ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_AI evaluates the Airy function Ai of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_AI, the Airy function Ai of X.
//
{
  static float aifcs[9] = {
   -0.03797135849666999750E+00,
    0.05919188853726363857E+00,
    0.00098629280577279975E+00,
    0.00000684884381907656E+00,
    0.00000002594202596219E+00,
    0.00000000006176612774E+00,
    0.00000000000010092454E+00,
    0.00000000000000012014E+00,
    0.00000000000000000010E+00 };
  static float aigcs[8] = {
    0.01815236558116127E+00,
    0.02157256316601076E+00,
    0.00025678356987483E+00,
    0.00000142652141197E+00,
    0.00000000457211492E+00,
    0.00000000000952517E+00,
    0.00000000000001392E+00,
    0.00000000000000001E+00 };
  static int naif = 0;
  static int naig = 0;
  float theta;
  float value;
  static float x3sml = 0.0;
  float xm;
  static float xmax = 0.0;
  float z;

  if ( naif == 0 )
  {
    naif = r4_inits ( aifcs, 9, 0.1E+00 * r4_mach ( 3 ) );
    naig = r4_inits ( aigcs, 8, 0.1E+00 * r4_mach ( 3 ) );
    x3sml = r4_power ( r4_mach ( 3 ), 0.3334E+00 );
    xmax = r4_power ( - 1.5E+00 * log ( r4_mach ( 1 ) ), 0.6667E+00 );
    xmax = xmax - xmax * log ( xmax ) 
      / ( 4.0E+00 * xmax * sqrt ( xmax ) + 1.0E+00 ) - 0.01E+00;
  }

  if ( x < - 1.0E+00 )
  {
    r4_aimp ( x, xm, theta );
    value = xm * cos ( theta );
  }
  else if ( r4_abs ( x ) <= x3sml )
  {
    z = 0.0E+00;
    value = 0.375E+00 + ( r4_csevl ( z, aifcs, naif ) 
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) );
  }
  else if ( x <= 1.0E+00 )
  {
    z = x * x * x;
    value = 0.375E+00 + ( r4_csevl ( z, aifcs, naif ) 
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) );
  }
  else if ( x <= xmax )
  {
    value = r4_aie ( x )
      * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  else
  {
    value = 0.0E+00;
  }
  return value;
}
//****************************************************************************80

float r4_aid ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_AID evaluates the derivative of the Airy function Ai of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_AID, the derivative of the Airy function Ai of X.
//
{
  static float aifcs[8] = {
     0.10527461226531408809E+00,
     0.01183613628152997844E+00,
     0.00012328104173225664E+00,
     0.00000062261225638140E+00,
     0.00000000185298887844E+00,
     0.00000000000363328873E+00,
     0.00000000000000504622E+00,
     0.00000000000000000522E+00 };
  static float aigcs[9] = {
     0.021233878150918666852E+00,
     0.086315930335214406752E+00,
     0.001797594720383231358E+00,
     0.000014265499875550693E+00,
     0.000000059437995283683E+00,
     0.000000000152403366479E+00,
     0.000000000000264587660E+00,
     0.000000000000000331562E+00,
     0.000000000000000000314E+00 };
  static int naif = 0;
  static int naig = 0;
  float phi;
  float value;
  float x2;
  static float x2sml = 0.0;
  float x3;
  static float x3sml = 0.0;
  float xn;

  if ( naif == 0 )
  {
    naif = r4_inits ( aifcs, 8, 0.1E+00 * r4_mach ( 3 ) );
    naig = r4_inits ( aigcs, 9, 0.1E+00 * r4_mach ( 3 ) );
    x3sml = r4_power ( r4_mach ( 3 ), 0.3334E+00 );
    x2sml = sqrt ( r4_mach ( 3 ) );
  }

  if ( x < - 1.0E+00 )
  {
    r4_admp ( x, xn, phi );
    value = xn * cos ( phi );
  }
  else if ( r4_abs ( x ) <= x2sml )
  {
    x2 = 0.0E+00;
    x3 = 0.0E+00;
    value = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) -
      r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00;
  }
  else if ( r4_abs ( x ) <= x3sml )
  {
    x2 = x * x;
    x3 = 0.0E+00;
    value = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) -
      r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00;
  }
  else if ( x <= 1.0E+00 )
  {
    x2 = x * x;
    x3 = x * x * x;
    value = ( x2 * ( 0.125E+00 + r4_csevl ( x3, aifcs, naif ) ) -
      r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00;
  }
  else
  {
    value = r4_aide ( x ) 
      * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  return value;
}
//****************************************************************************80

float r4_aide ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_AIDE: exponentially scaled derivative, Airy function Ai of an R4 argument.
//
//  Discussion:
//
//    if X <= 0,
//      R4_AIDE ( X ) = R4_AID ( X )
//    else
//      R4_AIDE ( X ) = R4_AID ( X ) * exp ( 2/3 * X^(3/2) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_AIDE, the exponentially scaled derivative of 
//    the Airy function Ai of X.
//
{
  static float aifcs[8] = {
     0.10527461226531408809E+00,
     0.01183613628152997844E+00,
     0.00012328104173225664E+00,
     0.00000062261225638140E+00,
     0.00000000185298887844E+00,
     0.00000000000363328873E+00,
     0.00000000000000504622E+00,
     0.00000000000000000522E+00 };
  static float aigcs[9] = {
     0.021233878150918666852E+00,
     0.086315930335214406752E+00,
     0.001797594720383231358E+00,
     0.000014265499875550693E+00,
     0.000000059437995283683E+00,
     0.000000000152403366479E+00,
     0.000000000000264587660E+00,
     0.000000000000000331562E+00,
     0.000000000000000000314E+00 };
  static float aip1cs[25] = {
     0.0358865097808301538E+00,
     0.0114668575627764899E+00,
    -0.0007592073583861400E+00,
     0.0000869517610893841E+00,
    -0.0000128237294298592E+00,
     0.0000022062695681038E+00,
    -0.0000004222295185921E+00,
     0.0000000874686415726E+00,
    -0.0000000192773588418E+00,
     0.0000000044668460054E+00,
    -0.0000000010790108052E+00,
     0.0000000002700029447E+00,
    -0.0000000000696480108E+00,
     0.0000000000184489907E+00,
    -0.0000000000050027817E+00,
     0.0000000000013852243E+00,
    -0.0000000000003908218E+00,
     0.0000000000001121536E+00,
    -0.0000000000000326862E+00,
     0.0000000000000096619E+00,
    -0.0000000000000028935E+00,
     0.0000000000000008770E+00,
    -0.0000000000000002688E+00,
     0.0000000000000000832E+00,
    -0.0000000000000000260E+00 };
  static float aip2cs[15] = {
     0.0065457691989713757E+00,
     0.0023833724120774592E+00,
    -0.0000430700770220586E+00,
     0.0000015629125858629E+00,
    -0.0000000815417186163E+00,
     0.0000000054103738057E+00,
    -0.0000000004284130883E+00,
     0.0000000000389497963E+00,
    -0.0000000000039623161E+00,
     0.0000000000004428184E+00,
    -0.0000000000000536297E+00,
     0.0000000000000069650E+00,
    -0.0000000000000009620E+00,
     0.0000000000000001403E+00,
    -0.0000000000000000215E+00 };
  float eta;
  static int naif = 0;
  static int naig = 0;
  static int naip1 = 0;
  static int naip2 = 0;
  float phi;
  float sqrtx;
  float value;
  float x2;
  static float x2sml = 0.0;
  float x3;
  static float x32sml = 0.0;
  static float x3sml = 0.0;
  float xn;
  static float xbig = 0.0;
  float z;

  if ( naif == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    naif = r4_inits ( aifcs, 8, eta);
    naig = r4_inits ( aigcs, 9, eta);
    naip1 = r4_inits ( aip1cs, 25, eta );
    naip2 = r4_inits ( aip2cs, 15, eta );
    x2sml = sqrt ( eta );
    x3sml = r4_power ( eta, 0.3333E+00 );
    x32sml = 1.3104E+00 * x3sml * x3sml;
    xbig = r4_power ( r4_mach ( 2 ), 0.6666E+00 );
  }

  if ( x < - 1.0E+00 )
  {
    r4_admp ( x, xn, phi );
    value = xn * cos ( phi );
  }
  else if ( r4_abs ( x ) <= x2sml )
  {
    x2 = 0.0E+00;
    x3 = 0.0E+00;
    value = ( x2 * ( 0.125E+00 
      + r4_csevl ( x3, aifcs, naif ) ) 
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00;
  }
  else if ( r4_abs ( x ) <= x3sml )
  {
    x2 = x * x;
    x3 = 0.0E+00;
    value = ( x2 * ( 0.125E+00 
      + r4_csevl ( x3, aifcs, naif ) ) 
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00;
  }
  else if ( r4_abs ( x ) <= x32sml )
  {
    x2 = x * x;
    x3 = x * x * x;
    value = ( x2 * ( 0.125E+00 
      + r4_csevl ( x3, aifcs, naif ) ) 
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00;
  }
  else if ( x <= 1.0E+00 )
  {
    x2 = x * x;
    x3 = x * x * x;
    value = ( x2 * ( 0.125E+00 
      + r4_csevl ( x3, aifcs, naif ) ) 
      - r4_csevl ( x3, aigcs, naig ) ) - 0.25E+00;
    value = value * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  else if ( x <= 4.0E+00 )
  {
    sqrtx = sqrt ( x );
    z = ( 16.0E+00 / ( x * sqrtx ) - 9.0E+00 ) / 7.0E+00;
    value = ( - 0.28125E+00 
      - r4_csevl ( z, aip1cs, naip1 ) ) * sqrt ( sqrtx );
  }
  else if ( x < xbig )
  {
    sqrtx = sqrt ( x );
    z = 16.0E+00 / ( x * sqrtx ) - 1.0E+00;
    value = ( - 0.28125E+00 
      - r4_csevl ( z, aip2cs, naip2 ) ) * sqrt ( sqrtx );
  }
  else
  {
    sqrtx = sqrt ( x );
    z = - 1.0E+00;
    value = ( - 0.28125E+00 
      - r4_csevl ( z, aip2cs, naip2 ) ) * sqrt ( sqrtx );
  }

  return value;
}
//****************************************************************************80

float r4_aie ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_AIE evaluates the exponential scaled Airy function Ai of an R4 argument.
//
//  Discussion:
//
//    If X <= 0
//      R4_AIE ( X ) = R4_AI ( X )
//    else
//      R4_AIE ( X ) = R4_AI ( X ) * exp ( 2/3 X^(3/2) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_AIE, the Airy function Ai of X.
//
{
  static float aifcs[9] = {
   -0.03797135849666999750E+00,
    0.05919188853726363857E+00,
    0.00098629280577279975E+00,
    0.00000684884381907656E+00,
    0.00000002594202596219E+00,
    0.00000000006176612774E+00,
    0.00000000000010092454E+00,
    0.00000000000000012014E+00,
    0.00000000000000000010E+00 };
  static float aigcs[8] = {
    0.01815236558116127E+00,
    0.02157256316601076E+00,
    0.00025678356987483E+00,
    0.00000142652141197E+00,
    0.00000000457211492E+00,
    0.00000000000952517E+00,
    0.00000000000001392E+00,
    0.00000000000000001E+00 };
  static float aipcs[34] = {
   -0.0187519297793868E+00,
   -0.0091443848250055E+00,
    0.0009010457337825E+00,
   -0.0001394184127221E+00,
    0.0000273815815785E+00,
   -0.0000062750421119E+00,
    0.0000016064844184E+00,
   -0.0000004476392158E+00,
    0.0000001334635874E+00,
   -0.0000000420735334E+00,
    0.0000000139021990E+00,
   -0.0000000047831848E+00,
    0.0000000017047897E+00,
   -0.0000000006268389E+00,
    0.0000000002369824E+00,
   -0.0000000000918641E+00,
    0.0000000000364278E+00,
   -0.0000000000147475E+00,
    0.0000000000060851E+00,
   -0.0000000000025552E+00,
    0.0000000000010906E+00,
   -0.0000000000004725E+00,
    0.0000000000002076E+00,
   -0.0000000000000924E+00,
    0.0000000000000417E+00,
   -0.0000000000000190E+00,
    0.0000000000000087E+00,
   -0.0000000000000040E+00,
    0.0000000000000019E+00,
   -0.0000000000000009E+00,
    0.0000000000000004E+00,
   -0.0000000000000002E+00,
    0.0000000000000001E+00,
   -0.0000000000000000E+00 };
  float eta;
  static int naif = 0;
  static int naig = 0;
  static int naip = 0;
  float sqrtx;
  float theta;
  float value;
  static float x32sml = 0.0;
  static float x3sml = 0.0;
  static float xbig = 0.0;
  float xm;
  float z;

  if ( naif == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    naif = r4_inits ( aifcs, 9, eta );
    naig = r4_inits ( aigcs, 8, eta );
    naip = r4_inits ( aipcs, 34, eta );
    x3sml = r4_power ( eta, 0.3333E+00 );
    x32sml = 1.3104E+00 * x3sml * x3sml;
    xbig = r4_power ( r4_mach ( 2 ), 0.6666E+00 );
  }

  if ( x < - 1.0E+00 )
  {
    r4_aimp ( x, xm, theta );
    value = xm * cos ( theta );
  }
  else if ( r4_abs ( x ) <= x32sml )
  {
    z = 0.0E+00;
    value = 0.375E+00 + ( r4_csevl ( z, aifcs, naif ) 
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) );
  }
  else if ( r4_abs ( x ) <= x3sml )
  {
    z = 0.0E+00;
    value = 0.375E+00 + ( r4_csevl ( z, aifcs, naif ) 
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) );
    value = value * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  else if ( x <= 1.0E+00 )
  {
    z = x * x * x;
    value = 0.375E+00 + ( r4_csevl (z, aifcs, naif ) 
      - x * ( 0.25E+00 + r4_csevl ( z, aigcs, naig ) ) );
    value = value * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  else if ( x < xbig )
  {
    sqrtx = sqrt ( x );
    z = 2.0E+00 / ( x * sqrtx ) - 1.0E+00;
    value = ( 0.28125E+00 
      + r4_csevl ( z, aipcs, naip ) ) / sqrt ( sqrtx );
  }
  else
  {
    sqrtx = sqrt ( x );
    z = - 1.0E+00;
    value = ( 0.28125E+00 
      + r4_csevl ( z, aipcs, naip ) ) / sqrt ( sqrtx );
  }

  return value;
}
//****************************************************************************80

void r4_aimp ( float x, float &ampl, float &theta )

//****************************************************************************80
//
//  Purpose:
//
//    R4_AIMP evaluates the modulus and phase of the Airy function.
//
//  Description:
//
//    This function must only be called when X <= -1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float &AMPL, &PHI, the modulus and phase of the 
//    Airy function.
//
{
  static float am21cs[40] = {
   +0.0065809191761485E+00,
   +0.0023675984685722E+00,
   +0.0001324741670371E+00,
   +0.0000157600904043E+00,
   +0.0000027529702663E+00,
   +0.0000006102679017E+00,
   +0.0000001595088468E+00,
   +0.0000000471033947E+00,
   +0.0000000152933871E+00,
   +0.0000000053590722E+00,
   +0.0000000020000910E+00,
   +0.0000000007872292E+00,
   +0.0000000003243103E+00,
   +0.0000000001390106E+00,
   +0.0000000000617011E+00,
   +0.0000000000282491E+00,
   +0.0000000000132979E+00,
   +0.0000000000064188E+00,
   +0.0000000000031697E+00,
   +0.0000000000015981E+00,
   +0.0000000000008213E+00,
   +0.0000000000004296E+00,
   +0.0000000000002284E+00,
   +0.0000000000001232E+00,
   +0.0000000000000675E+00,
   +0.0000000000000374E+00,
   +0.0000000000000210E+00,
   +0.0000000000000119E+00,
   +0.0000000000000068E+00,
   +0.0000000000000039E+00,
   +0.0000000000000023E+00,
   +0.0000000000000013E+00,
   +0.0000000000000008E+00,
   +0.0000000000000005E+00,
   +0.0000000000000003E+00,
   +0.0000000000000001E+00,
   +0.0000000000000001E+00,
   +0.0000000000000000E+00,
   +0.0000000000000000E+00,
   +0.0000000000000000E+00 };
  static float am22cs[33] = {
   -0.01562844480625341E+00,
   +0.00778336445239681E+00,
   +0.00086705777047718E+00,
   +0.00015696627315611E+00,
   +0.00003563962571432E+00,
   +0.00000924598335425E+00,
   +0.00000262110161850E+00,
   +0.00000079188221651E+00,
   +0.00000025104152792E+00,
   +0.00000008265223206E+00,
   +0.00000002805711662E+00,
   +0.00000000976821090E+00,
   +0.00000000347407923E+00,
   +0.00000000125828132E+00,
   +0.00000000046298826E+00,
   +0.00000000017272825E+00,
   +0.00000000006523192E+00,
   +0.00000000002490471E+00,
   +0.00000000000960156E+00,
   +0.00000000000373448E+00,
   +0.00000000000146417E+00,
   +0.00000000000057826E+00,
   +0.00000000000022991E+00,
   +0.00000000000009197E+00,
   +0.00000000000003700E+00,
   +0.00000000000001496E+00,
   +0.00000000000000608E+00,
   +0.00000000000000248E+00,
   +0.00000000000000101E+00,
   +0.00000000000000041E+00,
   +0.00000000000000017E+00,
   +0.00000000000000007E+00,
   +0.00000000000000002E+00 };
  static float ath1cs[36] = {
   -0.07125837815669365E+00,
   -0.00590471979831451E+00,
   -0.00012114544069499E+00,
   -0.00000988608542270E+00,
   -0.00000138084097352E+00,
   -0.00000026142640172E+00,
   -0.00000006050432589E+00,
   -0.00000001618436223E+00,
   -0.00000000483464911E+00,
   -0.00000000157655272E+00,
   -0.00000000055231518E+00,
   -0.00000000020545441E+00,
   -0.00000000008043412E+00,
   -0.00000000003291252E+00,
   -0.00000000001399875E+00,
   -0.00000000000616151E+00,
   -0.00000000000279614E+00,
   -0.00000000000130428E+00,
   -0.00000000000062373E+00,
   -0.00000000000030512E+00,
   -0.00000000000015239E+00,
   -0.00000000000007758E+00,
   -0.00000000000004020E+00,
   -0.00000000000002117E+00,
   -0.00000000000001132E+00,
   -0.00000000000000614E+00,
   -0.00000000000000337E+00,
   -0.00000000000000188E+00,
   -0.00000000000000105E+00,
   -0.00000000000000060E+00,
   -0.00000000000000034E+00,
   -0.00000000000000020E+00,
   -0.00000000000000011E+00,
   -0.00000000000000007E+00,
   -0.00000000000000004E+00,
   -0.00000000000000002E+00 };
  static float ath2cs[32] = {
   +0.00440527345871877E+00,
   -0.03042919452318455E+00,
   -0.00138565328377179E+00,
   -0.00018044439089549E+00,
   -0.00003380847108327E+00,
   -0.00000767818353522E+00,
   -0.00000196783944371E+00,
   -0.00000054837271158E+00,
   -0.00000016254615505E+00,
   -0.00000005053049981E+00,
   -0.00000001631580701E+00,
   -0.00000000543420411E+00,
   -0.00000000185739855E+00,
   -0.00000000064895120E+00,
   -0.00000000023105948E+00,
   -0.00000000008363282E+00,
   -0.00000000003071196E+00,
   -0.00000000001142367E+00,
   -0.00000000000429811E+00,
   -0.00000000000163389E+00,
   -0.00000000000062693E+00,
   -0.00000000000024260E+00,
   -0.00000000000009461E+00,
   -0.00000000000003716E+00,
   -0.00000000000001469E+00,
   -0.00000000000000584E+00,
   -0.00000000000000233E+00,
   -0.00000000000000093E+00,
   -0.00000000000000037E+00,
   -0.00000000000000015E+00,
   -0.00000000000000006E+00,
   -0.00000000000000002E+00 };
  float eta;
  static int nam21 = 0;
  static int nam22 = 0;
  static int nath1 = 0;
  static int nath2 = 0;
  static float pi4 = 0.78539816339744831E+00;
  float sqrtx;
  static float xsml = 0.0;
  float z;

  if ( nam21 == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    nam21 = r4_inits ( am21cs, 40, eta );
    nath1 = r4_inits ( ath1cs, 36, eta );
    nam22 = r4_inits ( am22cs, 33, eta );
    nath2 = r4_inits ( ath2cs, 32, eta );
    xsml = - r4_power ( 16.0E+00 / r4_mach ( 3 ), 0.3333E+00 );
  }

  if ( x <= xsml )
  {
    z = 1.0E+00;
    ampl = 0.3125E+00 + r4_csevl ( z, am21cs, nam21 );
    theta = - 0.625E+00 + r4_csevl ( z, ath1cs, nath1 );
  }
  else if ( x < - 2.0E+00 )
  {
    z = 16.0 / x / x / x + 1.0;
    ampl = 0.3125E+00 + r4_csevl ( z, am21cs, nam21 );
    theta = - 0.625E+00 + r4_csevl ( z, ath1cs, nath1 );
  }
  else if ( x <= - 1.0E+00 )
  {
    z = ( 16.0E+00 / x / x / x + 9.0E+00 ) / 7.0E+00;
    ampl = 0.3125E+00 + r4_csevl ( z, am22cs, nam22 );
    theta = - 0.625E+00 + r4_csevl ( z, ath2cs, nath2 );
  }
  else
  {
    cerr << "\n";
    cerr << "R4_AIMP - Fatal error!\n";
    cerr << "  - 1.0 < X.\n";
    exit ( 1 );
  }

  sqrtx = sqrt ( - x );
  ampl = sqrt ( ampl / sqrtx );
  theta = pi4 - x * sqrtx * theta;

  return;
}
//****************************************************************************80

float r4_aint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_AINT truncates an R4 argument to an integer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    1 September 2011
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_AINT, the truncated version of X.
//
{
  float value;

  if ( x < 0.0E+00 )
  {
    value = - ( float ) ( ( int ) ( r4_abs ( x ) ) );
  }
  else
  {
    value =   ( float ) ( ( int ) ( r4_abs ( x ) ) );
  }

  return value;
}
//****************************************************************************80

float r4_asin ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ASIN evaluates the arc-sine of an R4 argument.
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
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_ASIN, the arc-sine of X.
//
{
  static float asincs[20] = {
    0.10246391753227159E+00,
    0.054946487221245833E+00,
    0.004080630392544969E+00,
    0.000407890068546044E+00,
    0.000046985367432203E+00,
    0.000005880975813970E+00,
    0.000000777323124627E+00,
    0.000000106774233400E+00,
    0.000000015092399536E+00,
    0.000000002180972408E+00,
    0.000000000320759842E+00,
    0.000000000047855369E+00,
    0.000000000007225128E+00,
    0.000000000001101833E+00,
    0.000000000000169476E+00,
    0.000000000000026261E+00,
    0.000000000000004095E+00,
    0.000000000000000642E+00,
    0.000000000000000101E+00,
    0.000000000000000016E+00 };
  static int nterms = 0;
  static float pi2 = 1.57079632679489661923E+00;
  static float sqeps = 0.0;
  float value;
  float y;
  float z;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( asincs, 20, 0.1 * r4_mach ( 3 ) );
    sqeps = sqrt ( 6.0 * r4_mach ( 3 ) );
  }

  y = r4_abs ( x );

  if ( x < - 1.0 - sqeps )
  {
    cerr << "\n";
    cerr << "R4_ASIN - Fatal error!\n";
    cerr << "  X < - 1.0\n";
    exit ( 1 );
  }
  else if ( x < - 1.0 )
  {
    value = - pi2;
  }
  else if ( x < 1.0 )
  {
    z = 0.0;
    if ( sqeps < y )
    {
      z = y * y;
    }

    if ( z <= 0.5 )
    {
      value = x * ( 1.0 + r4_csevl ( 4.0 * z - 1.0, asincs, nterms ) );
    }
    else
    {
      value = pi2 - sqrt ( 1.0 - z ) * ( 1.0 + 
        r4_csevl ( 3.0 - 4.0 * z, asincs, nterms ) );
    }

    if ( x < 0.0 )
    {
      value = - r4_abs ( value );
    }
    else if ( 0.0 < x )
    {
      value = + r4_abs ( value );
    }
  }
  else if ( x < 1.0 + sqeps )
  {
    value = pi2;
  }
  else
  {
    cerr << "\n";
    cerr << "R4_ASIN - Fatal error!\n";
    cerr << "  1.0 < X\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

float r4_asinh ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ASINH evaluates the arc-sine of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_ASINH, the arc-hyperbolic sine of X.
//
{
  static float aln2 = 0.69314718055994530942E+00;
  static float asnhcs[20] = {
   -0.12820039911738186E+00,
   -0.058811761189951768E+00,
    0.004727465432212481E+00,
   -0.000493836316265361E+00,
    0.000058506207058557E+00,
   -0.000007466998328931E+00,
    0.000001001169358355E+00,
   -0.000000139035438587E+00,
    0.000000019823169483E+00,
   -0.000000002884746841E+00,
    0.000000000426729654E+00,
   -0.000000000063976084E+00,
    0.000000000009699168E+00,
   -0.000000000001484427E+00,
    0.000000000000229037E+00,
   -0.000000000000035588E+00,
    0.000000000000005563E+00,
   -0.000000000000000874E+00,
    0.000000000000000138E+00,
   -0.000000000000000021E+00 };
  static int nterms = 0;
  static float sqeps = 0.0;
  float value;
  static float xmax = 0.0;
  float y;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( asnhcs, 20, 0.1E+00 * r4_mach ( 3 ) );
    sqeps = sqrt ( r4_mach ( 3 ) );
    xmax = 1.0E+00 / sqeps;
  }

  y = r4_abs ( x );

  if ( y <= 1.0E+00 )
  {
    value = x;
    if ( sqeps < y )
    {
      value = x * ( 1.0E+00 
        + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, asnhcs, nterms ) );
    }
  }
  else
  {
    if ( y < xmax )
    {
      value = log ( y + sqrt ( y * y + 1.0E+00 ) );
    }
    else
    {
      value = aln2 + log ( y );
    }

    if ( x < 0.0E+00 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

float r4_atan ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ATAN evaluates the arc-tangent of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_ATAN, the arc-tangent of X.
//
{
  static float atancs[9] = {
     0.48690110349241406E+00,
    -0.006510831636717464E+00,
     0.000038345828265245E+00,
    -0.000000268722128762E+00,
     0.000000002050093098E+00,
    -0.000000000016450717E+00,
     0.000000000000136509E+00,
    -0.000000000000001160E+00,
     0.000000000000000010E+00 };
  static float conpi8[4] = {
    0.375E+00,
    0.75E+00,
    1.125E+00,
    1.5E+00 };
  int n;
  static int nterms = 0;
  float pi8[4] = {
    0.176990816987241548E-01,
    0.353981633974483096E-01,
    0.530972450961724644E-01,
    0.0707963267948966192E+00 };
  static float sqeps = 0.0;
  float t;
  static float tanp8[3] = {
    0.414213562373095048E+00,
    1.0E+00,
    2.41421356237309504E+00 };
  float value;
  static float xbig = 0.0;
  static float xbnd1 = +0.198912367379658006E+00;
  static float xbnd2 = +0.668178637919298919E+00;
  static float xbnd3 = +1.49660576266548901E+00;
  static float xbnd4 = +5.02733949212584810E+00;
  float y;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( atancs, 9, 0.1E+00 * r4_mach ( 3 ) );
    sqeps = sqrt ( 6.0E+00 * r4_mach ( 3 ) );
    xbig = 1.0E+00 / r4_mach ( 3 );
  }

  y = r4_abs ( x );

  if ( y <= xbnd1 )
  {
    value = x;
    if ( sqeps < y )
    {
      value = x * ( 0.75E+00 + r4_csevl ( 
        50.0E+00 * y * y - 1.0E+00, atancs, nterms ) );
    }
  }
  else if ( y <= xbnd4 )
  {
    if ( xbnd3 < y )
    {
      n = 3;
    }
    else if ( xbnd2 < y )
    {
      n = 2;
    }
    else
    {
      n = 1;
    }

    t = ( y - tanp8[n-1] ) / ( 1.0E+00 + y * tanp8[n-1] );

    value = conpi8[n-1] + ( pi8[n-1] + t * ( 0.75E+00 +
      r4_csevl ( 50.0E+00 * t * t - 1.0E+00, atancs, nterms ) ) );
  }
  else
  {
    value = conpi8[3] + pi8[3];

    if ( y < xbig )
    {
      value = conpi8[3] + ( pi8[3] - ( 0.75E+00 +
        r4_csevl ( 50.0E+00 / y / y - 1.0E+00, atancs, 
        nterms ) ) / y );
    }
  }

  if ( x < 0.0E+00 )
  {
    value = - r4_abs ( value );
  }
  else
  {
    value = + r4_abs ( value );
  }
  return value;
}
//****************************************************************************80

float r4_atan2 ( float sn, float cs )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ATAN2 evaluates the arc-tangent of two R4 arguments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double SN, CS, the Y and X coordinates of a point on the angle.
//
//    Output, double R4_ATAN2, the arc-tangent of the angle.
//
{
  float abscs;
  float abssn;
  static float big = 0.0;
  static float pi = 3.14159265358979323846E+00;
  static float sml = 0.0;
  float value;

  if ( sml == 0.0E+00 )
  {
    sml = r4_mach ( 1 );
    big = r4_mach ( 2 );
  }
//
//  We now make sure SN can be divided by CS.  It is painful.
//
  abssn = r4_abs ( sn );
  abscs = r4_abs ( cs );

  if ( abscs <= abssn )
  {
    if ( abscs < 1.0E+00 && abscs * big <= abssn )
    {
      if ( sn < 0.0E+00 )
      {
        value = - 0.5E+00 * pi;
      }
      else if ( sn == 0.0E+00 )
      {
        cerr << "\n";
        cerr << "R4_ATAN2 - Fatal error!\n";
        cerr << "  Both arguments are 0.\n";
        exit ( 1 );
      }
      else
      {
        value = 0.5E+00 * pi;
      }
      return value;
    }
  }
  else
  {
    if ( 1.0E+00 < abscs && abssn <= abscs * sml )
    {
      if ( 0.0E+00 <= cs )
      {
        value = 0.0E+00;
      }
      else
      {
        value = pi;
      }
      return value;
    }
  }

  value = atan ( sn / cs );

  if ( cs < 0.0E+00 )
  {
    value = value + pi;
  }

  if ( pi < value )
  {
    value = value - 2.0E+00 * pi;
  }
  return value;
}
//****************************************************************************80

float r4_atanh ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ATANH evaluates the arc-hyperbolic tangent of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_ATANH, the arc-hyperbolic tangent of X.
//
{
  static float atnhcs[15] = {
    0.094395102393195492E+00,
    0.049198437055786159E+00,
    0.002102593522455432E+00,
    0.000107355444977611E+00,
    0.000005978267249293E+00,
    0.000000350506203088E+00,
    0.000000021263743437E+00,
    0.000000001321694535E+00,
    0.000000000083658755E+00,
    0.000000000005370503E+00,
    0.000000000000348665E+00,
    0.000000000000022845E+00,
    0.000000000000001508E+00,
    0.000000000000000100E+00,
    0.000000000000000006E+00 };
  static float dxrel = 0.0;
  static int nterms = 0;
  static float sqeps = 0.0;
  float value;
  float y;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( atnhcs, 15, 0.1E+00 * r4_mach ( 3 ) );
    dxrel = sqrt ( r4_mach ( 4 ) );
    sqeps = sqrt ( 3.0E+00 * r4_mach ( 3 ) );
  }

  y = r4_abs ( x );

  if ( y <= sqeps )
  {
    value = x;
  }
  else if ( y <= 0.5E+00 )
  {
    value = x * ( 1.0E+00 
      + r4_csevl ( 8.0E+00 * x * x - 1.0E+00, atnhcs, nterms ) );
  }
  else if ( y < 1.0E+00 )
  {
    value = 0.5E+00 * log ( ( 1.0E+00 + x ) / ( 1.0E+00 - x ) );
  }
  else
  {
    cerr << "\n";
    cerr << "R4_ATANH - Fatal error!\n";
    cerr << "  1 <= |X|.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

float r4_besi0 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESI0 evaluates the Bessel function I of order 0 of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESI0, the Bessel function I of order 0 of X.
//
{
  static float bi0cs[12] = {
   -0.07660547252839144951E+00,
    1.927337953993808270E+00,
    0.2282644586920301339E+00,
    0.01304891466707290428E+00,
    0.00043442709008164874E+00,
    0.00000942265768600193E+00,
    0.00000014340062895106E+00,
    0.00000000161384906966E+00,
    0.00000000001396650044E+00,
    0.00000000000009579451E+00,
    0.00000000000000053339E+00,
    0.00000000000000000245E+00 };
  static int nti0 = 0;
  float value;
  static float xmax = 0.0;
  static float xsml = 0.0;
  float y;

  if ( nti0 == 0 )
  {
    nti0 = r4_inits ( bi0cs, 12, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) );
    xmax = log ( r4_mach ( 2 ) );
  }

  y = r4_abs ( x );

  if ( y <= xsml )
  {
    value = 1.0E+00;
  }
  else if ( y <= 3.0E+00 )
  {
    value = 2.75E+00 + r4_csevl ( y * y / 4.5E+00 - 1.0E+00, 
      bi0cs, nti0 );
  }
  else if ( y <= xmax )
  {
    value = exp ( y ) * r4_besi0e ( x );
  }
  else
  {
    cerr << "\n";
    cerr << "R4_BESI0 - Fatal error!\n";
    cerr << "  Result overflows.\n";
    cerr << "  |X| too large.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

float r4_besi0e ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESI0E evaluates the exponentially scaled Bessel function I0(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESI0E, the exponentially scaled Bessel 
//    function I0(X).
//
{
  static float ai02cs[22] = {
    0.05449041101410882E+00,
    0.00336911647825569E+00,
    0.00006889758346918E+00,
    0.00000289137052082E+00,
    0.00000020489185893E+00,
    0.00000002266668991E+00,
    0.00000000339623203E+00,
    0.00000000049406022E+00,
    0.00000000001188914E+00,
   -0.00000000003149915E+00,
   -0.00000000001321580E+00,
   -0.00000000000179419E+00,
    0.00000000000071801E+00,
    0.00000000000038529E+00,
    0.00000000000001539E+00,
   -0.00000000000004151E+00,
   -0.00000000000000954E+00,
    0.00000000000000382E+00,
    0.00000000000000176E+00,
   -0.00000000000000034E+00,
   -0.00000000000000027E+00,
    0.00000000000000003E+00 };
  static float ai0cs[21] = {
    0.07575994494023796E+00,
    0.00759138081082334E+00,
    0.00041531313389237E+00,
    0.00001070076463439E+00,
   -0.00000790117997921E+00,
   -0.00000078261435014E+00,
    0.00000027838499429E+00,
    0.00000000825247260E+00,
   -0.00000001204463945E+00,
    0.00000000155964859E+00,
    0.00000000022925563E+00,
   -0.00000000011916228E+00,
    0.00000000001757854E+00,
    0.00000000000112822E+00,
   -0.00000000000114684E+00,
    0.00000000000027155E+00,
   -0.00000000000002415E+00,
   -0.00000000000000608E+00,
    0.00000000000000314E+00,
   -0.00000000000000071E+00,
    0.00000000000000007E+00 };
  static float bi0cs[12] = {
   -0.07660547252839144951E+00,
    1.927337953993808270E+00,
    0.2282644586920301339E+00,
    0.01304891466707290428E+00,
    0.00043442709008164874E+00,
    0.00000942265768600193E+00,
    0.00000014340062895106E+00,
    0.00000000161384906966E+00,
    0.00000000001396650044E+00,
    0.00000000000009579451E+00,
    0.00000000000000053339E+00,
    0.00000000000000000245E+00 };
  static int ntai0 = 0;
  static int ntai02 = 0;
  static int nti0 = 0;
  float value;
  static float xsml = 0.0;
  float y;

  if ( nti0 == 0 )
  {
    nti0 = r4_inits ( bi0cs, 12, 0.1E+00 * r4_mach ( 3 ) );
    ntai0 = r4_inits ( ai0cs, 21, 0.1E+00 * r4_mach ( 3 ) );
    ntai02 = r4_inits ( ai02cs, 22, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) );
  }

  y = r4_abs ( x );

  if ( y <= xsml )
  {
    value = 1.0E+00;
  }
  else if ( y <= 3.0E+00 )
  {
    value = exp ( - y ) * ( 2.75E+00 +
      r4_csevl ( y * y / 4.5E+00 - 1.0E+00, bi0cs, nti0 ) );
  }
  else if ( y <= 8.0E+00 )
  {
    value = ( 0.375E+00 + r4_csevl 
      ( ( 48.0E+00 / y - 11.0E+00 ) / 5.0E+00, ai0cs, ntai0 ) ) 
      / sqrt ( y );
  }
  else
  {
    value = ( 0.375E+00 + r4_csevl 
      ( 16.0E+00 / y - 1.0E+00, ai02cs, ntai02 ) ) / sqrt ( y );
  }
  return value;
}
//****************************************************************************80

float r4_besi1 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESI1 evaluates the Bessel function I of order 1 of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESI1, the Bessel function I of order 1 of X.
//
{
  static float bi1cs[11] = {
   -0.001971713261099859E+00,
    0.40734887667546481E+00,
    0.034838994299959456E+00,
    0.001545394556300123E+00,
    0.000041888521098377E+00,
    0.000000764902676483E+00,
    0.000000010042493924E+00,
    0.000000000099322077E+00,
    0.000000000000766380E+00,
    0.000000000000004741E+00,
    0.000000000000000024E+00 };
  static int nti1 = 0;
  float value;
  static float xmax = 0.0;
  static float xmin = 0.0;
  static float xsml = 0.0;
  float y;

  if ( nti1 == 0 )
  {
    nti1 = r4_inits ( bi1cs, 11, 0.1E+00 * r4_mach ( 3 ) );
    xmin = 2.0E+00 * r4_mach ( 1 );
    xsml = sqrt ( 8.0E+00 * r4_mach ( 3 ) );
    xmax = log ( r4_mach ( 2 ) );
  }

  y = r4_abs ( x );

  if ( y <= xmin )
  {
    value = 0.0E+00;
  }
  else if ( y <= xsml )
  {
    value = 0.5E+00 * x;
  }
  else if ( y <= 3.0E+00 )
  {
    value = x * ( 0.875E+00 + r4_csevl 
      ( y * y / 4.5E+00 - 1.0E+00, bi1cs, nti1 ) );
  }
  else if ( y <= xmax )
  {
    value = exp ( y ) * r4_besi1e ( x );
  }
  else
  {
    cerr << "\n";
    cerr << "R4_BESI1 - Fatal error!\n";
    cerr << "  Result overflows.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

float r4_besi1e ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESI1E: exponentially scaled Bessel function I of order 1 of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESI1E, the exponentially scaled Bessel function I 
//    of order 1 of X.
//
{
  static float ai1cs[21] = {
   -0.02846744181881479E+00,
   -0.01922953231443221E+00,
   -0.00061151858579437E+00,
   -0.00002069971253350E+00,
    0.00000858561914581E+00,
    0.00000104949824671E+00,
   -0.00000029183389184E+00,
   -0.00000001559378146E+00,
    0.00000001318012367E+00,
   -0.00000000144842341E+00,
   -0.00000000029085122E+00,
    0.00000000012663889E+00,
   -0.00000000001664947E+00,
   -0.00000000000166665E+00,
    0.00000000000124260E+00,
   -0.00000000000027315E+00,
    0.00000000000002023E+00,
    0.00000000000000730E+00,
   -0.00000000000000333E+00,
    0.00000000000000071E+00,
   -0.00000000000000006E+00 };
  static float ai12cs[22] = {
    0.02857623501828014E+00,
   -0.00976109749136147E+00,
   -0.00011058893876263E+00,
   -0.00000388256480887E+00,
   -0.00000025122362377E+00,
   -0.00000002631468847E+00,
   -0.00000000383538039E+00,
   -0.00000000055897433E+00,
   -0.00000000001897495E+00,
    0.00000000003252602E+00,
    0.00000000001412580E+00,
    0.00000000000203564E+00,
   -0.00000000000071985E+00,
   -0.00000000000040836E+00,
   -0.00000000000002101E+00,
    0.00000000000004273E+00,
    0.00000000000001041E+00,
   -0.00000000000000382E+00,
   -0.00000000000000186E+00,
    0.00000000000000033E+00,
    0.00000000000000028E+00,
   -0.00000000000000003E+00 };
  static float bi1cs[11] = {
   -0.001971713261099859E+00,
    0.40734887667546481E+00,
    0.034838994299959456E+00,
    0.001545394556300123E+00,
    0.000041888521098377E+00,
    0.000000764902676483E+00,
    0.000000010042493924E+00,
    0.000000000099322077E+00,
    0.000000000000766380E+00,
    0.000000000000004741E+00,
    0.000000000000000024E+00 };
  static int ntai1 = 0;
  static int ntai12 = 0;
  static int nti1 = 0;
  float value;
  static float xmin = 0.0;
  static float xsml = 0.0;
  float y;

  if ( nti1 == 0 )
  {
    nti1 = r4_inits ( bi1cs, 11, 0.1E+00 * r4_mach ( 3 ) );
    ntai1 = r4_inits ( ai1cs, 21, 0.1E+00 * r4_mach ( 3 ) );
    ntai12 = r4_inits ( ai12cs, 22, 0.1E+00 * r4_mach ( 3 ) );
    xmin = 2.0E+00 * r4_mach ( 1 );
    xsml = sqrt ( 8.0E+00 * r4_mach ( 3 ) );
  }

  y = r4_abs ( x );

  if ( x == 0.0E+00 )
  {
    value = 0.0E+00;
  }
  else if ( y <= xmin )
  {
    value = 0.0E+00;
  }
  else if ( y <= xsml )
  {
    value = 0.5E+00 * x;
    value = exp ( - y ) * value;
  }
  else if ( y <= 3.0E+00 )
  {
    value = x * ( 0.875E+00 
      + r4_csevl ( y * y / 4.5E+00 - 1.0E+00, bi1cs, nti1 ) );
    value = exp ( - y ) * value;
  }
  else if ( y <= 8.0E+00 )
  {
    value = ( 0.375E+00 + r4_csevl ( ( 48.0E+00 / y - 11.0E+00 ) / 5.0E+00, 
      ai1cs, ntai1) ) / sqrt ( y );
    if ( x < 0.0E+00 )
    {
      value = - value;
    }
  }
  else
  {
    value = ( 0.375E+00 
      + r4_csevl ( 16.0E+00 / y - 1.0E+00, ai12cs, ntai12 ) ) / sqrt ( y );
    if ( x < 0.0E+00 )
    {
      value = - value;
    }
  }

  return value;
}
//****************************************************************************80

float r4_besj0 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESJ0 evaluates the Bessel function J of order 0 of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESJ0, the Bessel function J of order 0 of X.
//
{
  float ampl;
  static float bj0cs[13] = {
    0.100254161968939137E+00,
   -0.665223007764405132E+00,
    0.248983703498281314E+00,
   -0.0332527231700357697E+00,
    0.0023114179304694015E+00,
   -0.0000991127741995080E+00,
    0.0000028916708643998E+00,
   -0.0000000612108586630E+00,
    0.0000000009838650793E+00,
   -0.0000000000124235515E+00,
    0.0000000000001265433E+00,
   -0.0000000000000010619E+00,
    0.0000000000000000074E+00 };
  static float bm0cs[21] = {
    0.09284961637381644E+00,
   -0.00142987707403484E+00,
    0.00002830579271257E+00,
   -0.00000143300611424E+00,
    0.00000012028628046E+00,
   -0.00000001397113013E+00,
    0.00000000204076188E+00,
   -0.00000000035399669E+00,
    0.00000000007024759E+00,
   -0.00000000001554107E+00,
    0.00000000000376226E+00,
   -0.00000000000098282E+00,
    0.00000000000027408E+00,
   -0.00000000000008091E+00,
    0.00000000000002511E+00,
   -0.00000000000000814E+00,
    0.00000000000000275E+00,
   -0.00000000000000096E+00,
    0.00000000000000034E+00,
   -0.00000000000000012E+00,
    0.00000000000000004E+00 };
  static float bth0cs[24] = {
   -0.24639163774300119E+00,
    0.001737098307508963E+00,
   -0.000062183633402968E+00,
    0.000004368050165742E+00,
   -0.000000456093019869E+00,
    0.000000062197400101E+00,
   -0.000000010300442889E+00,
    0.000000001979526776E+00,
   -0.000000000428198396E+00,
    0.000000000102035840E+00,
   -0.000000000026363898E+00,
    0.000000000007297935E+00,
   -0.000000000002144188E+00,
    0.000000000000663693E+00,
   -0.000000000000215126E+00,
    0.000000000000072659E+00,
   -0.000000000000025465E+00,
    0.000000000000009229E+00,
   -0.000000000000003448E+00,
    0.000000000000001325E+00,
   -0.000000000000000522E+00,
    0.000000000000000210E+00,
   -0.000000000000000087E+00,
    0.000000000000000036E+00 };
  static int ntj0 = 0;
  static int ntm0 = 0;
  static int ntth0 = 0;
  static float pi4 = 0.78539816339744831E+00;
  float theta;
  float value;
  static float xmax = 0.0;
  static float xsml = 0.0;
  float y;
  float z;

  if ( ntj0 == 0 )
  {
    ntj0 = r4_inits ( bj0cs, 13, 0.1E+00 * r4_mach ( 3 ) );
    ntm0 = r4_inits ( bm0cs, 21, 0.1E+00 * r4_mach ( 3 ) );
    ntth0 = r4_inits ( bth0cs, 24, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) );
    xmax = 1.0E+00 / r4_mach ( 4 );
  }

  y = r4_abs ( x );
 
  if ( y <= xsml )
  {
    value = 1.0E+00;
  }
  else if ( y <= 4.0E+00 )
  {
    value = r4_csevl ( 0.125E+00 * y * y - 1.0E+00, bj0cs, ntj0 );
  }
  else
  {
    z = 32.0E+00 / y / y - 1.0E+00;
    ampl = ( 0.75E+00 + r4_csevl ( z, bm0cs, ntm0 ) ) / sqrt ( y );
    theta = y - pi4 + r4_csevl ( z, bth0cs, ntth0 ) / y;
    value = ampl * cos ( theta );
  }

  return value;
}
//****************************************************************************80

float r4_besj1 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESJ1 evaluates the Bessel function J of order 1 of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESJ1, the Bessel function J of order 1 of X.
//
{
  float ampl;
  static float bj1cs[12] = {
   -0.11726141513332787E+00,
   -0.25361521830790640E+00,
   +0.050127080984469569E+00,
   -0.004631514809625081E+00,
   +0.000247996229415914E+00,
   -0.000008678948686278E+00,
   +0.000000214293917143E+00,
   -0.000000003936093079E+00,
   +0.000000000055911823E+00,
   -0.000000000000632761E+00,
   +0.000000000000005840E+00,
   -0.000000000000000044E+00 };
  static float bm1cs[21] = {
   +0.1047362510931285E+00,
   +0.00442443893702345E+00,
   -0.00005661639504035E+00,
   +0.00000231349417339E+00,
   -0.00000017377182007E+00,
   +0.00000001893209930E+00,
   -0.00000000265416023E+00,
   +0.00000000044740209E+00,
   -0.00000000008691795E+00,
   +0.00000000001891492E+00,
   -0.00000000000451884E+00,
   +0.00000000000116765E+00,
   -0.00000000000032265E+00,
   +0.00000000000009450E+00,
   -0.00000000000002913E+00,
   +0.00000000000000939E+00,
   -0.00000000000000315E+00,
   +0.00000000000000109E+00,
   -0.00000000000000039E+00,
   +0.00000000000000014E+00,
   -0.00000000000000005E+00 };
  static float bth1cs[24] = {
   +0.74060141026313850E+00,
   -0.004571755659637690E+00,
   +0.000119818510964326E+00,
   -0.000006964561891648E+00,
   +0.000000655495621447E+00,
   -0.000000084066228945E+00,
   +0.000000013376886564E+00,
   -0.000000002499565654E+00,
   +0.000000000529495100E+00,
   -0.000000000124135944E+00,
   +0.000000000031656485E+00,
   -0.000000000008668640E+00,
   +0.000000000002523758E+00,
   -0.000000000000775085E+00,
   +0.000000000000249527E+00,
   -0.000000000000083773E+00,
   +0.000000000000029205E+00,
   -0.000000000000010534E+00,
   +0.000000000000003919E+00,
   -0.000000000000001500E+00,
   +0.000000000000000589E+00,
   -0.000000000000000237E+00,
   +0.000000000000000097E+00,
   -0.000000000000000040E+00 };
  static int ntj1 = 0;
  static int ntm1 = 0;
  static int ntth1 = 0;
  static float pi4 = 0.78539816339744831E+00;
  float theta;
  float value;
  static float xmax = 0.0;
  static float xmin = 0.0;
  static float xsml = 0.0;
  float y;
  float z;

  if ( ntj1 == 0 )
  {
    ntj1 = r4_inits ( bj1cs, 12, 0.1E+00 * r4_mach ( 3 ) );
    ntm1 = r4_inits ( bm1cs, 21, 0.1E+00 * r4_mach ( 3 ) );
    ntth1 = r4_inits ( bth1cs, 24, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( 8.0E+00 * r4_mach ( 3 ) );
    xmin = 2.0E+00 * r4_mach ( 1 );
    xmax = 1.0E+00 / r4_mach ( 4 );
  }

  y = r4_abs ( x );

  if ( y <= xmin )
  {
    value = 0.0E+00;
  }
  else if ( y <= xsml )
  {
    value = 0.5E+00 * x;
  }
  else if ( y <= 4.0E+00 )
  {
    value = x * ( 0.25E+00 
      + r4_csevl ( 0.125E+00 * y * y - 1.0E+00, bj1cs, ntj1 ) );
  }
  else
  {
    z = 32.0E+00 / y / y - 1.0E+00;
    ampl = ( 0.75E+00 + r4_csevl ( z, bm1cs, ntm1 ) ) / sqrt ( y );
    theta = y - 3.0E+00 * pi4 + r4_csevl ( z, bth1cs, ntth1 ) / y;
    if ( x < 0.0E+00 )
    {
      value = - ampl * cos ( theta );
    }
    else
    {
      value = + ampl * cos ( theta );
    }
  }
  return value;
}
//****************************************************************************80

float r4_besk0 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESK0 evaluates the Bessel function K of order 0 of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESK0, the Bessel function K of order 0 of X.
//
{
  static float bk0cs[11] = {
   -0.03532739323390276872E+00,
    0.3442898999246284869E+00,
    0.03597993651536150163E+00,
    0.00126461541144692592E+00,
    0.00002286212103119451E+00,
    0.00000025347910790261E+00,
    0.00000000190451637722E+00,
    0.00000000001034969525E+00,
    0.00000000000004259816E+00,
    0.00000000000000013744E+00,
    0.00000000000000000035E+00 };
  static int ntk0 = 0;
  float value;
  static float xmax = 0.0;
  static float xsml = 0.0;
  float y;

  if ( ntk0 == 0 )
  {
    ntk0 = r4_inits ( bk0cs, 11, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) );
    xmax = - log ( r4_mach ( 1 ) );
    xmax = xmax - 0.5E+00 * xmax * log ( xmax ) / ( xmax + 0.5E+00 ) - 0.01E+00;
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_BESK0 = Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0E+00;
    value = - log ( 0.5E+00 * x ) * r4_besi0 ( x ) 
      - 0.25E+00 + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk0cs, ntk0 );
  }
  else if ( x <= 2.0E+00 )
  {
    y = x * x;
    value = - log ( 0.5E+00 * x ) * r4_besi0 ( x ) 
      - 0.25E+00 + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk0cs, ntk0 );
  }
  else if ( x <= xmax )
  {
    value = exp ( - x ) * r4_besk0e ( x );
  }
  else
  {
    value = 0.0E+00;
  }

  return value;
}
//****************************************************************************80

float r4_besk0e ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESK0E evaluates the exponentially scaled Bessel function K0(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESK0E, the exponentially scaled Bessel 
//    function K0(X).
//
{
  static float ak02cs[14] = {
   -0.01201869826307592E+00,
   -0.00917485269102569E+00,
   +0.00014445509317750E+00,
   -0.00000401361417543E+00,
   +0.00000015678318108E+00,
   -0.00000000777011043E+00,
   +0.00000000046111825E+00,
   -0.00000000003158592E+00,
   +0.00000000000243501E+00,
   -0.00000000000020743E+00,
   +0.00000000000001925E+00,
   -0.00000000000000192E+00,
   +0.00000000000000020E+00,
   -0.00000000000000002E+00 };
  static float ak0cs[17] = {
   -0.07643947903327941E+00,
   -0.02235652605699819E+00,
   +0.00077341811546938E+00,
   -0.00004281006688886E+00,
   +0.00000308170017386E+00,
   -0.00000026393672220E+00,
   +0.00000002563713036E+00,
   -0.00000000274270554E+00,
   +0.00000000031694296E+00,
   -0.00000000003902353E+00,
   +0.00000000000506804E+00,
   -0.00000000000068895E+00,
   +0.00000000000009744E+00,
   -0.00000000000001427E+00,
   +0.00000000000000215E+00,
   -0.00000000000000033E+00,
   +0.00000000000000005E+00 };
  static float bk0cs[11] = {
   -0.03532739323390276872E+00,
   +0.3442898999246284869E+00,
   +0.03597993651536150163E+00,
   +0.00126461541144692592E+00,
   +0.00002286212103119451E+00,
   +0.00000025347910790261E+00,
   +0.00000000190451637722E+00,
   +0.00000000001034969525E+00,
   +0.00000000000004259816E+00,
   +0.00000000000000013744E+00,
   +0.00000000000000000035E+00 };
  static int ntak0 = 0;
  static int ntak02 = 0;
  static int ntk0 = 0;
  float value;
  static float xsml = 0.0;
  float y;

  if ( ntk0 == 0 )
  {
    ntk0 = r4_inits ( bk0cs, 11, 0.1E+00 * r4_mach ( 3 ) );
    ntak0 = r4_inits ( ak0cs, 17, 0.1E+00 * r4_mach ( 3 ) );
    ntak02 = r4_inits ( ak02cs, 14, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_BESK0E - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0E+00;
    value = exp ( x ) * ( - log ( 0.5E+00 * x ) * r4_besi0 ( x ) - 0.25E+00 
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk0cs, ntk0 ) );
  }
  else if ( x <= 2.0E+00 )
  {
    y = x * x;
    value = exp ( x ) * ( - log ( 0.5E+00 * x ) * r4_besi0 ( x ) - 0.25E+00 
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk0cs, ntk0 ) );
  }
  else if ( x <= 8.0E+00 )
  {
    value = ( 1.25E+00 + r4_csevl ( ( 16.0E+00 / x - 5.0E+00 ) / 3.0E+00, 
      ak0cs, ntak0 ) ) / sqrt ( x );
  }
  else
  {
    value = ( 1.25E+00 
      + r4_csevl ( 16.0E+00 / x - 1.0E+00, ak02cs, ntak02 ) ) / sqrt ( x );
  }
  return value;
}
//****************************************************************************80

float r4_besk1 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESK1 evaluates the Bessel function K of order 1 of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESK1, the Bessel function K of order 1 of X.
//
{
  static float bk1cs[11] = {
    0.0253002273389477705E+00,
   -0.353155960776544876E+00,
   -0.122611180822657148E+00,
   -0.0069757238596398643E+00,
   -0.0001730288957513052E+00,
   -0.0000024334061415659E+00,
   -0.0000000221338763073E+00,
   -0.0000000001411488392E+00,
   -0.0000000000006666901E+00,
   -0.0000000000000024274E+00,
   -0.0000000000000000070E+00 };
  static int ntk1 = 0;
  float value;
  static float xmax = 0.0;
  static float xmin = 0.0;
  static float xsml = 0.0;
  float y;

  if ( ntk1 == 0 )
  {
    ntk1 = r4_inits ( bk1cs, 11, 0.1E+00 * r4_mach ( 3 ) );
    xmin = exp ( r4_max ( log ( r4_mach ( 1 ) ), 
      - log ( r4_mach ( 2 ) ) ) + 0.01E+00 );
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) );
    xmax = - log ( r4_mach ( 1 ) );
    xmax = xmax - 0.5E+00 * xmax * log ( xmax )
      / ( xmax + 0.5E+00 ) - 0.01E+00;
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_BESK1 = Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0E+00;
    value = log ( 0.5E+00 * x ) * r4_besi1 ( x ) + ( 0.75E+00 
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk1cs, ntk1 ) ) / x;
  }
  else if ( x <= 2.0E+00 )
  {
    y = x * x;
    value = log ( 0.5E+00 * x ) * r4_besi1 ( x ) + ( 0.75E+00 
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk1cs, ntk1 ) ) / x;
  }
  else if ( x <= xmax )
  {
    value = exp ( - x ) * r4_besk1e ( x );
  }
  else
  {
    value = 0.0E+00;
  }
  return value;
}
//****************************************************************************80

float r4_besk1e ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESK1E evaluates the exponentially scaled Bessel function K1(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESK1E, the exponentially scaled Bessel 
//    function K1(X).
//
{
  static float ak12cs[14] = {
   +0.06379308343739001E+00,
   +0.02832887813049721E+00,
   -0.00024753706739052E+00,
   +0.00000577197245160E+00,
   -0.00000020689392195E+00,
   +0.00000000973998344E+00,
   -0.00000000055853361E+00,
   +0.00000000003732996E+00,
   -0.00000000000282505E+00,
   +0.00000000000023720E+00,
   -0.00000000000002176E+00,
   +0.00000000000000215E+00,
   -0.00000000000000022E+00,
   +0.00000000000000002E+00 };
  static float ak1cs[17] = {
   +0.2744313406973883E+00,
   +0.07571989953199368E+00,
   -0.00144105155647540E+00,
   +0.00006650116955125E+00,
   -0.00000436998470952E+00,
   +0.00000035402774997E+00,
   -0.00000003311163779E+00,
   +0.00000000344597758E+00,
   -0.00000000038989323E+00,
   +0.00000000004720819E+00,
   -0.00000000000604783E+00,
   +0.00000000000081284E+00,
   -0.00000000000011386E+00,
   +0.00000000000001654E+00,
   -0.00000000000000248E+00,
   +0.00000000000000038E+00,
   -0.00000000000000006E+00 };
  static float bk1cs[11] = {
   +0.0253002273389477705E+00,
   -0.353155960776544876E+00,
   -0.122611180822657148E+00,
   -0.0069757238596398643E+00,
   -0.0001730288957513052E+00,
   -0.0000024334061415659E+00,
   -0.0000000221338763073E+00,
   -0.0000000001411488392E+00,
   -0.0000000000006666901E+00,
   -0.0000000000000024274E+00,
   -0.0000000000000000070E+00 };
  static int ntak1 = 0;
  static int ntak12 = 0;
  static int ntk1 = 0;
  float value;
  static float xmin = 0.0;
  static float xsml = 0.0;
  float y;

  if ( ntk1 == 0 )
  {
    ntk1 = r4_inits ( bk1cs, 11, 0.1E+00 * r4_mach ( 3 ) );
    ntak1 = r4_inits ( ak1cs, 17, 0.1E+00 * r4_mach ( 3 ) );
    ntak12 = r4_inits ( ak12cs, 14, 0.1E+00 * r4_mach ( 3 ) );
    xmin = exp ( r4_max ( log ( r4_mach ( 1 ) ), 
      - log ( r4_mach ( 2 ) ) ) + 0.01E+00 );
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_BESK1E = Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0E+00;
    value = exp ( x ) * ( log ( 0.5E+00 * x ) * r4_besi1 ( x ) 
      + ( 0.75E+00 
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk1cs, ntk1 ) ) / x );
  }
  else if ( x <= 2.0E+00 )
  {
    y = x * x;
    value = exp ( x ) * ( log ( 0.5E+00 * x ) * r4_besi1 ( x ) + ( 0.75E+00 
      + r4_csevl ( 0.5E+00 * y - 1.0E+00, bk1cs, ntk1 ) ) / x );
  }
  else if ( x <= 8.0E+00 )
  {
    value = ( 1.25E+00 + r4_csevl ( ( 16.0E+00 / x - 5.0E+00 ) / 3.0E+00, 
      ak1cs, ntak1 ) ) / sqrt ( x );
  }
  else
  {
    value = ( 1.25E+00 
      + r4_csevl ( 16.E+00 / x - 1.0E+00, ak12cs, ntak12 ) ) / sqrt ( x );
  }
  return value;
}
//****************************************************************************80

float *r4_beskes ( float xnu, float x, int nin )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESKES evaluates a sequence of exponentially scaled K Bessel functions at X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float XNU, the order of the first function.
//    |XNU| < 1.
//
//    Input, float X, the argument.
//
//    Input, int NIN, the absolute value of NIN indicates the number of terms to compute.
//    If NIN < 0, successive values of NU count DOWN from XNU.
//    If NIN > 0, successive values of NU count UP from XNU.
//
//    Output, float R4_BESKES[abs(NIN)], the exponentially scaled K Bessel functions.
//
{
  static float alnbig = 0.0;
  float *bke;
  float bknu1;
  float direct;
  int i;
  int iswtch;
  int n;
  float v;
  float vend;
  float vincr;

  if ( alnbig == 0.0E+00 )
  {
    alnbig = log ( r4_mach ( 2 ) );
  }

  v = r4_abs ( xnu );
  n = abs ( nin );

  if ( 1.0E+00 <= v )
  {
    cerr << "\n";
    cerr << "R4_BESKES - Fatal error!\n";
    cerr << "  |XNU| must be less than 1.\n";
    exit ( 1 );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_BESKES - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    cerr << "\n";
    cerr << "R4_BESKES - Fatal error!\n";
    cerr << "  N = 0.\n";
    exit ( 1 );
  }

  bke = new float[abs(nin)];

  r4_knus ( v, x, bke[0], bknu1, iswtch );

  if ( n == 1 )
  {
    return bke;
  }

  if ( nin < 0 )
  {
    vincr = - 1.0E+00;
  }
  else
  {
    vincr = + 1.0E+00;
  }

  if ( xnu < 0.0E+00 )
  {
    direct = - vincr;
  }
  else
  {
    direct = vincr;
  }

  bke[1] = bknu1;

  if ( direct < 0.0E+00 )
  {
    r4_knus ( r4_abs ( xnu + vincr ), x, bke[1], bknu1, iswtch );
  }

  if ( n == 2 )
  {
    return bke;
  }

  vend = r4_abs ( xnu + ( float ) ( nin ) ) - 1.0E+00;

  v = xnu;
  for ( i = 2; i < n; i++ )
  {
    v = v + vincr;
    bke[i] = 2.0E+00 * v * bke[i-1] / x + bke[i-2];
  }
  return bke;
}
//****************************************************************************80

float *r4_besks ( float xnu, float x, int nin )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESKS evaluates a sequence of K Bessel functions at X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float XNU, the order of the first function.
//    |XNU| < 1.
//
//    Input, float X, the argument.
//
//    Input, int NIN, the absolute value of NIN indicates the number of terms to compute.
//    If NIN < 0, successive values of NU count DOWN from XNU.
//    If NIN > 0, successive values of NU count UP from XNU.
//
//    Output, float R4_BESKS[abs(NIN)], the K Bessel functions.
//
{
  float *bk;
  float expxi;
  int i;
  int n;
  static float xmax = 0.0;

  if ( xmax == 0.0E+00 )
  {
    xmax = - log ( r4_mach ( 1 ) );
    xmax = xmax + 0.5E+00 * log ( 3.14E+00 * 0.5E+00 / xmax );
  }

  bk = r4_beskes ( xnu, x, nin );

  expxi = exp ( - x );
  n = abs ( nin );

  for ( i = 0; i < n; i++ )
  {
    bk[i] = expxi * bk[i];
  }

  return bk;
}
//****************************************************************************80

float r4_besy0 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESY0 evaluates the Bessel function Y of order 0 of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESY0, the Bessel function Y of order 0 of X.
//
{
  static float alnhaf = -0.693147180559945309E+00;
  float ampl;
  static float bm0cs[21] = {
   +0.09284961637381644E+00,
   -0.00142987707403484E+00,
   +0.00002830579271257E+00,
   -0.00000143300611424E+00,
   +0.00000012028628046E+00,
   -0.00000001397113013E+00,
   +0.00000000204076188E+00,
   -0.00000000035399669E+00,
   +0.00000000007024759E+00,
   -0.00000000001554107E+00,
   +0.00000000000376226E+00,
   -0.00000000000098282E+00,
   +0.00000000000027408E+00,
   -0.00000000000008091E+00,
   +0.00000000000002511E+00,
   -0.00000000000000814E+00,
   +0.00000000000000275E+00,
   -0.00000000000000096E+00,
   +0.00000000000000034E+00,
   -0.00000000000000012E+00,
   +0.00000000000000004E+00 };
  static float bth0cs[24] = {
   -0.24639163774300119E+00,
   +0.001737098307508963E+00,
   -0.000062183633402968E+00,
   +0.000004368050165742E+00,
   -0.000000456093019869E+00,
   +0.000000062197400101E+00,
   -0.000000010300442889E+00,
   +0.000000001979526776E+00,
   -0.000000000428198396E+00,
   +0.000000000102035840E+00,
   -0.000000000026363898E+00,
   +0.000000000007297935E+00,
   -0.000000000002144188E+00,
   +0.000000000000663693E+00,
   -0.000000000000215126E+00,
   +0.000000000000072659E+00,
   -0.000000000000025465E+00,
   +0.000000000000009229E+00,
   -0.000000000000003448E+00,
   +0.000000000000001325E+00,
   -0.000000000000000522E+00,
   +0.000000000000000210E+00,
   -0.000000000000000087E+00,
   +0.000000000000000036E+00 };
  static float by0cs[13] = {
   -0.011277839392865573E+00,
   -0.12834523756042035E+00,
   -0.10437884799794249E+00,
   +0.023662749183969695E+00,
   -0.002090391647700486E+00,
   +0.000103975453939057E+00,
   -0.000003369747162423E+00,
   +0.000000077293842676E+00,
   -0.000000001324976772E+00,
   +0.000000000017648232E+00,
   -0.000000000000188105E+00,
   +0.000000000000001641E+00,
   -0.000000000000000011E+00 };
  static int ntm0 = 0;
  static int ntth0 = 0;
  static int nty0 = 0;
  static float pi4 = 0.78539816339744831E+00;
  float theta;
  static float twodpi = 0.63661977236758134E+00;
  float value;
  static float xmax = 0.0;
  static float xsml = 0.0;
  float y;
  float z;

  if ( nty0 == 0 )
  {
    nty0 = r4_inits ( by0cs, 13, 0.1E+00 * r4_mach ( 3 ) );
    ntm0 = r4_inits ( bm0cs, 21, 0.1E+00 * r4_mach ( 3 ) );
    ntth0 = r4_inits ( bth0cs, 24, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) );
    xmax = 1.0E+00 / r4_mach ( 4 );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_BESY0 - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0E+00;
    value = twodpi * ( alnhaf + log ( x ) ) 
      * r4_besj0 ( x ) + 0.375E+00 
      + r4_csevl ( 0.125E+00 * y - 1.0E+00, by0cs, nty0 );
  }
  else if ( x <= 4.0E+00 )
  {
    y = x * x;
    value = twodpi * ( alnhaf + log ( x ) ) 
      * r4_besj0 ( x ) + 0.375E+00 
      + r4_csevl ( 0.125E+00 * y - 1.0E+00, by0cs, nty0 );
  }
  else 
  {
    z = 32.0E+00 / x / x - 1.0E+00;
    ampl = ( 0.75E+00 + r4_csevl ( z, bm0cs, ntm0 ) ) / sqrt ( x );
    theta = x - pi4 + r4_csevl (z, bth0cs, ntth0 ) / x;
    value = ampl * sin ( theta );
  }
  return value;
}
//****************************************************************************80

float r4_besy1 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BESY1 evaluates the Bessel function Y of order 1 of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BESY1, the Bessel function Y of order 1 of X.
//
{
  static float ampl;
  static float bm1cs[21] = {
   +0.1047362510931285E+00,
   +0.00442443893702345E+00,
   -0.00005661639504035E+00,
   +0.00000231349417339E+00,
   -0.00000017377182007E+00,
   +0.00000001893209930E+00,
   -0.00000000265416023E+00,
   +0.00000000044740209E+00,
   -0.00000000008691795E+00,
   +0.00000000001891492E+00,
   -0.00000000000451884E+00,
   +0.00000000000116765E+00,
   -0.00000000000032265E+00,
   +0.00000000000009450E+00,
   -0.00000000000002913E+00,
   +0.00000000000000939E+00,
   -0.00000000000000315E+00,
   +0.00000000000000109E+00,
   -0.00000000000000039E+00,
   +0.00000000000000014E+00,
   -0.00000000000000005E+00 };
  static float bth1cs[24] = {
   +0.74060141026313850E+00,
   -0.004571755659637690E+00,
   +0.000119818510964326E+00,
   -0.000006964561891648E+00,
   +0.000000655495621447E+00,
   -0.000000084066228945E+00,
   +0.000000013376886564E+00,
   -0.000000002499565654E+00,
   +0.000000000529495100E+00,
   -0.000000000124135944E+00,
   +0.000000000031656485E+00,
   -0.000000000008668640E+00,
   +0.000000000002523758E+00,
   -0.000000000000775085E+00,
   +0.000000000000249527E+00,
   -0.000000000000083773E+00,
   +0.000000000000029205E+00,
   -0.000000000000010534E+00,
   +0.000000000000003919E+00,
   -0.000000000000001500E+00,
   +0.000000000000000589E+00,
   -0.000000000000000237E+00,
   +0.000000000000000097E+00,
   -0.000000000000000040E+00 };
  static float by1cs[14] = {
   +0.03208047100611908629E+00,
   +1.262707897433500450E+00,
   +0.00649996189992317500E+00,
   -0.08936164528860504117E+00,
   +0.01325088122175709545E+00,
   -0.00089790591196483523E+00,
   +0.00003647361487958306E+00,
   -0.00000100137438166600E+00,
   +0.00000001994539657390E+00,
   -0.00000000030230656018E+00,
   +0.00000000000360987815E+00,
   -0.00000000000003487488E+00,
   +0.00000000000000027838E+00,
   -0.00000000000000000186E+00 };
  static int ntm1 = 0;
  static int ntth1 = 0;
  static int nty1 = 0;
  static float pi4 = 0.78539816339744831E+00;
  float theta;
  static float twodpi = 0.63661977236758134E+00;
  float value;
  static float xmax = 0.0;
  static float xmin = 0.0;
  static float xsml = 0.0;
  float y;
  float z;

  if ( nty1 == 0 )
  {
    nty1 = r4_inits ( by1cs, 14, 0.1E+00 * r4_mach ( 3 ) );
    ntm1 = r4_inits ( bm1cs, 21, 0.1E+00 * r4_mach ( 3 ) );
    ntth1 = r4_inits ( bth1cs, 24, 0.1E+00 * r4_mach ( 3 ) );
    xmin = 1.571E+00 * exp ( r4_max ( log ( r4_mach ( 1 ) ), 
      - log ( r4_mach ( 2 ) ) ) + 0.01E+00 );
    xsml = sqrt ( 4.0E+00 * r4_mach ( 3 ) );
    xmax = 1.0E+00 / r4_mach ( 4 );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_BESY1 - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0E+00;
    value = twodpi * log ( 0.5E+00 * x ) * r4_besj1 ( x ) 
      + ( 0.5E+00 + r4_csevl ( 0.125E+00 * y - 1.0E+00, by1cs, 
      nty1 ) ) / x;
  }
  else if ( x <= 4.0E+00 )
  {
    y = x * x;
    value = twodpi * log ( 0.5E+00 * x ) * r4_besj1 ( x ) 
      + ( 0.5E+00 + r4_csevl ( 0.125E+00 * y - 1.0E+00, by1cs, 
      nty1 ) ) / x;
  }
  else
  {
    z = 32.0E+00 / x / x - 1.0E+00;
    ampl = ( 0.75E+00 + r4_csevl ( z, bm1cs, ntm1 ) ) / sqrt ( x );
    theta = x - 3.0E+00 * pi4 + r4_csevl ( z, bth1cs, ntth1 ) / x;
    value = ampl * sin ( theta );
  }
  return value;
}
//****************************************************************************80

float r4_beta ( float a, float b )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BETA evaluates the beta function of R4 arguments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, B, the arguments.
//
//    Output, float R4_BETA, the beta function of A and B.
//
{
  float alnsml = 0.0;
  float value;
  float xmax = 0.0;
  float xmin;

  if ( xmax == 0.0E+00 )
  {
    r4_gaml ( xmin, xmax );
    alnsml = log ( r4_mach ( 1 ) );
  }

  if ( a <= 0.0E+00 || b <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_BETA - Fatal error!\n";
    cerr << "  A and B must be greater than 0.\n";
    exit ( 1 );
  }

  if ( a + b < xmax )
  {
    value = r4_gamma ( a ) * r4_gamma ( b ) / r4_gamma ( a + b );
    return value;
  }

  value = r4_lbeta ( a, b );

  value = exp ( value );

  return value;
}
//****************************************************************************80

float r4_betai ( float x, float pin, float qin )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BETAI evaluates the incomplete beta ratio of R4 arguments.
//
//  Discussion:
//
//    The incomplete Beta function ratio is the probability that a
//    random variable from a beta distribution having parameters
//    P and Q will be less than or equal to X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Nancy Bosten, EL Battiste,
//    Remark on Algorithm 179: 
//    Incomplete Beta Ratio,
//    Communications of the ACM,
//    Volume 17, Number 3, March 1974, pages 156-157.
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the upper limit of integration.
//    0.0 <= X <= 1.0.
//
//    Input, float PIN, the first distribution parameter.
//    0.0 < PIN.
//
//    Input, float QIN, the second distribution parameter.
//    0.0 < QIN.
//
//    Output, float R4_BETAI, the incomplete beta function ratio.
//
{
  static float alneps = 0.0;
  static float alnsml = 0.0;
  float c;
  static float eps = 0.0;
  float finsum;
  int i;
  int ib;
  int n;
  float p;
  float p1;
  float ps;
  float q;
  static float sml = 0.0;
  float term;
  float value;
  float xb;
  float y;

  if ( eps == 0.0E+00 )
  {
    eps = r4_mach ( 3 );
    alneps = log ( eps );
    sml = r4_mach ( 1 );
    alnsml = log ( sml );
  }

  if ( x < 0.0E+00 || 1.0E+00 < x )
  {
    cerr << "\n";
    cerr << "R4_BETAI - Fatal error!\n";
    cerr << "  0 <= X <= 1 is required.\n";
    exit ( 1 );
  }

  if ( pin <= 0.0E+00 || qin <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_BETAI - Fatal error!\n";
    cerr << "  P or Q <= 0.0.\n";
    exit ( 1 );
  }

  y = x;
  p = pin;
  q = qin;

  if ( p < q || 0.8E+00 <= x )
  {
    if ( 0.2E+00 <= x )
    {
      y = 1.0E+00 - y;
      p = qin;
      q = pin;
    }
  }

  if ( ( p + q ) * y / ( p + 1.0E+00 ) < eps )
  {
    value = 0.0E+00;

    xb = p * log ( r4_max ( y, sml ) ) - log ( p ) - r4_lbeta ( p, q );

    if ( alnsml < xb && y != 0.0E+00 )
    {
      value = exp ( xb );
    }

    if ( y != x || p != pin )
    {
      value = 1.0E+00 - value;
    }
    return value;
  }
//
//  Evaluate the infinite sum first.
//  TERM will equal y**p/beta(ps,p) * (1.-ps)i * y**i / fac(i)
//
  ps = q - r4_aint ( q );
  if ( ps == 0.0E+00 ) 
  {
    ps = 1.0E+00;
  }

  xb = p * log ( y ) - r4_lbeta ( ps, p ) - log ( p );

  if ( xb < alnsml )
  {
    value = 0.0E+00;
  }
  else
  {
    value = exp ( xb );
    term = value * p;

    if ( ps != 1.0E+00 )
    {
      n = ( int ) ( r4_max ( alneps / log ( y ), 4.0E+00 ) );
      for ( i = 1; i <= n; i++ )
      {
        term = term * ( ( float ) ( i ) - ps ) * y / ( float ) ( i );
        value = value + term / ( p + ( float ) ( i ) );
      }
    }
  }
//
//  Now evaluate the finite sum.
//
  if  ( 1.0E+00 < q )
  {
    xb = p * log ( y ) + q * log ( 1.0E+00 - y ) 
      - r4_lbeta ( p, q ) - log ( q );
    ib = ( int ) ( r4_max ( xb / alnsml, 0.0E+00 ) );
    term = exp ( xb - ( float ) ( ib ) * alnsml );
    c = 1.0E+00 / ( 1.0E+00 - y );
    p1 = q * c / ( p + q - 1.0E+00 );

    finsum = 0.0E+00;
    n = ( int ) ( q );
    if ( q == ( float ) ( n ) )
    {
      n = n - 1;
    }

    for ( i = 1; i <= n; i++ )
    {
      if ( p1 <= 1.0E+00 && term / eps <= finsum )
      {
        break;
      }

      term = ( q - ( float ) ( i - 1 ) ) * c * term 
        / ( p + q - ( float ) ( i ) );

      if ( 1.0E+00 < term )
      {
        ib = ib - 1;
        term = term * sml;
      }

      if ( ib == 0 )
      {
        finsum = finsum + term;
      }
    }
    value = value + finsum;
  }

  if ( y != x || p != pin )
  {
    value = 1.0E+00 - value;
  }

  if ( value < 0.0E+00 )
  {
    value = 0.0E+00;
  }

  if ( 1.0E+00 < value )
  {
    value = 1.0E+00;
  }

  return value;
}
//****************************************************************************80

float r4_bi ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BI evaluates the Airy function Bi of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BI, the Airy function Bi of X.
//
{
  static float bif2cs[10] = {
    0.09984572693816041E+00,
    0.478624977863005538E+00,
    0.0251552119604330118E+00,
    0.0005820693885232645E+00,
    0.0000074997659644377E+00,
    0.0000000613460287034E+00,
    0.0000000003462753885E+00,
    0.0000000000014288910E+00,
    0.0000000000000044962E+00,
    0.0000000000000000111E+00 };
  static float bifcs[9] = {
   -0.01673021647198664948E+00,
    0.1025233583424944561E+00,
    0.00170830925073815165E+00,
    0.00001186254546774468E+00,
    0.00000004493290701779E+00,
    0.00000000010698207143E+00,
    0.00000000000017480643E+00,
    0.00000000000000020810E+00,
    0.00000000000000000018E+00 };
  static float big2cs[10] = {
    0.033305662145514340E+00,
    0.161309215123197068E+00,
    0.0063190073096134286E+00,
    0.0001187904568162517E+00,
    0.0000013045345886200E+00,
    0.0000000093741259955E+00,
    0.0000000000474580188E+00,
    0.0000000000001783107E+00,
    0.0000000000000005167E+00,
    0.0000000000000000011E+00 };
  static float bigcs[8] = {
    0.02246622324857452E+00,
    0.03736477545301955E+00,
    0.00044476218957212E+00,
    0.00000247080756363E+00,
    0.00000000791913533E+00,
    0.00000000001649807E+00,
    0.00000000000002411E+00,
    0.00000000000000002E+00 };
  float eta;
  static int nbif = 0;
  static int nbif2 = 0;
  static int nbig = 0;
  static int nbig2 = 0;
  float theta;
  float value;
  static float x3sml = 0.0;
  float xm;
  static float xmax = 0.0;
  float z;

  if ( nbif == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    nbif = r4_inits ( bifcs, 9, eta );
    nbig = r4_inits ( bigcs, 8, eta );
    nbif2 = r4_inits ( bif2cs, 10, eta );
    nbig2 = r4_inits ( big2cs, 10, eta );
    x3sml = r4_power ( eta, 0.3333E+00 );
    xmax = r4_power ( 1.5E+00 * log ( r4_mach ( 2 ) ), 0.6666E+00 );
  }

  if ( x <= - 1.0E+00 )
  {
    r4_aimp ( x, xm, theta );
    value = xm * sin ( theta );
  }
  else if ( r4_abs ( x ) <= x3sml )
  {
    z = 0.0E+00;

    value = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) 
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) );
  }
  else if ( x <= 1.0E+00 )
  {
    z = x * x * x;
    value = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) 
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) );
  }
  else if ( x <= 2.0E+00 )
  {
    z = ( 2.0E+00 * x * x * x - 9.0E+00 ) / 7.0E+00;
    value = 1.125E+00 + r4_csevl ( z, bif2cs, nbif2 ) 
      + x * ( 0.625E+00 + r4_csevl ( z, big2cs, nbig2 ) );
  }
  else
  {
    value = r4_bie ( x ) * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  return value;
}
//****************************************************************************80

float r4_bid ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BID evaluates the derivative of the Airy function Bi of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BID, the derivative of the Airy function Bi of X.
//
{
  static float bif2cs[10] = {
     0.323493987603522033521E+00,
     0.086297871535563559139E+00,
     0.002994025552655397426E+00,
     0.000051430528364661637E+00,
     0.000000525840250036811E+00,
     0.000000003561751373958E+00,
     0.000000000017146864007E+00,
     0.000000000000061663520E+00,
     0.000000000000000171911E+00,
     0.000000000000000000382E+00 };
  static float bifcs[8] = {
     0.1153536790828570243E+00,
     0.0205007894049192875E+00,
     0.0002135290278902876E+00,
     0.0000010783960614677E+00,
     0.0000000032094708833E+00,
     0.0000000000062930407E+00,
     0.0000000000000087403E+00,
     0.0000000000000000090E+00 };
  static float big2cs[10] = {
     1.6062999463621294578E+00,
     0.7449088819876088652E+00,
     0.0470138738610277380E+00,
     0.0012284422062548239E+00,
     0.0000173222412256624E+00,
     0.0000001521901652368E+00,
     0.0000000009113560249E+00,
     0.0000000000039547918E+00,
     0.0000000000000130017E+00,
     0.0000000000000000335E+00 };
  static float bigcs[9] = {
    -0.097196440416443537390E+00,
     0.149503576843167066571E+00,
     0.003113525387121326042E+00,
     0.000024708570579821297E+00,
     0.000000102949627731379E+00,
     0.000000000263970373987E+00,
     0.000000000000458279271E+00,
     0.000000000000000574283E+00,
     0.000000000000000000544E+00 };
  float eta;
  static int nbif = 0;
  static int nbif2 = 0;
  static int nbig = 0;
  static int nbig2 = 0;
  float phi;
  float value;
  float x2;
  static float x2sml = 0.0;
  float x3;
  static float x3sml = 0.0;
  static float xmax = 0.0;
  float xn;
  float z;

  if ( nbif == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    nbif = r4_inits ( bifcs, 8, eta );
    nbig = r4_inits ( bigcs, 9, eta );
    nbif2 = r4_inits ( bif2cs, 10, eta );
    nbig2 = r4_inits ( big2cs, 10, eta );
    x2sml = sqrt ( eta );
    x3sml = r4_power ( eta, 0.3333E+00 );
    xmax = r4_power ( 1.5E+00 * log ( r4_mach ( 2 ) ), 0.6666E+00 );
  }

  if ( x < - 1.0E+00 )
  {
    r4_admp ( x, xn, phi );
    value = xn * sin ( phi );
  }
  else if ( r4_abs ( x ) <= x2sml )
  {
    x2 = 0.0E+00;
    x3 = 0.0E+00;
    value = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) 
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00;
  }
  else if ( r4_abs ( x ) <= x3sml )
  {
    x2 = x * x;
    x3 = 0.0E+00;
    value = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) 
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00;
  }
  else if ( x <= 1.0E+00 )
  {
    x2 = x * x;
    x3 = x * x * x;
    value = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) 
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00;
  }
  else if ( x <= 2.0E+00 )
  {
    z = ( 2.0E+00 * x * x * x - 9.0E+00 ) / 7.0E+00;
    value = x * x * ( r4_csevl ( z, bif2cs, nbif2 ) + 0.25E+00 ) 
      + r4_csevl ( z, big2cs, nbig2 ) + 0.5E+00;
  }
  else
  {
    value = r4_bide ( x )
      * exp ( 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  return value;
}
//****************************************************************************80

float r4_bide ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BIDE: exponentially scaled derivative, Airy function Bi of an R4 argument.
//
//  Discussion:
//
//    if X < 0,
//      R4_BIDE ( X ) = R4_BID ( X )
//    else
//      R4_BIDE ( X ) = R4_BID ( X ) * exp ( - 2/3 * X**(3/2) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BIDE, the exponentially scaled derivative of 
//    the Airy function Bi of X.
//
{
  static float atr = 8.7506905708484345E+00;
  static float bif2cs[10] = {
     0.323493987603522033521E+00,
     0.086297871535563559139E+00,
     0.002994025552655397426E+00,
     0.000051430528364661637E+00,
     0.000000525840250036811E+00,
     0.000000003561751373958E+00,
     0.000000000017146864007E+00,
     0.000000000000061663520E+00,
     0.000000000000000171911E+00,
     0.000000000000000000382E+00 };
  static float bifcs[8] = {
     0.1153536790828570243E+00,
     0.0205007894049192875E+00,
     0.0002135290278902876E+00,
     0.0000010783960614677E+00,
     0.0000000032094708833E+00,
     0.0000000000062930407E+00,
     0.0000000000000087403E+00,
     0.0000000000000000090E+00 };
  static float big2cs[10] = {
     1.6062999463621294578E+00,
     0.7449088819876088652E+00,
     0.0470138738610277380E+00,
     0.0012284422062548239E+00,
     0.0000173222412256624E+00,
     0.0000001521901652368E+00,
     0.0000000009113560249E+00,
     0.0000000000039547918E+00,
     0.0000000000000130017E+00,
     0.0000000000000000335E+00 };
  static float bigcs[9] = {
    -0.097196440416443537390E+00,
     0.149503576843167066571E+00,
     0.003113525387121326042E+00,
     0.000024708570579821297E+00,
     0.000000102949627731379E+00,
     0.000000000263970373987E+00,
     0.000000000000458279271E+00,
     0.000000000000000574283E+00,
     0.000000000000000000544E+00 };
  static float bip1cs[24] = {
    -0.1729187351079553719E+00,
    -0.0149358492984694364E+00,
    -0.0005471104951678566E+00,
     0.0001537966292958408E+00,
     0.0000154353476192179E+00,
    -0.0000065434113851906E+00,
     0.0000003728082407879E+00,
     0.0000002072078388189E+00,
    -0.0000000658173336470E+00,
     0.0000000074926746354E+00,
     0.0000000011101336884E+00,
    -0.0000000007265140553E+00,
     0.0000000001782723560E+00,
    -0.0000000000217346352E+00,
    -0.0000000000020302035E+00,
     0.0000000000019311827E+00,
    -0.0000000000006044953E+00,
     0.0000000000001209450E+00,
    -0.0000000000000125109E+00,
    -0.0000000000000019917E+00,
     0.0000000000000015154E+00,
    -0.0000000000000004977E+00,
     0.0000000000000001155E+00,
    -0.0000000000000000186E+00 };
  static float bip2cs[29] = {
    -0.13269705443526630495E+00,
    -0.00568443626045977481E+00,
    -0.00015643601119611610E+00,
    -0.00001136737203679562E+00,
    -0.00000143464350991284E+00,
    -0.00000018098531185164E+00,
     0.00000000926177343611E+00,
     0.00000001710005490721E+00,
     0.00000000476698163504E+00,
    -0.00000000035195022023E+00,
    -0.00000000058890614316E+00,
    -0.00000000006678499608E+00,
     0.00000000006395565102E+00,
     0.00000000001554529427E+00,
    -0.00000000000792397000E+00,
    -0.00000000000258326243E+00,
     0.00000000000121655048E+00,
     0.00000000000038707207E+00,
    -0.00000000000022487045E+00,
    -0.00000000000004953477E+00,
     0.00000000000004563782E+00,
     0.00000000000000332998E+00,
    -0.00000000000000921750E+00,
     0.00000000000000094157E+00,
     0.00000000000000167154E+00,
    -0.00000000000000055134E+00,
    -0.00000000000000022369E+00,
     0.00000000000000017487E+00,
     0.00000000000000000207E+00 };
  static float btr = -2.0938363213560543E+00;
  float eta;
  static int nbif = 0;
  static int nbif2 = 0;
  static int nbig = 0;
  static int nbig2 = 0;
  static int nbip1 = 0;
  static int nbip2 = 0;
  float phi;
  float sqrtx;
  float value;
  float x2;
  static float x2sml = 0.0;
  float x3;
  static float x3sml = 0.0;
  static float x32sml = 0.0;
  static float xbig = 0.0;
  float xn;
  float z;

  if ( nbif == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    nbif = r4_inits ( bifcs, 8, eta );
    nbig = r4_inits ( bigcs, 9, eta );
    nbif2 = r4_inits ( bif2cs, 10, eta );
    nbig2 = r4_inits ( big2cs, 10, eta );
    nbip1 = r4_inits ( bip1cs, 24, eta );
    nbip2 = r4_inits ( bip2cs, 29, eta );
    x2sml = sqrt ( eta );
    x3sml = r4_power ( eta, 0.3333E+00 );
    x32sml = 1.3104 * x3sml * x3sml;
    xbig = r4_power ( r4_mach ( 2 ), 0.6666E+00 );
  }

  if ( x <= - 1.0E+00 )
  {
    r4_admp ( x, xn, phi );
    value = xn * sin ( phi );
  }
  else if ( 0.0E+00 <= x && x <= x32sml )
  {
    x2 = 0.0E+00;
    x3 = 0.0E+00;
    value = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) 
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00;
  }
  else if ( r4_abs ( x ) <= x2sml )
  {
    x2 = 0.0E+00;
    x3 = 0.0E+00;
    value = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) 
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00;
    value = value * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  else if ( x <= x3sml )
  {
    x2 = x * x;
    x3 = 0.0E+00;
    value = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) 
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00;
    value = value * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  else if ( x <= 1.0E+00 )
  {
    x2 = x * x;
    x3 = x * x * x;
    value = x2 * ( r4_csevl ( x3, bifcs, nbif ) + 0.25E+00 ) 
      + r4_csevl ( x3, bigcs, nbig ) + 0.5E+00;
    value = value * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  else if ( x <= 2.0E+00 )
  {
    z = ( 2.0E+00 * x * x * x - 9.0E+00 ) / 7.0E+00;
    value = exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 ) 
      * ( x * x * ( 0.25E+00 + r4_csevl ( z, bif2cs, nbif2 ) ) 
      + 0.5E+00 + r4_csevl ( z, big2cs, nbig2 ) );
  }
  else if ( x <= 4.0E+00 )
  {
    sqrtx = sqrt ( x );
    z = atr / ( x * sqrtx ) + btr;
    value = ( 0.625E+00 + r4_csevl ( z, bip1cs, nbip1 ) ) * sqrt ( sqrtx );
  }
  else if ( x < xbig )
  {
    sqrtx = sqrt ( x );
    z = 16.0E+00 / ( x * sqrtx ) - 1.0E+00;
    value = ( 0.625E+00 + r4_csevl ( z, bip2cs, nbip2 ) ) * sqrt ( sqrtx );
  }
  else
  {
    sqrtx = sqrt ( x );
    z = - 1.0E+00;
    value = ( 0.625E+00 + r4_csevl ( z, bip2cs, nbip2 ) ) * sqrt ( sqrtx );
  }
  return value;
}
//****************************************************************************80

float r4_bie ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BIE evaluates the exponentially scaled Airy function Bi of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_BIE, the exponentially scaled Airy function Bi of X.
//
{
  static float atr = 8.7506905708484345E+00;
  static float bif2cs[10] = {
    0.09984572693816041E+00,
    0.478624977863005538E+00,
    0.0251552119604330118E+00,
    0.0005820693885232645E+00,
    0.0000074997659644377E+00,
    0.0000000613460287034E+00,
    0.0000000003462753885E+00,
    0.0000000000014288910E+00,
    0.0000000000000044962E+00,
    0.0000000000000000111E+00 };
  static float bifcs[9] = {
   -0.01673021647198664948E+00,
    0.1025233583424944561E+00,
    0.00170830925073815165E+00,
    0.00001186254546774468E+00,
    0.00000004493290701779E+00,
    0.00000000010698207143E+00,
    0.00000000000017480643E+00,
    0.00000000000000020810E+00,
    0.00000000000000000018E+00 };
  static float big2cs[10] = {
    0.033305662145514340E+00,
    0.161309215123197068E+00,
    0.0063190073096134286E+00,
    0.0001187904568162517E+00,
    0.0000013045345886200E+00,
    0.0000000093741259955E+00,
    0.0000000000474580188E+00,
    0.0000000000001783107E+00,
    0.0000000000000005167E+00,
    0.0000000000000000011E+00 };
  static float bigcs[8] = {
    0.02246622324857452E+00,
    0.03736477545301955E+00,
    0.00044476218957212E+00,
    0.00000247080756363E+00,
    0.00000000791913533E+00,
    0.00000000001649807E+00,
    0.00000000000002411E+00,
    0.00000000000000002E+00 };
  static float bip2cs[29] = {
   -0.113596737585988679E+00,
    0.0041381473947881595E+00,
    0.0001353470622119332E+00,
    0.0000104273166530153E+00,
    0.0000013474954767849E+00,
    0.0000001696537405438E+00,
   -0.0000000100965008656E+00,
   -0.0000000167291194937E+00,
   -0.0000000045815364485E+00,
    0.0000000003736681366E+00,
    0.0000000005766930320E+00,
    0.0000000000621812650E+00,
   -0.0000000000632941202E+00,
   -0.0000000000149150479E+00,
    0.0000000000078896213E+00,
    0.0000000000024960513E+00,
   -0.0000000000012130075E+00,
   -0.0000000000003740493E+00,
    0.0000000000002237727E+00,
    0.0000000000000474902E+00,
   -0.0000000000000452616E+00,
   -0.0000000000000030172E+00,
    0.0000000000000091058E+00,
   -0.0000000000000009814E+00,
   -0.0000000000000016429E+00,
    0.0000000000000005533E+00,
    0.0000000000000002175E+00,
   -0.0000000000000001737E+00,
   -0.0000000000000000010E+00 };
  static float bipcs[24] = {
   -0.08322047477943447E+00,
    0.01146118927371174E+00,
    0.00042896440718911E+00,
   -0.00014906639379950E+00,
   -0.00001307659726787E+00,
    0.00000632759839610E+00,
   -0.00000042226696982E+00,
   -0.00000019147186298E+00,
    0.00000006453106284E+00,
   -0.00000000784485467E+00,
   -0.00000000096077216E+00,
    0.00000000070004713E+00,
   -0.00000000017731789E+00,
    0.00000000002272089E+00,
    0.00000000000165404E+00,
   -0.00000000000185171E+00,
    0.00000000000059576E+00,
   -0.00000000000012194E+00,
    0.00000000000001334E+00,
    0.00000000000000172E+00,
   -0.00000000000000145E+00,
    0.00000000000000049E+00,
   -0.00000000000000011E+00,
    0.00000000000000001E+00 };
  static float btr = -2.093836321356054E+00;
  float eta;
  static int nbif;
  static int nbif2;
  static int nbig;
  static int nbig2;
  static int nbip;
  static int nbip2;
  float sqrtx;
  float theta;
  float value;
  static float x32sml = 0.0;
  static float x3sml = 0.0;
  static float xbig = 0.0;
  float xm;
  float z;

  if ( nbif == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    nbif = r4_inits ( bifcs, 9, eta );
    nbig = r4_inits ( bigcs, 8, eta );
    nbif2 = r4_inits ( bif2cs, 10, eta );
    nbig2 = r4_inits ( big2cs, 10, eta );
    nbip  = r4_inits ( bipcs, 24, eta );
    nbip2 = r4_inits ( bip2cs, 29, eta );
    x3sml = r4_power ( eta, 0.3333E+00 );
    x32sml = 1.3104E+00 * x3sml * x3sml;
    xbig = r4_power ( r4_mach ( 2 ), 0.6666E+00 );
  }

  if ( x < -1.0E+00 )
  {
    r4_aimp ( x, xm, theta );
    value = xm * sin ( theta );
  }
  else if ( r4_abs ( x ) <= x32sml )
  {
    z = 0.0E+00;
    value = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) 
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) );
  }
  else if ( r4_abs ( x ) <= x3sml )
  {
    z = 0.0E+00;
    value = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) 
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) );
    value = value * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  else if ( x <= 1.0E+00 )
  {
    z = x * x * x;
    value = 0.625E+00 + r4_csevl ( z, bifcs, nbif ) 
      + x * ( 0.4375E+00 + r4_csevl ( z, bigcs, nbig ) );
    value = value * exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 );
  }
  else if ( x <= 2.0E+00 )
  {
    z = ( 2.0E+00 * x * x * x - 9.0E+00 ) / 7.0E+00;
    value = exp ( - 2.0E+00 * x * sqrt ( x ) / 3.0E+00 ) 
      * ( 1.125E+00 + r4_csevl ( z, bif2cs, nbif2 )
      + x * ( 0.625E+00 + r4_csevl ( z, big2cs, nbig2 ) ) );
  }
  else if ( x <= 4.0E+00 )
  {
    sqrtx = sqrt ( x );
    z = atr / ( x * sqrtx ) + btr;
    value = ( 0.625E+00 + r4_csevl ( z, bipcs, nbip ) ) / sqrt ( sqrtx );
  }
  else if ( x <= xbig )
  {
    sqrtx = sqrt ( x );
    z = 16.0E+00 / ( x * sqrtx ) - 1.0E+00;
    value = ( 0.625E+00 + r4_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx );
  }
  else
  {
    sqrtx = sqrt ( x );
    z = - 1.0E+00;
    value = ( 0.625E+00 + r4_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx );
  }
  return value;
}
//****************************************************************************80

float r4_binom ( int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    R4_BINOM evaluates the binomial coefficient using R4 arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, int N, M, the arguments.
//
//    Output, float R4_BINOM, the binomial coefficient.
//
{
  static float bilnmx = 0.0;
  float corr;
  static float fintmx = 0.0;
  int i;
  int k;
  static float sq2pil = 0.91893853320467274E+00;
  float value;
  float xk;
  float xn;
  float xnk;

  if ( bilnmx == 0.0E+00 );
  {
    bilnmx = log ( r4_mach ( 2 ) );
    fintmx = 0.9E+00 / r4_mach ( 3 );
  }

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "R4_BINOM - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  if ( m < 0 )
  {
    cerr << "\n";
    cerr << "R4_BINOM - Fatal error!\n";
    cerr << "  M < 0.\n";
    exit ( 1 );
  }

  if ( n < m )
  {
    cerr << "\n";
    cerr << "R4_BINOM - Fatal error!\n";
    cerr << "  N < M.\n";
    exit ( 1 );
  }

  k = i4_min ( m, n - m );

  if ( k <= 20 &&  ( float ) ( k ) * log ( ( float ) ( i4_max ( n, 1 ) ) ) <= bilnmx )
  {
    value = 1.0E+00;
    for ( i = 1; i <= k; i++ )
    {
      value = value * ( float ) ( n - i + 1 ) / ( float ) ( i );
    }
  }
  else
  {
    if ( k < 9 )
    {
      cerr << "\n";
      cerr << "R4_BINOM - Fatal error!\n";
      cerr << "  Result overflows.\n";
      cerr << "  N or M is too big.\n";
      exit ( 1 );
    }

    xn = ( float ) ( n + 1 );
    xk = ( float ) ( k + 1 );
    xnk = ( float ) ( n - k + 1 );

    corr = r4_lgmc ( xn ) - r4_lgmc ( xk ) - r4_lgmc ( xnk );

    value = xk * log ( xnk / xk ) 
      - xn * r4_lnrel ( - ( xk - 1.0E+00 ) / xn )
      - 0.5E+00 * log ( xn * xnk / xk ) + 1.0E+00 - sq2pil + corr;

    if ( bilnmx < value )
    {
      cerr << "\n";
      cerr << "R4_BINOM - Fatal error!\n";
      cerr << "  Result overflows.\n";
      cerr << "  N or M is too big.\n";
      exit ( 1 );
    }
    value = exp ( value );
  }

  if ( value < fintmx )
  {
    value = r4_aint ( value + 0.5E+00 );
  }

  return value;
}
//****************************************************************************80

float r4_cbrt ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CBRT computes the cube root of an R4.
//
//  Discussion:
//
//    The approximation is a generalized Chebyshev series converted
//    to polynomial form.  The approximation is nearly best in the 
//    sense of relative error with 4.085 digits accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2012
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the number whose square root is desired.
//
//    Output, float R4_CBRT, the cube root of X.
//
{
  static float cbrt2[5] = {
     0.62996052494743658E+00,
     0.79370052598409974E+00,
     1.0E+00,
     1.25992104989487316E+00,
     1.58740105196819947E+00 };
  int irem;
  int iter;
  int ixpnt;
  int n;
  static int niter = 0;
  float value;
  float vsq;
  float y;

  if ( niter == 0 )
  {
    niter = ( int ) ( 1.443E+00 * log ( -0.106E+00 
      * log ( 0.1E+00 * r4_mach ( 3 ) ) ) + 1.0E+00 );
  }

  value = 0.0E+00;

  if ( x != 0.0E+00 )
  {
    r4_upak ( r4_abs ( x ), y, n );
    ixpnt = n / 3;
    irem = n - 3 * ixpnt + 3;

    value = 0.439581E+00 + y * ( 
            0.928549E+00 + y * (
          - 0.512653E+00 + y *
            0.144586E+00 ) );

    for ( iter = 1; iter <= niter; iter++ )
    {
      vsq = value * value;
      value = value + ( y - value * vsq ) / ( 3.0E+00 * vsq );
    }

    if ( x < 0.0E+00 )
    {
      value = - r4_abs ( value );
    }
    else
    {
      value = + r4_abs ( value );
    }
    value = r4_pak ( cbrt2[irem-1] * value, ixpnt );
  }
  return value;
}
//****************************************************************************80

float r4_chi ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CHI evaluates the hyperbolic cosine integral of an R4 argument.
//
//  Discussion:
//
//    The hyperbolic cosine integral is defined by
//
//      CHI(X) = gamma + log ( x ) 
//        + integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
//
//    where gamma is Euler's constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_CHI, the hyperbolic cosine integral
//    evaluated at X.
//
{
  float value;

  value = 0.5E+00 * ( r4_ei ( x ) - r4_e1 ( x ) );

  return value;
}
//****************************************************************************80

float r4_chu ( float a, float b, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CHU evaluates the confluent hypergeometric function of R4 arguments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, B, the parameters.
//
//    Input, float X, the argument.
//
//    Output, float R4_CHU, the function value.
//
{
  float a0;
  float aintb;
  float alnx;
  float b0;
  float beps;
  float c0;
  static float eps = 0.0;
  float factor;
  float gamri1;
  float gamrni;
  int i;
  int istrt;
  int m;
  int n;
  static float pi = 3.14159265358979324E+00;
  float pch1ai;
  float pch1i;
  float pochai;
  float sum;
  float t;
  float value;
  float xeps1;
  float xi;
  float xi1;
  float xn;
  float xtoeps;

  if ( eps == 0.0E+00 )
  {
    eps = r4_mach ( 3 );
  }

  if ( x < 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_CHU - Fatal error!\n";
    cerr << "  X < 0.\n";
    exit ( 1 );
  }

  if ( x == 0.0E+00 )
  {
    if ( 1.0E+00 <= b )
    {
      cerr << "\n";
      cerr << "R4_CHU - Fatal error!\n";
      cerr << "  X = 0 and 1 <= B.\n";
      exit ( 1 );
    }
    value = r4_gamma ( 1.0E+00 - b ) / r4_gamma ( 1.0E+00 + a - b );
    return value;
  }

  if ( r4_max ( r4_abs ( a ), 1.0E+00 ) 
    * r4_max ( r4_abs ( 1.0E+00 + a - b ), 1.0E+00 ) 
    < 0.99E+00 * r4_abs ( x ) )
  {
    value = r4_power ( x, - a ) * r4_chu_scaled ( a, b, x );
    return value;
  }
//
//  The ascending series will be used, because the descending rational
//  approximation (which is based on the asymptotic series) is unstable.
//
  if ( b < 0.0E+00 )
  {
    aintb = r4_aint ( b - 0.5E+00 );
  }
  else
  {
    aintb = r4_aint ( b + 0.5E+00 );
  }
  beps = b - aintb;
  n = aintb;

  alnx = log ( x );
  xtoeps = exp ( - beps * alnx );
//
//  Evaluate the finite sum.
//
//  Consider the case b < 1.0 first.
//
  if ( n < 1 )
  {
    sum = 1.0E+00;
    t = 1.0E+00;
    m = - n;
    for ( i = 1; i <= m; i++ )
    {
      xi1 = ( float ) ( i - 1 );
      t = t * ( a + xi1 ) * x / ( ( b + xi1 ) * ( xi1 + 1.0E+00 ) );
      sum = sum + t;
    }
    sum = r4_poch ( 1.0E+00 + a - b, - a ) * sum;
  }
//
//  Now consider the case 1.-0 <= b
//
  else
  {
    sum = 0.0E+00;
    m = n - 2;

    if ( 0 <= m )
    {
      t = 1.0E+00;
      sum = 1.0E+00;

      for ( i = 1; i <= m; i++ )
      {
        xi = ( float ) ( i );
        t = t * ( a - b + xi ) * x / ( ( 1.0E+00 - b + xi ) * xi );
        sum = sum + t;
      }

      sum = r4_gamma ( b - 1.0E+00 ) * r4_gamr ( a ) 
        * r4_power ( x, ( 1.0 - ( float ) n ) ) * xtoeps * sum;
    }
  }
//
//  Now evaluate the infinite sum.
//
  if ( n < 1 )
  {
    istrt = 1 - n;
  }
  else
  {
    istrt = 0;
  }

  xi = ( float ) ( istrt );

  factor = ( float ) i4_pow ( - 1, n ) * r4_gamr ( 1.0E+00 + a - b ) 
    * r4_power ( x, ( float ) istrt );

  if ( beps != 0.0E+00 )
  {
    factor = factor * beps * pi / sin ( beps * pi );
  }

  pochai = r4_poch ( a, xi );
  gamri1 = r4_gamr ( xi + 1.0E+00 );
  gamrni = r4_gamr ( aintb + xi );
  b0 = factor * r4_poch ( a, xi - beps ) * gamrni 
    * r4_gamr ( xi + 1.0E+00 - beps );
//
//  x^(-beps) is close to 1.0, so we must be careful in evaluating
//  the differences.
//
  if ( r4_abs ( xtoeps - 1.0E+00 ) <= 0.5E+00 )
  {
    pch1ai = r4_poch1 ( a + xi, - beps );
    pch1i = r4_poch1 ( xi + 1.0E+00 - beps, beps );
    c0 = factor * pochai * gamrni * gamri1 * (
      - r4_poch1 ( b + xi, -beps ) + pch1ai 
      - pch1i + beps * pch1ai * pch1i );
//
//  xeps1 = (1.0 - x^(-beps)) / beps
//
    xeps1 = alnx * r4_exprel ( - beps * alnx );

    value = sum + c0 + xeps1 * b0;
    xn = ( float ) ( n );

    for ( i = 1; i <= 1000; i++ )
    {
      xi = ( float ) ( istrt + i );
      xi1 = ( float ) ( istrt + i - 1 );
      b0 = ( a + xi1 - beps ) * b0 * x 
        / ( ( xn + xi1 ) * ( xi - beps ) );
      c0 = ( a + xi1 ) * c0 * x / ( ( b + xi1 ) * xi ) 
        - ( ( a - 1.0E+00 ) * ( xn + 2.0E+00 * xi - 1.0E+00 )
        + xi * ( xi - beps ) ) * b0 
        / ( xi * ( b + xi1 ) * ( a + xi1 - beps ) );
      t = c0 + xeps1 * b0;
      value = value + t;
      if ( r4_abs ( t ) < eps * r4_abs ( value ) )
      {
        return value;
      }
    }
    cerr << "\n";
    cerr << "R4_CHU - Fatal error!\n";
    cerr << "  No convergence in 1000 terms.\n";
    exit ( 1 );
  }
//
//  x^(-beps) is very different from 1.0, so the straightforward
//  formulation is stable.
//
  a0 = factor * pochai * r4_gamr ( b + xi ) * gamri1 / beps;
  b0 = xtoeps * b0 / beps;

  value = sum + a0 - b0;

  for ( i = 1; i <= 1000; i++ )
  {
    xi = ( float ) ( istrt + i );
    xi1 = ( float ) ( istrt + i - 1 );
    a0 = ( a + xi1 ) * a0 * x / ( ( b + xi1 ) * xi );
    b0 = ( a + xi1 - beps ) * b0 * x 
      / ( ( aintb + xi1 ) * ( xi - beps ) );
    t = a0 - b0;
    value = value + t;
    if ( r4_abs ( t ) < eps * r4_abs ( value ) )
    {
      return value;
    }
  }

  cerr << "\n";
  cerr << "R4_CHU - Fatal error!\n";
  cerr << "  No convergence in 1000 terms.\n";
  exit ( 1 );
}
//****************************************************************************80

float r4_chu_scaled ( float a, float b, float z )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CHU_SCALED: scaled confluent hypergeometric function of R4 arguments.
//
//  Discussion:
//
//    Evaluate, for large z, z**a * u(a,b,z)  where U is the logarithmic
//    confluent hypergeometric function.  A rational approximation due to
//    Y L Luke is used.  When U is not in the asymptotic region, that is, when A
//    or B is large compared with Z, considerable significance loss occurs.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, B, the parameters.
//
//    Input, float Z, the argument.
//
//    Output, float R4_CHU_SCALED, the function value.
//
{
  float aa[4];
  float ab;
  float anbn;
  float bb[4];
  float bp;
  float c2;
  float ct1;
  float ct2;
  float ct3;
  float d1z;
  static float eps = 0.0;
  float g1;
  float g2;
  float g3;
  int i;
  int j;
  float sab;
  static float sqeps = 0.0;
  float value;
  float x2i1;

  if ( eps == 0.0E+00 )
  {
    eps = 4.0E+00 * r4_mach ( 4 );
    sqeps = sqrt ( r4_mach ( 4 ) );
  }

  bp = 1.0E+00 + a - b;
  ab = a * bp;
  ct2 = 2.0E+00 * ( z - ab );
  sab = a + bp;

  bb[0] = 1.0E+00;
  aa[0] = 1.0E+00;

  ct3 = sab + 1.0E+00 + ab;
  bb[1] = 1.0E+00 + 2.0E+00 * z / ct3;
  aa[1] = 1.0E+00 + ct2 / ct3;

  anbn = ct3 + sab + 3.0E+00;
  ct1 = 1.0E+00 + 2.0E+00 * z / anbn;
  bb[2] = 1.0E+00 + 6.0E+00 * ct1 * z / ct3;
  aa[2] = 1.0E+00 + 6.0E+00 * ab / anbn + 3.0E+00 * ct1 * ct2 / ct3;

  for ( i = 4; i <= 300; i++ )
  {
    x2i1 = ( float ) ( 2 * i - 3 );
    ct1 = x2i1 / ( x2i1 - 2.0 );
    anbn = anbn + x2i1 + sab;
    ct2 = ( x2i1 - 1.0E+00 ) / anbn;
    c2 = x2i1 * ct2 - 1.0E+00;
    d1z = x2i1 * 2.0E+00 * z / anbn;

    ct3 = sab * ct2;
    g1 = d1z + ct1 * ( c2 + ct3 );
    g2 = d1z - c2;
    g3 = ct1 * ( 1.0E+00 - ct3 - 2.0E+00 * ct2 );

    bb[3] = g1 * bb[2] + g2 * bb[1] + g3 * bb[0];
    aa[3] = g1 * aa[2] + g2 * aa[1] + g3 * aa[0];

    value = aa[3] / bb[3];

    if ( r4_abs ( value - aa[0] / bb[0] ) < eps * r4_abs ( value ) )
    {
      return value;
    }
//
//  If overflows or underflows prove to be a problem, the statements
//  below could be altered to incorporate a dynamically adjusted scale
//  factor.
//
    for ( j = 0; j < 3; j++ )
    {
      bb[j] = bb[j+1];
      aa[j] = aa[j+1];
    }
  }

  cerr << "\n";
  cerr << "R4_CHU_SCALED - Fatal error!\n";
  cerr << "  No convergence in 300 terms.\n";

  exit ( 1 );
}
//****************************************************************************80

float r4_ci ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CI evaluates the cosine integral Ci of an R4 argument.
//
//  Discussion:
//
//    The cosine integral is defined by
//
//      CI(X) = - integral ( X <= T < oo ) ( cos ( T ) ) / T  dT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_CI, the cosine integral Ci evaluated at X.
//
{
  static float cics[13] = {
    -0.34004281856055363156E+00,
    -1.03302166401177456807E+00,
     0.19388222659917082877E+00,
    -0.01918260436019865894E+00,
     0.00110789252584784967E+00,
    -0.00004157234558247209E+00,
     0.00000109278524300229E+00,
    -0.00000002123285954183E+00,
     0.00000000031733482164E+00,
    -0.00000000000376141548E+00,
     0.00000000000003622653E+00,
    -0.00000000000000028912E+00,
     0.00000000000000000194E+00 };
  float f;
  float g;
  static int nci = 0;
  float sinx;
  float value;
  static float xsml = 0.0;
  float y;

  if ( nci == 0 )
  {
    nci = r4_inits ( cics, 13, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( r4_mach ( 3 ) );
  }

  if ( x <= 0.0E+00 )
  { 
    cerr << "\n";
    cerr << "R4_CI - Fatal error!\n";
    cerr << "  X <= 0.0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = - 1.0E+00;
    value = log ( x ) - 0.5E+00 + r4_csevl ( y, cics, nci );
  }
  else if ( x <= 4.0E+00 )
  {
    y = ( x * x - 8.0E+00 ) * 0.125E+00;
    value = log ( x ) - 0.5E+00 + r4_csevl ( y, cics, nci );
  }
  else
  {
    r4_sifg ( x, f, g );
    sinx = sin ( x );
    value = f * sinx - g * cos ( x );
  }

  return value;
}
//****************************************************************************80

float r4_cin ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CIN evaluates the alternate cosine integral Cin of an R4 argument.
//
//  Discussion:
//
//    CIN(X) = gamma + log(X) 
//      + integral ( 0 <= T <= X ) ( cos ( T ) - 1 ) / T  dT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_CIN, the cosine integral Cin evaluated at X.
//
{
  float absx;
  static float cincs[11] = {
     0.3707450175090968874E+00,
    -0.0589357489636444683E+00,
     0.0053818964211356912E+00,
    -0.0002986005284196214E+00,
     0.0000109557257532162E+00,
    -0.0000002840545487735E+00,
     0.0000000054697399488E+00,
    -0.0000000000812418746E+00,
     0.0000000000009586859E+00,
    -0.0000000000000092027E+00,
     0.0000000000000000733E+00 };
  static float eul = 0.57721566490153286E+00;
  float f;
  float g;
  static int ncin = 0;
  float sinx;
  float value;
  static float xmin = 0.0;

  if ( ncin == 0 )
  {
    ncin = r4_inits ( cincs, 11, 0.1E+00 * r4_mach ( 3 ) );
    xmin = sqrt ( r4_mach ( 1 ) );
  }

  absx = r4_abs ( x );

  if ( absx <= xmin )
  {
    value = 0.0E+00;
  }
  else if ( x <= 4.0E+00 )
  {
    value = x * x * r4_csevl ( ( x * x - 8.0E+00 ) * 0.125E+00, cincs, ncin );
  }
  else
  {
    r4_sifg ( x, f, g );
    sinx = sin ( absx );
    value = - f * sinx + g * cos ( absx ) + log ( absx ) + eul;
  }
  return value;
}
//****************************************************************************80

float r4_cinh ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CINH evaluates the alternate hyperbolic cosine integral Cinh of an R4 argument.
//
//  Discussion:
//
//    Cinh ( x ) = Integral ( 0 <= t <= x ) ( cosh ( t ) - 1 ) dt / t
//
//    The original text of this program had a mistake:
//      y = x * x / 9.0E+00 - 1.0E+00
//    has been corrected to
//      y = x * x / 4.5E+00 - 1.0E+00
//    JVB, 27 March 2010
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_CINH, the hyperbolic cosine integral Cinh
//    evaluated at X.
//
{
  float absx;
  static float cinhcs[10] = {
     0.1093291636520734431E+00,
     0.0573928847550379676E+00,
     0.0028095756978830353E+00,
     0.0000828780840721357E+00,
     0.0000016278596173914E+00,
     0.0000000227809519256E+00,
     0.0000000002384484842E+00,
     0.0000000000019360830E+00,
     0.0000000000000125454E+00,
     0.0000000000000000664E+00 };
  static float eul = 0.57721566490153286E+00;
  static int ncinh = 0;
  float value;
  static float xmin = 0.0;
  static float xsml = 0.0;
  float y;

  if ( ncinh == 0 )
  {
    ncinh = r4_inits ( cinhcs, 10, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( r4_mach ( 3 ) );
    xmin = 2.0E+00 * sqrt ( r4_mach ( 1 ) );
  }

  absx = r4_abs ( x );

  if ( x == 0.0E+00 )
  {
    value = 0.0E+00;
  }
  else if ( absx <= xmin )
  {
    value = 0.0E+00;
  }
  else if ( x <= xsml )
  {
    y = - 1.0E+00;
    value = x * x * ( 0.25E+00 + r4_csevl ( y, cinhcs, ncinh ) );
  }
  else if ( x <= 3.0E+00 )
  {
    y = x * x / 4.5E+00 - 1.0E+00;
    value = x * x * ( 0.25E+00 + r4_csevl ( y, cinhcs, ncinh ) );
  }
  else
  {
    value = r4_chi ( absx ) - eul - log ( absx );
  }
  return value;
}
//****************************************************************************80

float r4_cos ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_COS evaluates the cosine of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_COS, the cosine of X.
//
{
  float absx;
  float f;
  int n2;
  static int ntsn = 0;
  static float pi2 = 1.570796326794896619E+00;
  static float pi2rec = 0.636619772367581343E+00;
  static float pihi = 3.140625E+00;
  static float pilo = 9.6765358979323846E-04;
  static float pirec = 0.31830988618379067E+00;
  static float sincs[10] = {
     -0.374991154955873175840E+00,
     -0.181603155237250201864E+00,
     +0.005804709274598633559E+00,
     -0.000086954311779340757E+00,
     +0.000000754370148088851E+00,
     -0.000000004267129665056E+00,
     +0.000000000016980422945E+00,
     -0.000000000000050120579E+00,
     +0.000000000000000114101E+00,
     -0.000000000000000000206E+00 };
  float value;
  static float xmax = 0.0;
  float xn;
  static float xsml = 0.0;
  static float xwarn = 0.0;
  float y;

  if ( ntsn == 0 )
  {
    ntsn = r4_inits ( sincs, 10, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( 2.0E+00 * r4_mach ( 3 ) );
    xmax = 1.0E+00 / r4_mach ( 4 );
    xwarn = sqrt ( xmax );
  }
  
  absx = r4_abs ( x );
  y = absx + pi2;

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R4_COS - Warning!\n";
    cerr << "  No precision because |X| is big.\n";
    value = 0.0E+00;
    return value;
  }

  if ( xwarn < y )
  {
    cerr << "\n";
    cerr << "R4_COS - Warning!\n";
    cerr << "  Answer < half precision because |X| is big.\n";
  }

  value = 1.0E+00;

  if ( absx < xsml )
  {
    return value;
  }

  xn = r4_aint ( y * pirec + 0.5E+00 );
  n2 = ( int ) ( r4_mod ( xn, 2.0E+00 ) + 0.5E+00 );
  xn = xn - 0.5E+00;
  f = ( absx - xn * pihi ) - xn * pilo;

  xn = 2.0E+00 * ( f * pi2rec ) * ( f * pi2rec ) - 1.0E+00;
  value = f + f * r4_csevl ( xn, sincs, ntsn );

  if ( n2 != 0 )
  {
    value = - value;
  }

  if ( value < - 1.0E+00 )
  {
    value = - 1.0E+00;
  }
  else if ( 1.0E+00 < value )
  {
    value = + 1.0E+00;
  }
  return value;
}
//****************************************************************************80

float r4_cos_deg ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_COS_DEG evaluates the cosine of an R4 argument in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument in degrees.
//
//    Output, float R4_COS_DEG, the cosine of X.
//
{
  int n;
  static float raddeg = 0.017453292519943296E+00;
  float value;

  value = cos ( raddeg * x );

  if ( r4_mod ( x, 90.0E+00 ) == 0.0E+00 )
  {
    n = ( int ) ( r4_abs ( x ) / 90.0E+00 + 0.5E+00 );
    n = ( n % 2 );

    if ( n == 1 )
    {
      value = 0.0E+00;
    }
    else if ( value < 0.0E+00 )
    {
      value = - 1.0E+00;;
    }
    else
    {
      value = + 1.0E+00;
    }
  }
  return value;
}
//****************************************************************************80

float r4_cosh ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_COSH evaluates the hyperbolic cosine of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    1 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_COSH, the hyperbolic cosine of X.
//
{
  float value;
  float y;
  static float ymax = 0.0;

  if ( ymax == 0.0E+00 )
  {
    ymax = 1.0E+00 / sqrt ( r4_mach ( 3 ) );
  }

  y = exp ( r4_abs ( x ) );

  value = 0.5E+00 * y;

  if ( y < ymax )
  {
    value = 0.5E+00 * ( y + 1.0E+00 / y );
  }

  return value;
}
//****************************************************************************80

float r4_cot ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_COT evaluates the cotangent of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    1 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_COT, the cotangent of X.
//
{
  float ainty;
  float ainty2;
  static float cotcs[8] = {
    0.24025916098295630E+00,
   -0.016533031601500228E+00,
   -0.000042998391931724E+00,
   -0.000000159283223327E+00,
   -0.000000000619109313E+00,
   -0.000000000002430197E+00,
   -0.000000000000009560E+00,
   -0.000000000000000037E+00 };
  int ifn;
  static int nterms = 0;
  static float pi2rec = 0.0116197723675813430E+00;
  float prodbg;
  static float sqeps = 0.0;
  float value;
  static float xmax = 0.0;
  static float xmin = 0.0;
  static float xsml = 0.0;
  float y;
  float yrem;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( cotcs, 8, 0.1E+00 * r4_mach ( 3 ) );
    xmax = 1.0E+00 / r4_mach ( 4 );
    xsml = sqrt ( 3.0E+00 * r4_mach ( 3 ) );
    xmin = exp ( r4_max ( log ( r4_mach ( 1 ) ), 
      - log ( r4_mach ( 2 ) ) ) + 0.01E+00 );
    sqeps = sqrt ( r4_mach ( 4 ) );
  }

  y = r4_abs ( x );

  if ( y < xmin )
  {
    cerr << "\n";
    cerr << "R4_COT - Fatal error!\n";
    cerr << "  |X| is too small.\n";
    exit ( 1 );
  }

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R4_COT - Fatal error!\n";
    cerr << "  |X| is too big.\n";
    exit ( 1 );
  }
//
//  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
//  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
//  = aint(.625*y) + aint(z) + rem(z)
//
  ainty = r4_aint ( y );
  yrem = y - ainty;
  prodbg = 0.625E+00 * ainty;
  ainty = r4_aint ( prodbg );
  y = ( prodbg - ainty ) + 0.625E+00 * yrem + y * pi2rec;
  ainty2 = r4_aint ( y );
  ainty = ainty + ainty2;
  y = y - ainty2;

  ifn = ( int ) r4_mod ( ainty, 2.0E+00 );

  if ( ifn == 1 )
  {
    y = 1.0E+00 - y;
  }

  if ( 0.5E+00 < r4_abs ( x ) && y < r4_abs ( x ) * sqeps )
  {
    cerr << "\n";
    cerr << "R4_COT - Warning!\n";
    cerr << "  Answer less than half precision.\n";
    cerr << "  |X| too big, or X nearly a nonzero multiple of pi.\n";
  }

  if ( y == 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_COT - Fatal error!\n";
    cerr << "  X is a multiple of pi.\n";
    exit ( 1 );
  }
  else if ( y <= xsml )
  {
    value = 1.0E+00 / y;
  }
  else if ( y <= 0.25E+00 )
  {
    value = ( 0.5E+00 
      + r4_csevl ( 32.0E+00 * y * y - 1.0E+00, cotcs, nterms ) ) / y;
  }
  else if ( y <= 0.5E+00 )
  {
    value = ( 0.5E+00 + r4_csevl ( 8.0E+00 * y * y - 1.0E+00, 
      cotcs, nterms ) ) / ( 0.5E+00 * y );

    value = ( value * value - 1.0E+00 ) * 0.5E+00 / value;
  }
  else
  {
    value = ( 0.5E+00 + r4_csevl ( 2.0E+00 * y * y - 1.0E+00, 
      cotcs, nterms ) ) / ( 0.25E+00 * y );
    value = ( value * value - 1.0E+00 ) * 0.5E+00 / value;
    value = ( value * value - 1.0E+00 ) * 0.5E+00 / value;
  }

  if ( x < 0.0E+00 )
  {
    value = - r4_abs ( value );
  }
  else
  {
    value = + r4_abs ( value );
  }

  if ( ifn == 1 )
  {
    value = - value;
  }

  return value;
}
//****************************************************************************80

 float r4_csevl ( float x, float a[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CSEVL evaluates a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, float X, the evaluation point.
//
//    Input, float CS[N], the Chebyshev coefficients.
//
//    Input, int N, the number of Chebyshev coefficients.
//
//    Output, float R4_CSEVL, the Chebyshev series evaluated at X.
//
{
  float b0;
  float b1;
  float b2;
  int i;
  float twox;
  float value;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R4_CSEVL - Fatal error!\n";
    cerr << "  Number of terms <= 0.\n";
    exit ( 1 );
  }

  if ( 1000 < n )
  {
    cerr << "\n";
    cerr << "R4_CSEVL - Fatal error!\n";
    cerr << "  Number of terms greater than 1000.\n";
    exit ( 1 );
 }

  if ( x < -1.1 || 1.1 < x )
  {
    cerr << "\n";
    cerr << "R4_CSEVL - Fatal error!\n";
    cerr << "  X outside (-1,+1).\n";
    exit ( 1 );
  }

  twox = 2.0 * x;
  b1 = 0.0;
  b0 = 0.0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[i];
  }

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

float r4_dawson ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_DAWSON evaluates Dawson's integral of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_DAWSON, the value of Dawson's integral at X.
//
{
  static float daw2cs[29] = {
   -0.056886544105215527E+00,
   -0.31811346996168131E+00,
    0.20873845413642237E+00,
   -0.12475409913779131E+00,
    0.067869305186676777E+00,
   -0.033659144895270940E+00,
    0.015260781271987972E+00,
   -0.006348370962596214E+00,
    0.002432674092074852E+00,
   -0.000862195414910650E+00,
    0.000283765733363216E+00,
   -0.000087057549874170E+00,
    0.000024986849985481E+00,
   -0.000006731928676416E+00,
    0.000001707857878557E+00,
   -0.000000409175512264E+00,
    0.000000092828292216E+00,
   -0.000000019991403610E+00,
    0.000000004096349064E+00,
   -0.000000000800324095E+00,
    0.000000000149385031E+00,
   -0.000000000026687999E+00,
    0.000000000004571221E+00,
   -0.000000000000751873E+00,
    0.000000000000118931E+00,
   -0.000000000000018116E+00,
    0.000000000000002661E+00,
   -0.000000000000000377E+00,
    0.000000000000000051E+00 };
  static float dawacs[26] = {
    0.01690485637765704E+00,
    0.00868325227840695E+00,
    0.00024248640424177E+00,
    0.00001261182399572E+00,
    0.00000106645331463E+00,
    0.00000013581597947E+00,
    0.00000002171042356E+00,
    0.00000000286701050E+00,
   -0.00000000019013363E+00,
   -0.00000000030977804E+00,
   -0.00000000010294148E+00,
   -0.00000000000626035E+00,
    0.00000000000856313E+00,
    0.00000000000303304E+00,
   -0.00000000000025236E+00,
   -0.00000000000042106E+00,
   -0.00000000000004431E+00,
    0.00000000000004911E+00,
    0.00000000000001235E+00,
   -0.00000000000000578E+00,
   -0.00000000000000228E+00,
    0.00000000000000076E+00,
    0.00000000000000038E+00,
   -0.00000000000000011E+00,
   -0.00000000000000006E+00,
    0.00000000000000002E+00 };
  static float dawcs[13] = {
   -0.006351734375145949E+00,
   -0.22940714796773869E+00,
    0.022130500939084764E+00,
   -0.001549265453892985E+00,
    0.000084973277156849E+00,
   -0.000003828266270972E+00,
    0.000000146285480625E+00,
   -0.000000004851982381E+00,
    0.000000000142146357E+00,
   -0.000000000003728836E+00,
    0.000000000000088549E+00,
   -0.000000000000001920E+00,
    0.000000000000000038E+00 };
  float eps;
  static int ntdaw = 0;
  static int ntdaw2 = 0;
  static int ntdawa = 0;
  float value;
  static float xbig = 0.0;
  static float xmax = 0.0;
  static float xsml = 0.0;
  float y;

  if ( ntdaw == 0 )
  {
    eps = r4_mach ( 3 );
    ntdaw  = r4_inits ( dawcs,  13, 0.1E+00 * eps );
    ntdaw2 = r4_inits ( daw2cs, 29, 0.1E+00 * eps );
    ntdawa = r4_inits ( dawacs, 26, 0.1E+00 * eps );

    xsml = sqrt ( 1.5E+00 * eps );
    xbig = sqrt ( 0.5E+00 / eps );
    xmax = exp ( r4_min ( - log ( 2.0E+00 * r4_mach ( 1 ) ), 
      log ( r4_mach ( 2 ) ) ) - 0.01E+00 );
  }

  y = r4_abs ( x );

  if ( y <= xsml )
  {
    value = x;
  }
  else if ( y <= 1.0E+00 )
  {
    value = x * ( 0.75E+00 
      + r4_csevl ( 2.0E+00 * y * y - 1.0E+00, dawcs, ntdaw ) );
  }
  else if ( y <= 4.0E+00 )
  {
    value = x * ( 0.25E+00 
      + r4_csevl ( 0.125E+00 * y * y - 1.0E+00, daw2cs, ntdaw2 ) );
  }
  else if ( y < xbig )
  {
    value = ( 0.5E+00 
      + r4_csevl ( 32.0E+00 / y / y - 1.0E+00, dawacs, ntdawa ) ) / x;
  }
  else if ( y <= xmax )
  {
    value = 0.5E+00 / x;
  }
  else
  {
    value = 0.0E+00;
  }
  return value;
}
//****************************************************************************80

float r4_e1 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_E1 evaluates the exponential integral E1 for an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_E1, the exponential integral E1 evaluated at X.
//
{
  static float ae11cs[39] = {
    0.12150323971606579E+00,
   -0.065088778513550150E+00,
    0.004897651357459670E+00,
   -0.000649237843027216E+00,
    0.000093840434587471E+00,
    0.000000420236380882E+00,
   -0.000008113374735904E+00,
    0.000002804247688663E+00,
    0.000000056487164441E+00,
   -0.000000344809174450E+00,
    0.000000058209273578E+00,
    0.000000038711426349E+00,
   -0.000000012453235014E+00,
   -0.000000005118504888E+00,
    0.000000002148771527E+00,
    0.000000000868459898E+00,
   -0.000000000343650105E+00,
   -0.000000000179796603E+00,
    0.000000000047442060E+00,
    0.000000000040423282E+00,
   -0.000000000003543928E+00,
   -0.000000000008853444E+00,
   -0.000000000000960151E+00,
    0.000000000001692921E+00,
    0.000000000000607990E+00,
   -0.000000000000224338E+00,
   -0.000000000000200327E+00,
   -0.000000000000006246E+00,
    0.000000000000045571E+00,
    0.000000000000016383E+00,
   -0.000000000000005561E+00,
   -0.000000000000006074E+00,
   -0.000000000000000862E+00,
    0.000000000000001223E+00,
    0.000000000000000716E+00,
   -0.000000000000000024E+00,
   -0.000000000000000201E+00,
   -0.000000000000000082E+00,
    0.000000000000000017E+00 };
  static float ae12cs[25] = {
    0.58241749513472674E+00,
   -0.15834885090578275E+00,
   -0.006764275590323141E+00,
    0.005125843950185725E+00,
    0.000435232492169391E+00,
   -0.000143613366305483E+00,
   -0.000041801320556301E+00,
   -0.000002713395758640E+00,
    0.000001151381913647E+00,
    0.000000420650022012E+00,
    0.000000066581901391E+00,
    0.000000000662143777E+00,
   -0.000000002844104870E+00,
   -0.000000000940724197E+00,
   -0.000000000177476602E+00,
   -0.000000000015830222E+00,
    0.000000000002905732E+00,
    0.000000000001769356E+00,
    0.000000000000492735E+00,
    0.000000000000093709E+00,
    0.000000000000010707E+00,
   -0.000000000000000537E+00,
   -0.000000000000000716E+00,
   -0.000000000000000244E+00,
   -0.000000000000000058E+00 };
  static float ae13cs[25] = {
   -0.60577324664060346E+00,
   -0.11253524348366090E+00,
    0.013432266247902779E+00,
   -0.001926845187381145E+00,
    0.000309118337720603E+00,
   -0.000053564132129618E+00,
    0.000009827812880247E+00,
   -0.000001885368984916E+00,
    0.000000374943193568E+00,
   -0.000000076823455870E+00,
    0.000000016143270567E+00,
   -0.000000003466802211E+00,
    0.000000000758754209E+00,
   -0.000000000168864333E+00,
    0.000000000038145706E+00,
   -0.000000000008733026E+00,
    0.000000000002023672E+00,
   -0.000000000000474132E+00,
    0.000000000000112211E+00,
   -0.000000000000026804E+00,
    0.000000000000006457E+00,
   -0.000000000000001568E+00,
    0.000000000000000383E+00,
   -0.000000000000000094E+00,
    0.000000000000000023E+00 };
  static float ae14cs[26] = {
   -0.1892918000753017E+00,
   -0.08648117855259871E+00,
    0.00722410154374659E+00,
   -0.00080975594575573E+00,
    0.00010999134432661E+00,
   -0.00001717332998937E+00,
    0.00000298562751447E+00,
   -0.00000056596491457E+00,
    0.00000011526808397E+00,
   -0.00000002495030440E+00,
    0.00000000569232420E+00,
   -0.00000000135995766E+00,
    0.00000000033846628E+00,
   -0.00000000008737853E+00,
    0.00000000002331588E+00,
   -0.00000000000641148E+00,
    0.00000000000181224E+00,
   -0.00000000000052538E+00,
    0.00000000000015592E+00,
   -0.00000000000004729E+00,
    0.00000000000001463E+00,
   -0.00000000000000461E+00,
    0.00000000000000148E+00,
   -0.00000000000000048E+00,
    0.00000000000000016E+00,
   -0.00000000000000005E+00 };
  static float e11cs[19] = {
  -16.113461655571494026E+00,
    7.7940727787426802769E+00,
   -1.9554058188631419507E+00,
    0.37337293866277945612E+00,
   -0.05692503191092901938E+00,
    0.00721107776966009185E+00,
   -0.00078104901449841593E+00,
    0.00007388093356262168E+00,
   -0.00000620286187580820E+00,
    0.00000046816002303176E+00,
   -0.00000003209288853329E+00,
    0.00000000201519974874E+00,
   -0.00000000011673686816E+00,
    0.00000000000627627066E+00,
   -0.00000000000031481541E+00,
    0.00000000000001479904E+00,
   -0.00000000000000065457E+00,
    0.00000000000000002733E+00,
   -0.00000000000000000108E+00 };
  static float e12cs[16] = {
    -0.037390214792202795E+00,
    0.042723986062209577E+00,
   -0.1303182079849700544E+00,
    0.01441912402469889073E+00,
   -0.00134617078051068022E+00,
    0.00010731029253063780E+00,
   -0.00000742999951611943E+00,
    0.00000045377325690753E+00,
   -0.00000002476417211390E+00,
    0.00000000122076581374E+00,
   -0.00000000005485141480E+00,
    0.00000000000226362142E+00,
   -0.00000000000008635897E+00,
    0.00000000000000306291E+00,
   -0.00000000000000010148E+00,
    0.00000000000000000315E+00 };
  float eta;
  static int ntae11 = 0;
  static int ntae12 = 0;
  static int ntae13 = 0;
  static int ntae14 = 0;
  static int nte11 = 0;
  static int nte12 = 0;
  float value;
  static float xmax = 0.0;

  if ( ntae11 == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    ntae11 = r4_inits ( ae11cs, 39, eta );
    ntae12 = r4_inits ( ae12cs, 25, eta );
    nte11 = r4_inits ( e11cs, 19, eta );
    nte12 = r4_inits ( e12cs, 16, eta );
    ntae13 = r4_inits ( ae13cs, 25, eta );
    ntae14 = r4_inits ( ae14cs, 26, eta );
    xmax = - log ( r4_mach ( 1 ) );
    xmax = xmax - log ( xmax );
  }

  if ( x <= - 10.0E+00 )
  {
    value = exp ( - x ) / x * ( 1.0E+00 
      + r4_csevl ( 20.0E+00 / x + 1.0E+00, ae11cs, ntae11 ) );
  }
  else if ( x <= - 4.0E+00 )
  {
    value = exp ( - x ) / x * ( 1.0E+00 + r4_csevl ( 
      ( 40.0E+00 / x + 7.0E+00 ) / 3.0E+00, ae12cs, ntae12 ) );
  }
  else if ( x <= - 1.0E+00 )
  {
    value = - log ( r4_abs ( x ) ) + r4_csevl ( 
      ( 2.0E+00 * x + 5.0E+00 ) / 3.0E+00, e11cs, nte11 );
  }
  else if ( x == 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_E1 - Fatal error!\n";
    cerr << "  X is zero.\n";
    exit ( 1 );
  }
  else if ( x <= 1.0E+00 )
  {
    value = ( - log ( r4_abs ( x ) ) - 0.6875E+00 + x ) 
      + r4_csevl ( x, e12cs, nte12 );
  }
  else if ( x <= 4.0E+00 )
  {
    value = exp ( - x ) / x * ( 1.0E+00 + r4_csevl ( 
      ( 8.0E+00 / x - 5.0E+00 ) / 3.0E+00, ae13cs, ntae13 ) );
  }
  else if ( x <= xmax )
  {
    value = exp ( - x ) / x * ( 1.0E+00 + r4_csevl (
      8.0E+00 / x - 1.0E+00, ae14cs, ntae14 ) );
  }
  else
  {
    value = 0.0E+00;
  }
  return value;
}
//****************************************************************************80

float r4_ei ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_EI evaluates the exponential integral Ei for an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_EI, the exponential integral Ei evaluated at X.
//
{
  float value;

  value = - r4_e1 ( - x );

  return value;
}
//****************************************************************************80

float r4_erf ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ERF evaluates the error function of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_ERF, the error function of X.
//
{
  static float erfcs[13] = {
   -0.049046121234691808E+00,
   -0.14226120510371364E+00,
    0.010035582187599796E+00,
   -0.000576876469976748E+00,
    0.000027419931252196E+00,
   -0.000001104317550734E+00,
    0.000000038488755420E+00,
   -0.000000001180858253E+00,
    0.000000000032334215E+00,
   -0.000000000000799101E+00,
    0.000000000000017990E+00,
   -0.000000000000000371E+00,
    0.000000000000000007E+00 };
  static float nterf = 0;
  static float sqeps = 0.0;
  static float sqrtpi = 1.7724538509055160E+00;
  float value;
  static float xbig = 0.0;
  float y;

  if ( nterf == 0 )
  {
    nterf = r4_inits ( erfcs, 13, 0.1E+00 * r4_mach ( 3 ) );
    xbig = sqrt ( - log ( sqrtpi * r4_mach ( 3 ) ) );
    sqeps = sqrt ( 2.0E+00 * r4_mach ( 3 ) );
  }

  y = r4_abs ( x );

  if ( y <= sqeps )
  {
    value = 2.0E+00 * x / sqrtpi;
  }
  else if ( y <= 1.0E+00 )
  {
    value = x * ( 1.0E+00 
      + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, erfcs, nterf ) );
  }
  else if ( y <= xbig )
  {
    value = 1.0E+00 - r4_erfc ( y );
    if ( x < 0.0E+00 )
    {
      value = - value;
    }
  }
  else
  {
    value = 1.0E+00;
    if ( x < 0.0E+00 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

float r4_erfc ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ERFC evaluates the co-error function of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_ERFC, the co-error function of X.
//
{
  static float erc2cs[23] = {
   -0.069601346602309501E+00,
   -0.041101339362620893E+00,
    0.003914495866689626E+00,
   -0.000490639565054897E+00,
    0.000071574790013770E+00,
   -0.000011530716341312E+00,
    0.000001994670590201E+00,
   -0.000000364266647159E+00,
    0.000000069443726100E+00,
   -0.000000013712209021E+00,
    0.000000002788389661E+00,
   -0.000000000581416472E+00,
    0.000000000123892049E+00,
   -0.000000000026906391E+00,
    0.000000000005942614E+00,
   -0.000000000001332386E+00,
    0.000000000000302804E+00,
   -0.000000000000069666E+00,
    0.000000000000016208E+00,
   -0.000000000000003809E+00,
    0.000000000000000904E+00,
   -0.000000000000000216E+00,
    0.000000000000000052E+00 };
  static float erfccs[24] = {
    0.0715179310202925E+00,
   -0.026532434337606719E+00,
    0.001711153977920853E+00,
   -0.000163751663458512E+00,
    0.000019871293500549E+00,
   -0.000002843712412769E+00,
    0.000000460616130901E+00,
   -0.000000082277530261E+00,
    0.000000015921418724E+00,
   -0.000000003295071356E+00,
    0.000000000722343973E+00,
   -0.000000000166485584E+00,
    0.000000000040103931E+00,
   -0.000000000010048164E+00,
    0.000000000002608272E+00,
   -0.000000000000699105E+00,
    0.000000000000192946E+00,
   -0.000000000000054704E+00,
    0.000000000000015901E+00,
   -0.000000000000004729E+00,
    0.000000000000001432E+00,
   -0.000000000000000439E+00,
    0.000000000000000138E+00,
   -0.000000000000000048E+00 };
  static float erfcs[13] = {
   -0.049046121234691808E+00,
   -0.14226120510371364E+00,
    0.010035582187599796E+00,
   -0.000576876469976748E+00,
    0.000027419931252196E+00,
   -0.000001104317550734E+00,
    0.000000038488755420E+00,
   -0.000000001180858253E+00,
    0.000000000032334215E+00,
   -0.000000000000799101E+00,
    0.000000000000017990E+00,
   -0.000000000000000371E+00,
    0.000000000000000007E+00 };
  float eta;
  static int nterc2 = 0;
  static int nterf = 0;
  static int nterfc = 0;
  static float sqeps;
  static float sqrtpi = 1.7724538509055160E+00;
  float value;
  static float xmax = 0.0;
  static float xsml = 0.0;
  float y;

  if ( nterf == 0 )
  {
    eta = 0.1E+00 * r4_mach ( 3 );
    nterf = r4_inits ( erfcs, 13, eta );
    nterfc = r4_inits ( erfccs, 24, eta );
    nterc2 = r4_inits ( erc2cs, 23, eta );

    xsml = - sqrt ( - log ( sqrtpi * r4_mach ( 3 ) ) );
    xmax = sqrt ( - log ( sqrtpi * r4_mach ( 1 ) ) );
    xmax = xmax - 0.5E+00 * log ( xmax ) / xmax - 0.01E+00;
    sqeps = sqrt ( 2.0E+00 * r4_mach ( 3 ) );
  }

  if ( x <= xsml )
  {
    value = 2.0E+00;
    return value;
  }

  if ( xmax < x )
  {
    cerr << "\n";
    cerr << "R4_ERFC - Warning!\n";
    cerr << "  X so big that ERFC underflows.\n";
    value = 0.0E+00;
    return value;
  }

  y = r4_abs ( x );

  if ( y < sqeps )
  {
    value = 1.0E+00 - 2.0E+00 * x / sqrtpi;
    return value;
  }
  else if ( y <= 1.0E+00 )
  {
    value = 1.0E+00 - x * ( 1.0E+00 
      + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, erfcs, nterf ) );
    return value;
  }

  y = y * y;

  if ( y <= 4.0E+00 )
  {
    value = exp ( - y ) / r4_abs ( x ) * ( 0.5E+00 
      + r4_csevl ( ( 8.0E+00 / y - 5.0E+00 ) / 3.0E+00, erc2cs, 
      nterc2 ) );
  }
  else
  {
    value = exp ( - y ) / r4_abs ( x ) * ( 0.5E+00 
      + r4_csevl ( 8.0E+00 / y - 1.0E+00, erfccs, nterfc ) );
  }

  if ( x < 0.0E+00 )
  {
    value = 2.0E+00 - value;
  }

  return value;
}
//****************************************************************************80

float r4_exp ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_EXP evaluates the exponential of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_EXP, the exponential of X.
//
{
  static float aln216 = 0.083120654223414518;
  static float expcs[8] = {
      0.086656949331498571,
      0.000938494869299839,
      0.000006776039709981,
      0.000000036693120039,
      0.000000000158959053,
      0.000000000000573859,
      0.000000000000001775,
      0.000000000000000004 };
  float f;
  int n;
  int n16;
  int ndx;
  static int nterms = 0;
  static float twon16[17] = {
     0.0E+00,
     0.44273782427413840E-01,
     0.90507732665257659E-01,
     0.13878863475669165E+00,
     0.18920711500272107E+00,
     0.24185781207348405E+00,
     0.29683955465100967E+00,
     0.35425554693689273E+00,
     0.41421356237309505E+00,
     0.47682614593949931E+00,
     0.54221082540794082E+00,
     0.61049033194925431E+00,
     0.68179283050742909E+00,
     0.75625216037329948E+00,
     0.83400808640934246E+00,
     0.91520656139714729E+00,
     1.0E+00 };
  float value;
  float xint;
  static float xmax = 0.0;
  static float xmin = 0.0;
  float y;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( expcs, 8, 0.1E+00 * r4_mach ( 3 ) );
    xmin = log ( r4_mach ( 1 ) ) + 0.01E+00;
    xmax = log ( r4_mach ( 2 ) ) - 0.001E+00;
  }

  if ( x < xmin )
  {
    cerr << "\n";
    cerr << "R4_EXP - Warning!\n";
    cerr << "  X so small that exp(X) underflows.\n";
    value = 0.0E+00;
  }
  else if ( x <= xmax )
  {
    xint = r4_aint ( x );
    y = x - xint;

    y = 23.0E+00 * y + x * aln216;
    n = ( int ) ( y );
    f = y - ( float ) ( n );
    n = ( int ) ( 23.0E+00 * xint + ( float ) ( n ) );
    n16 = n / 16;
    if ( n < 0 )
    {
      n16 = n16 - 1;
    }
    ndx = n - 16 * n16 + 1;

    value = 1.0E+00 + ( twon16[ndx-1] + f * ( 1.0E+00 + twon16[ndx-1] ) 
      * r4_csevl ( f, expcs, nterms ) );

    value = r4_pak ( value, n16 );
  }
  else
  {
    cerr << "\n";
    cerr << "R4_EXP - Fatal error!\n";
    cerr << "  X so large that exp(X) overflows.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

float r4_exprel ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_EXPREL evaluates the exponential relative error term of an R4 argument.
//
//  Discussion:
//
//    The relative error term is ( exp ( x ) - 1 ) / x.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_EXPREL, the exponential relative error term
//    at X.
//
{
  float absx;
  float alneps;
  int i;
  static int nterms = 0;
  float value;
  static float xbnd = 0.0;
  float xln;
  float xn;

  if ( nterms == 0 )
  {
    alneps = log ( r4_mach ( 3 ) );
    xn = 3.72E+00 - 0.3E+00 * alneps;
    xln = log ( ( xn + 1.0E+00 ) / 1.36E+00 );
    nterms = ( int ) ( xn - ( xn * xln + alneps ) / ( xln + 1.36E+00 ) 
      + 1.5E+00 );
    xbnd = r4_mach ( 3 );
  }

  absx = r4_abs ( x );

  if ( absx < xbnd )
  {
    value = 1.0E+00;
  }
  else if ( absx <= 0.5E+00 )
  {
    value = 0.0E+00;
    for ( i = 1; i <= nterms; i++ )
    {
      value = 1.0E+00 + value * x / ( float ) ( nterms + 2 - i );
    }
  }
  else
  {
    value = ( exp ( x ) - 1.0E+00 ) / x;
  }
  return value;
}
//****************************************************************************80

float r4_fac ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4_FAC evaluates the factorial of an I4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, int N, the argument.
//
//    Output, float R4_FAC, the factorial of N.
//
{
  static float facn[26] = {
    1.0E+00,
    1.0E+00,
    2.0E+00,
    6.0E+00,
    24.0E+00,
    120.0E+00,
    720.0E+00,
    5040.0E+00,
    40320.0E+00,
    362880.0E+00,
    3628800.0E+00,
    39916800.0E+00,
    479001600.0E+00,
    6227020800.0E+00,
    87178291200.0E+00,
    1307674368000.0E+00,
    20922789888000.0E+00,
    355687428096000.0E+00,
    6402373705728000.0E+00,
    0.12164510040883200E+18,
    0.24329020081766400E+19,
    0.51090942171709440E+20,
    0.11240007277776077E+22,
    0.25852016738884977E+23,
    0.62044840173323944E+24,
    0.15511210043330986E+26 };
  static int nmax = 0;
  static float sq2pil = 0.91893853320467274E+00;
  float value;
  float x;
  float xmax;
  float xmin;

  if ( nmax == 0 )
  {
    r4_gaml ( xmin, xmax );
    nmax = ( int ) ( xmax - 1.0E+00 );
  }

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "R4_FAC - Fatal error!\n";
    cerr << "  Input argument is negative.\n";
    exit ( 1 );
  }
  else if ( n <= 25 )
  {
    value = facn[n];
  }
  else if ( n <= nmax )
  {
    x = ( float ) ( n + 1 );
    value = exp ( ( x - 0.5E+00 ) * log ( x ) - x + sq2pil 
      + r4_lgmc ( x ) );
  }
  else
  {
    cerr << "\n";
    cerr << "R4_FAC - Fatal error!\n";
    cerr << "  Factorial overflows.\n";
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

float r4_gami ( float a, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_GAMI evaluates the incomplete gamma function for an R4 argument.
//
//  Discussion:
//
//    GAMI = Integral ( 0 <= T <= X ) exp ( - t ) * t^( a - 1 )  dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, the parameter.
//
//    Input, float X, the argument.
//
//    Output, float R4_GAMI, the value of the incomplete gamma function.
//
{
  float factor;
  float value;

  if ( a <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_GAMI - Fatal error!\n";
    cerr << "  A <= 0.\n";
    exit ( 1 );
  }

  if ( x < 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_GAMI - Fatal error!\n";
    cerr << "  X < 0.\n";
    exit ( 1 );
  }
  else if ( x == 0.0E+00 )
  {
    value = 0.0E+00;
  }
  else
  {
    factor = exp ( r4_lngam ( a ) + a * log ( x ) );

    value = factor * r4_gamit ( a, x );
  }
  return value;
}
//****************************************************************************80

float r4_gamic ( float a, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_GAMIC evaluates the complementary incomplete gamma function.
//
//  Discussion:
//
//    GAMIC = integral ( x <= t < oo ) exp(-t) * t^(a-1) dt
//
//    GAMIC is evaluated for arbitrary float values of A and non-negative
//    values X (even though GAMIC is defined for X < 0.0), except that
//    for X = 0 and A <= 0.0, GAMIC is undefined.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//    Walter Gautschi,
//    A Computational Procedure for Incomplete Gamma Functions,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 4, December 1979, pages 466-481.
//
//  Parameters:
//
//    Input, float A, the parameter.
//
//    Input, float X, the evaluation point.
//
//    Output, float R4_GAMIC, the value of the incomplete gamma function.
//
{
  float aeps;
  float algap1;
  static float alneps = 0.0;
  float alngs;
  float alx;
  static float bot = 0.0;
  float e;
  static float eps = 0.0;
  float fm;
  float gstar;
  float h;
  int izero;
  int ma;
  float sga;
  float sgng;
  float sgngam;
  float sgngs;
  static float sqeps = 0.0;
  float t;
  float value;

  if ( eps == 0.0E+00 )
  {
    eps = 0.5E+00 * r4_mach ( 3 );
    sqeps = sqrt ( r4_mach ( 4 ) );
    alneps = - log ( r4_mach ( 3 ) );
    bot = log ( r4_mach ( 1 ) );
  }

  if ( x < 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_GAMIC - Fatal error!\n";
    cerr << "  X < 0.\n";
    exit ( 1 );
  }

  if ( x == 0.0E+00 )
  {
    if ( a <= 0.0E+00 )
    {
      cerr << "\n";
      cerr << "R4_GAMIC - Fatal error!\n";
      cerr << "  X = 0 and A <= 0.\n";
      exit ( 1 );
    }

    value = exp ( r4_lngam ( a + 1.0E+00 ) - log ( a ) );
    return value;
  }

  alx = log ( x );
  if ( a < 0.0E+00 )
  {
    sga = - 1.0E+00;
  }
  else
  {
    sga = + 1.0E+00;
  }

  ma = ( int ) ( a + 0.5E+00 * sga );
  aeps = a - ( float ) ( ma );

  izero = 0;

  if ( x < 1.0E+00 )
  {
    if ( a <= 0.5E+00 && r4_abs ( aeps ) <= 0.001E+00 )
    {
      fm = - ( float ) ( ma );

      if ( fm <= 1.0E+00 )
      {
        e = 2.0E+00;
      }
      else
      {
        e = 2.0E+00 * ( fm + 2.0E+00 ) / ( fm * fm - 1.0E+00 );
      }

      e = e - alx * r4_power ( x, - 0.001E+00 );

      if ( e * r4_abs ( aeps ) <= eps )
      {
        value = r4_gmic ( a, x, alx );
        return value;
      }
    }

    r4_lgams ( a + 1.0E+00, algap1, sgngam );
    gstar = r4_gmit ( a, x, algap1, sgngam, alx );

    if ( gstar == 0.0E+00 )
    {
      izero = 1;
    }
    else
    {
      alngs = log ( r4_abs ( gstar ) );
      sgngs = r4_sign ( gstar );
    }
  }
  else
  {
    if ( a < x )
    {
      value = exp ( r4_lgic ( a, x, alx ) );
      return value;
    }
    sgngam = 1.0E+00;
    algap1 = r4_lngam ( a + 1.0E+00 );
    sgngs = 1.0E+00;
    alngs = r4_lgit ( a, x, algap1 );
  }

  h = 1.0E+00;

  if ( izero != 1 )
  {
    t = a * alx + alngs;

    if ( alneps < t )
    {
      sgng = - sgngs * sga * sgngam;
      t = t + algap1 - log ( r4_abs ( a ) );
      value = sgng * exp ( t );
      return value;
    }

    if ( - alneps < t )
    {
      h = 1.0E+00 - sgngs * exp ( t );
    }
  }
  sgng = r4_sign ( h ) * sga * sgngam;
  t = log ( r4_abs ( h ) ) + algap1 - log ( r4_abs ( a ) );
  value = sgng * exp ( t );
  return value;
}
//****************************************************************************80

float r4_gamit ( float a, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_GAMIT evaluates Tricomi's incomplete gamma function for an R4 argument.
//
//  Discussion:
//
//      GAMIT = x^(-a) / gamma(a) 
//        * Integral ( 0 <= t <= x ) exp(-t) * t^(a-1) dt
//
//    with analytic continuation for a <= 0.0.  Gamma(x) is the complete
//    gamma function of X.  GAMIT is evaluated for arbitrary float values of
//    A and for non-negative values of X (even though GAMIT is defined for
//    X < 0.0).
//
//    A slight deterioration of 2 or 3 digits accuracy will occur when
//    gamit is very large or very small in absolute value, because log-
//    arithmic variables are used.  Also, if the parameter A is very close
//    to a negative integer (but not a negative integer), there is a loss
//    of accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//    Walter Gautschi,
//    A Computational Procedure for Incomplete Gamma Functions,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 4, December 1979, pages 466-481.
//
//  Parameters:
//
//    Input, float A, the parameter.
//
//    Input, float X, the argument.
//
//    Output, float R4_GAMIT, the function value.
//
{
  float aeps;
  float ainta;
  float algap1;
  static float alneps = 0.0;
  float alng;
  float alx;
  static float bot = 0.0;
  float h;
  float sga;
  float sgngam;
  static float sqeps = 0.0;
  float t;
  float value;

  if ( alneps == 0.0E+00 )
  {
    alneps = - log ( r4_mach ( 3 ) );
    sqeps = sqrt ( r4_mach ( 4 ) );
    bot = log ( r4_mach ( 1 ) );
  }

  if ( x < 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_GAMIT - Fatal error!\n";
    cerr << "  X is negative.\n";
    exit ( 1 );
  }
  else if ( x == 0.0E+00 )
  {
    alx = 0.0E+00;
  }
  else
  {
    alx = log ( x );
  }

  if ( a < 0.0E+00 )
  {
    sga = - 1.0E+00;
  }
  else
  {
    sga = + 1.0E+00;
  }

  ainta = r4_aint ( a + 0.5E+00 * sga );
  aeps = a - ainta;

  if ( x == 0.0E+00 )
  {
    if ( 0.0E+00 < ainta || aeps != 0.0E+00 )
    {
      value = r4_gamr ( a + 1.0E+00 );
    }
    else
    {
      value = 0.0E+00;
    }
    return value;
  }

  if ( x <= 1.0E+00 )
  {
    if ( - 0.5E+00 <= a || aeps != 0.0E+00 )
    {
      r4_lgams ( a + 1.0E+00, algap1, sgngam );
    }
    value = r4_gmit ( a, x, algap1, sgngam, alx );
    return value;
  }

  if ( x <= a )
  {
    t = r4_lgit ( a, x, r4_lngam ( a + 1.0E+00 ) );
    value = exp ( t );
    return value;
  }

  alng = r4_lgic ( a, x, alx );
//
//  Evaluate in terms of log(value(a,x))
//
  h = 1.0E+00;

  if ( aeps != 0.0E+00 || 0.0E+00 < ainta )
  {
    r4_lgams ( a + 1.0E+00, algap1, sgngam );
    t = log ( r4_abs ( a ) ) + alng - algap1;

    if ( alneps < t )
    {
      t = t - a * alx;
      value = - sga * sgngam * exp ( t );
      return value;
    }

    if ( - alneps < t )
    {
      h = 1.0E+00 - sga * sgngam * exp ( t );
    }

  }

  t = - a * alx + log ( r4_abs ( h ) );

  if ( h < 0.0E+00 )
  {
    value = - exp ( t );
  }
  else
  {
    value = + exp ( t );
  }
  return value;
}
//****************************************************************************80

void r4_gaml ( float &xmin, float &xmax )

//****************************************************************************80
//
//  Purpose:
//
//    R4_GAML evaluates bounds for an R4 argument of the gamma function.
//
//  Discussion:
//
//    This function calculates the minimum and maximum legal bounds 
//    for X in the evaluation of GAMMA ( X ).
//
//    XMIN and XMAX are not the only bounds, but they are the only 
//    non-trivial ones to calculate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Output, float &XMIN, &XMAX, the bounds.
//
{
  float alnbig;
  float alnsml;
  int i;
  int j;
  float xln;
  float xold;

  alnsml = log ( r4_mach ( 1 ) );
  xmin = - alnsml;

  for ( i = 1; i <= 10; i++ )
  {
    xold = xmin;
    xln = log ( xmin );
    xmin = xmin - xmin * ( ( xmin + 0.5E+00 ) * xln - xmin 
      - 0.2258E+00 + alnsml ) / ( xmin * xln + 0.5E+00 );

    if ( r4_abs ( xmin - xold ) < 0.005E+00 )
    {
      xmin = - xmin + 0.01E+00;

      alnbig = log ( r4_mach ( 2 ) );
      xmax = alnbig;

      for ( j = 1; j <= 10; j++ )
      {
        xold = xmax;
        xln = log ( xmax );
        xmax = xmax - xmax * ( ( xmax - 0.5E+00 ) * xln - xmax 
          + 0.9189E+00 - alnbig ) / ( xmax * xln - 0.5E+00 );

        if ( r4_abs ( xmax - xold ) < 0.005E+00 )
        {
          xmax = xmax - 0.01E+00;
          xmin = r4_max ( xmin, - xmax + 1.0E+00 );
          return;
        }
      }
      cerr << "\n";
      cerr << "R4_GAML - Fatal error!\n";
      cerr << "  Unable to find XMAX.\n";
      exit ( 1 );
    }
  }
  cerr << "\n";
  cerr << "R4_GAML - Fatal error!\n";
  cerr << "  Unable to find XMIN.\n";
  exit ( 1 );
}
//****************************************************************************80

float r4_gamma ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_GAMMA evaluates the gamma function of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_GAMMA, the gamma function of X.
//
{
  static float dxrel = 0.0;
  static float gcs[23] = {
     0.008571195590989331E+00,
     0.004415381324841007E+00,
     0.05685043681599363E+00,
    -0.004219835396418561E+00,
     0.001326808181212460E+00,
    -0.0001893024529798880E+00,
     0.0000360692532744124E+00,
    -0.0000060567619044608E+00,
     0.0000010558295463022E+00,
    -0.0000001811967365542E+00,
     0.0000000311772496471E+00,
    -0.0000000053542196390E+00,
     0.0000000009193275519E+00,
    -0.0000000001577941280E+00,
     0.0000000000270798062E+00,
    -0.0000000000046468186E+00,
     0.0000000000007973350E+00,
    -0.0000000000001368078E+00,
     0.0000000000000234731E+00,
    -0.0000000000000040274E+00,
     0.0000000000000006910E+00,
    -0.0000000000000001185E+00,
     0.0000000000000000203E+00 };
  int i;
  int n;
  static int ngcs = 0;
  static float pi = 3.14159265358979324E+00;
  float sinpiy;
  static float sq2pil = 0.91893853320467274E+00;
  float value;
  static float xmax = 0.0;
  static float xmin = 0.0;
  static float xsml = 0.0;
  float y;

  if ( ngcs == 0 )
  {
    ngcs = r4_inits ( gcs, 23, 0.1E+00 * r4_mach ( 3 ) );
    r4_gaml ( xmin, xmax );
    xsml = exp ( r4_max ( log ( r4_mach ( 1 ) ),
      - log ( r4_mach ( 2 ) ) ) + 0.01E+00 );
    dxrel = sqrt ( r4_mach ( 4 ) );
  }

  y = r4_abs ( x );

  if ( y <= 10.0E+00 )
  {
    n = ( int ) ( x );
    if ( x < 0.0E+00 )
    {
      n = n - 1;
    }
    y = x - ( float ) ( n );
    n = n - 1;
    value = 0.9375E+00 + r4_csevl ( 2.0E+00 * y - 1.0E+00, gcs, ngcs );

    if ( n == 0 )
    {
      return value;
    }
    else if ( n < 0 )
    {
      n = - n;

      if ( x == 0.0E+00 )
      {
        cerr << "\n";
        cerr << "R4_GAMMA - Fatal error!\n";
        cerr << "  X is 0.\n";
        exit ( 1 );
      }

      if ( x < 0.0E+00 && x + ( float ) ( n - 2 ) == 0.0E+00 )
      {
        cerr << "\n";
        cerr << "R4_GAMMA - Fatal error!\n";
        cerr << "  X is a negative integer.\n";
        exit ( 1 );
      }

      if ( x < - 0.5E+00 && r4_abs ( ( x - r4_aint ( x - 0.5E+00 ) ) / x ) < dxrel )
      {
        cerr << "\n";
        cerr << "R4_GAMMA - Warning!\n";
        cerr << "  X too near a negative integer,\n";
        cerr << "  answer is half precision.\n";
      }

      if ( y < xsml )
      {
        cerr << "\n";
        cerr << "R4_GAMMA - Fatal error!\n";
        cerr << "  X is so close to zero that Gamma overflows.\n";
        cerr << "  X = " << x << "\n";
        exit ( 1 );
      }

      for ( i = 1; i <= n; i++ )
      {
        value = value / ( x + ( float ) ( i - 1 ) );
      }
    }
    else if ( n == 0 )
    {
    }
    else
    {
      for ( i = 1; i <= n; i++ )
      {
        value = ( y + ( float ) ( i ) ) * value;
      }
    }
  }
  else
  {
    if ( xmax < x )
    {
      cerr << "\n";
      cerr << "R4_GAMMA - Fatal error!\n";
      cerr << "  X so big that Gamma overflows.\n";
      exit ( 1 );
    }
//
//  Underflow.
//
    if ( x < xmin )
    {
      value = 0.0E+00;
      return value;
    }

    value = exp ( ( y - 0.5E+00 ) * log ( y ) - y + sq2pil + r4_lgmc ( y ) );

    if ( 0.0E+00 < x )
    {
      return value;
    }

    if ( r4_abs ( ( x - r4_aint ( x - 0.5E+00 ) ) / x ) < dxrel )
    {
      cerr << "\n";
      cerr << "R4_GAMMA - Warning!\n";
      cerr << "  X too near a negative integer,\n";
      cerr << "  answer is half precision.\n";
    }

    sinpiy = sin ( pi * y );

    if ( sinpiy == 0.0E+00 )
    {
      cerr << "\n";
      cerr << "R4_GAMMA - Fatal error!\n";
      cerr << "  X is a negative integer.\n";
      exit ( 1 );
    }
    value = - pi / ( y * sinpiy * value );
  }
  return value;
}
//****************************************************************************80

float r4_gamr ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_GAMR evaluates the reciprocal gamma function of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_GAMR, the value of the reciprocal gamma
//    function at X.
//
{
  float alngx;
  float sgngx;
  float value;

  if ( x <= 0.0E+00 && r4_aint ( x ) == x )
  {
    value = 0.0E+00;
  }
  else if ( r4_abs ( x ) <= 10.0E+00 )
  {
    value = 1.0E+00 / r4_gamma ( x );
  }
  else
  {
    r4_lgams ( x, alngx, sgngx );
    value = sgngx * exp ( - alngx );
  }

  return value;
 }
//****************************************************************************80

float r4_gmic ( float a, float x, float alx )

//****************************************************************************80
//
//  Purpose:
//
//    R4_GMIC: complementary incomplete gamma, small X, A near negative integer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, the parameter.
//
//    Input, float X, the argument.
//
//    Input, float ALX, the logarithm of X.
//
//    Output, float R4_GMIC, the complementary incomplete 
//    gamma function.
//
{
  float alng;
  static float bot = 0.0;
  bool converged;
  static float eps = 0.0;
  static float euler = 0.5772156649015329E+00;
  float fk;
  float fkp1;
  float fm;
  int k;
  int m;
  int ma;
  int mm1;
  float s;
  float sgng;
  float t;
  float te;
  float value;

  if ( eps == 0.0E+00 )
  {
    eps = 0.5E+00 * r4_mach ( 3 );
    bot = log ( r4_mach ( 1 ) );
  }

  if ( 0.0E+00 < a )
  {
    cerr << "\n";
    cerr << "R4_GMIC - Fatal error!\n";
    cerr << "  A must be near a negative integer.\n";
    exit ( 1 );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_GMIC - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  ma = ( int ) ( a - 0.5E+00 );
  fm = - ( float ) ( ma );
  m = - ma;

  te = 1.0E+00;
  t = 1.0E+00;
  s = t;
  converged = false;

  for ( k = 1; k <= 200; k++ )
  {
    fkp1 = ( float ) ( k + 1 );
    te = - x * te / ( fm + fkp1 );
    t = te / fkp1;
    s = s + t;
    if ( r4_abs ( t ) < eps * s )
    {
      converged = true;
      break;
    }
  }

  if ( !converged )
  {
    cerr << "\n";
    cerr << "R4_GMIC - Fatal error!\n";
    cerr << "  No convergence after 200 iterations.\n";
    exit ( 1 );
  }

  value = - alx - euler + x * s / ( fm + 1.0E+00 );

  if ( m == 0 )
  {
    return value;
  }
  else if ( m == 1 )
  {
    value = - value - 1.0E+00 + 1.0E+00 / x;
    return value;
  }

  te = fm;
  t = 1.0E+00;
  s = t;
  mm1 = m - 1;
  for ( k = 1; k <= mm1; k++ )
  {
    fk = ( float ) ( k );
    te = - x * te / fk;
    t = te / ( fm - fk );
    s = s + t;
    if ( r4_abs ( t ) < eps * r4_abs ( s ) )
    {
      break;
    }
  }

  for ( k = 1; k <= m; k++ )
  {
    value = value + 1.0E+00 / ( float ) ( k );
  }

  if ( ( m % 2 ) == 1 )
  {
    sgng = - 1.0E+00;
  }
  else
  {
    sgng = + 1.0E+00;
  }

  alng = log ( value ) - r4_lngam ( fm + 1.0E+00 );

  if ( bot < alng )
  {
    value = sgng * exp ( alng );
  }
  else
  {
    value = 0.0E+00;
  }

  if ( s != 0.0E+00 )
  {
    value = value + r4_sign ( s ) 
      * exp ( - fm * alx + log ( r4_abs ( s ) / fm ) );
  }

  return value;
}
//****************************************************************************80

float r4_gmit ( float a, float x, float algap1, float sgngam, float alx )

//****************************************************************************80
//
//  Purpose:
//
//    R4_GMIT: Tricomi's incomplete gamma function for small X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, the parameter.
//
//    Input, float X, the argument.
//
//    Input, float ALGAP1, the logarithm of Gamma ( A + 1 ).
//
//    Input, float SGNGAM, the sign of Gamma ( A + 1 ).
//
//    Input, float ALX, the logarithm of X.
//
//    Output, float R4_GMIT, the Tricomi incomplete gamma function.
//
{
  float ae;
  float aeps;
  float alg2;
  float algs;
  static float bot = 0.0;
  bool converged;
  static float eps = 0.0;
  float fk;
  int k;
  int m;
  int ma;
  float s;
  float sgng2;
  float t;
  float te;
  float value;

  if ( eps == 0.0E+00 )
  {
    eps = 0.5E+00 * r4_mach ( 3 );
    bot = log ( r4_mach ( 1 ) );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_GMIT - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  if ( a < 0.0E+00 )
  {
    ma = ( int ) ( a - 0.5E+00 );
  }
  else
  {
    ma = ( int ) ( a + 0.5E+00 );
  }

  aeps = a - ( float ) ( ma );

  if ( a < - 0.5E+00 )
  {
    ae = aeps;
  }
  else
  {
    ae = a;
  }

  t = 1.0E+00;
  te = ae;
  s = t;
  converged = false;
  for ( k = 1; k <= 200; k++ )
  {
    fk = ( float ) ( k );
    te = - x * te / fk;
    t = te / ( ae + fk );
    s = s + t;
    if ( r4_abs ( t ) < eps * r4_abs ( s ) )
    {
      converged = true;
      break;
    }
  }

  if ( !converged )
  {
    cerr << "\n";
    cerr << "R4_GMIT - Fatal error!\n";
    cerr << "  No convergence in 200 iterations.\n";
    exit ( 1 );
  }

  if ( - 0.5E+00 <= a )
  {
    algs = - algap1 + log ( s );
    value = exp ( algs );
    return value;
  }

  algs = - r4_lngam ( 1.0E+00 + aeps ) + log ( s );
  s = 1.0E+00;
  m = - ma - 1;
  t = 1.0E+00;

  for ( k = 1; k <= m; k++ )
  {
    t = x * t / ( aeps - ( float ) ( m + 1 - k ) );
    s = s + t;
    if ( r4_abs ( t ) < eps * r4_abs ( s ) )
    {
      break;
    }
  }

  value = 0.0E+00;
  algs = - ( float ) ( ma ) * log ( x ) + algs;

  if ( s == 0.0E+00 || aeps == 0.0E+00 )
  {
    value = exp ( algs );
    return value;
  }

  sgng2 = sgngam * r4_sign ( s );
  alg2 = - x - algap1 + log ( r4_abs ( s ) );

  if ( bot < alg2 )
  {
    value = sgng2 * exp ( alg2 );
  }

  if ( bot < algs )
  {
    value = value + exp ( algs );
  }

  return value;
}
//****************************************************************************80

int r4_inits ( float dos[], int nos, float eta )

//****************************************************************************80
//
//  Purpose:
//
//    R4_INITS initializes a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, float DOS[NOS], the Chebyshev coefficients.
//
//    Input, int NOS, the number of coefficients.
//
//    Input, float ETA, the desired accuracy.
//
//    Output, int R4_INITS, the number of terms of the series needed
//    to ensure the requested accuracy.
//
{
  float err;
  int i;
  int value;

  if ( nos < 1 )
  {
    cerr << "\n";
    cerr << "R4_INITS - Fatal error!\n";
    cerr << "  Number of coefficients < 1.\n";
    exit ( 1 );
  }

  err = 0.0;

  for ( i = nos - 1; 0 <= i; i-- )
  {
    err = err + r4_abs ( dos[i] );
    if ( eta < err )
    {
      value = i + 1;
      return value;
    }
  }

  value = i;
  cerr << "\n";
  cerr << "R4_INITS - Warning!\n";
  cerr << "  ETA may be too small.\n";

  return value;
}
//****************************************************************************80

float r4_int ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_INT returns the integer part of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    2 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_INT, the integer part of X.
//
{
  int i;
  int ibase;
  int ipart;
  int iscale;
  static int npart = 0;
  float part;
  static float scale = 0.0;
  float value;
  static float xbig = 0.0;
  static float xmax = 0.0;
  float xscl;

  if ( npart == 0 )
  {
    ibase = i4_mach ( 10 );
    xmax = 1.0E+00 / r4_mach ( 4 );
    xbig = r4_min ( ( float ) ( i4_mach ( 9 ) ), xmax );
    iscale = i4_pow ( ibase, ( int ) ( log ( xbig ) 
      / log ( ( float ) ( ibase ) ) - 0.5E+00 ) );
    scale = ( float ) iscale;
    npart = log ( xmax ) / log ( scale ) + 1.0E+00;
  }

  if ( x < - xmax )
  {
    value = x;
  }
  else if ( x < - xbig )
  {
    xscl = - x;

    for ( i = 1; i <= npart; i++ )
    {
      xscl = xscl / scale;
    }

    value = 0.0E+00;
    for ( i = 1; i <= npart; i++ )
    {
      xscl = xscl * scale;
      ipart = ( int ) ( xscl );
      part = ( float ) ( ipart );
      xscl = xscl - part;
      value = value * scale + part;
    }
    value = - value;
  }
  else if ( x < + xbig )
  {
    value = ( int ) ( x );
  }
  else if ( x < + xmax )
  {
    xscl = x;

    for ( i = 1; i <= npart; i++ )
    {
      xscl = xscl / scale;
    }

    value = 0.0E+00;
    for ( i = 1; i <= npart; i++ )
    {
      xscl = xscl * scale;
      ipart = ( int ) ( xscl );
      part = ( float ) ( ipart );
      xscl = xscl - part;
      value = value * scale + part;
    }
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

void r4_knus ( float xnu, float x, float &bknu, float &bknu1, int &iswtch )

//****************************************************************************80
//
//  Purpose:
//
//    R4_KNUS computes a sequence of K Bessel functions.
//
//  Discussion:
//
//    This routine computes Bessel functions 
//      exp(x) * k-sub-xnu (x)  
//    and
//      exp(x) * k-sub-xnu+1 (x) 
//    for 0.0 <= xnu < 1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float XNU, the order parameter.
//
//    Input, float X, the argument.
//
//    Output, float &BKNU, &BKNU1, the two K Bessel functions.
//
//    Output, integer &ISWTCH, ?
//
{
  float a[15];
  float a0;
  static float aln2 = 0.69314718055994531E+00;
  static float alnbig = 0.0;
  static float alneps = 0.0;
  static float alnsml = 0.0;
  float alnz;
  float alpha[15];
  float an;
  float b0;
  float beta[15];
  float bknu0;
  float bknud;
  float bn;
  float c0;
  static float c0kcs[16] = {
    0.060183057242626108E+00,
   -0.15364871433017286E+00,
   -0.011751176008210492E+00,
   -0.000852487888919795E+00,
   -0.000061329838767496E+00,
   -0.000004405228124551E+00,
   -0.000000316312467283E+00,
   -0.000000022710719382E+00,
   -0.000000001630564460E+00,
   -0.000000000117069392E+00,
   -0.000000000008405206E+00,
   -0.000000000000603466E+00,
   -0.000000000000043326E+00,
   -0.000000000000003110E+00,
   -0.000000000000000223E+00,
   -0.000000000000000016E+00 };
  static float euler = 0.57721566490153286E+00;
  float expx;
  int i;
  int ii;
  int inu;
  int n;
  static int ntc0k = 0;
  int nterms;
  static int ntznu1 = 0;
  float p1;
  float p2;
  float p3;
  float qq;
  float result;
  static float sqpi2 = 1.2533141373155003E+00;
  float sqrtx;
  float v;
  float vlnz;
  float x2n;
  float x2tov;
  float xi;
  float xmu;
  static float xnusml = 0.0;
  static float xsml = 0.0;
  float z;
  static float znu1cs[12] = {
    0.20330675699419173E+00,
    0.14007793341321977E+00,
    0.007916796961001613E+00,
    0.000339801182532104E+00,
    0.000011741975688989E+00,
    0.000000339357570612E+00,
    0.000000008425941769E+00,
    0.000000000183336677E+00,
    0.000000000003549698E+00,
    0.000000000000061903E+00,
    0.000000000000000981E+00,
    0.000000000000000014E+00 };
  float ztov;

  if ( ntc0k == 0 )
  {
    ntc0k = r4_inits ( c0kcs, 16, 0.1E+00 * r4_mach ( 3 ) );
    ntznu1 = r4_inits ( znu1cs, 12, 0.1E+00 * r4_mach ( 3 ) );
    xnusml = sqrt ( r4_mach ( 3 ) / 8.0E+00 );
    xsml = 0.1E+00 * r4_mach ( 3 );
    alnsml = log ( r4_mach ( 1 ) );
    alnbig = log ( r4_mach ( 2 ) );
    alneps = log ( 0.1E+00 * r4_mach ( 3 ) );
  }

  if ( xnu < 0.0E+00 || 1.0E+00 <= xnu )
  {
    cerr << "\n";
    cerr << "R4_KNUS - Fatal error!\n";
    cerr << "  XNU < 0 or. 1 <= XNU.\n";
    exit ( 1 );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_KNUS - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  iswtch = 0;
//
//  X is small.  Compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
//  then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-.5,+.5)
//  then to (0., .5), because k of negative order (-nu) = k of positive
//  order (+nu).
//
  if ( x <= 2.0E+00 )
  {
    if ( 0.5E+00 < xnu )
    {
      v = 1.0E+00 - xnu;
    }
    else
    {
      v = xnu;
    }
//
//  Carefully find (x/2)^xnu and z^xnu where z = x*x/4.
//
    alnz = 2.0E+00 * ( log ( x ) - aln2 );

    if ( x <= xnu )
    {
      if ( alnbig < - 0.5E+00 * xnu * alnz - aln2 - log ( xnu ) )
      {
        cerr << "\n";
        cerr << "R4_KNUS - Fatal error!\n";
        cerr << "  Small X causing overflow.\n";
        exit ( 1 );
      }
    }
    vlnz = v * alnz;
    x2tov = exp ( 0.5E+00 * vlnz );

    if ( alnsml < vlnz )
    {
      ztov = x2tov * x2tov;
    }
    else
    {
      ztov = 0.0E+00;
    }

    a0 = 0.5E+00 * r4_gamma ( 1.0E+00 + v );
    b0 = 0.5E+00 * r4_gamma ( 1.0E+00 - v );
    c0 = - euler;

    if ( 0.5E+00 < ztov && xnusml < v )
    {
      c0 = - 0.75E+00 + r4_csevl ( ( 8.0E+00 * v ) * v - 1.0E+00, c0kcs, ntc0k );
    }

    if ( ztov <= 0.5E+00 )
    {
      alpha[0] = ( a0 - ztov * b0 ) / v;
    }
    else
    {
      alpha[0] = c0 - alnz * ( 0.75E+00 +
        r4_csevl ( vlnz / 0.35E+00 + 1.0E+00, znu1cs, ntznu1 ) ) * b0;
    }

    beta[0] = - 0.5E+00 * ( a0 + ztov * b0 );

    if ( xsml < x )
    {
      z = 0.25E+00 * x * x;
    }
    else
    {
      z = 0.0E+00;
    }

    nterms = ( int ) r4_max ( 2.0E+00, 11.0E+00 
      + ( 8.0E+00 * alnz - 25.19E+00 - alneps ) / ( 4.28E+00 - alnz ) );

    for ( i = 2; i <= nterms; i++ )
    {
      xi = ( float ) ( i - 1 );
      a0 = a0 / ( xi * ( xi - v) );
      b0 = b0 / ( xi * ( xi + v) );
      alpha[i-1] = ( alpha[i-2] + 2.0E+00 * xi * a0 ) 
        / ( xi * ( xi + v ) );
      beta[i-1] = ( xi - 0.5E+00 * v ) * alpha[i-1] - ztov * b0;
    }

    bknu = alpha[nterms-1];
    bknud = beta[nterms-1];
    for ( ii = 2; ii <= nterms; ii++ )
    {
      i = nterms + 1 - ii;
      bknu = alpha[i-1] + bknu * z;
      bknud = beta[i-1] + bknud * z;
    }

    expx = exp ( x );
    bknu = expx * bknu / x2tov;

    if ( alnbig < - 0.5E+00 * ( xnu + 1.0E+00 ) * alnz - 2.0E+00 * aln2 )
    {
      iswtch = 1;
      return;
    }

    bknud = expx * bknud * 2.0E+00 / ( x2tov * x );

    if ( xnu <= 0.5E+00 )
    {
      bknu1 = v * bknu / x - bknud;
      return;
    }

    bknu0 = bknu;
    bknu = - v * bknu / x - bknud;
    bknu1 = 2.0E+00 * xnu * bknu / x + bknu0;
  }
//
//  X is large.  Find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke's
//  rational expansion.
//
  else
  {
    sqrtx = sqrt ( x );

    if ( 1.0E+00 / xsml < x )
    {
      bknu = sqpi2 / sqrtx;
      bknu1 = bknu;
      return;
    }

    an = - 1.56E+00 + 4.0E+00 / x;
    bn = - 0.29E+00 - 0.22E+00 / x;
    nterms = i4_min ( 15, ( int ) r4_max ( 3.0E+00, an + bn * alneps ) );

    for ( inu = 1; inu <= 2; inu++ )
    {
      if ( inu == 1 )
      {
        if ( xnusml < xnu )
        {
          xmu = ( 4.0E+00 * xnu ) * xnu;
        }
        else
        {
          xmu = 0.0E+00;
        }
      }
      else
      {
        xmu = 4.0E+00 * pow ( r4_abs ( xnu ) + 1.0E+00, 2 );
      }

      a[0] = 1.0E+00 - xmu;
      a[1] = 9.0E+00 - xmu;
      a[2] = 25.0E+00 - xmu;

      if ( a[1] == 0.0E+00 )
      {
        result = sqpi2 * ( 16.0E+00 * x + xmu + 7.0E+00 )
          / ( 16.0E+00 * x * sqrtx );
      }
      else
      {
        alpha[0] = 1.0E+00;
        alpha[1] = ( 16.0E+00 * x + a[1] ) / a[1];
        alpha[2] = ( ( 768.0E+00 * x + 48.0E+00 * a[2] ) * x 
          + a[1] * a[2] ) / ( a[1] * a[2] );

        beta[0] = 1.0E+00;
        beta[1] = ( 16.0E+00 * x + ( xmu + 7.0E+00 ) ) / a[1];
        beta[2] = ( ( 768.0E+00 * x + 48.0E+00 * ( xmu + 23.0E+00 ) ) * x 
          + ( ( xmu + 62.0E+00 ) * xmu + 129.0E+00 ) ) / ( a[1] * a[2] );

        for ( i = 4; i <= nterms; i++ )
        {
          n = i - 1;
          x2n = ( float ) ( 2 * n - 1 );

          a[i-1] = pow ( x2n + 2.0E+00, 2 ) - xmu;
          qq = 16.0E+00 * x2n / a[i-1];
          p1 = - x2n * ( ( float ) ( 12 * n * n - 20 * n ) 
            - a[0] ) / ( ( x2n - 2.0E+00 ) * a[i-1] ) - qq * x;
          p2 = ( ( float ) ( 12 * n * n - 28 * n + 8 ) - a[0] ) 
            / a[i-1] - qq * x;
          p3 = - x2n * a[i-4] / ( ( x2n - 2.0E+00 ) * a[i-1] );

          alpha[i-1] = - p1 * alpha[i-2] 
                       - p2 * alpha[i-3] 
                       - p3 * alpha[i-4];

          beta[i-1] =  - p1 * beta[i-2] 
                       - p2 * beta[i-3]  
                       - p3 * beta[i-4];

        }
        result = sqpi2 * beta[nterms-1] / ( sqrtx * alpha[nterms-1] );
      }

      if ( inu == 1 )
      {
        bknu = result;
      }
      else
      {
        bknu1 = result;
      }
    }
  }
  return;
}
//****************************************************************************80

float r4_lbeta ( float a, float b )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LBETA evaluates the logarithm of the beta function of R4 arguments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, B, the arguments.
//
//    Output, float R4_LBETA, the logarithm of the beta function of A
//    and B.
//
{
  float corr;
  float p;
  float q;
  static float sq2pil = 0.91893853320467274E+00;
  float value;

  p = r4_min ( a, b );
  q = r4_max ( a, b );

  if ( p <=  0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_LBETA - Fatal error!\n";
    cerr << "  Both arguments must be greater than 0.\n";
    exit ( 1 );
  }
  else if ( p < 10.0E+00 && q <= 10.0E+00 )
  {
    value = log ( r4_gamma ( p ) * ( r4_gamma ( q ) / r4_gamma ( p + q ) ) );
  }
  else if ( p < 10.0E+00 )
  {
    corr = r4_lgmc ( q ) - r4_lgmc ( p + q );

    value = r4_lngam ( p ) + corr + p - p * log ( p + q ) +
      ( q - 0.5E+00 ) * r4_lnrel ( - p / ( p + q ) );
  }
  else
  {
    corr = r4_lgmc ( p ) + r4_lgmc ( q ) - r4_lgmc ( p + q );

    value = - 0.5E+00 * log ( q ) + sq2pil + corr 
      + ( p - 0.5E+00 ) * log ( p / ( p + q ) ) 
      + q * r4_lnrel ( - p / ( p + q ) );
  }
  return value;
}
//****************************************************************************80

void r4_lgams ( float x, float &algam, float &sgngam )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LGAMS evaluates the log of |gamma(x)| and sign, for an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float &ALGAM, the logarithm of the absolute value of
//    gamma ( X ).
//
//    Output, float &SGNGAM, the sign (+1 or -1 ) of gamma ( X ).
//
{
  int k;

  algam = r4_lngam ( x );
  sgngam = 1.0E+00;

  if ( x <= 0.0E+00 )
  {
    k = ( int ) ( r4_mod ( - r4_aint ( x ), 2.0E+00 ) + 0.1E+00 );

    if ( k == 0 )
    {
      sgngam = - 1.0E+00;
    }
  }
  return;
}
//****************************************************************************80

float r4_lgic ( float a, float x, float alx )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LGIC evaluates the log complementary incomplete gamma function for large X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, the parameter.
//
//    Input, float X, the argument.
//
//    Input, float ALX, the logarithm of X.
//
//    Output, float R4_LGIC, the log complementary incomplete gamma function.
//
{
  static float eps = 0.0;
  float fk;
  int k;
  float p;
  float r;
  float s;
  float t;
  float value;
  float xma;
  float xpa;

  if ( eps == 0.0E+00 )
  {
    eps = 0.5E+00 * r4_mach ( 3 );
  }

  xpa = x + 1.0E+00 - a;
  xma = x - 1.0E+00 - a;

  r = 0.0E+00;
  p = 1.0E+00;
  s = p;
  for ( k = 1; k <= 200; k++ )
  {
    fk = ( float ) ( k );
    t = fk * ( a - fk ) * ( 1.0E+00 + r );
    r = - t / ( ( xma + 2.0E+00 * fk ) * ( xpa + 2.0E+00 * fk ) + t );
    p = r * p;
    s = s + p;
    if ( r4_abs ( p ) < eps * s )
    {
      value = a * alx - x + log ( s / xpa );
      return value;
    }
  }
  cerr << "\n";
  cerr << "R4_LGIC - Fatal error!\n";
  cerr << "  No convergence in 200 iterations.\n";

  exit ( 1 );
}
//****************************************************************************80

float r4_lgit ( float a, float x, float algap1 )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LGIT evaluates the log of Tricomi's incomplete gamma function.
//
//  Discussion:
//
//    Perron's continued fraction is used for large X and X <= A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, the parameter.
//
//    Input, float X, the argument.
//
//    Input, float ALGAP1, the logarithm of A+1.
//
//    Output, float R4_LGIT, the log of Tricomi's incomplete
//    gamma function.
//
{
  float a1x;
  float ax;
  static float eps = 0.0;
  float fk;
  float hstar;
  int k;
  float p;
  float r;
  float s;
  static float sqeps = 0.0;
  float t;
  float value;

  if ( eps == 0.0E+00 )
  {
    eps = 0.5E+00 * r4_mach ( 3 );
    sqeps = sqrt ( r4_mach ( 4 ) );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_LGIT - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  if ( a < x )
  {
    cerr << "\n";
    cerr << "R4_LGIT - Fatal error!\n";
    cerr << "  A < X.\n";
    exit ( 1 );
  }

  ax = a + x;
  a1x = ax + 1.0E+00;
  r = 0.0E+00;
  p = 1.0E+00;
  s = p;

  for ( k = 1; k <= 200; k++ )
  {
    fk = ( float ) k;
    t = ( a + fk ) * x * ( 1.0E+00 + r );
    r = t / ( ( ax + fk ) * ( a1x + fk ) - t );
    p = r * p;
    s = s + p;
    if ( r4_abs ( p ) < eps * s )
    {
      hstar = 1.0E+00 - x * s / a1x;
      value = - x - algap1 - log ( hstar );
      return value;
    }
  }

  cerr << "\n";
  cerr << "R4_LGIT - Fatal error!\n";
  cerr << "  No convergence after 200 iterations.\n";
  exit ( 1 );
}
//****************************************************************************80

float r4_lgmc ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LGMC evaluates the log gamma correction factor for an R4 argument.
//
//  Discussion:
//
//    For 10 <= X, compute the log gamma correction factor so that
//
//      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
//                          + ( x - 0.5 ) * log ( x ) - x 
//                          + value ( x )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_LGMC, the correction factor.
//
{
  static float algmcs[6] = {
    0.166638948045186E+00,
   -0.0000138494817606E+00,
    0.0000000098108256E+00,
   -0.0000000000180912E+00,
    0.0000000000000622E+00,
   -0.0000000000000003E+00 };
  static int nalgm = 0;
  float value;
  static float xbig = 0.0;
  static float xmax = 0.0;

  if ( nalgm == 0 )
  {
    nalgm = r4_inits ( algmcs, 6, r4_mach ( 3 ) );
    xbig = 1.0E+00 / sqrt ( r4_mach ( 3 ) );
    xmax = exp ( r4_min ( log ( r4_mach ( 2 ) / 12.0E+00 ), 
      - log ( 12.0E+00 * r4_mach ( 1 ) ) ) );
  }

  if ( x < 10.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_LGMC - Fatal error!\n";
    cerr << "  X must be at least 10.\n";
    exit ( 1 );
  }
  else if ( x < xbig )
  {
    value = r4_csevl ( 2.0E+00 * ( 10.0E+00 / x ) 
      * ( 10.0E+00 / x ) - 1.0E+00, algmcs, nalgm ) / x;
  }
  else if ( x < xmax )
  {
    value = 1.0E+00 / ( 12.0E+00 * x );
  }
  else
  {
    value = 0.0E+00;
  }
  return value;
}
//****************************************************************************80

float r4_li ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LI evaluates the logarithmic integral for an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_LI, the logarithmic integral evaluated at X.
//
{
  static float sqeps = 0.0;
  float value;

  if ( sqeps == 0.0E+00 )
  {
    sqeps = sqrt ( r4_mach ( 3 ) );
  }

  if ( x < 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_LI - Fatal error!\n";
    cerr << "  Function undefined for X <= 0.\n";
    exit ( 1 );
  }

  if ( x == 0.0E+00 )
  {
    value = 0.0E+00;
    return value;
  }

  if ( x == 1.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_LI - Fatal error!\n";
    cerr << "  Function undefined for X = 1.\n";
    exit ( 1 );
  }

  if ( r4_abs ( 1.0E+00 - x ) < sqeps )
  {
    cerr << "\n";
    cerr << "R4_LI - Warning!\n";
    cerr << "  Answer less than half precision.\n";
    cerr << "  X is too close to 1.\n";
  }

  value = r4_ei ( log ( x ) );

  return value;
}
//****************************************************************************80

float r4_lngam ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LNGAM evaluates the log of the absolute value of gamma of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_LNGAM, the logarithm of the absolute value of
//    the gamma function of X.
//
{
  static float dxrel = 0.0;
  static float pi = 3.14159265358979324E+00;
  float sinpiy;
  static float sq2pil = 0.91893853320467274E+00;
  static float sqpi2l = 0.22579135264472743E+00;
  float value;
  static float xmax = 0.0;
  float y;

  if ( xmax == 0.0E+00 )
  {
    xmax = r4_mach ( 2 ) / log ( r4_mach ( 2 ) );
    dxrel = sqrt ( r4_mach ( 4 ) );
  }

  y = r4_abs ( x );

  if ( y <= 10.0E+00 )
  {
    value = log ( r4_abs ( r4_gamma ( x ) ) );
    return value;
  }

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R4_LNGAM - Fatal error!\n";
    cerr << "  Result overflows, |X| too big.\n";
    exit ( 1 );
  }

  if ( 0.0E+00 < x )
  {
    value = sq2pil + ( x - 0.5E+00 ) * log ( x ) - x + r4_lgmc ( y );
    return value;
  }

  sinpiy = r4_abs ( sin ( pi * y ) );

  if ( sinpiy == 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_LNGAM - Fatal error!\n";
    cerr << "  X is a negative integer.\n";
    exit ( 1 );
  }

  value = sqpi2l + ( x - 0.5E+00 ) * log ( y ) - x 
    - log ( sinpiy ) - r4_lgmc ( y );

  if ( r4_abs ( ( x - r4_aint ( x - 0.5E+00 ) ) * value / x ) < dxrel )
  {
    cerr << "\n";
    cerr << "R4_LNGAM - Warning!\n";
    cerr << "  Result is half precision because\n";
    cerr << "  X is too near a negative integer.\n";
  }

  return value;
}
//****************************************************************************80

float r4_lnrel ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LNREL evaluates log ( 1 + X ) for an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_LNREL, the value of log ( 1 + X ).
//
{
  static float alnrcs[23] = {
    1.0378693562743770E+00,
   -0.13364301504908918E+00,
    0.019408249135520563E+00,
   -0.003010755112753577E+00,
    0.000486946147971548E+00,
   -0.000081054881893175E+00,
    0.000013778847799559E+00,
   -0.000002380221089435E+00,
    0.000000416404162138E+00,
   -0.000000073595828378E+00,
    0.000000013117611876E+00,
   -0.000000002354670931E+00,
    0.000000000425227732E+00,
   -0.000000000077190894E+00,
    0.000000000014075746E+00,
   -0.000000000002576907E+00,
    0.000000000000473424E+00,
   -0.000000000000087249E+00,
    0.000000000000016124E+00,
   -0.000000000000002987E+00,
    0.000000000000000554E+00,
   -0.000000000000000103E+00,
    0.000000000000000019E+00 };
  static int nlnrel = 0;
  float value;
  static float xmin = 0.0;

  if ( nlnrel == 0 )
  {
    nlnrel = r4_inits ( alnrcs, 23, 0.1E+00 * r4_mach ( 3 ) );
    xmin = - 1.0E+00 + sqrt ( r4_mach ( 4 ) );
  }

  if ( x <= - 1.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_LNREL - Fatal error!\n";
    cerr << "  X <= -1.\n";
    exit ( 1 );
  }
  else if ( x < xmin )
  {
    cerr << "\n";
    cerr << "R4_LNREL - Warning!\n";
    cerr << "  Result is less than half precision.\n";
    cerr << "  X is too close to - 1.\n";
  }

  if ( r4_abs ( x ) <= 0.375E+00 )
  {
    value = x * ( 1.0E+00 - x * r4_csevl ( x / 0.375E+00, alnrcs, nlnrel ) );
  }
  else
  {
    value = log ( 1.0E+00 + x );
  }

  return value;
}
//****************************************************************************80

float r4_log ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LOG evaluates the logarithm of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the evaluation point.
//
//    Output, float R4_LOG, the logarithm of X.
//
{
  static float aln2 = 0.068147180559945309E+00;
  static float alncen[5] = {
    0.0E+00,
    +0.223143551314209755E+00,
    +0.405465108108164381E+00,
    +0.559615787935422686E+00,
    +0.693147180559945309E+00 };
  static float alncs[6] = {
    1.3347199877973882E+00,
    0.000693756283284112E+00,
    0.000000429340390204E+00,
    0.000000000289338477E+00,
    0.000000000000205125E+00,
    0.000000000000000150E+00 };
  static float center[4] = {
    1.0E+00,
    1.25E+00,
    1.50E+00,
    1.75E+00 };
  int n;
  static int nterms = 0;
  int ntrval;
  float t;
  float t2;
  float value;
  float xn;
  float y;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( alncs, 6, 28.9E+00 * r4_mach ( 3 ) );
  }

  if ( x <= 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_LOG - Fatal error!\n";
    cerr << "  X <= 0.0\n";
    exit ( 1 );
  }

  r4_upak ( x, y, n );

  xn = ( float ) ( n - 1 );
  y = 2.0E+00 * y;
  ntrval = ( int ) ( 4.0E+00 * y - 2.5E+00 );

  if ( ntrval == 5 )
  {
    t = ( ( y - 1.0E+00 ) - 1.0E+00 ) / ( y + 2.0E+00 );
  }
  else if ( ntrval < 5 )
  {
    t = ( y - center[ntrval-1] ) / ( y + center[ntrval-1] );
  }

  t2 = t * t;

  value = 0.625E+00 * xn + ( aln2 * xn + alncen[ntrval-1] 
    + 2.0E+00 * t 
    + t * t2 * r4_csevl ( 578.0E+00 * t2 - 1.0, alncs, nterms ) );

  return value;
}
//****************************************************************************80

float r4_log10 ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_LOG10 evaluates the logarithm, base 10, of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the evaluation point.
//
//    Output, float R4_LOG10, the logarithm, base 10, of X.
//
{
  static float aloge = 0.43429448190325182765E+00;
  float value;

  value = aloge * log ( x );

  return value;
}
//****************************************************************************80

float r4_mach ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MACH returns single precision real machine constants.
//
//  Discussion:
//
//    Assume that single precision real numbers are stored with a mantissa 
//    of T digits in base B, with an exponent whose value must lie 
//    between EMIN and EMAX.  Then for values of I between 1 and 5, 
//    R4_MACH will return the following values:
//
//      R4_MACH(1) = B^(EMIN-1), the smallest positive magnitude.
//      R4_MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
//      R4_MACH(3) = B^(-T), the smallest relative spacing.
//      R4_MACH(4) = B^(1-T), the largest relative spacing.
//      R4_MACH(5) = log10(B)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, Norman Schryer,
//    Algorithm 528,
//    Framework for a Portable Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, Number 2, June 1978, page 176-188.
//
//  Parameters:
//
//    Input, int I, chooses the parameter to be returned.
//    1 <= I <= 5.
//
//    Output, float R4_MACH, the value of the chosen parameter.
//
{
  float value;

  if ( i == 1 )
  {
    value = 1.1754944E-38;
  }
  else if ( i == 2 )
  {
    value = 3.4028235E+38;
  }
  else if ( i == 3 )
  {
    value = 5.9604645E-08;
  }
  else if ( i == 4 )
  {
    value = 1.1920929E-07;
  }
  else if ( i == 5 )
  {
    value = 0.3010300E+00;
  }
  else
  {
    cerr << "\n";
    cerr << "R4_MACH - Fatal error!\n";
    cerr << "  The input argument I is out of bounds.\n";
    cerr << "  Legal values satisfy 1 <= I <= 5.\n";
    cerr << "  I = " << i << "\n";
    value = 0.0;
    exit ( 1 );
  }
 
  return value;
}
//****************************************************************************80

void r4_machar ( long int *ibeta, long int *it, long int *irnd, long int *ngrd,
  long int *machep, long int *negep, long int *iexp, long int *minexp,
  long int *maxexp, float *eps, float *epsneg, float *xmin, float *xmax )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MACHAR computes machine constants for R4 arithmetic.
//
//  Discussion:
//
//    This routine determines the parameters of the floating-point 
//    arithmetic system specified below.  The determination of the first 
//    three uses an extension of an algorithm due to Malcolm, 
//    incorporating some of the improvements suggested by Gentleman and 
//    Marovich.  
//
//    A FORTRAN version of this routine appeared as ACM algorithm 665.
//
//    This routine is a C translation of the FORTRAN code, and appeared
//    as part of ACM algorithm 722.
//
//    An earlier version of this program was published in Cody and Waite.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2006
//
//  Author:
//
//    Original FORTRAN77 version by William Cody.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
//    machine parameters,
//    ACM Transactions on Mathematical Software,
//    Volume 14, Number 4, pages 303-311, 1988.
//
//    William Cody and W Waite,
//    Software Manual for the Elementary Functions,
//    Prentice Hall, 1980.
//
//    M Gentleman and S Marovich,
//    Communications of the ACM,
//    Volume 17, pages 276-277, 1974.
//
//    M. Malcolm,
//    Communications of the ACM,
//    Volume 15, pages 949-951, 1972.
//
//  Parameters:
//
//    Output, long int* IBETA, the radix for the floating-point representation.
//
//    Output, long int* IT, the number of base IBETA digits in the floating-point
//    significand.
//
//    Output, long int* IRND:
//    0, if floating-point addition chops.
//    1, if floating-point addition rounds, but not in the IEEE style.
//    2, if floating-point addition rounds in the IEEE style.
//    3, if floating-point addition chops, and there is partial underflow.
//    4, if floating-point addition rounds, but not in the IEEE style, and 
//      there is partial underflow.
//    5, if floating-point addition rounds in the IEEE style, and there is 
//      partial underflow.
//
//    Output, long int* NGRD, the number of guard digits for multiplication with
//    truncating arithmetic.  It is
//    0, if floating-point arithmetic rounds, or if it truncates and only 
//      IT base IBETA digits participate in the post-normalization shift of the
//      floating-point significand in multiplication;
//    1, if floating-point arithmetic truncates and more than IT base IBETA
//      digits participate in the post-normalization shift of the floating-point
//      significand in multiplication.
//
//    Output, long int* MACHEP, the largest negative integer such that
//      1.0 + ( float ) IBETA ^ MACHEP != 1.0, 
//    except that MACHEP is bounded below by - ( IT + 3 ).
//
//    Output, long int* NEGEPS, the largest negative integer such that
//      1.0 - ( float ) IBETA ) ^ NEGEPS != 1.0, 
//    except that NEGEPS is bounded below by - ( IT + 3 ).
//
//    Output, long int* IEXP, the number of bits (decimal places if IBETA = 10)
//    reserved for the representation of the exponent (including the bias or
//    sign) of a floating-point number.
//
//    Output, long int* MINEXP, the largest in magnitude negative integer such 
//    that
//      ( float ) IBETA ^ MINEXP 
//    is positive and normalized.
//
//    Output, long int* MAXEXP, the smallest positive power of BETA that overflows.
// 
//    Output, float* EPS, the smallest positive floating-point number such
//    that  
//      1.0 + EPS != 1.0. 
//    in particular, if either IBETA = 2  or IRND = 0, 
//      EPS = ( float ) IBETA ^ MACHEP.
//    Otherwise,  
//      EPS = ( ( float ) IBETA ^ MACHEP ) / 2.
//
//    Output, float* EPSNEG, a small positive floating-point number such that
//      1.0 - EPSNEG != 1.0. 
//    In particular, if IBETA = 2 or IRND = 0, 
//    EPSNEG = ( float ) IBETA ^ NEGEPS.
//    Otherwise,  
//      EPSNEG = ( float ) IBETA ^ NEGEPS ) / 2.  
//    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
//    smallest number that can alter 1.0 by subtraction.
//
//    Output, float* XMIN, the smallest non-vanishing normalized floating-point
//    power of the radix:
//      XMIN = ( float ) IBETA ^ MINEXP
//
//    Output, float* XMAX, the largest finite floating-point number.  In
//    particular,
//      XMAX = ( 1.0 - EPSNEG ) * ( float ) IBETA ^ MAXEXP
//    On some machines, the computed value of XMAX will be only the second, 
//    or perhaps third, largest number, being too small by 1 or 2 units in 
//    the last digit of the significand.
//
{
  float a;
  float b;
  float beta;
  float betah;
  float betain;
  int i;
  int itmp;
  int iz;
  int j;
  int k;
  int mx;
  int nxres;
  float one;
  float t;
  float tmp;
  float tmp1;
  float tmpa;
  float two;
  float y;
  float z;
  float zero;

  (*irnd) = 1;
  one = (float) (*irnd);
  two = one + one;
  a = two;
  b = a;
  zero = 0.0e0;
//
//  Determine IBETA and BETA ala Malcolm.
//
  tmp = ( ( a + one ) - a ) - one;

  while ( tmp == zero )
  {
    a = a + a;
    tmp = a + one;
    tmp1 = tmp - a;
    tmp = tmp1 - one;
  }

  tmp = a + b;
  itmp = ( int ) ( tmp - a );

  while ( itmp == 0 )
  {
    b = b + b;
    tmp = a + b;
    itmp = ( int ) ( tmp - a );
  }

  *ibeta = itmp;
  beta = ( float ) ( *ibeta );
//
//  Determine IRND, IT.
//
  ( *it ) = 0;
  b = one;
  tmp = ( ( b + one ) - b ) - one;

  while ( tmp == zero )
  {
    *it = *it + 1;
    b = b * beta;
    tmp = b + one;
    tmp1 = tmp - b;
    tmp = tmp1 - one;
  }

  *irnd = 0;
  betah = beta / two;
  tmp = a + betah;
  tmp1 = tmp - a;

  if ( tmp1 != zero )
  {
    *irnd = 1;
  }

  tmpa = a + beta;
  tmp = tmpa + betah;

  if ( ( *irnd == 0 ) && ( tmp - tmpa != zero ) )
  {
    *irnd = 2;
  }
//
//  Determine NEGEP, EPSNEG.
//
  (*negep) = (*it) + 3;
  betain = one / beta;
  a = one;
 
  for ( i = 1; i <= (*negep); i++ )
  {
    a = a * betain;
  }
 
  b = a;
  tmp = ( one - a );
  tmp = tmp - one;

  while ( tmp == zero )
  {
    a = a * beta;
    *negep = *negep - 1;
    tmp1 = one - a;
    tmp = tmp1 - one;
  }

  (*negep) = -(*negep);
  (*epsneg) = a;
//
//  Determine MACHEP, EPS.
//

  (*machep) = -(*it) - 3;
  a = b;
  tmp = one + a;

  while ( tmp - one == zero)
  {
    a = a * beta;
    *machep = *machep + 1;
    tmp = one + a;
  }

  *eps = a;
//
//  Determine NGRD.
//
  (*ngrd) = 0;
  tmp = one + *eps;
  tmp = tmp * one;

  if ( ( (*irnd) == 0 ) && ( tmp - one ) != zero )
  {
    (*ngrd) = 1;
  }
//
//  Determine IEXP, MINEXP and XMIN.
//
//  Loop to determine largest I such that (1/BETA) ** (2**(I))
//  does not underflow.  Exit from loop is signaled by an underflow.
//

  i = 0;
  k = 1;
  z = betain;
  t = one + *eps;
  nxres = 0;

  for ( ; ; )
  {
    y = z;
    z = y * y;
//
//  Check for underflow
//

    a = z * one;
    tmp = z * t;

    if ( ( a + a == zero ) || ( r4_abs ( z ) > y ) )
    {
      break;
    }

    tmp1 = tmp * betain;

    if ( tmp1 * beta == z )
    {
      break;
    }

    i = i + 1;
    k = k + k;
  }
//
//  Determine K such that (1/BETA)**K does not underflow.
//  First set  K = 2 ** I.
//
  (*iexp) = i + 1;
  mx = k + k;

  if ( *ibeta == 10 )
  {
//
//  For decimal machines only
//

    (*iexp) = 2;
    iz = *ibeta;
    while ( k >= iz )
    {
      iz = iz * ( *ibeta );
      (*iexp) = (*iexp) + 1;
    }
    mx = iz + iz - 1;
  }
 
//
//  Loop to determine MINEXP, XMIN.
//  Exit from loop is signaled by an underflow.
//
  for ( ; ; )
  {
    (*xmin) = y;
    y = y * betain;
    a = y * one;
    tmp = y * t;
    tmp1 = a + a;

    if ( ( tmp1 == zero ) || ( r4_abs ( y ) >= ( *xmin ) ) )
    {
      break;
    }

    k = k + 1;
    tmp1 = tmp * betain;
    tmp1 = tmp1 * beta;

    if ( ( tmp1 == y ) && ( tmp != y ) )
    {
      nxres = 3;
      *xmin = y;
      break;
    }

  }

  (*minexp) = -k;

//
//  Determine MAXEXP, XMAX.
//
  if ( ( mx <= k + k - 3 ) && ( ( *ibeta ) != 10 ) )
  {
    mx = mx + mx;
    (*iexp) = (*iexp) + 1;
  }

  (*maxexp) = mx + (*minexp);
//
//  Adjust IRND to reflect partial underflow.
//
  (*irnd) = (*irnd) + nxres;
//
//  Adjust for IEEE style machines.
//
  if ( ( *irnd) >= 2 )
  {
    (*maxexp) = (*maxexp) - 2;
  }
//
//  Adjust for machines with implicit leading bit in binary
//  significand and machines with radix point at extreme
//  right of significand.
//
  i = (*maxexp) + (*minexp);

  if ( ( ( *ibeta ) == 2 ) && ( i == 0 ) )
  {
    (*maxexp) = (*maxexp) - 1;
  }

  if ( i > 20 )
  {
    (*maxexp) = (*maxexp) - 1;
  }

  if ( a != y )
  {
    (*maxexp) = (*maxexp) - 2;
  }

  (*xmax) = one - (*epsneg);
  tmp = (*xmax) * one;

  if ( tmp != (*xmax) )
  {
    (*xmax) = one - beta * (*epsneg);
  }

  (*xmax) = (*xmax) / ( beta * beta * beta * (*xmin) );
  i = (*maxexp) + (*minexp) + 3;

  if ( i > 0 )
  {
 
    for ( j = 1; j <= i; j++ )
    {
      if ( (*ibeta) == 2 )
      {
        (*xmax) = (*xmax) + (*xmax);
      }
      if ( (*ibeta) != 2 )
      {
        (*xmax) = (*xmax) * beta;
      }
    }

  }
  return;
}
//****************************************************************************80

float r4_max ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MAX returns the maximum of two R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, Y, the quantities to compare.
//
//    Output, float R4_MAX, the maximum of X and Y.
//
{
  float value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

float r4_min ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MIN returns the minimum of two R4's..
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, Y, the quantities to compare.
//
//    Output, float R_MIN, the minimum of X and Y.
//
{
  float value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

float r4_mod ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MOD returns the remainder of R4 division.
//
//  Discussion:
//
//    If
//      REM = R8_MOD ( X, Y )
//      RMULT = ( X - REM ) / Y
//    then
//      X = Y * RMULT + REM
//    where REM has the same sign as X, and abs ( REM ) < Y.
//
//  Example:
//
//        X         Y     R4_MOD   R4_MOD  Factorization
//
//      107        50       7     107 =  2 *  50 + 7
//      107       -50       7     107 = -2 * -50 + 7
//     -107        50      -7    -107 = -2 *  50 - 7
//     -107       -50      -7    -107 =  2 * -50 - 7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number to be divided.
//
//    Input, float Y, the number that divides X.
//
//    Output, float R4_MOD, the remainder when X is divided by Y.
//
{
  float value;

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R4_MOD - Fatal error!\n";
    cerr << "  R4_MOD ( X, Y ) called with Y = " << y << "\n";
    exit ( 1 );
  }

  value = x - ( ( float ) ( ( int ) ( x / y ) ) ) * y;

  if ( x < 0.0 && 0.0 < value )
  {
    value = value - r4_abs ( y );
  }
  else if ( 0.0 < x && value < 0.0 )
  {
    value = value + r4_abs ( y );
  }

  return value;
}
//****************************************************************************80

float r4_mop ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MOP returns the I-th power of -1 as an R4 value.
//
//  Discussion:
//
//    An R4 is a float value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, float R4_MOP, the I-th power of -1.
//
{
  float value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}
//****************************************************************************80

float r4_pak ( float y, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4_PAK packs a base 2 exponent into an R4.
//
//  Discussion:
//
//    This routine is almost the inverse of R4_UPAK.  It is not exactly 
//    the inverse, because abs(x) need not be between 0.5 and 1.0.  
//    If both R4_PAK and 2.0^n were known to be in range, we could compute
//    R4_PAK = x * 2.0^n .
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, float Y, the mantissa.
//
//    Input, integer N, the exponent.
//
//    Output, float R4_PAK, the packed value.
//
{
  static float aln210 = 3.321928094887362E+00;
  float aln2b;
  static int nmax = 0;
  static int nmin = 0;
  int nsum;
  int ny;
  float value;

  if ( nmin == 0 )
  {
    aln2b = 1.0E+00;
    if ( i4_mach ( 10 ) != 2 )
    {
      aln2b = r4_mach ( 5 ) * aln210;
    }
    nmin = aln2b * ( float ) ( i4_mach ( 12 ) );
    nmax = aln2b * ( float ) ( i4_mach ( 13 ) );
  }

  r4_upak ( y, value, ny );

  nsum = n + ny;

  if ( nsum < nmin )
  {
    cerr << "\n";
    cerr << "R4_PAK - Warning!\n";
    cerr << "  Packed number underflows.\n";
    value = 0.0E+00;
    return value;
  }

  if ( nmax < nsum )
  {
    cerr << "\n";
    cerr << "R4_PAK - Fatal error!\n";
    cerr << "  Packed number overflows.\n";
    exit ( 1 );
  }

  while ( nsum < 0 )
  {
    value = 0.5E+00 * value;
    nsum = nsum + 1;
  }

  while ( 0 < nsum )
  {
    value = 2.0E+00 * value;
    nsum = nsum - 1;
  }

  return value;
}
//****************************************************************************80

float r4_poch ( float a, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_POCH evaluates Pochhammer's function of R4 arguments.
//
//  Discussion:
//
//    POCH ( A, X ) = Gamma ( A + X ) / Gamma ( A ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, X, the arguments.
//
//    Output, float R4_POCH, the Pochhammer function of A and X.
//
{
  float absa;
  float absax;
  float alnga;
  float alngax;
  float ax;
  float b;
  float cospia;
  float cospix;
  float den;
  static float eps = 0.0;
  float err;
  float errpch;
  int i;
  int n;
  static float pi = 3.141592653589793238E+00;
  float value;
  float sgnga;
  float sgngax;
  float sinpia;
  float sinpix;
  static float sqeps = 0.0;

  if ( eps == 0.0E+00 )
  {
    eps = r4_mach ( 4 );
    sqeps = sqrt ( eps );
  }

  ax = a + x;

  if ( ax <= 0.0E+00 && r4_aint ( ax ) == ax )
  {
    if ( 0.0E+00 < a || r4_aint ( a ) != a )
    {
      cerr << "\n";
      cerr << "R4_POCH - Fatal error!\n";
      cerr << "  A + X is nonpositive integer,\n";
      cerr << "  but A is not.\n";
      exit ( 1 );
    }
//
//  We know here that both A+X and A are non-positive integers.
//
    if ( x == 0.0E+00 )
    {
      value = 1.0E+00;
    }
    else if ( - 20.0E+00 <= r4_min ( a + x, a ) )
    {
      n = ( int ) ( x );
      value = r4_mop ( n ) * r4_fac ( - ( int ) ( a ) ) 
        / r4_fac ( - ( int ) ( a ) - n );
    }
    else
    {
      n = ( int ) ( x );
      value = r4_mop ( n ) * exp ( ( a - 0.5E+00 ) 
        * r4_lnrel ( x / ( a - 1.0E+00 ) ) 
        + x * log ( - a + 1.0E+00 - x ) - x 
        + r4_lgmc ( - a + 1.0E+00 ) 
        - r4_lgmc ( - a - x + 1.0E+00 ) );
    }

    return value;
  }
//
//  Here we know A+X is not zero or a negative integer.
//
  if ( a <= 0.0E+00 && r4_aint ( a ) == a )
  {
    value = 0.0E+00;
    return value;
  }

  n = ( int ) r4_abs ( x );
//
//  X is a small non-positive integer, presummably a common case.
//
  if ( ( float ) ( n ) == x && n <= 20 )
  {
    value = 1.0E+00;
    for ( i = 1; i <= n; i++ )
    {
      value = value * ( a + ( float ) ( i - 1 ) );
    }
    return value;
  }

  absax = r4_abs ( a + x );
  absa = r4_abs ( a );

  if ( r4_max ( absax, absa ) <= 20.0E+00 )
  {
    value = r4_gamma ( a + x ) * r4_gamr ( a );
    return value;
  }

  if ( 0.5E+00 * absa < r4_abs ( x ) )
  {
    r4_lgams ( a + x, alngax, sgngax );
    r4_lgams ( a, alnga, sgnga );
    value = sgngax * sgnga * exp ( alngax - alnga );
    return value;
  }
//
//  Here abs(x) is small and both abs(a+x) and abs(a) are large.  Thus,
//  a+x and a must have the same sign.  For negative a, we use
//  gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) *
//  sin(pi*a)/sin(pi*(a+x))
//
  if ( a < 0.0E+00 )
  {
    b = - a - x + 1.0E+00;
  }
  else
  {
    b = a;
  }

  value = exp ( ( b - 0.5E+00 ) * r4_lnrel ( x / b ) 
    + x * log ( b + x ) - x + r4_lgmc ( b + x ) - r4_lgmc ( b ) );

  if ( 0.0E+00 <= a || value == 0.0E+00 )
  {
    return value;
  }

  cospix = cos ( pi * x );
  sinpix = sin ( pi * x );
  cospia = cos ( pi * a );
  sinpia = sin ( pi * a );

  errpch = r4_abs ( x ) * ( 1.0E+00 + log ( b ) );
  den = cospix + cospia * sinpix / sinpia;
  err = ( r4_abs ( x ) * ( r4_abs ( sinpix ) 
    + r4_abs ( cospia * cospix / sinpia ) )
    + r4_abs ( a * sinpix ) / ( sinpia * sinpia ) ) * pi;
  err = errpch + err / r4_abs ( den );

  value = value / den;

  return value;
}
//****************************************************************************80

float r4_poch1 ( float a, float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_POCH1 evaluates a quantity related to Pochhammer's symbol.
//
//  Discussion:
//
//    Evaluate a generalization of Pochhammer's symbol for special
//    situations that require especially accurate values when x is small in
//      poch1(a,x) = (poch(a,x)-1)/x
//                 = (gamma(a+x)/gamma(a) - 1.0)/x .
//    This specification is particularly suited for stably computing
//    expressions such as
//      (gamma(a+x)/gamma(a) - gamma(b+x)/gamma(b))/x
//           = poch1(a,x) - poch1(b,x)
//    Note that poch1(a,0.0) = psi(a)
//
//    When abs(x) is so small that substantial cancellation will occur if
//    the straightforward formula is used, we  use an expansion due
//    to fields and discussed by y. l. luke, the special functions and their
//    approximations, vol. 1, academic press, 1969, page 34.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float A, the parameter.
//
//    Input, float X, the evaluation point.
//
//    Output, float R4_POCH1, the value of the function.
//
{
  float absa;
  float absx;
  static float alneps = 0.0;
  float alnvar;
  float b;
  static float bern[9] = {
   0.83333333333333333E-01,
  -0.13888888888888889E-02,
   0.33068783068783069E-04,
  -0.82671957671957672E-06,
   0.20876756987868099E-07,
  -0.52841901386874932E-09,
   0.13382536530684679E-10,
  -0.33896802963225829E-12,
   0.85860620562778446E-14 };
  float binv;
  float bp;
  float gbern[10];
  float gbk;
  int i;
  int ii;
  int incr;
  int j;
  int k;
  int ndx;
  int nterms;
  static float pi = 3.14159265358979324E+00;
  float poly1;
  float q;
  float rho;
  float sinpx2;
  float sinpxx;
  static float sqtbig = 0.0;
  float term;
  float trig;
  float value;
  float var;
  float var2;

  if ( sqtbig == 0.0E+00 )
  {
    sqtbig = 1.0E+00 / sqrt ( 24.0E+00 * r4_mach ( 1 ) );
    alneps = log ( r4_mach ( 3 ) );
  }

  if ( x == 0.0E+00 )
  {
    value = r4_psi ( a );
    return value;
  }

  absx = r4_abs ( x );
  absa = r4_abs ( a );

  if ( 0.1E+00 * absa < absx || 0.1E+00 < absx * log ( r4_max ( absa, 2.0E+00 ) ) )
  {
    value = r4_poch ( a, x );
    value = ( value - 1.0E+00 ) / x;
    return value;
  }

  if ( a < - 0.5E+00 )
  {
    bp = 1.0E+00 - a - x;
  }
  else
  {
    bp = a;
  }

  if ( bp < 10.0E+00 )
  {
    incr = r4_aint ( 11.0E+00 - bp );
  }
  else
  {
    incr = 0;
  }

  b = bp + ( float ) ( incr );

  var = b + 0.5E+00 * ( x - 1.0E+00 );
  alnvar = log ( var );
  q = x * alnvar;
  poly1 = 0.0E+00;

  if ( var < sqtbig )
  {
    var2 = 1.0E+00 / var / var;

    rho = 0.5E+00 * ( x + 1.0E+00 );
    gbern[0] = 1.0E+00;
    gbern[1] = - rho / 12.0E+00;
    term = var2;
    poly1 = gbern[1] * term;

    nterms = ( int ) ( - 0.5E+00 * alneps / alnvar + 1.0E+00 );

    if ( 9 < nterms )
    {
      cerr << "\n";
      cerr << "R4_POCH1 - Fatal error!\n";
      cerr << "  9 < NTERMS.\n";
      exit ( 1 );
    } 

    for ( k = 2; k <= nterms; k++ )
    {
      gbk = 0.0E+00;
      for ( j = 1; j <= k; j++ )
      {
        ndx = k - j + 1;
        gbk = gbk + bern[ndx-1] * gbern[j-1];
      }
      gbern[k] = - rho * gbk / ( float ) ( k );
      term = term * ( ( float ) ( 2 * k - 2 ) - x ) 
        * ( ( float ) ( 2 * k - 1 ) - x ) * var2;
      poly1 = poly1 + gbern[k] * term;
    }
  }

  poly1 = ( x - 1.0E+00 ) * poly1;
  value = r4_exprel ( q ) * ( alnvar + q * poly1 ) + poly1;
//
//  We have poch1(b,x), but bp is small, so we use backwards recursion
//  to obtain poch1(bp,x).
//
  for ( ii = 1; ii <= incr; ii++ )
  {
    i = incr - ii;
    binv = 1.0E+00 / ( bp + ( float ) ( i ) );
    value = ( value - binv ) / ( 1.0E+00 + x * binv );
  }

  if ( bp == a )
  {
    return value;
  }
//
//  We have poch1(bp,x), but a is lt -0.5.  We therefore use a reflection
//  formula to obtain poch1(a,x).
//
  sinpxx = sin ( pi * x ) / x;
  sinpx2 = sin ( 0.5E+00 * pi * x );
  trig = sinpxx * r4_cot ( pi * b ) - 2.0E+00 * sinpx2 * ( sinpx2 / x );

  value = trig + ( 1.0E+00 + x * trig ) * value;

  return value;
}
//****************************************************************************80

float r4_power ( float a, float b )

//****************************************************************************80
//
//  Purpose:
//
//    R4_POWER evaluates A^B.
//
//  Discussion:
//
//    Defining this function obviates the idiotic compiler complaints
//    that arise when the arguments of POW are not both simple floats.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, the base.
//
//    Input, float B, the exponent.
//
//    Output, float R4_POWER, the value of A^B.
//
{
  float value;

  value = pow ( a, b );

  return value;
}
//****************************************************************************80

float r4_psi ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_PSI evaluates the psi function of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_PSI, the psi function of X.
//
 {
  static float apsics[16] = {
   -0.0204749044678185E+00,
   -0.0101801271534859E+00,
    0.0000559718725387E+00,
   -0.0000012917176570E+00,
    0.0000000572858606E+00,
   -0.0000000038213539E+00,
    0.0000000003397434E+00,
   -0.0000000000374838E+00,
    0.0000000000048990E+00,
   -0.0000000000007344E+00,
    0.0000000000001233E+00,
   -0.0000000000000228E+00,
    0.0000000000000045E+00,
   -0.0000000000000009E+00,
    0.0000000000000002E+00,
   -0.0000000000000000E+00 };
  float aux;
  static float dxrel = 0.0;
  int i;
  int n;
  static int ntapsi = 0;
  static int ntpsi = 0;
  static float pi = 3.14159265358979324E+00;
  static float psics[23] = {
   -0.038057080835217922E+00,
    0.49141539302938713E+00,
   -0.056815747821244730E+00,
    0.008357821225914313E+00,
   -0.001333232857994342E+00,
    0.000220313287069308E+00,
   -0.000037040238178456E+00,
    0.000006283793654854E+00,
   -0.000001071263908506E+00,
    0.000000183128394654E+00,
   -0.000000031353509361E+00,
    0.000000005372808776E+00,
   -0.000000000921168141E+00,
    0.000000000157981265E+00,
   -0.000000000027098646E+00,
    0.000000000004648722E+00,
   -0.000000000000797527E+00,
    0.000000000000136827E+00,
   -0.000000000000023475E+00,
    0.000000000000004027E+00,
   -0.000000000000000691E+00,
    0.000000000000000118E+00,
   -0.000000000000000020E+00 };
  float value;
  static float xbig = 0.0;
  float y;

  if ( ntpsi == 0 )
  {
    ntpsi = r4_inits ( psics, 23, 0.1E+00 * r4_mach ( 3 ) );
    ntapsi = r4_inits ( apsics, 16, 0.1E+00 * r4_mach ( 3 ) );
    xbig = 1.0E+00 / sqrt ( r4_mach ( 3 ) );
    dxrel = sqrt ( r4_mach ( 4 ) );
  }

  y = r4_abs ( x );

  if ( y < 2.0E+00 )
  {
    n = ( int ) ( x );
    if ( x < 0.0E+00 ) 
    {
      n = n - 1;
    }
    y = x - ( float ) ( n );
    n = n - 1;
    value = r4_csevl ( 2.0E+00 * y - 1.0E+00, psics, ntpsi );

    if ( n == 0 )
    {
      return value;
    }

    n = - n;

    if ( x == 0.0E+00 )
    {
      cerr << "\n";
      cerr << "R4_PSI - Fatal error!\n";
      cerr << "  X is zero.\n";
      exit ( 1 );
    }

    if ( x < 0.0E+00 && x + ( float ) ( n - 2 ) == 0.0E+00 )
    {
      cerr << "\n";
      cerr << "R4_PSI - Fatal error!\n";
      cerr << "  X is a negative integer.\n";
      exit ( 1 );
    }

    if ( x < - 0.5E+00 &&
      r4_abs ( ( x - r4_aint ( x - 0.5E+00 ) ) / x ) < dxrel )
    {
      cerr << "\n";
      cerr << "R4_PSI - Warning!\n";
      cerr << "  Answer is less than half precision\n";
      cerr << "  because X is near a negative integer.\n";
    }

    for ( i = 1; i <= n; i++ )
    {
      value = value - 1.0E+00 / ( x + ( float ) ( i - 1 ) );
    }
  }
  else
  {
    if ( y < xbig )
    {
      aux = r4_csevl ( 8.0E+00 / y / y - 1.0E+00, apsics, ntapsi );
    }
    else
    {
      aux = 0.0E+00;
    }

    if ( x < 0.0E+00 )
    {
      value = log ( r4_abs ( x ) ) - 0.5E+00 / x + aux 
        - pi * r4_cot ( pi * x );
    }
    else if ( 0.0E+00 < x )
    {
      value = log ( x ) - 0.5E+00 / x + aux;
    }
  }
  return value;
}
//****************************************************************************80

float r4_rand ( float r )

//****************************************************************************80
//
//  Purpose:
//
//    R4_RAND is a portable pseudorandom number generator.
//
//  Discussion:
//
//    This pseudo-random number generator is portable amoung a wide
//    variety of computers.  It is undoubtedly not as good as many
//    readily available installation dependent versions, and so this
//    routine is not recommended for widespread usage.  Its redeeming
//    feature is that the exact same random numbers (to within final round-
//    off error) can be generated from machine to machine.  Thus, programs
//    that make use of random numbers can be easily transported to and
//    checked in a new environment.
//
//    The random numbers are generated by the linear congruential
//    method described by Knuth in seminumerical methods (p.9),
//    addison-wesley, 1969.  Given the i-th number of a pseudo-random
//    sequence, the i+1 -st number is generated from
//      x(i+1) = (a*x(i) + c) mod m,
//    where here m = 2^22 = 4194304, c = 1731 and several suitable values
//    of the multiplier a are discussed below.  Both the multiplier a and
//    random number x are represented in double precision as two 11-bit
//    words.  The constants are chosen so that the period is the maximum
//    possible, 4194304.
//
//    In order that the same numbers be generated from machine to
//    machine, it is necessary that 23-bit ints be reducible modulo
//    2^11 exactly, that 23-bit ints be added exactly, and that 11-bit
//    ints be multiplied exactly.  Furthermore, if the restart option
//    is used (where r is between 0 and 1), then the product r*2^22 =
//    r*4194304 must be correct to the nearest int.
//
//    The first four random numbers should be 
//
//      0.0004127026,
//      0.6750836372, 
//      0.1614754200, 
//      0.9086198807.
//
//    The tenth random number is 
//
//      0.5527787209.
//
//    The hundredth random number is 
//
//      0.3600893021.  
//
//    The thousandth number should be 
//
//      0.2176990509.
//
//    In order to generate several effectively independent sequences
//    with the same generator, it is necessary to know the random number
//    for several widely spaced calls.  The I-th random number times 2^22,
//    where I=K*P/8 and P is the period of the sequence (P = 2^22), is
//    still of the form L*P/8.  In particular we find the I-th random
//    number multiplied by 2^22 is given by
//      I   =  0  1*p/8  2*p/8  3*p/8  4*p/8  5*p/8  6*p/8  7*p/8  8*p/8
//      RAND=  0  5*p/8  2*p/8  7*p/8  4*p/8  1*p/8  6*p/8  3*p/8  0
//    thus the 4*P/8 = 2097152 random number is 2097152/2^22.
//
//    Several multipliers have been subjected to the spectral test
//    (see Knuth, p. 82).  Four suitable multipliers roughly in order of
//    goodness according to the spectral test are
//      3146757 = 1536*2048 + 1029 = 2^21 + 2^20 + 2^10 + 5
//      2098181 = 1024*2048 + 1029 = 2^21 + 2^10 + 5
//      3146245 = 1536*2048 +  517 = 2^21 + 2^20 + 2^9 + 5
//      2776669 = 1355*2048 + 1629 = 5^9 + 7^7 + 1
//
//    In the table below log10(NU(I)) gives roughly the number of
//    random decimal digits in the random numbers considered I at a time.
//    C is the primary measure of goodness.  In both cases bigger is better.
//
//                     log10 nu(i)              c(i)
//         a       i=2  i=3  i=4  i=5    i=2  i=3  i=4  i=5
//
//      3146757    3.3  2.0  1.6  1.3    3.1  1.3  4.6  2.6
//      2098181    3.3  2.0  1.6  1.2    3.2  1.3  4.6  1.7
//      3146245    3.3  2.2  1.5  1.1    3.2  4.2  1.1  0.4
//      2776669    3.3  2.1  1.6  1.3    2.5  2.0  1.9  2.6
//     best
//      possible   3.3  2.3  1.7  1.4    3.6  5.9  9.7  14.9
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, real R, determines the action.
//    * R = 0.0, the next random number of the sequence is generated.
//    * R < 0.0, the last generated number will be returned for
//    possible use in a restart procedure.
//    * R > 0.0, the sequence of random numbers will start with the
//    seed ( R mod 1 ).  This seed is also returned as the value of
//    R4_RAND provided the arithmetic is done exactly.
//
//    Output, real R4_RAND, a pseudo-random number between 0.0 and 1.0.
//
{
  static int ia0 = 1029;
  static int ia1 = 1536;
  static int ia1ma0 = 507;
  static int ic = 1731;
  static int ix0 = 0;
  static int ix1 = 0;
  int iy0;
  int iy1;
  float value;

  if ( r == 0.0E+00 )
  {
    iy0 = ia0 * ix0;
    iy1 = ia1 * ix1 + ia1ma0 * ( ix0 - ix1 ) + iy0;
    iy0 = iy0 + ic;
    ix0 = ( iy0 % 2048 );
    iy1 = iy1 + ( iy0 - ix0 ) / 2048;
    ix1 = ( iy1 % 2048 );
  }

  if ( 0.0 < r )
  {
    ix1 = ( int ) ( r4_mod ( r, 1.0E+00 ) * 4194304.0 + 0.5E+00 );
    ix0 = ( ix1 % 2048 );
    ix1 = ( ix1 - ix0 ) / 2048;
  }

  value = ( float ) ( ix1 * 2048 + ix0 );
  value = value / 4194304.0E+00;
 
  return value;
}
//****************************************************************************80

float r4_randgs ( float xmean, float sd )

//****************************************************************************80
//
//  Purpose:
//
//    R4_RANDGS generates a normally distributed random number.
//
//  Discussion:
//
//    This function generate a normally distributed random number, that is, 
//    it generates random numbers with a Gaussian distribution.  These 
//    random numbers are not exceptionally good, especially in the tails 
//    of the distribution, but this implementation is simple and suitable 
//    for most applications.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Hamming,
//    Numerical Methods for Scientists and Engineers,
//    Dover, 1986,
//    ISBN: 0486652416,
//    LC: QA297.H28.
//
//  Parameters:
//
//    Input, float XMEAN, the mean of the Gaussian distribution.
//
//    Input, float SD, the standard deviation of the Gaussian function.
//
//    Output, float R4_RANDGS, a normally distributed random number.
//
{
  int i;
  float value;

  value = - 6.0E+00;
  for ( i = 1; i <= 12; i++ )
  {
    value = value + r4_rand ( 0.0E+00 );
  }
  value = xmean + sd * value;
  return value;
}
//****************************************************************************80

float r4_random ( float t[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4_RANDOM is a portable pseudorandom number generator.
//
//  Discussion:
//
//    This random number generator is portable amoung a wide variety of
//    computers.  It generates a random number between 0.0 and 1.0 
//    according to the algorithm presented by Bays and Durham.
//
//    The motivation for using this scheme, which resembles the
//    Maclaren-Marsaglia method, is to greatly increase the period of the
//    random sequence.  If the period of the basic generator is P,
//    then the expected mean period of the sequence generated by this
//    generator is given by
//
//      new mean P = sqrt ( pi * factorial ( N ) / ( 8 * P ) ),
//
//    where factorial ( N ) must be much greater than P in this 
//    asymptotic formula.  Generally, N should be 16 to maybe 32.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Carter Bays, Stephen Durham,
//    Improving a Poor Random Number Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 2, Number 1, March 1976, pages 59-64.
//
//  Parameters:
//
//    Input/output, float T[abs(N)+1], an array of random numbers from 
//    a previous invocation of this function.  Whenever N is positive 
//    and differs from the old N, the table is initialized.  The first 
//    iabs(N) numbers are the table discussed in the reference, and the 
//    last value is Y.  This array may be saved in order to restart a sequence.
//
//    Input, int N.  The absolute value of N is the number of random 
//    numbers in an auxiliary table.  Note though that iabs(N)+1 is the 
//    number of items in array T.  If N is positive and differs from its 
//    value in the previous invocation, then the table is initialized for 
//    the new value of N.  If N is negative, abs(N) is the number of items 
//    in an auxiliary table, but the tables are now assumed already to
//    be initialized.  This option enables the user to save the table T at 
//    the end of a long computer run and to restart with the same sequence.  
//    Normally, this function would be called at most once with negative N.  
//    Subsequent invocations would have N positive and of the correct magnitude.
//
//    Output, R4_RANDOM, a random number between 0.0 and 1.0.
//
{
  float dummy;
  static float floatn = -1.0;
  int i;
  int j;
  static int nold = -1;
  float value;

  if ( n != nold )
  {
    nold = abs ( n );
    floatn = ( float ) nold;
 
    if ( n < 0 )
    {
      dummy = r4_rand ( t[nold] );
    }
    else
    {
      for ( i = 1; i <= nold; i++ )
      {
        t[i-1] = r4_rand ( 0.0E+00 );
      }
      t[nold] = r4_rand ( 0.0E+00 );
    }
  }

  j = ( int ) ( t[nold] * floatn + 1.0E+00 );
  t[nold] = t[j-1];
  value = t[j-1];
  t[j-1] = r4_rand ( 0.0E+00 );

  return value;
}
//****************************************************************************80

float r4_ranf ( float sw )

//****************************************************************************80
//
//  Purpose:
//
//    R4_RANF is a driver for R4_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Carter Bays, Stephen Durham,
//    Improving a Poor Random Number Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 2, Number 1, March 1976, pages 59-64.
//
//  Parameters:
//
//    Input, float SW, chooses the action.
//    0.0 <= SW, compute and return the next random number.
//    0.0 > SW, print the internal table, and return the current (old)
//    random number.
//
//    Output, float R4_RANF, the random value.
//
{
  int i;
  static float ranold = 0.0;
  static float t[33];
  float value;

  if ( 0.0E+00 <= sw || ranold == 0.0E+00 )
  {
    value = r4_random ( t, 32 );
    ranold = value;
  }

  if ( sw < 0.0E+00 )
  {
    cout << "\n";
    cout << "  Current random number table:\n";
    cout << "\n:";
    for ( i = 0; i < 33; i++ )
    {
      cout << "  " << i
           << "  " << t[i] << "\n";
    }
    value = ranold;
  }
  return value;
}
//****************************************************************************80

float r4_ren ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_REN is a simple random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Malcolm Pike, David Hill,
//    Algorithm 266:
//    Pseudo-Random Numbers,
//    Communications of the ACM,
//    Volume 8, Number 10, October 1965, page 605.
//
//  Parameters:
//
//    Output, float R4_REN, the random value.
//
{
  static int iy = 100001;
  double value;

  iy = iy * 125;
  iy = iy - ( iy / 2796203 ) * 2796203;
  value = ( float ) ( iy ) / 2796203.0;

  return value;
}
//****************************************************************************80

float r4_shi ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SHI evaluates the hyperbolic sine integral Shi of an R4 argument.
//
//  Discussion:
//
//    Shi ( x ) = Integral ( 0 <= t <= x ) sinh ( t ) dt / t
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_SHI, the hyperbolic sine integral Shi evaluated at X.
//
{
  float absx;
  static int nshi = 0;
  static float shics[7] = {
     0.0078372685688900950695E+00,
     0.0039227664934234563973E+00,
     0.0000041346787887617267E+00,
     0.0000000024707480372883E+00,
     0.0000000000009379295591E+00,
     0.0000000000000002451817E+00,
     0.0000000000000000000467E+00 };
  float value;
  static float xsml = 0.0;

  if ( nshi == 0 )
  {
    nshi = r4_inits ( shics, 7, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( r4_mach ( 3 ) );
  }

  absx = r4_abs ( x );

  if ( absx <= xsml )
  {
    value = x;
  }
  else if ( absx <= 0.375E+00 )
  {
    value = x * ( 1.0E+00 
      + r4_csevl ( 128.0E+00 * x * x / 9.0E+00 - 1.0E+00, 
      shics, nshi ) );
  }
  else
  {
    value = 0.5E+00 * ( r4_ei ( x ) + r4_e1 ( x ) );
  }
  return value;
}
//****************************************************************************80

float r4_si ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SI evaluates the sine integral Si of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_SI, the sine integral Si evaluated at X.
//
{
  float absx;
  float cosx;
  float eps;
  float f;
  float g;
  static int nsi = 0;
  static float pi2 = 1.5707963267948966E+00;
  static float sics[12] = {
    -0.1315646598184841929E+00,
    -0.2776578526973601892E+00,
     0.0354414054866659180E+00,
    -0.0025631631447933978E+00,
     0.0001162365390497009E+00,
    -0.0000035904327241606E+00,
     0.0000000802342123706E+00,
    -0.0000000013562997693E+00,
     0.0000000000179440722E+00,
    -0.0000000000001908387E+00,
     0.0000000000000016670E+00,
    -0.0000000000000000122E+00 };
  float value;
  static float xsml = 0.0;

  if ( nsi == 0 )
  {
    nsi = r4_inits ( sics, 12, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( eps );
  }

  absx = r4_abs ( x );

  if ( absx < xsml )
  {
    value = x;
  }
  else if ( absx <= 4.0E+00 )
  {
    value = x * ( 0.75E+00 
      + r4_csevl ( ( x * x - 8.0E+00 ) * 0.125E+00, sics, nsi ) );
  }
  else
  {
    r4_sifg ( absx, f, g );
    cosx = cos ( absx );

    if ( x < 0.0E+00 )
    {
      value = - pi2 + f * cosx + g * sin ( x );
    }
    else
    {
      value = pi2 - f * cosx - g * sin ( x );
    }
  }

  return value;
}
//****************************************************************************80

void r4_sifg ( float x, float &f, float &g )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIFG is a utility routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float &F, &G.
//
{
  static float f1cs[20] = {
    -0.1191081969051363610E+00,
    -0.0247823144996236248E+00,
     0.0011910281453357821E+00,
    -0.0000927027714388562E+00,
     0.0000093373141568271E+00,
    -0.0000011058287820557E+00,
     0.0000001464772071460E+00,
    -0.0000000210694496288E+00,
     0.0000000032293492367E+00,
    -0.0000000005206529618E+00,
     0.0000000000874878885E+00,
    -0.0000000000152176187E+00,
     0.0000000000027257192E+00,
    -0.0000000000005007053E+00,
     0.0000000000000940241E+00,
    -0.0000000000000180014E+00,
     0.0000000000000035063E+00,
    -0.0000000000000006935E+00,
     0.0000000000000001391E+00,
    -0.0000000000000000282E+00 };
  static float f2cs[29] = {
    -0.0348409253897013234E+00,
    -0.0166842205677959686E+00,
     0.0006752901241237738E+00,
    -0.0000535066622544701E+00,
     0.0000062693421779007E+00,
    -0.0000009526638801991E+00,
     0.0000001745629224251E+00,
    -0.0000000368795403065E+00,
     0.0000000087202677705E+00,
    -0.0000000022601970392E+00,
     0.0000000006324624977E+00,
    -0.0000000001888911889E+00,
     0.0000000000596774674E+00,
    -0.0000000000198044313E+00,
     0.0000000000068641396E+00,
    -0.0000000000024731020E+00,
     0.0000000000009226360E+00,
    -0.0000000000003552364E+00,
     0.0000000000001407606E+00,
    -0.0000000000000572623E+00,
     0.0000000000000238654E+00,
    -0.0000000000000101714E+00,
     0.0000000000000044259E+00,
    -0.0000000000000019634E+00,
     0.0000000000000008868E+00,
    -0.0000000000000004074E+00,
     0.0000000000000001901E+00,
    -0.0000000000000000900E+00,
     0.0000000000000000432E+00 };
  static float g1cs[21] = {
    -0.3040578798253495954E+00,
    -0.0566890984597120588E+00,
     0.0039046158173275644E+00,
    -0.0003746075959202261E+00,
     0.0000435431556559844E+00,
    -0.0000057417294453025E+00,
     0.0000008282552104503E+00,
    -0.0000001278245892595E+00,
     0.0000000207978352949E+00,
    -0.0000000035313205922E+00,
     0.0000000006210824236E+00,
    -0.0000000001125215474E+00,
     0.0000000000209088918E+00,
    -0.0000000000039715832E+00,
     0.0000000000007690431E+00,
    -0.0000000000001514697E+00,
     0.0000000000000302892E+00,
    -0.0000000000000061400E+00,
     0.0000000000000012601E+00,
    -0.0000000000000002615E+00,
     0.0000000000000000548E+00 };
  static float g2cs[34] = {
    -0.0967329367532432218E+00,
    -0.0452077907957459871E+00,
     0.0028190005352706523E+00,
    -0.0002899167740759160E+00,
     0.0000407444664601121E+00,
    -0.0000071056382192354E+00,
     0.0000014534723163019E+00,
    -0.0000003364116512503E+00,
     0.0000000859774367886E+00,
    -0.0000000238437656302E+00,
     0.0000000070831906340E+00,
    -0.0000000022318068154E+00,
     0.0000000007401087359E+00,
    -0.0000000002567171162E+00,
     0.0000000000926707021E+00,
    -0.0000000000346693311E+00,
     0.0000000000133950573E+00,
    -0.0000000000053290754E+00,
     0.0000000000021775312E+00,
    -0.0000000000009118621E+00,
     0.0000000000003905864E+00,
    -0.0000000000001708459E+00,
     0.0000000000000762015E+00,
    -0.0000000000000346151E+00,
     0.0000000000000159996E+00,
    -0.0000000000000075213E+00,
     0.0000000000000035970E+00,
    -0.0000000000000017530E+00,
     0.0000000000000008738E+00,
    -0.0000000000000004487E+00,
     0.0000000000000002397E+00,
    -0.0000000000000001347E+00,
     0.0000000000000000801E+00,
    -0.0000000000000000501E+00 };
  static int nf1 = 0;
  static int nf2 = 0;
  static int ng1 = 0;
  static int ng2 = 0;
  float tol;
  static float xbig = 0.0;
  static float xbnd = 0.0;
  static float xmaxf = 0.0;
  static float xmaxg = 0.0;

  if ( nf1 == 0 )
  {
    tol = 0.1E+00 * r4_mach ( 3 );
    nf1 = r4_inits ( f1cs, 20, tol );
    nf2 = r4_inits ( f2cs, 29, tol );
    ng1 = r4_inits ( g1cs, 21, tol );
    ng2 = r4_inits ( g2cs, 34, tol );
    xbig = sqrt ( 1.0E+00 / r4_mach ( 3 ) );
    xmaxf = exp ( r4_min ( - log ( r4_mach ( 1 ) ), 
      log ( r4_mach ( 2 ) ) ) - 0.01E+00 );
    xmaxg = 1.0E+00 / sqrt ( r4_mach ( 1 ) );
    xbnd = sqrt ( 50.0E+00 );
  }

  if ( x < 4.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_SIFG - Fatal error!\n";
    cerr << "  Approximation invalid for X < 4.\n";
    exit ( 1 );
  }

  if ( x <= xbnd )
  {
    f = ( 1.0E+00 + r4_csevl ( ( 1.0E+00 / x / x - 0.04125E+00 )
      / 0.02125E+00, f1cs, nf1 ) ) / x;
    g = ( 1.0E+00 + r4_csevl ( ( 1.0E+00 / x / x - 0.04125E+00 )
      / 0.02125E+00, g1cs, ng1 ) ) / x / x;
  }
  else if ( x <= xbig )
  {
    f = ( 1.0E+00 + r4_csevl ( 100.0E+00 / x / x - 1.0E+00, 
      f2cs, nf2 ) ) / x;
    g = ( 1.0E+00 + r4_csevl ( 100.0E+00 / x / x - 1.0E+00, 
      g2cs, ng2) ) / x / x;
  }
  else
  {
    if ( x < xmaxf )
    {
      f = 1.0E+00 / x;
    }
    else
    {
      f = 0.0E+00;
    }

    if ( x < xmaxg )
    {
      g = 1.0E+00 / x / x;
    }
    else
    {
      g = 0.0E+00;
    }
  }
  return;
}
//****************************************************************************80

float r4_sign ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIGN returns the sign of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose sign is desired.
//
//    Output, float R4_SIGN, the sign of X.
//
{
  float value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

float r4_sin ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIN evaluates the sine of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_SIN, the sine of X.
//
{
  float f;
  int n2;
  static int ntsn = 0;
  static float pi2rec = 0.636619772367581343E+00;
  static float pihi = 3.140625E+00;
  static float pilo = 9.6765358979323846E-04;
  static float pirec = 0.31830988618379067E+00;
  static float sincs[10] = {
     -0.374991154955873175840E+00,
     -0.181603155237250201864E+00,
     +0.005804709274598633559E+00,
     -0.000086954311779340757E+00,
     +0.000000754370148088851E+00,
     -0.000000004267129665056E+00,
     +0.000000000016980422945E+00,
     -0.000000000000050120579E+00,
     +0.000000000000000114101E+00,
     -0.000000000000000000206E+00 };
  float sgn;
  float value;
  static float xmax = 0.0;
  float xn;
  static float xsml = 0.0;
  static float xwarn = 0.0;
  float y;
//
//  pihi + pilo = pi.  pihi is exactly representable on all machines
//  with at least 8 bits of precision.  whether it is exactly
//  represented depends on the compiler.  this routine is more
//  accurate if it is exactly represented.
//
  if ( ntsn == 0 )
  {
    ntsn = r4_inits ( sincs, 10, 0.1E+00 * r4_mach ( 3 ) );
    xsml = sqrt ( 6.0E+00 * r4_mach ( 3 ) );
    xmax = 1.0E+00 / r4_mach ( 4 );
    xwarn = sqrt ( xmax );
  }

  y = r4_abs ( x );

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R4_SIN - Warning!\n";
    cerr << "  No precision because |X| is big.\n";
    value = 0.0E+00;
    return value;
  }

  if ( xwarn < y )
  {
    cerr << "\n";
    cerr << "R4_SIN - Warning!\n";
    cerr << "  Answer < half precision because |X| is big.\n";
  }

  value = x;

  if ( y < xsml )
  {
    return value;
  }

  xn = r4_aint ( y * pirec + 0.5E+00 );
  n2 = ( int ) ( r4_mod ( xn, 2.0E+00 ) + 0.5E+00 );
  sgn = x;
  if ( n2 != 0 )
  {
    sgn = - sgn;
  }

  f = ( y - xn * pihi ) - xn * pilo;
  xn = 2.0E+00 * ( f * pi2rec ) * ( f * pi2rec ) - 1.0E+00;

  value = f + f * r4_csevl ( xn, sincs, ntsn );

  if ( sgn < 0.0E+00 )
  {
    value = - value;
  }

  if ( value < - 1.0E+00 )
  {
    value = - 1.0E+00;
  }
  else if ( 1.0E+00 < value )
  {
    value = 1.0E+00;
  }
  return value;
}
//****************************************************************************80

float r4_sin_deg ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIN_DEG evaluates the sine of an R4 argument in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument in degrees.
//
//    Output, float R4_SIN_DEG, the sine of X.
//
{
  int n;
  static float raddeg = 0.017453292519943296E+00;
  float value;

  value = sin ( raddeg * x );

  if ( r4_mod ( x, 90.0E+00 ) == 0.0E+00 )
  {
    n = ( int ) ( r4_abs ( x ) / 90.0E+00 + 0.5E+00 );
    n = ( n % 2 );

    if ( n == 0 )
    {
      value = 0.0E+00;
    }
    else if ( value < 0.0E+00 )
    {
      value = - 1.0E+00;
    }
    else
    {
      value = + 1.0E+00;
    }
  }
  return value;
}
//****************************************************************************80

float r4_sinh ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SINH evaluates the hyperbolic sine of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_SINH, the hyperbolic sine of X.
//
{
  int nterms = 0;
  static float sinhcs[8] = {
    0.1730421940471796E+00,
    0.08759422192276048E+00,
    0.00107947777456713E+00,
    0.00000637484926075E+00,
    0.00000002202366404E+00,
    0.00000000004987940E+00,
    0.00000000000007973E+00,
    0.00000000000000009E+00 };
  static float sqeps = 0.0;
  float value;
  float y;
  static float ymax = 0.0;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( sinhcs, 8, 0.1E+00 * r4_mach ( 3 ) );
    sqeps = sqrt ( 6.0E+00 * r4_mach ( 3 ) );
    ymax = 1.0E+00 / sqrt ( r4_mach ( 3 ) );
  }

  y = r4_abs ( x );

  if ( y <= sqeps )
  {
    value = x;
  }
  else if ( y <= 1.0E+00 )
  {
    value = x * ( 1.0E+00 
      + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, sinhcs, nterms ) );
  }
  else
  {
    y = exp ( y );

    if ( ymax <= y )
    {
      value = 0.5E+00 * y;
    }
    else
    {
      value = 0.5E+00 * ( y - 1.0E+00 / y );
    }

    if ( x < 0.0E+00 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

float r4_spence ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SPENCE evaluates a form of Spence's function for an R4 argument.
//
//  Discussion:
//
//    This function evaluates a form of Spence's function defined by
//
//      f(x) = Integral ( 0 <= y <= x ) - log ( abs ( 1 - y ) ) / y dy
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions, page 1004,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//    K Mitchell,
//    Tables of the function Integral ( 0 < y < x ) - log | 1 - y | dy / y
//    with an account of some properties of this and related functions,
//    Philosophical Magazine,
//    Volume 40, pages 351-368, 1949.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_SPENCE, Spence's function evaluated at X.
//
{
  float aln;
  static int nspenc = 0;
  static float pi26 = 1.644934066848226E+00;
  static float spencs[19] = {
    0.1527365598892406E+00,
    0.08169658058051014E+00,
    0.00581415714077873E+00,
    0.00053716198145415E+00,
    0.00005724704675185E+00,
    0.00000667454612164E+00,
    0.00000082764673397E+00,
    0.00000010733156730E+00,
    0.00000001440077294E+00,
    0.00000000198444202E+00,
    0.00000000027940058E+00,
    0.00000000004003991E+00,
    0.00000000000582346E+00,
    0.00000000000085767E+00,
    0.00000000000012768E+00,
    0.00000000000001918E+00,
    0.00000000000000290E+00,
    0.00000000000000044E+00,
    0.00000000000000006E+00 };
  float value;
  static float xbig = 0.0;

  if ( nspenc == 0 )
  {
    nspenc = r4_inits ( spencs, 19, 0.1E+00 * r4_mach ( 3 ) );
    xbig = 1.0E+00 / r4_mach ( 3 );
  }

  if ( x <= - xbig )
  {
    aln = log ( 1.0E+00 - x );
    value = - pi26 - 0.5E+00 * aln * ( 2.0E+00 * log ( - x ) - aln );
  }
  else if ( x <= - 1.0E+00 )
  {
    aln = log ( 1.0E+00 - x );

    value = - pi26 - 0.5E+00 * aln * ( 2.0E+00 
      * log ( - x ) - aln ) + ( 1.0E+00 + r4_csevl (
      4.0E+00 / ( 1.0E+00 - x ) - 1.0E+00, spencs, nspenc ) ) 
      / ( 1.0E+00 - x );
  }
  else if ( x <= 0.0E+00 )
  {
    value = - 0.5E+00 * log ( 1.0E+00 - x ) 
      * log ( 1.0E+00 - x ) - x * ( 1.0E+00 + r4_csevl ( 
      4.0E+00 * x / ( x - 1.0E+00 ) - 1.0E+00, spencs, nspenc ) ) 
      / ( x - 1.0E+00 );
  }
  else if ( x <= 0.5E+00 )
  {
    value = x * ( 1.0E+00 + r4_csevl ( 4.0E+00 * x - 1.0E+00, 
      spencs, nspenc ) );
  }
  else if ( x < 1.0E+00 )
  {
    value = pi26 - log ( x ) * log ( 1.0E+00 - x )
      - ( 1.0E+00 - x ) * ( 1.0E+00 + r4_csevl ( 4.0E+00 
      * ( 1.0E+00 - x ) - 1.0E+00, spencs, nspenc ) );
  }
  else if ( x == 1.0E+00 )
  {
    value = pi26;
  }
  else if ( x <= 2.0E+00 )
  {
    value = pi26 - 0.5E+00 * log ( x ) 
      * log ( ( x - 1.0E+00 ) * ( x - 1.0E+00 ) / x )
      + ( x - 1.0E+00 ) * ( 1.0E+00 + r4_csevl ( 4.0E+00 
      * ( x - 1.0E+00 ) / x - 1.0E+00, spencs, nspenc ) ) / x;
  }
  else if ( x < xbig )
  {
    value = 2.0E+00 * pi26 - 0.5E+00 * log ( x ) * log ( x )
      - ( 1.0E+00 + r4_csevl ( 4.0E+00 / x - 1.0E+00, spencs, 
      nspenc ) ) / x;
  }
  else
  {
    value = 2.0E+00 * pi26 - 0.5E+00 * log ( x ) * log ( x );
  }
  return value;
}
//****************************************************************************80

float r4_sqrt ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SQRT computes the square root of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the number whose square root is desired.
//
//    Output, float R4_SQRT, the square root of X.
//
{
  int irem;
  int iter;
  int ixpnt;
  int n;
  static int niter = 0;
  static float sqrt2[3] = {
     0.70710678118654752E+00,
     1.0E+00,
     1.41421356237309505E+00 };
  float value;
  float y;

  if ( niter == 0 )
  {
    niter = 1.443E+00 * r4_log ( - 0.104E+00 
      * r4_log ( 0.1E+00 * r4_mach ( 3 ) ) ) + 1.0E+00;
  }

  if ( x < 0.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_SQRT - Fatal error!\n";
    cerr << "  X < 0.0\n";
    exit ( 1 );
  }
  else if ( x == 0.0E+00 )
  {
    value = 0.0E+00;
  }
  else
  {
    r4_upak ( x, y, n );
    ixpnt = n / 2;
    irem = n - 2 * ixpnt + 2;
    value = 0.261599E+00 + y * ( 1.114292E+00 
      + y * ( -0.516888E+00 + y * 0.141067E+00 ) );

    for ( iter = 1; iter <= niter; iter++ )
    {
      value = value + 0.5E+00 * ( y - value * value ) / value;
    }
    value = r4_pak ( sqrt2[irem-1] * value, ixpnt );
  }
  return value;
}
//****************************************************************************80

float r4_tan ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_TAN evaluates the tangent of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_TAN, the tangent of X.
//
{
  float ainty;
  float ainty2;
  int ifn;
  static int nterms = 0;
  static float pi2rec = 0.0116197723675813430E+00;
  float prodbg;
  static float sqeps = 0.0;
  static float tancs[11] = {
     0.226279327631293578E+00,
     0.0430179131465489618E+00,
     0.0006854461068256508E+00,
     0.0000110453269475970E+00,
     0.0000001781747790392E+00,
     0.0000000028744968582E+00,
     0.0000000000463748541E+00,
     0.0000000000007481760E+00,
     0.0000000000000120704E+00,
     0.0000000000000001947E+00,
     0.0000000000000000031E+00 };
  float value;
  static float xmax = 0.0;
  static float xsml = 0.0;
  float y;
  float yrem;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( tancs, 11, 0.1E+00 * r4_mach ( 3 ) );
    xmax = 1.0E+00 / r4_mach ( 4 );
    xsml = sqrt ( 3.0E+00 * r4_mach ( 3 ) );
    sqeps = sqrt ( r4_mach ( 4 ) );
  }

  y = r4_abs ( x );

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R4_TAN - Warning\n";
    cerr << "  No precision because |X| is big.\n";
    value = 0.0E+00;
    return value;
  }
//
//  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
//  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
//  = aint(.625*y) + aint(z) + rem(z)
//
  ainty = r4_aint ( y );
  yrem = y - ainty;
  prodbg = 0.625E+00 * ainty;
  ainty = r4_aint ( prodbg );
  y = ( prodbg - ainty ) + 0.625E+00 * yrem + y * pi2rec;
  ainty2 = r4_aint ( y );
  ainty = ainty + ainty2;
  y = y - ainty2;

  ifn = ( int ) r4_mod ( ainty, 2.0E+00 );
  if ( ifn == 1 )
  {
    y = 1.0E+00 - y;
  }

  if ( 1.0E+00 - y < r4_abs ( x ) * sqeps )
  {
    cerr << "\n";
    cerr << "R4_TAN - Warning!\n";
    cerr << "  Answer < half precision.\n";
    cerr << "  |X| big or X near pi/2 or 3*pi/2.\n";
  }

  if ( y == 1.0E+00 )
  {
    cerr << "\n";
    cerr << "R4_TAN - Fatal error!\n";
    cerr << "  X is pi/2 or 3*pi/2.\n";
    exit ( 1 );
  }

  if ( y <= 0.25E+00 )
  {
    value = y;
    if ( xsml < y )
    {
      value = y * ( 1.5E+00 
        + r4_csevl ( 32.0E+00 * y * y - 1.0E+00, tancs, nterms ) );
    }
  }
  else if ( y <= 0.5E+00 )
  {
    value = 0.5E+00 * y * ( 1.5E+00 
      + r4_csevl ( 8.0E+00 * y * y - 1.0E+00, tancs, nterms ) );

    value = 2.0E+00 * value / ( 1.0E+00 - value * value );
  }
  else
  {
    value = 0.25E+00 * y * ( 1.5E+00 
      + r4_csevl ( 2.0E+00 * y * y - 1.0E+00, tancs, nterms ) );
    value = 2.0E+00 * value / ( 1.0E+00 - value * value );
    value = 2.0E+00 * value / ( 1.0E+00 - value * value );
  }

  if ( x < 0.0E+00 )
  {
    value = - r4_abs ( value );
 }
  else if ( 0.0E+00 < x )
  {
    value = + r4_abs ( value );
  }

  if ( ifn == 1 )
  {
    value = - value;
  }

  return value;
}
//****************************************************************************80

float r4_tanh ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_TANH evaluates the hyperbolic tangent of an R4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, float X, the argument.
//
//    Output, float R4_TANH, the hyperbolic tangent of X.
//
{
  static int nterms = 0;
  static float sqeps = 0.0;
  static float tanhcs[17] = {
   -0.25828756643634710E+00,
   -0.11836106330053497E+00,
   +0.009869442648006398E+00,
   -0.000835798662344582E+00,
   +0.000070904321198943E+00,
   -0.000006016424318120E+00,
   +0.000000510524190800E+00,
   -0.000000043320729077E+00,
   +0.000000003675999055E+00,
   -0.000000000311928496E+00,
   +0.000000000026468828E+00,
   -0.000000000002246023E+00,
   +0.000000000000190587E+00,
   -0.000000000000016172E+00,
   +0.000000000000001372E+00,
   -0.000000000000000116E+00,
   +0.000000000000000009E+00 };
  float value;
  static float xmax = 0.0;
  float y;

  if ( nterms == 0 )
  {
    nterms = r4_inits ( tanhcs, 17, 0.1E+00 * r4_mach ( 3 ) );
    sqeps = sqrt ( 3.0E+00 * r4_mach ( 3 ) );
    xmax = - 0.5E+00 * log ( r4_mach ( 3 ) );
  }

  y = r4_abs ( x );

  if ( y <= sqeps )
  {
    value = x;
  }
  else if ( y <= 1.0E+00 )
  {
    value = x * ( 1.0E+00 
      + r4_csevl ( 2.0E+00 * x * x - 1.0E+00, tanhcs, nterms ) );
  }
  else if ( y <= xmax )
  {
    y = exp ( y );
    value = ( y - 1.0E+00 / y ) / ( y + 1.0E+00 / y );
    if ( x < 0.0E+00 )
    {
      value = - value;
    }
  }
  else
  {
    if ( x < 0.0E+00 )
    {
      value = - 1.0E+00;
    }
    else
    {
      value = + 1.0E+00;
    }
  }
  return value;
}
//****************************************************************************80

void r4_upak ( float x, float &y, int &n )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UPAK unpacks an R4 into a mantissa and exponent.
//
//  Discussion:
//
//    This function unpacks a floating point number x so that
//
//      x = y * 2.0^n
//
//    where
//
//      0.5 <= abs ( y ) < 1.0 .
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, float X, the number to be unpacked.
//
//    Output, float &Y, the mantissa.
//
//    Output, int &N, the exponent.
//
{
  float absx;

  absx = r4_abs ( x );
  n = 0;
  y = 0.0;

  if ( x == 0.0 )
  {
    return;
  }

  while ( absx < 0.5 )
  {
    n = n - 1;
    absx = absx * 2.0;
  }

  while ( 1.0 <= absx )
  {
    n = n + 1;
    absx = absx * 0.5;
  }

  if ( x < 0.0 )
  {
    y = - absx;
  }
  else
  {
    y = + absx;
  }

  return;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_acos ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ACOS evaluates the arc-cosine of an R8 argument.
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
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ACOS, the arc-cosine of X.
//
{
  static double pi2 = 1.57079632679489661923132169163975;
  double value;

  value = pi2 - r8_asin ( x );

  return value;
}
//****************************************************************************80

double r8_acosh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ACOSH evaluates the arc-hyperbolic cosine of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ACOSH, the arc-hyperbolic cosine of X.
//
{
  static double dln2 = 0.69314718055994530941723212145818;
  double value;
  static double xmax = 0.0;

  if ( xmax == 0.0 )
  {
    xmax = 1.0 / sqrt ( r8_mach ( 3 ) );
  }

  if ( x < 1.0 )
  {
    cerr << "\n";
    cerr << "R8_ACOSH - Fatal error!\n";
    cerr << "  X < 1.0\n";
    exit ( 1 );
  }
  else if ( x < xmax )
  {
    value = log ( x + sqrt ( x * x - 1.0 ) );
  }
  else
  {
    value = dln2 + log ( x );
  }
  return value;
}
//****************************************************************************80

void r8_admp ( double x, double &ampl, double &phi )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ADMP: modulus and phase of the derivative of the Airy function.
//
//  Description:
//
//    This function must only be called when X <= -1.0.
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
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double &AMPL, &PHI, the modulus and phase of the 
//    derivative of the Airy function.
//
{
  static double an20cs[57] = {
      0.0126732217145738027154610751034240,
     -0.0005212847072615621184780942309478,
     -0.0000052672111140370429809074052969,
     -0.0000001628202185026483752632460680,
     -0.0000000090991442687371386325973075,
     -0.0000000007438647126242192890685403,
     -0.0000000000795494751591469486122822,
     -0.0000000000104050944288303742803960,
     -0.0000000000015932425598414551523990,
     -0.0000000000002770648272341913946674,
     -0.0000000000000535342629237606295104,
     -0.0000000000000113061541781728314051,
     -0.0000000000000025772190078943167788,
     -0.0000000000000006278033116032485076,
     -0.0000000000000001621295400189939757,
     -0.0000000000000000440992985240675353,
     -0.0000000000000000125655516553258972,
     -0.0000000000000000037336906988015204,
     -0.0000000000000000011524626926724671,
     -0.0000000000000000003683081499099144,
     -0.0000000000000000001215206965331797,
     -0.0000000000000000000412916177724016,
     -0.0000000000000000000144177364239347,
     -0.0000000000000000000051631842875864,
     -0.0000000000000000000018931242668250,
     -0.0000000000000000000007096054668569,
     -0.0000000000000000000002715406646904,
     -0.0000000000000000000001059486979400,
     -0.0000000000000000000000421030035685,
     -0.0000000000000000000000170233781664,
     -0.0000000000000000000000069966677028,
     -0.0000000000000000000000029206643813,
     -0.0000000000000000000000012373128203,
     -0.0000000000000000000000005315871095,
     -0.0000000000000000000000002314622618,
     -0.0000000000000000000000001020779922,
     -0.0000000000000000000000000455706227,
     -0.0000000000000000000000000205831071,
     -0.0000000000000000000000000094015189,
     -0.0000000000000000000000000043405874,
     -0.0000000000000000000000000020247792,
     -0.0000000000000000000000000009539214,
     -0.0000000000000000000000000004537234,
     -0.0000000000000000000000000002178016,
     -0.0000000000000000000000000001054823,
     -0.0000000000000000000000000000515242,
     -0.0000000000000000000000000000253763,
     -0.0000000000000000000000000000125983,
     -0.0000000000000000000000000000063030,
     -0.0000000000000000000000000000031771,
     -0.0000000000000000000000000000016131,
     -0.0000000000000000000000000000008248,
     -0.0000000000000000000000000000004246,
     -0.0000000000000000000000000000002200,
     -0.0000000000000000000000000000001147,
     -0.0000000000000000000000000000000602,
     -0.0000000000000000000000000000000318 };
  static double an21cs[60] = {
      0.0198313155263169394420342483165643,
     -0.0029376249067087533460593745594484,
     -0.0001136260695958195549872611137182,
     -0.0000100554451087156009750981645918,
     -0.0000013048787116563250421785598252,
     -0.0000002123881993150664830666079609,
     -0.0000000402270833384269040347850109,
     -0.0000000084996745953161799142201792,
     -0.0000000019514839426178614099532934,
     -0.0000000004783865343840384282992480,
     -0.0000000001236733992099450501137105,
     -0.0000000000334137486398754232219789,
     -0.0000000000093702823540766329897780,
     -0.0000000000027130128156139564687240,
     -0.0000000000008075953800583479535949,
     -0.0000000000002463214304700125252160,
     -0.0000000000000767655689109321564410,
     -0.0000000000000243882598807354919791,
     -0.0000000000000078831466358760308462,
     -0.0000000000000025882400995585864077,
     -0.0000000000000008619457862945690828,
     -0.0000000000000002907994739663128534,
     -0.0000000000000000992846796122890484,
     -0.0000000000000000342720229187774480,
     -0.0000000000000000119511048205515026,
     -0.0000000000000000042069729043678359,
     -0.0000000000000000014939697762818400,
     -0.0000000000000000005348981161589517,
     -0.0000000000000000001929877577826238,
     -0.0000000000000000000701313701018203,
     -0.0000000000000000000256585738509682,
     -0.0000000000000000000094475894562734,
     -0.0000000000000000000034996401941465,
     -0.0000000000000000000013037622466397,
     -0.0000000000000000000004883334163346,
     -0.0000000000000000000001838477586152,
     -0.0000000000000000000000695527324058,
     -0.0000000000000000000000264351910209,
     -0.0000000000000000000000100918094655,
     -0.0000000000000000000000038688924289,
     -0.0000000000000000000000014892036525,
     -0.0000000000000000000000005754342426,
     -0.0000000000000000000000002231725971,
     -0.0000000000000000000000000868607480,
     -0.0000000000000000000000000339220403,
     -0.0000000000000000000000000132910128,
     -0.0000000000000000000000000052239309,
     -0.0000000000000000000000000020594383,
     -0.0000000000000000000000000008142614,
     -0.0000000000000000000000000003228473,
     -0.0000000000000000000000000001283529,
     -0.0000000000000000000000000000511622,
     -0.0000000000000000000000000000204451,
     -0.0000000000000000000000000000081901,
     -0.0000000000000000000000000000032886,
     -0.0000000000000000000000000000013235,
     -0.0000000000000000000000000000005338,
     -0.0000000000000000000000000000002158,
     -0.0000000000000000000000000000000874,
     -0.0000000000000000000000000000000355 };
  static double an22cs[74] = {
      0.0537418629629794329091103360917783,
     -0.0126661435859883193466312085036450,
     -0.0011924334106593006840848916913681,
     -0.0002032327627275654552687155176363,
     -0.0000446468963075163979516164905945,
     -0.0000113359036053123490416997893086,
     -0.0000031641352378546107356671355827,
     -0.0000009446708886148939120888532442,
     -0.0000002966562236471765527900905456,
     -0.0000000969118892024367799908661433,
     -0.0000000326822538653274091533072559,
     -0.0000000113144618963583865900447294,
     -0.0000000040042691001741501738278050,
     -0.0000000014440333683907423778522199,
     -0.0000000005292853746152611585663541,
     -0.0000000001967763373707889528245726,
     -0.0000000000740800095755849858816731,
     -0.0000000000282016314294661982842740,
     -0.0000000000108440066463128331337590,
     -0.0000000000042074800682644236920617,
     -0.0000000000016459149670634819724739,
     -0.0000000000006486826705121018896077,
     -0.0000000000002574095003354105832300,
     -0.0000000000001027889029407822132143,
     -0.0000000000000412845827195222720128,
     -0.0000000000000166711029332862509726,
     -0.0000000000000067656696165608023403,
     -0.0000000000000027585448232693576823,
     -0.0000000000000011296397915297168938,
     -0.0000000000000004644848225457314333,
     -0.0000000000000001917198035033912928,
     -0.0000000000000000794197570111893530,
     -0.0000000000000000330116492300368930,
     -0.0000000000000000137658057726549714,
     -0.0000000000000000057578093720012791,
     -0.0000000000000000024152700858632017,
     -0.0000000000000000010159301700933666,
     -0.0000000000000000004284434955330055,
     -0.0000000000000000001811344052168016,
     -0.0000000000000000000767602045619422,
     -0.0000000000000000000326026346758614,
     -0.0000000000000000000138773806682627,
     -0.0000000000000000000059191627103729,
     -0.0000000000000000000025297256431944,
     -0.0000000000000000000010832077293819,
     -0.0000000000000000000004646674880404,
     -0.0000000000000000000001996797783865,
     -0.0000000000000000000000859524108705,
     -0.0000000000000000000000370584152073,
     -0.0000000000000000000000160027503479,
     -0.0000000000000000000000069208124999,
     -0.0000000000000000000000029974448994,
     -0.0000000000000000000000013000356362,
     -0.0000000000000000000000005646100942,
     -0.0000000000000000000000002455341103,
     -0.0000000000000000000000001069119686,
     -0.0000000000000000000000000466095090,
     -0.0000000000000000000000000203441579,
     -0.0000000000000000000000000088900866,
     -0.0000000000000000000000000038891813,
     -0.0000000000000000000000000017032637,
     -0.0000000000000000000000000007467295,
     -0.0000000000000000000000000003277097,
     -0.0000000000000000000000000001439618,
     -0.0000000000000000000000000000633031,
     -0.0000000000000000000000000000278620,
     -0.0000000000000000000000000000122743,
     -0.0000000000000000000000000000054121,
     -0.0000000000000000000000000000023884,
     -0.0000000000000000000000000000010549,
     -0.0000000000000000000000000000004663,
     -0.0000000000000000000000000000002063,
     -0.0000000000000000000000000000000913,
     -0.0000000000000000000000000000000405 };
  static double aph0cs[53] = {
     -0.0855849241130933256920124260179491,
      0.0011214378867065260735786722471124,
      0.0000042721029353664113951573742015,
      0.0000000817607381483243644018062323,
      0.0000000033907645000492724207816418,
      0.0000000002253264422619113939845276,
      0.0000000000206284209229015251256990,
      0.0000000000023858762828130887627258,
      0.0000000000003301618105886705480628,
      0.0000000000000527009648508328581123,
      0.0000000000000094555482203813492868,
      0.0000000000000018709426951344836908,
      0.0000000000000004023980041825392741,
      0.0000000000000000930192879258983167,
      0.0000000000000000229038635402379945,
      0.0000000000000000059634359822083386,
      0.0000000000000000016320279659403399,
      0.0000000000000000004671145658861339,
      0.0000000000000000001392334415363502,
      0.0000000000000000000430642670285155,
      0.0000000000000000000137781416318755,
      0.0000000000000000000045476710480396,
      0.0000000000000000000015448420203026,
      0.0000000000000000000005389770551212,
      0.0000000000000000000001927726737155,
      0.0000000000000000000000705659320166,
      0.0000000000000000000000263985084827,
      0.0000000000000000000000100791301805,
      0.0000000000000000000000039228928481,
      0.0000000000000000000000015547422955,
      0.0000000000000000000000006268306372,
      0.0000000000000000000000002568563962,
      0.0000000000000000000000001068858883,
      0.0000000000000000000000000451347253,
      0.0000000000000000000000000193267262,
      0.0000000000000000000000000083865369,
      0.0000000000000000000000000036857386,
      0.0000000000000000000000000016396202,
      0.0000000000000000000000000007379298,
      0.0000000000000000000000000003358392,
      0.0000000000000000000000000001544891,
      0.0000000000000000000000000000718013,
      0.0000000000000000000000000000337026,
      0.0000000000000000000000000000159710,
      0.0000000000000000000000000000076382,
      0.0000000000000000000000000000036855,
      0.0000000000000000000000000000017935,
      0.0000000000000000000000000000008800,
      0.0000000000000000000000000000004353,
      0.0000000000000000000000000000002170,
      0.0000000000000000000000000000001090,
      0.0000000000000000000000000000000551,
      0.0000000000000000000000000000000281 };
  static double aph1cs[58] = {
     -0.1024172908077571694021123321813917,
      0.0071697275146591248047211649144704,
      0.0001209959363122328589813856491397,
      0.0000073361512841219912080297845684,
      0.0000007535382954271607069982903869,
      0.0000001041478171741301926885109155,
      0.0000000174358728518545691858907606,
      0.0000000033399795033346451660184961,
      0.0000000007073075174363527083399508,
      0.0000000001619187515189773266792272,
      0.0000000000394539981881954889879668,
      0.0000000000101192281734227133292631,
      0.0000000000027092778259520332198030,
      0.0000000000007523806418422548885854,
      0.0000000000002156368733008966357328,
      0.0000000000000635282777126068410174,
      0.0000000000000191756972641501729345,
      0.0000000000000059143072446464891558,
      0.0000000000000018597128517275028357,
      0.0000000000000005950444923946103668,
      0.0000000000000001934229956430180252,
      0.0000000000000000637843021489504324,
      0.0000000000000000213127290087312393,
      0.0000000000000000072081380656728500,
      0.0000000000000000024652494144769247,
      0.0000000000000000008519110570266154,
      0.0000000000000000002972384468491170,
      0.0000000000000000001046426648811446,
      0.0000000000000000000371493036347327,
      0.0000000000000000000132923247793472,
      0.0000000000000000000047912837925909,
      0.0000000000000000000017390619859336,
      0.0000000000000000000006353585173501,
      0.0000000000000000000002335643614263,
      0.0000000000000000000000863643881606,
      0.0000000000000000000000321123006944,
      0.0000000000000000000000120031540983,
      0.0000000000000000000000045091488699,
      0.0000000000000000000000017020228580,
      0.0000000000000000000000006453744630,
      0.0000000000000000000000002457788564,
      0.0000000000000000000000000939897684,
      0.0000000000000000000000000360863150,
      0.0000000000000000000000000139077884,
      0.0000000000000000000000000053797184,
      0.0000000000000000000000000020882551,
      0.0000000000000000000000000008133371,
      0.0000000000000000000000000003178080,
      0.0000000000000000000000000001245700,
      0.0000000000000000000000000000489742,
      0.0000000000000000000000000000193099,
      0.0000000000000000000000000000076349,
      0.0000000000000000000000000000030269,
      0.0000000000000000000000000000012032,
      0.0000000000000000000000000000004795,
      0.0000000000000000000000000000001915,
      0.0000000000000000000000000000000767,
      0.0000000000000000000000000000000308 };
  static double aph2cs[72] = {
     -0.2057088719781465106973648665602125,
      0.0422196961357771921673114980369460,
      0.0020482560511207275042660577813334,
      0.0002607800735165005631187879922652,
      0.0000474824268004728875381750519293,
      0.0000105102756431611743473630026955,
      0.0000026353534014667945109314041983,
      0.0000007208824863499147299790783731,
      0.0000002103236664473352859749477082,
      0.0000000644975634555295598437362273,
      0.0000000205802377264368507978116888,
      0.0000000067836273920906428963513918,
      0.0000000022974015284009400168343792,
      0.0000000007961306765491187534883226,
      0.0000000002813860609741591719003632,
      0.0000000001011749056931973922841793,
      0.0000000000369306737952476559097060,
      0.0000000000136615066127098031778842,
      0.0000000000051142751416045045119388,
      0.0000000000019351688931706516247975,
      0.0000000000007393606916493224217271,
      0.0000000000002849792219222743597555,
      0.0000000000001107280782459648335733,
      0.0000000000000433412199370134633169,
      0.0000000000000170800825265670367471,
      0.0000000000000067733080195631114673,
      0.0000000000000027016904789262414108,
      0.0000000000000010834720751810782141,
      0.0000000000000004367060312970286167,
      0.0000000000000001768511738053366608,
      0.0000000000000000719359213093645717,
      0.0000000000000000293823610002933154,
      0.0000000000000000120482811525848357,
      0.0000000000000000049586659491091389,
      0.0000000000000000020479438315847217,
      0.0000000000000000008486019944410629,
      0.0000000000000000003527351765384506,
      0.0000000000000000001470563996804903,
      0.0000000000000000000614817826902188,
      0.0000000000000000000257737706565077,
      0.0000000000000000000108323903590042,
      0.0000000000000000000045638898024998,
      0.0000000000000000000019273635403662,
      0.0000000000000000000008157668569775,
      0.0000000000000000000003460202828346,
      0.0000000000000000000001470726482427,
      0.0000000000000000000000626356074088,
      0.0000000000000000000000267261292780,
      0.0000000000000000000000114246948763,
      0.0000000000000000000000048923460516,
      0.0000000000000000000000020985807810,
      0.0000000000000000000000009016618807,
      0.0000000000000000000000003880129464,
      0.0000000000000000000000001672282170,
      0.0000000000000000000000000721790800,
      0.0000000000000000000000000311982573,
      0.0000000000000000000000000135035015,
      0.0000000000000000000000000058524861,
      0.0000000000000000000000000025397686,
      0.0000000000000000000000000011035457,
      0.0000000000000000000000000004800788,
      0.0000000000000000000000000002090956,
      0.0000000000000000000000000000911743,
      0.0000000000000000000000000000397998,
      0.0000000000000000000000000000173923,
      0.0000000000000000000000000000076083,
      0.0000000000000000000000000000033316,
      0.0000000000000000000000000000014604,
      0.0000000000000000000000000000006407,
      0.0000000000000000000000000000002814,
      0.0000000000000000000000000000001237,
      0.0000000000000000000000000000000544 };
  double eta;
  static int nan20 = 0;
  static int nan21 = 0;
  static int nan22 = 0;
  static int naph0 = 0;
  static int naph1 = 0;
  static int naph2 = 0;
  static double pi34 = 2.35619449019234492884698253745962716313;
  double sqrtx;
  static double xsml = 0.0;
  double z;

  if ( nan20 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nan20 = r8_inits ( an20cs, 57, eta );
    nan21 = r8_inits ( an21cs, 60, eta );
    nan22 = r8_inits ( an22cs, 74, eta );
    naph0 = r8_inits ( aph0cs, 53, eta );
    naph1 = r8_inits ( aph1cs, 58, eta );
    naph2 = r8_inits ( aph2cs, 72, eta );
    xsml = - r8_power ( 128.0 / r8_mach ( 3 ), 0.3333 );
  }

  if ( x < xsml )
  {
    z = 1.0;
    ampl = 0.3125 + r8_csevl ( z, an20cs, nan20 );
    phi = - 0.625 + r8_csevl ( z, aph0cs, naph0 );
  }
  else if ( x < - 4.0 )
  {
    z = 128.0 / x / x / x + 1.0;
    ampl = 0.3125 + r8_csevl ( z, an20cs, nan20 );
    phi = - 0.625 + r8_csevl ( z, aph0cs, naph0 );
  }
  else if ( x < - 2.0 )
  {
    z = ( 128. / x / x / x + 9.0 ) / 7.0;
    ampl = 0.3125 + r8_csevl ( z, an21cs, nan21 );
    phi = - 0.625 + r8_csevl ( z, aph1cs, naph1 );
  }
  else if ( x <= - 1.0 )
  {
    z = ( 16.0 / x / x / x + 9.0 ) / 7.0;
    ampl = 0.3125 + r8_csevl ( z, an22cs, nan22 );
    phi = - 0.625 + r8_csevl ( z, aph2cs, naph2 );
  }
  else
  {
    cerr << "\n";
    cerr << "R8_ADMP - Fatal error!\n";
    cerr << "  - 1.0 < X.\n";
    exit ( 1 );
  }

  sqrtx = sqrt ( - x );
  ampl = sqrt ( ampl * sqrtx );
  phi = pi34 - x * sqrtx * phi;

  return;
}
//****************************************************************************80

double r8_ai ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_AI evaluates the Airy function Ai of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_AI, the Airy function Ai of X.
//
{
  static double aifcs[13] = {
    -0.37971358496669997496197089469414E-01,
    +0.59191888537263638574319728013777E-01,
    +0.98629280577279975365603891044060E-03,
    +0.68488438190765667554854830182412E-05,
    +0.25942025962194713019489279081403E-07,
    +0.61766127740813750329445749697236E-10,
    +0.10092454172466117901429556224601E-12,
    +0.12014792511179938141288033225333E-15,
    +0.10882945588716991878525295466666E-18,
    +0.77513772196684887039238400000000E-22,
    +0.44548112037175638391466666666666E-25,
    +0.21092845231692343466666666666666E-28,
    +0.83701735910741333333333333333333E-32 };
  static double aigcs[13] = {
    +0.18152365581161273011556209957864E-01,
    +0.21572563166010755534030638819968E-01,
    +0.25678356987483249659052428090133E-03,
    +0.14265214119792403898829496921721E-05,
    +0.45721149200180426070434097558191E-08,
    +0.95251708435647098607392278840592E-11,
    +0.13925634605771399051150420686190E-13,
    +0.15070999142762379592306991138666E-16,
    +0.12559148312567778822703205333333E-19,
    +0.83063073770821340343829333333333E-23,
    +0.44657538493718567445333333333333E-26,
    +0.19900855034518869333333333333333E-29,
    +0.74702885256533333333333333333333E-33 };
  static int naif = 0;
  static int naig = 0;
  double theta;
  double value;
  static double x3sml = 0.0;
  double xm;
  static double xmax = 0.0;
  double z;

  if ( naif == 0 )
  {
    naif = r8_inits ( aifcs, 13, 0.1 * r8_mach ( 3 ) );
    naig = r8_inits ( aigcs, 13, 0.1 * r8_mach ( 3 ) );
    x3sml = r8_power ( r8_mach ( 3 ), 0.3334 );
    xmax = r8_power ( - 1.5 * log ( r8_mach ( 1 ) ), 0.6667 );
    xmax = xmax - xmax * log ( xmax ) / 
      ( 4.0 * xmax * sqrt ( xmax ) + 1.0 ) - 0.01;
  }

  if ( x < - 1.0 )
  {
    r8_aimp ( x, xm, theta );
    value = xm * cos ( theta );
  }
  else if ( r8_abs ( x ) <= x3sml )
  {
    z = 0.0;
    value = 0.375 + ( r8_csevl ( z, aifcs, naif ) 
      - x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
  }
  else if ( x <= 1.0 )
  {
    z = x * x * x;
    value = 0.375 + ( r8_csevl ( z, aifcs, naif ) 
      - x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
  }
  else if ( x <= xmax )
  {
    value = r8_aie ( x ) * exp ( - 2.0 * x * sqrt ( x ) / 3.0 );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double r8_aid ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_AID evaluates the derivative of the Airy function Ai of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_AID, the derivative of the Airy function 
//    Ai of X.
//
{
  static double aifcs[13] = {
     0.105274612265314088088970057325134114,
     0.011836136281529978442889292583980840,
     0.000123281041732256643051689242469164,
     0.000000622612256381399016825658693579,
     0.000000001852988878441452950548140821,
     0.000000000003633288725904357915995625,
     0.000000000000005046217040440664768330,
     0.000000000000000005223816555471480985,
     0.000000000000000000004185745090748989,
     0.000000000000000000000002672887324883,
     0.000000000000000000000000001392128006,
     0.000000000000000000000000000000602653,
     0.000000000000000000000000000000000220 };
  static double aigcs[13] = {
     0.0212338781509186668523122276848937,
     0.0863159303352144067524942809461604,
     0.0017975947203832313578033963225230,
     0.0000142654998755506932526620687495,
     0.0000000594379952836832010488787064,
     0.0000000001524033664794478945214786,
     0.0000000000002645876603490435305100,
     0.0000000000000003315624296815020591,
     0.0000000000000000003139789757594792,
     0.0000000000000000000002325767379040,
     0.0000000000000000000000001384384231,
     0.0000000000000000000000000000676629,
     0.0000000000000000000000000000000276 };
  double eta;
  static int naif = 0;
  static int naig = 0;
  double phi;
  double value;
  double x2;
  static double x2sml = 0.0;
  double x3;
  static double x3sml = 0.0;
  double xn;

  if ( naif == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    naif = r8_inits ( aifcs, 13, eta );
    naig = r8_inits ( aigcs, 13, eta );
    x3sml = r8_power ( r8_mach ( 3 ), 0.3334 );
    x2sml = sqrt ( r8_mach ( 3 ) );
  }

  if ( x < - 1.0 )
  {
    r8_admp ( x, xn, phi );
    value = xn * cos ( phi );
  }
  else if ( r8_abs ( x ) <= x2sml )
  {
    x2 = 0.0;
    x3 = 0.0;
    value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) ) 
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25;
  }
  else if ( r8_abs ( x ) <= x3sml )
  {
    x2 = x * x;
    x3 = 0.0;
    value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) ) 
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25;
  }
  else if ( x <= 1.0 )
  {
    x2 = x * x;
    x3 = x * x * x;
    value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) ) 
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25;
  }
  else
  {
    value = r8_aide ( x ) * exp ( - 2.0 * x * sqrt ( x ) / 3.0 );
  }
  return value;
}
//****************************************************************************80

double r8_aide ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_AIDE: exponentially scaled derivative, Airy function Ai of an R8 argument.
//
//  Discussion:
//
//    if X <= 0,
//      R8_AIDE ( X ) = R8_AID ( X )
//    else
//      R8_AIDE ( X ) = R8_AID ( X ) * exp ( 2/3 * X^(3/2) )
//
//    Thanks to Aleksandra Piper for pointing out a correction involving 
//    the computation of Z in the last two cases, 02 February 2012.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 February 2012
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_AIDE, the exponentially scaled derivative of 
//    the Airy function Ai of X.
//
{
  double aifcs[13] = {
     0.105274612265314088088970057325134114,
      0.011836136281529978442889292583980840,
      0.000123281041732256643051689242469164,
      0.000000622612256381399016825658693579,
      0.000000001852988878441452950548140821,
      0.000000000003633288725904357915995625,
      0.000000000000005046217040440664768330,
      0.000000000000000005223816555471480985,
      0.000000000000000000004185745090748989,
      0.000000000000000000000002672887324883,
      0.000000000000000000000000001392128006,
      0.000000000000000000000000000000602653,
      0.000000000000000000000000000000000220 };
  double aigcs[13] = {
      0.0212338781509186668523122276848937,
      0.0863159303352144067524942809461604,
      0.0017975947203832313578033963225230,
      0.0000142654998755506932526620687495,
      0.0000000594379952836832010488787064,
      0.0000000001524033664794478945214786,
      0.0000000000002645876603490435305100,
      0.0000000000000003315624296815020591,
      0.0000000000000000003139789757594792,
      0.0000000000000000000002325767379040,
      0.0000000000000000000000001384384231,
      0.0000000000000000000000000000676629,
      0.0000000000000000000000000000000276 };
  double aip1cs[57] = {
     0.0358865097808301537956710489261688,
     0.0114668575627764898572700883121766,
    -0.0007592073583861400301335647601603,
     0.0000869517610893841271948619434021,
    -0.0000128237294298591691789607600486,
     0.0000022062695681038336934376250420,
    -0.0000004222295185920749486945988432,
     0.0000000874686415726348479356130376,
    -0.0000000192773588418365388625693417,
     0.0000000044668460054492719699777137,
    -0.0000000010790108051948168015747466,
     0.0000000002700029446696248083071434,
    -0.0000000000696480108007915257318929,
     0.0000000000184489907003246687076806,
    -0.0000000000050027817358071698301149,
     0.0000000000013852243366012168297298,
    -0.0000000000003908218466657048253473,
     0.0000000000001121536072524563451273,
    -0.0000000000000326861522579502522443,
     0.0000000000000096619179010090805752,
    -0.0000000000000028934767442698434271,
     0.0000000000000008770086661150897069,
    -0.0000000000000002688046261195853754,
     0.0000000000000000832498823872342992,
    -0.0000000000000000260343254786947057,
     0.0000000000000000082159528142686287,
    -0.0000000000000000026150406704984940,
     0.0000000000000000008390563463261051,
    -0.0000000000000000002712685618629660,
     0.0000000000000000000883333375271942,
    -0.0000000000000000000289603206822333,
     0.0000000000000000000095562185928676,
    -0.0000000000000000000031727463569051,
     0.0000000000000000000010595576960768,
    -0.0000000000000000000003558253765402,
     0.0000000000000000000001201334680517,
    -0.0000000000000000000000407666883800,
     0.0000000000000000000000139016944446,
    -0.0000000000000000000000047628165730,
     0.0000000000000000000000016391265551,
    -0.0000000000000000000000005665491354,
     0.0000000000000000000000001966381969,
    -0.0000000000000000000000000685230229,
     0.0000000000000000000000000239706939,
    -0.0000000000000000000000000084166831,
     0.0000000000000000000000000029659364,
    -0.0000000000000000000000000010487947,
     0.0000000000000000000000000003721150,
    -0.0000000000000000000000000001324570,
     0.0000000000000000000000000000472976,
    -0.0000000000000000000000000000169405,
     0.0000000000000000000000000000060855,
    -0.0000000000000000000000000000021924,
     0.0000000000000000000000000000007920,
    -0.0000000000000000000000000000002869,
     0.0000000000000000000000000000001042,
    -0.0000000000000000000000000000000379 };
  double aip2cs[37] = {
     0.0065457691989713756794276979067064,
     0.0023833724120774591992772552886923,
    -0.0000430700770220585862775012110584,
     0.0000015629125858629202330785369063,
    -0.0000000815417186162706965112501015,
     0.0000000054103738056935918208008783,
    -0.0000000004284130882614696528766222,
     0.0000000000389497962832286424862198,
    -0.0000000000039623161264979257658071,
     0.0000000000004428184214405989602353,
    -0.0000000000000536296527150689675318,
     0.0000000000000069649872139936028200,
    -0.0000000000000009619636286095319210,
     0.0000000000000001403454967784808032,
    -0.0000000000000000215097136525875715,
     0.0000000000000000034471230632678283,
    -0.0000000000000000005753907621819442,
     0.0000000000000000000997001165824168,
    -0.0000000000000000000178811436021458,
     0.0000000000000000000033110307923551,
    -0.0000000000000000000006315885529506,
     0.0000000000000000000001238666952364,
    -0.0000000000000000000000249324053394,
     0.0000000000000000000000051426030999,
    -0.0000000000000000000000010854236402,
     0.0000000000000000000000002341316852,
    -0.0000000000000000000000000515542099,
     0.0000000000000000000000000115758841,
    -0.0000000000000000000000000026479669,
     0.0000000000000000000000000006165328,
    -0.0000000000000000000000000001459931,
     0.0000000000000000000000000000351331,
    -0.0000000000000000000000000000085863,
     0.0000000000000000000000000000021297,
    -0.0000000000000000000000000000005358,
     0.0000000000000000000000000000001367,
    -0.0000000000000000000000000000000353 };
  double eta;
  static int naif = 0;
  static int naig = 0;
  static int naip1 = 0;
  static int naip2 = 0;
  double phi;
  double sqrtx;
  double value;
  double x2;
  static double x2sml = 0.0;
  double x3;
  static double x32sml = 0.0;
  static double x3sml = 0.0;
  static double xbig = 0.0;
  double xn;
  double z;

  if ( naif == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    naif = r8_inits ( aifcs, 13, eta );
    naig = r8_inits ( aigcs, 13, eta );
    naip1 = r8_inits ( aip1cs, 57, eta );
    naip2 = r8_inits ( aip2cs, 37, eta );
    x2sml = sqrt ( eta );
    x3sml = r8_power ( eta, 0.3333 );
    x32sml = 1.3104 * x3sml * x3sml;
    xbig = r8_power ( r8_mach ( 2 ), 0.6666 );
  }

  if ( x < - 1.0 )
  {
    r8_admp ( x, xn, phi );
    value = xn * cos ( phi );
  }
  else if ( r8_abs ( x ) < x2sml )
  {
    x2 = 0.0;
    x3 = 0.0;
    value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) ) 
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25;
    if ( x32sml < x )
    {
      value = value * exp ( 2.0 * x * sqrt ( x ) / 3.0 );
    }
  }
  else if ( r8_abs ( x ) < x3sml )
  {
    x2 = x * x;
    x3 = 0.0;
    value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) ) 
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25;
    if ( x32sml < x )
    {
      value = value * exp ( 2.0 * x * sqrt ( x ) / 3.0 );
    }
  }
  else if ( x <= 1.0 )
  {
    x2 = x * x;
    x3 = x * x;
    value = ( x2 * ( 0.125 + r8_csevl ( x3, aifcs, naif ) ) 
      - r8_csevl ( x3, aigcs, naig ) ) - 0.25;
    if ( x32sml < x )
    {
      value = value * exp ( 2.0 * x * sqrt ( x ) / 3.0 );
    }
  }
  else if ( x <= 4.0 )
  {
    sqrtx = sqrt ( x );
    z = ( 16.0  / ( x * sqrtx ) - 9.0 ) / 7.0;
    value = ( - 0.28125 - r8_csevl ( z, aip1cs, naip1 ) ) * sqrt ( sqrtx );
  }
  else if ( x < xbig )
  {
    sqrtx = sqrt ( x );
    z = 16.0  / ( x * sqrtx ) - 1.0;
    value = ( - 0.28125 - r8_csevl ( z, aip2cs, naip2 ) ) * sqrt ( sqrtx );
  }
  else
  {
    sqrtx = sqrt ( x );
    z = - 1.0;
    value = ( - 0.28125 - r8_csevl ( z, aip2cs, naip2 ) ) * sqrt ( sqrtx );
  }
  return value;
}
//****************************************************************************80

double r8_aie ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_AIE evaluates the exponentially scaled Airy function Ai of an R8 argument.
//
//  Discussion:
//
//    if X <= 0,
//      R8_AIE ( X ) = R8_AI ( X )
//    else
//      R8_AIE ( X ) = R8_AI ( X ) * exp ( 2/3 * X^(3/2) )
//
//    Thanks to Aleksandra Piper for pointing out a correction involving a
//    missing assignment to SQRTX, 27 January 2012.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2012
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_AIE, the exponentially scaled Airy function Ai of X.
//
{
  static double aifcs[13] = {
    -0.37971358496669997496197089469414E-01,
    +0.59191888537263638574319728013777E-01,
    +0.98629280577279975365603891044060E-03,
    +0.68488438190765667554854830182412E-05,
    +0.25942025962194713019489279081403E-07,
    +0.61766127740813750329445749697236E-10,
    +0.10092454172466117901429556224601E-12,
    +0.12014792511179938141288033225333E-15,
    +0.10882945588716991878525295466666E-18,
    +0.77513772196684887039238400000000E-22,
    +0.44548112037175638391466666666666E-25,
    +0.21092845231692343466666666666666E-28,
    +0.83701735910741333333333333333333E-32 };
  static double aigcs[13] = {
    +0.18152365581161273011556209957864E-01,
    +0.21572563166010755534030638819968E-01,
    +0.25678356987483249659052428090133E-03,
    +0.14265214119792403898829496921721E-05,
    +0.45721149200180426070434097558191E-08,
    +0.95251708435647098607392278840592E-11,
    +0.13925634605771399051150420686190E-13,
    +0.15070999142762379592306991138666E-16,
    +0.12559148312567778822703205333333E-19,
    +0.83063073770821340343829333333333E-23,
    +0.44657538493718567445333333333333E-26,
    +0.19900855034518869333333333333333E-29,
    +0.74702885256533333333333333333333E-33 };
  static double aip1cs[57] = {
    -0.2146951858910538455460863467778E-01,
    -0.7535382535043301166219720865565E-02,
    +0.5971527949026380852035388881994E-03,
    -0.7283251254207610648502368291548E-04,
    +0.1110297130739299666517381821140E-04,
    -0.1950386152284405710346930314033E-05,
    +0.3786973885159515193885319670057E-06,
    -0.7929675297350978279039072879154E-07,
    +0.1762247638674256075568420122202E-07,
    -0.4110767539667195045029896593893E-08,
    +0.9984770057857892247183414107544E-09,
    -0.2510093251387122211349867730034E-09,
    +0.6500501929860695409272038601725E-10,
    -0.1727818405393616515478877107366E-10,
    +0.4699378842824512578362292872307E-11,
    -0.1304675656297743914491241246272E-11,
    +0.3689698478462678810473948382282E-12,
    -0.1061087206646806173650359679035E-12,
    +0.3098414384878187438660210070110E-13,
    -0.9174908079824139307833423547851E-14,
    +0.2752049140347210895693579062271E-14,
    -0.8353750115922046558091393301880E-15,
    +0.2563931129357934947568636168612E-15,
    -0.7950633762598854983273747289822E-16,
    +0.2489283634603069977437281175644E-16,
    -0.7864326933928735569664626221296E-17,
    +0.2505687311439975672324470645019E-17,
    -0.8047420364163909524537958682241E-18,
    +0.2604097118952053964443401104392E-18,
    -0.8486954164056412259482488834184E-19,
    +0.2784706882142337843359429186027E-19,
    -0.9195858953498612913687224151354E-20,
    +0.3055304318374238742247668225583E-20,
    -0.1021035455479477875902177048439E-20,
    +0.3431118190743757844000555680836E-21,
    -0.1159129341797749513376922463109E-21,
    +0.3935772844200255610836268229154E-22,
    -0.1342880980296717611956718989038E-22,
    +0.4603287883520002741659190305314E-23,
    -0.1585043927004064227810772499387E-23,
    +0.5481275667729675908925523755008E-24,
    -0.1903349371855047259064017948945E-24,
    +0.6635682302374008716777612115968E-25,
    -0.2322311650026314307975200986453E-25,
    +0.8157640113429179313142743695359E-26,
    -0.2875824240632900490057489929557E-26,
    +0.1017329450942901435079714319018E-26,
    -0.3610879108742216446575703490559E-27,
    +0.1285788540363993421256640342698E-27,
    -0.4592901037378547425160693022719E-28,
    +0.1645597033820713725812102485333E-28,
    -0.5913421299843501842087920271360E-29,
    +0.2131057006604993303479369509546E-29,
    -0.7701158157787598216982761745066E-30,
    +0.2790533307968930417581783777280E-30,
    -0.1013807715111284006452241367039E-30,
    +0.3692580158719624093658286216533E-31 };
  static double aip2cs[37] = {
    -0.174314496929375513390355844011E-02,
    -0.167893854325541671632190613480E-02,
    +0.359653403352166035885983858114E-04,
    -0.138081860273922835457399383100E-05,
    +0.741122807731505298848699095233E-07,
    -0.500238203900133013130422866325E-08,
    +0.400693917417184240675446866355E-09,
    -0.367331242795905044199318496207E-10,
    +0.376034439592373852439592002918E-11,
    -0.422321332718747538026564938968E-12,
    +0.513509454033657070919618754120E-13,
    -0.669095850390477595651681356676E-14,
    +0.926667545641290648239550724382E-15,
    -0.135514382416070576333397356591E-15,
    +0.208115496312830995299006549335E-16,
    -0.334116499159176856871277570256E-17,
    +0.558578584585924316868032946585E-18,
    -0.969219040152365247518658209109E-19,
    +0.174045700128893206465696557738E-19,
    -0.322640979731130400247846333098E-20,
    +0.616074471106625258533259618986E-21,
    -0.120936347982490059076420676266E-21,
    +0.243632763310138108261570095786E-22,
    -0.502914221497457468943403144533E-23,
    +0.106224175543635689495470626133E-23,
    -0.229284284895989241509856324266E-24,
    +0.505181733929503744986884778666E-25,
    -0.113498123714412404979793920000E-25,
    +0.259765565985606980698374144000E-26,
    -0.605124621542939506172231679999E-27,
    +0.143359777966772800720295253333E-27,
    -0.345147757060899986280721066666E-28,
    +0.843875190213646740427025066666E-29,
    -0.209396142298188169434453333333E-29,
    +0.527008873478945503182848000000E-30,
    -0.134457433014553385789030399999E-30,
    +0.347570964526601147340117333333E-31 };
  double eta;
  static int naif = 0;
  static int naig =  0;
  static int naip1 = 0;
  static int naip2 = 0;
  double sqrtx;
  double theta;
  double value;
  static double x32sml = 0.0;
  static double x3sml = 0.0;
  static double xbig = 0.0;
  double xm;
  double z;

  if ( naif == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    naif = r8_inits ( aifcs, 13, eta );
    naig = r8_inits ( aigcs, 13, eta );
    naip1 = r8_inits ( aip1cs, 57, eta );
    naip2 = r8_inits ( aip2cs, 37, eta );
    x3sml = r8_power ( eta, 0.3333 );
    x32sml = 1.3104 * x3sml * x3sml;
    xbig = r8_power ( r8_mach ( 2 ), 0.6666 );
  }

  if ( x < - 1.0 )
  {
    r8_aimp ( x, xm, theta );
    value = xm * cos ( theta );
  }
  else if ( 0.0 <= x && x <= x32sml )
  {
    z = 0.0;
    value = 0.3750 + ( r8_csevl ( z, aifcs, naif ) 
      - x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
  }
  else if ( r8_abs ( x ) <= x3sml )
  {
    z = 0.0;
    value = 0.3750 + ( r8_csevl ( z, aifcs, naif ) 
      - x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
    value = value * exp ( 2.0 * x * sqrt ( x ) / 3.0 );
  }
  else if ( x <= 1.0 )
  {
    z = x * x * x;
    value = 0.3750 + ( r8_csevl ( z, aifcs, naif ) 
      - x * ( 0.25 + r8_csevl ( z, aigcs, naig ) ) );
    value = value * exp ( 2.0 * x * sqrt ( x ) / 3.0 );
  }
  else if ( x <= 4.0 )
  {
    sqrtx = sqrt ( x );
    z = ( 16.0 / ( x * sqrtx ) - 9.0 ) / 7.0;
    value = ( 0.28125 + r8_csevl ( z, aip1cs, naip1 ) ) / sqrt ( sqrtx );
  }
  else if ( x < xbig )
  {
    sqrtx = sqrt ( x );
    z = 16.0 / ( x * sqrtx ) - 1.0;
    value = ( 0.28125 + r8_csevl ( z, aip2cs, naip2 ) ) / sqrt ( sqrtx );
  }
  else
  {
    sqrtx = sqrt ( x );
    z = - 1.0;
    value = ( 0.28125 + r8_csevl ( z, aip2cs, naip2 ) ) / sqrt ( sqrtx );
  }
  return value;
}
//****************************************************************************80

void r8_aimp ( double x, double &ampl, double &theta )

//****************************************************************************80
//
//  Purpose:
//
//    R8_AIMP evaluates the modulus and phase of the Airy function.
//
//  Description:
//
//    This function must only be called when X <= -1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double &AMPL, &PHI, the modulus and phase of the 
//    Airy function.
//
{
  static double am20cs[57] = {
    +0.108716749086561856615730588125E-01,
    +0.369489228982663555091728665146E-03,
    +0.440680100484689563667507001327E-05,
    +0.143686762361911153929183952833E-06,
    +0.824275552390078308670628855353E-08,
    +0.684426758893661606173927278180E-09,
    +0.739566697282739287731004740213E-10,
    +0.974595633696825017638702600847E-11,
    +0.150076885829405775650973119497E-11,
    +0.262147910221527634206252854802E-12,
    +0.508354111376487180357278966914E-13,
    +0.107684753358811440492985997070E-13,
    +0.246091286618433429335914062617E-14,
    +0.600786380358656418436110373550E-15,
    +0.155449156102388071150651388384E-15,
    +0.423535125035576604426382780182E-16,
    +0.120862166289299840154401109189E-16,
    +0.359609651214658240861499706423E-17,
    +0.111134218386395638261774604677E-17,
    +0.355559532432366609893680289225E-18,
    +0.117433021600139309998766947387E-18,
    +0.399397454661077561389162200966E-19,
    +0.139576671528916310425606325640E-19,
    +0.500240055309236041393459280716E-20,
    +0.183552760958132679184834866457E-20,
    +0.688490998179202743197790112404E-21,
    +0.263631035611417012359996885105E-21,
    +0.102924890237338360287153563785E-21,
    +0.409246966671594885489762960571E-22,
    +0.165558573406734651039727903828E-22,
    +0.680797467063033356116599685727E-23,
    +0.284326559934079832419751134476E-23,
    +0.120507398348965255097287818819E-23,
    +0.517961243287505217976613610424E-24,
    +0.225622613427562816303268640887E-24,
    +0.995418801147745168832117078246E-25,
    +0.444551696397342424308280582053E-25,
    +0.200865195461501101425916097338E-25,
    +0.917786344151775165973885645402E-26,
    +0.423872958105589240661672197948E-26,
    +0.197789272007846092370846251490E-26,
    +0.932116351284620665680435253373E-27,
    +0.443482133249918099955611379722E-27,
    +0.212945672365573895594589552837E-27,
    +0.103158569651075977552209344907E-27,
    +0.504023773022591199157904590029E-28,
    +0.248301304570155945304046541005E-28,
    +0.123301783128562196054198238560E-28,
    +0.617033449920521746121976730507E-29,
    +0.311092617415918897233869792213E-29,
    +0.157983085201706173015269071503E-29,
    +0.807931987538283607678121339092E-30,
    +0.415997394138667562722951360052E-30,
    +0.215610934097716900471935862504E-30,
    +0.112468857265869178296752823613E-30,
    +0.590331560632838091123040811797E-31,
    +0.311735667692928562046280505333E-31 };
  static double am21cs[60] = {
    +0.592790266721309588375717482814E-02,
    +0.200569405393165186428695217690E-02,
    +0.911081850262275893553072526291E-04,
    +0.849894306372047155633172107475E-05,
    +0.113297908976913076637929215494E-05,
    +0.187517946100666496180950627804E-06,
    +0.359306519018245832699035211192E-07,
    +0.765757714071683864039093517470E-08,
    +0.176999967168039173925953460744E-08,
    +0.436259555654598932720546585535E-09,
    +0.113291641337853230035520085219E-09,
    +0.307257690982419244137868398126E-10,
    +0.864482416482201075541200465766E-11,
    +0.251015250060924402115104562212E-11,
    +0.749102496764440371601802227751E-12,
    +0.228996928487994073089565214432E-12,
    +0.715113658927987694949327491175E-13,
    +0.227607924959566841946395165061E-13,
    +0.736942142760886513969953227782E-14,
    +0.242328675267827490463991742006E-14,
    +0.808153774548239869283406558403E-15,
    +0.273008079804356086659174563386E-15,
    +0.933236070891385318473519474326E-16,
    +0.322508099681084622213867546973E-16,
    +0.112581932346444541217757573416E-16,
    +0.396699463986938821660259459530E-17,
    +0.141006567944319504660865034527E-17,
    +0.505302086537851213375537393032E-18,
    +0.182461523215945141197999102789E-18,
    +0.663584568262130466928029121642E-19,
    +0.242963731631276179741747455826E-19,
    +0.895238915123687802013669922963E-20,
    +0.331845289350050791260229250755E-20,
    +0.123706196188658315384437905922E-20,
    +0.463636677012390840306767734243E-21,
    +0.174653135947764475469758765989E-21,
    +0.661116810234991176307910643111E-22,
    +0.251409918994072486176125666459E-22,
    +0.960274995571732568694034386998E-23,
    +0.368324952289296395686436898078E-23,
    +0.141843138269159136145535939553E-23,
    +0.548342674276935830106345800990E-24,
    +0.212761054623118806650372562616E-24,
    +0.828443700849418591487734760953E-25,
    +0.323670563926127001421028600927E-25,
    +0.126868882963286057355055062493E-25,
    +0.498843818992121626935068934362E-26,
    +0.196734584467649390967119381790E-26,
    +0.778135971020326957713212064836E-27,
    +0.308633941498911152919192968451E-27,
    +0.122744647045453119789338037234E-27,
    +0.489431279134292205885241216204E-28,
    +0.195646879802909821175925099724E-28,
    +0.783988952922426171166311492266E-29,
    +0.314896914002484223748298978099E-29,
    +0.126769763137250681307067842559E-29,
    +0.511470691906900141641632107724E-30,
    +0.206801709795538770250900316706E-30,
    +0.837891344768519001325996867583E-31,
    +0.340168991971489802052339079577E-31 };
  static double am22cs[74] = {
    -0.156284448062534112753545828583E-01,
    +0.778336445239681307018943100334E-02,
    +0.867057770477189528406072812110E-03,
    +0.156966273156113719469953482266E-03,
    +0.356396257143286511324100666302E-04,
    +0.924598335425043154495080090994E-05,
    +0.262110161850422389523194982066E-05,
    +0.791882216516012561489469982263E-06,
    +0.251041527921011847803162690862E-06,
    +0.826522320665407734472997712940E-07,
    +0.280571166281305264396384290014E-07,
    +0.976821090484680786674631273890E-08,
    +0.347407923227710343287279035573E-08,
    +0.125828132169836914219092738164E-08,
    +0.462988260641895264497330784625E-09,
    +0.172728258813604072468143128696E-09,
    +0.652319200131154135148574124970E-10,
    +0.249047168520982056019881087112E-10,
    +0.960156820553765948078189890126E-11,
    +0.373448002067726856974776596757E-11,
    +0.146417565032053391722216189678E-11,
    +0.578265471168512825475827881553E-12,
    +0.229915407244706118560254184494E-12,
    +0.919780711231997257150883662365E-13,
    +0.370060068813090065807504045556E-13,
    +0.149675761698672987823326345205E-13,
    +0.608361194938461148720451399443E-14,
    +0.248404087115121397635425326873E-14,
    +0.101862476526769080727914465339E-14,
    +0.419383856352753989429640310957E-15,
    +0.173318901762930756149702493501E-15,
    +0.718821902388508517820445406811E-16,
    +0.299123633598403607712470896113E-16,
    +0.124868990433238627855713110880E-16,
    +0.522829344609483661928651193632E-17,
    +0.219532961724713396595998454359E-17,
    +0.924298325229777281154410024332E-18,
    +0.390157708236091407825543197309E-18,
    +0.165093892693863707213759030367E-18,
    +0.700221815715994367565716554487E-19,
    +0.297651833616786915573214963506E-19,
    +0.126796539086902072571134261229E-19,
    +0.541243400697077628687581725061E-20,
    +0.231487350218155252296382133283E-20,
    +0.991920288386566563462623851167E-21,
    +0.425803015323732357158897608174E-21,
    +0.183101842973024501678402003088E-21,
    +0.788678712311075375564526811022E-22,
    +0.340254607386229874956582997235E-22,
    +0.147020881405712530791860892535E-22,
    +0.636211018324916957733348071767E-23,
    +0.275707050680980721919395987768E-23,
    +0.119645858090104071356261780457E-23,
    +0.519912545729242147981768210567E-24,
    +0.226217674847104475260575286850E-24,
    +0.985526113754431819448565068283E-25,
    +0.429870630332508717223681286187E-25,
    +0.187723641661580639829657670189E-25,
    +0.820721941772842137268801052115E-26,
    +0.359214665604615507812767944463E-26,
    +0.157390594612773315611458940587E-26,
    +0.690329781039333834965319153586E-27,
    +0.303092079078968534607859331415E-27,
    +0.133204934160481219185689121944E-27,
    +0.585978836851523490117937981442E-28,
    +0.258016868489487806338425080457E-28,
    +0.113712433637283667223632182863E-28,
    +0.501592557226068509236430548549E-29,
    +0.221445829395509373322569708484E-29,
    +0.978470283886507289984691416411E-30,
    +0.432695414934180170112000952983E-30,
    +0.191497288193994570612929860440E-30,
    +0.848164622402392354171298331562E-31,
    +0.375947065173955919947455052934E-31 };
  static double ath0cs[53] = {
    -0.8172601764161634499840208700543E-01,
    -0.8004012824788273287596481113068E-03,
    -0.3186525268782113203795553628242E-05,
    -0.6688388266477509330741698865033E-07,
    -0.2931759284994564516506822463184E-08,
    -0.2011263760883621669049030307186E-09,
    -0.1877522678055973426074008166652E-10,
    -0.2199637137704601251899002199848E-11,
    -0.3071616682592272449025746605586E-12,
    -0.4936140553673418361025600985389E-13,
    -0.8902833722583660416935236969866E-14,
    -0.1768987764615272613656814199467E-14,
    -0.3817868689032277014678199609600E-15,
    -0.8851159014819947594156286509984E-16,
    -0.2184818181414365953149677679568E-16,
    -0.5700849046986452380599442295119E-17,
    -0.1563121122177875392516031795495E-17,
    -0.4481437996768995067906688776353E-18,
    -0.1337794883736188022044566044098E-18,
    -0.4143340036874114453776852445442E-19,
    -0.1327263385718805025080481164652E-19,
    -0.4385728589128440522215756835955E-20,
    -0.1491360695952818067686201743956E-20,
    -0.5208104738630711377154238188773E-21,
    -0.1864382222390498923872526604979E-21,
    -0.6830263751167969012975435381881E-22,
    -0.2557117058029329629296207591347E-22,
    -0.9770158640254300218246907254046E-23,
    -0.3805161433416679084068428254886E-23,
    -0.1509022750737054063493926482995E-23,
    -0.6087551341242424929005568014525E-24,
    -0.2495879513809711495425982124058E-24,
    -0.1039157654581920948909588084274E-24,
    -0.4390235913976846536974594969051E-25,
    -0.1880790678447990211675826820582E-25,
    -0.8165070764199462948863022205753E-26,
    -0.3589944503749750514266435585041E-26,
    -0.1597658126632132872981291608708E-26,
    -0.7193250175703823969113802835305E-27,
    -0.3274943012727856506209351132721E-27,
    -0.1507042445783690665816975047272E-27,
    -0.7006624198319904717843967949140E-28,
    -0.3289907402983718226528815678356E-28,
    -0.1559518084365146526445322711496E-28,
    -0.7460690508208254582833851119721E-29,
    -0.3600877034824662020563277249431E-29,
    -0.1752851437473772257350402219197E-29,
    -0.8603275775188512909623778628724E-30,
    -0.4256432603226946534668039480105E-30,
    -0.2122161865044262927723650698206E-30,
    -0.1065996156704879052472060798561E-30,
    -0.5393568608816949116410688086892E-31,
    -0.2748174851043954822278496517870E-31 };
  static double ath1cs[58] = {
    -0.6972849916208883845888148415037E-01,
    -0.5108722790650044987073448077961E-02,
    -0.8644335996989755094525334749512E-04,
    -0.5604720044235263542188698916125E-05,
    -0.6045735125623897409156376640077E-06,
    -0.8639802632488334393219721138499E-07,
    -0.1480809484309927157147782480780E-07,
    -0.2885809334577236039999449908712E-08,
    -0.6191631975665699609309191231800E-09,
    -0.1431992808860957830931365259879E-09,
    -0.3518141102137214721504616874321E-10,
    -0.9084761919955078290070339808051E-11,
    -0.2446171672688598449343283664767E-11,
    -0.6826083203213446240828996710264E-12,
    -0.1964579931194940171278546257802E-12,
    -0.5808933227139693164009191265856E-13,
    -0.1759042249527441992795400959024E-13,
    -0.5440902932714896613632538945319E-14,
    -0.1715247407486806802622358519451E-14,
    -0.5500929233576991546871101847161E-15,
    -0.1791878287739317259495152638754E-15,
    -0.5920372520086694197778411062231E-16,
    -0.1981713027876483962470972206590E-16,
    -0.6713232347016352262049984343790E-17,
    -0.2299450243658281116122358619832E-17,
    -0.7957300928236376595304637145634E-18,
    -0.2779994027291784157172290233739E-18,
    -0.9798924361326985224406795480814E-19,
    -0.3482717006061574386702645565849E-19,
    -0.1247489122558599057173300058084E-19,
    -0.4501210041478228113487751824452E-20,
    -0.1635346244013352135596114164667E-20,
    -0.5980102897780336268098762265941E-21,
    -0.2200246286286123454028196295475E-21,
    -0.8142463073515085897408205291519E-22,
    -0.3029924773660042537432330709674E-22,
    -0.1133390098574623537722943969689E-22,
    -0.4260766024749295719283049889791E-23,
    -0.1609363396278189718797500634453E-23,
    -0.6106377190825026293045330444287E-24,
    -0.2326954318021694061836577887573E-24,
    -0.8903987877472252604474129558186E-25,
    -0.3420558530005675024117914752341E-25,
    -0.1319026715257272659017212100607E-25,
    -0.5104899493612043091316191177386E-26,
    -0.1982599478474547451242444663466E-26,
    -0.7725702356880830535636111851519E-27,
    -0.3020234733664680100815776863573E-27,
    -0.1184379739074169993712946380800E-27,
    -0.4658430227922308520573252840106E-28,
    -0.1837554188100384647157502006613E-28,
    -0.7268566894427990953321876684800E-29,
    -0.2882863120391468135527089875626E-29,
    -0.1146374629459906350417591664639E-29,
    -0.4570031437748533058179991688533E-30,
    -0.1826276602045346104809934028799E-30,
    -0.7315349993385250469111066350933E-31,
    -0.2936925599971429781637815773866E-31 };
  static double ath2cs[72] = {
    +0.4405273458718778997061127057775E-02,
    -0.3042919452318454608483844239873E-01,
    -0.1385653283771793791602692842653E-02,
    -0.1804443908954952302670486910952E-03,
    -0.3380847108327308671057465323618E-04,
    -0.7678183535229023055257676817765E-05,
    -0.1967839443716035324690935417077E-05,
    -0.5483727115877700361586143659281E-06,
    -0.1625461550532612452712696212258E-06,
    -0.5053049981268895015277637842078E-07,
    -0.1631580701124066881183851715617E-07,
    -0.5434204112348517507963436694817E-08,
    -0.1857398556409900325763850109630E-08,
    -0.6489512033326108816213513640676E-09,
    -0.2310594885800944720482995987079E-09,
    -0.8363282183204411682819329546745E-10,
    -0.3071196844890191462660661303891E-10,
    -0.1142367142432716819409514579892E-10,
    -0.4298116066345803065822470108971E-11,
    -0.1633898699596715440601646086632E-11,
    -0.6269328620016619432123443754076E-12,
    -0.2426052694816257357356159203991E-12,
    -0.9461198321624039090742527765052E-13,
    -0.3716060313411504806847798281269E-13,
    -0.1469155684097526763170138810309E-13,
    -0.5843694726140911944556401363094E-14,
    -0.2337502595591951298832675034934E-14,
    -0.9399231371171435401160167358411E-15,
    -0.3798014669372894500076335263715E-15,
    -0.1541731043984972524883443681775E-15,
    -0.6285287079535307162925662365202E-16,
    -0.2572731812811455424755383992774E-16,
    -0.1057098119354017809340974866555E-16,
    -0.4359080267402696966695992699964E-17,
    -0.1803634315959978013953176945540E-17,
    -0.7486838064380536821719431676914E-18,
    -0.3117261367347604656799597209985E-18,
    -0.1301687980927700734792871620696E-18,
    -0.5450527587519522468973883909909E-19,
    -0.2288293490114231872268635931903E-19,
    -0.9631059503829538655655060440088E-20,
    -0.4063281001524614089092195416434E-20,
    -0.1718203980908026763900413858510E-20,
    -0.7281574619892536367415322473328E-21,
    -0.3092352652680643127960680345790E-21,
    -0.1315917855965440490383417023254E-21,
    -0.5610606786087055512664907412668E-22,
    -0.2396621894086355206020304337895E-22,
    -0.1025574332390581200832954423924E-22,
    -0.4396264138143656476403607323663E-23,
    -0.1887652998372577373342508719450E-23,
    -0.8118140359576807603579433230445E-24,
    -0.3496734274366286856375952089214E-24,
    -0.1508402925156873215171751475867E-24,
    -0.6516268284778671059787773834341E-25,
    -0.2818945797529207424505942114583E-25,
    -0.1221127596512262744598094464505E-25,
    -0.5296674341169867168620011705073E-26,
    -0.2300359270773673431358870971744E-26,
    -0.1000279482355367494781220348930E-26,
    -0.4354760404180879394806893162179E-27,
    -0.1898056134741477522515482827030E-27,
    -0.8282111868712974697554009309315E-28,
    -0.3617815493066569006586213484374E-28,
    -0.1582018896178003654858941843636E-28,
    -0.6925068597802270011772820383247E-29,
    -0.3034390239778629128908629727335E-29,
    -0.1330889568166725224761977446509E-29,
    -0.5842848522173090120487606971706E-30,
    -0.2567488423238302631121274357678E-30,
    -0.1129232322268882185791505819151E-30,
    -0.4970947029753336916550570105023E-31 };
  double eta;
  static int nam20 = 0;
  static int nam21 = 0;
  static int nam22 = 0;
  static int nath0 = 0;
  static int nath1 = 0;
  static int nath2 = 0;
  static double pi4 = 0.78539816339744830961566084581988;
  double sqrtx;
  static double xsml = 0.0;
  double z;

  if ( nam20 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nam20 = r8_inits ( am20cs, 57, eta );
    nath0 = r8_inits ( ath0cs, 53, eta );
    nam21 = r8_inits ( am21cs, 60, eta );
    nath1 = r8_inits ( ath1cs, 58, eta );
    nam22 = r8_inits ( am22cs, 74, eta );
    nath2 = r8_inits ( ath2cs, 72, eta );
    xsml = - r8_power ( 128.0 / r8_mach ( 3 ), 0.3333 );
  }

  if ( x <= xsml )
  {
    z = 1.0;
    ampl = 0.3125 + r8_csevl ( z, am20cs, nam20 );
    theta = - 0.625 + r8_csevl ( z, ath0cs, nath0 );
  }
  else if ( x < - 4.0 )
  {
    z = 128.0 / x / x / x + 1.0;
    ampl = 0.3125 + r8_csevl ( z, am20cs, nam20 );
    theta = - 0.625 + r8_csevl ( z, ath0cs, nath0 );
  }
  else if ( x < - 2.0 )
  {
    z = ( 128.0 / x / x / x + 9.0 ) / 7.0;
    ampl = 0.3125 + r8_csevl ( z, am21cs, nam21 );
    theta = - 0.625 + r8_csevl ( z, ath1cs, nath1 );
  }
  else if ( x <= - 1.0 )
  {
    z = ( 16.0 / x / x / x + 9.0 ) / 7.0;
    ampl = 0.3125 + r8_csevl ( z, am22cs, nam22 );
    theta = - 0.625 + r8_csevl ( z, ath2cs, nath2 );
  }
  else
  {
    cerr << "\n";
    cerr << "R8_AIMP - Fatal error!\n";
    cerr << "  -1.0 < X.\n";
    exit ( 1 );
  }

  sqrtx = sqrt ( - x );
  ampl = sqrt ( ampl / sqrtx );
  theta = pi4 - x * sqrtx * theta;

  return;
}
//****************************************************************************80

double r8_aint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_AINT truncates an R8 argument to an integer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    1 September 2011
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_AINT, the truncated version of X.
//
{
  double value;

  if ( x < 0.0E+00 )
  {
    value = - ( double ) ( ( int ) ( r8_abs ( x ) ) );
  }
  else
  {
    value =   ( double ) ( ( int ) ( r8_abs ( x ) ) );
  }

  return value;
}
//****************************************************************************80

double r8_asin ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ASIN evaluates the arc-sine of an R8 argument.
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
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ASIN, the arc-sine of X.
//
{
  static double asincs[39] = {
    +0.10246391753227159336573148305785E+00,
    +0.54946487221245833306011195902924E-01,
    +0.40806303925449692851307056149246E-02,
    +0.40789006854604435455598823905612E-03,
    +0.46985367432203691616048530136218E-04,
    +0.58809758139708058986454385552074E-05,
    +0.77732312462777632750557528163795E-06,
    +0.10677423340082039235047504956587E-06,
    +0.15092399536022808262386434401064E-07,
    +0.21809724080055385496609614713930E-08,
    +0.32075984262789614433261959667376E-09,
    +0.47855369646781034461493133918953E-10,
    +0.72251287362910432263848754537112E-11,
    +0.11018334742255783705372701334987E-11,
    +0.16947632539203354877423745651078E-12,
    +0.26261558667348224162283241502416E-13,
    +0.40958299813281178408828069291110E-14,
    +0.64244793108803655891727944887091E-15,
    +0.10128142198228221693973361222041E-15,
    +0.16039221897380787560050597464746E-16,
    +0.25503501355807141715298789676373E-17,
    +0.40701403797862382855487165672106E-18,
    +0.65172671712881144437889267575466E-19,
    +0.10467453037096796954244891716266E-19,
    +0.16858725563380328094989095185066E-20,
    +0.27221936305040227625164341247999E-21,
    +0.44059293900347550617126830079999E-22,
    +0.71466685243375937853063168000000E-23,
    +0.11615793343859516051798971733333E-23,
    +0.18915234552354685801184187733333E-24,
    +0.30855772044244342399827968000000E-25,
    +0.50416366022162453412970495999999E-26,
    +0.82502725502400865081753600000000E-27,
    +0.13520032631020947208055466666666E-27,
    +0.22184326876541720216644266666666E-28,
    +0.36442494054085079212578133333333E-29,
    +0.59920218558643813307733333333333E-30,
    +0.98584812059573785810261333333333E-31,
    +0.16222501166399014393173333333333E-31 };
  static int nterms = 0;
  static double pi2 = 1.57079632679489661923132169163975;
  static double sqeps = 0.0;
  double value;
  double y;
  double z;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( asincs, 39, 0.1 * r8_mach ( 3 ) );
    sqeps = sqrt ( 6.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( x < - 1.0 - sqeps )
  {
    cerr << "\n";
    cerr << "R8_ASIN - Fatal error!\n";
    cerr << "  X < - 1.0\n";
    exit ( 1 );
  }
  else if ( x < - 1.0 )
  {
    value = - pi2;
  }
  else if ( x < 1.0 )
  {
    z = 0.0;
    if ( sqeps < y )
    {
      z = y * y;
    }

    if ( z <= 0.5 )
    {
      value = x * ( 1.0 + r8_csevl ( 4.0 * z - 1.0, asincs, nterms ) );
    }
    else
    {
      value = pi2 - sqrt ( 1.0 - z ) * ( 1.0 + 
        r8_csevl ( 3.0 - 4.0 * z, asincs, nterms ) );
    }

    if ( x < 0.0 )
    {
      value = - r8_abs ( value );
    }
    else if ( 0.0 < x )
    {
      value = + r8_abs ( value );
    }
  }
  else if ( x < 1.0 + sqeps )
  {
    value = pi2;
  }
  else
  {
    cerr << "\n";
    cerr << "R8_ASIN - Fatal error!\n";
    cerr << "  1.0 < X\n";
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

double r8_asinh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ASINH evaluates the arc-sine of an R8 argument.
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
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ASINH, the arc-hyperbolic sine of X.
//
{
  static double aln2 = 0.69314718055994530941723212145818;
  static double asnhcs[39] = {
    -0.12820039911738186343372127359268E+00,
    -0.58811761189951767565211757138362E-01,
    +0.47274654322124815640725249756029E-02,
    -0.49383631626536172101360174790273E-03,
    +0.58506207058557412287494835259321E-04,
    -0.74669983289313681354755069217188E-05,
    +0.10011693583558199265966192015812E-05,
    -0.13903543858708333608616472258886E-06,
    +0.19823169483172793547317360237148E-07,
    -0.28847468417848843612747272800317E-08,
    +0.42672965467159937953457514995907E-09,
    -0.63976084654366357868752632309681E-10,
    +0.96991686089064704147878293131179E-11,
    -0.14844276972043770830246658365696E-11,
    +0.22903737939027447988040184378983E-12,
    -0.35588395132732645159978942651310E-13,
    +0.55639694080056789953374539088554E-14,
    -0.87462509599624678045666593520162E-15,
    +0.13815248844526692155868802298129E-15,
    -0.21916688282900363984955142264149E-16,
    +0.34904658524827565638313923706880E-17,
    -0.55785788400895742439630157032106E-18,
    +0.89445146617134012551050882798933E-19,
    -0.14383426346571317305551845239466E-19,
    +0.23191811872169963036326144682666E-20,
    -0.37487007953314343674570604543999E-21,
    +0.60732109822064279404549242880000E-22,
    -0.98599402764633583177370173440000E-23,
    +0.16039217452788496315232638293333E-23,
    -0.26138847350287686596716134399999E-24,
    +0.42670849606857390833358165333333E-25,
    -0.69770217039185243299730773333333E-26,
    +0.11425088336806858659812693333333E-26,
    -0.18735292078860968933021013333333E-27,
    +0.30763584414464922794065920000000E-28,
    -0.50577364031639824787046399999999E-29,
    +0.83250754712689142224213333333333E-30,
    -0.13718457282501044163925333333333E-30,
    +0.22629868426552784104106666666666E-31 };
  static int nterms = 0;
  static double sqeps = 0.0;
  double value;
  static double xmax = 0.0;
  double y;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( asnhcs, 39, 0.1 * r8_mach ( 3 ) );
    sqeps = sqrt ( r8_mach ( 3 ) );
    xmax = 1.0 / sqeps;
  }

  y = r8_abs ( x );

  if ( y <= sqeps )
  {
    value = x;
  }
  else if ( y <= 1.0 )
  {
    value = x * ( 1.0 + r8_csevl ( 2.0 * x * x - 1.0, asnhcs, nterms ) );
  }
  else if ( y < xmax )
  {
    value = log ( y + sqrt ( y * y + 1.0 ) );
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  else
  {
    value = aln2 + log ( y );
    if ( x < 0.0 );
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

double r8_atan ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ATAN evaluates the arc-tangent of an R8 argument.
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
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ATAN, the arc-tangent of X.
//
{
  static double atancs[16] = {
    +0.48690110349241406474636915902891E+00,
    -0.65108316367174641818869794945974E-02,
    +0.38345828265245177653569992430456E-04,
    -0.26872212876223146539595410518788E-06,
    +0.20500930985824269846636514686688E-08,
    -0.16450717395484269455734135285348E-10,
    +0.13650975274390773423813528484428E-12,
    -0.11601779591998246322891309834666E-14,
    +0.10038333943866273835797657402666E-16,
    -0.88072747152163859327073696000000E-19,
    +0.78136321005661722180580266666666E-21,
    -0.69954535148267456086613333333333E-23,
    +0.63105905713702136004266666666666E-25,
    -0.57296075370213874346666666666666E-27,
    +0.52274796280602282666666666666666E-29,
    -0.48327903911608320000000000000000E-31 };
  static double conpi8[4] = {
    0.375,
    0.75,
    1.125,
    1.5 };
  int n;
  static int nterms = 0;
  static double pi8[4] = {
    +0.17699081698724154807830422909937E-01,
    +0.35398163397448309615660845819875E-01,
    +0.53097245096172464423491268729813E-01,
    +0.70796326794896619231321691639751E-01 };
  static double sqeps = 0.0;
  double t;
  static double tanp8[3] = {
    +0.41421356237309504880168872420969,
    +1.0,
    +2.4142135623730950488016887242096 };
  double value;
  static double xbig = 0.0;
  static double xbnd1 = +0.19891236737965800691159762264467;
  static double xbnd2 = +0.66817863791929891999775768652308;
  static double xbnd3 = +1.4966057626654890176011351349424;
  static double xbnd4 = +5.0273394921258481045149750710640;
  double y;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( atancs, 16, 0.1 * r8_mach ( 3 ) );
    sqeps = sqrt ( 6.0 * r8_mach ( 3 ) );
    xbig = 1.0 / r8_mach ( 3 );
  }

  y = r8_abs ( x );

  if ( y <= xbnd1 )
  {
    value = x;
    if ( sqeps < y )
    {
      value = x * ( 0.75 + r8_csevl ( 50.0 * y * y - 1.0, atancs, nterms ) );
    }
  }
  else if ( y <= xbnd4 )
  {
    if ( xbnd3 < y )
    {
      n = 3;
    }
    else if ( xbnd2 < y )
    {
      n = 2;
    }
    else
    {
      n = 1;
    }

    t = ( y - tanp8[n-1] ) / ( 1.0 + y * tanp8[n-1] );

    value = conpi8[n-1] + ( pi8[n-1] + t * ( 0.75 +
      r8_csevl ( 50.0 * t * t - 1.0, atancs, nterms ) ) );
  }
  else
  {
    value = conpi8[3] + pi8[3];

    if ( y < xbig )
    {
      value = conpi8[3] + ( pi8[3] - ( 0.75 +
        r8_csevl ( 50.0 / y / y - 1.0, atancs, nterms ) ) / y );
    }
  }

  if ( x < 0.0 )
  {
    value = - r8_abs ( value );
  }
  else
  {
    value = + r8_abs ( value );
  }

  return value;
}
//****************************************************************************80

double r8_atan2 ( double sn, double cs )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ATAN2 evaluates the arc-tangent of two R8 arguments.
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
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double SN, CS, the Y and X coordinates of a 
//    point on the angle.
//
//    Output, double R8_ATAN2, the arc-tangent of the angle.
//
{
  double abscs;
  double abssn;
  static double big = 0.0;
  static double pi = 3.14159265358979323846264338327950;
  static double sml = 0.0;
  double value;

  if ( sml == 0.0 )
  {
    sml = r8_mach ( 1 );
    big = r8_mach ( 2 );
  }
//
//  We now make sure SN can be divided by CS.  It is painful.
//
  abssn = r8_abs ( sn );
  abscs = r8_abs ( cs );

  if ( abscs <= abssn )
  {
    if ( abscs < 1.0 && abscs * big <= abssn )
    {
      if ( sn < 0.0 )
      {
        value = - 0.5 * pi;
      }
      else if ( sn == 0.0 )
      {
        cerr << "\n";
        cerr << "R8_ATAN2 - Fatal error!\n";
        cerr << "  Both arguments are 0.\n";
        exit ( 1 );
      }
      else
      {
        value = 0.5 * pi;
      }
      return value;
    }
  }
  else
  {
    if ( 1.0 < abscs && abssn <= abscs * sml )
    {
      if ( 0.0 <= cs )
      {
        value = 0.0;
      }
      else
      {
        value = pi;
      }
      return value;
    }
  }
  value = atan ( sn / cs );

  if ( cs < 0.0 )
  {
    value = value + pi;
  }

  if ( pi < value )
  {
    value = value - 2.0 * pi;
  }
  return value;
}
//****************************************************************************80

double r8_atanh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ATANH evaluates the arc-hyperbolic tangent of an R8 argument.
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
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ATANH, the arc-hyperbolic tangent of X.
//
{
  static double atnhcs[27] = {
    +0.9439510239319549230842892218633E-01,
    +0.4919843705578615947200034576668E-01,
    +0.2102593522455432763479327331752E-02,
    +0.1073554449776116584640731045276E-03,
    +0.5978267249293031478642787517872E-05,
    +0.3505062030889134845966834886200E-06,
    +0.2126374343765340350896219314431E-07,
    +0.1321694535715527192129801723055E-08,
    +0.8365875501178070364623604052959E-10,
    +0.5370503749311002163881434587772E-11,
    +0.3486659470157107922971245784290E-12,
    +0.2284549509603433015524024119722E-13,
    +0.1508407105944793044874229067558E-14,
    +0.1002418816804109126136995722837E-15,
    +0.6698674738165069539715526882986E-17,
    +0.4497954546494931083083327624533E-18,
    +0.3032954474279453541682367146666E-19,
    +0.2052702064190936826463861418666E-20,
    +0.1393848977053837713193014613333E-21,
    +0.9492580637224576971958954666666E-23,
    +0.6481915448242307604982442666666E-24,
    +0.4436730205723615272632320000000E-25,
    +0.3043465618543161638912000000000E-26,
    +0.2091881298792393474047999999999E-27,
    +0.1440445411234050561365333333333E-28,
    +0.9935374683141640465066666666666E-30,
    +0.6863462444358260053333333333333E-31 };
  static double dxrel = 0.0;
  static int nterms = 0;
  static double sqeps = 0.0;
  double value;
  double y;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( atnhcs, 27, 0.1 * r8_mach ( 3 ) );
    dxrel = sqrt ( r8_mach ( 4 ) );
    sqeps = sqrt ( 3.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= sqeps )
  {
    value = x;
  }
  else if ( y <= 0.5 )
  {
    value = x * ( 1.0 + r8_csevl ( 8.0 * x * x - 1.0, atnhcs, nterms ) );
  }
  else if ( y < 1.0 )
  {
    value = 0.5 * log ( ( 1.0 + x ) / ( 1.0 - x ) );
    if ( 1.0 - y < dxrel )
    {
      cerr << "\n";
      cerr << "R8_ATANH - Warning:\n";
      cerr << "  Answer lt half precision because |X| too near 1.\n";
    }
  }
  else
  {
    cerr << "\n";
    cerr << "R8_ATANH - Fatal error!\n";
    cerr << "  1 <= |X|.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

void r8_b0mp ( double x, double &ampl, double &theta )

//****************************************************************************80
//
//  Purpose:
//
//    R8_B0MP evaluates the modulus and phase for the Bessel J0 and Y0 functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double &AMPL, &THETA, the modulus and phase.
//
{
  static double bm0cs[37] = {
    +0.9211656246827742712573767730182E-01,
    -0.1050590997271905102480716371755E-02,
    +0.1470159840768759754056392850952E-04,
    -0.5058557606038554223347929327702E-06,
    +0.2787254538632444176630356137881E-07,
    -0.2062363611780914802618841018973E-08,
    +0.1870214313138879675138172596261E-09,
    -0.1969330971135636200241730777825E-10,
    +0.2325973793999275444012508818052E-11,
    -0.3009520344938250272851224734482E-12,
    +0.4194521333850669181471206768646E-13,
    -0.6219449312188445825973267429564E-14,
    +0.9718260411336068469601765885269E-15,
    -0.1588478585701075207366635966937E-15,
    +0.2700072193671308890086217324458E-16,
    -0.4750092365234008992477504786773E-17,
    +0.8615128162604370873191703746560E-18,
    -0.1605608686956144815745602703359E-18,
    +0.3066513987314482975188539801599E-19,
    -0.5987764223193956430696505617066E-20,
    +0.1192971253748248306489069841066E-20,
    -0.2420969142044805489484682581333E-21,
    +0.4996751760510616453371002879999E-22,
    -0.1047493639351158510095040511999E-22,
    +0.2227786843797468101048183466666E-23,
    -0.4801813239398162862370542933333E-24,
    +0.1047962723470959956476996266666E-24,
    -0.2313858165678615325101260800000E-25,
    +0.5164823088462674211635199999999E-26,
    -0.1164691191850065389525401599999E-26,
    +0.2651788486043319282958336000000E-27,
    -0.6092559503825728497691306666666E-28,
    +0.1411804686144259308038826666666E-28,
    -0.3298094961231737245750613333333E-29,
    +0.7763931143074065031714133333333E-30,
    -0.1841031343661458478421333333333E-30,
    +0.4395880138594310737100799999999E-31 };
  static double bm02cs[40] = {
    +0.9500415145228381369330861335560E-01,
    -0.3801864682365670991748081566851E-03,
    +0.2258339301031481192951829927224E-05,
    -0.3895725802372228764730621412605E-07,
    +0.1246886416512081697930990529725E-08,
    -0.6065949022102503779803835058387E-10,
    +0.4008461651421746991015275971045E-11,
    -0.3350998183398094218467298794574E-12,
    +0.3377119716517417367063264341996E-13,
    -0.3964585901635012700569356295823E-14,
    +0.5286111503883857217387939744735E-15,
    -0.7852519083450852313654640243493E-16,
    +0.1280300573386682201011634073449E-16,
    -0.2263996296391429776287099244884E-17,
    +0.4300496929656790388646410290477E-18,
    -0.8705749805132587079747535451455E-19,
    +0.1865862713962095141181442772050E-19,
    -0.4210482486093065457345086972301E-20,
    +0.9956676964228400991581627417842E-21,
    -0.2457357442805313359605921478547E-21,
    +0.6307692160762031568087353707059E-22,
    -0.1678773691440740142693331172388E-22,
    +0.4620259064673904433770878136087E-23,
    -0.1311782266860308732237693402496E-23,
    +0.3834087564116302827747922440276E-24,
    -0.1151459324077741271072613293576E-24,
    +0.3547210007523338523076971345213E-25,
    -0.1119218385815004646264355942176E-25,
    +0.3611879427629837831698404994257E-26,
    -0.1190687765913333150092641762463E-26,
    +0.4005094059403968131802476449536E-27,
    -0.1373169422452212390595193916017E-27,
    +0.4794199088742531585996491526437E-28,
    -0.1702965627624109584006994476452E-28,
    +0.6149512428936330071503575161324E-29,
    -0.2255766896581828349944300237242E-29,
    +0.8399707509294299486061658353200E-30,
    -0.3172997595562602355567423936152E-30,
    +0.1215205298881298554583333026514E-30,
    -0.4715852749754438693013210568045E-31 };
  static double bt02cs[39] = {
    -0.24548295213424597462050467249324,
    +0.12544121039084615780785331778299E-02,
    -0.31253950414871522854973446709571E-04,
    +0.14709778249940831164453426969314E-05,
    -0.99543488937950033643468850351158E-07,
    +0.85493166733203041247578711397751E-08,
    -0.86989759526554334557985512179192E-09,
    +0.10052099533559791084540101082153E-09,
    -0.12828230601708892903483623685544E-10,
    +0.17731700781805131705655750451023E-11,
    -0.26174574569485577488636284180925E-12,
    +0.40828351389972059621966481221103E-13,
    -0.66751668239742720054606749554261E-14,
    +0.11365761393071629448392469549951E-14,
    -0.20051189620647160250559266412117E-15,
    +0.36497978794766269635720591464106E-16,
    -0.68309637564582303169355843788800E-17,
    +0.13107583145670756620057104267946E-17,
    -0.25723363101850607778757130649599E-18,
    +0.51521657441863959925267780949333E-19,
    -0.10513017563758802637940741461333E-19,
    +0.21820381991194813847301084501333E-20,
    -0.46004701210362160577225905493333E-21,
    +0.98407006925466818520953651199999E-22,
    -0.21334038035728375844735986346666E-22,
    +0.46831036423973365296066286933333E-23,
    -0.10400213691985747236513382399999E-23,
    +0.23349105677301510051777740800000E-24,
    -0.52956825323318615788049749333333E-25,
    +0.12126341952959756829196287999999E-25,
    -0.28018897082289428760275626666666E-26,
    +0.65292678987012873342593706666666E-27,
    -0.15337980061873346427835733333333E-27,
    +0.36305884306364536682359466666666E-28,
    -0.86560755713629122479172266666666E-29,
    +0.20779909972536284571238399999999E-29,
    -0.50211170221417221674325333333333E-30,
    +0.12208360279441714184191999999999E-30,
    -0.29860056267039913454250666666666E-31 };
  static double bth0cs[44] = {
    -0.24901780862128936717709793789967,
    +0.48550299609623749241048615535485E-03,
    -0.54511837345017204950656273563505E-05,
    +0.13558673059405964054377445929903E-06,
    -0.55691398902227626227583218414920E-08,
    +0.32609031824994335304004205719468E-09,
    -0.24918807862461341125237903877993E-10,
    +0.23449377420882520554352413564891E-11,
    -0.26096534444310387762177574766136E-12,
    +0.33353140420097395105869955014923E-13,
    -0.47890000440572684646750770557409E-14,
    +0.75956178436192215972642568545248E-15,
    -0.13131556016891440382773397487633E-15,
    +0.24483618345240857495426820738355E-16,
    -0.48805729810618777683256761918331E-17,
    +0.10327285029786316149223756361204E-17,
    -0.23057633815057217157004744527025E-18,
    +0.54044443001892693993017108483765E-19,
    -0.13240695194366572724155032882385E-19,
    +0.33780795621371970203424792124722E-20,
    -0.89457629157111779003026926292299E-21,
    +0.24519906889219317090899908651405E-21,
    -0.69388422876866318680139933157657E-22,
    +0.20228278714890138392946303337791E-22,
    -0.60628500002335483105794195371764E-23,
    +0.18649748964037635381823788396270E-23,
    -0.58783732384849894560245036530867E-24,
    +0.18958591447999563485531179503513E-24,
    -0.62481979372258858959291620728565E-25,
    +0.21017901684551024686638633529074E-25,
    -0.72084300935209253690813933992446E-26,
    +0.25181363892474240867156405976746E-26,
    -0.89518042258785778806143945953643E-27,
    +0.32357237479762298533256235868587E-27,
    -0.11883010519855353657047144113796E-27,
    +0.44306286907358104820579231941731E-28,
    -0.16761009648834829495792010135681E-28,
    +0.64292946921207466972532393966088E-29,
    -0.24992261166978652421207213682763E-29,
    +0.98399794299521955672828260355318E-30,
    -0.39220375242408016397989131626158E-30,
    +0.15818107030056522138590618845692E-30,
    -0.64525506144890715944344098365426E-31,
    +0.26611111369199356137177018346367E-31 };
  double eta;
  static int nbm0 = 0;
  static int nbm02 = 0;
  static int nbt02 = 0;
  static int nbth0 = 0;
  static double pi4 = 0.785398163397448309615660845819876;
  static double xmax = 0.0;
  double z;

  if ( nbm0 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nbm0 = r8_inits ( bm0cs, 37, eta );
    nbt02 = r8_inits ( bt02cs, 39, eta );
    nbm02 = r8_inits ( bm02cs, 40, eta );
    nbth0 = r8_inits ( bth0cs, 44, eta );
    xmax = 1.0 / r8_mach ( 4 );
  }

  if ( x < 4.0 )
  {
    cerr << "\n";
    cerr << "R8_B0MP - Fatal error!\n";
    cerr << "  X < 4.\n";
    exit ( 1 );
  }
  else if ( x <= 8.0 )
  {
    z = ( 128.0 / x / x - 5.0 ) / 3.0;
    ampl = ( 0.75 + r8_csevl ( z, bm0cs, nbm0 ) ) / sqrt ( x );
    theta = x - pi4 + r8_csevl ( z, bt02cs, nbt02 ) / x;
  }
  else
  {
    z = 128.0 / x / x - 1.0;
    ampl = ( 0.75 + r8_csevl ( z, bm02cs, nbm02) ) / sqrt ( x );
    theta = x - pi4 + r8_csevl ( z, bth0cs, nbth0 ) / x;
  }
  return;
}
//****************************************************************************80

void r8_b1mp ( double x, double &ampl, double &theta )

//****************************************************************************80
//
//  Purpose:
//
//    R8_B1MP evaluates the modulus and phase for the Bessel J1 and Y1 functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double &AMPL, &THETA, the modulus and phase.
//
{
  static double bm12cs[40] = {
    +0.9807979156233050027272093546937E-01,
    +0.1150961189504685306175483484602E-02,
    -0.4312482164338205409889358097732E-05,
    +0.5951839610088816307813029801832E-07,
    -0.1704844019826909857400701586478E-08,
    +0.7798265413611109508658173827401E-10,
    -0.4958986126766415809491754951865E-11,
    +0.4038432416421141516838202265144E-12,
    -0.3993046163725175445765483846645E-13,
    +0.4619886183118966494313342432775E-14,
    -0.6089208019095383301345472619333E-15,
    +0.8960930916433876482157048041249E-16,
    -0.1449629423942023122916518918925E-16,
    +0.2546463158537776056165149648068E-17,
    -0.4809472874647836444259263718620E-18,
    +0.9687684668292599049087275839124E-19,
    -0.2067213372277966023245038117551E-19,
    +0.4646651559150384731802767809590E-20,
    -0.1094966128848334138241351328339E-20,
    +0.2693892797288682860905707612785E-21,
    -0.6894992910930374477818970026857E-22,
    +0.1830268262752062909890668554740E-22,
    -0.5025064246351916428156113553224E-23,
    +0.1423545194454806039631693634194E-23,
    -0.4152191203616450388068886769801E-24,
    +0.1244609201503979325882330076547E-24,
    -0.3827336370569304299431918661286E-25,
    +0.1205591357815617535374723981835E-25,
    -0.3884536246376488076431859361124E-26,
    +0.1278689528720409721904895283461E-26,
    -0.4295146689447946272061936915912E-27,
    +0.1470689117829070886456802707983E-27,
    -0.5128315665106073128180374017796E-28,
    +0.1819509585471169385481437373286E-28,
    -0.6563031314841980867618635050373E-29,
    +0.2404898976919960653198914875834E-29,
    -0.8945966744690612473234958242979E-30,
    +0.3376085160657231026637148978240E-30,
    -0.1291791454620656360913099916966E-30,
    +0.5008634462958810520684951501254E-31 };
  static double bm1cs[37] = {
    +0.1069845452618063014969985308538,
    +0.3274915039715964900729055143445E-02,
    -0.2987783266831698592030445777938E-04,
    +0.8331237177991974531393222669023E-06,
    -0.4112665690302007304896381725498E-07,
    +0.2855344228789215220719757663161E-08,
    -0.2485408305415623878060026596055E-09,
    +0.2543393338072582442742484397174E-10,
    -0.2941045772822967523489750827909E-11,
    +0.3743392025493903309265056153626E-12,
    -0.5149118293821167218720548243527E-13,
    +0.7552535949865143908034040764199E-14,
    -0.1169409706828846444166290622464E-14,
    +0.1896562449434791571721824605060E-15,
    -0.3201955368693286420664775316394E-16,
    +0.5599548399316204114484169905493E-17,
    -0.1010215894730432443119390444544E-17,
    +0.1873844985727562983302042719573E-18,
    -0.3563537470328580219274301439999E-19,
    +0.6931283819971238330422763519999E-20,
    -0.1376059453406500152251408930133E-20,
    +0.2783430784107080220599779327999E-21,
    -0.5727595364320561689348669439999E-22,
    +0.1197361445918892672535756799999E-22,
    -0.2539928509891871976641440426666E-23,
    +0.5461378289657295973069619199999E-24,
    -0.1189211341773320288986289493333E-24,
    +0.2620150977340081594957824000000E-25,
    -0.5836810774255685901920938666666E-26,
    +0.1313743500080595773423615999999E-26,
    -0.2985814622510380355332778666666E-27,
    +0.6848390471334604937625599999999E-28,
    -0.1584401568222476721192960000000E-28,
    +0.3695641006570938054301013333333E-29,
    -0.8687115921144668243012266666666E-30,
    +0.2057080846158763462929066666666E-30,
    -0.4905225761116225518523733333333E-31 };
  static double bt12cs[39] = {
    +0.73823860128742974662620839792764,
    -0.33361113174483906384470147681189E-02,
    +0.61463454888046964698514899420186E-04,
    -0.24024585161602374264977635469568E-05,
    +0.14663555577509746153210591997204E-06,
    -0.11841917305589180567005147504983E-07,
    +0.11574198963919197052125466303055E-08,
    -0.13001161129439187449366007794571E-09,
    +0.16245391141361731937742166273667E-10,
    -0.22089636821403188752155441770128E-11,
    +0.32180304258553177090474358653778E-12,
    -0.49653147932768480785552021135381E-13,
    +0.80438900432847825985558882639317E-14,
    -0.13589121310161291384694712682282E-14,
    +0.23810504397147214869676529605973E-15,
    -0.43081466363849106724471241420799E-16,
    +0.80202544032771002434993512550400E-17,
    -0.15316310642462311864230027468799E-17,
    +0.29928606352715568924073040554666E-18,
    -0.59709964658085443393815636650666E-19,
    +0.12140289669415185024160852650666E-19,
    -0.25115114696612948901006977706666E-20,
    +0.52790567170328744850738380799999E-21,
    -0.11260509227550498324361161386666E-21,
    +0.24348277359576326659663462400000E-22,
    -0.53317261236931800130038442666666E-23,
    +0.11813615059707121039205990399999E-23,
    -0.26465368283353523514856789333333E-24,
    +0.59903394041361503945577813333333E-25,
    -0.13690854630829503109136383999999E-25,
    +0.31576790154380228326413653333333E-26,
    -0.73457915082084356491400533333333E-27,
    +0.17228081480722747930705920000000E-27,
    -0.40716907961286507941068800000000E-28,
    +0.96934745136779622700373333333333E-29,
    -0.23237636337765716765354666666666E-29,
    +0.56074510673522029406890666666666E-30,
    -0.13616465391539005860522666666666E-30,
    +0.33263109233894654388906666666666E-31 };
  static double bth1cs[44] = {
    +0.74749957203587276055443483969695,
    -0.12400777144651711252545777541384E-02,
    +0.99252442404424527376641497689592E-05,
    -0.20303690737159711052419375375608E-06,
    +0.75359617705690885712184017583629E-08,
    -0.41661612715343550107630023856228E-09,
    +0.30701618070834890481245102091216E-10,
    -0.28178499637605213992324008883924E-11,
    +0.30790696739040295476028146821647E-12,
    -0.38803300262803434112787347554781E-13,
    +0.55096039608630904934561726208562E-14,
    -0.86590060768383779940103398953994E-15,
    +0.14856049141536749003423689060683E-15,
    -0.27519529815904085805371212125009E-16,
    +0.54550796090481089625036223640923E-17,
    -0.11486534501983642749543631027177E-17,
    +0.25535213377973900223199052533522E-18,
    -0.59621490197413450395768287907849E-19,
    +0.14556622902372718620288302005833E-19,
    -0.37022185422450538201579776019593E-20,
    +0.97763074125345357664168434517924E-21,
    -0.26726821639668488468723775393052E-21,
    +0.75453300384983271794038190655764E-22,
    -0.21947899919802744897892383371647E-22,
    +0.65648394623955262178906999817493E-23,
    -0.20155604298370207570784076869519E-23,
    +0.63417768556776143492144667185670E-24,
    -0.20419277885337895634813769955591E-24,
    +0.67191464220720567486658980018551E-25,
    -0.22569079110207573595709003687336E-25,
    +0.77297719892989706370926959871929E-26,
    -0.26967444512294640913211424080920E-26,
    +0.95749344518502698072295521933627E-27,
    -0.34569168448890113000175680827627E-27,
    +0.12681234817398436504211986238374E-27,
    -0.47232536630722639860464993713445E-28,
    +0.17850008478186376177858619796417E-28,
    -0.68404361004510395406215223566746E-29,
    +0.26566028671720419358293422672212E-29,
    -0.10450402527914452917714161484670E-29,
    +0.41618290825377144306861917197064E-30,
    -0.16771639203643714856501347882887E-30,
    +0.68361997776664389173535928028528E-31,
    -0.28172247861233641166739574622810E-31 };
  double eta;
  static int nbm1 = 0;
  static int nbm12 = 0;
  static int nbt12 = 0;
  static int nbth1 = 0;
  static double pi4 = 0.785398163397448309615660845819876;
  static double xmax = 0.0;
  double z;

  if ( nbm1 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nbm1 = r8_inits ( bm1cs, 37, eta );
    nbt12 = r8_inits ( bt12cs, 39, eta );
    nbm12 = r8_inits ( bm12cs, 40, eta );
    nbth1 = r8_inits ( bth1cs, 44, eta );
    xmax = 1.0 / r8_mach ( 4 );
  }

  if ( x < 4.0 )
  {
    cerr << "\n";
    cerr << "R8_B1MP - Fatal error!\n";
    cerr << "  X < 4.\n";
    exit ( 1 );
  }
  else if ( x <= 8.0 )
  {
    z = ( 128.0 / x / x - 5.0 ) / 3.0;
    ampl = ( 0.75 + r8_csevl ( z, bm1cs, nbm1 ) ) / sqrt ( x );
    theta = x - 3.0 * pi4 + r8_csevl ( z, bt12cs, nbt12 ) / x;
  }
  else
  {
    z = 128.0 / x / x - 1.0;
    ampl = ( 0.75 + r8_csevl ( z, bm12cs, nbm12 ) ) / sqrt ( x );
    theta = x - 3.0 * pi4 + r8_csevl ( z, bth1cs, nbth1 ) / x;
  }
  return;
}
//****************************************************************************80

double r8_besi0 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESI0 evaluates the Bessel function I of order 0 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESI0, the Bessel function I of order 0 of X.
//
{
  static double bi0cs[18] = {
    -0.7660547252839144951081894976243285E-01,
    +0.1927337953993808269952408750881196E+01,
    +0.2282644586920301338937029292330415,
    +0.1304891466707290428079334210691888E-01,
    +0.4344270900816487451378682681026107E-03,
    +0.9422657686001934663923171744118766E-05,
    +0.1434006289510691079962091878179957E-06,
    +0.1613849069661749069915419719994611E-08,
    +0.1396650044535669699495092708142522E-10,
    +0.9579451725505445344627523171893333E-13,
    +0.5333981859862502131015107744000000E-15,
    +0.2458716088437470774696785919999999E-17,
    +0.9535680890248770026944341333333333E-20,
    +0.3154382039721427336789333333333333E-22,
    +0.9004564101094637431466666666666666E-25,
    +0.2240647369123670016000000000000000E-27,
    +0.4903034603242837333333333333333333E-30,
    +0.9508172606122666666666666666666666E-33 };
  static int nti0 = 0;
  double value;
  static double xmax = 0.0;
  static double xsml = 0.0;
  double y;

  if ( nti0 == 0 )
  {
    nti0 = r8_inits ( bi0cs, 18, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( 8.0 * r8_mach ( 3 ) );
    xmax = log ( r8_mach ( 2 ) );
  }

  y = r8_abs ( x );

  if ( y <= xsml )
  {
    value = 1.0;
  }
  else if ( y <= 3.0 )
  {
    value = 2.75 + r8_csevl ( y * y / 4.5 - 1.0, bi0cs, nti0 );
  }
  else if ( y <= xmax )
  {
    value = exp ( y ) * r8_besi0e ( x );
  }
  else
  {
    cerr << "\n";
    cerr << "R8_BESI0 - Fatal error!\n";
    cerr << "  |X| too large.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

double r8_besi0e ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESI0E evaluates the exponentially scaled Bessel function I0(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESI0E, the exponentially scaled Bessel 
//    function I0(X).
//
{
  static double ai02cs[69] = {
    +0.5449041101410883160789609622680E-01,
    +0.3369116478255694089897856629799E-02,
    +0.6889758346916823984262639143011E-04,
    +0.2891370520834756482966924023232E-05,
    +0.2048918589469063741827605340931E-06,
    +0.2266668990498178064593277431361E-07,
    +0.3396232025708386345150843969523E-08,
    +0.4940602388224969589104824497835E-09,
    +0.1188914710784643834240845251963E-10,
    -0.3149916527963241364538648629619E-10,
    -0.1321581184044771311875407399267E-10,
    -0.1794178531506806117779435740269E-11,
    +0.7180124451383666233671064293469E-12,
    +0.3852778382742142701140898017776E-12,
    +0.1540086217521409826913258233397E-13,
    -0.4150569347287222086626899720156E-13,
    -0.9554846698828307648702144943125E-14,
    +0.3811680669352622420746055355118E-14,
    +0.1772560133056526383604932666758E-14,
    -0.3425485619677219134619247903282E-15,
    -0.2827623980516583484942055937594E-15,
    +0.3461222867697461093097062508134E-16,
    +0.4465621420296759999010420542843E-16,
    -0.4830504485944182071255254037954E-17,
    -0.7233180487874753954562272409245E-17,
    +0.9921475412173698598880460939810E-18,
    +0.1193650890845982085504399499242E-17,
    -0.2488709837150807235720544916602E-18,
    -0.1938426454160905928984697811326E-18,
    +0.6444656697373443868783019493949E-19,
    +0.2886051596289224326481713830734E-19,
    -0.1601954907174971807061671562007E-19,
    -0.3270815010592314720891935674859E-20,
    +0.3686932283826409181146007239393E-20,
    +0.1268297648030950153013595297109E-22,
    -0.7549825019377273907696366644101E-21,
    +0.1502133571377835349637127890534E-21,
    +0.1265195883509648534932087992483E-21,
    -0.6100998370083680708629408916002E-22,
    -0.1268809629260128264368720959242E-22,
    +0.1661016099890741457840384874905E-22,
    -0.1585194335765885579379705048814E-23,
    -0.3302645405968217800953817667556E-23,
    +0.1313580902839239781740396231174E-23,
    +0.3689040246671156793314256372804E-24,
    -0.4210141910461689149219782472499E-24,
    +0.4791954591082865780631714013730E-25,
    +0.8459470390221821795299717074124E-25,
    -0.4039800940872832493146079371810E-25,
    -0.6434714653650431347301008504695E-26,
    +0.1225743398875665990344647369905E-25,
    -0.2934391316025708923198798211754E-26,
    -0.1961311309194982926203712057289E-26,
    +0.1503520374822193424162299003098E-26,
    -0.9588720515744826552033863882069E-28,
    -0.3483339380817045486394411085114E-27,
    +0.1690903610263043673062449607256E-27,
    +0.1982866538735603043894001157188E-28,
    -0.5317498081491816214575830025284E-28,
    +0.1803306629888392946235014503901E-28,
    +0.6213093341454893175884053112422E-29,
    -0.7692189292772161863200728066730E-29,
    +0.1858252826111702542625560165963E-29,
    +0.1237585142281395724899271545541E-29,
    -0.1102259120409223803217794787792E-29,
    +0.1886287118039704490077874479431E-30,
    +0.2160196872243658913149031414060E-30,
    -0.1605454124919743200584465949655E-30,
    +0.1965352984594290603938848073318E-31 };
  static double ai0cs[46] = {
    +0.7575994494023795942729872037438E-01,
    +0.7591380810823345507292978733204E-02,
    +0.4153131338923750501863197491382E-03,
    +0.1070076463439073073582429702170E-04,
    -0.7901179979212894660750319485730E-05,
    -0.7826143501438752269788989806909E-06,
    +0.2783849942948870806381185389857E-06,
    +0.8252472600612027191966829133198E-08,
    -0.1204463945520199179054960891103E-07,
    +0.1559648598506076443612287527928E-08,
    +0.2292556367103316543477254802857E-09,
    -0.1191622884279064603677774234478E-09,
    +0.1757854916032409830218331247743E-10,
    +0.1128224463218900517144411356824E-11,
    -0.1146848625927298877729633876982E-11,
    +0.2715592054803662872643651921606E-12,
    -0.2415874666562687838442475720281E-13,
    -0.6084469888255125064606099639224E-14,
    +0.3145705077175477293708360267303E-14,
    -0.7172212924871187717962175059176E-15,
    +0.7874493403454103396083909603327E-16,
    +0.1004802753009462402345244571839E-16,
    -0.7566895365350534853428435888810E-17,
    +0.2150380106876119887812051287845E-17,
    -0.3754858341830874429151584452608E-18,
    +0.2354065842226992576900757105322E-19,
    +0.1114667612047928530226373355110E-19,
    -0.5398891884396990378696779322709E-20,
    +0.1439598792240752677042858404522E-20,
    -0.2591916360111093406460818401962E-21,
    +0.2238133183998583907434092298240E-22,
    +0.5250672575364771172772216831999E-23,
    -0.3249904138533230784173432285866E-23,
    +0.9924214103205037927857284710400E-24,
    -0.2164992254244669523146554299733E-24,
    +0.3233609471943594083973332991999E-25,
    -0.1184620207396742489824733866666E-26,
    -0.1281671853950498650548338687999E-26,
    +0.5827015182279390511605568853333E-27,
    -0.1668222326026109719364501503999E-27,
    +0.3625309510541569975700684800000E-28,
    -0.5733627999055713589945958399999E-29,
    +0.3736796722063098229642581333333E-30,
    +0.1602073983156851963365512533333E-30,
    -0.8700424864057229884522495999999E-31,
    +0.2741320937937481145603413333333E-31 };
  static double bi0cs[18] = {
    -0.7660547252839144951081894976243285E-01,
    +0.1927337953993808269952408750881196E+01,
    +0.2282644586920301338937029292330415,
    +0.1304891466707290428079334210691888E-01,
    +0.4344270900816487451378682681026107E-03,
    +0.9422657686001934663923171744118766E-05,
    +0.1434006289510691079962091878179957E-06,
    +0.1613849069661749069915419719994611E-08,
    +0.1396650044535669699495092708142522E-10,
    +0.9579451725505445344627523171893333E-13,
    +0.5333981859862502131015107744000000E-15,
    +0.2458716088437470774696785919999999E-17,
    +0.9535680890248770026944341333333333E-20,
    +0.3154382039721427336789333333333333E-22,
    +0.9004564101094637431466666666666666E-25,
    +0.2240647369123670016000000000000000E-27,
    +0.4903034603242837333333333333333333E-30,
    +0.9508172606122666666666666666666666E-33 };
  double eta;
  static int ntai02 = 0;
  static int ntai0 = 0;
  static int nti0 = 0;
  double value;
  static double xsml = 0.0;
  double y;

  if ( nti0 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nti0 = r8_inits ( bi0cs, 18, eta );
    ntai0 = r8_inits ( ai0cs, 46, eta );
    ntai02 = r8_inits ( ai02cs, 69, eta );
    xsml = sqrt ( 8.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= xsml )
  {
    value = 1.0;
  }
  else if ( y <= 3.0 )
  {
    value = exp ( - y ) * ( 2.75 
      + r8_csevl ( y * y / 4.5 - 1.0, bi0cs, nti0 ) );
  }
  else if ( y <= 8.0 )
  {
    value = ( 0.375
      + r8_csevl ( ( 48.0 / y - 11.0 ) / 5.0, ai0cs, ntai0 ) ) / sqrt ( y );
  }
  else
  {
    value = ( 0.375 
      + r8_csevl ( 16.0 / y - 1.0, ai02cs, ntai02 ) )  / sqrt ( y );
  }
  return value;
}
//****************************************************************************80

double r8_besi1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESI1 evaluates the Bessel function I of order 1 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESI1, the Bessel function I of order 1 of X.
//
{
  static double bi1cs[17] = {
    -0.19717132610998597316138503218149E-02,
    +0.40734887667546480608155393652014,
    +0.34838994299959455866245037783787E-01,
    +0.15453945563001236038598401058489E-02,
    +0.41888521098377784129458832004120E-04,
    +0.76490267648362114741959703966069E-06,
    +0.10042493924741178689179808037238E-07,
    +0.99322077919238106481371298054863E-10,
    +0.76638017918447637275200171681349E-12,
    +0.47414189238167394980388091948160E-14,
    +0.24041144040745181799863172032000E-16,
    +0.10171505007093713649121100799999E-18,
    +0.36450935657866949458491733333333E-21,
    +0.11205749502562039344810666666666E-23,
    +0.29875441934468088832000000000000E-26,
    +0.69732310939194709333333333333333E-29,
    +0.14367948220620800000000000000000E-31 };
  static int nti1 = 0;
  double  value;
  static double xmax = 0.0;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( nti1 == 0 )
  {
    nti1 = r8_inits ( bi1cs, 17, 0.1 * r8_mach ( 3 ) );
    xmin = 2.0 * r8_mach ( 1 );
    xsml = sqrt ( 8.0 * r8_mach ( 3 ) );
    xmax = log ( r8_mach ( 2 ) );
  }

  y = r8_abs ( x );

  if ( y <= xmin )
  {
    value = 0.0;
  }
  else if ( y <= xsml )
  {
    value = 0.5 * x;
  }
  else if ( y <= 3.0 )
  {
    value = x * ( 0.875 + r8_csevl ( y * y / 4.5 - 1.0, bi1cs, nti1 ) );
  }
  else if ( y <= xmax )
  {
    value = exp ( y ) * r8_besi1e ( x );
  }
  else
  {
    cerr << "\n";
    cerr << "R8_BESI1 - Fatal error!\n";
    cerr << "  Result overflows.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

double r8_besi1e ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESI1E evaluates the exponentially scaled Bessel function I1(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESI1E, the exponentially scaled Bessel 
//    function I1(X).
//
{
  static double ai12cs[69] = {
    +0.2857623501828012047449845948469E-01,
    -0.9761097491361468407765164457302E-02,
    -0.1105889387626237162912569212775E-03,
    -0.3882564808877690393456544776274E-05,
    -0.2512236237870208925294520022121E-06,
    -0.2631468846889519506837052365232E-07,
    -0.3835380385964237022045006787968E-08,
    -0.5589743462196583806868112522229E-09,
    -0.1897495812350541234498925033238E-10,
    +0.3252603583015488238555080679949E-10,
    +0.1412580743661378133163366332846E-10,
    +0.2035628544147089507224526136840E-11,
    -0.7198551776245908512092589890446E-12,
    -0.4083551111092197318228499639691E-12,
    -0.2101541842772664313019845727462E-13,
    +0.4272440016711951354297788336997E-13,
    +0.1042027698412880276417414499948E-13,
    -0.3814403072437007804767072535396E-14,
    -0.1880354775510782448512734533963E-14,
    +0.3308202310920928282731903352405E-15,
    +0.2962628997645950139068546542052E-15,
    -0.3209525921993423958778373532887E-16,
    -0.4650305368489358325571282818979E-16,
    +0.4414348323071707949946113759641E-17,
    +0.7517296310842104805425458080295E-17,
    -0.9314178867326883375684847845157E-18,
    -0.1242193275194890956116784488697E-17,
    +0.2414276719454848469005153902176E-18,
    +0.2026944384053285178971922860692E-18,
    -0.6394267188269097787043919886811E-19,
    -0.3049812452373095896084884503571E-19,
    +0.1612841851651480225134622307691E-19,
    +0.3560913964309925054510270904620E-20,
    -0.3752017947936439079666828003246E-20,
    -0.5787037427074799345951982310741E-22,
    +0.7759997511648161961982369632092E-21,
    -0.1452790897202233394064459874085E-21,
    -0.1318225286739036702121922753374E-21,
    +0.6116654862903070701879991331717E-22,
    +0.1376279762427126427730243383634E-22,
    -0.1690837689959347884919839382306E-22,
    +0.1430596088595433153987201085385E-23,
    +0.3409557828090594020405367729902E-23,
    -0.1309457666270760227845738726424E-23,
    -0.3940706411240257436093521417557E-24,
    +0.4277137426980876580806166797352E-24,
    -0.4424634830982606881900283123029E-25,
    -0.8734113196230714972115309788747E-25,
    +0.4045401335683533392143404142428E-25,
    +0.7067100658094689465651607717806E-26,
    -0.1249463344565105223002864518605E-25,
    +0.2867392244403437032979483391426E-26,
    +0.2044292892504292670281779574210E-26,
    -0.1518636633820462568371346802911E-26,
    +0.8110181098187575886132279107037E-28,
    +0.3580379354773586091127173703270E-27,
    -0.1692929018927902509593057175448E-27,
    -0.2222902499702427639067758527774E-28,
    +0.5424535127145969655048600401128E-28,
    -0.1787068401578018688764912993304E-28,
    -0.6565479068722814938823929437880E-29,
    +0.7807013165061145280922067706839E-29,
    -0.1816595260668979717379333152221E-29,
    -0.1287704952660084820376875598959E-29,
    +0.1114548172988164547413709273694E-29,
    -0.1808343145039336939159368876687E-30,
    -0.2231677718203771952232448228939E-30,
    +0.1619029596080341510617909803614E-30,
    -0.1834079908804941413901308439210E-31 };
  static double ai1cs[46] = {
    -0.2846744181881478674100372468307E-01,
    -0.1922953231443220651044448774979E-01,
    -0.6115185857943788982256249917785E-03,
    -0.2069971253350227708882823777979E-04,
    +0.8585619145810725565536944673138E-05,
    +0.1049498246711590862517453997860E-05,
    -0.2918338918447902202093432326697E-06,
    -0.1559378146631739000160680969077E-07,
    +0.1318012367144944705525302873909E-07,
    -0.1448423418183078317639134467815E-08,
    -0.2908512243993142094825040993010E-09,
    +0.1266388917875382387311159690403E-09,
    -0.1664947772919220670624178398580E-10,
    -0.1666653644609432976095937154999E-11,
    +0.1242602414290768265232168472017E-11,
    -0.2731549379672432397251461428633E-12,
    +0.2023947881645803780700262688981E-13,
    +0.7307950018116883636198698126123E-14,
    -0.3332905634404674943813778617133E-14,
    +0.7175346558512953743542254665670E-15,
    -0.6982530324796256355850629223656E-16,
    -0.1299944201562760760060446080587E-16,
    +0.8120942864242798892054678342860E-17,
    -0.2194016207410736898156266643783E-17,
    +0.3630516170029654848279860932334E-18,
    -0.1695139772439104166306866790399E-19,
    -0.1288184829897907807116882538222E-19,
    +0.5694428604967052780109991073109E-20,
    -0.1459597009090480056545509900287E-20,
    +0.2514546010675717314084691334485E-21,
    -0.1844758883139124818160400029013E-22,
    -0.6339760596227948641928609791999E-23,
    +0.3461441102031011111108146626560E-23,
    -0.1017062335371393547596541023573E-23,
    +0.2149877147090431445962500778666E-24,
    -0.3045252425238676401746206173866E-25,
    +0.5238082144721285982177634986666E-27,
    +0.1443583107089382446416789503999E-26,
    -0.6121302074890042733200670719999E-27,
    +0.1700011117467818418349189802666E-27,
    -0.3596589107984244158535215786666E-28,
    +0.5448178578948418576650513066666E-29,
    -0.2731831789689084989162564266666E-30,
    -0.1858905021708600715771903999999E-30,
    +0.9212682974513933441127765333333E-31,
    -0.2813835155653561106370833066666E-31 };
  static double bi1cs[17] = {
    -0.19717132610998597316138503218149E-02,
    +0.40734887667546480608155393652014,
    +0.34838994299959455866245037783787E-01,
    +0.15453945563001236038598401058489E-02,
    +0.41888521098377784129458832004120E-04,
    +0.76490267648362114741959703966069E-06,
    +0.10042493924741178689179808037238E-07,
    +0.99322077919238106481371298054863E-10,
    +0.76638017918447637275200171681349E-12,
    +0.47414189238167394980388091948160E-14,
    +0.24041144040745181799863172032000E-16,
    +0.10171505007093713649121100799999E-18,
    +0.36450935657866949458491733333333E-21,
    +0.11205749502562039344810666666666E-23,
    +0.29875441934468088832000000000000E-26,
    +0.69732310939194709333333333333333E-29,
    +0.14367948220620800000000000000000E-31 };
  double eta;
  static int ntai1 = 0;
  static int ntai12 = 0;
  static int nti1 = 0;
  double value;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( nti1 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nti1 = r8_inits ( bi1cs, 17, eta );
    ntai1 = r8_inits ( ai1cs, 46, eta );
    ntai12 = r8_inits ( ai12cs, 69, eta );
    xmin = 2.0 * r8_mach ( 1 );
    xsml = sqrt ( 8.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= xmin )
  {
    value = 0.0;
  }
  else if ( y <= xsml )
  {
    value = 0.5 * x * exp ( - y );
  }
  else if ( y <= 3.0 )
  {
    value = x * ( 0.875 + r8_csevl ( y * y / 4.5 - 1.0, bi1cs, nti1 ) )
      * exp ( - y );
  }
  else if ( y <= 8.0 )
  {
    value = ( 0.375 + r8_csevl ( ( 48.0 / y - 11.0) / 5.0, 
      ai1cs, ntai1 ) ) / sqrt ( y );
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  else
  {
    value = ( 0.375 + r8_csevl ( 16.0 / y - 1.0, ai12cs, ntai12 ) ) / sqrt ( y );
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

double r8_besj0 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESJ0 evaluates the Bessel function J of order 0 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESJ0, the Bessel function J of order 0 of X.
//
{
  double ampl;
  static double bj0cs[19] = {
    +0.10025416196893913701073127264074,
    -0.66522300776440513177678757831124,
    +0.24898370349828131370460468726680,
    -0.33252723170035769653884341503854E-01,
    +0.23114179304694015462904924117729E-02,
    -0.99112774199508092339048519336549E-04,
    +0.28916708643998808884733903747078E-05,
    -0.61210858663032635057818407481516E-07,
    +0.98386507938567841324768748636415E-09,
    -0.12423551597301765145515897006836E-10,
    +0.12654336302559045797915827210363E-12,
    -0.10619456495287244546914817512959E-14,
    +0.74706210758024567437098915584000E-17,
    -0.44697032274412780547627007999999E-19,
    +0.23024281584337436200523093333333E-21,
    -0.10319144794166698148522666666666E-23,
    +0.40608178274873322700800000000000E-26,
    -0.14143836005240913919999999999999E-28,
    +0.43910905496698880000000000000000E-31 };
  static int ntj0 = 0;
  double theta;
  double value;
  static double xsml = 0.0;
  double y;

  if ( ntj0 == 0 )
  {
    ntj0 = r8_inits ( bj0cs, 19, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= xsml )
  {
    value = 1.0;
  }
  else if ( y <= 4.0 )
  {
    value = r8_csevl ( 0.125 * y * y - 1.0, bj0cs, ntj0 );
  }
  else
  {
    r8_b0mp ( y, ampl, theta );
    value = ampl * cos ( theta );
  }
  return value;
}
//****************************************************************************80

double r8_besj1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESJ1 evaluates the Bessel function J of order 1 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESJ1, the Bessel function J of order 1 of X.
//
{
  double ampl;
  static double bj1cs[19] = {
    -0.117261415133327865606240574524003,
    -0.253615218307906395623030884554698,
    +0.501270809844695685053656363203743E-01,
    -0.463151480962508191842619728789772E-02,
    +0.247996229415914024539124064592364E-03,
    -0.867894868627882584521246435176416E-05,
    +0.214293917143793691502766250991292E-06,
    -0.393609307918317979229322764073061E-08,
    +0.559118231794688004018248059864032E-10,
    -0.632761640466139302477695274014880E-12,
    +0.584099161085724700326945563268266E-14,
    -0.448253381870125819039135059199999E-16,
    +0.290538449262502466306018688000000E-18,
    -0.161173219784144165412118186666666E-20,
    +0.773947881939274637298346666666666E-23,
    -0.324869378211199841143466666666666E-25,
    +0.120223767722741022720000000000000E-27,
    -0.395201221265134933333333333333333E-30,
    +0.116167808226645333333333333333333E-32 };
  static int ntj1 = 0;
  double theta;
  double value;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ntj1 == 0 )
  {
    ntj1 = r8_inits ( bj1cs, 19, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
    xmin = 2.0 * r8_mach ( 1 );
  }

  y = r8_abs ( x );

  if ( y <= xmin )
  {
    value = 0.0;
  }
  else if ( y <= xsml )
  {
    value = 0.5 * x;
  }
  else if ( y <= 4.0 )
  {
    value = x * ( 0.25 + r8_csevl ( 0.125 * y * y - 1.0, bj1cs, ntj1 ) );
  }
  else
  {
    r8_b1mp ( y, ampl, theta );
    if ( x < 0.0 )
    {
      value = - ampl * cos ( theta );
    }
    else
    {
      value = + ampl * cos ( theta );
    }
  }
  return value;
}
//****************************************************************************80

double r8_besk0 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESK0 evaluates the Bessel function K of order 0 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESK0, the Bessel function K of order 0 of X.
//
{
  static double bk0cs[16] = {
    -0.353273932339027687201140060063153E-01,
    +0.344289899924628486886344927529213,
    +0.359799365153615016265721303687231E-01,
    +0.126461541144692592338479508673447E-02,
    +0.228621210311945178608269830297585E-04,
    +0.253479107902614945730790013428354E-06,
    +0.190451637722020885897214059381366E-08,
    +0.103496952576336245851008317853089E-10,
    +0.425981614279108257652445327170133E-13,
    +0.137446543588075089694238325440000E-15,
    +0.357089652850837359099688597333333E-18,
    +0.763164366011643737667498666666666E-21,
    +0.136542498844078185908053333333333E-23,
    +0.207527526690666808319999999999999E-26,
    +0.271281421807298560000000000000000E-29,
    +0.308259388791466666666666666666666E-32 };
  static int ntk0 = 0;
  double value;
  static double xmax = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ntk0 == 0 )
  {
    ntk0 = r8_inits (bk0cs, 16, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
    xmax = - log ( r8_mach ( 1 ) );
    xmax = xmax - 0.5 * xmax * log ( xmax ) / ( xmax + 0.5 );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESK0 = Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0;
    value = - log ( 0.5 * x ) * r8_besi0 ( x )  
      - 0.25 + r8_csevl ( 0.5 * y - 1.0, bk0cs, ntk0 );
  }
  else if ( x <= 2.0 )
  {
    y = x * x;
    value = - log ( 0.5 * x ) * r8_besi0 ( x )  
      - 0.25 + r8_csevl ( 0.5 * y - 1.0, bk0cs, ntk0 );
  }
  else if ( x <= xmax )
  {
    value = exp ( - x ) * r8_besk0e ( x );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double r8_besk0e ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESK0E evaluates the exponentially scaled Bessel function K0(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESK0E, the exponentially scaled Bessel 
//    function K0(X).
//
{
  static double ak02cs[33] = {
    -0.1201869826307592239839346212452E-01,
    -0.9174852691025695310652561075713E-02,
    +0.1444550931775005821048843878057E-03,
    -0.4013614175435709728671021077879E-05,
    +0.1567831810852310672590348990333E-06,
    -0.7770110438521737710315799754460E-08,
    +0.4611182576179717882533130529586E-09,
    -0.3158592997860565770526665803309E-10,
    +0.2435018039365041127835887814329E-11,
    -0.2074331387398347897709853373506E-12,
    +0.1925787280589917084742736504693E-13,
    -0.1927554805838956103600347182218E-14,
    +0.2062198029197818278285237869644E-15,
    -0.2341685117579242402603640195071E-16,
    +0.2805902810643042246815178828458E-17,
    -0.3530507631161807945815482463573E-18,
    +0.4645295422935108267424216337066E-19,
    -0.6368625941344266473922053461333E-20,
    +0.9069521310986515567622348800000E-21,
    -0.1337974785423690739845005311999E-21,
    +0.2039836021859952315522088960000E-22,
    -0.3207027481367840500060869973333E-23,
    +0.5189744413662309963626359466666E-24,
    -0.8629501497540572192964607999999E-25,
    +0.1472161183102559855208038400000E-25,
    -0.2573069023867011283812351999999E-26,
    +0.4601774086643516587376640000000E-27,
    -0.8411555324201093737130666666666E-28,
    +0.1569806306635368939301546666666E-28,
    -0.2988226453005757788979199999999E-29,
    +0.5796831375216836520618666666666E-30,
    -0.1145035994347681332155733333333E-30,
    +0.2301266594249682802005333333333E-31 };
  static double ak0cs[38] = {
    -0.7643947903327941424082978270088E-01,
    -0.2235652605699819052023095550791E-01,
    +0.7734181154693858235300618174047E-03,
    -0.4281006688886099464452146435416E-04,
    +0.3081700173862974743650014826660E-05,
    -0.2639367222009664974067448892723E-06,
    +0.2563713036403469206294088265742E-07,
    -0.2742705549900201263857211915244E-08,
    +0.3169429658097499592080832873403E-09,
    -0.3902353286962184141601065717962E-10,
    +0.5068040698188575402050092127286E-11,
    -0.6889574741007870679541713557984E-12,
    +0.9744978497825917691388201336831E-13,
    -0.1427332841884548505389855340122E-13,
    +0.2156412571021463039558062976527E-14,
    -0.3349654255149562772188782058530E-15,
    +0.5335260216952911692145280392601E-16,
    -0.8693669980890753807639622378837E-17,
    +0.1446404347862212227887763442346E-17,
    -0.2452889825500129682404678751573E-18,
    +0.4233754526232171572821706342400E-19,
    -0.7427946526454464195695341294933E-20,
    +0.1323150529392666866277967462400E-20,
    -0.2390587164739649451335981465599E-21,
    +0.4376827585923226140165712554666E-22,
    -0.8113700607345118059339011413333E-23,
    +0.1521819913832172958310378154666E-23,
    -0.2886041941483397770235958613333E-24,
    +0.5530620667054717979992610133333E-25,
    -0.1070377329249898728591633066666E-25,
    +0.2091086893142384300296328533333E-26,
    -0.4121713723646203827410261333333E-27,
    +0.8193483971121307640135680000000E-28,
    -0.1642000275459297726780757333333E-28,
    +0.3316143281480227195890346666666E-29,
    -0.6746863644145295941085866666666E-30,
    +0.1382429146318424677635413333333E-30,
    -0.2851874167359832570811733333333E-31 };
  static double bk0cs[16] = {
    -0.353273932339027687201140060063153E-01,
    +0.344289899924628486886344927529213,
    +0.359799365153615016265721303687231E-01,
    +0.126461541144692592338479508673447E-02,
    +0.228621210311945178608269830297585E-04,
    +0.253479107902614945730790013428354E-06,
    +0.190451637722020885897214059381366E-08,
    +0.103496952576336245851008317853089E-10,
    +0.425981614279108257652445327170133E-13,
    +0.137446543588075089694238325440000E-15,
    +0.357089652850837359099688597333333E-18,
    +0.763164366011643737667498666666666E-21,
    +0.136542498844078185908053333333333E-23,
    +0.207527526690666808319999999999999E-26,
    +0.271281421807298560000000000000000E-29,
    +0.308259388791466666666666666666666E-32 };
  double eta;
  static int ntak0 = 0;
  static int ntak02 = 0;
  static int ntk0 = 0;
  double value;
  static double xsml = 0.0;
  double y;

  if ( ntk0 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    ntk0 = r8_inits ( bk0cs, 16, eta );
    ntak0 = r8_inits ( ak0cs, 38, eta );
    ntak02 = r8_inits ( ak02cs, 33, eta );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESK0E = Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0;
    value = exp ( x ) * ( - log ( 0.5 * x ) * r8_besi0 ( x ) - 0.25 
      + r8_csevl ( 0.5 * y - 1.0, bk0cs, ntk0 ) );
  }
  else if ( x <= 2.0 )
  {
    y = x * x;
    value = exp ( x ) * ( - log ( 0.5 * x ) * r8_besi0 ( x ) - 0.25 
      + r8_csevl ( 0.5 * y - 1.0, bk0cs, ntk0 ) );
  }
  else if ( x <= 8.0 )
  {
    value = ( 1.25 + r8_csevl ( ( 16.0 / x - 5.0 ) / 3.0, ak0cs, 
      ntak0 ) ) / sqrt ( x );
  }
  else
  {
    value = ( 1.25 + r8_csevl ( 16.0 / x - 1.0, ak02cs, ntak02 ) ) / sqrt ( x );
  } 
  return value;
}
//****************************************************************************80

double r8_besk1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESK1 evaluates the Bessel function K of order 1 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESK1, the Bessel function K of order 1 of X.
//
{
  static double bk1cs[16] = {
    +0.25300227338947770532531120868533E-01,
    -0.35315596077654487566723831691801,
    -0.12261118082265714823479067930042,
    -0.69757238596398643501812920296083E-02,
    -0.17302889575130520630176507368979E-03,
    -0.24334061415659682349600735030164E-05,
    -0.22133876307347258558315252545126E-07,
    -0.14114883926335277610958330212608E-09,
    -0.66669016941993290060853751264373E-12,
    -0.24274498505193659339263196864853E-14,
    -0.70238634793862875971783797120000E-17,
    -0.16543275155100994675491029333333E-19,
    -0.32338347459944491991893333333333E-22,
    -0.53312750529265274999466666666666E-25,
    -0.75130407162157226666666666666666E-28,
    -0.91550857176541866666666666666666E-31 };
  static int ntk1 = 0;
  double value;
  static double xmax = 0.0;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ntk1 == 0 )
  {
    ntk1 = r8_inits ( bk1cs, 16, 0.1 * r8_mach ( 3 ) );
    xmin = exp ( r8_max ( log ( r8_mach ( 1 ) ), 
      - log ( r8_mach ( 2 ) ) ) + 0.01 );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
    xmax = - log ( r8_mach ( 1 ) );
    xmax = xmax - 0.5 * xmax * log ( xmax ) 
      / ( xmax + 0.5 ) - 0.01;
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESK1 = Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0;
    value = log ( 0.5 * x ) * r8_besi1 ( x ) + ( 0.75 
      + r8_csevl ( 0.5 * y - 1.0, bk1cs, ntk1 ) ) / x;
  }
  else if ( x <= 2.0 )
  {
    y = x * x;
    value = log ( 0.5 * x ) * r8_besi1 ( x ) + ( 0.75 
      + r8_csevl ( 0.5 * y - 1.0, bk1cs, ntk1 ) ) / x;
  }
  else if ( x <= xmax )
  {
    value = exp ( - x ) * r8_besk1e ( x );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double r8_besk1e ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESK1E evaluates the exponentially scaled Bessel function K1(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESK1E, the exponentially scaled Bessel 
//    function K1(X).
//
{
  static double ak12cs[33] = {
    +0.6379308343739001036600488534102E-01,
    +0.2832887813049720935835030284708E-01,
    -0.2475370673905250345414545566732E-03,
    +0.5771972451607248820470976625763E-05,
    -0.2068939219536548302745533196552E-06,
    +0.9739983441381804180309213097887E-08,
    -0.5585336140380624984688895511129E-09,
    +0.3732996634046185240221212854731E-10,
    -0.2825051961023225445135065754928E-11,
    +0.2372019002484144173643496955486E-12,
    -0.2176677387991753979268301667938E-13,
    +0.2157914161616032453939562689706E-14,
    -0.2290196930718269275991551338154E-15,
    +0.2582885729823274961919939565226E-16,
    -0.3076752641268463187621098173440E-17,
    +0.3851487721280491597094896844799E-18,
    -0.5044794897641528977117282508800E-19,
    +0.6888673850418544237018292223999E-20,
    -0.9775041541950118303002132480000E-21,
    +0.1437416218523836461001659733333E-21,
    -0.2185059497344347373499733333333E-22,
    +0.3426245621809220631645388800000E-23,
    -0.5531064394246408232501248000000E-24,
    +0.9176601505685995403782826666666E-25,
    -0.1562287203618024911448746666666E-25,
    +0.2725419375484333132349439999999E-26,
    -0.4865674910074827992378026666666E-27,
    +0.8879388552723502587357866666666E-28,
    -0.1654585918039257548936533333333E-28,
    +0.3145111321357848674303999999999E-29,
    -0.6092998312193127612416000000000E-30,
    +0.1202021939369815834623999999999E-30,
    -0.2412930801459408841386666666666E-31 };
  static double ak1cs[38] = {
    +0.27443134069738829695257666227266,
    +0.75719899531993678170892378149290E-01,
    -0.14410515564754061229853116175625E-02,
    +0.66501169551257479394251385477036E-04,
    -0.43699847095201407660580845089167E-05,
    +0.35402774997630526799417139008534E-06,
    -0.33111637792932920208982688245704E-07,
    +0.34459775819010534532311499770992E-08,
    -0.38989323474754271048981937492758E-09,
    +0.47208197504658356400947449339005E-10,
    -0.60478356628753562345373591562890E-11,
    +0.81284948748658747888193837985663E-12,
    -0.11386945747147891428923915951042E-12,
    +0.16540358408462282325972948205090E-13,
    -0.24809025677068848221516010440533E-14,
    +0.38292378907024096948429227299157E-15,
    -0.60647341040012418187768210377386E-16,
    +0.98324256232648616038194004650666E-17,
    -0.16284168738284380035666620115626E-17,
    +0.27501536496752623718284120337066E-18,
    -0.47289666463953250924281069568000E-19,
    +0.82681500028109932722392050346666E-20,
    -0.14681405136624956337193964885333E-20,
    +0.26447639269208245978085894826666E-21,
    -0.48290157564856387897969868800000E-22,
    +0.89293020743610130180656332799999E-23,
    -0.16708397168972517176997751466666E-23,
    +0.31616456034040694931368618666666E-24,
    -0.60462055312274989106506410666666E-25,
    +0.11678798942042732700718421333333E-25,
    -0.22773741582653996232867840000000E-26,
    +0.44811097300773675795305813333333E-27,
    -0.88932884769020194062336000000000E-28,
    +0.17794680018850275131392000000000E-28,
    -0.35884555967329095821994666666666E-29,
    +0.72906290492694257991679999999999E-30,
    -0.14918449845546227073024000000000E-30,
    +0.30736573872934276300799999999999E-31 };
  static double bk1cs[16] = {
    +0.25300227338947770532531120868533E-01,
    -0.35315596077654487566723831691801,
    -0.12261118082265714823479067930042,
    -0.69757238596398643501812920296083E-02,
    -0.17302889575130520630176507368979E-03,
    -0.24334061415659682349600735030164E-05,
    -0.22133876307347258558315252545126E-07,
    -0.14114883926335277610958330212608E-09,
    -0.66669016941993290060853751264373E-12,
    -0.24274498505193659339263196864853E-14,
    -0.70238634793862875971783797120000E-17,
    -0.16543275155100994675491029333333E-19,
    -0.32338347459944491991893333333333E-22,
    -0.53312750529265274999466666666666E-25,
    -0.75130407162157226666666666666666E-28,
    -0.91550857176541866666666666666666E-31 };
  double eta;
  static int ntak1 = 0;
  static int ntak12 = 0;
  static int ntk1 = 0;
  double value;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ntk1 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    ntk1 = r8_inits ( bk1cs, 16, eta );
    ntak1 = r8_inits ( ak1cs, 38, eta );
    ntak12 = r8_inits ( ak12cs, 33, eta );
    xmin = exp ( r8_max ( log ( r8_mach ( 1 ) ), 
      - log ( r8_mach ( 2 ) ) ) + 0.01 );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESK1E = Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0;
    value = exp ( x ) * ( log ( 0.5 * x ) * r8_besi1 ( x )
      + ( 0.75 + r8_csevl ( 0.5 * y - 1.0, bk1cs, ntk1 ) ) / x );
  }
  else if ( x <= 2.0 )
  {
    y = x * x;
    value = exp ( x ) * ( log ( 0.5 * x ) * r8_besi1 ( x )
      + ( 0.75 + r8_csevl ( 0.5 * y - 1.0, bk1cs, ntk1 ) ) / x );
  }
  else if ( x <= 8.0 )
  {
    value = ( 1.25 
      + r8_csevl ( ( 16.0 / x - 5.0 ) / 3.0, ak1cs, ntak1 ) ) / sqrt ( x );
  }
  else
  {
    value = ( 1.25 +
      r8_csevl ( 16.0 / x - 1.0, ak12cs, ntak12 ) ) / sqrt ( x );
  }
  return value;
}
//****************************************************************************80

void r8_beskes ( double xnu, double x, int nin, double bke[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESKES: a sequence of exponentially scaled K Bessel functions at X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double XNU, ?
//    |XNU| < 1.
//
//    Input, double X, the argument.
//
//    Input, int NIN, indicates the number of terms to compute.
//
//    Output, double BKE(abs(NIN)), the exponentially scaled 
//    K Bessel functions.
//
{
  static double alnbig = 0.0;
  double bknu1;
  double direct;
  int i;
  int iswtch;
  int n;
  double v;
  double vend;
  double vincr;

  if ( alnbig == 0.0 )
  {
    alnbig = log ( r8_mach ( 2 ) );
  }

  v = r8_abs ( xnu );
  n = i4_abs ( nin );

  if ( 1.0 <= v )
  {
    cerr << "\n";
    cerr << "R8_BESKES - Fatal error!\n";
    cerr << "  |XNU| must be less than 1.\n";
    exit ( 1 );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESKES - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    cerr << "\n";
    cerr << "R8_BESKES - Fatal error!\n";
    cerr << "  N = 0.\n";
    exit ( 1 );
  }

  r8_knus ( v, x, bke[0], bknu1, iswtch );

  if ( n == 1 )
  {
    return;
  }

  if ( nin < 0 )
  {
    vincr = - 1.0;
  }
  else
  {
    vincr = + 1.0;
  }

  if ( xnu < 0.0 )
  {
    direct = - vincr;
  }
  else
  {
    direct = vincr;
  }

  bke[1] = bknu1;

  if ( direct < 0.0 )
  {
    r8_knus ( r8_abs ( xnu + vincr ), x, bke[1], bknu1, iswtch );
  }

  if ( n == 2 )
  {
    return;
  }

  vend = r8_abs ( xnu + ( double ) ( nin ) ) - 1.0;

  v = xnu;
  for ( i = 3; i <= n; i++ )
  {
    v = v + vincr;
    bke[i-1] = 2.0 * v * bke[i-2] / x + bke[i-3];
  }
  return;
}
//****************************************************************************80

void r8_besks ( double xnu, double x, int nin, double bk[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESKS evaluates a sequence of K Bessel functions at X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double XNU, ?
//    |XNU| < 1.
//
//    Input, double X, the argument.
//
//    Input, int NIN, indicates the number of terms to compute.
//
//    Output, double BK(abs(NIN)), the K Bessel functions.
//
{
  double expxi;
  int i;
  int n;
  static double xmax = 0.0;

  if ( xmax == 0.0 )
  {
    xmax = - log ( r8_mach ( 1 ) );
    xmax = xmax + 0.5 * log ( 3.14 * 0.5 / xmax );
  }

  r8_beskes ( xnu, x, nin, bk );

  expxi = exp ( - x );
  n = i4_abs ( nin );

  for ( i = 0; i < n; i++ )
  {
    bk[i] = expxi * bk[i];
  }
  return;
}
//****************************************************************************80

double r8_besy0 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESY0 evaluates the Bessel function Y of order 0 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESY0, the Bessel function Y of order 0 of X.
//
{
  static double alnhaf = -0.69314718055994530941723212145818;
  double ampl;
  static double by0cs[19] = {
    -0.1127783939286557321793980546028E-01,
    -0.1283452375604203460480884531838,
    -0.1043788479979424936581762276618,
    +0.2366274918396969540924159264613E-01,
    -0.2090391647700486239196223950342E-02,
    +0.1039754539390572520999246576381E-03,
    -0.3369747162423972096718775345037E-05,
    +0.7729384267670667158521367216371E-07,
    -0.1324976772664259591443476068964E-08,
    +0.1764823261540452792100389363158E-10,
    -0.1881055071580196200602823012069E-12,
    +0.1641865485366149502792237185749E-14,
    -0.1195659438604606085745991006720E-16,
    +0.7377296297440185842494112426666E-19,
    -0.3906843476710437330740906666666E-21,
    +0.1795503664436157949829120000000E-23,
    -0.7229627125448010478933333333333E-26,
    +0.2571727931635168597333333333333E-28,
    -0.8141268814163694933333333333333E-31 };
  static int nty0 = 0;
  double theta;
  static double twodpi = 0.636619772367581343075535053490057;
  double value;
  static double xsml = 0.0;
  double y;

  if ( nty0 == 0 )
  {
    nty0 = r8_inits ( by0cs, 19, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESY0 - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0;
    value = twodpi * ( alnhaf + log ( x ) ) * r8_besj0 ( x ) 
      + 0.375 + r8_csevl ( 0.125 * y - 1.0, by0cs, nty0 );
  }
  else if ( x <= 4.0 )
  {
    y = x * x;
    value = twodpi * ( alnhaf + log ( x ) ) * r8_besj0 ( x ) + 0.375 
      + r8_csevl ( 0.125 * y - 1.0, by0cs, nty0 );
  }
  else
  {
    r8_b0mp ( x, ampl, theta );
    value = ampl * sin ( theta );
  }
  return value;
}
//****************************************************************************80

double r8_besy1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESY1 evaluates the Bessel function Y of order 1 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESY1, the Bessel function Y of order 1 of X.
//
{
  double ampl;
  double static by1cs[20] = {
    +0.320804710061190862932352018628015E-01,
    +0.126270789743350044953431725999727E+01,
    +0.649996189992317500097490637314144E-02,
    -0.893616452886050411653144160009712E-01,
    +0.132508812217570954512375510370043E-01,
    -0.897905911964835237753039508298105E-03,
    +0.364736148795830678242287368165349E-04,
    -0.100137438166600055549075523845295E-05,
    +0.199453965739017397031159372421243E-07,
    -0.302306560180338167284799332520743E-09,
    +0.360987815694781196116252914242474E-11,
    -0.348748829728758242414552947409066E-13,
    +0.278387897155917665813507698517333E-15,
    -0.186787096861948768766825352533333E-17,
    +0.106853153391168259757070336000000E-19,
    -0.527472195668448228943872000000000E-22,
    +0.227019940315566414370133333333333E-24,
    -0.859539035394523108693333333333333E-27,
    +0.288540437983379456000000000000000E-29,
    -0.864754113893717333333333333333333E-32 };
  static int nty1 = 0;
  double theta;
  static double twodpi = 0.636619772367581343075535053490057;
  double value;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( nty1 == 0 )
  {
    nty1 = r8_inits ( by1cs, 20, 0.1 * r8_mach ( 3 ) );
    xmin = 1.571 * exp ( r8_max ( log ( r8_mach ( 1 ) ), 
      - log ( r8_mach ( 2 ) ) ) + 0.01 );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESY1 - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xmin )
  {
    y = 0.0;
    value = twodpi * log ( 0.5 * x ) * r8_besj1 ( x ) 
      + ( 0.5 + r8_csevl ( 0.125 * y - 1.0, by1cs, nty1 ) ) / x;
  }
  else if ( x <= 4.0 )
  {
    y = x * x;
    value = twodpi * log ( 0.5 * x ) * r8_besj1 ( x ) 
      + ( 0.5 + r8_csevl ( 0.125 * y - 1.0, by1cs, nty1 ) ) / x;
  }
  else
  {
    r8_b1mp ( x, ampl, theta );
    value = ampl * sin ( theta );
  }
  return value;
}
//****************************************************************************80

double r8_beta ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BETA evaluates the beta function of R8 arguments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, B, the arguments.
//
//    Output, double R8_BETA, the beta function of A and B.
//
{
  static double alnsml = 0.0;
  double value;
  static double xmax = 0.0;
  double xmin;

  if ( xmax == 0.0 )
  {
    r8_gaml ( xmin, xmax );
    alnsml = log ( r8_mach ( 1 ) );
  }

  if ( a <= 0.0 || b <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BETA - Fatal error!\n";
    cerr << "  A and B must be greater than 0.\n";
    exit ( 1 );
  }

  if ( a + b < xmax )
  {
    value = r8_gamma ( a ) * r8_gamma ( b ) / r8_gamma ( a + b );
    return value;
  }

  value = r8_lbeta ( a, b );

  value = exp ( value );

  return value;
}
//****************************************************************************80

double r8_betai ( double x, double pin, double qin )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BETAI evaluates the incomplete beta ratio of R8 arguments.
//
//  Discussion:
//
//    The incomplete Beta function ratio is the probability that a
//    random variable from a beta distribution having parameters
//    P and Q will be less than or equal to X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Nancy Bosten, EL Battiste,
//    Remark on Algorithm 179: 
//    Incomplete Beta Ratio,
//    Communications of the ACM,
//    Volume 17, Number 3, March 1974, pages 156-157.
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the upper limit of integration.
//    0.0 <= X <= 1.0.
//
//    Input, double PIN, the first distribution parameter.
//    0.0 < PIN.
//
//    Input, double QIN, the second distribution parameter.
//    0.0 < QIN.
//
//    Output, double R8_BETAI, the incomplete beta function ratio.
//
{
  static double alneps = 0.0;
  static double alnsml = 0.0;
  double c;
  static double eps = 0.0;
  double finsum;
  int i;
  int ib;
  int n;
  double p;
  double p1;
  double ps;
  double q;
  static double sml = 0.0;
  double term;
  double value;
  double xb;
  double xi;
  double y;

  if ( eps == 0.0 )
  {
    eps = r8_mach ( 3 );
    alneps = log ( eps );
    sml = r8_mach ( 1 );
    alnsml = log ( sml );
  }

  if ( x < 0.0 || 1.0 < x )
  {
    cerr << "\n";
    cerr << "R8_BETAI - Fatal error!\n";
    cerr << "  0 <= X <= 1 is required.\n";
    exit ( 1 );
  }

  if ( pin <= 0.0 || qin <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BETAI - Fatal error!\n";
    cerr << "  P or Q <= 0.0.\n";
    exit ( 1 );
  }

  y = x;
  p = pin;
  q = qin;

  if ( p < q || 0.8 <= x )
  {
    if ( 0.2 <= x )
    {
      y = 1.0 - y;
      p = qin;
      q = pin;
    }
  }

  if ( ( p + q ) * y / ( p + 1.0 ) < eps )
  {
    value = 0.0;
    xb = p * log ( r8_max ( y, sml ) ) - log ( p ) - r8_lbeta ( p, q );
    if ( alnsml < xb && y != 0.0 )
    {
      value = exp ( xb );
    }
    if ( y != x || p != pin )
    {
      value = 1.0 - value;
    }
    return value;
  }

  ps = q - r8_aint ( q );
  if ( ps == 0.0 )
  {
    ps = 1.0;
  }

  xb = p * log ( y ) - r8_lbeta ( ps, p ) - log ( p );

  if ( xb < alnsml )
  {
    value = 0.0;
  }
  else
  {
    value = exp ( xb );
    term = value * p;
    if ( ps != 1.0 )
    {
      n = ( int ) ( r8_max ( alneps / log ( y ), 4.0 ) );
      for ( i = 1; i <= n; i++ )
      {
        xi = ( double ) ( i );
        term = term * ( xi - ps ) * y / xi;
        value = value + term / ( p + xi );
      }
    }
  }

  if ( 1.0 < q )
  {
    xb = p * log ( y ) + q * log ( 1.0 - y ) 
      - r8_lbeta ( p, q ) - log ( q );
    ib = ( int ) ( r8_max ( xb / alnsml, 0.0 ) );
    term = exp ( xb - ( double ) ( ib ) * alnsml );
    c = 1.0 / ( 1.0 - y );
    p1 = q * c / ( p + q - 1.0 );

    finsum = 0.0;
    n = ( int ) ( q );
    if ( q == ( double ) ( n ) )
    {
      n = n - 1;
    }

    for ( i = 1; i <= n; i++ )
    {
      if ( p1 <= 1.0 && term / eps <= finsum )
      {
        break;
      }

      xi = ( double ) ( i );
      term = ( q - xi + 1.0 ) * c * term / ( p + q - xi );

      if ( 1.0 < term )
      {
        ib = ib - 1;
        term = term * sml;
      }

      if ( ib == 0 )
      {
        finsum = finsum + term;
      }
    }
    value = value + finsum;
  }

  if ( y != x || p != pin )
  {
    value = 1.0 - value;
  }

  if ( value < 0.0 )
  {
    value =  0.0;
  }

  if ( 1.0 < value )
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

double r8_bi ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BI evaluates the Airy function Bi of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BI, the Airy function Bi of X.
//
{
  static double bifcs[13] = {
    -0.16730216471986649483537423928176E-01,
    +0.10252335834249445611426362777757,
    +0.17083092507381516539429650242013E-02,
    +0.11862545467744681179216459210040E-04,
    +0.44932907017792133694531887927242E-07,
    +0.10698207143387889067567767663628E-09,
    +0.17480643399771824706010517628573E-12,
    +0.20810231071761711025881891834399E-15,
    +0.18849814695665416509927971733333E-18,
    +0.13425779173097804625882666666666E-21,
    +0.77159593429658887893333333333333E-25,
    +0.36533879617478566399999999999999E-28,
    +0.14497565927953066666666666666666E-31 };
  static double bif2cs[15] = {
    +0.0998457269381604104468284257993,
    +0.47862497786300553772211467318231,
    +0.25155211960433011771324415436675E-01,
    +0.58206938852326456396515697872216E-03,
    +0.74997659644377865943861457378217E-05,
    +0.61346028703493836681403010356474E-07,
    +0.34627538851480632900434268733359E-09,
    +0.14288910080270254287770846748931E-11,
    +0.44962704298334641895056472179200E-14,
    +0.11142323065833011708428300106666E-16,
    +0.22304791066175002081517866666666E-19,
    +0.36815778736393142842922666666666E-22,
    +0.50960868449338261333333333333333E-25,
    +0.60003386926288554666666666666666E-28,
    +0.60827497446570666666666666666666E-31 };
  static double bigcs[13] = {
    +0.22466223248574522283468220139024E-01,
    +0.37364775453019545441727561666752E-01,
    +0.44476218957212285696215294326639E-03,
    +0.24708075636329384245494591948882E-05,
    +0.79191353395149635134862426285596E-08,
    +0.16498079851827779880887872402706E-10,
    +0.24119906664835455909247501122841E-13,
    +0.26103736236091436985184781269333E-16,
    +0.21753082977160323853123792000000E-19,
    +0.14386946400390433219483733333333E-22,
    +0.77349125612083468629333333333333E-26,
    +0.34469292033849002666666666666666E-29,
    +0.12938919273216000000000000000000E-32 };
  static double big2cs[15] = {
    +0.033305662145514340465176188111647,
    +0.161309215123197067613287532084943,
    +0.631900730961342869121615634921173E-02,
    +0.118790456816251736389780192304567E-03,
    +0.130453458862002656147116485012843E-05,
    +0.937412599553521729546809615508936E-08,
    +0.474580188674725153788510169834595E-10,
    +0.178310726509481399800065667560946E-12,
    +0.516759192784958180374276356640000E-15,
    +0.119004508386827125129496251733333E-17,
    +0.222982880666403517277063466666666E-20,
    +0.346551923027689419722666666666666E-23,
    +0.453926336320504514133333333333333E-26,
    +0.507884996513522346666666666666666E-29,
    +0.491020674696533333333333333333333E-32 };
  double eta;
  static int nbif = 0;
  static int nbif2 = 0;
  static int nbig = 0;
  static int nbig2 = 0;
  double theta;
  double value;
  static double x3sml = 0.0;
  double xm;
  static double xmax = 0.0;
  double z;

  if ( nbif == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nbif = r8_inits ( bifcs, 13, eta );
    nbig = r8_inits ( bigcs, 13, eta );
    nbif2 = r8_inits ( bif2cs, 15, eta );
    nbig2 = r8_inits ( big2cs, 15, eta );
    x3sml = r8_power ( eta, 0.3333 );
    xmax = r8_power ( 1.5 * log ( r8_mach ( 2 ) ), 0.6666 );
  }

  if ( x < - 1.0 )
  {
    r8_aimp ( x, xm, theta );
    value = xm * sin ( theta );
  }
  else if ( r8_abs ( x ) <= x3sml )
  {
    z = 0.0;
    value = 0.625 + r8_csevl ( z, bifcs, nbif ) 
      + x * ( 0.4375 + r8_csevl ( z, bigcs, nbig ) );
  }
  else if ( x <= 1.0 )
  {
    z = x * x * x;
    value = 0.625 + r8_csevl ( z, bifcs, nbif ) 
      + x * ( 0.4375 + r8_csevl ( z, bigcs, nbig ) );
  }
  else if ( x <= 2.0 )
  {
    z = ( 2.0 * x * x * x - 9.0 ) / 7.0;
    value = 1.125 + r8_csevl ( z, bif2cs, nbif2 ) 
      + x * ( 0.625 + r8_csevl ( z, big2cs, nbig2 ) );
  }
  else
  {
    value = r8_bie ( x ) * exp ( 2.0 * x * sqrt ( x ) / 3.0 );
  }
  return value;
}
//****************************************************************************80

double r8_bid ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BID evaluates the derivative of the Airy function Bi of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BID, the derivative of the Airy function Bi of X.
//
{
  static double bif2cs[15] = {
     0.32349398760352203352119193596266015,
     0.08629787153556355913888835323811100,
     0.00299402555265539742613821050727155,
     0.00005143052836466163720464316950821,
     0.00000052584025003681146026033098613,
     0.00000000356175137395770028102730600,
     0.00000000001714686400714584830518308,
     0.00000000000006166351969232555406693,
     0.00000000000000017191082154315985806,
     0.00000000000000000038236889518803943,
     0.00000000000000000000069424173624884,
     0.00000000000000000000000104833932510,
     0.00000000000000000000000000133721972,
     0.00000000000000000000000000000145986,
      0.00000000000000000000000000000000138 };
  static double bifcs[13] = {
     0.115353679082857024267474446284908879,
     0.020500789404919287530357789445940252,
     0.000213529027890287581892679619451158,
     0.000001078396061467683042209155523569,
     0.000000003209470883320666783353670420,
     0.000000000006293040671833540390213316,
     0.000000000000008740304300063083340121,
     0.000000000000000009047915683496049529,
     0.000000000000000000007249923164709251,
     0.000000000000000000000004629576649604,
     0.000000000000000000000000002411236436,
     0.000000000000000000000000000001043825,
     0.000000000000000000000000000000000382 };
  static double big2cs[16] = {
     1.606299946362129457759284537862622883,
     0.744908881987608865201476685194753972,
     0.047013873861027737964095177635353019,
     0.001228442206254823907016188785848091,
     0.000017322241225662362670987355613727,
     0.000000152190165236801893711508366563,
     0.000000000911356024911957704145528786,
     0.000000000003954791842356644201722554,
     0.000000000000013001737033862320007309,
     0.000000000000000033493506858269079763,
     0.000000000000000000069419094403694057,
     0.000000000000000000000118248256604581,
     0.000000000000000000000000168462493472,
     0.000000000000000000000000000203684674,
     0.000000000000000000000000000000211619,
     0.000000000000000000000000000000000191 };
  static double bigcs[13] = {
    -0.0971964404164435373897790974606802,
     0.1495035768431670665710843445326264,
     0.0031135253871213260419419176839631,
     0.0000247085705798212967777021920569,
     0.0000001029496277313786081987324295,
     0.0000000002639703739869432892676778,
     0.0000000000004582792707803206608181,
     0.0000000000000005742829740893447321,
     0.0000000000000000005438275385238549,
     0.0000000000000000000004028347267083,
     0.0000000000000000000000002397823826,
     0.0000000000000000000000000001171956,
     0.0000000000000000000000000000000479 };
  double eta;
  static int nbif = 0;
  static int nbif2 = 0;
  static int nbig = 0;
  static int nbig2 = 0;
  double phi;
  double value;
  double x2;
  static double x2sml = 0.0;
  double x3;
  static double x3sml = 0.0;
  static double xmax = 0.0;
  double xn;
  double z;

  if ( nbif == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nbif = r8_inits ( bifcs, 13, eta );
    nbig = r8_inits ( bigcs, 13, eta );
    nbif2 = r8_inits ( bif2cs, 15, eta );
    nbig2 = r8_inits ( big2cs, 16, eta );
    x2sml = sqrt ( eta );
    x3sml = r8_power ( eta, 0.3333 );
    xmax = r8_power ( 1.5 * log ( r8_mach ( 2 ) ), 0.6666 );
  }

  if ( x < - 1.0 )
  {
    r8_admp ( x, xn, phi );
    value = xn * sin ( phi );
  }
  else if ( r8_abs ( x ) <= x2sml )
  {
    x2 = 0.0;
    x3 = 0.0;
    value = x2 * ( r8_csevl ( x3, bifcs, nbif ) + 0.25 ) 
      + r8_csevl ( x3, bigcs, nbig ) + 0.5;
  }
  else if ( r8_abs ( x ) <= x3sml )
  {
    x2 = x * x;
    x3 = 0.0;
    value = x2 * ( r8_csevl ( x3, bifcs, nbif ) + 0.25 ) 
      + r8_csevl ( x3, bigcs, nbig ) + 0.5;
  }
  else if ( x <= 1.0 )
  {
    x2 = x * x;
    x3 = x * x * x;
    value = x2 * ( r8_csevl ( x3, bifcs, nbif ) + 0.25 ) 
      + r8_csevl ( x3, bigcs, nbig ) + 0.5;
  }
  else if ( x <= 2.0 )
  {
    z = ( 2.0 * x * x * x - 9.0 ) / 7.0;
    value = x * x * ( r8_csevl ( z, bif2cs, nbif2 ) + 0.25 ) 
      + r8_csevl ( z, big2cs, nbig2 ) + 0.5;
  }
  else
  {
    value = r8_bide ( x ) * exp ( 2.0 * x * sqrt ( x ) / 3.0 );
  }
  return value;
}
//****************************************************************************80

double r8_bide ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BIDE: exponentially scaled derivative, Airy function Bi of an R8 argument.
//
//  Discussion:
//
//    if X < 0,
//      R8_BIDE ( X ) = R8_BID ( X )
//    else
//      R8_BIDE ( X ) = R8_BID ( X ) * exp ( - 2/3 * X^(3/2) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BIDE, the exponentially scaled derivative of 
//    the Airy function Bi of X.
//
{
  static double atr = 8.75069057084843450880771988210148;
  static double bif2cs[15] = {
     0.32349398760352203352119193596266015,
     0.08629787153556355913888835323811100,
     0.00299402555265539742613821050727155,
     0.00005143052836466163720464316950821,
     0.00000052584025003681146026033098613,
     0.00000000356175137395770028102730600,
     0.00000000001714686400714584830518308,
     0.00000000000006166351969232555406693,
     0.00000000000000017191082154315985806,
     0.00000000000000000038236889518803943,
     0.00000000000000000000069424173624884,
     0.00000000000000000000000104833932510,
     0.00000000000000000000000000133721972,
     0.00000000000000000000000000000145986,
     0.00000000000000000000000000000000138 };
  static double bifcs[13] = {
     0.115353679082857024267474446284908879,
     0.020500789404919287530357789445940252,
     0.000213529027890287581892679619451158,
     0.000001078396061467683042209155523569,
     0.000000003209470883320666783353670420,
     0.000000000006293040671833540390213316,
     0.000000000000008740304300063083340121,
     0.000000000000000009047915683496049529,
     0.000000000000000000007249923164709251,
     0.000000000000000000000004629576649604,
     0.000000000000000000000000002411236436,
     0.000000000000000000000000000001043825,
     0.000000000000000000000000000000000382 };
  static double big2cs[16] = {
     1.606299946362129457759284537862622883,
     0.744908881987608865201476685194753972,
     0.047013873861027737964095177635353019,
     0.001228442206254823907016188785848091,
     0.000017322241225662362670987355613727,
     0.000000152190165236801893711508366563,
     0.000000000911356024911957704145528786,
     0.000000000003954791842356644201722554,
     0.000000000000013001737033862320007309,
     0.000000000000000033493506858269079763,
     0.000000000000000000069419094403694057,
     0.000000000000000000000118248256604581,
     0.000000000000000000000000168462493472,
     0.000000000000000000000000000203684674,
     0.000000000000000000000000000000211619,
     0.000000000000000000000000000000000191 };
  static double bigcs[13] = {
    -0.0971964404164435373897790974606802,
     0.1495035768431670665710843445326264,
     0.0031135253871213260419419176839631,
     0.0000247085705798212967777021920569,
     0.0000001029496277313786081987324295,
     0.0000000002639703739869432892676778,
     0.0000000000004582792707803206608181,
     0.0000000000000005742829740893447321,
     0.0000000000000000005438275385238549,
     0.0000000000000000000004028347267083,
     0.0000000000000000000000002397823826,
     0.0000000000000000000000000001171956,
     0.0000000000000000000000000000000479 };
  static double bip1cs[47] = {
    -0.17291873510795537186124679823741003,
    -0.01493584929846943639486231021818675,
    -0.00054711049516785663990658697874460,
     0.00015379662929584083449573727856666,
     0.00001543534761921794131028948022869,
    -0.00000654341138519060129226087106765,
     0.00000037280824078787032232152275240,
     0.00000020720783881887480080810710514,
    -0.00000006581733364696191689495883922,
     0.00000000749267463539288212986048985,
     0.00000000111013368840707147698890101,
    -0.00000000072651405529159512323880794,
     0.00000000017827235598470153962165668,
    -0.00000000002173463524809506269656807,
    -0.00000000000203020349653882594017049,
     0.00000000000193118272294077519319859,
    -0.00000000000060449525048290296023117,
     0.00000000000012094496248933664277802,
    -0.00000000000001251088360074479784619,
    -0.00000000000000199173832424881344036,
     0.00000000000000151540816342864303038,
    -0.00000000000000049768927059816240250,
     0.00000000000000011545959731810501403,
    -0.00000000000000001863286862907983871,
     0.00000000000000000099330392344759104,
     0.00000000000000000068182083667412417,
    -0.00000000000000000034854456479650551,
     0.00000000000000000010860382134235961,
    -0.00000000000000000002599290185240166,
     0.00000000000000000000476895370459000,
    -0.00000000000000000000051946940777177,
    -0.00000000000000000000005925575044912,
     0.00000000000000000000005746008970972,
    -0.00000000000000000000002186119806494,
     0.00000000000000000000000624124294738,
    -0.00000000000000000000000146003421785,
     0.00000000000000000000000027493893904,
    -0.00000000000000000000000003474678018,
    -0.00000000000000000000000000109303694,
     0.00000000000000000000000000261972744,
    -0.00000000000000000000000000112365018,
     0.00000000000000000000000000035152059,
    -0.00000000000000000000000000009167601,
     0.00000000000000000000000000002040203,
    -0.00000000000000000000000000000373038,
     0.00000000000000000000000000000046070,
     0.00000000000000000000000000000001748 };
  static double bip2cs[88] = {
    -0.13269705443526630494937031210217135,
    -0.00568443626045977481306046339037428,
    -0.00015643601119611609623698471216660,
    -0.00001136737203679562267336053207940,
    -0.00000143464350991283669643136951338,
    -0.00000018098531185164131868746481700,
     0.00000000926177343610865546229511422,
     0.00000001710005490720592181887296162,
     0.00000000476698163503781708252686849,
    -0.00000000035195022023163141945397159,
    -0.00000000058890614315886871574147635,
    -0.00000000006678499607795537597612089,
     0.00000000006395565101720391190697713,
     0.00000000001554529427064394106403245,
    -0.00000000000792396999744612971684001,
    -0.00000000000258326242689717798947525,
     0.00000000000121655047787849117978773,
     0.00000000000038707207172899985942258,
    -0.00000000000022487045479618229130656,
    -0.00000000000004953476515684046293493,
     0.00000000000004563781601526912756017,
     0.00000000000000332998314345014118494,
    -0.00000000000000921750185832874202719,
     0.00000000000000094156670658958205765,
     0.00000000000000167153952640716157721,
    -0.00000000000000055134268782182410852,
    -0.00000000000000022368651572006617795,
     0.00000000000000017486948976520089209,
     0.00000000000000000206518666352329750,
    -0.00000000000000003973060018130712479,
     0.00000000000000001154836935724892335,
     0.00000000000000000553906053678276421,
    -0.00000000000000000457174427396478267,
     0.00000000000000000026567111858284432,
     0.00000000000000000101599148154167823,
    -0.00000000000000000044821231272196246,
    -0.00000000000000000007959149661617295,
     0.00000000000000000014583615616165794,
    -0.00000000000000000004015127893061405,
    -0.00000000000000000002079152963743616,
     0.00000000000000000001972630449634388,
    -0.00000000000000000000336033404001683,
    -0.00000000000000000000376504832685507,
     0.00000000000000000000269935508825595,
    -0.00000000000000000000026985946069808,
    -0.00000000000000000000061794011788222,
     0.00000000000000000000038782693311711,
    -0.00000000000000000000002420094005071,
    -0.00000000000000000000009844051058925,
     0.00000000000000000000005954353358494,
    -0.00000000000000000000000361274446366,
    -0.00000000000000000000001552634578088,
     0.00000000000000000000000977819380304,
    -0.00000000000000000000000092239447509,
    -0.00000000000000000000000241545903934,
     0.00000000000000000000000169558652255,
    -0.00000000000000000000000026762408641,
    -0.00000000000000000000000036188116265,
     0.00000000000000000000000030372404951,
    -0.00000000000000000000000007422876903,
    -0.00000000000000000000000004930678544,
     0.00000000000000000000000005468790028,
    -0.00000000000000000000000001920315188,
    -0.00000000000000000000000000516335154,
     0.00000000000000000000000000957723167,
    -0.00000000000000000000000000463659079,
    -0.00000000000000000000000000004509226,
     0.00000000000000000000000000155617519,
    -0.00000000000000000000000000104156509,
     0.00000000000000000000000000019565323,
     0.00000000000000000000000000021335380,
    -0.00000000000000000000000000021461958,
     0.00000000000000000000000000007875791,
     0.00000000000000000000000000001713768,
    -0.00000000000000000000000000003917137,
     0.00000000000000000000000000002233559,
    -0.00000000000000000000000000000269383,
    -0.00000000000000000000000000000577764,
     0.00000000000000000000000000000519650,
    -0.00000000000000000000000000000183361,
    -0.00000000000000000000000000000045763,
     0.00000000000000000000000000000099235,
    -0.00000000000000000000000000000058938,
     0.00000000000000000000000000000009568,
     0.00000000000000000000000000000013758,
    -0.00000000000000000000000000000014066,
     0.00000000000000000000000000000005964,
     0.00000000000000000000000000000000437 };
  static double btr = -2.09383632135605431360096498526268;
  double eta;
  static int nbif = 0;
  static int nbif2 = 0;
  static int nbig = 0;
  static int nbig2 = 0;
  static int nbip1 = 0;
  static int nbip2 = 0;
  double phi;
  double sqrtx;
  double value;
  double x2;
  static double x2sml = 0.0;
  double x3;
  static double x3sml = 0.0;
  static double x32sml = 0.0;
  static double xbig = 0.0;
  double xn;
  double z;

  if ( nbif == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nbif = r8_inits ( bifcs, 13, eta );
    nbig = r8_inits ( bigcs, 13, eta );
    nbif2 = r8_inits ( bif2cs, 15, eta );
    nbig2 = r8_inits ( big2cs, 16, eta );
    nbip1 = r8_inits ( bip1cs, 47, eta );
    nbip2 = r8_inits ( bip2cs, 88, eta );
    x2sml = sqrt ( eta );
    x3sml = r8_power ( eta, 0.3333 );
    x32sml = 1.3104 * x3sml * x3sml;
    xbig = r8_power ( r8_mach ( 2 ), 0.6666 );
  }

  if ( x < -1.0 )
  {
    r8_admp ( x, xn, phi );
    value = xn * sin ( phi );
  }
  else if ( r8_abs ( x ) <= x2sml )
  {
    x2 = 0.0;
    x3 = 0.0;
    value = x2 * ( r8_csevl ( x3, bifcs, nbif ) 
      + 0.25 ) + r8_csevl ( x3, bigcs, nbig ) + 0.5;
    if ( x32sml < x )
    {
      value = value * exp ( - 2.0 * x * sqrt ( x ) / 3.0 );
    }
  }
  else if ( r8_abs ( x ) <= x3sml )
  {
    x2 = x * x;
    x3 = 0.0;
    value = x2 * ( r8_csevl ( x3, bifcs, nbif ) 
      + 0.25 ) + r8_csevl ( x3, bigcs, nbig ) + 0.5;
    if ( x32sml < x )
    {
      value = value * exp ( - 2.0 * x * sqrt ( x ) / 3.0 );
    }
  }
  else if ( x <= 1.0 )
  {
    x2 = x * x;
    x3 = x * x * x;
    value = x2 * ( r8_csevl ( x3, bifcs, nbif ) 
      + 0.25 ) + r8_csevl ( x3, bigcs, nbig ) + 0.5;
    if ( x32sml < x )
    {
      value = value * exp ( - 2.0 * x * sqrt ( x ) / 3.0 );
    }
  }
  else if ( x <= 2.0 )
  {
    z = ( 2.0 * x * x * x - 9.0 ) / 7.0;
    value = exp ( - 2.0 * x * sqrt ( x ) / 3.0 ) 
      * ( x * x * ( 0.25 + r8_csevl ( z, bif2cs, nbif2 ) ) 
      + 0.5 + r8_csevl ( z, big2cs, nbig2 ) );
  }
  else if ( x <= 4.0 )
  {
    sqrtx = sqrt ( x );
    z = atr / x / sqrtx + btr;
    value = ( 0.625 + r8_csevl ( z, bip1cs, nbip1 ) ) * sqrt ( sqrtx );
  }
  else if ( x <= xbig )
  {
    sqrtx = sqrt ( x );
    z = 16.0 / x / sqrtx - 1.0;
    value = ( 0.625 + r8_csevl ( z, bip2cs, nbip2 ) ) * sqrt ( sqrtx );
  }
  else
  {
    sqrtx = sqrt ( x );
    z = -1.0;
    value = ( 0.625 + r8_csevl ( z, bip2cs, nbip2 ) ) * sqrt ( sqrtx );
  }
  return value;
}
//****************************************************************************80

double r8_bie ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BIE evaluates the exponentially scaled Airy function Bi of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BIE, the exponentially scaled Airy function Bi of X.
//
{
  static double atr = 8.75069057084843450880771988210148;
  static double bif2cs[15] = {
    +0.0998457269381604104468284257993,
    +0.47862497786300553772211467318231,
    +0.25155211960433011771324415436675E-01,
    +0.58206938852326456396515697872216E-03,
    +0.74997659644377865943861457378217E-05,
    +0.61346028703493836681403010356474E-07,
    +0.34627538851480632900434268733359E-09,
    +0.14288910080270254287770846748931E-11,
    +0.44962704298334641895056472179200E-14,
    +0.11142323065833011708428300106666E-16,
    +0.22304791066175002081517866666666E-19,
    +0.36815778736393142842922666666666E-22,
    +0.50960868449338261333333333333333E-25,
    +0.60003386926288554666666666666666E-28,
    +0.60827497446570666666666666666666E-31 };
  static double bifcs[13] = {
    -0.16730216471986649483537423928176E-01,
    +0.10252335834249445611426362777757,
    +0.17083092507381516539429650242013E-02,
    +0.11862545467744681179216459210040E-04,
    +0.44932907017792133694531887927242E-07,
    +0.10698207143387889067567767663628E-09,
    +0.17480643399771824706010517628573E-12,
    +0.20810231071761711025881891834399E-15,
    +0.18849814695665416509927971733333E-18,
    +0.13425779173097804625882666666666E-21,
    +0.77159593429658887893333333333333E-25,
    +0.36533879617478566399999999999999E-28,
    +0.14497565927953066666666666666666E-31 };
  static double big2cs[15] = {
    +0.033305662145514340465176188111647,
    +0.161309215123197067613287532084943,
    +0.631900730961342869121615634921173E-02,
    +0.118790456816251736389780192304567E-03,
    +0.130453458862002656147116485012843E-05,
    +0.937412599553521729546809615508936E-08,
    +0.474580188674725153788510169834595E-10,
    +0.178310726509481399800065667560946E-12,
    +0.516759192784958180374276356640000E-15,
    +0.119004508386827125129496251733333E-17,
    +0.222982880666403517277063466666666E-20,
    +0.346551923027689419722666666666666E-23,
    +0.453926336320504514133333333333333E-26,
    +0.507884996513522346666666666666666E-29,
    +0.491020674696533333333333333333333E-32 };
  static double bigcs[13] = {
    +0.22466223248574522283468220139024E-01,
    +0.37364775453019545441727561666752E-01,
    +0.44476218957212285696215294326639E-03,
    +0.24708075636329384245494591948882E-05,
    +0.79191353395149635134862426285596E-08,
    +0.16498079851827779880887872402706E-10,
    +0.24119906664835455909247501122841E-13,
    +0.26103736236091436985184781269333E-16,
    +0.21753082977160323853123792000000E-19,
    +0.14386946400390433219483733333333E-22,
    +0.77349125612083468629333333333333E-26,
    +0.34469292033849002666666666666666E-29,
    +0.12938919273216000000000000000000E-32 };
  static double bip1cs[47] = {
    -0.83220474779434474687471864707973E-01,
    +0.11461189273711742889920226128031E-01,
    +0.42896440718911509494134472566635E-03,
    -0.14906639379950514017847677732954E-03,
    -0.13076597267876290663136340998881E-04,
    +0.63275983961030344754535716032494E-05,
    -0.42226696982681924884778515889433E-06,
    -0.19147186298654689632835494181277E-06,
    +0.64531062845583173611038157880934E-07,
    -0.78448546771397719289748310448628E-08,
    -0.96077216623785085879198533565432E-09,
    +0.70004713316443966339006074402068E-09,
    -0.17731789132814932022083128056698E-09,
    +0.22720894783465236347282126389311E-10,
    +0.16540456313972049847032860681891E-11,
    -0.18517125559292316390755369896693E-11,
    +0.59576312477117290165680715534277E-12,
    -0.12194348147346564781055769498986E-12,
    +0.13347869253513048815386347813597E-13,
    +0.17278311524339746664384792889731E-14,
    -0.14590732013016720735268871713166E-14,
    +0.49010319927115819978994989520104E-15,
    -0.11556545519261548129262972762521E-15,
    +0.19098807367072411430671732441524E-16,
    -0.11768966854492179886913995957862E-17,
    -0.63271925149530064474537459677047E-18,
    +0.33861838880715361614130191322316E-18,
    -0.10725825321758625254992162219622E-18,
    +0.25995709605617169284786933115562E-19,
    -0.48477583571081193660962309494101E-20,
    +0.55298913982121625361505513198933E-21,
    +0.49421660826069471371748197444266E-22,
    -0.55162121924145707458069720814933E-22,
    +0.21437560417632550086631884499626E-22,
    -0.61910313387655605798785061137066E-23,
    +0.14629362707391245659830967336959E-23,
    -0.27918484471059005576177866069333E-24,
    +0.36455703168570246150906795349333E-25,
    +0.58511821906188711839382459733333E-27,
    -0.24946950487566510969745047551999E-26,
    +0.10979323980338380977919579477333E-26,
    -0.34743388345961115015034088106666E-27,
    +0.91373402635349697363171082240000E-28,
    -0.20510352728210629186247720959999E-28,
    +0.37976985698546461748651622399999E-29,
    -0.48479458497755565887848448000000E-30,
    -0.10558306941230714314205866666666E-31 };
  static double bip2cs[88] = {
    -0.11359673758598867913797310895527,
    +0.41381473947881595760052081171444E-02,
    +0.13534706221193329857696921727508E-03,
    +0.10427316653015353405887183456780E-04,
    +0.13474954767849907889589911958925E-05,
    +0.16965374054383983356062511163756E-06,
    -0.10096500865641624301366228396373E-07,
    -0.16729119493778475127836973095943E-07,
    -0.45815364485068383217152795613391E-08,
    +0.37366813665655477274064749384284E-09,
    +0.57669303201452448119584643502111E-09,
    +0.62181265087850324095393408792371E-10,
    -0.63294120282743068241589177281354E-10,
    -0.14915047908598767633999091989487E-10,
    +0.78896213942486771938172394294891E-11,
    +0.24960513721857797984888064000127E-11,
    -0.12130075287291659477746664734814E-11,
    -0.37404939108727277887343460402716E-12,
    +0.22377278140321476798783446931091E-12,
    +0.47490296312192466341986077472514E-13,
    -0.45261607991821224810605655831294E-13,
    -0.30172271841986072645112245876020E-14,
    +0.91058603558754058327592683478908E-14,
    -0.98149238033807062926643864207709E-15,
    -0.16429400647889465253601245251589E-14,
    +0.55334834214274215451182114635164E-15,
    +0.21750479864482655984374381998156E-15,
    -0.17379236200220656971287029558087E-15,
    -0.10470023471443714959283909313604E-17,
    +0.39219145986056386925441403311462E-16,
    -0.11621293686345196925824005665910E-16,
    -0.54027474491754245533735411307773E-17,
    +0.45441582123884610882675428553304E-17,
    -0.28775599625221075729427585480086E-18,
    -0.10017340927225341243596162960440E-17,
    +0.44823931215068369856332561906313E-18,
    +0.76135968654908942328948982366775E-19,
    -0.14448324094881347238956060145422E-18,
    +0.40460859449205362251624847392112E-19,
    +0.20321085700338446891325190707277E-19,
    -0.19602795471446798718272758041962E-19,
    +0.34273038443944824263518958211738E-20,
    +0.37023705853905135480024651593154E-20,
    -0.26879595172041591131400332966712E-20,
    +0.28121678463531712209714454683364E-21,
    +0.60933963636177797173271119680329E-21,
    -0.38666621897150844994172977893413E-21,
    +0.25989331253566943450895651927228E-22,
    +0.97194393622938503767281175216084E-22,
    -0.59392817834375098415630478204591E-22,
    +0.38864949977113015409591960439444E-23,
    +0.15334307393617272869721512868769E-22,
    -0.97513555209762624036336521409724E-23,
    +0.96340644440489471424741339383726E-24,
    +0.23841999400208880109946748792454E-23,
    -0.16896986315019706184848044205207E-23,
    +0.27352715888928361222578444801478E-24,
    +0.35660016185409578960111685025730E-24,
    -0.30234026608258827249534280666954E-24,
    +0.75002041605973930653144204823232E-25,
    +0.48403287575851388827455319838748E-25,
    -0.54364137654447888432698010297766E-25,
    +0.19281214470820962653345978809756E-25,
    +0.50116355020532656659611814172172E-26,
    -0.95040744582693253786034620869972E-26,
    +0.46372646157101975948696332245611E-26,
    +0.21177170704466954163768170577046E-28,
    -0.15404850268168594303692204548726E-26,
    +0.10387944293201213662047889194441E-26,
    -0.19890078156915416751316728235153E-27,
    -0.21022173878658495471177044522532E-27,
    +0.21353099724525793150633356670491E-27,
    -0.79040810747961342319023537632627E-28,
    -0.16575359960435585049973741763592E-28,
    +0.38868342850124112587625586496537E-28,
    -0.22309237330896866182621562424717E-28,
    +0.27777244420176260265625977404382E-29,
    +0.57078543472657725368712433782772E-29,
    -0.51743084445303852800173371555280E-29,
    +0.18413280751095837198450927071569E-29,
    +0.44422562390957094598544071068647E-30,
    -0.98504142639629801547464958226943E-30,
    +0.58857201353585104884754198881995E-30,
    -0.97636075440429787961402312628595E-31,
    -0.13581011996074695047063597884122E-30,
    +0.13999743518492413270568048380345E-30,
    -0.59754904545248477620884562981118E-31,
    -0.40391653875428313641045327529856E-32 };
  static double btr = -2.09383632135605431360096498526268;
  double eta;
  static int nbif = 0;
  static int nbif2 = 0;
  static int nbig = 0;
  static int nbig2 = 0;
  static int nbip1 = 0;
  static int nbip2 = 0;
  double sqrtx;
  double theta;
  double value;
  static double x32sml;
  static double x3sml;
  static double xbig;
  double xm;
  double z;

  if ( nbif == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nbif = r8_inits ( bifcs, 13, eta );
    nbig = r8_inits ( bigcs, 13, eta );
    nbif2 = r8_inits ( bif2cs, 15, eta );
    nbig2 = r8_inits ( big2cs, 15, eta );
    nbip1 = r8_inits ( bip1cs, 47, eta );
    nbip2 = r8_inits ( bip2cs, 88, eta );
    x3sml = r8_power ( eta, 0.3333 );
    x32sml = 1.3104 * x3sml * x3sml;
    xbig = r8_power ( r8_mach ( 2 ), 0.6666 );
  }

  if ( x < - 1.0 )
  {
    r8_aimp ( x, xm, theta );
    value = xm * sin ( theta );
  }
  else if ( r8_abs ( x ) <= x3sml )
  {
    z = 0.0;
    value = 0.625 + r8_csevl ( z, bifcs, nbif ) 
      + x * ( 0.4375 + r8_csevl ( z, bigcs, nbig ) );
    if (  x32sml <= x )
    {
      value = value * exp ( - 2.0 * x * sqrt ( x ) / 3.0 );
    }
  }
  else if ( x <= 1.0 )
  {
    z = x * x * x;
    value = 0.625 + r8_csevl ( z, bifcs, nbif ) 
      + x * ( 0.4375 + r8_csevl ( z, bigcs, nbig ) );
    if (  x32sml <= x )
    {
      value = value * exp ( - 2.0 * x * sqrt ( x ) / 3.0 );
    }
  }
  else if ( x <= 2.0 )
  {
    z = ( 2.0 * x * x * x - 9.0 ) / 7.0;
    value = exp ( - 2.0 * x * sqrt ( x ) / 3.0 ) 
      * ( 1.125 + r8_csevl ( z, bif2cs, nbif2 ) 
      + x * ( 0.625 + r8_csevl ( z, big2cs, nbig2 ) ) );
  }
  else if ( x <= 4.0 )
  {
    sqrtx = sqrt ( x );
    z = atr / x / sqrtx + btr;
    value = ( 0.625 + r8_csevl ( z, bip1cs, nbip1 ) ) / sqrt ( sqrtx );
  }
  else if ( x < xbig )
  {
    sqrtx = sqrt ( x );
    z = 16.0 / ( x * sqrtx ) - 1.0;
    value = ( 0.625 + r8_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx );
  }
  else
  {
    sqrtx = sqrt ( x );
    z = - 1.0;
    value = ( 0.625 + r8_csevl ( z, bip2cs, nbip2 ) ) / sqrt ( sqrtx );
  }
  return value;
}
//****************************************************************************80

double r8_binom ( int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BINOM evaluates the binomial coefficient using R8 arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, int N, M, the arguments.
//
//    Output, double R8_BINOM, the binomial coefficient.
//
{
  static double bilnmx = 0.0;
  double corr;
  static double fintmx = 0.0;
  int i;
  int k;
  static double sq2pil = 0.91893853320467274178032973640562;
  double value;
  double xk;
  double xn;
  double xnk;

  if ( bilnmx == 0.0 )
  {
    bilnmx = log ( r8_mach ( 2 ) ) - 0.0001;
    fintmx = 0.9 / r8_mach ( 3 );
  }

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "R8_BINOM - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  if ( m < 0 )
  {
    cerr << "\n";
    cerr << "R8_BINOM - Fatal error!\n";
    cerr << "  M < 0.\n";
    exit ( 1 );
  }

  if ( n < m )
  {
    cerr << "\n";
    cerr << "R8_BINOM - Fatal error!\n";
    cerr << "  N < M.\n";
    exit ( 1 );
  }

  k = i4_min ( m, n - m );

  if ( k <= 20 && 
    ( double ) ( k ) * log ( ( double ) ( i4_max ( n, 1 ) ) ) <= bilnmx )
  {
    value = 1.0;
    for ( i = 1; i <= k; i++ )
    {
      value = value * ( double ) ( n - i + 1 ) / ( double ) ( i );
    }
  }
  else
  {
    if ( k < 9 )
    {
      cerr << "\n";
      cerr << "R8_BINOM - Fatal error!\n";
      cerr << "  Result overflows.\n";
      cerr << "  N or M is too big.\n";
      exit ( 1 );
    }

    xn = ( double ) ( n + 1 );
    xk = ( double ) ( k + 1 );
    xnk = ( double ) ( n - k + 1 );

    corr = r8_lgmc ( xn ) - r8_lgmc ( xk ) - r8_lgmc ( xnk );

    value = xk * log ( xnk / xk ) 
      - xn * r8_lnrel ( - ( xk - 1.0 ) / xn )
      - 0.5 * log ( xn * xnk / xk ) + 1.0 - sq2pil + corr;

    if ( bilnmx < value )
    {
      cerr << "\n";
      cerr << "R8_BINOM - Fatal error!\n";
      cerr << "  Result overflows.\n";
      cerr << "  N or M is too big.\n";
      exit ( 1 );
    }
    value = exp ( value );
  }

  if ( value < fintmx )
  {
    value = r8_aint ( value + 0.5 );
  }
  return value;
}
//****************************************************************************80

double r8_cbrt ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CBRT computes the cube root of an R8.
//
//  Discussion:
//
//    The approximation is a generalized Chebyshev series converted
//    to polynomial form.  The approximation is nearly best in the 
//    sense of relative error with 4.085 digits accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the number whose square root is desired.
//
//    Output, double R8_CBRT, the cube root of X.
//
{
  static double cbrt2[5] = {
    0.62996052494743658238360530363911,
    0.79370052598409973737585281963615,
    1.0,
    1.25992104989487316476721060727823,
    1.58740105196819947475170563927231 };
  int irem;
  int iter;
  int ixpnt;
  int n;
  static int niter = 0;
  double value;
  double vsq;
  double y;

  if ( niter == 0 )
  {
    niter = ( int ) ( 1.443 * log ( - 0.106 * log ( 0.1 * r8_mach ( 3 ) ) ) + 1.0 );
  }

  value = 0.0;

  if ( x != 0.0 )
  {
    r8_upak ( r8_abs ( x ), y, n );
    ixpnt = n / 3;
    irem = n - 3 * ixpnt + 3;

    value = 0.439581 + y * ( 
            0.928549 + y * (
          - 0.512653 + y *
            0.144586 ) );

    for ( iter = 1; iter <= niter; iter++ )
    {
      vsq = value * value;
      value = value + ( y - value * vsq ) / ( 3.0 * vsq );
    }

    if ( x < 0.0 )
    {
      value = - r8_abs ( value );
    }
    else
    {
      value = + r8_abs ( value );
    }
    value = r8_pak ( cbrt2[irem-1] * value, ixpnt );
  }
  return value;
}
//****************************************************************************80

double r8_chi ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHI evaluates the hyperbolic cosine integral of an R8 argument.
//
//  Discussion:
//
//    The hyperbolic cosine integral is defined by
//
//      CHI(X) = gamma + log ( x ) 
//        + integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
//
//    where gamma is Euler's constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_CHI, the hyperbolic cosine integral
//    evaluated at X.
//
{
  double value;

  value = 0.5 * ( r8_ei ( x ) - r8_e1 ( x ) );

  return value;
}
//****************************************************************************80

double r8_chu ( double a, double b, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHU evaluates the confluent hypergeometric function of R8 arguments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, B, the parameters.
//
//    Input, double X, the argument.
//
//    Output, double R8_CHU, the function value.
//
{
  double a0;
  double aintb;
  double alnx;
  double b0;
  double beps;
  double c0;
  static double eps = 0.0;
  double factor;
  double gamri1;
  double gamrni;
  int i;
  int istrt;
  int m;
  int n;
  double pch1ai;
  double pch1i;
  static double pi = 3.141592653589793238462643383279503;
  double pochai;
  double sum;
  double t;
  double value;
  double xeps1;
  double xi;
  double xi1;
  double xn;
  double xtoeps;

  if ( eps == 0.0 )
  {
    eps = r8_mach ( 3 );
  }

  if ( x < 0.0 )
  {
    cerr << "\n";
    cerr << "R8_CHU - Fatal error!\n";
    cerr << "  X < 0.\n";
    exit ( 1 );
  }

  if ( x == 0.0 )
  {
    if ( 1.0 <= b )
    {
      cerr << "\n";
      cerr << "R8_CHU - Fatal error!\n";
      cerr << "  X = 0 and 1 <= B.\n";
      exit ( 1 );
    }
    value = r8_gamma ( 1.0 - b ) / r8_gamma ( 1.0 + a - b );
    return value;
  }

  if ( r8_max ( r8_abs ( a ), 1.0 ) 
    * r8_max ( r8_abs ( 1.0 + a - b ), 1.0 ) < 0.99 * r8_abs ( x ) )
  {
    value = r8_power ( x, - a ) * r8_chu_scaled ( a, b, x );
    return value;
  }
//
//  The ascending series will be used, because the descending rational
//  approximation (which is based on the asymptotic series) is unstable.
//
  if ( 0.0 <= b )
  {
    aintb = r8_aint ( b + 0.5 );
  }
  else
  {
    aintb = r8_aint ( b - 0.5 );
  }
  beps = b - aintb;
  n = ( int ) aintb;
  alnx = log ( x );
  xtoeps = exp ( - beps * alnx );
//
//  Evaluate the finite sum.
//
//  Consider the case b < 1.0 first.
//
  if ( n < 1 )
  {
    sum = 1.0; 
    t = 1.0;
    m = - n;
    for ( i = 1; i <= m; i++ )
    {
      xi1 = ( double ) ( i - 1 );
      t = t * ( a + xi1 ) * x / ( ( b + xi1 ) * ( xi1 + 1.0 ) );
      sum = sum + t;
    }
    sum = r8_poch ( 1.0 + a - b, - a ) * sum;
  }
//
//  Now consider the case 1 <= b.
//
  else
  {
    sum = 0.0;
    m = n - 2;

    if ( 0 <= m )
    {
      t = 1.0;
      sum = 1.0;

      for ( i = 1; i <= m; i++ )
      {
        xi = ( double ) ( i );
        t = t * ( a - b + xi ) * x / ( ( 1.0 - b + xi ) * xi );
        sum = sum + t;
      }

      sum = r8_gamma ( b - 1.0 ) * r8_gamr ( a )
        * r8_power ( x, ( double ) ( 1 - n ) ) * xtoeps * sum;
    }
  }
//
//  Next evaluate the infinite sum.
//
  if ( n < 1 )
  {
    istrt = 1 - n;
  }
  else
  {
    istrt = 0;
  }

  xi = ( double ) ( istrt );

  factor = r8_mop ( n ) * r8_gamr ( 1.0 + a - b ) * r8_power ( x, xi );

  if ( beps != 0.0 )
  {
    factor = factor * beps * pi / sin ( beps * pi );
  }

  pochai = r8_poch ( a, xi );
  gamri1 = r8_gamr ( xi + 1.0 );
  gamrni = r8_gamr ( aintb + xi );
  b0 = factor * r8_poch ( a, xi - beps ) 
    * gamrni * r8_gamr ( xi + 1.0 - beps );
//
//  x^(-beps) is close to 1.0, so we must be careful in evaluating the
//  differences.
//
  if ( r8_abs ( xtoeps - 1.0 ) <= 0.5 )
  {
    pch1ai = r8_poch1 ( a + xi, -beps );
    pch1i = r8_poch1 ( xi + 1.0 - beps, beps );
    c0 = factor * pochai * gamrni * gamri1 * (
      - r8_poch1 ( b + xi,- beps ) + pch1ai 
      - pch1i + beps * pch1ai * pch1i );
//
//  xeps1 = (1.0 - x^(-beps))/beps = (x^(-beps) - 1.0)/(-beps)
//
    xeps1 = alnx* r8_exprel ( - beps * alnx );

    value = sum + c0 + xeps1 * b0;
    xn = ( double ) ( n );

    for ( i = 1; i <= 1000; i++ )
    {
      xi = ( double ) ( istrt + i );
      xi1 = ( double ) ( istrt + i - 1 );
      b0 = ( a + xi1 - beps ) * b0 * x 
        / ( ( xn + xi1 ) * ( xi - beps ) );
      c0 = ( a + xi1 ) * c0 * x / ( ( b + xi1) * xi )
        - ( ( a - 1.0 ) * ( xn + 2.0 * xi - 1.0 ) 
        + xi * ( xi - beps ) ) * b0 
        / ( xi * ( b + xi1 ) * ( a + xi1 - beps ) );
      t = c0 + xeps1 * b0;
      value = value + t;
      if ( r8_abs ( t ) < eps * r8_abs ( value ) )
      {
        return value;
      }
    }

    cerr << "\n";
    cerr << "R8_CHU - Fatal error!\n";
    cerr << "  No convergence in 1000 terms.\n";
    exit ( 1 );
  }
//
//  x^(-beps) is very different from 1.0, so the straightforward
//  formulation is stable.
//
  a0 = factor * pochai * r8_gamr ( b + xi ) * gamri1 / beps;
  b0 = xtoeps * b0 / beps;

  value = sum + a0 - b0;

  for ( i = 1; i <= 1000; i++ )
  {
    xi = ( double ) ( istrt + i );
    xi1 = ( double ) ( istrt + i - 1 );
    a0 = ( a + xi1 ) * a0 * x / ( ( b + xi1 ) * xi );
    b0 = ( a + xi1 - beps ) * b0 * x
      / ( ( aintb + xi1 ) * ( xi - beps ) );
    t = a0 - b0;
    value = value + t;
    if ( r8_abs ( t ) < eps * r8_abs ( value ) )
    {
      return value;
    }
  }

  cerr << "\n";
  cerr << "R8_CHU - Fatal error!\n";
  cerr << "  No convergence in 1000 terms.\n";
  exit ( 1 );
}
//****************************************************************************80

double r8_chu_scaled ( double a, double b, double z )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHU_SCALED: scaled confluent hypergeometric function of R8 arguments.
//
//  Discussion:
//
//    Evaluate, for large z, z^a * u(a,b,z)  where U is the logarithmic
//    confluent hypergeometric function.  A rational approximation due to
//    Y L Luke is used.  When U is not in the asymptotic region, that is, when A
//    or B is large compared with Z, considerable significance loss occurs.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, B, the parameters.
//
//    Input, double Z, the argument.
//
//    Output, double R8_CHU_SCALED, the function value.
//
{
  double aa[4];
  double ab;
  double anbn;
  double bb[4];
  double bp;
  double c2;
  double ct1;
  double ct2;
  double ct3;
  double d1z;
  static double eps = 0.0;
  double g1;
  double g2;
  double g3;
  int i;
  int j;
  double sab;
  static double sqeps = 0.0;
  double value;
  double x2i1;

  if ( eps == 0.0 )
  {
    eps = 4.0 * r8_mach ( 4 );
    sqeps = sqrt ( r8_mach ( 4 ) );
  }

  bp = 1.0 + a - b;
  ab = a * bp;
  ct2 = 2.0 * ( z - ab );
  sab = a + bp;

  bb[0] = 1.0;
  aa[0] = 1.0;

  ct3 = sab + 1.0 + ab;
  bb[1] = 1.0 + 2.0 * z / ct3;
  aa[1] = 1.0 + ct2 / ct3;

  anbn = ct3 + sab + 3.0;
  ct1 = 1.0 + 2.0 * z / anbn;
  bb[2] = 1.0 + 6.0 * ct1 * z / ct3;
  aa[2] = 1.0 + 6.0 * ab / anbn + 3.0 * ct1 * ct2 / ct3;

  for ( i = 4; i <= 300; i++ )
  {
    x2i1 = ( double ) ( 2 * i - 3 );
    ct1 = x2i1 / ( x2i1 - 2.0 );
    anbn = anbn + x2i1 + sab;
    ct2 = ( x2i1 - 1.0 ) /anbn;
    c2 = x2i1 * ct2 - 1.0;
    d1z = x2i1 * 2.0 * z / anbn;

    ct3 = sab * ct2;
    g1 = d1z + ct1 * ( c2 + ct3 );
    g2 = d1z - c2;
    g3 = ct1 * ( 1.0 - ct3 - 2.0 * ct2 );

    bb[3] = g1 * bb[2] + g2 * bb[1] + g3 * bb[0];
    aa[3] = g1 * aa[2] + g2 * aa[1] + g3 * aa[0];

    value = aa[3] / bb[3];

    if ( r8_abs ( value - aa[0] / bb[0] ) < eps * r8_abs ( value ) )
    {
      return value;
    }

    for ( j = 0; j < 3; j++ )
    {
      aa[j] = aa[j+1];
      bb[j] = bb[j+1];
    }
  }

  cerr << "\n";
  cerr << "R8_CHU_SCALED - Fatal error!\n";
  cerr << "  No convergence after 300 terms.\n";
  exit ( 1 );
}
//****************************************************************************80

double r8_ci ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CI evaluates the cosine integral Ci of an R8 argument.
//
//  Discussion:
//
//    The cosine integral is defined by
//
//      CI(X) = - integral ( X <= T < Infinity ) ( cos ( T ) ) / T  dT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_CI, the cosine integral Ci evaluated at X.
//
{
  static double cics[19] = {
    -0.34004281856055363156281076633129873,
    -1.03302166401177456807159271040163751,
     0.19388222659917082876715874606081709,
    -0.01918260436019865893946346270175301,
     0.00110789252584784967184098099266118,
    -0.00004157234558247208803840231814601,
     0.00000109278524300228715295578966285,
    -0.00000002123285954183465219601280329,
     0.00000000031733482164348544865129873,
    -0.00000000000376141547987683699381798,
     0.00000000000003622653488483964336956,
    -0.00000000000000028911528493651852433,
     0.00000000000000000194327860676494420,
    -0.00000000000000000001115183182650184,
     0.00000000000000000000005527858887706,
    -0.00000000000000000000000023907013943,
     0.00000000000000000000000000091001612,
    -0.00000000000000000000000000000307233,
     0.00000000000000000000000000000000926 };
  double f;
  double g;
  static int nci = 0;
  double sinx;
  double value;
  static double xsml = 0.0;
  double y;

  if ( nci == 0 )
  {
    nci = r8_inits ( cics, 19, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( r8_mach ( 3 ) );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_CI - Fatal error!\n";
    cerr << "  X <= 0.0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = - 1.0;
    value = log ( x ) - 0.5 + r8_csevl ( y, cics, nci );
  }
  else if ( x <= 4.0 )
  {
    y = ( x * x - 8.0 ) * 0.125;
    value = log ( x ) - 0.5 + r8_csevl ( y, cics, nci );
  }
  else
  {
    r8_sifg ( x, f, g );
    sinx = sin ( x );
    value = f * sinx - g * cos ( x );
  }
  return value;
}
//****************************************************************************80

double r8_cin ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CIN evaluates the alternate cosine integral Cin of an R8 argument.
//
//  Discussion:
//
//    CIN(X) = gamma + log(X) 
//      + integral ( 0 <= T <= X ) ( cos ( T ) - 1 ) / T  dT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_CIN, the cosine integral Cin evaluated at X.
//
{
  double absx;
  static double cincs[18] = {
     0.37074501750909688741654801228564992,
    -0.05893574896364446831956864397363697,
     0.00538189642113569124048745326203340,
    -0.00029860052841962135319594906563410,
     0.00001095572575321620077031054467306,
    -0.00000028405454877346630491727187731,
     0.00000000546973994875384912457861806,
    -0.00000000008124187461318157083277452,
     0.00000000000095868593117706609013181,
    -0.00000000000000920266004392351031377,
     0.00000000000000007325887999017895024,
    -0.00000000000000000049143726675842909,
     0.00000000000000000000281577746753902,
    -0.00000000000000000000001393986788501,
     0.00000000000000000000000006022485646,
    -0.00000000000000000000000000022904717,
     0.00000000000000000000000000000077273,
    -0.00000000000000000000000000000000233 };
  static double eul = 0.57721566490153286060651209008240;
  double f;
  double g;
  static int ncin = 0;
  double sinx;
  double value;
  static double xmin = 0.0;

  if ( ncin == 0 )
  {
    ncin = r8_inits ( cincs, 18, 0.1 * r8_mach ( 3 ) );
    xmin = sqrt ( r8_mach ( 1 ) );
  }

  absx = r8_abs ( x );

  if ( absx <= xmin )
  {
    value = 0.0;
  }
  else if ( absx <= 4.0 )
  {
    value = r8_csevl ( ( x * x - 8.0 ) * 0.125, cincs, ncin ) * x * x;
  }
  else
  {
    r8_sifg ( absx, f, g );
    sinx = sin ( absx );
    value = - f * sinx + g * cos ( absx ) + log ( absx ) + eul;
  }
  return value;
}
//****************************************************************************80

double r8_cinh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CINH: alternate hyperbolic cosine integral Cinh of an R8 argument.
//
//  Discussion:
//
//    Cinh ( x ) = Integral ( 0 <= t <= x ) ( cosh ( t ) - 1 ) dt / t
//
//    The original text of this program had a mistake:
//      y = x * x / 9.0 - 1.0
//    has been corrected to
//      y = x * x / 4.5 - 1.0
//    JVB, 27 March 2010
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_CINH, the hyperbolic cosine integral Cinh
//    evaluated at X.
//
{
  double absx;
  static double cinhcs[16] = {
     0.1093291636520734431407425199795917,
     0.0573928847550379676445323429825108,
     0.0028095756978830353416404208940774,
     0.0000828780840721356655731765069792,
     0.0000016278596173914185577726018815,
     0.0000000227809519255856619859083591,
     0.0000000002384484842463059257284002,
     0.0000000000019360829780781957471028,
     0.0000000000000125453698328172559683,
     0.0000000000000000663637449497262300,
     0.0000000000000000002919639263594744,
     0.0000000000000000000010849123956107,
     0.0000000000000000000000034499080805,
     0.0000000000000000000000000094936664,
     0.0000000000000000000000000000228291,
     0.0000000000000000000000000000000484 };
  static double eul = 0.57721566490153286060651209008240;
  static int ncinh = 0;
  double value;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ncinh == 0 )
  {
    ncinh = r8_inits ( cinhcs, 16, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( r8_mach ( 3 ) );
    xmin = 2.0 * sqrt ( r8_mach ( 1 ) );
  }

  absx = r8_abs ( x );

  if ( x == 0.0 )
  {
    value = 0.0;
  }
  else if ( absx <= xmin )
  {
    value = 0.0;
  }
  else if ( x <= xsml )
  {
    y = - 1.0;
    value = x * x * ( 0.25 + r8_csevl ( y, cinhcs, ncinh ) );
  }
  else if ( x <= 3.0 )
  {
    y = x * x / 4.5 - 1.0;
    value = x * x * ( 0.25 + r8_csevl ( y, cinhcs, ncinh ) );
  }
  else
  {
    value = r8_chi ( absx ) - eul - log ( absx );
  }
  return value;
}
//****************************************************************************80

double r8_cos ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_COS evaluates the cosine of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_COS, the cosine of X.
//
{
  double absx;
  double f;
  int n2;
  static int ntsn = 0;
  double pi2 = 1.57079632679489661923132169163975;
  static double pi2rec = 0.63661977236758134307553505349006;
  static double pihi = 3.140625;
  static double pilo = 9.6765358979323846264338327950288E-04;
  static double pirec = 0.31830988618379067153776752674503;
  static double sincs[15] = {
   -0.374991154955873175839919279977323464,
   -0.181603155237250201863830316158004754,
    0.005804709274598633559427341722857921,
   -0.000086954311779340757113212316353178,
    0.000000754370148088851481006839927030,
   -0.000000004267129665055961107126829906,
    0.000000000016980422945488168181824792,
   -0.000000000000050120578889961870929524,
    0.000000000000000114101026680010675628,
   -0.000000000000000000206437504424783134,
    0.000000000000000000000303969595918706,
   -0.000000000000000000000000371357734157,
    0.000000000000000000000000000382486123,
   -0.000000000000000000000000000000336623,
    0.000000000000000000000000000000000256 };
  double value;
  static double xmax = 0.0;
  double xn;
  static double xsml = 0.0;
  static double xwarn = 0.0;
  double y;

  if ( ntsn == 0 )
  {
    ntsn = r8_inits ( sincs, 15, 0.1 * r8_mach ( 3 ) );
    xsml = r8_sqrt ( 2.0 * r8_mach ( 3 ) );
    xmax = 1.0 / r8_mach ( 4 );
    xwarn = r8_sqrt ( xmax );
  }

  absx = r8_abs ( x );
  y = absx + pi2;

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R8_COS - Warning!\n";
    cerr << "  No precision because |X| is big.\n";
    value = 0.0;
    return value;
  }

  if ( xwarn < y )
  {
    cerr << "\n";
    cerr << "R8_COS - Warning!\n";
    cerr << "  Answer < half precision because |X| is big.\n";
  }

  value = 1.0;

  if ( absx < xsml )
  {
    return value;
  }

  xn = ( double ) ( ( int ) ( y * pirec + 0.5 ) );
  n2 = ( int ) ( r8_mod ( xn, 2.0 ) + 0.5 );
  xn = xn - 0.5;
  f = ( absx - xn * pihi ) - xn * pilo;

  xn = 2.0 * ( f * pi2rec ) * ( f * pi2rec ) - 1.0;
  value = f + f * r8_csevl ( xn, sincs, ntsn );

  if ( n2 != 0 )
  {
    value = - value;
  }

  if ( value < -1.0 )
  {
    value = -1.0;
  }
  else if ( 1.0 < value )
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

double r8_cos_deg ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_COS_DEG evaluates the cosine of an R8 argument in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument in degrees.
//
//    Output, double R8_COS_DEG, the cosine of X.
//
{
  int n;
  static double raddeg = 0.017453292519943295769236907684886;
  double value;

  value = cos ( raddeg * x );

  if ( fmod ( x, 90.0 ) == 0.0 )
  {
    n = ( int ) ( r8_abs ( x ) / 90.0 + 0.5 );
    n = ( n % 2 );

    if ( n == 1 )
    {
      value = 0.0;
    }
    else if ( value < 0.0 )
    {
      value = - 1.0;
    }
    else
    {
      value = + 1.0;
    }
  }
  return value;
}
//****************************************************************************80

double r8_cosh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_COSH evaluates the hyperbolic cosine of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_COSH, the hyperbolic cosine of X.
//
{
  double value;
  double y;
  static double ymax = 0.0;

  if ( ymax == 0.0 )
  {
    ymax = 1.0 / sqrt ( r8_mach ( 3 ) );
  }

  y = exp ( r8_abs ( x ) );

  if ( y < ymax )
  {
    value = 0.5 * ( y + 1.0 / y );
  }
  else
  {
    value = 0.5 * y;
  }
  return value;
}
//****************************************************************************80

double r8_cot ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_COT evaluates the cotangent of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_COT, the cotangent of X.
//
{
  double ainty;
  double ainty2;
  static double cotcs[15] = {
    +0.240259160982956302509553617744970,
    -0.165330316015002278454746025255758E-01,
    -0.429983919317240189356476228239895E-04,
    -0.159283223327541046023490851122445E-06,
    -0.619109313512934872588620579343187E-09,
    -0.243019741507264604331702590579575E-11,
    -0.956093675880008098427062083100000E-14,
    -0.376353798194580580416291539706666E-16,
    -0.148166574646746578852176794666666E-18,
    -0.583335658903666579477984000000000E-21,
    -0.229662646964645773928533333333333E-23,
    -0.904197057307483326719999999999999E-26,
    -0.355988551920600064000000000000000E-28,
    -0.140155139824298666666666666666666E-30,
    -0.551800436872533333333333333333333E-33 };
  int ifn;
  static int nterms = 0;
  static double pi2rec = 0.011619772367581343075535053490057;
  double prodbg;
  static double sqeps = 0.0;
  double value;
  static double xmax = 0.0;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;
  double yrem;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( cotcs, 15, 0.1 * r8_mach ( 3 ) );
    xmax = 1.0 / r8_mach ( 4 );
    xsml = sqrt ( 3.0 * r8_mach ( 3 ) );
    xmin = exp ( r8_max ( log ( r8_mach ( 1 ) ), 
      - log ( r8_mach ( 2 ) ) )  + 0.01 );
    sqeps = sqrt ( r8_mach ( 4 ) );
  }

  y = r8_abs ( x );

  if ( y < xmin )
  {
    cerr << "\n";
    cerr << "R8_COT - Fatal error!\n";
    cerr << "  |X| is too small.\n";
    exit ( 1 );
  }

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R8_COT - Fatal error!\n";
    cerr << "  |X| is too big.\n";
    exit ( 1 );
  }
//
//  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
//  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
//  = aint(.625*y) + aint(z) + rem(z)
//
  ainty = r8_aint ( y );
  yrem = y - ainty;
  prodbg = 0.625 * ainty;
  ainty = r8_aint ( prodbg );
  y = ( prodbg - ainty ) + 0.625 * yrem + y * pi2rec;
  ainty2 = r8_aint ( y );
  ainty = ainty + ainty2;
  y = y - ainty2;

  ifn = ( int ) fmod ( ainty, 2.0 );
  if ( ifn == 1 )
  {
    y = 1.0 - y;
  }

  if ( 0.5 < r8_abs ( x ) && y < r8_abs ( x ) * sqeps )
  {
    cerr << "\n";
    cerr << "R8_COT - Warning!\n";
    cerr << "  Answer less than half precision.\n";
    cerr << "  |X| too big, or X nearly a nonzero multiple of pi.\n";
    exit ( 1 );
  }

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R8_COT - Fatal error!\n";
    cerr << "  X is a multiple of pi.\n";
    exit ( 1 );
  }
  else if ( y <= xsml )
  {
    value = 1.0 / y;
  }
  else if ( y <= 0.25 )
  {
    value = ( 0.5 + r8_csevl ( 32.0 * y * y - 1.0, cotcs, nterms ) ) / y;
  }
  else if ( y <= 0.5 )
  {
    value = ( 0.5 + r8_csevl ( 8.0 * y * y - 1.0, 
      cotcs, nterms ) ) / ( 0.5 * y );

    value = ( value * value - 1.0 ) * 0.5 / value;
  }
  else
  {
    value = ( 0.5 + r8_csevl ( 2.0 * y * y - 1.0, 
      cotcs, nterms ) ) / ( 0.25 * y );
    value = ( value * value - 1.0 ) * 0.5 / value;
    value = ( value * value - 1.0 ) * 0.5 / value;
  }

  if ( x < 0.0 )
  {
    value = - r8_abs ( value );
  }
  else
  {
    value = + r8_abs ( value );
  }

  if ( ifn == 1 )
  {
    value = - value;
  }
  return value;
}
//****************************************************************************80

 double r8_csevl ( double x, double a[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CSEVL evaluates a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double CS[N], the Chebyshev coefficients.
//
//    Input, int N, the number of Chebyshev coefficients.
//
//    Output, double R8_CSEVL, the Chebyshev series evaluated at X.
//
{
  double b0;
  double b1;
  double b2;
  int i;
  double twox;
  double value;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R8_CSEVL - Fatal error!\n";
    cerr << "  Number of terms <= 0.\n";
    exit ( 1 );
  }

  if ( 1000 < n )
  {
    cerr << "\n";
    cerr << "R8_CSEVL - Fatal error!\n";
    cerr << "  Number of terms greater than 1000.\n";
    exit ( 1 );
 }

  if ( x < -1.1 || 1.1 < x )
  {
    cerr << "\n";
    cerr << "R8_CSEVL - Fatal error!\n";
    cerr << "  X outside (-1,+1).\n";
    exit ( 1 );
  }

  twox = 2.0 * x;
  b1 = 0.0;
  b0 = 0.0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[i];
  }

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double r8_dawson ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_DAWSON evaluates Dawson's integral of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_DAWSON, the value of Dawson's integral at X.
//
{
  static double daw2cs[45] = {
    -0.56886544105215527114160533733674E-01,
    -0.31811346996168131279322878048822,
    +0.20873845413642236789741580198858,
    -0.12475409913779131214073498314784,
    +0.67869305186676777092847516423676E-01,
    -0.33659144895270939503068230966587E-01,
    +0.15260781271987971743682460381640E-01,
    -0.63483709625962148230586094788535E-02,
    +0.24326740920748520596865966109343E-02,
    -0.86219541491065032038526983549637E-03,
    +0.28376573336321625302857636538295E-03,
    -0.87057549874170423699396581464335E-04,
    +0.24986849985481658331800044137276E-04,
    -0.67319286764160294344603050339520E-05,
    +0.17078578785573543710504524047844E-05,
    -0.40917551226475381271896592490038E-06,
    +0.92828292216755773260751785312273E-07,
    -0.19991403610147617829845096332198E-07,
    +0.40963490644082195241210487868917E-08,
    -0.80032409540993168075706781753561E-09,
    +0.14938503128761465059143225550110E-09,
    -0.26687999885622329284924651063339E-10,
    +0.45712216985159458151405617724103E-11,
    -0.75187305222043565872243727326771E-12,
    +0.11893100052629681879029828987302E-12,
    -0.18116907933852346973490318263084E-13,
    +0.26611733684358969193001612199626E-14,
    -0.37738863052129419795444109905930E-15,
    +0.51727953789087172679680082229329E-16,
    -0.68603684084077500979419564670102E-17,
    +0.88123751354161071806469337321745E-18,
    -0.10974248249996606292106299624652E-18,
    +0.13261199326367178513595545891635E-19,
    -0.15562732768137380785488776571562E-20,
    +0.17751425583655720607833415570773E-21,
    -0.19695006967006578384953608765439E-22,
    +0.21270074896998699661924010120533E-23,
    -0.22375398124627973794182113962666E-24,
    +0.22942768578582348946971383125333E-25,
    -0.22943788846552928693329592319999E-26,
    +0.22391702100592453618342297600000E-27,
    -0.21338230616608897703678225066666E-28,
    +0.19866196585123531518028458666666E-29,
    -0.18079295866694391771955199999999E-30,
    +0.16090686015283030305450666666666E-31 };
  static double dawacs[75] = {
    +0.1690485637765703755422637438849E-01,
    +0.8683252278406957990536107850768E-02,
    +0.2424864042417715453277703459889E-03,
    +0.1261182399572690001651949240377E-04,
    +0.1066453314636176955705691125906E-05,
    +0.1358159794790727611348424505728E-06,
    +0.2171042356577298398904312744743E-07,
    +0.2867010501805295270343676804813E-08,
    -0.1901336393035820112282492378024E-09,
    -0.3097780484395201125532065774268E-09,
    -0.1029414876057509247398132286413E-09,
    -0.6260356459459576150417587283121E-11,
    +0.8563132497446451216262303166276E-11,
    +0.3033045148075659292976266276257E-11,
    -0.2523618306809291372630886938826E-12,
    -0.4210604795440664513175461934510E-12,
    -0.4431140826646238312143429452036E-13,
    +0.4911210272841205205940037065117E-13,
    +0.1235856242283903407076477954739E-13,
    -0.5788733199016569246955765071069E-14,
    -0.2282723294807358620978183957030E-14,
    +0.7637149411014126476312362917590E-15,
    +0.3851546883566811728777594002095E-15,
    -0.1199932056928290592803237283045E-15,
    -0.6313439150094572347334270285250E-16,
    +0.2239559965972975375254912790237E-16,
    +0.9987925830076495995132891200749E-17,
    -0.4681068274322495334536246507252E-17,
    -0.1436303644349721337241628751534E-17,
    +0.1020822731410541112977908032130E-17,
    +0.1538908873136092072837389822372E-18,
    -0.2189157877645793888894790926056E-18,
    +0.2156879197938651750392359152517E-20,
    +0.4370219827442449851134792557395E-19,
    -0.8234581460977207241098927905177E-20,
    -0.7498648721256466222903202835420E-20,
    +0.3282536720735671610957612930039E-20,
    +0.8858064309503921116076561515151E-21,
    -0.9185087111727002988094460531485E-21,
    +0.2978962223788748988314166045791E-22,
    +0.1972132136618471883159505468041E-21,
    -0.5974775596362906638089584995117E-22,
    -0.2834410031503850965443825182441E-22,
    +0.2209560791131554514777150489012E-22,
    -0.5439955741897144300079480307711E-25,
    -0.5213549243294848668017136696470E-23,
    +0.1702350556813114199065671499076E-23,
    +0.6917400860836148343022185660197E-24,
    -0.6540941793002752512239445125802E-24,
    +0.6093576580439328960371824654636E-25,
    +0.1408070432905187461501945080272E-24,
    -0.6785886121054846331167674943755E-25,
    -0.9799732036214295711741583102225E-26,
    +0.2121244903099041332598960939160E-25,
    -0.5954455022548790938238802154487E-26,
    -0.3093088861875470177838847232049E-26,
    +0.2854389216344524682400691986104E-26,
    -0.3951289447379305566023477271811E-27,
    -0.5906000648607628478116840894453E-27,
    +0.3670236964668687003647889980609E-27,
    -0.4839958238042276256598303038941E-29,
    -0.9799265984210443869597404017022E-28,
    +0.4684773732612130606158908804300E-28,
    +0.5030877696993461051647667603155E-29,
    -0.1547395051706028239247552068295E-28,
    +0.6112180185086419243976005662714E-29,
    +0.1357913399124811650343602736158E-29,
    -0.2417687752768673088385304299044E-29,
    +0.8369074582074298945292887587291E-30,
    +0.2665413042788979165838319401566E-30,
    -0.3811653692354890336935691003712E-30,
    +0.1230054721884951464371706872585E-30,
    +0.4622506399041493508805536929983E-31,
    -0.6120087296881677722911435593001E-31,
    +0.1966024640193164686956230217896E-31 };
  static double dawcs[21] = {
    -0.6351734375145949201065127736293E-02,
    -0.2294071479677386939899824125866,
    +0.2213050093908476441683979161786E-01,
    -0.1549265453892985046743057753375E-02,
    +0.8497327715684917456777542948066E-04,
    -0.3828266270972014924994099521309E-05,
    +0.1462854806250163197757148949539E-06,
    -0.4851982381825991798846715425114E-08,
    +0.1421463577759139790347568183304E-09,
    -0.3728836087920596525335493054088E-11,
    +0.8854942961778203370194565231369E-13,
    -0.1920757131350206355421648417493E-14,
    +0.3834325867246327588241074439253E-16,
    -0.7089154168175881633584099327999E-18,
    +0.1220552135889457674416901120000E-19,
    -0.1966204826605348760299451733333E-21,
    +0.2975845541376597189113173333333E-23,
    -0.4247069514800596951039999999999E-25,
    +0.5734270767391742798506666666666E-27,
    -0.7345836823178450261333333333333E-29,
    +0.8951937667516552533333333333333E-31 };
  double eps;
  static int ntdaw = 0;
  static int ntdaw2 = 0;
  static int ntdawa = 0;
  double value;
  static double xbig = 0.0;
  static double xmax = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ntdaw == 0 )
  {
    eps = r8_mach ( 3 );
    ntdaw  = r8_inits ( dawcs,  21, 0.1 * eps );
    ntdaw2 = r8_inits ( daw2cs, 45, 0.1 * eps );
    ntdawa = r8_inits ( dawacs, 75, 0.1 * eps );
    xsml = sqrt ( 1.5 * eps );
    xbig = sqrt ( 0.5 / eps );
    xmax = exp ( r8_min ( - log ( 2.0 * r8_mach ( 1 ) ), 
      log ( r8_mach ( 2 ) ) ) - 0.01 );
  }

  y = r8_abs ( x );

  if ( y <= xsml )
  {
    value = x;
  }
  else if ( y <= 1.0 )
  {
    value = x * ( 0.75 + r8_csevl ( 2.0 * y * y - 1.0, dawcs, ntdaw ) );
  }
  else if ( y <= 4.0 )
  {
    value = x * ( 0.25 + r8_csevl ( 0.125 * y * y - 1.0, daw2cs, ntdaw2 ) );
  }
  else if ( y < xbig )
  {
    value = ( 0.5 + r8_csevl ( 32.0 / y / y - 1.0, dawacs, ntdawa ) ) / x;
  }
  else if ( y <= xmax )
  {
    value = 0.5 / x;
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double r8_e1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_E1 evaluates the exponential integral E1 for an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_E1, the exponential integral E1 evaluated at X.
//
{
  static double ae10cs[50] = {
    +0.3284394579616699087873844201881E-01,
    -0.1669920452031362851476184343387E-01,
    +0.2845284724361346807424899853252E-03,
    -0.7563944358516206489487866938533E-05,
    +0.2798971289450859157504843180879E-06,
    -0.1357901828534531069525563926255E-07,
    +0.8343596202040469255856102904906E-09,
    -0.6370971727640248438275242988532E-10,
    +0.6007247608811861235760831561584E-11,
    -0.7022876174679773590750626150088E-12,
    +0.1018302673703687693096652346883E-12,
    -0.1761812903430880040406309966422E-13,
    +0.3250828614235360694244030353877E-14,
    -0.5071770025505818678824872259044E-15,
    +0.1665177387043294298172486084156E-16,
    +0.3166753890797514400677003536555E-16,
    -0.1588403763664141515133118343538E-16,
    +0.4175513256138018833003034618484E-17,
    -0.2892347749707141906710714478852E-18,
    -0.2800625903396608103506340589669E-18,
    +0.1322938639539270903707580023781E-18,
    -0.1804447444177301627283887833557E-19,
    -0.7905384086522616076291644817604E-20,
    +0.4435711366369570103946235838027E-20,
    -0.4264103994978120868865309206555E-21,
    -0.3920101766937117541553713162048E-21,
    +0.1527378051343994266343752326971E-21,
    +0.1024849527049372339310308783117E-22,
    -0.2134907874771433576262711405882E-22,
    +0.3239139475160028267061694700366E-23,
    +0.2142183762299889954762643168296E-23,
    -0.8234609419601018414700348082312E-24,
    -0.1524652829645809479613694401140E-24,
    +0.1378208282460639134668480364325E-24,
    +0.2131311202833947879523224999253E-26,
    -0.2012649651526484121817466763127E-25,
    +0.1995535662263358016106311782673E-26,
    +0.2798995808984003464948686520319E-26,
    -0.5534511845389626637640819277823E-27,
    -0.3884995396159968861682544026146E-27,
    +0.1121304434507359382850680354679E-27,
    +0.5566568152423740948256563833514E-28,
    -0.2045482929810499700448533938176E-28,
    -0.8453813992712336233411457493674E-29,
    +0.3565758433431291562816111116287E-29,
    +0.1383653872125634705539949098871E-29,
    -0.6062167864451372436584533764778E-30,
    -0.2447198043989313267437655119189E-30,
    +0.1006850640933998348011548180480E-30,
    +0.4623685555014869015664341461674E-31 };
  static double ae11cs[60] = {
    +0.20263150647078889499401236517381,
    -0.73655140991203130439536898728034E-01,
    +0.63909349118361915862753283840020E-02,
    -0.60797252705247911780653153363999E-03,
    -0.73706498620176629330681411493484E-04,
    +0.48732857449450183453464992488076E-04,
    -0.23837064840448290766588489460235E-05,
    -0.30518612628561521027027332246121E-05,
    +0.17050331572564559009688032992907E-06,
    +0.23834204527487747258601598136403E-06,
    +0.10781772556163166562596872364020E-07,
    -0.17955692847399102653642691446599E-07,
    -0.41284072341950457727912394640436E-08,
    +0.68622148588631968618346844526664E-09,
    +0.53130183120506356147602009675961E-09,
    +0.78796880261490694831305022893515E-10,
    -0.26261762329356522290341675271232E-10,
    -0.15483687636308261963125756294100E-10,
    -0.25818962377261390492802405122591E-11,
    +0.59542879191591072658903529959352E-12,
    +0.46451400387681525833784919321405E-12,
    +0.11557855023255861496288006203731E-12,
    -0.10475236870835799012317547189670E-14,
    -0.11896653502709004368104489260929E-13,
    -0.47749077490261778752643019349950E-14,
    -0.81077649615772777976249734754135E-15,
    +0.13435569250031554199376987998178E-15,
    +0.14134530022913106260248873881287E-15,
    +0.49451592573953173115520663232883E-16,
    +0.79884048480080665648858587399367E-17,
    -0.14008632188089809829248711935393E-17,
    -0.14814246958417372107722804001680E-17,
    -0.55826173646025601904010693937113E-18,
    -0.11442074542191647264783072544598E-18,
    +0.25371823879566853500524018479923E-20,
    +0.13205328154805359813278863389097E-19,
    +0.62930261081586809166287426789485E-20,
    +0.17688270424882713734999261332548E-20,
    +0.23266187985146045209674296887432E-21,
    -0.67803060811125233043773831844113E-22,
    -0.59440876959676373802874150531891E-22,
    -0.23618214531184415968532592503466E-22,
    -0.60214499724601478214168478744576E-23,
    -0.65517906474348299071370444144639E-24,
    +0.29388755297497724587042038699349E-24,
    +0.22601606200642115173215728758510E-24,
    +0.89534369245958628745091206873087E-25,
    +0.24015923471098457555772067457706E-25,
    +0.34118376888907172955666423043413E-26,
    -0.71617071694630342052355013345279E-27,
    -0.75620390659281725157928651980799E-27,
    -0.33774612157467324637952920780800E-27,
    -0.10479325703300941711526430332245E-27,
    -0.21654550252170342240854880201386E-28,
    -0.75297125745288269994689298432000E-30,
    +0.19103179392798935768638084000426E-29,
    +0.11492104966530338547790728833706E-29,
    +0.43896970582661751514410359193600E-30,
    +0.12320883239205686471647157725866E-30,
    +0.22220174457553175317538581162666E-31 };
  static double ae12cs[41] = {
    +0.63629589796747038767129887806803,
    -0.13081168675067634385812671121135,
    -0.84367410213053930014487662129752E-02,
    +0.26568491531006685413029428068906E-02,
    +0.32822721781658133778792170142517E-03,
    -0.23783447771430248269579807851050E-04,
    -0.11439804308100055514447076797047E-04,
    -0.14405943433238338455239717699323E-05,
    +0.52415956651148829963772818061664E-08,
    +0.38407306407844323480979203059716E-07,
    +0.85880244860267195879660515759344E-08,
    +0.10219226625855003286339969553911E-08,
    +0.21749132323289724542821339805992E-10,
    -0.22090238142623144809523503811741E-10,
    -0.63457533544928753294383622208801E-11,
    -0.10837746566857661115340539732919E-11,
    -0.11909822872222586730262200440277E-12,
    -0.28438682389265590299508766008661E-14,
    +0.25080327026686769668587195487546E-14,
    +0.78729641528559842431597726421265E-15,
    +0.15475066347785217148484334637329E-15,
    +0.22575322831665075055272608197290E-16,
    +0.22233352867266608760281380836693E-17,
    +0.16967819563544153513464194662399E-19,
    -0.57608316255947682105310087304533E-19,
    -0.17591235774646878055625369408853E-19,
    -0.36286056375103174394755328682666E-20,
    -0.59235569797328991652558143488000E-21,
    -0.76030380926310191114429136895999E-22,
    -0.62547843521711763842641428479999E-23,
    +0.25483360759307648606037606400000E-24,
    +0.25598615731739857020168874666666E-24,
    +0.71376239357899318800207052800000E-25,
    +0.14703759939567568181578956800000E-25,
    +0.25105524765386733555198634666666E-26,
    +0.35886666387790890886583637333333E-27,
    +0.39886035156771301763317759999999E-28,
    +0.21763676947356220478805333333333E-29,
    -0.46146998487618942367607466666666E-30,
    -0.20713517877481987707153066666666E-30,
    -0.51890378563534371596970666666666E-31 };
  static double ae13cs[50] = {
    -0.60577324664060345999319382737747,
    -0.11253524348366090030649768852718,
    +0.13432266247902779492487859329414E-01,
    -0.19268451873811457249246838991303E-02,
    +0.30911833772060318335586737475368E-03,
    -0.53564132129618418776393559795147E-04,
    +0.98278128802474923952491882717237E-05,
    -0.18853689849165182826902891938910E-05,
    +0.37494319356894735406964042190531E-06,
    -0.76823455870552639273733465680556E-07,
    +0.16143270567198777552956300060868E-07,
    -0.34668022114907354566309060226027E-08,
    +0.75875420919036277572889747054114E-09,
    -0.16886433329881412573514526636703E-09,
    +0.38145706749552265682804250927272E-10,
    -0.87330266324446292706851718272334E-11,
    +0.20236728645867960961794311064330E-11,
    -0.47413283039555834655210340820160E-12,
    +0.11221172048389864324731799928920E-12,
    -0.26804225434840309912826809093395E-13,
    +0.64578514417716530343580369067212E-14,
    -0.15682760501666478830305702849194E-14,
    +0.38367865399315404861821516441408E-15,
    -0.94517173027579130478871048932556E-16,
    +0.23434812288949573293896666439133E-16,
    -0.58458661580214714576123194419882E-17,
    +0.14666229867947778605873617419195E-17,
    -0.36993923476444472706592538274474E-18,
    +0.93790159936721242136014291817813E-19,
    -0.23893673221937873136308224087381E-19,
    +0.61150624629497608051934223837866E-20,
    -0.15718585327554025507719853288106E-20,
    +0.40572387285585397769519294491306E-21,
    -0.10514026554738034990566367122773E-21,
    +0.27349664930638667785806003131733E-22,
    -0.71401604080205796099355574271999E-23,
    +0.18705552432235079986756924211199E-23,
    -0.49167468166870480520478020949333E-24,
    +0.12964988119684031730916087125333E-24,
    -0.34292515688362864461623940437333E-25,
    +0.90972241643887034329104820906666E-26,
    -0.24202112314316856489934847999999E-26,
    +0.64563612934639510757670475093333E-27,
    -0.17269132735340541122315987626666E-27,
    +0.46308611659151500715194231466666E-28,
    -0.12448703637214131241755170133333E-28,
    +0.33544574090520678532907007999999E-29,
    -0.90598868521070774437543935999999E-30,
    +0.24524147051474238587273216000000E-30,
    -0.66528178733552062817107967999999E-31 };
  static double ae14cs[64] = {
    -0.1892918000753016825495679942820,
    -0.8648117855259871489968817056824E-01,
    +0.7224101543746594747021514839184E-02,
    -0.8097559457557386197159655610181E-03,
    +0.1099913443266138867179251157002E-03,
    -0.1717332998937767371495358814487E-04,
    +0.2985627514479283322825342495003E-05,
    -0.5659649145771930056560167267155E-06,
    +0.1152680839714140019226583501663E-06,
    -0.2495030440269338228842128765065E-07,
    +0.5692324201833754367039370368140E-08,
    -0.1359957664805600338490030939176E-08,
    +0.3384662888760884590184512925859E-09,
    -0.8737853904474681952350849316580E-10,
    +0.2331588663222659718612613400470E-10,
    -0.6411481049213785969753165196326E-11,
    +0.1812246980204816433384359484682E-11,
    -0.5253831761558460688819403840466E-12,
    +0.1559218272591925698855028609825E-12,
    -0.4729168297080398718476429369466E-13,
    +0.1463761864393243502076199493808E-13,
    -0.4617388988712924102232173623604E-14,
    +0.1482710348289369323789239660371E-14,
    -0.4841672496239229146973165734417E-15,
    +0.1606215575700290408116571966188E-15,
    -0.5408917538957170947895023784252E-16,
    +0.1847470159346897881370231402310E-16,
    -0.6395830792759094470500610425050E-17,
    +0.2242780721699759457250233276170E-17,
    -0.7961369173983947552744555308646E-18,
    +0.2859308111540197459808619929272E-18,
    -0.1038450244701137145900697137446E-18,
    +0.3812040607097975780866841008319E-19,
    -0.1413795417717200768717562723696E-19,
    +0.5295367865182740958305442594815E-20,
    -0.2002264245026825902137211131439E-20,
    +0.7640262751275196014736848610918E-21,
    -0.2941119006868787883311263523362E-21,
    +0.1141823539078927193037691483586E-21,
    -0.4469308475955298425247020718489E-22,
    +0.1763262410571750770630491408520E-22,
    -0.7009968187925902356351518262340E-23,
    +0.2807573556558378922287757507515E-23,
    -0.1132560944981086432141888891562E-23,
    +0.4600574684375017946156764233727E-24,
    -0.1881448598976133459864609148108E-24,
    +0.7744916111507730845444328478037E-25,
    -0.3208512760585368926702703826261E-25,
    +0.1337445542910839760619930421384E-25,
    -0.5608671881802217048894771735210E-26,
    +0.2365839716528537483710069473279E-26,
    -0.1003656195025305334065834526856E-26,
    +0.4281490878094161131286642556927E-27,
    -0.1836345261815318199691326958250E-27,
    +0.7917798231349540000097468678144E-28,
    -0.3431542358742220361025015775231E-28,
    +0.1494705493897103237475066008917E-28,
    -0.6542620279865705439739042420053E-29,
    +0.2877581395199171114340487353685E-29,
    -0.1271557211796024711027981200042E-29,
    +0.5644615555648722522388044622506E-30,
    -0.2516994994284095106080616830293E-30,
    +0.1127259818927510206370368804181E-30,
    -0.5069814875800460855562584719360E-31 };
  static double e11cs[29] = {
    -0.16113461655571494025720663927566180E+02,
    +0.77940727787426802769272245891741497E+01,
    -0.19554058188631419507127283812814491E+01,
    +0.37337293866277945611517190865690209,
    -0.56925031910929019385263892220051166E-01,
    +0.72110777696600918537847724812635813E-02,
    -0.78104901449841593997715184089064148E-03,
    +0.73880933562621681878974881366177858E-04,
    -0.62028618758082045134358133607909712E-05,
    +0.46816002303176735524405823868362657E-06,
    -0.32092888533298649524072553027228719E-07,
    +0.20151997487404533394826262213019548E-08,
    -0.11673686816697793105356271695015419E-09,
    +0.62762706672039943397788748379615573E-11,
    -0.31481541672275441045246781802393600E-12,
    +0.14799041744493474210894472251733333E-13,
    -0.65457091583979673774263401588053333E-15,
    +0.27336872223137291142508012748799999E-16,
    -0.10813524349754406876721727624533333E-17,
    +0.40628328040434303295300348586666666E-19,
    -0.14535539358960455858914372266666666E-20,
    +0.49632746181648636830198442666666666E-22,
    -0.16208612696636044604866560000000000E-23,
    +0.50721448038607422226431999999999999E-25,
    -0.15235811133372207813973333333333333E-26,
    +0.44001511256103618696533333333333333E-28,
    -0.12236141945416231594666666666666666E-29,
    +0.32809216661066001066666666666666666E-31,
    -0.84933452268306432000000000000000000E-33 };
  static double e12cs[25] = {
    -0.3739021479220279511668698204827E-01,
    +0.4272398606220957726049179176528E-01,
    -0.130318207984970054415392055219726,
    +0.144191240246988907341095893982137E-01,
    -0.134617078051068022116121527983553E-02,
    +0.107310292530637799976115850970073E-03,
    -0.742999951611943649610283062223163E-05,
    +0.453773256907537139386383211511827E-06,
    -0.247641721139060131846547423802912E-07,
    +0.122076581374590953700228167846102E-08,
    -0.548514148064092393821357398028261E-10,
    +0.226362142130078799293688162377002E-11,
    -0.863589727169800979404172916282240E-13,
    +0.306291553669332997581032894881279E-14,
    -0.101485718855944147557128906734933E-15,
    +0.315482174034069877546855328426666E-17,
    -0.923604240769240954484015923200000E-19,
    +0.255504267970814002440435029333333E-20,
    -0.669912805684566847217882453333333E-22,
    +0.166925405435387319431987199999999E-23,
    -0.396254925184379641856000000000000E-25,
    +0.898135896598511332010666666666666E-27,
    -0.194763366993016433322666666666666E-28,
    +0.404836019024630033066666666666666E-30,
    -0.807981567699845120000000000000000E-32 };
  double eta;
  static int ntae10 = 0;
  static int ntae11 = 0;
  static int ntae12 = 0;
  static int ntae13 = 0;
  static int ntae14 = 0;
  static int nte11 = 0;
  static int nte12 = 0;
  double value;
  static double xmax = 0.0;

  if ( ntae10 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    ntae10 = r8_inits ( ae10cs, 50, eta );
    ntae11 = r8_inits ( ae11cs, 60, eta );
    ntae12 = r8_inits ( ae12cs, 41, eta );
    nte11 = r8_inits ( e11cs, 29, eta );
    nte12 = r8_inits ( e12cs, 25, eta );
    ntae13 = r8_inits ( ae13cs, 50, eta );
    ntae14 = r8_inits ( ae14cs, 64, eta );
    xmax = - log ( r8_mach ( 1 ) );
    xmax = xmax - log ( xmax );
  }

  if ( x <= - 32.0 )
  {
    value = exp ( - x ) / x * ( 1.0 
      + r8_csevl ( 64.0 / x + 1.0, ae10cs, ntae10 ) );
  }
  else if ( x <= - 8.0 )
  {
    value = exp ( - x ) / x * ( 1.0 
      + r8_csevl ( ( 64.0 / x + 5.0 ) / 3.0, ae11cs, ntae11 ) );
  }
  else if ( x <= - 4.0 )
  {
    value = exp ( - x ) / x * (1.0 
      + r8_csevl ( 16.0 / x + 3.0, ae12cs, ntae12 ) );
  }
  else if ( x <= - 1.0 )
  {
    value = - log ( - x ) 
      + r8_csevl ( ( 2.0 * x + 5.0 ) / 3.0, e11cs, nte11 );
  }
  else if ( x == 0.0 )
  {
    cerr << "\n";
    cerr << "R8_E1 - Fatal error!\n";
    cerr << "  X is zero.\n";
    exit ( 1 );
  }
  else if ( x <= 1.0 )
  {
    value = ( - log ( r8_abs ( x ) ) - 0.6875 + x ) 
      + r8_csevl ( x, e12cs, nte12 );
  }
  else if ( x <= 4.0 )
  {
    value = exp ( - x ) / x * ( 1.0 
      + r8_csevl ( ( 8.0 / x - 5.0 ) / 3.0, ae13cs, ntae13 ) );
  }
  else if ( x <= xmax )
  {
    value = exp ( - x ) / x * ( 1.0 
      + r8_csevl ( 8.0 / x - 1.0, ae14cs, ntae14 ) );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double r8_ei ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EI evaluates the exponential integral Ei for an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_EI, the exponential integral Ei evaluated at X.
//
{
  double value;

  value = - r8_e1 ( - x );

  return value;
}
//****************************************************************************80

double r8_erf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERF evaluates the error function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ERF, the error function of X.
//
{
  static double erfcs[21] = {
    -0.49046121234691808039984544033376E-01,
    -0.14226120510371364237824741899631,
    +0.10035582187599795575754676712933E-01,
    -0.57687646997674847650827025509167E-03,
    +0.27419931252196061034422160791471E-04,
    -0.11043175507344507604135381295905E-05,
    +0.38488755420345036949961311498174E-07,
    -0.11808582533875466969631751801581E-08,
    +0.32334215826050909646402930953354E-10,
    -0.79910159470045487581607374708595E-12,
    +0.17990725113961455611967245486634E-13,
    -0.37186354878186926382316828209493E-15,
    +0.71035990037142529711689908394666E-17,
    -0.12612455119155225832495424853333E-18,
    +0.20916406941769294369170500266666E-20,
    -0.32539731029314072982364160000000E-22,
    +0.47668672097976748332373333333333E-24,
    -0.65980120782851343155199999999999E-26,
    +0.86550114699637626197333333333333E-28,
    -0.10788925177498064213333333333333E-29,
    +0.12811883993017002666666666666666E-31 };
  static int nterf = 0;
  static double sqeps = 0.0;
  static double sqrtpi = 1.77245385090551602729816748334115;
  double value;
  static double xbig = 0.0;
  double y;

  if ( nterf == 0 )
  {
    nterf = r8_inits ( erfcs, 21, 0.1 * r8_mach ( 3 ) );
    xbig = sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) );
    sqeps = sqrt ( 2.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= sqeps )
  {
    value = 2.0 * x / sqrtpi;
  }
  else if ( y <= 1.0 )
  {
    value = x * ( 1.0 + r8_csevl ( 2.0 * x * x - 1.0, erfcs, nterf ) );
  }
  else if ( y <= xbig )
  {
    value = 1.0 - r8_erfc ( y );
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  else
  {
    value = 1.0;
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

double r8_erfc ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERFC evaluates the co-error function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ERFC, the co-error function of X.
//
{
  static double erc2cs[49] = {
    -0.6960134660230950112739150826197E-01,
    -0.4110133936262089348982212084666E-01,
    +0.3914495866689626881561143705244E-02,
    -0.4906395650548979161280935450774E-03,
    +0.7157479001377036380760894141825E-04,
    -0.1153071634131232833808232847912E-04,
    +0.1994670590201997635052314867709E-05,
    -0.3642666471599222873936118430711E-06,
    +0.6944372610005012589931277214633E-07,
    -0.1371220902104366019534605141210E-07,
    +0.2788389661007137131963860348087E-08,
    -0.5814164724331161551864791050316E-09,
    +0.1238920491752753181180168817950E-09,
    -0.2690639145306743432390424937889E-10,
    +0.5942614350847910982444709683840E-11,
    -0.1332386735758119579287754420570E-11,
    +0.3028046806177132017173697243304E-12,
    -0.6966648814941032588795867588954E-13,
    +0.1620854541053922969812893227628E-13,
    -0.3809934465250491999876913057729E-14,
    +0.9040487815978831149368971012975E-15,
    -0.2164006195089607347809812047003E-15,
    +0.5222102233995854984607980244172E-16,
    -0.1269729602364555336372415527780E-16,
    +0.3109145504276197583836227412951E-17,
    -0.7663762920320385524009566714811E-18,
    +0.1900819251362745202536929733290E-18,
    -0.4742207279069039545225655999965E-19,
    +0.1189649200076528382880683078451E-19,
    -0.3000035590325780256845271313066E-20,
    +0.7602993453043246173019385277098E-21,
    -0.1935909447606872881569811049130E-21,
    +0.4951399124773337881000042386773E-22,
    -0.1271807481336371879608621989888E-22,
    +0.3280049600469513043315841652053E-23,
    -0.8492320176822896568924792422399E-24,
    +0.2206917892807560223519879987199E-24,
    -0.5755617245696528498312819507199E-25,
    +0.1506191533639234250354144051199E-25,
    -0.3954502959018796953104285695999E-26,
    +0.1041529704151500979984645051733E-26,
    -0.2751487795278765079450178901333E-27,
    +0.7290058205497557408997703680000E-28,
    -0.1936939645915947804077501098666E-28,
    +0.5160357112051487298370054826666E-29,
    -0.1378419322193094099389644800000E-29,
    +0.3691326793107069042251093333333E-30,
    -0.9909389590624365420653226666666E-31,
    +0.2666491705195388413323946666666E-31 };
  static double erfccs[59] = {
    +0.715179310202924774503697709496E-01,
    -0.265324343376067157558893386681E-01,
    +0.171115397792085588332699194606E-02,
    -0.163751663458517884163746404749E-03,
    +0.198712935005520364995974806758E-04,
    -0.284371241276655508750175183152E-05,
    +0.460616130896313036969379968464E-06,
    -0.822775302587920842057766536366E-07,
    +0.159214187277090112989358340826E-07,
    -0.329507136225284321486631665072E-08,
    +0.722343976040055546581261153890E-09,
    -0.166485581339872959344695966886E-09,
    +0.401039258823766482077671768814E-10,
    -0.100481621442573113272170176283E-10,
    +0.260827591330033380859341009439E-11,
    -0.699111056040402486557697812476E-12,
    +0.192949233326170708624205749803E-12,
    -0.547013118875433106490125085271E-13,
    +0.158966330976269744839084032762E-13,
    -0.472689398019755483920369584290E-14,
    +0.143587337678498478672873997840E-14,
    -0.444951056181735839417250062829E-15,
    +0.140481088476823343737305537466E-15,
    -0.451381838776421089625963281623E-16,
    +0.147452154104513307787018713262E-16,
    -0.489262140694577615436841552532E-17,
    +0.164761214141064673895301522827E-17,
    -0.562681717632940809299928521323E-18,
    +0.194744338223207851429197867821E-18,
    -0.682630564294842072956664144723E-19,
    +0.242198888729864924018301125438E-19,
    -0.869341413350307042563800861857E-20,
    +0.315518034622808557122363401262E-20,
    -0.115737232404960874261239486742E-20,
    +0.428894716160565394623737097442E-21,
    -0.160503074205761685005737770964E-21,
    +0.606329875745380264495069923027E-22,
    -0.231140425169795849098840801367E-22,
    +0.888877854066188552554702955697E-23,
    -0.344726057665137652230718495566E-23,
    +0.134786546020696506827582774181E-23,
    -0.531179407112502173645873201807E-24,
    +0.210934105861978316828954734537E-24,
    -0.843836558792378911598133256738E-25,
    +0.339998252494520890627359576337E-25,
    -0.137945238807324209002238377110E-25,
    +0.563449031183325261513392634811E-26,
    -0.231649043447706544823427752700E-26,
    +0.958446284460181015263158381226E-27,
    -0.399072288033010972624224850193E-27,
    +0.167212922594447736017228709669E-27,
    -0.704599152276601385638803782587E-28,
    +0.297976840286420635412357989444E-28,
    -0.126252246646061929722422632994E-28,
    +0.539543870454248793985299653154E-29,
    -0.238099288253145918675346190062E-29,
    +0.109905283010276157359726683750E-29,
    -0.486771374164496572732518677435E-30,
    +0.152587726411035756763200828211E-30 };
  static double erfcs[21] = {
    -0.49046121234691808039984544033376E-01,
    -0.14226120510371364237824741899631,
    +0.10035582187599795575754676712933E-01,
    -0.57687646997674847650827025509167E-03,
    +0.27419931252196061034422160791471E-04,
    -0.11043175507344507604135381295905E-05,
    +0.38488755420345036949961311498174E-07,
    -0.11808582533875466969631751801581E-08,
    +0.32334215826050909646402930953354E-10,
    -0.79910159470045487581607374708595E-12,
    +0.17990725113961455611967245486634E-13,
    -0.37186354878186926382316828209493E-15,
    +0.71035990037142529711689908394666E-17,
    -0.12612455119155225832495424853333E-18,
    +0.20916406941769294369170500266666E-20,
    -0.32539731029314072982364160000000E-22,
    +0.47668672097976748332373333333333E-24,
    -0.65980120782851343155199999999999E-26,
    +0.86550114699637626197333333333333E-28,
    -0.10788925177498064213333333333333E-29,
    +0.12811883993017002666666666666666E-31 };
  double eta;
  static int nterc2 = 0;
  static int nterf = 0;
  static int nterfc = 0;
  static double sqeps = 0.0;
  static double sqrtpi = 1.77245385090551602729816748334115;
  double value;
  static double xmax = 0.0;
  static double xsml = 0.0;
  double y;

  if ( nterf == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nterf = r8_inits ( erfcs, 21, eta );
    nterfc = r8_inits ( erfccs, 59, eta );
    nterc2 = r8_inits ( erc2cs, 49, eta );

    xsml = - sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) );
    xmax = sqrt (- log ( sqrtpi * r8_mach ( 1 ) ) );
    xmax = xmax - 0.5 * log ( xmax ) / xmax - 0.01;
    sqeps = sqrt ( 2.0 * r8_mach ( 3 ) );
  }

  if ( x <= xsml )
  {
    value = 2.0;
    return value;
  }

  if ( xmax < x )
  {
    cerr << "\n";
    cerr << "R8_ERFC - Warning!\n";
    cerr << "  X so big that ERFC underflows.\n";
    value = 0.0;
    return value;
  }

  y = r8_abs ( x );

  if ( y < sqeps )
  {
    value = 1.0 - 2.0 * x / sqrtpi;
    return value;
  }
  else if ( y <= 1.0 )
  {
    value = 1.0 - x * ( 1.0 
      + r8_csevl ( 2.0 * x * x - 1.0, erfcs, nterf ) );
    return value;
  }

  y = y * y;

  if ( y <= 4.0 )
  {
    value = exp ( - y ) / r8_abs ( x ) * ( 0.5 
      + r8_csevl ( ( 8.0 / y - 5.0 ) / 3.0, erc2cs, nterc2 ) );
  }
  else 
  {
    value = exp ( - y ) / r8_abs ( x ) * ( 0.5 
      + r8_csevl ( 8.0 / y - 1.0, erfccs, nterfc ) );
  }

  if ( x < 0.0 )
  {
    value = 2.0 - value;
  }

  return value;
}
//****************************************************************************80

double r8_exp ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXP evaluates the exponential of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_EXP, the exponential of X.
//
{
  static double aln216 = +0.83120654223414517758794896030274E-01;
  static double expcs[14] = {
    +0.866569493314985712733404647266231E-01,
    +0.938494869299839561896336579701203E-03,
    +0.677603970998168264074353014653601E-05,
    +0.366931200393805927801891250687610E-07,
    +0.158959053617461844641928517821508E-09,
    +0.573859878630206601252990815262106E-12,
    +0.177574448591421511802306980226000E-14,
    +0.480799166842372422675950244533333E-17,
    +0.115716376881828572809260000000000E-19,
    +0.250650610255497719932458666666666E-22,
    +0.493571708140495828480000000000000E-25,
    +0.890929572740634240000000000000000E-28,
    +0.148448062907997866666666666666666E-30,
    +0.229678916630186666666666666666666E-33 };
  double f;
  int n;
  int n16;
  int ndx;
  static int nterms = 0;
  static double twon16[17] = {
    +0.0,
    +0.44273782427413840321966478739929E-01,
    +0.90507732665257659207010655760707E-01,
    +0.13878863475669165370383028384151,
    +0.18920711500272106671749997056047,
    +0.24185781207348404859367746872659,
    +0.29683955465100966593375411779245,
    +0.35425554693689272829801474014070,
    +0.41421356237309504880168872420969,
    +0.47682614593949931138690748037404,
    +0.54221082540794082361229186209073,
    +0.61049033194925430817952066735740,
    +0.68179283050742908606225095246642,
    +0.75625216037329948311216061937531,
    +0.83400808640934246348708318958828,
    +0.91520656139714729387261127029583,
    +1.0 };
  double value;
  double xint;
  static double xmax;
  static double xmin;
  double y;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( expcs, 14, 0.1 * r8_mach ( 3 ) );
    xmin = log ( r8_mach ( 1 ) ) + 0.01;
    xmax = log ( r8_mach ( 2 ) ) - 0.001;
  }

  if ( x < xmin )
  {
    cerr << "\n";
    cerr << "R8_EXP - Warning!\n";
    cerr << "  X so small that exp(X) underflows.\n";
    value = 0.0;
  }
  else if ( x <= xmax )
  {
    xint = r8_aint ( x );
    y = x - xint;

    y = 23.0 * y + x * aln216;
    n = ( int ) ( y );
    f = y - ( double ) ( n );
    n = 23.0 * xint + ( double ) ( n );
    n16 = n / 16;
    if ( n < 0 )
    {
      n16 = n16 - 1;
    }
    ndx = n - 16 * n16 + 1;

    value = 1.0 + ( twon16[ndx-1] + f * ( 1.0 + twon16[ndx-1] ) 
      * r8_csevl ( f, expcs, nterms ) );

    value = r8_pak ( value, n16 );
  }
  else
  {
    cerr << "\n";
    cerr << "R8_EXP - Fatal error!\n";
    cerr << "  X so large that exp(X) overflows.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

double r8_exprel ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EXPREL evaluates the exponential relative error term of an R8 argument.
//
//  Discussion:
//
//    The relative error term is ( exp ( x ) - 1 ) / x.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_EXPREL, the exponential relative error term
//    at X.
//
{
  double absx;
  double alneps;
  int i;
  static int nterms = 0;
  double value;
  static double xbnd = 0.0;
  double xln;
  double xn;

  if ( nterms == 0 )
  {
    alneps = log ( r8_mach ( 3 ) );
    xn = 3.72 - 0.3 * alneps;
    xln = log ( ( xn + 1.0 ) / 1.36 );
    nterms = ( int ) ( xn - ( xn * xln + alneps ) / ( xln + 1.36 ) + 1.5 );
    xbnd = r8_mach ( 3 );
  }

  absx = r8_abs ( x );

  if ( absx < xbnd )
  {
    value = 1.0;
  }
  else if ( absx <= 0.5 )
  {
    value = 0.0;
    for ( i = 1; i <= nterms; i++ )
    {
      value = 1.0 + value * x / ( double ) ( nterms + 2 - i );
    }
  }
  else
  {
    value = ( exp ( x ) - 1.0 ) / x;
  }
  return value;
}
//****************************************************************************80

double r8_fac ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FAC evaluates the factorial of an I4 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, int N, the argument.
//
//    Output, double R8_FAC, the factorial of N.
//
{
  static double facn[31] = {
    +0.100000000000000000000000000000000E+01,
    +0.100000000000000000000000000000000E+01,
    +0.200000000000000000000000000000000E+01,
    +0.600000000000000000000000000000000E+01,
    +0.240000000000000000000000000000000E+02,
    +0.120000000000000000000000000000000E+03,
    +0.720000000000000000000000000000000E+03,
    +0.504000000000000000000000000000000E+04,
    +0.403200000000000000000000000000000E+05,
    +0.362880000000000000000000000000000E+06,
    +0.362880000000000000000000000000000E+07,
    +0.399168000000000000000000000000000E+08,
    +0.479001600000000000000000000000000E+09,
    +0.622702080000000000000000000000000E+10,
    +0.871782912000000000000000000000000E+11,
    +0.130767436800000000000000000000000E+13,
    +0.209227898880000000000000000000000E+14,
    +0.355687428096000000000000000000000E+15,
    +0.640237370572800000000000000000000E+16,
    +0.121645100408832000000000000000000E+18,
    +0.243290200817664000000000000000000E+19,
    +0.510909421717094400000000000000000E+20,
    +0.112400072777760768000000000000000E+22,
    +0.258520167388849766400000000000000E+23,
    +0.620448401733239439360000000000000E+24,
    +0.155112100433309859840000000000000E+26,
    +0.403291461126605635584000000000000E+27,
    +0.108888694504183521607680000000000E+29,
    +0.304888344611713860501504000000000E+30,
    +0.884176199373970195454361600000000E+31,
    +0.265252859812191058636308480000000E+33 };
  static int nmax = 0;
  static double sq2pil = 0.91893853320467274178032973640562;
  double value;
  double x;
  double xmax;
  double xmin;

  if ( nmax == 0 )
  {
    r8_gaml ( xmin, xmax );
    nmax = ( int ) ( xmax - 1.0 );
  }

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "R8_FAC - Fatal error!\n";
    cerr << "  Input argument is negative.\n";
    exit ( 1 );
  }
  else if ( n <= 30 )
  {
    value = facn[n];
  }
  else if ( n <= nmax )
  {
    x = ( double ) ( n + 1 );
    value = exp ( ( x - 0.5 ) * log ( x ) - x + sq2pil + r8_lgmc ( x ) );
  }
  else
  {
    cerr << "\n";
    cerr << "R8_FAC - Fatal error!\n";
    cerr << "  Factorial overflows.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

double r8_gami ( double a, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMI evaluates the incomplete gamma function for an R8 argument.
//
//  Discussion:
//
//    GAMI = Integral ( 0 <= T <= X ) exp ( - t ) * t^( a - 1 )  dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, the parameter.
//
//    Input, double X, the argument.
//
//    Output, double R8_GAMI, the value of the incomplete gamma function.
//
{
  double factor;
  double value;

  if ( a <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_GAMI - Fatal error!\n";
    cerr << "  A <= 0.\n";
    exit ( 1 );
  }

  if ( x < 0.0 )
  {
    cerr << "\n";
    cerr << "R8_GAMI - Fatal error!\n";
    cerr << "  X < 0.\n";
    exit ( 1 );
  }
  else if ( x == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    factor = exp ( r8_lngam ( a ) + a * log ( x ) );
    value = factor * r8_gamit ( a, x );
  }
  return value;
}
//****************************************************************************80

double r8_gamic ( double a, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMIC evaluates the complementary incomplete gamma function.
//
//  Discussion:
//
//    GAMIC = integral ( x <= t < oo ) exp(-t) * t^(a-1) dt
//
//    GAMIC is evaluated for arbitrary real values of A and non-negative
//    values X (even though GAMIC is defined for X < 0.0), except that
//    for X = 0 and A <= 0.0, GAMIC is undefined.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//    Walter Gautschi,
//    A Computational Procedure for Incomplete Gamma Functions,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 4, December 1979, pages 466-481.
//
//  Parameters:
//
//    Input, double A, the parameter.
//
//    Input, double X, the evaluation point.
//
//    Output, double R8_GAMIC, the value of the incomplete 
//    gamma function.
//
{
  double aeps;
  double ainta;
  double algap1;
  static double alneps = 0.0;
  double alngs;
  double alx;
  static double bot = 0.0;
  double e;
  static double eps = 0.0;
  double gstar;
  double h;
  int izero;
  int ma;
  double sga;
  double sgng;
  double sgngam;
  double sgngs;
  static double sqeps = 0.0;
  double t;
  double value;

  if ( eps == 0.0 )
  {
    eps = 0.5 * r8_mach ( 3 );
    sqeps = sqrt ( r8_mach ( 4 ) );
    alneps = - log ( r8_mach ( 3 ) );
    bot = log ( r8_mach ( 1 ) );
  }

  if ( x < 0.0 )
  {
    cerr << "\n";
    cerr << "R8_GAMIC - Fatal error!\n";
    cerr << "  X < 0.\n";
    exit ( 1 );
  }

  if ( x == 0.0 )
  {
    if ( a <= 0.0 )
    {
      cerr << "\n";
      cerr << "R8_GAMIC - Fatal error!\n";
      cerr << "  X = 0 and A <= 0.\n";
      exit ( 1 );
    }
    value = exp ( r8_lngam ( a + 1.0 ) - log ( a ) );

    return value;
  }

  alx = log ( x );
  if ( a < 0.0 )
  {
    sga = - 1.0;
  }
  else
  {
    sga = + 1.0;
  }

  ainta = r8_aint ( a + 0.5 * sga );
  aeps = a - ainta;

  izero = 0;

  if ( x < 1.0 )
  {
    if ( a <= 0.5 && r8_abs ( aeps ) <= 0.001 )
    {
      if ( - ainta <= 1.0 )
      {
        e = 2.0;
      }
      else
      {
        e = 2.0 * ( - ainta + 2.0 ) / ( ainta * ainta - 1.0 );
      }
      e = e - alx * r8_power ( x, - 0.001 );

      if ( e * r8_abs ( aeps ) <= eps )
      {
        value = r8_gmic ( a, x, alx );
        return value;
      }
    }

    r8_lgams ( a + 1.0, algap1, sgngam );
    gstar = r8_gmit ( a, x, algap1, sgngam, alx );

    if ( gstar == 0.0 )
    {
      izero = 1;
    }
    else
    {
      alngs = log ( r8_abs ( gstar ) );
      sgngs = r8_sign ( gstar );
    }
  }
  else
  {
    if ( a < x )
    {
      value = exp ( r8_lgic ( a, x, alx ) );
      return value;
    }

    sgngam = 1.0;
    algap1 = r8_lngam ( a + 1.0 );
    sgngs = 1.0;
    alngs = r8_lgit ( a, x, algap1 );
  }

  h = 1.0;

  if ( izero != 1 )
  {
    t = a * alx + alngs;

    if ( alneps < t )
    {
      sgng = - sgngs * sga * sgngam;
      t = t + algap1 - log ( r8_abs ( a ) );
      value = sgng * exp ( t );
      return value;
    }

    if ( - alneps < t )
    {
      h = 1.0 - sgngs * exp ( t );
    }
  }
  sgng = r8_sign ( h ) * sga * sgngam;
  t = log ( r8_abs ( h ) ) + algap1 - log ( r8_abs ( a ) );
  value = sgng * exp ( t );

  return value;
}
//****************************************************************************80

double r8_gamit ( double a, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMIT evaluates Tricomi's incomplete gamma function for an R8 argument.
//
//  Discussion:
//
//      GAMIT = x^(-a) / gamma(a) 
//        * Integral ( 0 <= t <= x ) exp(-t) * t^(a-1) dt
//
//    with analytic continuation for a <= 0.0.  Gamma(x) is the complete
//    gamma function of X.  GAMIT is evaluated for arbitrary real values of
//    A and for non-negative values of X (even though GAMIT is defined for
//    X < 0.0).
//
//    A slight deterioration of 2 or 3 digits accuracy will occur when
//    gamit is very large or very small in absolute value, because log-
//    arithmic variables are used.  Also, if the parameter A is very close
//    to a negative integer (but not a negative integer), there is a loss
//    of accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//    Walter Gautschi,
//    A Computational Procedure for Incomplete Gamma Functions,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 4, December 1979, pages 466-481.
//
//  Parameters:
//
//    Input, double A, the parameter.
//
//    Input, double X, the argument.
//
//    Output, double R8_GAMIT, the function value.
//
{
  double aeps;
  double ainta;
  double algap1;
  static double alneps = 0.0;
  double alng;
  double alx;
  static double bot = 0.0;
  double h;
  double sga;
  double sgngam;
  static double sqeps = 0.0;
  double t;
  double value;
 
  if ( alneps == 0.0 )
  {
    alneps = - log ( r8_mach ( 3 ) );
    sqeps = sqrt ( r8_mach ( 4 ) );
    bot = log ( r8_mach ( 1 ) );
  }

  if ( x < 0.0 )
  {
    cerr << "\n";
    cerr << "R8_GAMIT - Fatal error!\n";
    cerr << "  X is negative.\n";
    exit ( 1 );
  }
  else if ( x == 0.0 )
  {
    alx = 0.0;
  }
  else
  {
    alx = log ( x );
  }

  if ( a < 0.0 )
  {
    sga = - 1.0;
  }
  else
  {
    sga = + 1.0;
  }

  ainta = r8_aint ( a + 0.5 * sga );
  aeps = a - ainta;

  if ( x == 0.0 )
  {
    if ( 0.0 < ainta || aeps != 0.0 )
    {
      value = r8_gamr ( a + 1.0 );
    }
    else
    {
      value = 0.0;
    }
    return value;
  }

  if ( x <= 1.0 )
  {
    if ( - 0.5 <= a || aeps != 0.0 )
    {
      r8_lgams ( a + 1.0, algap1, sgngam );
    }
    value = r8_gmit ( a, x, algap1, sgngam, alx );
    return value;
  }

  if ( x <= a )
  {
    t = r8_lgit (a, x, r8_lngam ( a + 1.0 ) );
    value = exp ( t );
    return value;
  }

  alng = r8_lgic ( a, x, alx );
//
//  Evaluate in terms of log (r8_gamic (a, x))
//
  h = 1.0;

  if ( aeps != 0.0 || 0.0 < ainta )
  {
    r8_lgams ( a + 1.0, algap1, sgngam );
    t = log ( r8_abs ( a ) ) + alng - algap1;

    if ( alneps < t )
    {
      t = t - a * alx;
      value = - sga * sgngam * exp ( t );
      return value;
    }

    if ( - alneps < t )
    {
      h = 1.0 - sga * sgngam * exp ( t );
    }
  }
  t = - a * alx + log ( r8_abs ( h ) );

  if ( h < 0.0 )
  {
    value = - exp ( t );
  }
  else
  {
    value = + exp ( t );
  }
  return value;
}
//****************************************************************************80

void r8_gaml ( double &xmin, double &xmax )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAML evaluates bounds for an R8 argument of the gamma function.
//
//  Discussion:
//
//    This function calculates the minimum and maximum legal bounds 
//    for X in the evaluation of GAMMA ( X ).
//
//    XMIN and XMAX are not the only bounds, but they are the only 
//    non-trivial ones to calculate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Output, double &XMIN, &XMAX, the bounds.
//
{
  double alnbig;
  double alnsml;
  int i;
  int j;
  double xln;
  double xold;

  alnsml = log ( r8_mach ( 1 ) );
  xmin = - alnsml;

  for ( i = 1; i <= 10; i++ )
  {
    xold = xmin;
    xln = log ( xmin );
    xmin = xmin - xmin * ( ( xmin + 0.5 ) * xln - xmin 
      - 0.2258 + alnsml ) / ( xmin * xln + 0.5 );

    if ( r8_abs ( xmin - xold ) < 0.005 )
    {
      xmin = - xmin + 0.01;

      alnbig = log ( r8_mach ( 2 ) );
      xmax = alnbig;

      for ( j = 1; j <= 10; j++ )
      {
        xold = xmax;
        xln = log ( xmax );
        xmax = xmax - xmax * ( ( xmax - 0.5 ) * xln - xmax 
          + 0.9189 - alnbig ) / ( xmax * xln - 0.5 );

        if ( r8_abs ( xmax - xold ) < 0.005 )
        {
          xmax = xmax - 0.01;
          xmin = r8_max ( xmin, - xmax + 1.0 );
          return;
        }
      }
      cerr << "\n";
      cerr << "R8_GAML - Fatal error!\n";
      cerr << "  Unable to find XMAX.\n";
      exit ( 1 );
    }
  }
  cerr << "\n";
  cerr << "R8_GAML - Fatal error!\n";
  cerr << "  Unable to find XMIN.\n";
  exit ( 1 );
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates the gamma function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_GAMMA, the gamma function of X.
//
{
  static double dxrel = 0.0;
  static double gcs[42] = {
    +0.8571195590989331421920062399942E-02,
    +0.4415381324841006757191315771652E-02,
    +0.5685043681599363378632664588789E-01,
    -0.4219835396418560501012500186624E-02,
    +0.1326808181212460220584006796352E-02,
    -0.1893024529798880432523947023886E-03,
    +0.3606925327441245256578082217225E-04,
    -0.6056761904460864218485548290365E-05,
    +0.1055829546302283344731823509093E-05,
    -0.1811967365542384048291855891166E-06,
    +0.3117724964715322277790254593169E-07,
    -0.5354219639019687140874081024347E-08,
    +0.9193275519859588946887786825940E-09,
    -0.1577941280288339761767423273953E-09,
    +0.2707980622934954543266540433089E-10,
    -0.4646818653825730144081661058933E-11,
    +0.7973350192007419656460767175359E-12,
    -0.1368078209830916025799499172309E-12,
    +0.2347319486563800657233471771688E-13,
    -0.4027432614949066932766570534699E-14,
    +0.6910051747372100912138336975257E-15,
    -0.1185584500221992907052387126192E-15,
    +0.2034148542496373955201026051932E-16,
    -0.3490054341717405849274012949108E-17,
    +0.5987993856485305567135051066026E-18,
    -0.1027378057872228074490069778431E-18,
    +0.1762702816060529824942759660748E-19,
    -0.3024320653735306260958772112042E-20,
    +0.5188914660218397839717833550506E-21,
    -0.8902770842456576692449251601066E-22,
    +0.1527474068493342602274596891306E-22,
    -0.2620731256187362900257328332799E-23,
    +0.4496464047830538670331046570666E-24,
    -0.7714712731336877911703901525333E-25,
    +0.1323635453126044036486572714666E-25,
    -0.2270999412942928816702313813333E-26,
    +0.3896418998003991449320816639999E-27,
    -0.6685198115125953327792127999999E-28,
    +0.1146998663140024384347613866666E-28,
    -0.1967938586345134677295103999999E-29,
    +0.3376448816585338090334890666666E-30,
    -0.5793070335782135784625493333333E-31 };
  int i;
  int n;
  static int ngcs = 0;
  static double pi = 3.14159265358979323846264338327950;
  double sinpiy;
  static double sq2pil = 0.91893853320467274178032973640562;
  double value;
  static double xmax = 0.0;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ngcs == 0 )
  {
    ngcs = r8_inits ( gcs, 42, 0.1 * r8_mach ( 3 ) );
    r8_gaml ( xmin, xmax );
    xsml = exp ( r8_max ( log ( r8_mach ( 1 ) ),
      - log ( r8_mach ( 2 ) ) ) + 0.01 );
    dxrel = sqrt ( r8_mach ( 4 ) );
  }

  y = r8_abs ( x );

  if ( y <= 10.0 )
  {
    n = ( int ) ( x );
    if ( x < 0.0 )
    {
      n = n - 1;
    }
    y = x - ( double ) ( n );
    n = n - 1;
    value = 0.9375 + r8_csevl ( 2.0 * y - 1.0, gcs, ngcs );

    if ( n == 0 )
    {
      return value;
    }
    else if ( n < 0 )
    {
      n = - n;

      if ( x == 0.0 )
      {
        cerr << "\n";
        cerr << "R8_GAMMA - Fatal error!\n";
        cerr << "  X is 0.\n";
        exit ( 1 );
      }

      if ( x < 0.0 && x + ( double ) ( n - 2 ) == 0.0 )
      {
        cerr << "\n";
        cerr << "R8_GAMMA - Fatal error!\n";
        cerr << "  X is a negative int.\n";
        exit ( 1 );
      }

      if ( x < - 0.5 && r8_abs ( ( x - r8_aint ( x - 0.5 ) ) / x ) < dxrel )
      {
        cerr << "\n";
        cerr << "R8_GAMMA - Warning!\n";
        cerr << "  X too near a negative int,\n";
        cerr << "  answer is half precision.\n";
      }

      if ( y < xsml )
      {
        cerr << "\n";
        cerr << "R8_GAMMA - Fatal error!\n";
        cerr << "  X is so close to zero that Gamma overflows.\n";
        exit ( 1 );
      }

      for ( i = 1; i <= n; i++ )
      {
        value = value / ( x + ( double ) ( i - 1 ) );
      }

    }
    else if ( n == 0 )
    {
    }
    else
    {
      for ( i = 1; i <= n; i++ )
      {
        value = ( y + ( double ) ( i ) ) * value;
      }
    }
  }
  else
  {
    if ( xmax < x )
    {
      cerr << "\n";
      cerr << "R8_GAMMA - Fatal error!\n";
      cerr << "  X so big that Gamma overflows.\n";
      exit ( 1 );
    }
//
//  Underflow.
//
    if ( x < xmin )
    {
      value = 0.0;
      return value;
    }

    value = exp ( ( y - 0.5 ) * log ( y ) - y + sq2pil + r8_lgmc ( y ) );

    if ( 0.0 < x )
    {
      return value;
    }

    if ( r8_abs ( ( x - r8_aint ( x - 0.5 ) ) / x ) < dxrel )
    {
      cerr << "\n";
      cerr << "R8_GAMMA - Warning!\n";
      cerr << "  X too near a negative int,\n";
      cerr << "  answer is half precision.\n";
    }

    sinpiy = sin ( pi * y );

    if ( sinpiy == 0.0 )
    {
      cerr << "\n";
      cerr << "R8_GAMMA - Fatal error!\n";
      cerr << "  X is a negative int.\n";
      exit ( 1 );
    }
    value = - pi / ( y * sinpiy * value );
  }
  return value;
}
//****************************************************************************80

double r8_gamr ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMR evaluates the reciprocal gamma function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_GAMR, the value of the reciprocal gamma
//    function at X.
//
{
  double alngx;
  double value;
  double sgngx;

  if ( x <= 0.0 && r8_aint ( x ) == x )
  {
    value = 0.0;
  } 
  else if ( r8_abs ( x ) <= 10.0 )
  {
    value = 1.0 / r8_gamma ( x );
  }
  else
  {
    r8_lgams ( x, alngx, sgngx );
    value = sgngx * exp ( - alngx );
  }
  return value;
}
//****************************************************************************80

double r8_gmic ( double a, double x, double alx )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GMIC: complementary incomplete gamma, small X, A near negative int.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, the parameter.
//
//    Input, double X, the argument.
//
//    Input, double ALX, the logarithm of X.
//
//    Output, double R8_GMIC, the complementary incomplete 
//    gamma function.
//
{
  double alng;
  static double bot = 0.0;
  bool converged;
  static double eps = 0.0;
  static double euler = 0.57721566490153286060651209008240;
  double fk;
  double fkp1;
  double fm;
  int k;
  int m;
  int ma;
  int mm1;
  double s;
  double sgng;
  double t;
  double te;
  double value;

  if ( eps == 0.0 )
  {
    eps = 0.5 * r8_mach ( 3 );
    bot = log ( r8_mach ( 1 ) );
  }

  if ( 0.0 < a )
  {
    cerr << "\n";
    cerr << "R8_GMIC - Fatal error!\n";
    cerr << "  A must be near a negative int.\n";
    exit ( 1 );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_GMIC - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  m = - ( int ) ( a - 0.5 );
  fm = ( double ) ( m );

  te = 1.0;
  t = 1.0;
  s = t;
  converged = false;

  for ( k = 1; k <= 200; k++ )
  {
    fkp1 = ( double ) ( k + 1 );
    te = - x * te / ( fm + fkp1 );
    t = te / fkp1;
    s = s + t;
    if ( r8_abs ( t ) < eps * s )
    {
      converged = true;
      break;
    }
  }

  if ( !converged )
  {
    cerr << "\n";
    cerr << "R8_GMIC - Fatal error!\n";
    cerr << "  No convergence after 200 iterations.\n";
    exit ( 1 );
  }

  value = - alx - euler + x * s / ( fm + 1.0 );

  if ( m == 0 )
  {
    return value;
  }
  else if ( m == 1 )
  {
    value = - value - 1.0 + 1.0 / x;
    return value;
  }

  te = fm;
  t = 1.0;
  s = t;
  mm1 = m - 1;

  for ( k = 1; k <= mm1; k++ )
  {
    fk = ( double ) ( k );
    te = - x * te / fk;
    t = te / ( fm - fk );
    s = s + t;
    if ( r8_abs ( t ) < eps * r8_abs ( s ) )
    {
      break;
    }
  }

  for ( k = 1; k <= m; k++ )
  {
    value = value + 1.0 / ( double ) ( k );
  }

  if ( ( m % 2 ) == 1 )
  {
    sgng = - 1.0;
  }
  else
  {
    sgng = + 1.0;
  }

  alng = log ( value ) - r8_lngam ( fm + 1.0 );

  if ( bot < alng )
  {
    value = sgng * exp ( alng );
  }
  else
  {
    value = 0.0;
  }

  if ( s != 0.0 )
  {
    value = value
      + r8_sign ( s ) * exp ( - fm * alx + log ( r8_abs ( s ) / fm ) );
  }

  return value;
}
//****************************************************************************80

double r8_gmit ( double a, double x, double algap1, double sgngam, double alx )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GMIT: Tricomi's incomplete gamma function for small X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, the parameter.
//
//    Input, double X, the argument.
//
//    Input, double ALGAP1, the logarithm of Gamma ( A + 1 ).
//
//    Input, double SGNGAM, the sign of Gamma ( A + 1 ).
//
//    Input, double ALX, the logarithm of X.
//
//    Output, double R8_GMIT, the Tricomi incomplete gamma function.
//
{
  double ae;
  double aeps;
  double alg2;
  double algs;
  static double bot = 0.0;
  bool converged;
  static double eps = 0.0;
  double fk;
  int k;
  int m;
  int ma;
  double s;
  double sgng2;
  double t;
  double te;
  double value;

  if ( eps == 0.0 )
  {
    eps = 0.5 * r8_mach ( 3 );
    bot = log ( r8_mach ( 1 ) );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_GMIT - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  if ( a < 0.0 )
  {
    ma = ( int ) ( a - 0.5 );
  }
  else
  {
    ma = ( int ) ( a + 0.5 );
  }

  aeps = a - ( double ) ( ma );

  if ( a < - 0.5 )
  {
    ae = aeps;
  }
  else
  {
    ae = a;
  }

  t = 1.0;
  te = ae;
  s = t;
  converged = false;

  for ( k = 1; k <= 200; k++ )
  {
    fk = ( double ) ( k );
    te = - x * te / fk;
    t = te / ( ae + fk );
    s = s + t;
    if ( r8_abs ( t ) < eps * r8_abs ( s ) )
    {
      converged = true;
      break;
    }
  }

  if ( !converged )
  {
    cerr << "\n";
    cerr << "R8_GMIT - Fatal error!\n";
    cerr << "  No convergence in 200 iterations.\n";
    exit ( 1 );
  }

  if ( - 0.5 <= a )
  {
    algs = - algap1 + log ( s );
    value = exp ( algs );
    return value;
  }

  algs = - r8_lngam ( 1.0 + aeps ) + log ( s );
  s = 1.0;
  m = - ma - 1;
  t = 1.0;

  for ( k = 1; k <= m; k++ )
  {
    t = x * t / ( aeps - ( double ) ( m + 1 - k ) );
    s = s + t;
    if ( r8_abs ( t ) < eps * r8_abs ( s ) )
    {
      break;
    }
  }

  value = 0.0;
  algs = - ( double ) ( ma ) * log ( x ) + algs;

  if ( s == 0.0 || aeps == 0.0 )
  {
    value = exp ( algs );
    return value;
  }

  sgng2 = sgngam * r8_sign ( s );
  alg2 = - x - algap1 + log ( r8_abs ( s ) );

  if ( bot < alg2 )
  {
    value = sgng2 * exp ( alg2 );
  }

  if ( bot < algs )
  {
    value = value + exp ( algs );
  }

  return value;
}
//****************************************************************************80

int r8_inits ( double dos[], int nos, double eta )

//****************************************************************************80
//
//  Purpose:
//
//    R8_INITS initializes a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double DOS[NOS], the Chebyshev coefficients.
//
//    Input, int NOS, the number of coefficients.
//
//    Input, double ETA, the desired accuracy.
//
//    Output, int R8_INITS, the number of terms of the series needed
//    to ensure the requested accuracy.
//
{
  double err;
  int i;
  int value;

  if ( nos < 1 )
  {
    cerr << "\n";
    cerr << "R8_INITS - Fatal error!\n";
    cerr << "  Number of coefficients < 1.\n";
    exit ( 1 );
  }

  err = 0.0;

  for ( i = nos - 1; 0 <= i; i-- )
  {
    err = err + r8_abs ( dos[i] );
    if ( eta < err )
    {
      value = i + 1;
      return value;
    }
  }

  value = i;
  cerr << "\n";
  cerr << "R8_INITS - Warning!\n";
  cerr << "  ETA may be too small.\n";

  return value;
}
//****************************************************************************80

double r8_int ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_INT returns the integer part of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_INT, the integer part of X.
//
{
  int i;
  int ibase;
  int ipart;
  static int npart = 0;
  double part;
  static double scale = 0.0;
  double value;
  static double xbig = 0.0;
  static double xmax = 0.0;
  double xscl;

  if ( npart == 0 )
  {
    ibase = i4_mach ( 10 );
    xmax = 1.0 / r8_mach ( 4 );
    xbig = r8_min ( ( double ) ( i4_mach ( 9 ) ), 1.0 / r8_mach ( 4 ) );
    scale = ( double ) i4_pow ( ibase,
      ( int ) ( log ( xbig ) / log ( ( double ) ( ibase ) ) - 0.5 ) );
    npart = log ( xmax ) / log ( scale ) + 1.0;
  }
//
//  X may be too small.
//
  if ( x < - xmax )
  {
    value = x;
  }
  else if ( x < - xbig )
  {
    xscl = - x;

    for ( i = 1; i <= npart; i++ )
    {
      xscl = xscl / scale;
    }

    value = 0.0;
    for ( i = 1; i <= npart; i++ )
    {
      xscl = xscl * scale;
      ipart = ( int ) ( xscl );
      part = ( double ) ( ipart );
      xscl = xscl - part;
      value = value * scale + part;
    }
    value = - value;
  }
  else if ( x <= xbig )
  {
    value = ( int ) ( x );
  }
  else if ( x <= xmax )
  {
    xscl = x;

    for ( i = 1; i <= npart; i++ )
    {
      xscl = xscl / scale;
    }

    value = 0.0;
    for ( i = 1; i <= npart; i++ )
    {
      xscl = xscl * scale;
      ipart = ( int ) ( xscl );
      part = ( double ) ( ipart );
      xscl = xscl - part;
      value = value * scale + part;
    }
  }
//
//  X may be too large.
//
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

void r8_knus ( double xnu, double x, double &bknu, double &bknu1, int &iswtch )

//****************************************************************************80
//
//  Purpose:
//
//    R8_KNUS computes a sequence of K Bessel functions.
//
//  Discussion:
//
//    This routine computes Bessel functions 
//      exp(x) * k-sub-xnu (x)  
//    and
//      exp(x) * k-sub-xnu+1 (x) 
//    for 0.0 <= xnu < 1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double XNU, the order parameter.
//
//    Input, double X, the argument.
//
//    Output, double &BKNU, &BKNU1, the two K Bessel functions.
//
//    Output, int &ISWTCH, ?
//
{
  double a[32];
  double a0;
  static double aln2 = 0.69314718055994530941723212145818;
  static double alnbig = 0;
  static double alneps = 0;
  static double alnsml = 0;
  double alnz;
  double alpha[32];
  double an;
  double b0;
  double beta[32];
  double bknu0;
  double bknud;
  double bn;
  double c0;
  static double c0kcs[29] = {
    +0.60183057242626108387577445180329E-01,
    -0.15364871433017286092959755943124,
    -0.11751176008210492040068229226213E-01,
    -0.85248788891979509827048401550987E-03,
    -0.61329838767496791874098176922111E-04,
    -0.44052281245510444562679889548505E-05,
    -0.31631246728384488192915445892199E-06,
    -0.22710719382899588330673771793396E-07,
    -0.16305644608077609552274620515360E-08,
    -0.11706939299414776568756044043130E-09,
    -0.84052063786464437174546593413792E-11,
    -0.60346670118979991487096050737198E-12,
    -0.43326960335681371952045997366903E-13,
    -0.31107358030203546214634697772237E-14,
    -0.22334078226736982254486133409840E-15,
    -0.16035146716864226300635791528610E-16,
    -0.11512717363666556196035697705305E-17,
    -0.82657591746836959105169479089258E-19,
    -0.59345480806383948172333436695984E-20,
    -0.42608138196467143926499613023976E-21,
    -0.30591266864812876299263698370542E-22,
    -0.21963541426734575224975501815516E-23,
    -0.15769113261495836071105750684760E-24,
    -0.11321713935950320948757731048056E-25,
    -0.81286248834598404082792349714433E-27,
    -0.58360900893453226552829349315949E-28,
    -0.41901241623610922519452337780905E-29,
    -0.30083737960206435069530504212862E-30,
    -0.21599152067808647728342168089832E-31 };
  double eta;
  static double euler = 0.57721566490153286060651209008240;
  double expx;
  int i;
  int ii;
  int inu;
  int n;
  static int ntc0k = 0;
  int nterms;
  static int ntznu1 = 0;
  double p1;
  double p2;
  double p3;
  double qq;
  double result;
  static double sqpi2 = +1.2533141373155002512078826424055;
  double sqrtx;
  double v;
  double vlnz;
  double x2n;
  double x2tov;
  double xi;
  double xmu;
  static double xnusml = 0.0;
  static double xsml = 0.0;
  double z;
  static double znu1cs[20] = {
    +0.203306756994191729674444001216911,
    +0.140077933413219771062943670790563,
    +0.791679696100161352840972241972320E-02,
    +0.339801182532104045352930092205750E-03,
    +0.117419756889893366664507228352690E-04,
    +0.339357570612261680333825865475121E-06,
    +0.842594176976219910194629891264803E-08,
    +0.183336677024850089184748150900090E-09,
    +0.354969844704416310863007064469557E-11,
    +0.619032496469887332205244342078407E-13,
    +0.981964535680439424960346115456527E-15,
    +0.142851314396490474211473563005985E-16,
    +0.191894921887825298966162467488436E-18,
    +0.239430979739498914162313140597128E-20,
    +0.278890246815347354835870465474995E-22,
    +0.304606650633033442582845214092865E-24,
    +0.313173237042191815771564260932089E-26,
    +0.304133098987854951645174908005034E-28,
    +0.279840384636833084343185097659733E-30,
    +0.244637186274497596485238794922666E-32 };
  double ztov;

  if ( ntc0k == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    ntc0k = r8_inits ( c0kcs, 29, eta );
    ntznu1 = r8_inits ( znu1cs, 20, eta );
    xnusml = sqrt ( r8_mach ( 3 ) / 8.0 );
    xsml = 0.1 * r8_mach ( 3 );
    alnsml = log ( r8_mach ( 1 ) );
    alnbig = log ( r8_mach ( 2 ) );
    alneps = log ( 0.1 * r8_mach ( 3 ) );
  }

  if ( xnu < 0.0 || 1.0 <= xnu )
  {
    cerr << "\n";
    cerr << "R8_KNUS - Fatal error!\n";
    cerr << "  XNU < 0 or 1 <= XNU.\n";
    exit ( 1 );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_KNUS - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  iswtch = 0;
//
//  X is small.  Compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
//  then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-0.5,+0.5)
//  then to (0., .5), because k of negative order (-nu) = k of positive
//  order (+nu).
//
  if ( x <= 2.0 )
  {
    if ( xnu <= 0.5 )
    {
      v = xnu;
    }
    else
    {
      v = 1.0 - xnu;
    }
//
//  carefully find (x/2)^xnu and z^xnu where z = x*x/4.
//
    alnz = 2.0 * ( log ( x ) - aln2 );

    if ( x <= xnu )
    {
      if ( alnbig < - 0.5 * xnu * alnz - aln2 - log ( xnu ) )
      {
        cerr << "\n";
        cerr << "R8_KNUS - Fatal error!\n";
        cerr << "  Small X causing overflow.\n";
        exit ( 1 );
      }
    }

    vlnz = v * alnz;
    x2tov = exp ( 0.5 * vlnz );

    if ( vlnz <= alnsml )
    {
      ztov = 0.0;
    }
    else
    {
      ztov = x2tov * x2tov;
    }

    a0 = 0.5 * r8_gamma ( 1.0 + v );
    b0 = 0.5 * r8_gamma ( 1.0 - v );
    c0 = - euler;

    if ( 0.5 <= ztov && xnusml < v )
    {
      c0 = - 0.75 + r8_csevl ( ( 8.0 * v ) * v - 1.0, c0kcs, ntc0k );
    }

    if ( ztov <= 0.5 )
    {
      alpha[0] = ( a0 - ztov * b0 ) / v;
    }
    else
    {
      alpha[0] = c0 - alnz * ( 0.75 +
        r8_csevl ( vlnz / 0.35 + 1.0, znu1cs, ntznu1 ) ) * b0;
    }

    beta[0] = - 0.5 * ( a0 + ztov * b0 );

    if ( x <= xsml )
    {
      z = 0.0;
    }
    else
    {
      z = 0.25 * x * x;
    }

    nterms = i4_max ( 2, ( int ) ( 11.0 
      + ( 8.0 * alnz - 25.19 - alneps ) / ( 4.28 - alnz ) ) );

    for ( i = 2; i <= nterms; i++ )
    {
      xi = ( double ) ( i - 1 );
      a0 = a0 / ( xi * ( xi - v ) );
      b0 = b0 / ( xi * ( xi + v ) );
      alpha[i-1] = ( alpha[i-2] + 2.0 * xi * a0 ) 
        / ( xi * ( xi + v ) );
      beta[i-1] = ( xi - 0.5 * v ) * alpha[i-1] - ztov * b0;
    }

    bknu = alpha[nterms-1];
    bknud = beta[nterms-1];
    for ( ii = 2; ii <= nterms; ii++ )
    {
      i = nterms + 1 - ii;
      bknu = alpha[i-1] + bknu * z;
      bknud = beta[i-1] + bknud * z;
    }

    expx = exp ( x );
    bknu = expx * bknu / x2tov;

    if ( alnbig < - 0.5 * ( xnu + 1.0 ) * alnz - 2.0 * aln2 )
    {
      iswtch = 1;
      return;
    }

    bknud = expx * bknud * 2.0 / ( x2tov * x );

    if ( xnu <= 0.5 )
    {
      bknu1 = v * bknu / x - bknud;
      return;
    }
    bknu0 = bknu;
    bknu = - v * bknu / x - bknud;
    bknu1 = 2.0 * xnu * bknu / x + bknu0;
  }
//
//  x is large.  find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke-s
//  rational expansion.
//
  else
  {
    sqrtx = sqrt ( x );

    if ( 1.0 / xsml < x )
    {
      bknu = sqpi2 / sqrtx;
      bknu1 = bknu;
      return;
    }

    an = - 0.60 - 1.02 / x;
    bn = - 0.27 - 0.53 / x;
    nterms = i4_min ( 32, i4_max ( 3, ( int ) ( an + bn * alneps ) ) );

    for ( inu = 1; inu <= 2; inu++ )
    {
      if ( inu == 1 )
      {
        if ( xnu <= xnusml )
        {
          xmu = 0.0;
        }
        else
        {
          xmu = ( 4.0 * xnu ) * xnu;
        }
      }
      else
      {
        xmu = 4.0 * ( r8_abs ( xnu ) + 1.0 ) * ( r8_abs ( xnu ) + 1.0 );
      }

      a[0] = 1.0 - xmu;
      a[1] = 9.0 - xmu;
      a[2] = 25.0 - xmu;

      if ( a[1] == 0.0 )
      {
        result = sqpi2 * ( 16.0 * x + xmu + 7.0 ) / ( 16.0 * x * sqrtx );
      }
      else
      {
        alpha[0] = 1.0;
        alpha[1] = ( 16.0 * x + a[1] ) / a[1];
        alpha[2] = ( ( 768.0 * x + 48.0 * a[2] ) * x 
          + a[1] * a[2] ) / ( a[1] * a[2] );

        beta[0] = 1.0;
        beta[1] = ( 16.0 * x + ( xmu + 7.0 ) ) / a[1];
        beta[2] = ( ( 768.0 * x + 48.0 * ( xmu + 23.0 ) ) * x +
          ( ( xmu + 62.0 ) * xmu + 129.0 ) ) / ( a[1] * a[2] );

        for ( i = 4; i <= nterms; i++ )
        {
          n = i - 1;
          x2n = ( double ) ( 2 * n - 1 );

          a[i-1] = ( x2n + 2.0 ) * ( x2n + 2.0 ) - xmu;
          qq = 16.0 * x2n / a[i-1];
          p1 = - x2n * ( ( double ) ( 12 * n * n - 20 * n ) - a[0] ) 
            / ( ( x2n - 2.0 ) * a[i-1] ) - qq * x;
          p2 = ( ( double ) ( 12 * n * n - 28 * n + 8 ) - a[0] ) 
            / a[i-1] - qq * x;
          p3 = - x2n * a[i-4] / ( ( x2n - 2.0 ) * a[i-1] );

          alpha[i-1] = - p1 * alpha[i-2]
                       - p2 * alpha[i-3] 
                       - p3 * alpha[i-4];

          beta[i-1] =  - p1 * beta[i-2]
                       - p2 * beta[i-3] 
                       - p3 * beta[i-4];

        }
        result = sqpi2 * beta[nterms-1] / ( sqrtx * alpha[nterms-1] );
      }

      if ( inu == 1 )
      {
        bknu = result;
      }
      else
      {
        bknu1 = result;
      }
    }
  }
  return;
}
//****************************************************************************80

double r8_lbeta ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LBETA evaluates the logarithm of the beta function of R8 arguments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, B, the arguments.
//
//    Output, double R8_LBETA, the logarithm of the beta function of A
//    and B.
//
{
  double corr;
  double p;
  double q;
  static double sq2pil = 0.91893853320467274178032973640562;
  double value;

  p = r8_min ( a, b );
  q = r8_max ( a, b );

  if ( p <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_LBETA - Fatal error!\n";
    cerr << "  Both arguments must be greater than 0.\n";
    exit ( 1 );
  }
  else if ( p < 10.0 && q <= 10.0 )
  {
    value = log ( r8_gamma ( p ) * ( r8_gamma ( q ) / r8_gamma ( p + q ) ) );
  }
  else if ( p < 10.0 )
  {
    corr = r8_lgmc ( q ) - r8_lgmc ( p + q );

    value = r8_lngam ( p ) + corr + p - p * log ( p + q ) +
      ( q - 0.5 ) * r8_lnrel ( - p / ( p + q ) );
  }
  else
  {
    corr = r8_lgmc ( p ) + r8_lgmc ( q ) - r8_lgmc ( p + q );

    value = - 0.5 * log ( q ) + sq2pil + corr 
      + ( p - 0.5 ) * log ( p / ( p + q ) ) 
      + q * r8_lnrel ( - p / ( p + q ) );
  }
  return value;
}

//****************************************************************************80

void r8_lgams ( double x, double &algam, double &sgngam )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LGAMS evaluates the log of |gamma(x)| and sign, for an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double &ALGAM, the logarithm of the absolute value of
//    gamma ( X ).
//
//    Output, double &SGNGAM, the sign (+1 or -1) of gamma ( X ).
//
{
  int k;

  algam = r8_lngam ( x );
  sgngam = 1.0;

  if ( x <= 0.0 )
  {
    k = ( int ) ( fmod ( - r8_aint ( x ), 2.0 ) + 0.1 );

    if ( k == 0 )
    {
      sgngam = - 1.0;
    }
  }
  return;
}
//****************************************************************************80

double r8_lgic ( double a, double x, double alx )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LGIC evaluates the log complementary incomplete gamma function for large X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, the parameter.
//
//    Input, double X, the argument.
//
//    Input, double ALX, the logarithm of X.
//
//    Output, double R8_LGIC, the log complementary incomplete 
//    gamma function.
//
{
  static double eps = 0.0;
  double fk;
  int k;
  double p;
  double r;
  double s;
  double t;
  double value;
  double xma;
  double xpa;

  if ( eps == 0.0 )
  {
    eps = 0.5 * r8_mach ( 3 );
  }

  xpa = x + 1.0 - a;
  xma = x - 1.0 - a;

  r = 0.0;
  p = 1.0;
  s = p;
  for ( k = 1; k <= 300; k++ )
  {
    fk = ( double ) ( k );
    t = fk * ( a - fk ) * ( 1.0 + r );
    r = - t / ( ( xma + 2.0 * fk ) * ( xpa + 2.0 * fk ) + t );
    p = r * p;
    s = s + p;
    if ( r8_abs ( p ) < eps * s )
    {
      value = a * alx - x + log ( s / xpa );
      return value;
    }
  }

  cerr << "\n";
  cerr << "R8_LGIC - Fatal error!\n";
  cerr << "  No convergence in 300 iterations.\n";

  exit ( 1 );
}
//****************************************************************************80

double r8_lgit ( double a, double x, double algap1 )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LGIT evaluates the log of Tricomi's incomplete gamma function.
//
//  Discussion:
//
//    Perron's continued fraction is used for large X and X <= A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, the parameter.
//
//    Input, double X, the argument.
//
//    Input, double ALGAP1, the logarithm of A+1.
//
//    Output, double R8_LGIT, the log of Tricomi's incomplete
//    gamma function.
//
{
  double a1x;
  double ax;
  static double eps = 0.0;
  double fk;
  double hstar;
  int k;
  double p;
  double r;
  double s;
  static double sqeps = 0.0;
  double t;
  double value;

  if ( eps == 0.0 )
  {
    eps = 0.5 * r8_mach ( 3 );
    sqeps = sqrt ( r8_mach ( 4 ) );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_LGIT - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  if ( a < x )
  {
    cerr << "\n";
    cerr << "R8_LGIT - Fatal error!\n";
    cerr << "  A < X.\n";
    exit ( 1 );
  }

  ax = a + x;
  a1x = ax + 1.0;
  r = 0.0;
  p = 1.0;
  s = p;

  for ( k = 1; k <= 200; k++ )
  {
    fk = ( double ) ( k );
    t = ( a + fk ) * x * ( 1.0 + r );
    r = t / ( ( ax + fk ) * ( a1x + fk ) - t );
    p = r * p;
    s = s + p;
    if ( r8_abs ( p ) < eps * s )
    {
      hstar = 1.0 - x * s / a1x;
      value = - x - algap1 - log ( hstar );
      return value;
    }
  }

  cerr << "\n";
  cerr << "R8_LGIT - Fatal error!\n";
  cerr << "  No convergence after 200 iterations.\n";
  exit ( 1 );
}
//****************************************************************************80

double r8_lgmc ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LGMC evaluates the log gamma correction factor for an R8 argument.
//
//  Discussion:
//
//    For 10 <= X, compute the log gamma correction factor so that
//
//      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
//                          + ( x - 0.5 ) * log ( x ) - x 
//                          + r8_lgmc ( x )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_LGMC, the correction factor.
//
{
  static double algmcs[15] = {
    +0.1666389480451863247205729650822,
    -0.1384948176067563840732986059135E-04,
    +0.9810825646924729426157171547487E-08,
    -0.1809129475572494194263306266719E-10,
    +0.6221098041892605227126015543416E-13,
    -0.3399615005417721944303330599666E-15,
    +0.2683181998482698748957538846666E-17,
    -0.2868042435334643284144622399999E-19,
    +0.3962837061046434803679306666666E-21,
    -0.6831888753985766870111999999999E-23,
    +0.1429227355942498147573333333333E-24,
    -0.3547598158101070547199999999999E-26,
    +0.1025680058010470912000000000000E-27,
    -0.3401102254316748799999999999999E-29,
    +0.1276642195630062933333333333333E-30 };
  static int nalgm = 0;
  double value;
  static double xbig = 0.0;
  static double xmax = 0.0;

  if ( nalgm == 0 )
  {
    nalgm = r8_inits ( algmcs, 15, r8_mach ( 3 ) );
    xbig = 1.0 / sqrt ( r8_mach ( 3 ) );
    xmax = exp ( r8_min ( log ( r8_mach ( 2 ) / 12.0 ), 
      - log ( 12.0 * r8_mach ( 1 ) ) ) );
  }

  if ( x < 10.0 )
  {
    cerr << "\n";
    cerr << "R8_LGMC - Fatal error!\n";
    cerr << "  X must be at least 10.\n";
    exit ( 1 );
  }
  else if ( x < xbig )
  {
    value = r8_csevl ( 2.0 * ( 10.0 / x ) 
      * ( 10.0 / x ) - 1.0, algmcs, nalgm ) / x;
  }
  else if ( x < xmax )
  {
    value = 1.0 / ( 12.0 * x );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double r8_li ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LI evaluates the logarithmic integral for an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_LI, the logarithmic integral evaluated at X.
//
{
  static double sqeps = 0.0;
  double value;

  if ( sqeps == 0.0 )
  {
    sqeps = sqrt ( r8_mach ( 3 ) );
  }

  if ( x < 0.0 )
  {
    cerr << "\n";
    cerr << "R8_LI - Fatal error!\n";
    cerr << "  Function undefined for X <= 0.\n";
    exit ( 1 );
  }

  if ( x == 0.0 )
  {
    value = 0.0;
    return value;
  }

  if ( x == 1.0 )
  {
    cerr << "\n";
    cerr << "R8_LI - Fatal error!\n";
    cerr << "  Function undefined for X = 1.\n";
    exit ( 1 );
  }

  if ( r8_abs ( 1.0 - x ) < sqeps )
  {
    cerr << "\n";
    cerr << "R8_LI - Warning!\n";
    cerr << "  Answer less than half precision.\n";
    cerr << "  X is too close to 1.\n";
  }

  value = r8_ei ( log ( x ) );

  return value;
}
//****************************************************************************80

double r8_lngam ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LNGAM: log of the absolute value of gamma of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_LNGAM, the logarithm of the absolute value of
//    the gamma function of X.
//
{
  static double dxrel = 0.0;
  static double pi = 3.14159265358979323846264338327950;
  double sinpiy;
  static double sq2pil = 0.91893853320467274178032973640562;
  static double sqpi2l = +0.225791352644727432363097614947441;
  double value;
  static double xmax = 0.0;
  double y;

  if ( xmax == 0.0 )
  {
    xmax = r8_mach ( 2 ) / log ( r8_mach ( 2 ) );
    dxrel = sqrt ( r8_mach ( 4 ) );
  }

  y = r8_abs ( x );

  if ( y <= 10.0 )
  {
    value = log ( r8_abs ( r8_gamma ( x ) ) );
    return value;
  }

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R8_LNGAM - Fatal error!\n";
    cerr << "  Result overflows, |X| too big.\n";
    exit ( 1 );
  }

  if ( 0.0 < x )
  {
    value = sq2pil + ( x - 0.5 ) * log ( x ) - x + r8_lgmc ( y );
    return value;
  }

  sinpiy = r8_abs ( sin ( pi * y ) );

  if ( sinpiy == 0.0 )
  {
    cerr << "\n";
    cerr << "R8_LNGAM - Fatal error!\n";
    cerr << "  X is a negative int.\n";
    exit ( 1 );
  }

  value = sqpi2l + ( x - 0.5 ) * log ( y ) - x - log ( sinpiy ) - r8_lgmc ( y );

  if ( r8_abs ( ( x - r8_aint ( x - 0.5 ) ) * value / x ) < dxrel )
  {
    cerr << "\n";
    cerr << "R8_LNGAM - Warning!\n";
    cerr << "  Result is half precision because\n";
    cerr << "  X is too near a negative int.\n";
  }

  return value;
}
//****************************************************************************80

double r8_lnrel ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LNREL evaluates log ( 1 + X ) for an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_LNREL, the value of LOG ( 1 + X ).
//
{
  static double alnrcs[43] = {
    +0.10378693562743769800686267719098E+01,
    -0.13364301504908918098766041553133,
    +0.19408249135520563357926199374750E-01,
    -0.30107551127535777690376537776592E-02,
    +0.48694614797154850090456366509137E-03,
    -0.81054881893175356066809943008622E-04,
    +0.13778847799559524782938251496059E-04,
    -0.23802210894358970251369992914935E-05,
    +0.41640416213865183476391859901989E-06,
    -0.73595828378075994984266837031998E-07,
    +0.13117611876241674949152294345011E-07,
    -0.23546709317742425136696092330175E-08,
    +0.42522773276034997775638052962567E-09,
    -0.77190894134840796826108107493300E-10,
    +0.14075746481359069909215356472191E-10,
    -0.25769072058024680627537078627584E-11,
    +0.47342406666294421849154395005938E-12,
    -0.87249012674742641745301263292675E-13,
    +0.16124614902740551465739833119115E-13,
    -0.29875652015665773006710792416815E-14,
    +0.55480701209082887983041321697279E-15,
    -0.10324619158271569595141333961932E-15,
    +0.19250239203049851177878503244868E-16,
    -0.35955073465265150011189707844266E-17,
    +0.67264542537876857892194574226773E-18,
    -0.12602624168735219252082425637546E-18,
    +0.23644884408606210044916158955519E-19,
    -0.44419377050807936898878389179733E-20,
    +0.83546594464034259016241293994666E-21,
    -0.15731559416479562574899253521066E-21,
    +0.29653128740247422686154369706666E-22,
    -0.55949583481815947292156013226666E-23,
    +0.10566354268835681048187284138666E-23,
    -0.19972483680670204548314999466666E-24,
    +0.37782977818839361421049855999999E-25,
    -0.71531586889081740345038165333333E-26,
    +0.13552488463674213646502024533333E-26,
    -0.25694673048487567430079829333333E-27,
    +0.48747756066216949076459519999999E-28,
    -0.92542112530849715321132373333333E-29,
    +0.17578597841760239233269760000000E-29,
    -0.33410026677731010351377066666666E-30,
    +0.63533936180236187354180266666666E-31 };
  static int nlnrel = 0;
  double value;
  static double xmin = 0.0;

  if ( nlnrel == 0 )
  {
    nlnrel = r8_inits ( alnrcs, 43, 0.1 * r8_mach ( 3 ) );
    xmin = - 1.0 + sqrt ( r8_mach ( 4 ) );
  }

  if ( x <= - 1.0 )
  {
    cerr << "\n";
    cerr << "R8_LNREL - Fatal error!\n";
    cerr << "  X <= -1.\n";
    exit ( 1 );
  }
  else if ( x < xmin )
  {
    cerr << "\n";
    cerr << "R8_LNREL - Warning!\n";
    cerr << "  Result is less than half precision.\n";
    cerr << "  X is too close to - 1.\n";
  }

  if ( r8_abs ( x ) <= 0.375 )
  {
    value = x * ( 1.0 - x * r8_csevl ( x / 0.375, alnrcs, nlnrel ) );
  }
  else
  {
    value = log ( 1.0 + x );
  }

  return value;
}
//****************************************************************************80

double r8_log ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LOG evaluates the logarithm of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double R8_LOG, the logarithm of X.
//
{
  static double aln2 = 0.06814718055994530941723212145818;
  static double alncen[5] = {
     0.0,
    +0.22314355131420975576629509030983,
    +0.40546510810816438197801311546434,
    +0.55961578793542268627088850052682,
    +0.69314718055994530941723212145817 };
  static double alncs[11] = {
    +0.13347199877973881561689386047187E+01,
    +0.69375628328411286281372438354225E-03,
    +0.42934039020450834506559210803662E-06,
    +0.28933847795432594580466440387587E-09,
    +0.20512517530340580901741813447726E-12,
    +0.15039717055497386574615153319999E-15,
    +0.11294540695636464284521613333333E-18,
    +0.86355788671171868881946666666666E-22,
    +0.66952990534350370613333333333333E-25,
    +0.52491557448151466666666666666666E-28,
    +0.41530540680362666666666666666666E-31 };
  static double center[4] = {
    1.0,
    1.25,
    1.50,
    1.75 };
  int n;
  static int nterms = 0;
  int ntrval;
  double t;
  double t2;
  double value;
  double xn;
  double y;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( alncs, 11, 28.9 * r8_mach ( 3 ) );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_LOG - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  r8_upak ( x, y, n );

  xn = ( double ) ( n - 1 );
  y = 2.0 * y;
  ntrval = ( int ) ( 4.0 * y - 2.5 );

  if ( ntrval == 5 )
  {
    t = ( ( y - 1.0 ) - 1.0 ) / ( y + 2.0 );
  }
  else if ( ntrval < 5 )
  {
    t = ( y - center[ntrval-1] ) / ( y + center[ntrval-1] );
  }

  t2 = t * t;
  value = 0.625 * xn + ( aln2 * xn + alncen[ntrval-1] + 2.0 * t 
    + t * t2 * r8_csevl ( 578.0 * t2 - 1.0, alncs, nterms) );

  return value;
}
//****************************************************************************80

double r8_log10 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LOG10 evaluates the logarithm, base 10, of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double R8_LOG10, the logarithm, base 10, of X.
//
{
  static double aloge = 0.43429448190325182765112891891661;
  double value;

  value = aloge * log ( x );

  return value;
}
//****************************************************************************80

double r8_mach ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MACH returns double precision real machine constants.
//
//  Discussion:
//
//    Assuming that the internal representation of a double precision real
//    number is in base B, with T the number of base-B digits in the mantissa,
//    and EMIN the smallest possible exponent and EMAX the largest possible 
//    exponent, then
//
//      R8_MACH(1) = B^(EMIN-1), the smallest positive magnitude.
//      R8_MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
//      R8_MACH(3) = B^(-T), the smallest relative spacing.
//      R8_MACH(4) = B^(1-T), the largest relative spacing.
//      R8_MACH(5) = log10(B).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, Norman Schryer,
//    Algorithm 528:
//    Framework for a Portable Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, Number 2, June 1978, page 176-188.
//
//  Parameters:
//
//    Input, int I, chooses the parameter to be returned.
//    1 <= I <= 5.
//
//    Output, double R8_MACH, the value of the chosen parameter.
//
{
  double value;

  if ( i == 1 )
  {
    value = 4.450147717014403E-308;
  }
  else if ( i == 2 )
  {
    value = 8.988465674311579E+307;
  }
  else if ( i == 3 )
  {
    value = 1.110223024625157E-016;
  }
  else if ( i == 4 )
  {
    value = 2.220446049250313E-016;
  }
  else if ( i == 5 )
  {
    value = 0.301029995663981E+000;
  }
  else
  {
    cerr << "\n";
    cerr << "R8_MACH - Fatal error!\n";
    cerr << "  The input argument I is out of bounds.\n";
    cerr << "  Legal values satisfy 1 <= I <= 5.\n";
    cerr << "  I = " << i << "\n";
    value = 0.0;
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

void r8_machar ( long int *ibeta, long int *it, long int *irnd, long int *ngrd,
  long int *machep, long int *negep, long int *iexp, long int *minexp,
  long int *maxexp, double *eps, double *epsneg, double *xmin, double *xmax )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MACHAR computes machine constants for R8 arithmetic.
//
//  Discussion:
//
//    This routine determines the parameters of the floating-point 
//    arithmetic system specified below.  The determination of the first 
//    three uses an extension of an algorithm due to Malcolm, 
//    incorporating some of the improvements suggested by Gentleman and 
//    Marovich.  
//
//    A FORTRAN version of this routine appeared as ACM algorithm 665.
//
//    This routine is a C translation of the FORTRAN code, and appeared
//    as part of ACM algorithm 722.
//
//    An earlier version of this program was published in Cody and Waite.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2006
//
//  Author:
//
//    Original FORTRAN77 version by William Cody.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
//      machine parameters,
//    ACM Transactions on Mathematical Software,
//    Volume 14, Number 4, pages 303-311, 1988.
//
//    William Cody and W Waite,
//    Software Manual for the Elementary Functions,
//    Prentice Hall, 1980.
//
//    M Gentleman and S Marovich,
//    Communications of the ACM,
//    Volume 17, pages 276-277, 1974.
//
//    M. Malcolm,
//    Communications of the ACM,
//    Volume 15, pages 949-951, 1972.
//
//  Parameters:
//
//    Output, long int* ibeta, the radix for the floating-point representation.
//
//    Output, long int* it, the number of base IBETA digits in the floating-point
//    significand.
//
//    Output, long int* irnd:
//    0, if floating-point addition chops.
//    1, if floating-point addition rounds, but not in the IEEE style.
//    2, if floating-point addition rounds in the IEEE style.
//    3, if floating-point addition chops, and there is partial underflow.
//    4, if floating-point addition rounds, but not in the IEEE style, and 
//      there is partial underflow.
//    5, if floating-point addition rounds in the IEEE style, and there is 
//      partial underflow.
//
//    Output, long int* ngrd, the number of guard digits for multiplication with
//    truncating arithmetic.  It is
//    0, if floating-point arithmetic rounds, or if it truncates and only 
//      IT base IBETA digits participate in the post-normalization shift of the
//      floating-point significand in multiplication;
//   1, if floating-point arithmetic truncates and more than IT base IBETA
//      digits participate in the post-normalization shift of the floating-point
//      significand in multiplication.
//
//    Output, long int* MACHEP, the largest negative integer such that
//      1.0 + ( double ) IBETA ^ MACHEP != 1.0, 
//    except that MACHEP is bounded below by - ( IT + 3 ).
//
//    Output, long int* NEGEPS, the largest negative integer such that
//      1.0 - ( double ) IBETA ) ^ NEGEPS != 1.0, 
//    except that NEGEPS is bounded below by - ( IT + 3 ).
//
//    Output, long int* IEXP, the number of bits (decimal places if IBETA = 10)
//    reserved for the representation of the exponent (including the bias or
//    sign) of a floating-point number.
//
//    Output, long int* MINEXP, the largest in magnitude negative integer such 
//    that
//      ( double ) IBETA ^ MINEXP 
//    is positive and normalized.
//
//    Output, long int* MAXEXP, the smallest positive power of BETA that overflows.
// 
//    Output, double* EPS, the smallest positive floating-point number such
//    that  
//      1.0 + EPS != 1.0. 
//    in particular, if either IBETA = 2  or IRND = 0, 
//      EPS = ( double ) IBETA ^ MACHEP.
//    Otherwise,  
//      EPS = ( ( double ) IBETA ^ MACHEP ) / 2.
//
//    Output, double* EPSNEG, a small positive floating-point number such that
//      1.0 - EPSNEG != 1.0. 
//    In particular, if IBETA = 2 or IRND = 0, 
//      EPSNEG = ( double ) IBETA ^ NEGEPS.
//    Otherwise,  
//      EPSNEG = ( double ) IBETA ^ NEGEPS ) / 2.  
//    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
//    smallest number that can alter 1.0 by subtraction.
//
//    Output, double* XMIN, the smallest non-vanishing normalized floating-point
//    power of the radix:
//      XMIN = ( double ) IBETA ^ MINEXP
//
//    Output, float* XMAX, the largest finite floating-point number.  In
//    particular,
//      XMAX = ( 1.0 - EPSNEG ) * ( double ) IBETA ^ MAXEXP
//    On some machines, the computed value of XMAX will be only the second, 
//    or perhaps third, largest number, being too small by 1 or 2 units in 
//    the last digit of the significand.
//
{
  double a;
  double b;
  double beta;
  double betah;
  double betain;
  int i;
  int itmp;
  int iz;
  int j;
  int k;
  int mx;
  int nxres;
  double one;
  double t;
  double tmp;
  double tmp1;
  double tmpa;
  double two;
  double y;
  double z;
  double zero;

  (*irnd) = 1;
  one = (double) (*irnd);
  two = one + one;
  a = two;
  b = a;
  zero = 0.0e0;
//
//  Determine IBETA and BETA ala Malcolm.
//
  tmp = ( ( a + one ) - a ) - one;

  while ( tmp == zero )
  {
    a = a + a;
    tmp = a + one;
    tmp1 = tmp - a;
    tmp = tmp1 - one;
  }

  tmp = a + b;
  itmp = ( int ) ( tmp - a );

  while ( itmp == 0 )
  {
    b = b + b;
    tmp = a + b;
    itmp = ( int ) ( tmp - a );
  }

  *ibeta = itmp;
  beta = ( double ) ( *ibeta );
//
//  Determine IRND, IT.
//
  ( *it ) = 0;
  b = one;
  tmp = ( ( b + one ) - b ) - one;

  while ( tmp == zero )
  {
    *it = *it + 1;
    b = b * beta;
    tmp = b + one;
    tmp1 = tmp - b;
    tmp = tmp1 - one;
  }

  *irnd = 0;
  betah = beta / two;
  tmp = a + betah;
  tmp1 = tmp - a;

  if ( tmp1 != zero )
  {
    *irnd = 1;
  }

  tmpa = a + beta;
  tmp = tmpa + betah;

  if ( ( *irnd == 0 ) && ( tmp - tmpa != zero ) )
  {
    *irnd = 2;
  }
//
//  Determine NEGEP, EPSNEG.
//
  (*negep) = (*it) + 3;
  betain = one / beta;
  a = one;
 
  for ( i = 1; i <= (*negep); i++ )
  {
    a = a * betain;
  }
 
  b = a;
  tmp = ( one - a );
  tmp = tmp - one;

  while ( tmp == zero )
  {
    a = a * beta;
    *negep = *negep - 1;
    tmp1 = one - a;
    tmp = tmp1 - one;
  }

  (*negep) = -(*negep);
  (*epsneg) = a;
//
//  Determine MACHEP, EPS.
//

  (*machep) = -(*it) - 3;
  a = b;
  tmp = one + a;

  while ( tmp - one == zero)
  {
    a = a * beta;
    *machep = *machep + 1;
    tmp = one + a;
  }

  *eps = a;
//
//  Determine NGRD.
//
  (*ngrd) = 0;
  tmp = one + *eps;
  tmp = tmp * one;

  if ( ( (*irnd) == 0 ) && ( tmp - one ) != zero )
  {
    (*ngrd) = 1;
  }
//
//  Determine IEXP, MINEXP and XMIN.
//
//  Loop to determine largest I such that (1/BETA) ** (2**(I))
//  does not underflow.  Exit from loop is signaled by an underflow.
//

  i = 0;
  k = 1;
  z = betain;
  t = one + *eps;
  nxres = 0;

  for ( ; ; )
  {
    y = z;
    z = y * y;
//
//  Check for underflow
//

    a = z * one;
    tmp = z * t;

    if ( ( a + a == zero ) || ( r8_abs ( z ) > y ) )
    {
      break;
    }

    tmp1 = tmp * betain;

    if ( tmp1 * beta == z )
    {
      break;
    }

    i = i + 1;
    k = k + k;
  }
//
//  Determine K such that (1/BETA)**K does not underflow.
//  First set  K = 2 ^ I.
//
  (*iexp) = i + 1;
  mx = k + k;
//
//  For decimal machines only
//
  if ( *ibeta == 10 )
  {
    (*iexp) = 2;
    iz = *ibeta;
    while ( iz <= k )
    {
      iz = iz * ( *ibeta );
      (*iexp) = (*iexp) + 1;
    }
    mx = iz + iz - 1;
  } 
//
//  Loop to determine MINEXP, XMIN.
//  Exit from loop is signaled by an underflow.
//
  for ( ; ; )
  {
    (*xmin) = y;
    y = y * betain;
    a = y * one;
    tmp = y * t;
    tmp1 = a + a;

    if ( ( tmp1 == zero ) || ( r8_abs ( y ) >= ( *xmin ) ) )
    {
      break;
    }

    k = k + 1;
    tmp1 = tmp * betain;
    tmp1 = tmp1 * beta;

    if ( ( tmp1 == y ) && ( tmp != y ) )
    {
      nxres = 3;
      *xmin = y;
      break;
    }

  }

  (*minexp) = -k;
//
//  Determine MAXEXP, XMAX.
//
  if ( ( mx <= k + k - 3 ) && ( ( *ibeta ) != 10 ) )
  {
    mx = mx + mx;
    (*iexp) = (*iexp) + 1;
  }

  (*maxexp) = mx + (*minexp);
//
//  Adjust IRND to reflect partial underflow.
//
  (*irnd) = (*irnd) + nxres;
//
//  Adjust for IEEE style machines.
//
  if ( ( *irnd) >= 2 )
  {
    (*maxexp) = (*maxexp) - 2;
  }
//
//  Adjust for machines with implicit leading bit in binary
//  significand and machines with radix point at extreme
//  right of significand.
//
  i = (*maxexp) + (*minexp);

  if ( ( ( *ibeta ) == 2 ) && ( i == 0 ) )
  {
    (*maxexp) = (*maxexp) - 1;
  }

  if ( i > 20 )
  {
    (*maxexp) = (*maxexp) - 1;
  }

  if ( a != y )
  {
    (*maxexp) = (*maxexp) - 2;
  }

  (*xmax) = one - (*epsneg);
  tmp = (*xmax) * one;

  if ( tmp != (*xmax) )
  {
    (*xmax) = one - beta * (*epsneg);
  }

  (*xmax) = (*xmax) / ( beta * beta * beta * (*xmin) );
  i = (*maxexp) + (*minexp) + 3;

  if ( i > 0 )
  {
 
    for ( j = 1; j <= i; j++ )
    {
      if ( (*ibeta) == 2 )
      {
        (*xmax) = (*xmax) + (*xmax);
      }
      if ( (*ibeta) != 2 )
      {
        (*xmax) = (*xmax) * beta;
      }
    }

  }
  return;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double r8_mod ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOD returns the remainder of R8 division.
//
//  Discussion:
//
//    If
//      REM = R8_MOD ( X, Y )
//      RMULT = ( X - REM ) / Y
//    then
//      X = Y * RMULT + REM
//    where REM has the same sign as X, and abs ( REM ) < Y.
//
//  Example:
//
//        X         Y     R8_MOD   R8_MOD  Factorization
//
//      107        50       7     107 =  2 *  50 + 7
//      107       -50       7     107 = -2 * -50 + 7
//     -107        50      -7    -107 = -2 *  50 - 7
//     -107       -50      -7    -107 =  2 * -50 - 7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number to be divided.
//
//    Input, double Y, the number that divides X.
//
//    Output, double R8_MOD, the remainder when X is divided by Y.
//
{
  double value;

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R8_MOD - Fatal error!\n";
    cerr << "  R8_MOD ( X, Y ) called with Y = " << y << "\n";
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( x < 0.0 && 0.0 < value )
  {
    value = value - r8_abs ( y );
  }
  else if ( 0.0 < x && value < 0.0 )
  {
    value = value + r8_abs ( y );
  }

  return value;
}
//****************************************************************************80

double r8_mop ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOP returns the I-th power of -1 as an R8 value.
//
//  Discussion:
//
//    An R8 is an double value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, double R8_MOP, the I-th power of -1.
//
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}
//****************************************************************************80

double r8_pak ( double y, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PAK packs a base 2 exponent into an R8.
//
//  Discussion:
//
//    This routine is almost the inverse of R8_UPAK.  It is not exactly 
//    the inverse, because abs(x) need not be between 0.5 and 1.0.  
//    If both R8_PAK and 2.0^n were known to be in range, we could compute
//    R8_PAK = x * 2.0^n .
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double Y, the mantissa.
//
//    Input, int N, the exponent.
//
//    Output, double R8_PAK, the packed value.
//
{
  static double aln210 = 3.321928094887362347870319429489;
  double aln2b;
  static int nmax = 0;
  static int nmin = 0;
  int nsum;
  int ny;
  double value;

  if ( nmin == 0 )
  {
    aln2b = 1.0;
    if ( i4_mach ( 10 ) != 2 );
    {
      aln2b = r8_mach ( 5 ) * aln210;
    }
    nmin = aln2b * ( double ) ( i4_mach ( 15 ) );
    nmax = aln2b * ( double ) ( i4_mach ( 16 ) );
  }

  r8_upak ( y, value, ny );

  nsum = n + ny;

  if ( nsum < nmin )
  {
    cerr << "\n";
    cerr << "R8_PAK - Warning!\n";
    cerr << "  Packed number underflows.\n";
    value = 0.0;
    return value;
  }

  if ( nmax < nsum )
  {
    cerr << "\n";
    cerr << "R8_PAK - Fatal error!\n";
    cerr << "  Packed number overflows.\n";
    exit ( 1 );
  }

  while ( nsum < 0 )
  {
    value = 0.5 * value;
    nsum = nsum + 1;
  }

  while ( 0 < nsum )
  {
    value = 2.0 * value;
    nsum = nsum - 1;
  }
  return value;
}
//****************************************************************************80

double r8_poch ( double a, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_POCH evaluates Pochhammer's function of R8 arguments.
//
//  Discussion:
//
//    POCH ( A, X ) = Gamma ( A + X ) / Gamma ( A ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, X, the arguments.
//
//    Output, double R8_POCH, the Pochhammer function of A and X.
//
{
  double absa;
  double absax;
  double alnga;
  double alngax;
  double ax;
  double b;
  double cospia;
  double cospix;
  double den;
  static double eps = 0.0;
  double err;
  double errpch;
  int i;
  int ia;
  int n;
  static double pi = 3.141592653589793238462643383279503;
  double sgnga;
  double sgngax;
  double sinpia;
  double sinpix;
  static double sqeps = 0.0;
  double value;

  if ( eps == 0.0 )
  {
    eps = r8_mach ( 4 );
    sqeps = sqrt ( eps );
  }

  ax = a + x;

  if ( ax <= 0.0 && r8_aint ( ax ) == ax )
  {
    if ( 0.0 < a || r8_aint ( a ) != a )
    {
      cerr << "\n";
      cerr << "R8_POCH - Fatal error!\n";
      cerr << "  A + X is nonpositive int,\n";
      cerr << "  but A is not.\n";
      exit ( 1 );
    }
//
//  We know here that both A+X and A are non-positive integers.
//
    if ( x == 0.0 )
    {
      value = 1.0;
    }
    else if ( - 20.0 < r8_min ( a + x, a ) )
    {
      n = ( int ) ( x );
      ia = ( int ) ( a );
      value = r8_mop ( n ) * r8_fac ( - ia ) / r8_fac ( - ia - n );
    }
    else
    {
      n = ( int ) ( x );
      value = r8_mop ( n ) * exp ( ( a - 0.5 ) 
        * r8_lnrel ( x / ( a - 1.0 ) ) 
        + x * log ( - a + 1.0 - x ) - x 
        + r8_lgmc ( - a + 1.0 ) 
        - r8_lgmc ( - a - x + 1.0 ) );
    }
    return value;
  }
//
//  A + X is not zero or a negative integer.
//
  if ( a <= 0.0 && r8_aint ( a ) == a )
  {
    value = 0.0;
    return value;
  }

  n = r8_abs ( x );
//
//  x is a small non-positive integer, presummably a common case.
//
  if ( ( double ) ( n ) == x && n <= 20 )
  {
    value = 1.0;
    for ( i = 1; i <= n; i++ )
    {
      value = value * ( a + ( double ) ( i - 1 ) );
    }
    return value;
  }

  absax = r8_abs ( a + x );
  absa = r8_abs ( a );

  if ( r8_max ( absax, absa ) <= 20.0 )
  {
    value = r8_gamma ( a + x ) * r8_gamr ( a );
    return value;
  }

  if ( 0.5 * absa < r8_abs ( x ) )
  {
    r8_lgams ( a + x, alngax, sgngax );
    r8_lgams ( a, alnga, sgnga );
    value = sgngax * sgnga * exp ( alngax - alnga );
    return value;
  }
//
//  abs(x) is small and both abs(a+x) and abs(a) are large.  thus,
//  a+x and a must have the same sign.  for negative a, we use
//  gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) *
//  sin(pi*a)/sin(pi*(a+x))
//
  if ( a < 0.0 )
  {
    b = - a - x + 1.0;
  }
  else
  {
    b = a;
  }

  value = exp ( ( b - 0.5 ) * r8_lnrel ( x / b ) 
    + x * log ( b + x ) - x + r8_lgmc ( b + x ) - r8_lgmc ( b ) );

  if ( 0.0 <= a || value == 0.0 )
  {
    return value;
  }

  cospix = cos ( pi * x );
  sinpix = sin ( pi * x );
  cospia = cos ( pi * a );
  sinpia = sin ( pi * a );

  errpch = r8_abs ( x ) * ( 1.0 + log ( b ) );
  den = cospix + cospia * sinpix / sinpia;
  err = ( r8_abs ( x ) * ( r8_abs ( sinpix ) 
    + r8_abs ( cospia * cospix / sinpia ) )
    + r8_abs ( a * sinpix ) / sinpia / sinpia ) * pi;
  err = errpch + err / r8_abs ( den );

  value = value / den;

  return value;
}
//****************************************************************************80

double r8_poch1 ( double a, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_POCH1 evaluates a quantity related to Pochhammer's symbol.
//
//  Discussion:
//
//    Evaluate a generalization of Pochhammer's symbol for special
//    situations that require especially accurate values when x is small in
//      poch1(a,x) = (poch(a,x)-1)/x
//                 = (gamma(a+x)/gamma(a) - 1.0)/x .
//    This specification is particularly suited for stably computing
//    expressions such as
//      (gamma(a+x)/gamma(a) - gamma(b+x)/gamma(b))/x
//           = poch1(a,x) - poch1(b,x)
//    Note that poch1(a,0.0) = psi(a)
//
//    When abs(x) is so small that substantial cancellation will occur if
//    the straightforward formula is used, we  use an expansion due
//    to fields and discussed by y. l. luke, the special functions and their
//    approximations, vol. 1, academic press, 1969, page 34.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double A, the parameter.
//
//    Input, double X, the evaluation point.
//
//    Output, double R8_POCH1, the value of the function.
//
{
  double absa;
  double absx;
  static double alneps = 0.0;
  double alnvar;
  double b;
  static double bern[20] = {
    +0.833333333333333333333333333333333E-01,
    -0.138888888888888888888888888888888E-02,
    +0.330687830687830687830687830687830E-04,
    -0.826719576719576719576719576719576E-06,
    +0.208767569878680989792100903212014E-07,
    -0.528419013868749318484768220217955E-09,
    +0.133825365306846788328269809751291E-10,
    -0.338968029632258286683019539124944E-12,
    +0.858606205627784456413590545042562E-14,
    -0.217486869855806187304151642386591E-15,
    +0.550900282836022951520265260890225E-17,
    -0.139544646858125233407076862640635E-18,
    +0.353470703962946747169322997780379E-20,
    -0.895351742703754685040261131811274E-22,
    +0.226795245233768306031095073886816E-23,
    -0.574472439520264523834847971943400E-24,
    +0.145517247561486490186626486727132E-26,
    -0.368599494066531017818178247990866E-28,
    +0.933673425709504467203255515278562E-30,
    -0.236502241570062993455963519636983E-31 };
  double binv;
  double bp;
  double gbern[21];
  double gbk;
  int i;
  int ii;
  int incr;
  int j;
  int k;
  int ndx;
  int nterms;
  static double pi = 3.141592653589793238462643383279503;
  double poly1;
  double q;
  double rho;
  double sinpxx;
  double sinpx2;
  static double sqtbig = 0.0;
  double term;
  double trig;
  double value;
  double var;
  double var2;

  if ( sqtbig == 0.0 )
  {
    sqtbig = 1.0 / sqrt ( 24.0 * r8_mach ( 1 ) );
    alneps = log ( r8_mach ( 3 ) );
  }

  if ( x == 0.0 )
  {
    value = r8_psi ( a );
    return value;
  }

  absx = r8_abs ( x );
  absa = r8_abs ( a );

  if ( 0.1 * absa < absx || 0.1 < absx * log ( r8_max ( absa, 2.0 ) ) )
  {
    value = r8_poch ( a, x );
    value = ( value - 1.0 ) / x;
    return value;
  }

  if ( a < - 0.5 )
  {
    bp = 1.0 - a - x;
  }
  else
  {
    bp = a;
  }

  if ( bp < 10.0 )
  {
    incr = r8_aint ( 11.0 - bp );
  }
  else
  {
    incr = 0;
  }

  b = bp + ( double ) ( incr );

  var = b + 0.5 * ( x - 1.0 );
  alnvar = log ( var );
  q = x * alnvar;
  poly1 = 0.0;

  if ( var < sqtbig )
  {
    var2 = 1.0 / var / var;

    rho = 0.5 * ( x + 1.0 );
    gbern[0] = 1.0;
    gbern[1] = - rho / 12.0;
    term = var2;
    poly1 = gbern[1] * term;

    nterms = ( int ) ( - 0.5 * alneps / alnvar + 1.0 );

    if ( 20 < nterms )
    {
      cerr << "\n";
      cerr << "R8_POCH1 - Fatal error!\n";
      cerr << " 20 < NTERMS.\n";
      exit ( 1 );
    } 

    for ( k = 2; k <= nterms; k++ )
    {
      gbk = 0.0;
      for ( j = 1; j <= k; j++ )
      {
        ndx = k - j + 1;
        gbk = gbk + bern[ndx-1] * gbern[j-1];
      }
      gbern[k] = - rho * gbk / ( double ) ( k );
      term = term * ( ( double ) ( 2 * k - 2 ) - x )
        * ( ( double ) ( 2 * k - 1 ) - x ) * var2;
      poly1 = poly1 + gbern[k] * term;
    }
  }
  poly1 = ( x - 1.0 ) * poly1;
  value = r8_exprel ( q ) * ( alnvar + q * poly1 ) + poly1;
//
//  we have r8_poch1(b,x), but bp is small, so we use backwards recursion
//  to obtain r8_poch1(bp,x).
//
  for ( ii = 1; ii <= incr; ii++ )
  {
    i = incr - ii;
    binv = 1.0 / ( bp + ( double ) ( i ) );
    value = ( value - binv ) / ( 1.0 + x * binv );
  }

  if ( bp == a )
  {
    return value;
  }
//
//  we have r8_poch1(bp,x), but a is lt -0.5.  we therefore use a reflection
//  formula to obtain r8_poch1(a,x).
//
  sinpxx = sin ( pi * x ) / x;
  sinpx2 = sin ( 0.5 * pi * x );
  trig = sinpxx * r8_cot ( pi * b ) - 2.0 * sinpx2 * ( sinpx2 / x );

  value = trig + ( 1.0 + x * trig ) * value;

  return value;
}
//****************************************************************************80

double r8_power ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    R8_POWER evaluates A^B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the base.
//
//    Input, double B, the exponent.
//
//    Output, double R8_POWER, the value of A^B.
//
{
  double value;

  value = pow ( a, b );

  return value;
}
//****************************************************************************80

double r8_psi ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PSI evaluates the psi function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_PSI, the psi function of X.
//
{
  static double apsics[16] = {
    -0.832710791069290760174456932269E-03,
    -0.416251842192739352821627121990E-03,
    +0.103431560978741291174463193961E-06,
    -0.121468184135904152987299556365E-09,
    +0.311369431998356155521240278178E-12,
    -0.136461337193177041776516100945E-14,
    +0.902051751315416565130837974000E-17,
    -0.831542997421591464829933635466E-19,
    +0.101224257073907254188479482666E-20,
    -0.156270249435622507620478933333E-22,
    +0.296542716808903896133226666666E-24,
    -0.674686886765702163741866666666E-26,
    +0.180345311697189904213333333333E-27,
    -0.556901618245983607466666666666E-29,
    +0.195867922607736251733333333333E-30,
    -0.775195892523335680000000000000E-32 };
  double aux;
  static double dxrel = 0.0;
  int i;
  int n;
  static int ntapsi = 0;
  static int ntpsi = 0;
  static double pi = 3.14159265358979323846264338327950;
  static double psics[42] = {
    -0.38057080835217921520437677667039E-01,
    +0.49141539302938712748204699654277,
    -0.56815747821244730242892064734081E-01,
    +0.83578212259143131362775650747862E-02,
    -0.13332328579943425998079274172393E-02,
    +0.22031328706930824892872397979521E-03,
    -0.37040238178456883592889086949229E-04,
    +0.62837936548549898933651418717690E-05,
    -0.10712639085061849855283541747074E-05,
    +0.18312839465484165805731589810378E-06,
    -0.31353509361808509869005779796885E-07,
    +0.53728087762007766260471919143615E-08,
    -0.92116814159784275717880632624730E-09,
    +0.15798126521481822782252884032823E-09,
    -0.27098646132380443065440589409707E-10,
    +0.46487228599096834872947319529549E-11,
    -0.79752725638303689726504797772737E-12,
    +0.13682723857476992249251053892838E-12,
    -0.23475156060658972717320677980719E-13,
    +0.40276307155603541107907925006281E-14,
    -0.69102518531179037846547422974771E-15,
    +0.11856047138863349552929139525768E-15,
    -0.20341689616261559308154210484223E-16,
    +0.34900749686463043850374232932351E-17,
    -0.59880146934976711003011081393493E-18,
    +0.10273801628080588258398005712213E-18,
    -0.17627049424561071368359260105386E-19,
    +0.30243228018156920457454035490133E-20,
    -0.51889168302092313774286088874666E-21,
    +0.89027730345845713905005887487999E-22,
    -0.15274742899426728392894971904000E-22,
    +0.26207314798962083136358318079999E-23,
    -0.44964642738220696772598388053333E-24,
    +0.77147129596345107028919364266666E-25,
    -0.13236354761887702968102638933333E-25,
    +0.22709994362408300091277311999999E-26,
    -0.38964190215374115954491391999999E-27,
    +0.66851981388855302310679893333333E-28,
    -0.11469986654920864872529919999999E-28,
    +0.19679385886541405920515413333333E-29,
    -0.33764488189750979801907200000000E-30,
    +0.57930703193214159246677333333333E-31 };
  double value;
  static double xbig = 0.0;
  double y;

  if ( ntpsi == 0 )
  {
    ntpsi = r8_inits ( psics, 42, 0.1 * r8_mach ( 3 ) );
    ntapsi = r8_inits ( apsics, 16, 0.1 * r8_mach ( 3 ) );
    xbig = 1.0 / sqrt ( r8_mach ( 3 ) );
    dxrel = sqrt ( r8_mach ( 4 ) );
  }

  y = r8_abs ( x );

  if ( y < 10.0 )
  {
    n = ( int ) ( x );
    if ( x < 0.0 )
    {
      n = n - 1;
    }
    y = x - ( double ) ( n );
    n = n - 1;
    value = r8_csevl ( 2.0 * y - 1.0, psics, ntpsi );

    if ( n == 0 )
    {
      return value;
    }
    else if ( n < 0 )
    {
      n = - n;

      if ( x == 0.0 )
      {
        cerr << "\n";
        cerr << "R8_PSI - Fatal error!\n";
        cerr << "  X is zero.\n";
        exit ( 1 );
      }

      if ( x < 0.0 && x + ( double ) ( n - 2 ) == 0.0 )
      {
        cerr << "\n";
        cerr << "R8_PSI - Fatal error!\n";
        cerr << "  X is a negative int.\n";
        exit ( 1 );
      }

      if ( x < - 0.5 && r8_abs ( ( x - r8_aint ( x - 0.5 ) ) / x ) < dxrel )
      {
        cerr << "\n";
        cerr << "R8_PSI - Warning!\n";
        cerr << "  Answer is less than half precision\n";
        cerr << "  because X is near a negative int.\n";
      }

      for ( i = 1; i <= n; i++ )
      {
        value = value - 1.0 / ( x + ( double ) ( i - 1 ) );
      }
    }
    else if ( 0 < n )
    {
      for ( i = 1; i <= n; i++ )
      {
        value = value + 1.0 / ( y + ( double ) ( i ) );
      }

    }
  }
  else
  {
    if ( y < xbig )
    {
      aux = r8_csevl ( 8.0 / y / y - 1.0, apsics, ntapsi );
    }
    else
    {
      aux = 0.0;
    }

    if ( x < 0.0 )
    {
      value = log ( r8_abs ( x ) ) - 0.5 / x + aux 
        - pi * r8_cot ( pi * x );
    }
    else if ( 0.0 < x )
    {
      value = log ( x ) - 0.5 / x + aux;
    }
  }
  return value;
}
//****************************************************************************80

double r8_ren ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_REN is a simple random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Malcolm Pike, David Hill,
//    Algorithm 266:
//    Pseudo-Random Numbers,
//    Communications of the ACM,
//    Volume 8, Number 10, October 1965, page 605.
//
//  Parameters:
//
//    Output, double R8_REN, the random value.
//
{
  static int iy = 100001;
  double value;

  iy = iy * 125;
  iy = iy - ( iy / 2796203 ) * 2796203;
  value = ( double ) ( iy ) / 2796203.0;

  return value;
}
//****************************************************************************80

double r8_shi ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SHI evaluates the hyperbolic sine integral Shi of an R8 argument.
//
//  Discussion:
//
//    Shi ( x ) = Integral ( 0 <= t <= x ) sinh ( t ) dt / t
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_SHI, the hyperbolic sine integral 
//    Shi evaluated at X.
//
{
  double absx;
  static int nshi = 0;
  static double shics[10] = {
     0.0078372685688900950695200984317332E+00,
     0.0039227664934234563972697574427225E+00,
     0.0000041346787887617266746747908275E+00,
     0.0000000024707480372882742135145302E+00,
     0.0000000000009379295590763630457157E+00,
     0.0000000000000002451817019520867353E+00,
     0.0000000000000000000467416155257592E+00,
     0.0000000000000000000000067803072389E+00,
     0.0000000000000000000000000007731289E+00,
     0.0000000000000000000000000000000711E+00 };
  double value;
  static double xsml = 0.0;

  if ( nshi == 0 )
  {
    nshi = r8_inits ( shics, 10, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( r8_mach ( 3 ) );
  }

  absx = r8_abs ( x );

  if ( absx <= xsml )
  {
    value = x;
  }
  else if ( absx <= 0.375 )
  {
    value = x * ( 1.0 
      + r8_csevl ( 128.0 * x * x / 9.0 - 1.0, shics, nshi ) );
  }
  else
  {
    value = 0.5 * ( r8_ei ( x ) + r8_e1 ( x ) );
  }
  return value;
}
//****************************************************************************80

double r8_si ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SI evaluates the sine integral Si of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_SI, the sine integral Si evaluated at X.
//
{
  double absx;
  double cosx;
  double f;
  double g;
  static int nsi = 0;
  static double pi2 = 1.57079632679489661923132169163975;
  static double sics[18] = {
    -0.1315646598184841928904275173000457,
    -0.2776578526973601892048287660157299,
     0.0354414054866659179749135464710086,
    -0.0025631631447933977658752788361530,
     0.0001162365390497009281264921482985,
    -0.0000035904327241606042670004347148,
     0.0000000802342123705710162308652976,
    -0.0000000013562997692540250649931846,
     0.0000000000179440721599736775567759,
    -0.0000000000001908387343087145490737,
     0.0000000000000016669989586824330853,
    -0.0000000000000000121730988368503042,
     0.0000000000000000000754181866993865,
    -0.0000000000000000000004014178842446,
     0.0000000000000000000000018553690716,
    -0.0000000000000000000000000075166966,
     0.0000000000000000000000000000269113,
    -0.0000000000000000000000000000000858 };
  double value;
  static double xsml = 0.0;

  if ( nsi == 0 )
  {
    nsi = r8_inits ( sics, 18, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( r8_mach ( 3 ) );
  }

  absx = r8_abs ( x );

  if ( absx < xsml )
  {
    value = x;
  }
  else if ( absx <= 4.0 )
  {
    value = x * ( 0.75 + r8_csevl ( ( x * x - 8.0 ) * 0.125, sics, nsi ) );
  }
  else
  {
    r8_sifg ( absx, f, g );
    cosx = cos ( absx );
    value = pi2 - f * cosx - g * sin ( x );
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

void r8_sifg ( double x, double &f, double &g )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIFG is a utility routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double &F, &G.
//
{
  static double f1cs[43] = {
-0.1191081969051363610348201965828918,
-0.0247823144996236247590074150823133,
 0.0011910281453357821268120363054457,
-0.0000927027714388561748308600360706,
 0.0000093373141568270996868204582766,
-0.0000011058287820557143938979426306,
 0.0000001464772071460162169336550799,
-0.0000000210694496287689532601227548,
 0.0000000032293492366848236382857374,
-0.0000000005206529617529375828014986,
 0.0000000000874878884570278750268316,
-0.0000000000152176187056123668294574,
 0.0000000000027257192405419573900583,
-0.0000000000005007053075968556290255,
 0.0000000000000940240902726068511779,
-0.0000000000000180014444791803678336,
 0.0000000000000035062621432741785826,
-0.0000000000000006935282926769149709,
 0.0000000000000001390925136454216568,
-0.0000000000000000282486885074170585,
 0.0000000000000000058031305693579081,
-0.0000000000000000012046901573375820,
 0.0000000000000000002525052443655940,
-0.0000000000000000000533980268805594,
 0.0000000000000000000113855786274122,
-0.0000000000000000000024462861505259,
 0.0000000000000000000005293659320439,
-0.0000000000000000000001153184940277,
 0.0000000000000000000000252786568318,
-0.0000000000000000000000055738645378,
 0.0000000000000000000000012358245621,
-0.0000000000000000000000002754350842,
 0.0000000000000000000000000616906808,
-0.0000000000000000000000000138817443,
 0.0000000000000000000000000031375329,
-0.0000000000000000000000000007121249,
 0.0000000000000000000000000001622778,
-0.0000000000000000000000000000371206,
 0.0000000000000000000000000000085221,
-0.0000000000000000000000000000019633,
 0.0000000000000000000000000000004538,
-0.0000000000000000000000000000001052,
 0.0000000000000000000000000000000245 };
  static double f2cs[99] = {
-0.03484092538970132330836049733745577,
-0.01668422056779596873246786312278676,
 0.00067529012412377385045207859239727,
-0.00005350666225447013628785577557429,
 0.00000626934217790075267050759431626,
-0.00000095266388019916680677790414293,
 0.00000017456292242509880425504427666,
-0.00000003687954030653093307097646628,
 0.00000000872026777051395264075816938,
-0.00000000226019703919738748530423167,
 0.00000000063246249765250612520444877,
-0.00000000018889118884717869240911480,
 0.00000000005967746729997813372620472,
-0.00000000001980443117372239011196007,
 0.00000000000686413954772103383713264,
-0.00000000000247310193070199106074890,
 0.00000000000092263594549941404196042,
-0.00000000000035523634999261784497297,
 0.00000000000014076049625351591461820,
-0.00000000000005726228499747652794311,
 0.00000000000002386537545413171810106,
-0.00000000000001017141890764597142232,
 0.00000000000000442594531078364424968,
-0.00000000000000196344933049189761979,
 0.00000000000000088688748314810461024,
-0.00000000000000040743345027311546948,
 0.00000000000000019016837215675339859,
-0.00000000000000009009707297478042442,
 0.00000000000000004329211274095668667,
-0.00000000000000002108144465322479526,
 0.00000000000000001039637907026452274,
-0.00000000000000000518891007948931936,
 0.00000000000000000261955324869899371,
-0.00000000000000000133690399951301570,
 0.00000000000000000068941057702931664,
-0.00000000000000000035905362610437250,
 0.00000000000000000018878077255791706,
-0.00000000000000000010016125265594380,
 0.00000000000000000005360725691578228,
-0.00000000000000000002893198974944827,
 0.00000000000000000001574065100202625,
-0.00000000000000000000863027106431206,
 0.00000000000000000000476715602862288,
-0.00000000000000000000265222739998504,
 0.00000000000000000000148582865063866,
-0.00000000000000000000083797235923135,
 0.00000000000000000000047565916422711,
-0.00000000000000000000027169073353112,
 0.00000000000000000000015612738881686,
-0.00000000000000000000009024555078347,
 0.00000000000000000000005246097049119,
-0.00000000000000000000003066450818697,
 0.00000000000000000000001801996250957,
-0.00000000000000000000001064443050752,
 0.00000000000000000000000631942158881,
-0.00000000000000000000000377013812246,
 0.00000000000000000000000225997542918,
-0.00000000000000000000000136100844814,
 0.00000000000000000000000082333232003,
-0.00000000000000000000000050025986091,
 0.00000000000000000000000030526245684,
-0.00000000000000000000000018705164021,
 0.00000000000000000000000011508404393,
-0.00000000000000000000000007108714611,
 0.00000000000000000000000004408065533,
-0.00000000000000000000000002743760867,
 0.00000000000000000000000001714144851,
-0.00000000000000000000000001074768860,
 0.00000000000000000000000000676259777,
-0.00000000000000000000000000426981348,
 0.00000000000000000000000000270500637,
-0.00000000000000000000000000171933331,
 0.00000000000000000000000000109636138,
-0.00000000000000000000000000070132573,
 0.00000000000000000000000000045001784,
-0.00000000000000000000000000028963835,
 0.00000000000000000000000000018697009,
-0.00000000000000000000000000012104646,
 0.00000000000000000000000000007859065,
-0.00000000000000000000000000005116867,
 0.00000000000000000000000000003340627,
-0.00000000000000000000000000002186851,
 0.00000000000000000000000000001435340,
-0.00000000000000000000000000000944523,
 0.00000000000000000000000000000623117,
-0.00000000000000000000000000000412101,
 0.00000000000000000000000000000273208,
-0.00000000000000000000000000000181558,
 0.00000000000000000000000000000120934,
-0.00000000000000000000000000000080737,
 0.00000000000000000000000000000054022,
-0.00000000000000000000000000000036227,
 0.00000000000000000000000000000024348,
-0.00000000000000000000000000000016401,
 0.00000000000000000000000000000011074,
-0.00000000000000000000000000000007497,
 0.00000000000000000000000000000005091,
-0.00000000000000000000000000000003470,
 0.00000000000000000000000000000002377 };
  static double g1cs[44] = {
-0.3040578798253495954499726682091083,
-0.0566890984597120587731339156118269,
 0.0039046158173275643919984071554082,
-0.0003746075959202260618619339867489,
 0.0000435431556559843679552220840065,
-0.0000057417294453025046561970723475,
 0.0000008282552104502629741937616492,
-0.0000001278245892594642727883913223,
 0.0000000207978352948687884439257529,
-0.0000000035313205921990798042032682,
 0.0000000006210824236308951068631449,
-0.0000000001125215474446292649336987,
 0.0000000000209088917684421605267019,
-0.0000000000039715831737681727689158,
 0.0000000000007690431314272089939005,
-0.0000000000001514696742731613519826,
 0.0000000000000302892146552359684119,
-0.0000000000000061399703834708825400,
 0.0000000000000012600605829510933553,
-0.0000000000000002615029250939483683,
 0.0000000000000000548278844891796821,
-0.0000000000000000116038182129526571,
 0.0000000000000000024771654107129795,
-0.0000000000000000005330672753223389,
 0.0000000000000000001155666075598465,
-0.0000000000000000000252280547744957,
 0.0000000000000000000055429038550786,
-0.0000000000000000000012252208421297,
 0.0000000000000000000002723664318684,
-0.0000000000000000000000608707831422,
 0.0000000000000000000000136724874476,
-0.0000000000000000000000030856626806,
 0.0000000000000000000000006995212319,
-0.0000000000000000000000001592587569,
 0.0000000000000000000000000364051056,
-0.0000000000000000000000000083539465,
 0.0000000000000000000000000019240303,
-0.0000000000000000000000000004446816,
 0.0000000000000000000000000001031182,
-0.0000000000000000000000000000239887,
 0.0000000000000000000000000000055976,
-0.0000000000000000000000000000013100,
 0.0000000000000000000000000000003074,
-0.0000000000000000000000000000000723 };
  static double g2cs[44] = {
-0.1211802894731646263541834046858267,
-0.0316761386394950286701407923505610,
 0.0013383199778862680163819429492182,
-0.0000895511011392252425531905069518,
 0.0000079155562961718213115249467924,
-0.0000008438793322241520181418982080,
 0.0000001029980425677530146647227274,
-0.0000000139295750605183835795834444,
 0.0000000020422703959875980400677594,
-0.0000000003196534694206427035434752,
 0.0000000000528147832657267698615312,
-0.0000000000091339554672671033735289,
 0.0000000000016426251238967760444819,
-0.0000000000003055897039322660002410,
 0.0000000000000585655825785779717892,
-0.0000000000000115229197730940120563,
 0.0000000000000023209469119988537310,
-0.0000000000000004774355834177535025,
 0.0000000000000001000996765800180573,
-0.0000000000000000213533778082256704,
 0.0000000000000000046277190777367671,
-0.0000000000000000010175807410227657,
 0.0000000000000000002267657399884672,
-0.0000000000000000000511630776076426,
 0.0000000000000000000116767014913108,
-0.0000000000000000000026935427672470,
 0.0000000000000000000006275665841146,
-0.0000000000000000000001475880557531,
 0.0000000000000000000000350145314739,
-0.0000000000000000000000083757732152,
 0.0000000000000000000000020191815152,
-0.0000000000000000000000004903567705,
 0.0000000000000000000000001199123348,
-0.0000000000000000000000000295170610,
 0.0000000000000000000000000073113112,
-0.0000000000000000000000000018217843,
 0.0000000000000000000000000004565148,
-0.0000000000000000000000000001150151,
 0.0000000000000000000000000000291267,
-0.0000000000000000000000000000074125,
 0.0000000000000000000000000000018953,
-0.0000000000000000000000000000004868,
 0.0000000000000000000000000000001256,
-0.0000000000000000000000000000000325 };
  static double g3cs[56] = {
-0.0280574367809472928402815264335299,
-0.0137271597162236975409100508089556,
 0.0002894032638760296027448941273751,
-0.0000114129239391197145908743622517,
 0.0000006813965590726242997720207302,
-0.0000000547952289604652363669058052,
 0.0000000055207429918212529109406521,
-0.0000000006641464199322920022491428,
 0.0000000000922373663487041108564960,
-0.0000000000144299088886682862611718,
 0.0000000000024963904892030710248705,
-0.0000000000004708240675875244722971,
 0.0000000000000957217659216759988140,
-0.0000000000000207889966095809030537,
 0.0000000000000047875099970877431627,
-0.0000000000000011619070583377173759,
 0.0000000000000002956508969267836974,
-0.0000000000000000785294988256492025,
 0.0000000000000000216922264368256612,
-0.0000000000000000062113515831676342,
 0.0000000000000000018384568838450977,
-0.0000000000000000005610887482137276,
 0.0000000000000000001761862805280062,
-0.0000000000000000000568111050541451,
 0.0000000000000000000187786279582313,
-0.0000000000000000000063531694151124,
 0.0000000000000000000021968802368238,
-0.0000000000000000000007754666550395,
 0.0000000000000000000002791018356581,
-0.0000000000000000000001023178525247,
 0.0000000000000000000000381693403919,
-0.0000000000000000000000144767895606,
 0.0000000000000000000000055779512634,
-0.0000000000000000000000021817239071,
 0.0000000000000000000000008656646309,
-0.0000000000000000000000003482157895,
 0.0000000000000000000000001419188130,
-0.0000000000000000000000000585714314,
 0.0000000000000000000000000244660482,
-0.0000000000000000000000000103387099,
 0.0000000000000000000000000044177299,
-0.0000000000000000000000000019080079,
 0.0000000000000000000000000008326038,
-0.0000000000000000000000000003669553,
 0.0000000000000000000000000001632875,
-0.0000000000000000000000000000733357,
 0.0000000000000000000000000000332327,
-0.0000000000000000000000000000151906,
 0.0000000000000000000000000000070020,
-0.0000000000000000000000000000032539,
 0.0000000000000000000000000000015240,
-0.0000000000000000000000000000007193,
 0.0000000000000000000000000000003420,
-0.0000000000000000000000000000001638,
 0.0000000000000000000000000000000790,
-0.0000000000000000000000000000000383 };
  static int nf1 = 0;
  static int nf2 = 0;
  static int ng1 = 0;
  static int ng2 = 0;
  static int ng3 = 0;
  double tol;
  static double xbig = 0.0;
  static double xbnd = 0.0;
  static double xbndg = 0.0;
  static double xmaxf = 0.0;
  static double xmaxg = 0.0;

  if ( nf1 == 0 )
  {
    tol = 0.1 * r8_mach ( 3 );
    nf1 = r8_inits ( f1cs, 43, tol );
    nf2 = r8_inits ( f2cs, 99, tol );
    ng1 = r8_inits ( g1cs, 44, tol );
    ng2 = r8_inits ( g2cs, 44, tol );
    ng3 = r8_inits ( g3cs, 56, tol );
    xbig = sqrt ( 1.0 / r8_mach ( 3 ) );
    xmaxf = exp ( r8_min ( - log ( r8_mach ( 1 ) ), 
      log ( r8_mach ( 2 ) ) ) - 0.01 );
    xmaxg = 1.0 / sqrt ( r8_mach ( 1 ) );
    xbnd = sqrt ( 50.0 );
    xbndg = sqrt ( 200.0 );
  }

  if ( x < 4.0 )
  {
    cerr << "\n";
    cerr << "R8_SIFG - Fatal error!\n";
    cerr << "  Approximation invalid for X < 4.\n";
    exit ( 1 );
  }
  else if ( x <= xbnd )
  {
    f = ( 1.0 + r8_csevl ( ( 1.0 / x / x - 0.04125 )
      / 0.02125, f1cs, nf1 ) ) / x;
    g = ( 1.0 + r8_csevl ( ( 1.0 / x / x - 0.04125 )
      / 0.02125, g1cs, ng1 ) ) / x / x;
  }
  else if ( x <= xbig )
  {
    f = ( 1.0 + r8_csevl ( 100. / x / x - 1.0, f2cs, nf2 ) ) / x;
    if ( x <= xbndg )
    {
      g = ( 1.0 + r8_csevl ( ( 10000.0 / x / x - 125.0 ) 
        / 75.0, g2cs, ng2 ) ) / x / x;
    }
    else
    {
      g = ( 1.0 + r8_csevl ( 400.0 / x / x - 1.0, g3cs, ng3 ) ) / x / x;
    }
  }
  else
  {
    if ( x < xmaxf )
    {
      f = 1.0 / x;
    }
    else
    {
      f = 0.0;
    }
    if ( x < xmaxg )
    {
      g = 1.0 / x / x;
    }
    else
    {
      g = 0.0;
    }
  }
  return;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

double r8_sin ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIN evaluates the sine of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_SIN, the sine of X.
//
{
  double f;
  int n2;
  static int ntsn = 0;
  static double pi2rec = 0.63661977236758134307553505349006;
  static double pihi = 3.140625;
  static double pilo = 9.6765358979323846264338327950288E-04;
  static double pirec = 0.31830988618379067153776752674503;
  double sgn;
  static double sincs[15] = {
    -0.374991154955873175839919279977323464,
    -0.181603155237250201863830316158004754,
     0.005804709274598633559427341722857921,
    -0.000086954311779340757113212316353178,
     0.000000754370148088851481006839927030,
    -0.000000004267129665055961107126829906,
     0.000000000016980422945488168181824792,
    -0.000000000000050120578889961870929524,
     0.000000000000000114101026680010675628,
    -0.000000000000000000206437504424783134,
     0.000000000000000000000303969595918706,
    -0.000000000000000000000000371357734157,
     0.000000000000000000000000000382486123,
    -0.000000000000000000000000000000336623,
     0.000000000000000000000000000000000256 };
  double value;
  static double xmax = 0.0;
  double xn;
  static double xsml = 0.0;
  static double xwarn = 0.0;
  double y;

  if ( ntsn == 0 )
  {
    ntsn = r8_inits ( sincs, 15, 0.1 * r8_mach ( 3 ) );
    xsml = r8_sqrt ( 2.0 * r8_mach ( 3 ) );
    xmax = 1.0 / r8_mach( 4 );
    xwarn = r8_sqrt ( xmax );
  }

  y = r8_abs ( x );

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R8_SIN - Warning!\n";
    cerr << "  No precision because |X| is big.\n";
    value = 0.0;
    return value;
  }

  if ( xwarn < y )
  {
    cerr << "\n";
    cerr << "R8_SIN - Warning!\n";
    cerr << "  Answer < half precision because |X| is big.\n";
  }

  value = x;
  if ( y < xsml )
  {
    return value;
  }

  xn = ( double ) ( ( int ) ( y * pirec + 0.5 ) );
  n2 = ( int ) ( r8_mod ( xn, 2.0 ) + 0.5 );

  sgn = x;
  if ( n2 != 0 )
  {
    sgn = - sgn;
  }

  f = ( y - xn * pihi ) - xn * pilo;
  xn = 2.0 * ( f * pi2rec ) * ( f * pi2rec ) - 1.0;
  value = f + f * r8_csevl ( xn, sincs, ntsn );

  if ( sgn < 0.0 )
  {
    value = - value;
  }

  if ( value < - 1.0 )
  {
    value = - 1.0;
  }
  else if ( 1.0 < value )
  {
    value = + 1.0;
  }

  return value;
}
//****************************************************************************80

double r8_sin_deg ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIN_DEG evaluates the sine of an R8 argument in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument in degrees.
//
//    Output, double R8_SIN_DEG, the sine of X.
//
{
  int n;
  static double raddeg = 0.017453292519943295769236907684886E+00;
  double value;

  value = sin ( raddeg * x );

  if ( r8_mod ( x, 90.0E+00 ) == 0.0E+00 )
  {
    n = ( int ) ( r8_abs ( x ) / 90.0E+00 + 0.5E+00 );
    n = ( n % 2 );

    if ( n == 0 )
    {
      value = 0.0E+00;
    }
    else if ( value < 0.0E+00 )
    {
      value = - 1.0E+00;
    }
    else
    {
      value = + 1.0E+00;
    }
  }
  return value;
}
//****************************************************************************80

double r8_sinh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SINH evaluates the hyperbolic sine of an R8 argument.
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
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_SINH, the hyperbolic sine of X.
//
{
  static int nterms = 0;
  static double sinhcs[13] = {
    +0.17304219404717963167588384698501E+00,
    +0.87594221922760477154900263454440E-01,
    +0.10794777745671327502427270651579E-02,
    +0.63748492607547504815685554571850E-05,
    +0.22023664049230530159190496019502E-07,
    +0.49879401804158493149425807203661E-10,
    +0.79730535541157304814411480441186E-13,
    +0.94731587130725443342927317226666E-16,
    +0.86934920504480078871023098666666E-19,
    +0.63469394403318040457397333333333E-22,
    +0.37740337870858485738666666666666E-25,
    +0.18630213719570056533333333333333E-28,
    +0.77568437166506666666666666666666E-32 };
  static double sqeps = 0.0;
  double value;
  double y;
  static double ymax = 0.0;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( sinhcs, 13, 0.1 * r8_mach ( 3 ) );
    sqeps = sqrt ( 6.0 * r8_mach ( 3 ) );
    ymax = 1.0 / sqrt ( r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= sqeps )
  {
    value = x;
  }
  else if ( y <= 1.0 )
  {
    value = x * ( 1.0 + r8_csevl ( 2.0 * x * x - 1.0, sinhcs, nterms ) );
  }
  else
  {
    y = exp ( y );

    if ( ymax <= y )
    {
      value = 0.5 * y;
    }
    else
    {
      value = 0.5 * ( y - 1.0 / y );
    }

    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

double r8_spence ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SPENCE evaluates a form of Spence's function for an R8 argument.
//
//  Discussion:
//
//    This function evaluates a form of Spence's function defined by
//
//      f(x) = Integral ( 0 <= y <= x ) - log ( abs ( 1 - y ) ) / y dy
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions, page 1004,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//    K Mitchell,
//    Tables of the function Integral ( 0 < y < x ) - log | 1 - y | dy / y
//    with an account of some properties of this and related functions,
//    Philosophical Magazine,
//    Volume 40, pages 351-368, 1949.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_SPENCE, Spence's function evaluated at X.
//
{
  double aln;
  static int nspenc = 0;
  static double pi26 = +1.644934066848226436472415166646025189219;
  double spencs[38] = {
    +0.1527365598892405872946684910028,
    +0.8169658058051014403501838185271E-01,
    +0.5814157140778730872977350641182E-02,
    +0.5371619814541527542247889005319E-03,
    +0.5724704675185826233210603054782E-04,
    +0.6674546121649336343607835438589E-05,
    +0.8276467339715676981584391689011E-06,
    +0.1073315673030678951270005873354E-06,
    +0.1440077294303239402334590331513E-07,
    +0.1984442029965906367898877139608E-08,
    +0.2794005822163638720201994821615E-09,
    +0.4003991310883311823072580445908E-10,
    +0.5823462892044638471368135835757E-11,
    +0.8576708692638689278097914771224E-12,
    +0.1276862586280193045989483033433E-12,
    +0.1918826209042517081162380416062E-13,
    +0.2907319206977138177795799719673E-14,
    +0.4437112685276780462557473641745E-15,
    +0.6815727787414599527867359135607E-16,
    +0.1053017386015574429547019416644E-16,
    +0.1635389806752377100051821734570E-17,
    +0.2551852874940463932310901642581E-18,
    +0.3999020621999360112770470379519E-19,
    +0.6291501645216811876514149171199E-20,
    +0.9933827435675677643803887752533E-21,
    +0.1573679570749964816721763805866E-21,
    +0.2500595316849476129369270954666E-22,
    +0.3984740918383811139210663253333E-23,
    +0.6366473210082843892691326293333E-24,
    +0.1019674287239678367077061973333E-24,
    +0.1636881058913518841111074133333E-25,
    +0.2633310439417650117345279999999E-26,
    +0.4244811560123976817224362666666E-27,
    +0.6855411983680052916824746666666E-28,
    +0.1109122433438056434018986666666E-28,
    +0.1797431304999891457365333333333E-29,
    +0.2917505845976095173290666666666E-30,
    +0.4742646808928671061333333333333E-31 };
  double value;
  static double xbig = 0.0;

  if ( nspenc == 0 )
  {
    nspenc = r8_inits ( spencs, 38, 0.1 * r8_mach ( 3 ) );
    xbig = 1.0 / r8_mach ( 3 );
  }

  if ( x <= - xbig )
  {
    aln = log ( 1.0 - x );
    value = - pi26  - 0.5 * aln * ( 2.0 * log ( - x ) - aln );
  }
  else if ( x <= - 1.0 )
  {
    aln = log ( 1.0 - x );

    value = - pi26 - 0.5 * aln * ( 2.0 
      * log ( - x ) - aln ) + ( 1.0 + r8_csevl (
      4.0 / ( 1.0 - x ) - 1.0, spencs, nspenc ) ) / ( 1.0 - x );
  }
  else if ( x <= 0.0 )
  {
    value = - 0.5 * log ( 1.0 - x ) 
      * log ( 1.0 - x ) - x * ( 1.0 + r8_csevl ( 
      4.0 * x / ( x - 1.0 ) - 1.0, spencs, nspenc ) ) / ( x - 1.0 );
  }
  else if ( x <= 0.5 )
  {
    value = x * ( 1.0 + r8_csevl ( 4.0 * x - 1.0, spencs, nspenc ) );
  }
  else if ( x < 1.0 )
  {
    value = pi26 - log ( x ) * log ( 1.0 - x )
      - ( 1.0 - x ) * ( 1.0 + r8_csevl ( 4.0 
      * ( 1.0 - x ) - 1.0, spencs, nspenc ) );
  }
  else if ( x == 1.0 )
  {
    value = pi26;
  }
  else if ( x <= 2.0 )
  {
    value = pi26 - 0.5 * log ( x ) 
      * log ( ( x - 1.0 ) * ( x - 1.0 ) / x )
      + ( x - 1.0 ) * ( 1.0 + r8_csevl ( 4.0 
      * ( x - 1.0 ) / x - 1.0, spencs, nspenc ) ) / x;
  }
  else if ( x < xbig )
  {
    value = 2.0 * pi26 - 0.5 * log ( x ) * log ( x )
      - ( 1.0 + r8_csevl ( 4.0 / x - 1.0, spencs, nspenc ) ) / x;
  }
  else
  {
    value = 2.0 * pi26 - 0.5 * log ( x ) * log ( x );
  }
  return value;
}
//****************************************************************************80

double r8_sqrt ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SQRT computes the square root of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the number whose square root is desired.
//
//    Output, double R8_SQRT, the square root of X.
//
{
  int irem;
  int iter;
  int ixpnt;
  int n;
  static int niter = 0;
  static double sqrt2[3] = {
    0.70710678118654752440084436210485,
    1.0,
    1.41421356237309504880168872420970 };
  double value;
  double y;

  if ( niter == 0 )
  {
    niter = 1.443 * r8_log ( - 0.104 * r8_log ( 0.1 * r8_mach ( 3 ) ) ) + 1.0;
  }

  if ( x < 0.0 )
  {
    cerr << "\n";
    cerr << "R8_SQRT - Fatal error!\n";
    cerr << "  X is negative.\n";
    exit ( 1 );
  }
  else if ( x == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    r8_upak ( x, y, n );
    ixpnt = n / 2;
    irem = n - 2 * ixpnt + 2;
    value = 0.261599 + y * ( 1.114292 + y * ( -0.516888 + y * 0.141067 ) );

    for ( iter = 1; iter <= niter; iter++ )
    {
      value = value + 0.5 * ( y - value * value ) / value;
    }
    value = r8_pak ( sqrt2[irem-1] * value, ixpnt );
  }
  return value;
}
//****************************************************************************80

double r8_tan ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TAN evaluates the tangent of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_TAN, the tangent of X.
//
{
  double ainty;
  double ainty2;
  int ifn;
  static int nterms = 0;
  static double pi2rec = 0.011619772367581343075535053490057;
  double prodbg;
  static double sqeps = 0.0;
  static double tancs[19] = {
    +0.22627932763129357846578636531752,
    +0.43017913146548961775583410748067E-01,
    +0.68544610682565088756929473623461E-03,
    +0.11045326947597098383578849369696E-04,
    +0.17817477903926312943238512588940E-06,
    +0.28744968582365265947529646832471E-08,
    +0.46374854195902995494137478234363E-10,
    +0.74817609041556138502341633308215E-12,
    +0.12070497002957544801644516947824E-13,
    +0.19473610812823019305513858584533E-15,
    +0.31417224874732446504614586026666E-17,
    +0.50686132555800153941904891733333E-19,
    +0.81773105159836540043979946666666E-21,
    +0.13192643412147384408951466666666E-22,
    +0.21283995497042377309866666666666E-24,
    +0.34337960192345945292800000000000E-26,
    +0.55398222121173811200000000000000E-28,
    +0.89375227794352810666666666666666E-30,
    +0.14419111371369130666666666666666E-31 };
  double value;
  static double xmax = 0.0;
  static double xsml = 0.0;
  double y;
  double yrem;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( tancs, 19, 0.1 * r8_mach ( 3 ) );
    xmax = 1.0 / r8_mach ( 4 );
    xsml = sqrt ( 3.0 * r8_mach ( 3 ) );
    sqeps = sqrt ( r8_mach ( 4 ) );
  }

  y = r8_abs ( x );

  if ( xmax < y )
  {
    cerr << "\n";
    cerr << "R8_TAN - Warning\n";
    cerr << "  No precision because |X| is big.\n";
    value = 0.0;
    return value;
  }
//
//  Carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
//  = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
//  = aint(.625*y) + aint(z) + rem(z)
//
  ainty = r8_aint ( y );
  yrem = y - ainty;
  prodbg = 0.625 * ainty;
  ainty = r8_aint ( prodbg );
  y = ( prodbg - ainty ) + 0.625 * yrem + pi2rec * y;
  ainty2 = r8_aint ( y );
  ainty = ainty + ainty2;
  y = y - ainty2;

  ifn = ( int ) fmod ( ainty, 2.0 );

  if ( ifn == 1 )
  {
    y = 1.0 - y;
  }

  if ( 1.0 - y < r8_abs ( x ) * sqeps )
  {
    cerr << "\n";
    cerr << "R8_TAN - Warning!\n";
    cerr << "  Answer < half precision.\n";
    cerr << "  |X| big or X near pi/2 or 3*pi/2.\n";
  }

  if ( y == 1.0 )
  {
    cerr << "\n";
    cerr << "R8_TAN - Fatal error!\n";
    cerr << "  X is pi/2 or 3*pi/2.\n";
    exit ( 1 );
  }

  if ( y <= 0.25 )
  {
    value = y;
    if ( xsml < y )
    {
      value = y * ( 1.5 + r8_csevl ( 32.0 * y * y - 1.0, tancs, nterms ) );
    }
  }
  else if ( y <= 0.5 )
  {
    value = 0.5 * y * ( 1.5 + r8_csevl ( 
      8.0 * y * y - 1.0, tancs, nterms ) );
    value = 2.0 * value / ( 1.0 - value * value );
  }
  else
  {
    value = 0.25 * y * ( 1.5 + r8_csevl ( 
      2.0 * y * y - 1.0, tancs, nterms ) );
    value = 2.0 * value / ( 1.0 - value * value );
    value = 2.0 * value / ( 1.0 - value * value );
  }

  if ( x < 0.0 )
  {
    value = - r8_abs ( value );
  }
  else if ( 0.0 < x )
  {
    value = + r8_abs ( value );
  }

  if ( ifn == 1 )
  {
    value = - value;
  }
  return value;
}
//****************************************************************************80

double r8_tanh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TANH evaluates the hyperbolic tangent of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_TANH, the hyperbolic tangent of X.
//
{
  static int nterms = 0;
  static double sqeps = 0.0;
  static double tanhcs[31] = {
    -0.25828756643634710438338151450605,
    -0.11836106330053496535383671940204,
    +0.98694426480063988762827307999681E-02,
    -0.83579866234458257836163690398638E-03,
    +0.70904321198943582626778034363413E-04,
    -0.60164243181207040390743479001010E-05,
    +0.51052419080064402965136297723411E-06,
    -0.43320729077584087216545467387192E-07,
    +0.36759990553445306144930076233714E-08,
    -0.31192849612492011117215651480953E-09,
    +0.26468828199718962579377758445381E-10,
    -0.22460239307504140621870997006196E-11,
    +0.19058733768288196054319468396139E-12,
    -0.16172371446432292391330769279701E-13,
    +0.13723136142294289632897761289386E-14,
    -0.11644826870554194634439647293781E-15,
    +0.98812684971669738285540514338133E-17,
    -0.83847933677744865122269229055999E-18,
    +0.71149528869124351310723506176000E-19,
    -0.60374242229442045413288837119999E-20,
    +0.51230825877768084883404663466666E-21,
    -0.43472140157782110106047829333333E-22,
    +0.36888473639031328479423146666666E-23,
    -0.31301874774939399883325439999999E-24,
    +0.26561342006551994468488533333333E-25,
    -0.22538742304145029883494399999999E-26,
    +0.19125347827973995102208000000000E-27,
    -0.16228897096543663117653333333333E-28,
    +0.13771101229854738786986666666666E-29,
    -0.11685527840188950118399999999999E-30,
    +0.99158055384640389120000000000000E-32 };
  double value;
  static double xmax = 0.0;
  double y;
  double yrec;

  if ( nterms == 0 )
  {
    nterms = r8_inits ( tanhcs, 31, 0.1 * r8_mach ( 3 ) );
    sqeps = sqrt ( 3.0 * r8_mach ( 3 ) );
    xmax = - 0.5 * log ( r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= sqeps )
  {
    value = x;
  }
  else if ( y <= 1.0 )
  {
    value = x * ( 1.0 + r8_csevl ( 2.0 * x * x - 1.0, tanhcs, nterms ) );
  }
  else if ( y <= xmax )
  {
    y = exp ( y );
    yrec = 1.0 / y;
    value = ( y - yrec ) / ( y + yrec );

    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  else
  {
    if ( x < 0.0 )
    {
      value = - 1.0;
    }
    else
    {
      value = + 1.0;
    }
  }
  return value;
}
//****************************************************************************80

void r8_upak ( double x, double &y, int &n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UPAK unpacks an R8 into a mantissa and exponent.
//
//  Discussion:
//
//    This function unpacks a floating point number x so that
//
//      x = y * 2.0^n
//
//    where
//
//      0.5 <= abs ( y ) < 1.0 .
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the number to be unpacked.
//
//    Output, double &Y, the mantissa.
//
//    Output, int &N, the exponent.
//
{
  double absx;

  absx = r8_abs ( x );
  n = 0;
  y = 0.0;

  if ( x == 0.0 )
  {
    return;
  }

  while ( absx < 0.5 )
  {
    n = n - 1;
    absx = absx * 2.0;
  }

  while ( 1.0 <= absx )
  {
    n = n + 1;
    absx = absx * 0.5;
  }

  if ( x < 0.0 )
  {
    y = - absx;
  }
  else
  {
    y = + absx;
  }

  return;
}
