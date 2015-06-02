# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>
# include <string>

using namespace std;

# include "c4lib.hpp"
# include "r4lib.hpp"

//****************************************************************************80

float c4_abs ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ABS returns the absolute value of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_ABS, the magnitude of X.
//
{
  float value;

  value = sqrt ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) );

  return value;
}
//****************************************************************************80

complex <float> c4_acos ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ACOS evaluates the inverse cosine of a C4.
//
//  Discussion:
//
//    Here we use the relationship:
//
//      C4_ACOS ( Z ) = pi/2 - C4_ASIN ( Z ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_ACOS, the function value.
//
{
  complex <float> c2;
  float c2_imag;
  float c2_real;
  float r4_pi_half = 1.57079632679489661923;

  c2 = c4_asin ( c1 );

  c2_real = r4_pi_half - real ( c2 );
  c2_imag =            - imag ( c2 );

  c2 = complex <float> ( c2_real, c2_imag );

  return c2;
}
//****************************************************************************80

complex <float> c4_acosh ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ACOSH evaluates the inverse hyperbolic cosine of a C4.
//
//  Discussion:
//
//    Here we use the relationship:
//
//      C4_ACOSH ( Z ) = i * C4_ACOS ( Z ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_ACOSH, the function value.
//
{
  complex <float> c2;

  c2 = c4_i ( ) * c4_acos ( c1 );
  
  return c2;
}
//****************************************************************************80

complex <float> c4_add ( complex <float> c1, complex <float> c2 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ADD adds two C4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, C2, the arguments.
//
//    Output, complex <float> C4_ADD, the sum of C1 and C2.
//
{
  complex <float> c3;

  c3 = c1 + c2;

  return c3;
}
//****************************************************************************80

float c4_arg ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ARG returns the argument of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose argument is desired.
//
//    Output, float C4_ARG, the argument of X.
//
{
  float value;

  if ( imag ( x ) == 0.0 && real ( x ) == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = atan2 ( imag ( x ), real ( x ) );
  }

  return value;
}
//****************************************************************************80

complex <float> c4_asin ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ASIN evaluates the inverse sine of a C4.
//
//  Discussion:
//
//    Here we use the relationship:
//
//      C4_ASIN ( Z ) = - i * log ( i * z + sqrt ( 1 - z * z ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_ASIN, the function value.
//
{
  complex <float> c2;
  complex <float> c3;
  complex <float> c4;
  complex <float> c5;
  complex <float> ce;

  c2 = c4_i ( );
  c5 = c4_one ( ) - c1 * c1;
  c3 = c4_sqrt ( c5 );
  c4 = c4_log ( c3 + c2 * c1 );
  ce = - c2 * c4;

  return ce;
}
//****************************************************************************80

complex <float> c4_asinh ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ASINH evaluates the inverse hyperbolic sine of a C4.
//
//  Discussion:
//
//    Here we use the relationship:
//
//      C4_ASINH ( Z ) = - i * C4_ASIN ( i * Z ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_ASINH, the function value.
//
{
  complex <float> c2;
  complex <float> c3;
  complex <float> c4;
  complex <float> c5;
  complex <float> c6;

  c2 = c4_i ( );
  c3 = c2 * c1;
  c4 = c4_asin ( c3 );
  c5 = c2 * c4;
  c6 = - c5;

  return c6;
}
//****************************************************************************80

complex <float> c4_atan ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ATAN evaluates the inverse tangent of a C4.
//
//  Discussion:
//
//    Here we use the relationship:
//
//      C4_ATAN ( Z ) = ( i / 2 ) * log ( ( 1 - i * z ) / ( 1 + i * z ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_ATAN, the function value.
//
{
  complex <float> c2;
  complex <float> c3;
  complex <float> c4;
  complex <float> c5;
  complex <float> c6;
  complex <float> c7;
  complex <float> c8;
  complex <float> c9;
  complex <float> cx;

  c2 = c4_i ( );
  c3 = c4_one ( );
  c4 = c4_mul ( c2, c1 );
  c5 = c4_sub ( c3, c4 );
  c6 = c4_add ( c3, c4 );
  c7 = c4_div ( c5, c6 );

  c8 = c4_log ( c7 );
  c9 = c4_mul ( c2, c8 );
  cx = c9 / complex <float> ( 2.0, 0.0 );

  return cx;
}
//****************************************************************************80

complex <float> c4_atanh ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ATANH evaluates the inverse hyperbolic tangent of a C4.
//
//  Discussion:
//
//    Here we use the relationship:
//
//      C4_ATANH ( Z ) = - i * C4_ATAN ( i * Z ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_ATANH, the function value.
//
{
  complex <float> c2;
  complex <float> c3;
  complex <float> c4;
  complex <float> c5;
  complex <float> c6;

  c2 = c4_i ( );

  c3 = c4_mul ( c2, c1 );
  c4 = c4_atan ( c3 );
  c5 = c4_mul ( c2, c4 );
  c6 = c4_neg ( c5 );

  return c6;
}
//****************************************************************************80

complex <float> c4_conj ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_CONJ conjugates a C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_CONJ, the function value.
//
{
  complex <float> c2;

  c2 = conj ( c1 );

  return c2;
}
//****************************************************************************80

void c4_copy ( complex <float> c1, complex <float> c2 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_COPY copies a C4.
//
//  Discussion:
//
//    The order of the arguments may seem unnatural, but it is arranged so
//    that the call
//
//      c4_copy ( c1, c2 )
//
//    mimics the assignment
//
//      c1 = c2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, complex <float> C1, the copy of C2.
//
//    Input, complex <float> C2, the value to be copied.
//
{
  c1 = c2;

  return;
}
//****************************************************************************80

complex <float> c4_cos ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_COS evaluates the cosine of a C4.
//
//  Discussion:
//
//    We use the relationship:
//
//      C4_COS ( C ) = ( C4_EXP ( i * C ) + C4_EXP ( - i * C ) ) / 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_COS, the function value.
//
{
  complex <float> c2;

  c2 = ( exp ( c1 * c4_i ( ) ) + exp ( - c1 * c4_i ( ) ) ) 
    / complex <float> ( 2.0, 0.0 );

  return c2;
}
//****************************************************************************80

complex <float> c4_cosh ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_COSH evaluates the hyperbolic cosine of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_COSH, the function value.
//
{
  complex <float> c2;
  complex <float> c3;
  complex <float> c4;
  complex <float> c5;
  complex <float> c6;

  c2 = c4_exp ( c1 );

  c3 = c4_neg ( c1 );
  c4 = c4_exp ( c3 );

  c5 = c4_add ( c2, c4 );
  c6 = c4_div_r4 ( c5, 2.0 );

  return c6;
}
//****************************************************************************80

complex <float> c4_cube_root ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_CUBE_ROOT returns the principal cube root of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the number whose cube root is desired.
//
//    Output, complex <float> C4_CUBE_ROOT, the cube root of X.
//
{
  float argument;
  float magnitude;
  complex <float> value;

  argument = c4_arg ( x );
  magnitude = c4_mag ( x );

  if ( magnitude == 0.0 )
  {
    value = complex <float> ( 0.0, 0.0 );
  }
  else
  {
    value = pow ( magnitude, ( float ) ( 1.0 / 3.0 ) ) 
      * complex <float> ( cos ( argument / 3.0 ), sin ( argument / 3.0 ) );
  }

  return value;
}
//****************************************************************************80

complex <float> c4_div ( complex <float> c1, complex <float> c2 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_DIV divides two C4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, C2, the arguments.
//
//    Output, complex <float> C4_DIV, the function value.
//
{
  float c2_norm;
  complex <float> c3;
  float c3_imag;
  float c3_real;

  c2_norm = c4_abs ( c2 );

  c3_real = ( real ( c1 ) * real ( c2 ) 
            + imag ( c1 ) * imag ( c2 ) ) / c2_norm / c2_norm;

  c3_imag = ( imag ( c1 ) * real ( c2 ) 
            - real ( c1 ) * imag ( c2 ) ) / c2_norm / c2_norm;

  c3 = complex <float> ( c3_real, c3_imag );

  return c3;
}
//****************************************************************************80

complex <float> c4_div_r4 ( complex <float> c1, float r )

//****************************************************************************80
//
//  Purpose:
//
//    C4_DIV_R4 divides a C4 by an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the value to be divided.
//
//    Input, float R, the divisor.
//
//    Output, complex <float> C4_DIV_R4, the function value.
//
{
  complex <float> c2;

  c2 = c1 / r;

  return c2;
}
//****************************************************************************80

complex <float> c4_exp ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_EXP exponentiates a C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_EXP, the function value.
//
{
  complex <float> c2;

  c2 = exp ( c1 );

  return c2;
}
//****************************************************************************80

complex <float> c4_i ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4_I returns the value of the imaginary unit, i as a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, complex <float> C4_I, the value of complex i.
//
{
  complex <float> value;

  value = complex <float> ( 0.0, 1.0 );

  return value;
}
//****************************************************************************80

float c4_imag ( complex <float> c )

//****************************************************************************80
//
//  Purpose:
//
//    C4_IMAG returns the imaginary part of a C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C, the argument.
//
//    Output, float C4_IMAG, the function value.
//
{
  float value;

  value = imag ( c );

  return value;
}
//****************************************************************************80

complex <float> c4_inv ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_INV inverts a C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_INV, the function value;
//
{
  complex <float> c2;

  c2 = c4_one ( ) / c1;

  return c2;
}
//****************************************************************************80

bool c4_le_l1 ( complex <float> x, complex <float> y )

//****************************************************************************80
//
//  Purpose:
//
//    C4_LE_L1 := X <= Y for C4 values, and the L1 norm.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    The L1 norm can be defined here as:
//
//      C4_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, Y, the values to be compared.
//
//    Output, bool C4_LE_L1, is TRUE if X <= Y.
//
{
  bool value;

  if ( r4_abs ( real ( x ) ) + r4_abs ( imag ( x ) ) <= 
       r4_abs ( real ( y ) ) + r4_abs ( imag ( y ) ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

bool c4_le_l2 ( complex <float> x, complex <float> y )

//****************************************************************************80
//
//  Purpose:
//
//    C4_LE_L2 := X <= Y for C4 values, and the L2 norm.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    The L2 norm can be defined here as:
//
//      C4_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag (X) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, Y, the values to be compared.
//
//    Output, bool C4_LE_L2, is TRUE if X <= Y.
//
{
  bool value;

  if ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) <= 
       pow ( real ( y ), 2 ) + pow ( imag ( y ), 2 ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

bool c4_le_li ( complex <float> x, complex <float> y )

//****************************************************************************80
//
//  Purpose:
//
//    C4_LE_LI := X <= Y for C4 values, and the L-oo norm.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    The L-oo norm can be defined here as:
//
//      C4_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, Y, the values to be compared.
//
//    Output, bool C4_LE_LI, is TRUE if X <= Y.
//
{
  bool value;

  if ( r4_max ( r4_abs ( real ( x ) ), r4_abs ( imag ( x ) ) ) <= 
       r4_max ( r4_abs ( real ( y ) ), r4_abs ( imag ( y ) ) ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

complex <float> c4_log ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_LOG evaluates the logarithm of a C4.
//
//  Discussion:
//
//    Here we use the relationship:
//
//      C4_LOG ( Z ) = LOG ( MAG ( Z ) ) + i * ARG ( Z )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_LOG, the function value.
//
{
  float arg;
  complex <float> c2;
  float mag;

  arg = c4_arg ( c1 );
  mag = c4_mag ( c1 );

  c2 = complex <float> ( log ( mag ), arg );

  return c2;
}
//****************************************************************************80

float c4_mag ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_MAG returns the magnitude of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_MAG, the magnitude of X.
//
{
  float value;

  value = sqrt ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) );

  return value;
}
//****************************************************************************80

complex <float> c4_mul ( complex <float> c1, complex <float> c2 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_MUL multiplies two C4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, C2, the arguments.
//
//    Output, complex <float> C4_MUL, the function value.
//
{
  complex <float> c3;

  c3 = c1 * c2;

  return c3;
}
//****************************************************************************80

complex <float> c4_neg ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NEG negates a C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_NEG, the function value.
//
{
  complex <float> c2;

  c2 = - c1;

  return c2;
}
//****************************************************************************80

complex <float> c4_nint ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NINT returns the nearest complex integer of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the value to be NINT'ed.
//
//    Output, complex <float> C4_NINT, the NINT'ed value.
//
{
  float r;
  float r_min;
  float x;
  float x_min;
  float xc;
  float y;
  float y_min;
  float yc;
  complex <float> value;

  xc = real ( c1 );
  yc = imag ( c1 );
//
//  Lower left.
//
  x = r4_floor ( real ( c1 ) );
  y = r4_floor ( imag ( c1 ) );
  r = pow ( x - xc, 2 ) + pow ( y - yc, 2 );
  r_min = r;
  x_min = x;
  y_min = y;
//
//  Lower right.
//
  x = r4_floor ( real ( c1 ) ) + 1.0;
  y = r4_floor ( imag ( c1 ) );
  r = pow ( x - xc, 2 ) + pow ( y - yc, 2 );
  if ( r < r_min )
  {
    r_min = r;
    x_min = x;
    y_min = y;
  }
//
//  Upper right.
//
  x = r4_floor ( real ( c1 ) ) + 1.0;
  y = r4_floor ( imag ( c1 ) ) + 1.0;
  r = pow ( x - xc, 2 ) + pow ( y - yc, 2 );
  if ( r < r_min )
  {
    r_min = r;
    x_min = x;
    y_min = y;
  }
//
//  Upper left.
//
  x = r4_floor ( real ( c1 ) );
  y = r4_floor ( imag ( c1 ) ) + 1.0;
  r = pow ( x - xc, 2 ) + pow ( y - yc, 2 );
  if ( r < r_min )
  {
    r_min = r;
    x_min = x;
    y_min = y;
  }

  value = complex <float> ( x_min, y_min );

  return value;
}
//****************************************************************************80

float c4_norm_l1 ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NORM_L1 evaluates the L1 norm of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    Numbers of equal norm lie along diamonds centered at (0,0).
//
//    The L1 norm can be defined here as:
//
//      C4_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_NORM_L1, the norm of X.
//
{
  float value;

  value = r4_abs ( real ( x ) ) + r4_abs ( imag ( x ) );

  return value;
}
//****************************************************************************80

float c4_norm_l2 ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NORM_L2 evaluates the L2 norm of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    Numbers of equal norm lie on circles centered at (0,0).
//
//    The L2 norm can be defined here as:
//
//      C4_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag ( X ) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_NORM_L2, the 2-norm of X.
//
{
  float value;

  value = sqrt ( pow ( real ( x ), 2 )
               + pow ( imag ( x ), 2 ) );

  return value;
}
//****************************************************************************80

float c4_norm_li ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NORM_LI evaluates the L-oo norm of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//    Numbers of equal norm lie along squares whose centers are at (0,0).
//
//    The L-oo norm can be defined here as:
//
//      C4_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the value whose norm is desired.
//
//    Output, float C4_NORM_LI, the L-oo norm of X.
//
{
  float value;

  value = r4_max ( r4_abs ( real ( x ) ), r4_abs ( imag ( x ) ) );

  return value;
}
//****************************************************************************80

complex <float> c4_normal_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4_NORMAL_01 returns a unit pseudonormal C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, complex <float> C4_NORMAL_01, a unit pseudonormal value.
//
{
  float pi = 3.141592653589793;
  float v1;
  float v2;
  complex <float> value;
  float x_c;
  float x_r;

  v1 = r4_uniform_01 ( seed );
  v2 = r4_uniform_01 ( seed );

  x_r = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * pi * v2 );
  x_c = sqrt ( - 2.0 * log ( v1 ) ) * sin ( 2.0 * pi * v2 );

  value = complex <float> ( x_r, x_c );

  return value;
}
//****************************************************************************80

complex <float> c4_one ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ONE returns the value of complex 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, complex <float> C4_ONE, the value of complex 1.
//
{
  complex <float> value;

  value = complex <float> ( 1.0, 0.0);

  return value;
}
//****************************************************************************80

void c4_print ( complex <float> a, string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4_PRINT prints a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> A, the value to be printed.
//
//    Input, string TITLE, a title.
//
{
  cout << title
       << "  ( " << setw(14) << real ( a )
       << ", "   << setw(14) << imag ( a ) << " )\n";

  return;
}
//****************************************************************************80

float c4_real ( complex <float> c )

//****************************************************************************80
//
//  Purpose:
//
//    C4_REAL returns the real part of a C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C, the complex number.
//
//    Output, float C4_REAL, the function value.
//
{
  float value;

  value = real ( c );

  return value;
}
//****************************************************************************80

complex <float> c4_sin ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_SIN evaluates the sine of a C4.
//
//  Discussion:
//
//    We use the relationship:
//
//      C4_SIN ( C ) = - i * ( C4_EXP ( i * C ) - C4_EXP ( - i * C ) ) / 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_SIN, the function value.
//
{
  complex <float> c2;
  complex <float> c3;
  complex <float> c4;
  complex <float> c5;
  complex <float> c6;
  complex <float> c7;
  complex <float> c8;
  complex <float> c9;
  complex <float> cx;
  float r;

  c2 = c4_i ( );

  c3 = c4_mul ( c2, c1 );
  c4 = c4_exp ( c3 );

  c5 = c4_neg ( c3 );
  c6 = c4_exp ( c5 );

  c7 = c4_sub ( c4, c6 );

  r = 2.0;
  c8 = c4_div_r4 ( c7, r );
  c9 = c4_mul ( c8, c2 );
  cx = c4_neg ( c9 );

  return cx;
}
//****************************************************************************80

complex <float> c4_sinh ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_SINH evaluates the hyperbolic sine of a C4.
//
//  Discussion:
//
//    We use the relationship:
//
//      C4_SINH ( C ) = ( C4_EXP ( C ) - C4_EXP ( - i * C ) ) / 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_SINH, the function value.
//
{
  complex <float> c2;
  complex <float> c3;
  complex <float> c4;
  complex <float> c5;
  complex <float> c6;
  float r;

  c2 = c4_exp ( c1 );

  c3 = c4_neg ( c1 );
  c4 = c4_exp ( c3 );

  c5 = c4_sub ( c2, c4 );

  r = 2.0;
  c6 = c4_div_r4 ( c5, r );

  return c6;
}
//****************************************************************************80

complex <float> c4_sqrt ( complex <float> x )

//****************************************************************************80
//
//  Purpose:
//
//    C4_SQRT returns the principal square root of a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> X, the number whose square root is desired.
//
//    Output, complex <float> C4_SQRT, the square root of X.
//
{
  float argument;
  float magnitude;
  complex <float> value;

  argument = c4_arg ( x );
  magnitude = c4_mag ( x );

  if ( magnitude == 0.0 )
  {
    value = complex <float> ( 0.0, 0.0 );
  }
  else
  {
    value = sqrt ( magnitude ) 
      * complex <float> ( cos ( argument / 2.0 ), sin ( argument / 2.0 ) );
  }

  return value;
}
//****************************************************************************80

complex <float> c4_sub ( complex <float> c1, complex <float> c2 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_SUB subtracts two C4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, C2, the arguments.
//
//    Output, complex <float> C4_SUB, the function value.
//
{
  complex <float> c3;

  c3 = c1 - c2;

  return c3;
}
//****************************************************************************80

void c4_swap ( complex <float> *x, complex <float> *y )

//****************************************************************************80
//
//  Purpose:
//
//    C4_SWAP swaps two C4's.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, complex <float> *X, *Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  complex <float> z;

   z = *x; 
  *x = *y;
  *y =  z;

  return;
}
//****************************************************************************80

complex <float> c4_tan ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_TAN evaluates the tangent of a C4.
//
//  Discussion:
//
//    We use the relationship:
//
//      C4_TAN ( C ) = - i * ( C4_EXP ( i * C ) - C4_EXP ( - i * C ) ) 
//                         / ( C4_EXP ( I * C ) + C4_EXP ( - i * C ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_TAN, the function value.
//
{
  complex <float> c2;
  complex <float> c3;
  complex <float> c4;
  complex <float> c5;
  complex <float> c6;
  complex <float> c7;
  complex <float> c8;
  complex <float> c9;
  complex <float> cx;
  complex <float> ce;

  c2 = c4_i ( );
  c3 = c4_mul ( c2, c1 );
  c4 = c4_neg ( c3 );
  
  c5 = c4_exp ( c3 );
  c6 = c4_exp ( c4 );

  c7 = c4_sub ( c5, c6 );
  c8 = c4_add ( c5, c6 );

  c9 = c4_div ( c7, c8 );
  cx = c4_mul ( c2, c9 );
  ce = c4_neg ( cx );

  return ce;
}
//****************************************************************************80

complex <float> c4_tanh ( complex <float> c1 )

//****************************************************************************80
//
//  Purpose:
//
//    C4_TANH evaluates the hyperbolic tangent of a C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C1, the argument.
//
//    Output, complex <float> C4_TANH, the function value.
//
{
  complex <float> c2;
  complex <float> c3;
  complex <float> c4;
  complex <float> c5;
  complex <float> c6;
  complex <float> c7;

  c2 = c4_exp ( c1 );

  c3 = c4_neg ( c1 );
  c4 = c4_exp ( c3 );

  c5 = c4_sub ( c2, c4 );
  c6 = c4_add ( c2, c4 );

  c7 = c4_div ( c5, c6 );

  return c7;
}
//****************************************************************************80

void c4_to_cartesian ( complex <float> c, float *x, float *y )

//****************************************************************************80
//
//  Purpose:
//
//    C4_TO_CARTESIAN converts a C4 to Cartesian form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C, the argument.
//
//    Output, float *X, *Y, the Cartesian form.
//
{
  *x = real ( c );
  *y = imag ( c );

  return;
}
//****************************************************************************80

void c4_to_polar ( complex <float> c, float *r, float *theta )

//****************************************************************************80
//
//  Purpose:
//
//    C4_TO_POLAR converts a C4 to polar form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> C, the argument.
//
//    Output, float *R, *THETA, the polar form.
//
{
  *r = c4_abs ( c );
  *theta = c4_arg ( c );

  return;
}
//****************************************************************************80

complex <float> c4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4_UNIFORM_01 returns a unit pseudorandom C4.
//
//  Discussion:
//
//    The angle should be uniformly distributed between 0 and 2 * PI,
//    the square root of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4_UNIFORM_01, a pseudorandom complex value.
//
{
  int i4_huge = 2147483647;
  int k;
  float pi = 3.1415926;
  float r;
  float theta;
  complex <float> value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  theta = 2.0 * pi *
    ( ( float ) ( *seed ) * 4.656612875E-10 );

  value = complex <float> ( r * cos ( theta ), r * sin ( theta ) );

  return value;
}
//****************************************************************************80

complex <float> c4_zero ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4_ZERO returns the value of 0 as a C4.
//
//  Discussion:
//
//    A C4 is a complex <float> value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, complex <float> C4_ZERO, the value of complex 0.
//
{
  complex <float> value;

  value = complex <float> ( 0.0, 0.0 );

  return value;
}
//****************************************************************************80

void c4mat_add ( int m, int n, complex <float> alpha, complex <float> a[],
  complex <float> beta, complex <float> b[], complex <float> c[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_ADD combines two C4MAT's using complex scalar factors.
//
//  Discussion:
//
//    An C4MAT is a doubly dimensioned array of complex single precision values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, complex <float> ALPHA, the first scale factor.
//
//    Input, complex <float> A[M*N], the first matrix.
//
//    Input, complex <float> BETA, the second scale factor.
//
//    Input, complex <float> B[M*N], the second matrix.
//
//    Output, complex <float> C[M*N], the result.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
    }
  }
  return;
}
//****************************************************************************80

void c4mat_add_r4 ( int m, int n, float alpha, complex <float> a[],
  float beta, complex <float> b[], complex <float> c[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_ADD_R4 combines two C4MAT's using real scalar factors.
//
//  Discussion:
//
//    An C4MAT is a doubly dimensioned array of complex float precision values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float ALPHA, the first scale factor.
//
//    Input, complex <float> A[M*N], the first matrix.
//
//    Input, float BETA, the second scale factor.
//
//    Input, complex <float> B[M*N], the second matrix.
//
//    Output, complex <float> C[M*N], the result.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
    }
  }
  return;
}
//****************************************************************************80

void c4mat_copy ( int m, int n, complex <float> a1[], complex <float> a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_COPY copies one C4MAT to another.
//
//  Discussion:
//
//    An C4MAT is a doubly dimensioned array of complex <float> values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, complex <float> A1[M*N], the matrix to be copied.
//
//    Output, complex <float> A2[M*N], the copy of A1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return;
}
//****************************************************************************80

complex <float> *c4mat_copy_new ( int m, int n, complex <float> a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_COPY_NEW copies one C4MAT to a "new" C4MAT.
//
//  Discussion:
//
//    An C4MAT is a doubly dimensioned array of complex <float> values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, complex <float> A1[M*N], the matrix to be copied.
//
//    Output, complex <float> C4MAT_COPY_NEW[M*N], the copy of A1.
//
{
  complex <float> *a2;
  int i;
  int j;

  a2 = new complex <float>[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
}
//****************************************************************************80

void c4mat_fss ( int n, complex <float> a[], int nb, complex <float> x[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_FSS factors and solves a system with multiple right hand sides.
//
//  Discussion:
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, complex <float> A[N*N].
//    On input, A is the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input, int NB, the number of right hand sides.
//
//    Input/output, complex <float> X[N*NB], on input, the right hand sides of the
//    linear systems.  On output, the solutions of the linear systems.
//
{
  int i;
  int ipiv;
  int j;
  int jcol;
  float piv;
  complex <float> t;

  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = c4_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < c4_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = c4_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cout << "\n";
      cout << "C4MAT_FSS - Fatal error!\n";
      cout << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
//
//  Scale the pivot row.
//
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != c4_zero ( ) )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return;
}
//****************************************************************************80

complex <float> *c4mat_fss_new ( int n, complex <float> a[], int nb, 
  complex <float> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_FSS_NEW factors and solves a system with multiple right hand sides.
//
//  Discussion:
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, complex <float> A[N*N].
//    On input, A is the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input, int NB, the number of right hand sides.
//
//    Input, complex <float> B[N*NB], the right hand sides of the linear systems.
//
//    Output, complex <float> C4MAT_FSS_NEW[N*NB], the solutions of the 
//    linear systems.
//
{
  int i;
  int ipiv;
  int j;
  int jcol;
  float piv;
  complex <float> t;
  complex <float> *x;

  x = new complex <float>[n*nb];

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }
  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = c4_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol + 1; i <= n; i++ )
    {
      if ( piv < c4_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = c4_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "C4MAT_FSS_NEW - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
//
//  Scale the pivot row.
//
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol + 1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != c4_zero ( ) )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return x;
}
//****************************************************************************80

complex <float> *c4mat_identity_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_IDENTITY_NEW sets a C4MAT to the identity.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <float> C4MAT_IDENTITY_NEW[N*N], the matrix.
//
{
  complex <float> *a;
  int i;
  int j;

  a = new complex <float> [n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = complex <float> ( 1.0, 0.0 );
      }
      else
      {
        a[i+j*n] = complex <float> ( 0.0, 0.0 );
      }
    }
  }
  return a;
}
//****************************************************************************80

complex <float> *c4mat_indicator_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_INDICATOR_NEW returns the C4MAT indicator matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, complex <float> C4MAT_INDICATOR_NEW[M*N], the matrix.
//
{
  complex <float> *a;
  int i;
  int j;

  a = new complex <float> [m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = complex <float> ( i, j );
    }
  }
  return a;
}
//****************************************************************************80

void c4mat_minvm ( int n1, int n2, complex <float> a[], 
  complex <float> b[], complex <float> c[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_MINVM returns inverse(A) * B for C4MAT's.
//
//  Discussion:
//
//    A C4MAT is an array of C4 values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the matrices.
//
//    Input, complex <float> A[N1*N1], B[N1*N2], the matrices.
//
//    Output, complex <float> C[N1*N2], the result, 
//    C = inverse(A) * B.
//
{
  complex <float> *alu;

  alu = c4mat_copy_new ( n1, n1, a );

  c4mat_copy ( n1, n2, b, c );

  c4mat_fss ( n1, alu, n2, c );
 
  delete [] alu;

  return;
}
//****************************************************************************80

complex <float> *c4mat_minvm_new ( int n1, int n2, complex <float> a[], 
  complex <float> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_MINVM_NEW returns inverse(A) * B for C4MAT's.
//
//  Discussion:
//
//    A C4MAT is an array of C4 values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the order of the matrices.
//
//    Input, complex <float> A[N1*N1], B[N1*N2], the matrices.
//
//    Output, complex <float> C4MAT_MINVM_NEW[N1*N2], the result, 
//    C = inverse(A) * B.
//
{
  complex <float> *alu;
  complex <float> *c;

  alu = c4mat_copy_new ( n1, n1, a );
  c = c4mat_fss_new ( n1, alu, n2, b );
 
  delete [] alu;

  return c;
}
//****************************************************************************80

void c4mat_mm ( int n1, int n2, int n3, complex <float> a[], 
  complex <float> b[], complex <float> c[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_MM multiplies two matrices.
//
//  Discussion:
//
//    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, complex <float> A[N1*N2], complex <float> B[N2*N3], 
//    the matrices to multiply.
//
//    Output, complex <float> C[N1*N3], the product matrix C = A * B.
//
{
  complex <float> *c1;
  int i;
  int j;
  int k;

  c1 = new complex <float> [n1*n3];

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c1[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c1[i+j*n1] = c1[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  c4mat_copy ( n1, n3, c1, c );

  delete [] c1;

  return;
}
//****************************************************************************80

complex <float> *c4mat_mm_new ( int n1, int n2, int n3, complex <float> a[], 
  complex <float> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    A C4MAT is a doubly dimensioned array of C4 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, complex <float> A[N1*N2], complex <float> B[N2*N3], 
//    the matrices to multiply.
//
//    Output, complex <float> C4MAT_MM_NEW[N1*N3], the product matrix C = A * B.
//
{
  complex <float> *c;
  int i;
  int j;
  int k;

  c = new complex <float> [n1*n3];

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }
  return c;
}
//****************************************************************************80

void c4mat_nint ( int m, int n, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_NINT rounds the entries of a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of A.
//
//    Input/output, complex <float> A[M*N], the matrix to be NINT'ed.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = complex <float> ( 
        r4_nint ( real ( a[i+j*m] ) ), 
        r4_nint ( imag ( a[i+j*m] ) ) );
    }
  }
  return;
}
//****************************************************************************80

float c4mat_norm_fro ( int m, int n, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_NORM_FRO returns the Frobenius norm of a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of C4 values.
//
//    The Frobenius norm is defined as
//
//      C4MAT_NORM_FRO = sqrt (
//        sum ( 1 <= I <= M ) Sum ( 1 <= J <= N ) |A(I,J)| )
//
//    The matrix Frobenius-norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      c4vec_norm_l2 ( A*x ) <= c4mat_norm_fro ( A ) * c4vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the order of the matrix.
//
//    Input, complex <float> A[M*N], the matrix.
//
//    Output, float C4MAT_NORM_FRO, the Frobenius norm of A.
//
{
  int i;
  int j;
  float value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( real ( a[i+j*m] ), 2 )
                    + pow ( imag ( a[i+j*m] ), 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

float c4mat_norm_l1 ( int m, int n, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_NORM_L1 returns the matrix L1 norm of a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an MxN array of C4's, stored by (I,J) -> [I+J*M].
//
//    The matrix L1 norm is defined as:
//
//      C4MAT_NORM_L1 = max ( 1 <= J <= N )
//        sum ( 1 <= I <= M ) abs ( A(I,J) ).
//
//    The matrix L1 norm is derived from the vector L1 norm, and
//    satisifies:
//
//      c4vec_norm_l1 ( A * x ) <= c4mat_norm_l1 ( A ) * c4vec_norm_l1 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, complex <float> A[M*N], the matrix whose L1 norm is desired.
//
//    Output, float C4MAT_NORM_L1, the L1 norm of A.
//
{
  float col_sum;
  int i;
  int j;
  float value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    col_sum = 0.0;
    for ( i = 0; i < m; i++ )
    {
      col_sum = col_sum + c4_abs ( a[i+j*m] );
    }
    value = r4_max ( value, col_sum );
  }

  return value;
}
//****************************************************************************80

float c4mat_norm_li ( int m, int n, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_NORM_LI returns the matrix L-oo norm of a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of C4 values.
//
//    The matrix L-oo norm is defined as:
//
//      C4MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
//
//    The matrix L-oo norm is derived from the vector L-oo norm,
//    and satisifies:
//
//      c4vec_norm_li ( A * x ) <= c4mat_norm_li ( A ) * c4vec_norm_li ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, complex <float> A[M*N], the matrix whose L-oo
//    norm is desired.
//
//    Output, float C4MAT_NORM_LI, the L-oo norm of A.
//
{
  int i;
  int j;
  float row_sum;
  float value;

  value = 0.0;

  for ( i = 0; i < m; i++ )
  {
    row_sum = 0.0;
    for ( j = 0; j < n; j++ )
    {
      row_sum = row_sum + c4_abs ( a[i+j*m] );
    }
    value = r4_max ( value, row_sum );
  }
  return value;
}
//****************************************************************************80

void c4mat_print ( int m, int n, complex <float> a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_PRINT prints a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <float> A[M*N], the matrix.
//
//    Input, string TITLE, a title.
//
{
  c4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void c4mat_print_some ( int m, int n, complex <float> a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_PRINT_SOME prints some of a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <float> A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
  complex <float> c;
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int incx = 4;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( j2lo = jlo; j2lo <= i4_min ( jhi, n ); j2lo = j2lo + incx )
  {
    j2hi = j2lo + incx - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      cout << "     " << setw(10) << j << "     ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = ilo;
    if ( i2lo < 1 )
    {
      i2lo = 1;
    }
    i2hi = ihi;
    if ( m < i2hi )
    {
      i2hi = m;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(5) << i << ":";
//
//  Print out (up to) INCX entries in row I, that lie in the current strip.
//
      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;
        c = a[i-1+(j-1)*m];
        cout << "  " << setw(8) << real ( c )
             << "  " << setw(8) << imag ( c );
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

void c4mat_scale ( int m, int n, complex <float> alpha, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_SCALE scales a C4MAT by a complex scalar factor.
//
//  Discussion:
//
//    An C4MAT is a doubly dimensioned array of complex float precision values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, complex <float> ALPHA, the scale factor.
//
//    Input/output, complex <float> A[M*N], the matrix to be scaled.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] * alpha;
    }
  }
  return;
}
//****************************************************************************80

void c4mat_scale_r4 ( int m, int n, float alpha, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_SCALE_R4 scales a C4MAT by a real scalar factor.
//
//  Discussion:
//
//    An C4MAT is a doubly dimensioned array of complex float precision values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, float ALPHA, the scale factor.
//
//    Input/output, complex <float> A[M*N], the matrix to be scaled.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] * alpha;
    }
  }
  return;
}
//****************************************************************************80

void c4mat_uniform_01 ( int m, int n, int *seed, complex <float> c[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C[M*N], the pseudorandom complex matrix.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  float r;
  int k;
  float pi = 3.1415926;
  float theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * complex <float> ( cos ( theta ), sin ( theta ) );
    }
  }

  return;
}
//****************************************************************************80

complex <float> *c4mat_uniform_01_new ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_UNIFORM_01_NEW returns a new unit pseudorandom C4MAT.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4MAT_UNIFORM_01[M*N], the pseudorandom complex matrix.
//
{
  complex <float> *c;
  int i;
  int i4_huge = 2147483647;
  int j;
  float r;
  int k;
  float pi = 3.1415926;
  float theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4MAT_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  c = new complex <float>[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

      c[i+j*m] = r * complex <float> ( cos ( theta ), sin ( theta ) );
    }
  }

  return c;
}
//****************************************************************************80

complex <float> *c4mat_zero_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_ZERO_NEW returns a new zeroed C4MAT.
//
//  Discussion:
//
//    An C4MAT is a doubly dimensioned array of complex float precision values, 
//    which may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, complex <float> C4MAT_ZERO_NEW[M*N], the zeroed matrix.
//
{
  complex <float> *a;
  int i;
  int j;

  a = new complex <float>[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

void c4vec_copy ( int n, complex <float> a1[], complex <float> a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_COPY copies a C4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, complex <float> A1[N], the vector to be copied.
//
//    Output, complex <float> A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

complex <float> *c4vec_copy_new ( int n, complex <float> a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_COPY_NEW copies a C4VEC to a "new" C4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, complex <float> A1[N], the vector to be copied.
//
//    Output, complex <float> C4VEC_COPY_NEW[N], the copy of A1.
//
{
  complex <float> *a2;
  int i;

  a2 = new complex <float>[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

complex <float> *c4vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_INDICATOR_NEW sets a C4VEC to the indicator vector.
//
//  Discussion:
//
//    A C4VEC is a vector of complex <float> values.
//
//    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, complex <float> C4VEC_INDICATOR_NEW[N], the array.
//
{
  complex <float> *a;
  int i;

  a = new complex <float> [n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = complex <float> ( i+1, -i-1 );
  }

  return a;
}
//****************************************************************************80

void c4vec_nint ( int n, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_NINT rounds the entries of a C4VEC.
//
//  Discussion:
//
//    A C4VEC is a vector of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input/output, complex <float> A[N], the vector to be nint'ed.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = c4_nint ( a[i] );
  }

  return;
}
//****************************************************************************80

float c4vec_norm_l2 ( int n, complex <float> a[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_NORM_L2 returns the L2 norm of a C4VEC.
//
//  Discussion:
//
//    The vector L2 norm is defined as:
//
//      value = sqrt ( sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, complex <float> A[N], the vector whose L2 norm is desired.
//
//    Output, float C4VEC_NORM_L2, the L2 norm of A.
//
{
  int i;
  float value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value 
          + real ( a[i] ) * real ( a[i] ) 
          + imag ( a[i] ) * imag ( a[i] );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

void c4vec_print ( int n, complex <float> a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_PRINT prints a C4VEC.
//
//  Discussion:
//
//    A C4VEC is a vector of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, complex <float> A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";

  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8) << i
         << ": " << real ( a[i] )
         << "  " << imag ( a[i] ) << "\n";
  }

  return;
}
//****************************************************************************80

void c4vec_print_part ( int n, complex <float> a[], int max_print, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_PRINT_PART prints "part" of a C4VEC.
//
//  Discussion:
//
//    The user specifies MAX_PRINT, the maximum number of lines to print.
//
//    If N, the size of the vector, is no more than MAX_PRINT, then
//    the entire vector is printed, one entry per line.
//
//    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
//    followed by a line of periods suggesting an omission,
//    and the last entry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, complex <float> A[N], the vector to be printed.
//
//    Input, int MAX_PRINT, the maximum number of lines
//    to print.
//
//    Input, string TITLE, a title.
//
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << real ( a[i] ) 
           << "  " << setw(14) << imag ( a[i] ) << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << real ( a[i] ) 
           << "  " << setw(14) << imag ( a[i] ) << "\n";
    }
    cout << "  ........  ..............  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
           << "  " << setw(14) << real ( a[i] ) 
           << "  " << setw(14) << imag ( a[i] ) << "\n";
  }
  else
  {
    for ( i= 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << real ( a[i] ) 
           << "  " << setw(14) << imag ( a[i] ) << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << "  " << setw(14) << real ( a[i] ) 
         << "  " << setw(14) << imag ( a[i] )
         << "  " << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

void c4vec_print_some ( int n, complex <float> a[], int i_lo, int i_hi, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_PRINT_SOME prints some of a C4VEC.
//
//  Discussion:
//
//    A C4VEC is a vector of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, complex <float> A[N], the vector to be printed.
//
//    Input, int I_LO, I_HI, the first and last entries to print.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = i4_max ( 0, i_lo ); i <= i4_min ( i_hi, n - 1 ); i++ )
  {
    cout << "  " << setw(6) << i
         << ": " << real ( a[i] )
         << "  " << imag ( a[i] ) << "\n";
  }

  return;
}
//****************************************************************************80

void c4vec_sort_a_l2 ( int n, complex <float> x[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_SORT_A_L2 ascending sorts a C4VEC by L2 norm.
//
//  Discussion:
//
//    The L2 norm of A+Bi is sqrt ( A * A + B * B ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input/output, complex <float> X[N].
//    On input, an unsorted array.
//    On output, X has been sorted.
//
{
  int i;
  int indx;
  int isgn;
  int j;
  float normsq_i;
  float normsq_j;
  complex <float> temp;

  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;

  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );

    if ( 0 < indx )
    {
      temp = x[i-1];
      x[i-1] = x[j-1];
      x[j-1] = temp;
    }
    else if ( indx < 0 )
    {
      normsq_i = pow ( real ( x[i-1] ), 2 )
               + pow ( imag ( x[i-1] ), 2 );

      normsq_j = pow ( real ( x[j-1] ), 2 )
               + pow ( imag ( x[j-1] ), 2 );

      if ( normsq_i < normsq_j )
      {
        isgn = -1;
      }
      else
      {
        isgn = +1;
      }
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

complex <float> *c4vec_spiral ( int n, int m, complex <float> c1, 
  complex <float> c2 )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_SPIRAL returns N points on a spiral between C1 and C2.
//
//  Discussion:
//
//    A C4VEC is a vector of C8's.
//
//    Let the polar form of C1 be ( R1, T1 ) and the polar form of C2 
//    be ( R2, T2 ) where, if necessary, we increase T2 by 2*PI so that T1 <= T2.
//    
//    Then the polar form of the I-th point C(I) is:
//
//      R(I) = ( ( N - I     ) * R1 
//             + (     I - 1 ) * R2 ) 
//              / ( N    - 1 )
//
//      T(I) = ( ( N - I     ) * T1 
//             + (     I - 1 ) * ( T2 + M * 2 * PI ) ) 
//             / ( N     - 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points on the spiral.
//
//    Input, int M, the number of full circuits the 
//    spiral makes.
//
//    Input, complex <float> C1, C2, the first and last points 
//    on the spiral.
//
//    Output, complex <float> C4VEC_SPIRAL_NEW[N], the points.
//
{
  complex <float> *c;
  int i;
  float r1;
  float r2;
  float ri;
  float r4_pi = 3.141592653589793;
  float t1;
  float t2;
  float ti;

  c = new complex <float>[n];

  r1 = c4_abs ( c1 );
  r2 = c4_abs ( c2 );

  t1 = c4_arg ( c1 );
  t2 = c4_arg ( c2 );

  if ( m == 0 )
  {
    if ( t2 < t1 )
    {
      t2 = t2 + 2.0 * r4_pi;
    }
  }
  else if ( 0 < m )
  {
    if ( t2 < t1 )
    {
      t2 = t2 + 2.0 * r4_pi;
    }
    t2 = t2 + ( float ) ( m ) * 2.0 * r4_pi;
  }
  else if ( m < 0 )
  {
    if ( t1 < t2 )
    {
      t2 = t2 - 2.0 * r4_pi;
    }
    t2 = t2 - ( float ) ( m ) * 2.0 * r4_pi;
  }

  for ( i = 0; i < n; i++ )
  {
    ri = ( ( float ) ( n - i - 1 ) * r1
         + ( float ) (     i     ) * r2 )
         / ( float ) ( n     - 1 );

    ti = ( ( float ) ( n - i - 1 ) * t1
         + ( float ) (     i     ) * t2 )
         / ( float ) ( n     - 1 );

    c[i] = polar_to_c4 ( ri, ti );
  }

  return c;
}
//****************************************************************************80

void c4vec_uniform_01 ( int n, int *seed, complex <float> c[] )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of values to compute.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C[N], the pseudorandom 
//    complex vector.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  float pi = 3.1415926;
  float r;
  float theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * complex <float> ( cos ( theta ), sin ( theta ) );
  }

  return;
}
//****************************************************************************80

complex <float> *c4vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNIFORM_01_NEW returns a unit pseudorandom C4VEC.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of values to compute.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4VEC_UNIFORM_01_NEW[N], the pseudorandom 
//    complex vector.
//
{
  complex <float> *c;
  int i;
  int i4_huge = 2147483647;
  int k;
  float pi = 3.1415926;
  float r;
  float theta;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "C4VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  c = new complex <float> [n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r = sqrt ( ( float ) ( *seed ) * 4.656612875E-10 );

    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    theta = 2.0 * pi * ( ( float ) ( *seed ) * 4.656612875E-10 );

    c[i] = r * complex <float> ( cos ( theta ), sin ( theta ) );
  }

  return c;
}
//****************************************************************************80

complex <float> *c4vec_unity ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNITY returns the N roots of unity in a C4VEC.
//
//  Discussion:
//
//    A C4VEC is a vector of complex <float> values.
//
//    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
//
//    X(1:N)^N = ( (1,0), (1,0), ..., (1,0) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, complex <float> C4VEC_UNITY[N], the N roots of unity.
//
{
  complex <float> *a;
  int i;
  float pi = 3.141592653589793;
  float theta;

  a = new complex <float> [n];

  for ( i = 0; i < n; i++ )
  {
    theta = pi * ( float ) ( 2 * i ) / ( float ) ( n );
    a[i] = complex <float> ( cos ( theta ), sin ( theta ) );
  }

  return a;
}
//****************************************************************************80

complex <float> cartesian_to_c4 ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    CARTESIAN_TO_C4 converts a Cartesian form to a C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, Y, the Cartesian form.
//
//    Output, complex <float> CARTESIAN_TO_C4, the complex number.
//
{
  complex <float> c;

  c = complex <float> ( x, y );

  return c;
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

complex <float> polar_to_c4 ( float r, float theta )

//****************************************************************************80
//
//  Purpose:
//
//    POLAR_TO_C4 converts a polar form to a C4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R, THETA, the polar form.
//
//    Output, complex <float> POLAR_TO_C4, the complex number.
//
{
  complex <float> c;

  c = complex <float> ( r * cos ( theta ), r * sin ( theta ) );

  return c;
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

complex <float> r4_csqrt ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CSQRT returns the complex square root of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose square root is desired.
//
//    Output, complex <float> R4_CSQRT, the square root of X:
//
{
  float argument;
  float magnitude;
  float pi = 3.141592653589793;
  complex <float> value;

  if ( 0.0 < x )
  {
    magnitude = x;
    argument = 0.0;
  }
  else if ( 0.0 == x )
  {
    magnitude = 0.0;
    argument = 0.0;
  }
  else if ( x < 0.0 )
  {
    magnitude = -x;
    argument = pi;
  }

  magnitude = sqrt ( magnitude );
  argument = argument / 2.0;

  value = magnitude * complex <float> ( cos ( argument ), sin ( argument ) );

  return value;
}
//****************************************************************************80

float r4_floor ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_FLOOR rounds an R4 "down" (towards -oo) to the next integer.
//
//  Example:
//
//    X        R4_FLOOR(X)
//
//   -1.1      -2.0
//   -1.0      -1.0
//   -0.9      -1.0
//   -0.1      -1.0
//    0.0       0.0
//    0.1       0.0
//    0.9       0.0
//    1.0       1.0
//    1.1       1.0
//    2.9       2.0
//    3.0       3.0
//    3.14159   3.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose floor is desired.
//
//    Output, float R4_FLOOR, the floor of X.
//
{
  float value;

  value = ( int ) x;

  if ( x < value )
  {
    value = value - 1.0;
  }

  return value;
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

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
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
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
  }
  return value;
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

float r4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_01 returns a unit pseudorandom R4.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      r4_uniform_01 = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R4_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  float r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( float ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r4poly2_root ( float a, float b, float c, complex <float> *r1,
  complex <float> *r2 )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY2_ROOT returns the two roots of a quadratic polynomial.
//
//  Discussion:
//
//    The polynomial has the form:
//
//      A * X^2 + B * X + C = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 October 2005
//
//  Parameters:
//
//    Input, float A, B, C, the coefficients of the polynomial.
//    A must not be zero.
//
//    Output, complex <float> *R1, *R2, the roots of the polynomial, which
//    might be real and distinct, real and equal, or complex conjugates.
//
{
  float disc;
  complex <float> q;
  float t;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R4POLY2_ROOT - Fatal error!\n";
    cerr << "  The coefficient A is zero.\n";
    exit ( 1 );
  }

  disc = b * b - 4.0 * a * c;
  t = - 0.5 * ( b + r4_sign ( b ) );
  q =  complex <float> ( t, 0.0 ) * r4_csqrt ( disc );
  *r1 = q / a;
  *r2 = c / q;

  return;
}
//****************************************************************************80

void r4poly3_root ( float a, float b, float c, float d,
  complex <float> *r1, complex <float> *r2, complex <float> *r3 )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY3_ROOT returns the three roots of a cubic polynomial.
//
//  Discussion:
//
//    The polynomial has the form
//
//      A * X^3 + B * X^2 + C * X + D = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2005
//
//  Parameters:
//
//    Input, float A, B, C, D, the coefficients of the polynomial.
//    A must not be zero.
//
//    Output, complex <float> *R1, *R2, *R3, the roots of the polynomial, which
//    will include at least one real root.
//
{
  complex <float> i;
  float pi = 3.141592653589793;
  float q;
  float r;
  float s1;
  float s2;
  complex <float> t;
  float temp;
  float theta;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R4POLY3_ROOT - Fatal error!\n";
    cerr << "  A must not be zero.\n";
    exit ( 1 );
  }

  i = complex <float> ( 0.0, 1.0 );

  q = ( pow ( b / a, 2 ) - 3.0 * ( c / a ) ) / 9.0;

  r = ( 2.0 * pow ( b / a, 3 ) - 9.0 * ( b / a ) * ( c / a )
      + 27.0 * ( d / a ) ) / 54.0;

  if ( r * r < q * q * q )
  {
    theta = acos ( r / sqrt ( pow ( q, 3 ) ) );
    *r1 = -2.0 * sqrt ( q ) * cos (   theta              / 3.0 );
    *r2 = -2.0 * sqrt ( q ) * cos ( ( theta + 2.0 * pi ) / 3.0 );
    *r3 = -2.0 * sqrt ( q ) * cos ( ( theta + 4.0 * pi ) / 3.0 );
  }
  else if ( q * q * q <= r * r )
  {
    temp = -r + sqrt ( r * r - q * q * q );
    s1 = r4_sign ( temp ) * pow ( r4_abs ( temp ), 1.0 / 3.0 );

    temp = -r - sqrt ( r * r - q * q * q );
    s2 = r4_sign ( temp ) * pow ( r4_abs ( temp ), 1.0 / 3.0 );

    *r1 = s1 + s2;
    *r2 = complex <float> ( -0.5 * ( s1 + s2 ),   0.5 * sqrt ( 3.0 ) * ( s1 - s2 ) );
    *r3 = complex <float> ( -0.5 * ( s1 + s2 ), - 0.5 * sqrt ( 3.0 ) * ( s1 - s2 ) );
  }

  t = complex <float> ( b / ( 3.0 * a ), 0.0 );

  *r1 = *r1 - t;
  *r2 = *r2 - t;
  *r3 = *r3 - t;

  return;
}
//****************************************************************************80

void r4poly4_root ( float a, float b, float c, float d, float e,
  complex <float> *r1, complex <float> *r2, complex <float> *r3,
  complex <float> *r4 )

//****************************************************************************80
//
//  Purpose:
//
//    R4POLY4_ROOT returns the four roots of a quartic polynomial.
//
//  Discussion:
//
//    The polynomial has the form:
//
//      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2005
//
//  Parameters:
//
//    Input, float A, B, C, D, the coefficients of the polynomial.
//    A must not be zero.
//
//    Output, complex <float> *R1, *R2, *R3, *R4, the roots of the polynomial.
//
{
  float a3;
  float a4;
  float b3;
  float b4;
  float c3;
  float c4;
  float d3;
  float d4;
  complex <float> p;
  complex <float> q;
  complex <float> r;
  complex <float> t;
  complex <float> zero;

  zero = 0.0;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "R4POLY4_ROOT - Fatal error!\n";
    cerr << "  A must not be zero.\n";
    exit ( 1 );
  }

  a4 = b / a;
  b4 = c / a;
  c4 = d / a;
  d4 = e / a;
//
//  Set the coefficients of the resolvent cubic equation.
//
  a3 = 1.0;
  b3 = -b4;
  c3 = a4 * c4 - 4.0 * d4;
  d3 = -a4 * a4 * d4 + 4.0 * b4 * d4 - c4 * c4;
//
//  Find the roots of the resolvent cubic.
//
  r4poly3_root ( a3, b3, c3, d3, r1, r2, r3 );
//
//  Choose one root of the cubic, here R1.
//
//  Set R = sqrt ( 0.25 * A4^2 - B4 + R1 )
//
  t = complex <float> ( 0.25 * a4 * a4 - b4, 0.0 );
  r = c4_sqrt (  t + ( *r1 ) );

  if ( real ( r ) != 0.0 || imag ( r ) != 0.0 )
  {
    p = c4_sqrt ( - r * r + complex <float> ( 0.75 * a4 * a4 - 2.0 * b4, 0.0 )
        + complex <float> ( 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) ) / r );

    q = c4_sqrt ( - r * r + complex <float> ( 0.75 * a4 * a4 - 2.0 * b4, 0.0 )
        - complex <float> ( 0.25 * ( 4.0 * a4 * b4 - 8.0 * c4 - a4 * a4 * a4 ) ) / r );
  }
  else
  {
    p = c4_sqrt ( complex <float> ( 0.75 * a4 * a4 - 2.0 * b4, 0.0 )
      + complex <float> ( 2.0, 0.0 ) 
      * c4_sqrt ( (*r1) * (*r1) - complex <float> ( 4.0 * d4, 0.0 ) ) );

    q = c4_sqrt ( complex <float> ( 0.75 * a4 * a4 - 2.0 * b4, 0.0 )
      - complex <float> ( 2.0, 0.0 )
      * c4_sqrt ( (*r1) * (*r1) - complex <float> ( 4.0 * d4, 0.0 ) ) );
  }
//
//  Set the roots.
//
  t = complex <float> ( 0.5, 0.0 );
  *r1 = complex <float> ( -0.25 * a4, 0.0 ) + t * r + t * p;
  *r2 = complex <float> ( -0.25 * a4, 0.0 ) + t * r - t * p;
  *r3 = complex <float> ( -0.25 * a4, 0.0 ) - t * r + t * q;
  *r4 = complex <float> ( -0.25 * a4, 0.0 ) - t * r - t * q;

  return;
}
//****************************************************************************80

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 ) 
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
