# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>
# include <cstring>

using namespace std;

# include "subpak.hpp"

//****************************************************************************80

double angle_shift ( double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_SHIFT shifts angle ALPHA to lie between BETA and BETA+2PI.
//
//  Discussion:
//
//    The input angle ALPHA is shifted by multiples of 2 * PI to lie
//    between BETA and BETA+2*PI.
//
//    The resulting angle GAMMA has all the same trigonometric function
//    values as ALPHA.
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
//    Input, double ALPHA, the angle to be shifted.
//
//    Input, double BETA, defines the lower endpoint of
//    the angle range.
//
//    Output, double ANGLE_SHIFT, the shifted angle.
//
{
  double gamma;
  double pi = 3.141592653589793;

  if ( alpha < beta )
  {
    gamma = beta - fmod (  beta - alpha, 2.0 * pi ) + 2.0 * pi;
  }
  else
  {
    gamma = beta + fmod ( alpha - beta, 2.0 * pi );
  }

  return gamma;
}
//****************************************************************************80

double angle_shift_deg ( double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_SHIFT_DEG shifts angle ALPHA to lie between BETA and BETA+360.
//
//  Discussion:
//
//    The input angle ALPHA is shifted by multiples of 360 to lie
//    between BETA and BETA+360.
//
//    The resulting angle GAMMA has all the same trigonometric function
//    values as ALPHA.
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
//    Input, double ALPHA, the angle to be shifted.
//
//    Input, double BETA, defines the lower endpoint of
//    the angle range.
//
//    Output, double ANGLE_SHIFT, the shifted angle.
//
{
  double gamma;

  if ( alpha < beta )
  {
    gamma = beta - fmod ( beta - alpha, 360.0 ) + 360.0;
  }
  else
  {
    gamma = beta + fmod ( alpha - beta, 360.0 );
  }

  return gamma;
}
//****************************************************************************80

double *angle_to_rgb ( double angle )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, the angle in the color hexagon.  The sextants are
//    defined by the following points:
//        0 degrees, 1, 0, 0, red;
//       60 degrees, 1, 1, 0, yellow;
//      120 degrees, 0, 1, 0, green;
//      180 degrees, 0, 1, 1, cyan;
//      240 degrees, 0, 0, 1, blue;
//      300 degrees, 1, 0, 1, magenta.
//
//    Output, double ANGLE_TO_RGB[3], the RGB specifications for the 
//    color that lies at the given angle, on the perimeter of the 
//    color hexagon.  One value will be 1, and one value will be 0.
//
{
# define DEGREES_TO_RADIANS ( 3.141592653589793 / 180.0 )

  double angle2;
  double *rgb;

  rgb = new double[3];

  angle = r8_modp ( angle, 360.0 );

  if ( angle <= 60.0 )
  {
    angle2 = DEGREES_TO_RADIANS * 3.0 * angle / 4.0;
    rgb[0] = 1.0;
    rgb[1] = tan ( angle2 );
    rgb[2] = 0.0;
  }
  else if ( angle <= 120.0 )
  {
    angle2 = DEGREES_TO_RADIANS * 3.0 * angle / 4.0;
    rgb[0] = cos ( angle2 ) / sin ( angle2 );
    rgb[1] = 1.0;
    rgb[2] = 0.0;
  }
  else if ( angle <= 180.0 )
  {
    angle2 = DEGREES_TO_RADIANS * 3.0 * ( angle - 120.0 ) / 4.0;
    rgb[0] = 0.0;
    rgb[1] = 1.0;
    rgb[2] = tan ( angle2 );
  }
  else if ( angle <= 240.0 )
  {
    angle2 = DEGREES_TO_RADIANS * 3.0 * ( angle - 120.0 ) / 4.0;
    rgb[0] = 0.0;
    rgb[1] = cos ( angle2 ) / sin ( angle2 );
    rgb[2] = 1.0;
  }
  else if ( angle <= 300.0 )
  {
    angle2 = DEGREES_TO_RADIANS * 3.0 * ( angle - 240.0 ) / 4.0;
    rgb[0] = tan ( angle2 );
    rgb[1] = 0.0;
    rgb[2] = 1.0;
  }
  else if ( angle <= 360.0 )
  {
    angle2 = DEGREES_TO_RADIANS * 3.0 * ( angle - 240.0 ) / 4.0;
    rgb[0] = 1.0;
    rgb[1] = 0.0;
    rgb[2] = cos ( angle2 ) / sin ( angle2 );
  }

  return rgb;
# undef DEGREES_TO_RADIANS
}
//****************************************************************************80

void axis_limits ( double xmin, double xmax, int ndivs, double *pxmin, 
  double *pxmax, double *pxdiv, int *nticks )

//****************************************************************************80
//
//  Purpose:
//
//    AXIS_LIMITS returns "nice" axis limits for a plot.
//
//  Discussion:
//
//    The routine is given information about the range of a variable, and
//    the number of divisions desired.  It returns suggestions for
//    labeling a plotting axis for the variable, including the
//    starting and ending points, the length of a single division,
//    and a suggested tick marking for the axis.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XMIN, XMAX, the lower and upper values that must be
//    included on the axis.  XMIN must be less than XMAX.
//
//    Input, int NDIVS, the number of divisions desired along
//    the axis.
//
//    Output, double *PXMIN, *PXMAX, the recommended lower and upper axis
//    bounds.  It will be the case that PXMIN <= XMIN < XMAX <= PXMAX.
//
//    Output, double *PXDIV, the recommended size of a single division.
//
//    Output, int *NTICKS, a suggested number of ticks to use,
//    if subdividing each of the NDIVS divisions of the axis.
//
{
# define NSTEPS 5

  double best;
  double good;
  int i;
  int ihi;
  int ilo;
  int intlog;
  int iticks[NSTEPS] = { 5, 4, 4, 5, 5 };
  int ival;
  int j;
  double pxmax2;
  double pxmin2;
  double pxdiv2;
  double reldif;
  double steps[NSTEPS] = { 1.0,  2.0,  4.0,  5.0, 10.0 };
  double temp;

  if ( xmin == xmax )
  {
    xmin = xmin - 0.5;
    xmax = xmax + 0.5;
  }
  else if ( xmax < xmin )
  {
    temp = xmin;
    xmin = xmax;
    xmax = temp;
  }

  if ( ndivs <= 0 )
  {
    ndivs = 5;
  }
//
//  Set RELDIF, the size of the X interval divided by the largest X.
//
  if ( xmax != xmin )
  {
    reldif = ( xmax - xmin ) / r8_max ( r8_abs ( xmax ), r8_abs ( xmin ) );
  }
  else
  {
    reldif = 0.0;
  }
//
//  If RELDIF tells us that XMIN and XMAX are extremely close,
//  do some simple things.
//
  if ( reldif < 0.00001 )
  {
    if ( xmax == 0.0 )
    {
      *pxdiv = 1.0;
    }
    else
    {
      intlog = ( int ) ( r8_log_10 ( xmax ) );

      if ( intlog < 0 )
      {
        intlog = intlog - 1;
      }

      *pxdiv = pow ( 10.0, intlog );

      if ( 1.0 < *pxdiv )
      {
        *pxdiv = 1.0;
      }

    }

    *nticks = 5;
    *pxmin = xmax - ( double ) ( ndivs / 2 ) * (*pxdiv);
    *pxmax = xmax + ( double ) ( ndivs - ( ndivs / 2 ) ) * (*pxdiv);
  }
//
//  But now handle the more general case, when XMIN and XMAX
//  are relatively far apart.
//
  else
  {
    best = -999.0;
//
//  On second loop, increase INTLOG by 1.
//
    for ( j = 1; j <= 2; j++ )
    {
//
//  Compute INTLOG, roughly the logarithm base 10 of the range
//  divided by the number of divisions.
//
      intlog = ( int ) ( r8_log_10 ( ( xmax - xmin ) / ( double ) ( ndivs ) ) ) 
        + ( j - 1 );

      if ( xmax - xmin  < ( double ) ( ndivs ) )
      {
        intlog = intlog - 1;
      }
//
//  Now consider taking 1, 2, 4, 5 or 10 steps of size 10**INTLOG:
//
      for ( i = 1; i <= NSTEPS; i++ )
      {
//
//  Compute the size of each step.
//
        pxdiv2 = steps[i-1] * pow ( 10.0, intlog );
//
//  Make sure NDIVS steps can reach from XMIN to XMAX, at least.
//
        if ( xmax <= xmin + ndivs * pxdiv2 )
        {
//
//  Now decide where to start the axis.
//  Start the axis at PXMIN2, to the left of XMIN, and
//  representing a whole number of steps of size PXDIV2.
//
          if ( 0.0 <= xmin )
          {
            ival = ( int ) ( xmin / pxdiv2 );
          }
          else
          {
            ival = ( int ) ( xmin / pxdiv2 ) - 1;
          }

          pxmin2 = ival * pxdiv2;
//
//  PXMAX2 is, of course, NDIVS steps above PXMIN2.
//
          pxmax2 = pxmin2 + ndivs * pxdiv2;
//
//  Only consider going on if PXMAX2 is at least XMAX.
//
          if ( xmax <= pxmax2 )
          {
//
//  Now judge this grid by the relative amount of wasted axis length.
//
            good = ( xmax - xmin ) / ( pxmax2 - pxmin2 );

            if ( best < good )
            {
              best = good;
              *pxmax = pxmax2;
              *pxmin = pxmin2;
              *pxdiv = pxdiv2;
              *nticks = iticks[i-1];
            }
          }
        }
      }
    }
  }
//
//  If necessary, adjust the locations of PXMIN and PXMAX so that the
//  interval is more symmetric in containing XMIN through XMAX.
//
  for ( ; ; )
  {
    ilo = ( int ) ( ( xmin - *pxmin ) / (*pxdiv) );
    ihi = ( int ) ( ( *pxmax - xmax ) / (*pxdiv) );

    if ( ihi < ilo + 2 )
    {
      break;
    }

    *pxmin = *pxmin - (*pxdiv);
    *pxmax = *pxmax - (*pxdiv);

  }

  return;
# undef NSTEPS
}
//****************************************************************************80

int bar_check ( int digit[12] )

//****************************************************************************80
//
//  Purpose:
//
//    BAR_CHECK computes the check digit for a barcode.
//
//  Discussion:
//
//    CHECK = SUM ( I = 1, 11, by 2's ) DIGIT(I)
//       + 3 * SUM ( I = 2, 10, by 2's ) DIGIT(I)
//
//    CHECK = MOD ( 10 - MOD ( CHECK, 10 ), 10 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIGIT[12], entries 1 through 11 of DIGIT contain
//    the digits of the bar code.  Each entry must be between 0 and 9.
//    The 12th digit should be the check digit.
//
//    Output, int BAR_CHECK, the correct check digit.  If the bar code
//    is correct, then DIGIT(12) should equal BAR_CHECK.
//
{
  int check;

  check =  
          ( digit[0] + digit[2] + digit[4] + digit[6] + digit[8] + digit[10] )
    + 3 * ( digit[1] + digit[3] + digit[5] + digit[7] + digit[9] );

  check =  ( 10 - ( check % 10 ) ) % 10;

  return check;
}
//****************************************************************************80

char *bar_code ( int digit[] )

//****************************************************************************80
//
//  Purpose:
//
//    BAR_CODE constructs the 113 character barcode from 11 digits.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int DIGIT(12).
//    On input, the first 11 entries of DIGIT contain a code to be
//    turned into a barcode.
//    On output, the 12-th entry of DIGIT is a check digit.
//
//    Output, char BAR_CODE[114], the 113 character bar code corresponding 
//    to the digit information.
//
{
  char *bar;
  int check;
  char *codel;
  char *coder;
  int i;

  bar = new char[114];
//
//  9 character quiet zone.
//
  strcpy ( bar, "000000000" );
//
//  3 character guard pattern.
//
  strcpy ( bar+9, "101" );
//
//  7 character product category.
//
  codel = bar_digit_code_left ( digit[0] );
  strcpy ( bar+12, codel );
  delete [] codel;
//
//  35 characters contain the 5 digit manufacturer code.
//
  for ( i = 0; i < 5; i++ )
  {
    codel = bar_digit_code_left ( digit[i+1] );
    strcpy ( bar+19+i*7, codel );
    delete [] codel;
  }
//
//  Center guard pattern.
//
  strcpy ( bar+54, "01010" );
//
//  35 characters contain the 5 digit product code.
//
  for ( i = 0; i < 5; i++ )
  {
    coder = bar_digit_code_right ( digit[i+6] );
    strcpy ( bar+59+i*7, coder );
    delete [] coder;
  }
//
//  Compute the check digit.
//
  check = bar_check ( digit );
  digit[11] = check;

  coder = bar_digit_code_right ( check );
  strcpy ( bar+94, coder );
//
//  Guard pattern.
//
  strcpy ( bar+101, "101" );
//
//  Quiet zone.
//
  strcpy ( bar+104, "000000000" );
  bar[113] = '\0';

  return bar;
}
//****************************************************************************80

char *bar_digit_code_left ( int digit )

//****************************************************************************80
//
//  Purpose:
//
//    BAR_DIGIT_CODE_LEFT returns the 7 character left bar code for a digit.
//
//  Example:
//
//    DIGIT = 3
//    CODEL = '0111101'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIGIT, the digit, between 0 and 9.
//
//    Output, char BAR_CODE_DIGIT_LEFT[8], the 7 character left code for the digit.
//
{
  char *codel;

  codel = new char[8];

  if ( digit == 0 )
  {
    strcpy ( codel, "0001101" );
  }
  else if ( digit == 1 )
  {
    strcpy ( codel, "0011001" );
  }
  else if ( digit == 2 )
  {
    strcpy ( codel, "0010011" );
  }
  else if ( digit == 3 )
  {
    strcpy ( codel, "0111101" );
  }
  else if ( digit == 4 )
  {
    strcpy ( codel, "0100011" );
  }
  else if ( digit == 5 )
  {
    strcpy ( codel, "0110001" );
  }
  else if ( digit == 6 )
  {
    strcpy ( codel, "0101111" );
  }
  else if ( digit == 7 )
  {
    strcpy ( codel, "0111011" );
  }
  else if ( digit == 8 )
  {
    strcpy ( codel, "0110111" );
  }
  else if ( digit == 9 )
  {
    strcpy ( codel, "0001011" );
  }
  else
  {
    strcpy ( codel, "???????" );
  }

  return codel;
}
//****************************************************************************80

char *bar_digit_code_right ( int digit )

//****************************************************************************80
//
//  Purpose:
//
//    BAR_DIGIT_CODE_RIGHT returns the 7 character right bar code for a digit.
//
//  Example:
//
//    DIGIT = 3
//    CODER = "1000010"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIGIT, the digit, between 0 and 9.
//
//    Output, char BAR_DIGIT_CODE_RIGHT[8], the 7 character right code.
//
{
  char *coder;

  coder = new char[8];

  if ( digit == 0 )
  {
    strcpy ( coder, "1110010" );
  }
  else if ( digit == 1 )
  {
    strcpy ( coder, "1100110" );
  }
  else if ( digit == 2 )
  {
    strcpy ( coder, "1101100" );
  }
  else if ( digit == 3 )
  {
    strcpy ( coder, "1000010" );
  }
  else if ( digit == 4 )
  {
    strcpy ( coder, "1011100" );
  }
  else if ( digit == 5 )
  {
    strcpy ( coder, "1001110" );
  }
  else if ( digit == 6 )
  {
    strcpy ( coder, "1010000" );
  }
  else if ( digit == 7 )
  {
    strcpy ( coder, "1000100" );
  }
  else if ( digit == 8 )
  {
    strcpy ( coder, "1001000" );
  }
  else if ( digit == 9 )
  {
    strcpy ( coder, "1110100" );
  }
  else
  {
    strcpy ( coder, "???????" );
  }

  return coder;
}
//****************************************************************************80

double bmi_english ( double w_lb, double h_ft, double h_in )

//****************************************************************************80
//
//  Purpose:
//
//    BMI_ENGLISH computes the body mass index given English measurements.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double W_LB, the body weight in pounds.
//
//    Input, double H_FT, H_IN, the body height in feet and inches
//
//    Output, double BMI_ENGLISH, the body mass index.
//
{
  double h_m;
  double value;
  double w_kg;

  w_kg = pounds_to_kilograms ( w_lb );

  h_m = feet_to_meters ( h_ft + ( h_in / 12.0 ) );

  value = bmi_metric ( w_kg, h_m );

  return value;
}
//****************************************************************************80

double bmi_metric ( double w_kg, double h_m )

//****************************************************************************80
//
//  Purpose:
//
//    BMI_METRIC computes the body mass index given metric measurements.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double W_KG, the body weight in kilograms.
//
//    Input, double H_M, the body height in meters.
//
//    Output, double BMI_METRIC, the body mass index.
//
{
  double value;

  value = w_kg / ( h_m * h_m );

  return value;
}
//****************************************************************************80

double c8_argument ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_ARGUMENT returns the argument of a C8.
//
//  Discussion:
//
//    A C8 is a double precision complex value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> X, the value whose argument is desired.
//
//    Output, double C8_ARGUMENT, the argument of X.
//
{
  double value;

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

double c8_magnitude ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_MAGNITUDE returns the magnitude of a C8.
//
//  Discussion:
//
//    A C8 is a double precision complex value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> X, the value whose norm is desired.
//
//    Output, double C8_MAGNITUDE, the magnitude of X.
//
{
  double magnitude;

  magnitude = sqrt ( pow ( real ( x ), 2 ) + pow ( imag ( x ), 2 ) );

  return magnitude;
}
//****************************************************************************80

complex <double> c8_sqrt ( complex <double> x )

//****************************************************************************80
//
//  Purpose:
//
//    C8_SQRT returns the principal square root of a C8.
//
//  Discussion:
//
//    A C8 is a double precision complex value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> X, the number whose square root is desired.
//
//    Output, complex <double> C8_SQRT, the square root of X.
//
{
  double argument;
  double magnitude;
  complex <double> value;

  argument = c8_argument ( x );
  magnitude = c8_magnitude ( x );

  if ( magnitude == 0.0 )
  {
    value = complex <double> ( 0.0, 0.0 );
  }
  else
  {
    value = sqrt ( magnitude ) 
      * complex <double> ( cos ( argument / 2.0 ), sin ( argument / 2.0 ) );
  }

  return value;
}
//****************************************************************************80

char ch_cap ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= c && c <= 122 ) 
  {
    c = c - 32;
  }   

  return c;
}
//****************************************************************************80

bool ch_eqi ( char ch1, char ch2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH1, CH2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  bool value;

  if ( 97 <= ch1 && ch1 <= 122 ) 
  {
    ch1 = ch1 - 32;
  } 
  if ( 97 <= ch2 && ch2 <= 122 ) 
  {
    ch2 = ch2 - 32;
  }     

  value = ( ch1 == ch2 );

  return value;
}
//****************************************************************************80

bool ch_is_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_DIGIT returns TRUE if a character is a decimal digit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to be analyzed.
//
//    Output, bool CH_IS_DIGIT, is TRUE if C is a digit.
//
{
  bool value;

  if ( '0' <= c && c <= '9' )
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

int ch_to_digit ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     CH  DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the 
//    character was 'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

int charstar_len_trim ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    CHARSTAR_LEN_TRIM returns the length of a CHAR* to the last nonblank.
//
//  Discussion:
//
//    This function used to be called S_LEN_TRIM.  However, it seems preferable
//    to use the STRING class rather than CHAR*.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int CHARSTAR_LEN_TRIM, the length of the string to the last nonblank.
//    If CHARSTAR_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n ) 
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

double degrees_to_radians ( double degrees )

//****************************************************************************80
//
//  Purpose:
//
//    DEGREES_TO_RADIANS converts an angle measure from degrees to radians.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double DEGREES, the angle measure in degrees.
//
//    Output, double DEGREES_TO_RADIANS, the angle measure in radians.
//
{
  double value;

  value = ( degrees / 180.0 ) * 3.141592653589793;

  return value;
}
//****************************************************************************80

double e_constant ( )

//****************************************************************************80
//
//  Purpose:
//
//    E_CONSTANT returns the value of E.
//
//  Discussion:
//
//    "E" was named in honor of Euler, but is known as Napier's constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double E_CONSTANT, the base of the natural 
//    logarithm system.
//
{
  double value = 2.718281828459045;

  return value;
}
//****************************************************************************80

double euler_constant ( )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
//
//  Discussion:
//
//    The Euler-Mascheroni constant is often denoted by a lower-case
//    Gamma.  Gamma is defined as
//
//      Gamma = limit ( M -> oo ) ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EULER_CONSTANT, the value of the 
//    Euler-Mascheroni constant.
//
{
  double value = 0.5772156649015328;

  return value;
}
//****************************************************************************80

void fac_div ( int prime_num, int npower1[], int npower2[], int npower3[] )

//****************************************************************************80
//
//  Purpose:
//
//    FAC_DIV divides two quantities represented as prime factors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PRIME_NUM, the index of the highest prime number
//    used in the representations.
//
//    Input, int NPOWER1[PRIME_NUM], the powers of primes
//    in the representation of the first quantity.
//
//    Input, int NPOWER2[PRIME_NUM], the powers of primes
//    in the representation of the second quantity.
//
//    Output, int NPOWER3[PRIME_NUM], the powers of primes
//    in the representation of the quotient.
//
{
  int i;

  for ( i = 0; i < prime_num; i++ )
  {
    npower3[i] = npower1[i] - npower2[i];
  }

  return;
}
//****************************************************************************80

void fac_gcd ( int prime_num, int npower1[], int npower2[], int npower3[] )

//****************************************************************************80
//
//  Purpose:
//
//    FAC_GCD finds the GCD of two products of prime factors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PRIME_NUM, the index of the highest prime number
//    used in the representations.
//
//    Input, int NPOWER1[PRIME_NUM], the powers of primes
//    in the representation of the first quantity.  All the powers
//    must be nonnegative.
//
//    Input, int NPOWER2[PRIME_NUM], the powers of primes
//    in the representation of the second quantity.  All the powers
//    must be nonnegative.
//
//    Output, int NPOWER3[PRIME_NUM], the powers of primes
//    in the representation of the GCD.
//
{
  int i;

  for ( i = 0; i < prime_num; i++ )
  {
    if ( npower1[i] < 0 )
    {
      cerr << "\n";
      cerr << "FAC_GCD - Fatal error!\n";
      cerr << "  One of the powers is negative.\n";
      exit ( 1 );
    }

    if ( npower2[i] < 0 )
    {
      cerr << "\n";
      cerr << "FAC_GCD - Fatal error!\n";
      cerr << "  One of the powers is negative.\n";
      exit ( 1 );
    }

    npower3[i] = i4_min ( npower1[i], npower2[i] );
  }

  return;
}
//****************************************************************************80

void fac_lcm ( int prime_num, int npower1[], int npower2[], int npower3[] )

//****************************************************************************80
//
//  Purpose:
//
//    FAC_LCM finds the LCM of two products of prime factors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PRIME_NUM, the index of the highest prime number
//    used in the representations.
//
//    Input, int NPOWER1[PRIME_NUM], the powers of primes
//    in the representation of the first quantity.
//
//    Input, int NPOWER2[PRIME_NUM], the powers of primes
//    in the representation of the second quantity.
//
//    Output, int NPOWER3[PRIME_NUM], the powers of primes
//    in the representation of the LCM.
//
{
  int i;

  for ( i = 0; i < prime_num; i++ )
  {
    if ( npower1[i] < 0 )
    {
      cerr << "\n";
      cerr << "FAC_LCM - Fatal error!\n";
      cerr << "  One of the powers is negative.\n";
      exit ( 1 );
    }

    if ( npower2[i] < 0 )
    {
      cerr << "\n";
      cerr << "FAC_LCM - Fatal error!\n";
      cerr << "  One of the powers is negative.\n";
      exit ( 1 );
    }

    npower3[i] = i4_max ( npower1[i], npower2[i] );
  }

  return;
}
//****************************************************************************80

void fac_mul ( int prime_num, int npower1[], int npower2[], int npower3[] )

//****************************************************************************80
//
//  Purpose:
//
//    FAC_MUL multiplies two quantities represented as prime factors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PRIME_NUM, the index of the highest prime number
//    used in the representations.
//
//    Input, int NPOWER1[PRIME_NUM], the powers of primes
//    in the representation of the first quantity.
//
//    Input, int NPOWER2[PRIME_NUM], the powers of primes
//    in the representation of the second quantity.
//
//    Output, int NPOWER3[PRIME_NUM], the powers of primes
//    in the representation of the product.
//
{
  int i;

  for ( i = 0; i < prime_num; i++ )
  {
    npower3[i] = npower1[i] + npower2[i];
  }
  return;
}
//****************************************************************************80

void fac_print ( int prime_num, int npower[] )

//****************************************************************************80
//
//  Purpose:
//
//    FAC_PRINT prints a product of prime factors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PRIME_NUM, the index of the highest prime number
//    used in the representations.
//
//    Input, int NPOWER[PRIME_NUM], the powers of primes
//    in the representation of the quantity.
//
{
  int i;

  cout << "\n";
  cout << "   Prime     Power\n";
  cout << "\n";

  for ( i = 0; i < prime_num; i++ )
  {
    if ( npower[i] != 0 )
    {
      cout << "  " << setw(8) << prime ( i+1 )
           << "  " << setw(8) << npower[i] << "\n";
    }
  }

  return;
}
//****************************************************************************80

int fac_to_i4 ( int prime_num, int npower[] )

//****************************************************************************80
//
//  Purpose:
//
//    FAC_TO_I4 converts a product of prime factors into an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PRIME_NUM, the index of the highest prime number
//    used in the representations.
//
//    Input, int NPOWER[PRIME_NUM], the powers of primes
//    in the representation of the quantity.  If any of these powers
//    are negative, then INTVAL will be set to 0.
//
//    Output, int FAC_TO_I4, the integer represented by the product of the
//    prime factors.
//
{
  int factor;
  int i;
  int intval;
  int j;

  intval = 1;

  for ( i = 0; i < prime_num; i++ )
  {
    if ( npower[i] < 0 )
    {
      cerr << "\n";
      cerr << "FAC_TO_I4 - Fatal error!\n";
      cerr << "  One of the powers is negative.\n";
      exit ( 1 );
    }
    factor = prime ( i+1 );

    for ( j = 1; j <= npower[i]; j++ )
    {
      intval = intval * factor;
    }
  }

  return intval;
}
//****************************************************************************80

void fac_to_rat ( int prime_num, int npower[], int *top, int *bot )

//****************************************************************************80
//
//  Purpose:
//
//    FAC_TO_RAT converts a prime factorization into a rational value.
//
//  Example:
//
//    Start with the prime factorization representation:
//
//      40/9 = 2**3 * 3**(-2) * 5
//
//    Input:
//
//      NPOWER = ( 3, -2, 1 )
//
//    Output:
//
//      TOP = 40 ( = 2**3 * 5**1 = PRIME(1)**3                 * PRIME(3)**1 )
//      BOT = 9  ( = 3**2        =               PRIME(2)**2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PRIME_NUM, the index of the highest prime number
//    used in the representations.
//
//    Input, int NPOWER[PRIME_NUM].  NPOWER(I) is the power of
//    the I-th prime in the prime factorization.  NPOWER(I) may
//    be positive or negative.
//
//    Output, int TOP, BOT, the top and bottom of a rational value.
//
{
  int i;

  *top = 1;
  *bot = 1;

  for ( i = 0; i < prime_num; i++ )
  {
    if ( 0 < npower[i] )
    {
      *top = *top * i4_power ( prime ( i + 1 ), npower[i] );
    }
    else if ( npower[i] < 0 )
    {
      *bot = *bot * i4_power ( prime ( i + 1 ), -npower[i] );
    }
  }

  return;
}
//****************************************************************************80

double feet_to_meters ( double ft )

//****************************************************************************80
//
//  Purpose:
//
//    FEET_TO_METERS converts a measurement in feet to meters.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double FT, the length in feet.
//
//    Output, double FEET_TO_METERS, the corresponding length in meters.
//
{
  double value;

  value = 0.0254 * 12.0 * ft;

  return value;
}
//****************************************************************************80

double gauss_sum ( int ndim, int n, double amplitude[], double center[], 
  double width[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    GAUSS_SUM evaluates a function that is the sum of Gaussians.
//
//  Discussion:
//
//    Gauss_Sum(X) = Sum ( 1 <= J <= Ngauss ) Amplitude(I) * exp ( -Arg )
//
//    where
//
//      Arg = sum ( 1 <= I <= NDIM ) ( ( ( X(I) - Center(I,J) ) / Width(J) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDIM, the spatial dimension.
//
//    Input, int N, the number of component Gaussian functions.
//
//    Input, double AMPLITUDE[N], CENTER[NDIM*N], WIDTH[N],
//    the amplitude, center and width for the component Gaussian functions.
//
//    Input, double X[NDIM], the point at which the function
//    is to be evaluated.
//
//    Output, double GAUSS_SUM, the value of the function.
//
{
  double arg;
  int i;
  int j;
  double value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    arg = 0.0;
    for ( i = 0; i < ndim; i++ )
    {
      arg = arg + pow ( ( x[i] - center[i+j*ndim] ) / width[j], 2 );
    }
    value = value + amplitude[j] * exp ( -arg );
  }

  return value;
}
//****************************************************************************80

unsigned long get_seed ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, unsigned long GET_SEED, a random seed value.
//
{
# define UNSIGNED_LONG_MAX 4294967295UL

  time_t clock;
  int hours;
  int minutes;
  int seconds;
  struct tm *lt;
  static unsigned long seed = 0;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  if ( seed == 0 )
  {
    clock = time ( &tloc );
    lt = localtime ( &clock );
//
//  Extract HOURS.
//
    hours = lt->tm_hour;
//
//  In case of 24 hour clocks, shift so that HOURS is between 1 and 12.
//
    if ( 12 < hours )
    {
      hours = hours - 12;
    }
//
//  Move HOURS to 0, 1, ..., 11
//
    hours = hours - 1;

    minutes = lt->tm_min;

    seconds = lt->tm_sec;

    seed = seconds + 60 * ( minutes + 60 * hours );
//
//  We want values in [1,43200], not [0,43199].
//
    seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,UNSIGNED_LONG_MAX].
//
    seed = ( unsigned long ) 
      ( ( ( double ) seed )
      * ( ( double ) UNSIGNED_LONG_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
  }
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;

# undef UNSIGNED_LONG_MAX
}
//****************************************************************************80

double *grid1 ( int dim_num, int nstep, double x1[], double x2[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID1 finds grid points between X1 and X2 in N dimensions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the points X1 and X2.
//
//    Input, int NSTEP, the number of points to be generated.
//    NSTEP must be at least 2.
//
//    Input, double X1[DIM_NUM], X2[DIM_NUM], the first and last
//    points, between which the equally spaced points are
//    to be computed.
//
//    Output, double X[DIM_NUM*NSTEP], the set of equally spaced
//    points.  Each column of X represents one point, with X[*,1] = X1
//    and X[*,NSTEP] = X2.
//
{
  int i;
  int j;
  double *x;

  if ( dim_num < 1 )
  {
    cerr << "\n";
    cerr << "GRID1 - Fatal error!\n";
    cerr << "  DIM_NUM < 1.\n";
    cerr << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  if ( nstep < 2 )
  {
    cerr << "\n";
    cerr << "GRID1 - Fatal error!\n";
    cerr << "  NSTEP < 2.\n";
    cerr << "  NSTEP = " << nstep << "\n";
    exit ( 1 );
  }

  x = new double[dim_num*nstep];

  for ( j = 1; j <= nstep; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x[i+(j-1)*dim_num] =
        ( ( double ) ( nstep - j     ) * x1[i]   
        + ( double ) (         j - 1 ) * x2[i] ) 
        / ( double ) ( nstep     - 1 );
    }
  }
  return x;
}
//****************************************************************************80

double *grid1n ( int j, int dim_num, int nstep, double x1[], double x2[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID1N finds the I-th grid point between X1 and X2 in N dimensions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int J, the number of the desired point.
//    Normally J would be between 1 and NSTEP, but that is
//    not necessary.  Note that J = 1 returns X1 and J = NSTEP
//    returns X2.
//
//    Input, int DIM_NUM, the dimension of the points X, X1 and X2.
//
//    Input, int NSTEP, this is the number of equally
//    spaced points that are between X1 and X2.  NSTEP must
//    be at least 2, because X1 and X2 are always included
//    in the set of points.
//
//    Input, double X1[DIM_NUM], X2[DIM_NUM], the first and last
//    points, between which the equally spaced points lie.
//
//    Output, double GRID1N[DIM_NUM], the J-th grid point between X1
//    and X2.
//
{
  int i;
  double *x;

  if ( dim_num < 1 )
  {
    cerr << "\n";
    cerr << "GRID1N - Fatal error!\n";
    cerr << "  DIM_NUM < 1.\n";
    cerr << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  if ( nstep < 2 )
  {
    cerr << "\n";
    cerr << "GRID1N - Fatal error!\n";
    cerr << "  NSTEP < 2.\n";
    cerr << "  NSTEP = " << nstep << "\n";
    exit ( 1 );
  }

  x = new double[dim_num];

  for ( i = 0; i < dim_num; i++ )
  {
    x[i] = ( ( double ) ( nstep - j     ) * x1[i]
           + ( double ) (         j - 1 ) * x2[i] )
           / ( double ) ( nstep     - 1 );
  }

  return x;
}
//****************************************************************************80

double *grid2 ( int j1, int j2, int dim_num, int nstep, double x1[], 
  double x2[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID2 computes grid points between X1 and X2 in N dimensions.
//
//  Discussion:
//
//    GRID2 computes grid points between X1 and X2 in N dimensions.
//
//    However, X1 need not be the first point computed, nor X2 the last.
//    The user must specify the steps on which X1 and X2 are passed
//    through.  These steps may even be outside the range of 1 through NSTEP.
//
//    We assume that a set of equally spaced points have
//    been drawn on the line through X1 and X2, and that
//    they have been numbered, with X1 labeled J1 and X2
//    labeled J2.  J1 or J2 may be between 1 and NSTEP,
//    in which case X1 or X2 will actually be returned in the
//    X array, but there is no requirement that J1 or J2
//    satisfy this condition.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int J1, J2.  J1 specifies the step on which
//    X1 would be computed, and similarly for J2.  
//    J1 and J2 must be distinct.
//
//    Input, int DIM_NUM, the dimension of the points X1 and X2.
//
//    Input, int NSTEP, this is the number of equally
//    spaced points that are to be generated.
//    NSTEP should be at least 1.
//
//    Input, double X1[DIM_NUM], X2[DIM_NUM], the points that define
//    the line along which the equally spaced points are generated, and
//    which may or may not be included in the set of computed points.
//
//    Output, double GRID2[DIM_NUM*NSTEP], the set of equally spaced
//    points.  Each column of X represents one point.
//    If 1 <= J1 <= NSTEP, then X(*,J1) = X1, and similarly for J2.
//
{
  int i;
  int j;
  double *x;

  if ( dim_num < 1 )
  {
    cerr << "\n";
    cerr << "GRID2 - Fatal error!\n";
    cerr << "  DIM_NUM < 1.\n";
    cerr << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  if ( j1 == j2 )
  {
    cerr << "\n";
    cerr << "GRID2 - Fatal error!\n";
    cerr << "  J1 = J2, leading to zero denominator.\n";
    cerr << "  J1 = " << j1 << "\n";
    cerr << "  J2 = " << j2 << "\n";
    exit ( 1 );
  }

  x = new double[nstep*dim_num];

  for ( j = 1; j <= nstep; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x[i+(j-1)*dim_num] = ( ( double ) ( j2 - j      ) * x1[i]   
                           + ( double ) (      j - j1 ) * x2[i] ) 
                           / ( double ) ( j2     - j1 );
    }
  }

  return x;
}
//****************************************************************************80

double *grid2n ( int j, int j1, int j2, int dim_num, double x1[], double x2[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID2N computes one grid point between X1 and X2 in N dimensions.
//
//  Discussion:
//
//    However, X1 need not be the first point computed, nor X2 the last.
//    The user must specify the steps on which X1 and X2 are passed through.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int J, the J coordinate of the desired point.
//    Note that if J = J1, X will be returned as X1, and if
//    J = J2, X will be returned as X2.
//
//    Input, int J1, J2.  J1 specifies the step on which
//    X1 would be computed, and similarly for J2.  That is,
//    we assume that a set of equally spaced points have
//    been drawn on the line through X1 and X2, and that
//    they have been numbered, with X1 labeled J1 and X2
//    labeled J2.  J1 and J2 must be distinct.
//
//    Input, int DIM_NUM, the dimension of the points X1 and X2.
//
//    Input, double X1[DIM_NUM], X2[DIM_NUM], the points that define
//    the line along which the equally spaced points are
//    generated, and which may or may not be included in the
//    set of computed points.
//
//    Output, double GRID_2N[DIM_NUM].  X(I) is the J-th point from the
//    set of equally spaced points.
//
{
  int i;
  double *x;

  x = new double[dim_num];

  if ( dim_num < 1 )
  {
    cerr << "\n";
    cerr << "GRID2N - Fatal error!\n";
    cerr << "  DIM_NUM < 1.\n";
    cerr << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  if ( j1 == j2 )
  {
    cerr << "\n";
    cerr << "GRID2N - Fatal error!\n";
    cerr << "  J1 = J2, leading to zero denominator.\n";
    cerr << "  J1 = " << j1 << "\n";
    cerr << "  J2 = " << j2 << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < dim_num; i++ )
  {
    x[i] = ( ( double ) ( j2 - j      ) * x1[i]   
           + ( double ) (      j - j1 ) * x2[i] )
           / ( double ) ( j2     - j1 );
  }

  return x;
}
//****************************************************************************80

double *grid3 ( int dim_num, int nstep1, int nstep2, double x1[], double x2[],
  double x3[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID3 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
//
//  Discussion:
//
//    The line between X1 and X2 will have NSTEP1 points generated along 
//    it, and the line between X1 and X3 will have NSTEP2 points generated
//    along it.
//
//    Fixing the second and third indices of X represents one point, with
//    the following special values:
//
//      X(*,1,1)      = X1
//      X(*,NSTEP1,1) = X2
//      X(*,1,NSTEP2) = X3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the points X1, X2 and X3.
//
//    Input, int NSTEP1, NSTEP2.  These are the number of
//    equally spaced points to generate in the first and second
//    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
//    X3 are always included in the set of points.
//
//    Input, double X1[DIM_NUM], X2[DIM_NUM], X3[DIM_NUM], the points
//    which define three corners of the parallelogram on
//    which the grid will be generated.
//
//    Output, double GRID3[DIM_NUM*NSTEP1*NSTEP2], the set of equally
//    spaced points.  
//
{
  int i;
  int j;
  int k;
  double psi1;
  double psi2;
  double psi3;
  double *x;

  if ( dim_num < 1 )
  {
    cerr << "\n";
    cerr << "GRID3 - Fatal error!\n";
    cerr << "  DIM_NUM < 1.\n";
    cerr << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  if ( nstep1 < 2 )
  {
    cerr << "\n";
    cerr << "GRID3 - Fatal error!\n";
    cerr << "  NSTEP1 < 2.\n";
    cerr << "  NSTEP1 = " << nstep1 << "\n";
    exit ( 1 );
  }

  if ( nstep2 < 2 )
  {
    cerr << "\n";
    cerr << "GRID3 - Fatal error!\n";
    cerr << "  NSTEP2 < 2.\n";
    cerr << "  NSTEP2 = " << nstep2 << "\n";
    exit ( 1 );
  }

  x = new double[dim_num*nstep1*nstep2];

  for ( j = 1; j <= nstep1; j++ )
  {
    psi2 = ( double ) ( j      - 1 ) 
         / ( double ) ( nstep1 - 1 );

    for ( k = 1; k <= nstep2; k++ )
    {
      psi3 = ( double ) (      k - 1 ) 
           / ( double ) ( nstep2 - 1 );

      psi1 = 1.0 - psi2 - psi3;

      for ( i = 0; i < dim_num; i++ )
      {
        x[i+(j-1)*dim_num+(k-1)*dim_num*nstep1] = 
            psi1 * x1[i]
          + psi2 * x2[i]
          + psi3 * x3[i];
      }
    }
  }

  return x;
}
//****************************************************************************80

double *grid3n ( int j, int k, int dim_num, int nstep1, int nstep2, 
  double x1[], double x2[], double x3[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID3N computes a parallelogram grid on 3 points in N dimensions.
//
//  Discussion:
//
//    The line between X1 and X2 will have NSTEP1
//    points generated along it, and the line between X1 and
//    X3 will have NSTEP2 points generated along it.
//
//    The following special values are:
//
//      J       K         X
//
//      1       1         X1
//      NSTEP1  1         X2
//      1       NSTEP2    X3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int J, K, the parallelogram coordinates
//    of the point.  J measures steps from X1 to X2, and
//    K measures steps from X1 to X3.  Normally, J would
//    be between 1 and NSTEP1, K between 1 and NSTEP2,
//    but this is not necessary.
//
//    Input, int DIM_NUM, the dimension of the points X1, X2 and X3.
//
//    Input, int NSTEP1, NSTEP2.  These are the number of
//    equally spaced points to generate in the first and second
//    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
//    X3 are always included in the set of points.
//
//    Input, double X1[DIM_NUM], X2[DIM_NUM], X3[DIM_NUM], the points
//    which define three corners of the parallelogram on
//    which the grid will be generated.
//
//    Output, double GRID3N[DIM_NUM], the point with coordinates (J,K)
//    from the the set of equally spaced points.  
//
{
  int i;
  double psi1;
  double psi2;
  double psi3;
  double *x;

  if ( dim_num < 1 )
  {
    cerr << "\n";
    cerr << "GRID3N - Fatal error!\n";
    cerr << "  DIM_NUM < 1.\n";
    cerr << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  if ( nstep1 < 2 )
  {
    cerr << "\n";
    cerr << "GRID3N - Fatal error!\n";
    cerr << "  NSTEP1 < 2.\n";
    cerr << "  NSTEP1 = " << nstep1 << "\n";
    exit ( 1 );
  }

  if ( nstep2 < 2 )
  {
    cerr << "\n";
    cerr << "GRID3N - Fatal error!\n";
    cerr << "  NSTEP2 < 2.\n";
    cerr << "  NSTEP2 = " << nstep2 << "\n";
    exit ( 1 );
  }

  x = new double[dim_num];

  psi2 = ( double ) ( j - 1  ) / ( double ) ( nstep1 - 1 );

  psi3 = ( double ) ( k - 1  ) / ( double ) ( nstep2 - 1 );

  psi1 = 1.0 - psi2 - psi3;

  for ( i = 0; i < dim_num; i++ )
  {
    x[i] = psi1 * x1[i] + psi2 * x2[i] + psi3 * x3[i];
  }

  return x;
}
//****************************************************************************80

double *grid4 ( int j1, int j2, int k1, int k2, int dim_num, int nstep1, 
  int nstep2, double x1[], double x2[], double x3[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID4 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
//
//  Discussion:
//
//    Unlike GRID3, GRID4 does not necessarily place X1 at the
//    "origin" of the parallelogram, with X2 and X3 set at the
//    extreme J and K coordinates.  Instead, the user is free
//    to specify the J and K coordinates of the points, although
//    they are required to lie on a subparallelogram of the
//    larger one.
//
//    The line through X1 and X2 will have NSTEP1
//    points generated along it, and the line through X1 and
//    X3 will have NSTEP2 points generated along it.
//
//    If we imagine that the
//    main parallelogram is drawn first, with coordinate
//    ranges 1 <= J <= NSTEP1 and 1 <= K <= NSTEP2, then
//    these indices determine the (J,K) coordinates of the
//    three points, namely:
//
//      X1 : (J1,K1)
//      X2 : (J2,K1)
//      X3 : (J1,K2)
//
//    Of course, we actually start with the points X1, X2,
//    and X3, and they define a parallelogram and a (J,K)
//    coordinate system over the plane containing them.  We
//    then are free to consider the parallelogram defined
//    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
//    which may or may not contain any of the points X1, X2
//    and X3.
//
//    Assuming that the indices J1, J2, K1 and K2 are "within
//    bounds", the following special values will be computed:
//
//      X(*,J1,K1) = X1
//      X(*,J2,K1) = X2
//      X(*,J1,K2) = X3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int J1, J2, K1, K2, the indices.  
//
//    Input, int DIM_NUM, the dimension of the points X1, X2 and X3.
//
//    Input, int NSTEP1, NSTEP2.  These are the number of
//    equally spaced points to generate in the first and second
//    directions.  NSTEP1 and NSTEP2 should be at least 1.
//
//    Input, double X1[DIM_NUM], X2[DIM_NUM], X3[DIM_NUM], the points
//    which define three corners of the parallelogram on
//    which the grid will be generated.
//
//    Output, double X[DIM_NUM*NSTEP1*NSTEP2], the set of equally
//    spaced points.  Fixing the second and third indices
//    of X represents one point.
//
{
  int i;
  int j;
  int k;
  double psi1;
  double psi2;
  double psi3;
  double *x;

  if ( dim_num < 1 )
  {
    cerr << "\n";
    cerr << "GRID4 - Fatal error!\n";
    cerr << "  DIM_NUM < 1.\n";
    cerr << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  if ( nstep1 < 2 )
  {
    cerr << "\n";
    cerr << "GRID4 - Fatal error!\n";
    cerr << "  NSTEP1 < 2.\n";
    cerr << "  NSTEP1 = " << nstep1 << "\n";
    exit ( 1 );
  }

  if ( nstep2 < 2 )
  {
    cerr << "\n";
    cerr << "GRID4 - Fatal error!\n";
    cerr << "  NSTEP2 < 2.\n";
    cerr << "  NSTEP2 = " << nstep2 << "\n";
    exit ( 1 );
  }

  if ( j1 == j2 )
  {
    cerr << "\n";
    cerr << "GRID4 - Fatal error!\n";
    cerr << "  J1 = J2, leading to zero denominator.\n";
    cerr << "  J1 = " << j1 << "\n";
    cerr << "  J2 = " << j2 << "\n";
    exit ( 1 );
  }

  if ( k1 == k2 )
  {
    cerr << "\n";
    cerr << "GRID4 - Fatal error!\n";
    cerr << "  K1 = K2, leading to zero denominator.\n";
    cerr << "  K1 = " << k1 << "\n";
    cerr << "  K2 = " << k2 << "\n";
    exit ( 1 );
  }

  x = new double[dim_num*nstep1*nstep2];

  for ( j = 1; j <= nstep1; j++ )
  {
    psi2 = ( double ) (  j - j1 ) 
         / ( double ) ( j2 - j1 );

    for ( k = 1; k <= nstep2; k++ )
    {
      psi3 = ( double ) (  k - k1 ) 
           / ( double ) ( k2 - k1 );

      psi1 = 1.0 - psi2 - psi3;

      for ( i = 0; i < dim_num; i++ )
      {
        x[i+(j-1)*dim_num+(k-1)*dim_num*nstep1] =
            psi1 * x1[i]
          + psi2 * x2[i]
          + psi3 * x3[i];
      }
    }
  }

  return x;
}
//****************************************************************************80

double *grid4n ( int j, int j1, int j2, int k, int k1, int k2, int dim_num, 
  int nstep1, int nstep2, double x1[], double x2[], double x3[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID4N computes a single point on a parallelogram grid in N space.
//
//  Discussion:
//
//    The computation is identical to that of GRID4, except that
//    only one point at a time is computed.
//
//    The line through X1 and X2 will have NSTEP1
//    points generated along it, and the line through X1 and
//    X3 will have NSTEP2 points generated along it.
//
//    The following special values will be computed:
//
//      J  K  X
//
//      J1 K1 X1
//      J2 K2 X2
//      J1 K2 X3
//
//    If we imagine that the main parallelogram is drawn first, with 
//    coordinate ranges 1 <= J <= NSTEP1 and 1 <= K <= NSTEP2, then
//    the indices J and K determine the (J,K) coordinates of the
//    three points X1, X2, and X3, namely:
//
//      X1 : (J1,K1)
//      X2 : (J2,K1)
//      X3 : (J1,K2)
//
//    Of course, we actually start with the points X1, X2,
//    and X3, and they define a parallelogram and an (J,K)
//    coordinate system over the plane containing them.  We
//    then are free to consider the parallelogram defined
//    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
//    which may or may not contain any of the points X1, X2
//    and X3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int J, the J coordinate of the point X.
//
//    Input, int J1, J2.  See discussion.
//
//    Input, int K, the K coordinate of the point X.
//
//    Input, int K1, K2.  See discussion.
//
//    Input, int DIM_NUM, the dimension of the points X, X1, X2 and X3.
//
//    Input, int NSTEP1, NSTEP2.  These are the number of
//    equally spaced points generated in the first and second
//    directions.
//    NSTEP1 and NSTEP2 should be at least 1.
//
//    Input, double X1[DIM_NUM], X2[DIM_NUM], X3[DIM_NUM], the points
//    which define three corners of the parallelogram on
//    which the grid will be generated.
//
//    Output, double GRID4N[DIM_NUM], the point whose parallelogram
//    coordinates are (J,K).
//
{
  int i;
  double psi1;
  double psi2;
  double psi3;
  double *x;

  if ( dim_num < 1 )
  {
    cerr << "\n";
    cerr << "GRID4N - Fatal error!\n";
    cerr << "  DIM_NUM < 1.\n";
    cerr << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  if ( nstep1 < 2 )
  {
    cerr << "\n";
    cerr << "GRID4N - Fatal error!\n";
    cerr << "  NSTEP1 < 2.\n";
    cerr << "  NSTEP1 = " << nstep1 << "\n";
    exit ( 1 );
  }

  if ( nstep2 < 2 )
  {
    cerr << "\n";
    cerr << "GRID4N - Fatal error!\n";
    cerr << "  NSTEP2 < 2.\n";
    cerr << "  NSTEP2 = " << nstep2 << "\n";
    exit ( 1 );
  }

  if ( j1 == j2 )
  {
    cerr << "\n";
    cerr << "GRID4N - Fatal error!\n";
    cerr << "  J1 = J2, leading to zero denominator.\n";
    cerr << "  J1 = " << j1 << "\n";
    cerr << "  J2 = " << j2 << "\n";
    exit ( 1 );
  }

  if ( k1 == k2 )
  {
    cerr << "\n";
    cerr << "GRID4N - Fatal error!\n";
    cerr << "  K1 = K2, leading to zero denominator.\n";
    cerr << "  K1 = " << k1 << "\n";
    cerr << "  K2 = " << k2 << "\n";
    exit ( 1 );
  }

  psi2 = ( double ) ( j  - j1 ) / ( double ) ( j2 - j1 );

  psi3 = ( double ) ( k  - k1 ) / ( double ) ( k2 - k1 );

  psi1 = 1.0 - psi2 - psi3;

  x = new double[dim_num];

  for ( i = 0; i < dim_num; i++ )
  {
    x[i] = psi1 * x1[i] + psi2 * x2[i] + psi3 * x3[i];
  }

  return x;
}
//****************************************************************************80

short int i2_reverse_bytes ( short int x )

//****************************************************************************80
//
//  Purpose:
//
//    I2_REVERSE_BYTES reverses the two bytes in an I2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, short int X, a value whose bytes are to be reversed.
//
//    Output, short int I2_REVERSE_BYTES, a value with
//    bytes in reverse order.
//
{
  char c;
  union 
  {
    short int yshortint;
    char ychar[2];
  } y;

  y.yshortint = x;
  
  c = y.ychar[0];
  y.ychar[0] = y.ychar[1];
  y.ychar[1] = c;

  return ( y.yshortint );
}
//****************************************************************************80

int i4_log_10 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
//
//  Example:
//
//        I  I4_LOG_10
//    -----  --------
//        0    0
//        1    0
//        2    0
//        9    0
//       10    1
//       11    1
//       99    1
//      100    2
//      101    2
//      999    2
//     1000    3
//     1001    3
//     9999    3
//    10000    4
//
//  Discussion:
//
//    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number whose logarithm base 10 is desired.
//
//    Output, int I4_LOG_10, the integer part of the logarithm base 10 of
//    the absolute value of X.
//
{
  int i_abs;
  int ten_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    ten_pow = 10;

    i_abs = abs ( i );

    while ( ten_pow <= i_abs )
    {
      value = value + 1;
      ten_pow = ten_pow * 10;
    }

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

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If 
//      NREM = I4_MODP ( I, J ) 
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
// 
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is 
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
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
//    Output, int I4_POWER, the value of I^J.
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
      cerr << "I4_POWER - Fatal error!\n";
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
      cerr << "I4_POWER - Fatal error!\n";
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

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
  int value;

  if ( i < 0 ) 
  {
    value = -1;
  }
  else
  {
    value = 1;
  }
  return value;
}
//****************************************************************************80

void i4_swap ( int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = *i;
  *i = *j;
  *j = k;
 
  return;
}
//****************************************************************************80

int *i4_to_digits_decimal ( int i, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_DIGITS_DECIMAL determines the last N decimal digits of an I4.
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
//    Input, int I, the integer to be analyzed.
//
//    Input, int N, the number of digits to determine.
//
//    Output, int I4_TO_DIGITS_DECIMAL[N], the last N decimal digits of I.
//    DIGIT[I-1] is the "coefficient" of 10**(I-1).
//
{
  int *digit;
  int j;

  digit = new int[n];

  i = abs ( i );

  for ( j = 1; j <= n; j++ )
  {
    digit[j-1] = i % 10;
    i = ( i - digit[j-1] ) / 10;
  }

  return digit;
}
//****************************************************************************80

int *i4_to_fac ( int i, int prime_num )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_FAC converts an I4 into a product of prime factors.
//
//  Discussion:
//
//    This routine will fail if the input integer is not positive,
//    or if PRIME_NUM is too small to account for the factors of the integer.
//
//    The formula is:
//
//      I = Product ( 1 <= J <= PRIME_NUM ) PRIME(J)**NPOWER(J).
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
//    Input, int I, the integer to be factored.
//
//    Input, int PRIME_NUM, the number of prime factors for
//    which storage has been allocated.
//
//    Output, int I4_TO_FAC[PRIME_NUM], the powers of the primes.
//
{
  int j;
  int *npower;
  int p;

  if ( i <= 0 )
  {
    cerr << "\n";
    cerr << "I4_TO_FAC - Fatal error!\n";
    cerr << "  Input integer I is not positive.\n";
    exit ( 1 );
  }

  npower = new int[prime_num];
//
//  Try dividing the remainder by each prime.
//
  for ( j = 1; j <= prime_num; j++ )
  {
    npower[j-1] = 0;

    p = prime ( j );

    while ( ( i % p ) == 0 )
    {
      npower[j-1] = npower[j-1] + 1;
      i = i / p;
    }

  }
  return npower;
}
//****************************************************************************80

char i4_to_isbn ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_ISBN converts an I4 to an ISBN digit.
//
//  Discussion:
//
//    Only the integers 0 through 10 can be input.  The representation
//    of 10 is 'X'.
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
//  Reference:
//
//    Book Industry Study Group,
//    The Evolution in Product Identification:
//    Sunrise 2005 and the ISBN-13,
//    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
//
//  Parameters:
//
//    Input, int I, an integer between 0 and 10.
//
//    Output, char I4_TO_ISBN, the ISBN character code of the integer.
//    If I is illegal, then I4_TO_ISBN is set to '?'.
//
{
       if ( i == 0 )
  {
    return '0';
  }
  else if ( i == 1 )
  {
    return '1';
  }
  else if ( i == 2 )
  {
    return '2';
  }
  else if ( i == 3 )
  {
    return '3';
  }
  else if ( i == 4 )
  {
    return '4';
  }
  else if ( i == 5 )
  {
    return '5';
  }
  else if ( i == 6 )
  {
    return '6';
  }
  else if ( i == 7 )
  {
    return '7';
  }
  else if ( i == 8 )
  {
    return '8';
  }
  else if ( i == 9 )
  {
    return '9';
  }
  else if ( i == 10 )
  {
    return 'X';
  }
  else
  {
    return '?';
  }
}
//****************************************************************************80

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2006
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
//    in Handbook of Simulation,
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
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 ) 
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

double i4int_to_r8int ( int imin, int imax, int i, double rmin, double rmax )

//****************************************************************************80
//
//  Purpose:
//
//    I4INT_TO_R8INT maps an I4 interval to an R8 interval.
//
//  Discussion:
//
//    The formula is
//
//      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IMIN, IMAX, the range.
//
//    Input, int I, the integer to be converted.
//
//    Input, double RMIN, RMAX, the range.
//
//    Output, double R, the corresponding value in [RMIN,RMAX].
//
{
  double r;

  if ( imax == imin )
  {
    r = 0.5 * ( rmin + rmax );
  }
  else
  {
    r = ( ( double ) ( imax - i        ) * rmin   
        + ( double ) (        i - imin ) * rmax ) 
        / ( double ) ( imax     - imin );
  }

  return r;
}
//****************************************************************************80

void i4vec_copy ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY copies an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int A2[N], the copy of A1.
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

int *i4vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR_NEW sets an I4VEC to the indicator vector.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int I4VEC_INDICATOR_NEW[N], the array.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }
  return a;
}
//****************************************************************************80

int i4vec_min ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MIN returns the value of the minimum element in an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MIN, the value of the minimum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] < value )
    {
      value = a[i];
    }
  }
  return value; 
}
//****************************************************************************80

void i4vec_permute ( int n, int p[], int base, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PERMUTE permutes an I4VEC in place.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    This routine permutes an array of integer "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   1,   3,   4,   0,   2 )
//      A = (   1,   2,   3,   4,   5 )
//
//    Output:
//
//      A    = (   2,   4,   5,   1,   3 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.  
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for 
//    a 1-based permutation.
//
//    Input/output, int A[N], the array to be permuted.
//
{
  int a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "I4VEC_PERMUTE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
//  So temporarily add 1-BASE to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "I4VEC_PERMUTE - Fatal error!\n";
          cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
//
//  Restore the base of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 + base;
  }

  return;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
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
         << ": " << setw(8) << a[i]  << "\n";
  }
  return;
}
//****************************************************************************80

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

void i4vec_zero ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}
//****************************************************************************80

int *i4vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO_NEW creates and zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
}
//****************************************************************************80

void ij_next ( int *i, int *j, int n )

//****************************************************************************80
//
//  Purpose:
//
//    IJ_NEXT returns the next matrix index.
//
//  Discussion:
//
//    For N = 3, the sequence of indices returned is:
//
//      (1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (3,1), (3,2), (3,3), (0,0).
//
//    Note that once the value (N,N) is returned, the next value returned
//    will be (0,0).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On input, the current pair of indices.
//    On output, the next pair of indices.  If either index is illegal on
//    input, the output value of (I,J) will be (1,1).
//
//    Input, int N, the maximum value for I and J.
//
{
  if ( n < 1 )
  {
    *i = 0;
    *j = 0;
    return;
  }

  if ( *i < 1 || n < *i || *j < 1 || n < *j )
  {
    *i = 1;
    *j = 1;
    return;
  }

  if ( *j < n )
  {
    *j = *j + 1;
  }
  else if ( *i < n )
  {
    *i = *i + 1;
    *j = 1;
  }
  else
  {
    *i = 0;
    *j = 0;
  }

  return;
}
//****************************************************************************80

void ij_next_gt ( int *i, int *j, int n )

//****************************************************************************80
//
//  Purpose:
//
//    IJ_NEXT_GT returns the next matrix index, with the constraint that I < J.
//
//  Discussion:
//
//    For N = 5, the sequence of indices returned is:
//
//      (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On input, the current pair of indices.
//    On output, the next pair of indices.  If either index is illegal on
//    input, the output value of (I,J) will be (1,2).
//
//    Input, int N, the maximum value for I and J.
//    A value of N less than 2 is nonsense.
//
{
  if ( n < 2 )
  {
    *i = 0;
    *j = 0;
    return;
  }

  if ( *i < 1 || n < *i || *j < 1 || n < *j || *j <= *i )
  {
    *i = 1;
    *j = 2;
    return;
  }

  if ( *j < n )
  {
    *j = *j + 1;
  }
  else if ( *i < n - 1 )
  {
    *i = *i + 1;
    *j = *i + 1;
  }
  else
  {
    *i = 0;
    *j = 0;
  }

  return;
}
//****************************************************************************80

void index_box2_next_2d ( int n1, int n2, int ic, int jc, int *i, int *j, 
  int *more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
//
//  Discussion:
//
//    The box has center at (IC,JC), and has half-widths N1 and N2.
//    The indices are exactly those which are between (IC-N1,JC-N2) and
//    (IC+N1,JC+N2) with the property that at least one of I and J
//    is an "extreme" value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the half-widths of the box, that is, the
//    maximum distance allowed between (IC,JC) and (I,J).
//
//    Input, int IC, JC, the central cell of the box.
//
//    Input/output, int *I, *J.  On input, the previous index set.
//    On output, the next index set.  On the first call, MORE should
//    be set to FALSE, and the input values of I and J are ignored.
//
//    Input/output, bool *MORE.
//    On the first call for a given box, the user should set MORE to FALSE.
//    On return, the routine sets MORE to TRUE.
//    When there are no more indices, the routine sets MORE to FALSE.
//
{
  if ( !(*more) )
  {
    *more = true;
    *i = ic - n1;
    *j = jc - n2;
    return;
  }

  if ( *i == ic + n1 && 
       *j == jc + n2 )
  {
    *more = false;
    return;
  }
//
//  Increment J.
//
  *j = *j + 1;
//
//  Check J.
//
  if ( jc + n2 < *j )
  {
    *j = jc - n2;
    *i = *i + 1;
  }
  else if ( *j < jc + n2 && ( *i == ic - n1 || *i == ic + n1 ) )
  {
    return;
  }
  else
  {
    *j = jc + n2;
    return;
  }

  return;
}
//****************************************************************************80

void index_box2_next_3d ( int n1, int n2, int n3, int ic, int jc, int kc, 
  int *i, int *j, int *k, bool *more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
//
//  Discussion:
//
//    The box has a central cell of (IC,JC,KC), with a half widths of
//    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3) and
//    (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J, and K
//    is an "extreme" value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the "half widths" of the box, that is, the
//    maximum distances from the central cell allowed for I, J and K.
//
//    Input, int IC, JC, KC, the central cell of the box.
//
//    Input/output, int *I, *J, *K.  On input, the previous index set.
//    On output, the next index set.  On the first call, MORE should
//    be set to FALSE, and the input values of I, J, and K are ignored.
//
//    Input/output, bool *MORE.
//    On the first call for a given box, the user should set MORE to FALSE.
//    On return, the routine sets MORE to TRUE.
//    When there are no more indices, the routine sets MORE to FALSE.
//
{
  if ( !(*more) )
  {
    *more = true;
    *i = ic - n1;
    *j = jc - n2;
    *k = kc - n3;
    return;
  }

  if ( *i == ic + n1 && 
       *j == jc + n2 && 
       *k == kc + n3 )
  {
    *more = false;
    return;
  }
//
//  Increment K.
//
  *k = *k + 1;
//
//  Check K.
//
  if ( kc + n3 < *k )
  {
    *k = kc - n3;
    *j = *j + 1;
  }
  else if ( *k < kc + n3 &&
    ( *i == ic - n1 || 
      *i == ic + n1 || 
      *j == jc - n2 || 
      *j == jc + n2 ) )
  {
    return;
  }
  else
  {
    *k = kc + n3;
    return;
  }
//
//  Check J.
//
  if ( jc + n2 < *j )
  {
    *j = jc - n2;
    *i = *i + 1;
  }
  else if ( *j < jc + n2 &&
    ( *i == ic - n1 || 
      *i == ic + n1 || 
      *k == kc - n3 || 
      *k == kc + n3 ) )
  {
    return;
  }
  else
  {
    *j = jc + n2;
    return;
  }

  return;
}
//****************************************************************************80

int index1_col ( int i_min, int i, int i_max, int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX1_COL indexes a 1D vector by columns.
//
//  Discussion:
//
//    This 1D routine is provided primarily for analogy.
//    Moreover, in 1D there is no difference between row and column indexing.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for the first index,
//    the minimum, the index, and the maximum.
//
//    Input, int INDEX_MIN, the index of element I_MIN.
//    Typically, this is 0 or 1.
//
//    Output, int INDEX1_COL, the index of element I.
//
{
  int value;

  value = index_min + ( i - i_min );

  return value;
}
//****************************************************************************80

int index1_row ( int i_min, int i, int i_max, int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX1_ROW indexes a 1D vector by rows.
//
//  Discussion:
//
//    This 1D routine is provided primarily for analogy.
//    Moreover, in 1D there is no difference between row and column indexing.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for the first index,
//    the minimum, the index, and the maximum.
//
//    Input, int INDEX_MIN, the index of element I_MIN.
//    Typically, this is 0 or 1.
//
//    Output, int INDEX1_ROW, the index of element I.
//
{
  int value;

  value = index_min + ( i - i_min );

  return value;
}
//****************************************************************************80

int index2_col ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX2_COL indexes a 2D array by columns.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
//    and increasing the row index first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Input, int INDEX_MIN, the index of element (I_MIN,J_MIN).
//    Typically, this is 0 or 1.
//
//    Output, int INDEX2_COL, the index of element (I,J).
//
{
  int value;

  value = index_min + ( i - i_min ) + ( j - j_min ) * ( i_max + 1 - i_min );

  return value;
}
//****************************************************************************80

int index2_row ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX2_ROW indexes a 2D array by row.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
//    and increasing the column index first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Input, int INDEX_MIN, the index of element (I_MIN,J_MIN).
//    Typically, this is 0 or 1.
//
//    Output, int INDEX2_ROW, the index of element (I,J).
//
{
  int value;

  value = index_min + ( j - j_min ) + ( i - i_min ) * ( j_max + 1 - j_min );

  return value;
}
//****************************************************************************80

int index3_col ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int k_min, int k, int k_max, int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX3_COL indexes a 3D array by columns.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry (I_MIN,J_MIN,K_MIN), 
//    and increasing the row index first, then the column index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Input, int K_MIN, K, K_MAX, for plane indices,
//    the minimum, the index, and the maximum.
//
//    Input, int INDEX_MIN, the index of (I_MIN,J_MIN,K_MIN).
//    Typically, this is 0 or 1.
//
//    Output, int INDEX3_COL, the index of element (I,J,K).
//
{
  int value;

  value = index_min 
             + ( i - i_min ) 
             + ( j - j_min ) * ( i_max + 1 - i_min ) 
             + ( k - k_min ) * ( j_max + 1 - j_min ) * ( i_max + 1 - i_min );

  return value;
}
//****************************************************************************80

int index3_row ( int i_min, int i, int i_max, int j_min, int j, int j_max, 
  int k_min, int k, int k_max, int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX3_ROW indexes a 3D array by rows.
//
//  Discussion:
//
//    When we say "by rows", we really just mean that entries of the array are 
//    indexed starting at entry (I_MIN,J_MIN,K_MIN), and the increasing the LAST
//    index first, then the next-to-the-last, and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I_MIN, I, I_MAX, for row indices,
//    the minimum, the index, and the maximum.
//
//    Input, int J_MIN, J, J_MAX, for column indices,
//    the minimum, the index, and the maximum.
//
//    Input, int K_MIN, K, K_MAX, for plane indices,
//    the minimum, the index, and the maximum.
//
//    Input, int INDEX_MIN, the index of (I_MIN,J_MIN,K_MIN).
//    Typically, this is 0 or 1.
//
//    Output, int INDEX3_ROW, the index of element (I,J,K).
//
{
  int value;

  value = index_min 
             + ( k - k_min ) 
             + ( j - j_min ) * ( k_max + 1 - k_min ) 
             + ( i - i_min ) * ( j_max + 1 - j_min ) * ( k_max + 1 - k_min );

  return value;
}
//****************************************************************************80

int index4_col ( int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max, 
  int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max, 
  int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX4_COL indexes a 4D array by columns.
//
//  Discussion:
//
//    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
//    and increasing the initial index first, then the second, third and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1_MIN, I1, I1_MAX, for index 1,
//    the minimum, the index, and the maximum.
//
//    Input, int I2_MIN, I2, I2_MAX, for index 2,
//    the minimum, the index, and the maximum.
//
//    Input, int I3_MIN, I3, I3_MAX, for index 3,
//    the minimum, the index, and the maximum.
//
//    Input, int I4_MIN, I4, I4_MAX, for index 4,
//    the minimum, the index, and the maximum.
//
//    Input, int INDEX_MIN, the index of 
//    (I1_MIN,I2_MIN,I3_MIN,I4_MIN).  Typically, this is 0 or 1.
//
//    Output, int INDEX4_COL, the index of element (I1,I2,I3,I4).
//
{
  int value;

  value = index_min 
    + ( i1 - i1_min ) 
    + ( i2 - i2_min ) * ( i1_max + 1 - i1_min ) 
    + ( i3 - i3_min ) * ( i2_max + 1 - i2_min ) * ( i1_max + 1 - i1_min ) 
    + ( i4 - i4_min ) * ( i3_max + 1 - i3_min ) * ( i2_max + 1 - i2_min ) 
      * ( i1_max + 1 - i1_min );

  return value;
}
//****************************************************************************80

int index4_row ( int i1_min, int i1, int i1_max, int i2_min, int i2, int i2_max, 
  int i3_min, int i3, int i3_max, int i4_min, int i4, int i4_max, 
  int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX4_ROW indexes a 4D array by rows.
//
//  Discussion:
//
//    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
//    and increasing the last index, then the next to last, and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1_MIN, I1, I1_MAX, for index 1,
//    the minimum, the index, and the maximum.
//
//    Input, int I2_MIN, I2, I2_MAX, for index 2,
//    the minimum, the index, and the maximum.
//
//    Input, int I3_MIN, I3, I3_MAX, for index 3,
//    the minimum, the index, and the maximum.
//
//    Input, int I4_MIN, I4, I4_MAX, for index 4,
//    the minimum, the index, and the maximum.
//
//    Input, int INDEX_MIN, the index of element 
//    (I1_MIN,I2_MIN,I3_MIN,I4_MIN).  Typically, this is 0 or 1.
//
//    Output, int INDEX4_ROW, the index of element (I1,I2,I3,I4).
//
{
  int value;

  value = index_min 
    + ( i4 - i4_min ) 
    + ( i3 - i3_min ) * ( i4_max + 1 - i4_min ) 
    + ( i2 - i2_min ) * ( i3_max + 1 - i3_min ) * ( i4_max + 1 - i4_min ) 
    + ( i1 - i1_min ) * ( i2_max + 1 - i2_min ) * ( i3_max + 1 - i3_min ) 
      * ( i4_max + 1 - i4_min );

  return value;
}
//****************************************************************************80

int indexn_col ( int n, int i_min[], int i[], int i_max[], int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEXN_COL indexes an ND array by columns.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry 
//      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
//    and increasing the first index up to I_MAX(1), 
//    then the second and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of indices.
//
//    Input, int I_MIN[N], the minimum indices.
//
//    Input, int I[N], the indices.
//
//    Input, int I_MAX[N], for maximum indices.
//
//    Input, int INDEX_MIN, the index of 
//    ( I_MIN[0], I_MIN[1],...,I_MIN[N-1] ).  Typically, this is 0 or 1.
//
//    Output, int INDEXN_COL, the index of element I.
//
{
  int j;
  int value;

  value = ( i[n-1] - i_min[n-1] );

  for ( j = n - 2; 0 <= j; j-- )
  {
    value = value * ( i_max[j] + 1 - i_min[j] ) + ( i[j] - i_min[j] );
  }
  value = value + index_min;

  return value;
}
//****************************************************************************80

int indexn_row ( int n, int i_min[], int i[], int i_max[], int index_min )

//****************************************************************************80
//
//  Purpose:
//
//    INDEXN_ROW indexes an ND array by rows.
//
//  Discussion:
//
//    Entries of the array are indexed starting at entry 
//      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
//    and increasing the last index up to I_MAX(N), 
//    then the next-to-last and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of indices.
//
//    Input, int I_MIN[N], the minimum indices.
//
//    Input, int I[N], the indices.
//
//    Input, int I_MAX[N], for maximum indices.
//
//    Input, int INDEX_MIN, the index of 
//    ( I_MIN[0], I_MIN[1],...,I_MIN[N-1] ).  Typically, this is 0 or 1.
//
//    Output, int INDEXN_ROW, the index of element I.
//
{
  int j;
  int value;

  value = ( i[0] - i_min[0] );

  for ( j = 1; j < n; j++ )
  {
    value = value * ( i_max[j] + 1 - i_min[j] ) + ( i[j] - i_min[j] );
  }
  value = value + index_min;

  return value;
}
//****************************************************************************80

int isbn_check ( string isbn )

//****************************************************************************80
//
//  Purpose:
//
//    ISBN_CHECK checks an ISBN code.
//
//  Discussion:
//
//    ISBN stands for International Standard Book Number.  A unique ISBN
//    is assigned to each new book.  The ISBN includes 10 digits.  There is
//    an initial digit, then a dash, then a set of digits which are a
//    code for the publisher, another digit, and then the check digit:
//
//      initial-publisher-book-check
//
//    The number of digits used for the publisher and book codes can vary,
//    but the check digit is always one digit, and the total number of
//    digits is always 10.
//
//    The check digit is interesting, because it is a way of trying to
//    make sure that an ISBN has not been incorrectly copied.  Specifically,
//    if the ISBN is correct, then its ten digits will satisfy
//
//       10 * A + 9 * B + 8 * C + 7 * D + 6 * E
//      + 5 * F * 4 * G * 3 * H + 2 * I +     J  = 0 mod 11.
//
//    Here, we've taken 'A' to represent the first digit and 'J' the
//    last (which is the check digit).  In order for the code to work,
//    the value of J must be allowed to be anything from 0 to 10.  In
//    cases where J works out to be 10, the special digit 'X' is used.
//    An 'X' digit can only occur in this last check-digit position
//    on an ISBN.
//
//  Example:
//
//    0-8493-9640-9
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Book Industry Study Group,
//    The Evolution in Product Identification:
//    Sunrise 2005 and the ISBN-13,
//    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
//
//  Parameters:
//
//    Input, string ISBN, an ISBN code.
//
//    Output, int ISBN_CHECK, the value of the ISBN check sum.
//    If CHECK is zero, the ISBN code is legitimate.
//    If CHECK is -1, then the ISBN code is not legitimate because it does
//    not contain exactly 10 digits.  If CHECK is between 1 and 10, then
//    the ISBN code has the right number of digits, but at least one of
//    the digits is incorrect.
//
{
  char c;
  int check;
  int digit[10];
  int i;
  int lenc;
  int num_digit;
//
//  Determine how many digits have been supplied.
//
  lenc = s_len_trim ( isbn );

  i = 0;
  num_digit = 0;

  for ( ; ; )
  {
    i = i + 1;

    if ( lenc < i )
    {
      break;
    }

    c = isbn[i-1];

    if ( ch_is_digit ( c ) )
    {
      digit[num_digit] = isbn_to_i4 ( c );
      num_digit = num_digit + 1;
    }
    else if ( ( num_digit == 9 && isbn[i-1] == 'X' ) ||
              ( num_digit == 9 && isbn[i-1] == 'x' ) )
    {
      digit[num_digit] = isbn_to_i4 ( c );
      num_digit = num_digit + 1;
    }

    if ( 10 <= num_digit )
    {
      break;
    }
  }
//
//  If we didn't get exactly 10 digits, return with an error.
//
  if ( num_digit != 10 )
  {
    check = -1;
    return check;
  }
//
//  Compute the checksum.
//
  check = 0;
  for ( i = 0; i < 10; i++ )
  {
    check = check + ( 11 - i ) * digit[i];
  }

  check = ( check % 11 );

  return check;
}
//****************************************************************************80

void isbn_fill ( string isbn )

//****************************************************************************80
//
//  Purpose:
//
//    ISBN_FILL fills in a missing digit in an ISBN code.
//
//  Example:
//
//    Input:
//
//      0-8493-9?40-9
//
//    Output:
//
//      0-8493-9640-9
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Book Industry Study Group,
//    The Evolution in Product Identification:
//    Sunrise 2005 and the ISBN-13,
//    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
//
//  Parameters:
//
//    Input/output, string ISBN, a partial ISBN code.  On input,
//    a single digit has been replaced by the character '?', signifying
//    that that digit is missing.  The routine replaces the question
//    mark by the correct digit.
//
{
  char c;
  int check;
  int digit[10];
  int digit_pos;
  int i;
  int isbn_pos;
  int j;
  int k;
  int lenc;
  int num_digit;

  lenc = s_len_trim ( isbn );

  i = 0;
  isbn_pos = -1;
  digit_pos = -1;
  num_digit = 0;

  for ( ; ; )
  {
    i = i + 1;

    if ( lenc < i )
    {
      break;
    }

    c = isbn[i-1];

    if ( ch_is_digit ( c ) )
    {
      num_digit = num_digit + 1;
      digit[num_digit-1] = isbn_to_i4 ( c );
    }
    else if ( ( num_digit == 9 && isbn[i-1] == 'X' ) || 
              ( num_digit == 9 && isbn[i-1] == 'x' ) )
    {
      num_digit = num_digit + 1;
      digit[num_digit-1] = isbn_to_i4 ( c );
    }
    else if ( c == '?' )
    {
      if ( isbn_pos == -1 )
      {
        num_digit = num_digit + 1;
        digit[num_digit-1] = 0;
        digit_pos = num_digit;
        isbn_pos = i;
      }
      else
      {
        cerr << "\n";
        cerr << "ISBN_FILL - Fatal error!\n";
        cerr << "  Only one question mark is allowed.\n";
        exit ( 1 );
      }
    }

    if ( 10 <= num_digit )
    {
      break;
    }
  }

  if ( num_digit != 10 )
  {
    cerr << "\n";
    cerr << "ISBN_FILL - Fatal error!\n";
    cerr << "  The input ISBN code did not have 10 digits.\n";
    exit ( 1 );
  }

  if ( isbn_pos == -1 )
  {
    cerr << "\n";
    cerr << "ISBN_FILL - Fatal error!\n";
    cerr << "  A question mark is required.\n";
    exit ( 1 );
  }

  check = 0;
  for ( i = 1; i <= 10; i++ )
  {
    check = check + ( 11 - i ) * digit[i-1];
  }

  check = ( check % 11 );

  if ( check == 0 )
  {
    k = 0;
  }
//
//  Need to solve the modular equation:
//
//    A * X = B mod C
//
//  Below is a stupid way.  One day I will come back and fix this up.
//
  else
  {
    for ( i = 1; i <= 10; i++ )
    {
      j = ( 11 - digit_pos ) * i + check;
      if ( ( j % 11 ) == 0 )
      {
        k = i;
      }
    }
  }

  isbn[isbn_pos-1] = i4_to_isbn ( k );

  return;
}
//****************************************************************************80

int isbn_to_i4 ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    ISBN_TO_I4 converts an ISBN character into an I4.
//
//  Discussion:
//
//    The characters '0' through '9' stand for themselves, but
//    the character 'X' or 'x' stands for 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Book Industry Study Group,
//    The Evolution in Product Identification:
//    Sunrise 2005 and the ISBN-13,
//    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
//
//  Parameters:
//
//    Input, char C, the ISBN character code to be converted.
//
//    Output, int ISBN_TO_I4, the numeric value of the character
//    code, between 0 and 10.  This value is returned as -1 if C is
//    not a valid character code.
//
{
  int value;

  if ( '0' <= c && c <= '9' )
  {
    value = c - '0';
  }
  else if ( c == 'X' || c == 'x' )
  {
    value = 10;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int iset2_compare ( int x1, int y1, int x2, int y2 )

//****************************************************************************80
//
//  Purpose:
//
//    ISET2_COMPARE compares two I2 sets.
//
//  Discussion:
//
//    The I2 set (X1,Y1) < (X2,Y2) if
//
//      min ( X1, Y1 ) < min ( X2, Y2 ) or
//      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) < max ( X2, Y2 )
//
//    The I2 set (X1,Y1) = (X2,Y2) if
//
//      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) = max ( X2, Y2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X1, Y1, the first I2 set.
//
//    Input, int X2, Y2, the second I2 set.
//
//    Output, int ISET2_COMPARE: 
//    -1, (X1,Y1) < (X2,Y2).
//     0, (X1,Y1) = (X2,Y2).
//    +1, (X1,Y1) > (X2,Y2).
//
{
  int a1;
  int a2;
  int b1;
  int b2;
  int value;

  a1 = i4_min ( x1, y1 );
  b1 = i4_max ( x1, y1 );

  a2 = i4_min ( x2, y2 );
  b2 = i4_max ( x2, y2 );

  if ( a1 < a2 )
  {
    value = -1;
  }
  else if ( a2 < a1 )
  {
    value = +1;
  }
  else if ( b1 < b2 )
  {
    value = -1;
  }
  else if ( b2 < b1 )
  {
    value = +1;
  }
  else
  {
    value = 0;
  }
  return value;
}
//****************************************************************************80

int lcm_12n ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    LCM_12N computes the least common multiple of the integers 1 through N.
//
//  Example:
//
//    N    LCM_12N
//
//    1          1
//    2          2
//    3          3
//    4         12
//    5         60
//    6         60
//    7        420
//    8        840
//    9       2520
//   10       2520
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the value of N.
//
//    Output, int LCM_12N, the least common multiple of the integers 1 to N.
//
{
  int i;
  int imult;
  int j;
  int value;

  value = 1;

  for ( i = 2; i <= n; i++ )
  {
    imult = i;
    for ( j = 1; j < i; j++ )
    {
      if ( ( imult % ( i - j ) ) == 0 )
      {
        imult = imult / ( i - j );
      }
    }
    value = value * imult;
  }

  return value;
}
//****************************************************************************80

void lmat_print ( int m, int n, bool a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    LMAT_PRINT prints an LMAT.
//
//  Discussion:
//
//    An LMAT is an array of L values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2011
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
//    Input, bool A[M*N], the matrix.
//
//    Input, string TITLE, a title.
//
{
  lmat_print_some ( m, n, a, 0, 0, m - 1, n - 1, title );

  return;
}
//****************************************************************************80

void lmat_print_some ( int m, int n, bool a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    LMAT_PRINT_SOME prints some of an LMAT.
//
//  Discussion:
//
//    An LMAT is an array of L values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, bool A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int incx = 35;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  for ( j2lo = i4_max ( jlo, 0 ); j2lo <= i4_min ( jhi, n - 1 ); j2lo = j2lo + incx )
  {
    j2hi = j2lo + incx - 1;
    j2hi = i4_min ( j2hi, n - 1 );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";

    if ( 100 <= j2hi )
    {
      cout << "      ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << " " << setw(1) << j / 100;
      }
      cout << "\n";
    }

    if ( 10 <= j2hi )
    {
      cout << "      ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << " " << setw(1) << ( j / 10 ) % 10;
      }
      cout << "\n";
    }

    cout << "  Col ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << " " << setw(1) << j % 10;
    }
    cout << "\n";

    cout << "  Row\n";
    cout << "\n";

    i2lo = i4_max ( ilo, 0 );
    i2hi = i4_min ( ihi, m - 1 );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(5) << i << ":";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << " " << setw(1) << a[i+j*m];
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void lmat_transpose_print ( int m, int n, bool a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    LMAT_TRANSPOSE_PRINT prints an LMAT, transposed.
//
//  Discussion:
//
//    An LMAT is an array of L values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, bool A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  lmat_transpose_print_some ( m, n, a, 0, 0, m - 1, n - 1, title );

  return;
}
//****************************************************************************80

void lmat_transpose_print_some ( int m, int n, bool a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    LMAT_TRANSPOSE_PRINT_SOME prints some of an LMAT, transposed.
//
//  Discussion:
//
//    An LMAT is an array of L values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, bool A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int incx = 35;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  for ( j2lo = i4_max ( ilo, 0 ); j2lo <= i4_min ( ihi, m - 1 ); j2lo = j2lo + incx )
  {
    i2hi = i2lo + incx - 1;
    i2hi = i4_min ( i2hi, m - 1 );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";

    if ( 100 <= i2hi )
    {
      cout << "      ";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << " " << setw(1) << i / 100;
      }
      cout << "\n";
    }

    if ( 10 <= i2hi )
    {
      cout << "      ";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << " " << setw(1) << ( i / 10 ) % 10;
      }
      cout << "\n";
    }

    cout << "  Row ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << " " << setw(1) << i % 10;
    }
    cout << "\n";

    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 0 );
    j2hi = i4_min ( jhi, n - 1 );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << ":";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << " " << setw(1) << a[i+j*m];
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

int luhn_check ( int digit_num, int digit[] ) 

//****************************************************************************80
//
//  Purpose:
//
//    LUHN_CHECK computes the Luhn checksum for a string of digits.
//
//  Discussion:
//
//    To compute the Luhn checksum, begin at the end of the string, and double
//    every other digit.  If a doubled digit is greater than 9, subtract 9.
//    Then sum the digits to get CHECK_SUM.
//
//    If mod ( CHECK_SUM, 10 ) = 0 the digit sequence is accepted.
//
//    The Luhn check sum will detect any single digit error, as well as
//    most errors involving a transposition of adjacent digits; it cannot
//    detect the transposition of the adjacent digits 0 and 9, however.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIGIT_NUM, the number of digits.
//
//    Input, int DIGIT[DIGIT_NUM], the string of digits.
//    Normally, these are true digits, that is, values between 0 and 9.
//
//    Output, int LUHN_CHECK, the check sum, which
//    should be divisible by 10 if the digit sequence is correct.
//
{
  int check_sum;
  int *digit_copy;
  int i;

  digit_copy = new int[digit_num];

  for ( i = 0; i < digit_num; i++ )
  {
    digit_copy[i] = digit[i];
  }

  for ( i = digit_num - 2; 0 <= i; i = i - 2 )
  {
    digit_copy[i] = 2 * digit_copy[i];
  }

  for ( i = 0; i < digit_num; i++ )
  {
    if ( 9 < digit_copy[i] )
    {
      digit_copy[i] = digit_copy[i] - 9;
    }
  }

  check_sum = i4vec_sum ( digit_num, digit_copy );

  delete [] digit_copy;

  return check_sum;
}
//****************************************************************************80

void lvec_print ( int n, bool a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    LVEC_PRINT prints a logical vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, bool A[N], the vector to be printed.
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
         << ": " << setw(1) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

int pause_input ( )

//****************************************************************************80
//
//  Purpose:
//
//    PAUSE_INPUT waits until an input character is entered.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
// 
//    13 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int PAUSE_INPUT, the character read from STDIN.
//
{
  int value;

  cout << "Press RETURN to continue.\n";

  value = getc ( stdin );

  return value;
}
//****************************************************************************80

void pause_seconds ( int seconds )

//****************************************************************************80
//
//  Purpose:
//
//    PAUSE_SECONDS waits a specified number of seconds.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
// 
//    13 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SECONDS, the number of seconds to pause.
//
{
  time_t t1;
  time_t t2;

  t1 = time ( NULL );
  t2 = t1;
  while ( t2 - t1 < seconds )
  {
    t2 = time ( NULL );
  }

  return;
}
//****************************************************************************80

bool perm_check ( int n, int p[], int base )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from BASE to
//    to BASE+N-1 occurs among the N entries of the permutation.
//
//    Set the input quantity BASE to 0, if P is a 0-based permutation,
//    or to 1 if P is a 1-based permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Input, int BASE, the index base.
//
//    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
//
{
  bool found;
  int i;
  int seek;

  for ( seek = base; seek < base + n; seek++ )
  {
    found = false;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = true;
        break;
      }
    }

    if ( !found )
    {
      return false;
    }

  }

  return true;
}
//****************************************************************************80

void perm_cycle ( int n, int p[], int *isgn, int *ncycle, int iopt )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CYCLE analyzes a permutation.
//
//  Discussion:
//
//    The routine will count cycles, find the sign of a permutation,
//    and tag a permutation.
//
//  Example:
//
//    Input:
//
//      N = 9
//      IOPT = 1
//      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
//
//    Output:
//
//      NCYCLE = 3
//      ISGN = +1
//      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
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
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N].  On input, P describes a
//    permutation, in the sense that entry I is to be moved to P[I].
//    If IOPT = 0, then P will not be changed by this routine.
//    If IOPT = 1, then on output, P will be "tagged".  That is,
//    one element of every cycle in P will be negated.  In this way,
//    a user can traverse a cycle by starting at any entry I1 of P
//    which is negative, moving to I2 = ABS(P[I1]), then to
//    P[I2], and so on, until returning to I1.
//
//    Output, int *ISGN, the "sign" of the permutation, which is
//    +1 if the permutation is even, -1 if odd.  Every permutation
//    may be produced by a certain number of pairwise switches.
//    If the number of switches is even, the permutation itself is
//    called even.
//
//    Output, int *NCYCLE, the number of cycles in the permutation.
//
//    Input, int IOPT, requests tagging.
//    0, the permutation will not be tagged.
//    1, the permutation will be tagged.
//
{
  int base = 1;
  int i;
  int i1;
  int i2;
  int is;

  if ( !perm_check ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "PERM_CYCLE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }

  is = 1;
  *ncycle = n;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p[i-1];

    while ( i < i1 )
    {
      *ncycle = *ncycle - 1;
      i2 = p[i1-1];
      p[i1-1] = -i2;
      i1 = i2;
    }

    if ( iopt != 0 )
    {
      is = - i4_sign ( p[i-1] );
    }
    p[i-1] = abs ( p[i-1] ) * i4_sign ( is );
  }

  *isgn = 1 - 2 * ( ( n - *ncycle ) % 2 );

  return;
}
//****************************************************************************80

int *perm_free ( int npart, int ipart[], int nfree )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_FREE reports the number of unused items in a partial permutation.
//
//  Discussion:
//
//    It is assumed that the N objects being permuted are the integers
//    from 1 to N, and that IPART contains a "partial" permutation, that
//    is, the NPART entries of IPART represent the beginning of a
//    permutation of all N items.
//
//    The routine returns in IFREE the items that have not been used yet.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NPART, the number of entries in IPART.  NPART may be 0.
//
//    Input, int IPART(NPART), the partial permutation, which should
//    contain, at most once, some of the integers between 1 and
//    NPART+NFREE.
//
//    Input, int NFREE, the number of integers that have not been
//    used in IPART.  This is simply N - NPART.  NFREE may be zero.
//
//    Output, int PERM_FREE[NFREE], the integers between 1 and NPART+NFREE
//    that were not used in IPART.
//
{
  int i;
  int *ifree;
  int j;
  int k;
  int match;
  int n;

  n = npart + nfree;

  if ( npart < 0 )
  {
    cerr << "\n";
    cerr << "PERM_FREE - Fatal error!\n";
    cerr << "  NPART < 0.\n";
    cerr << "  NPART = " << npart << "\n";
    exit ( 1 );
  }
  else if ( npart == 0 )
  {
    ifree = i4vec_indicator_new ( n );
    return ifree;
  }
  else if ( nfree < 0 )
  {
    cerr << "\n";
    cerr << "PERM_FREE - Fatal error!\n";
    cerr << "  NFREE < 0.\n";
    cerr << "  NFREE =  << nfree << \n";
    exit ( 1 );
  }
  else if ( nfree == 0 )
  {
    return NULL;
  }
  else
  {
    ifree = new int[nfree];

    k = 0;

    for ( i = 1; i <= n; i++ )
    {
      match = 0;

      for ( j = 1; j <= npart; j++ )
      {
        if ( ipart[j-1] == i )
        {
          match = j;
          break;
        }
      }

      if ( match == 0 )
      {
        k = k + 1;
        if ( nfree < k )
        {
          cerr << "\n";
          cerr << "PERM_FREE - Fatal error!\n";
          cerr << "  The partial permutation is illegal.\n";
          cerr << "  It should contain, at most once, some of\n";
          cerr << "  the integers between 1 and N = " << n << "\n";
          cerr << "  The number of integers that have not\n";
          cerr << "  been used is at least K = " << k << "\n";
          cerr << "  This should be exactly NFREE = " << nfree << "\n";
          i4vec_print ( npart, ipart, "  The partial permutation:" );
          exit ( 1 );
        }
        ifree[k-1] = i;
      }
    }
  }

  return ifree;
}
//****************************************************************************80

void perm_inverse ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE inverts a permutation "in place".
//
//  Discussion:
//
//    This algorithm assumes that the entries in the permutation vector are
//    strictly positive.  In particular, the value 0 must not occur.
//
//    When necessary, this function shifts the data temporarily so that
//    this requirement is satisfied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On output, P describes the inverse permutation
//
{
  int base;
  int i;
  int i0;
  int i1;
  int i2;
  int is;
  int p_min;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE - Fatal error!\n";
    cerr << "  Input value of N = " << n << "\n";
    exit ( 1 );
  }
//
//  Find the least value, and shift data so it begins at 1.
//
  p_min = i4vec_min ( n, p );
  base = 1;

  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - p_min + base;
  }
//
//  Now we can safely check the permutation.
//
  if ( !perm_check ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  Now we can invert the permutation.
//
  is = 1;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p[i-1];

    while ( i < i1 )
    {
      i2 = p[i1-1];
      p[i1-1] = -i2;
      i1 = i2;
    }

    is = - i4_sign ( p[i-1] );
    p[i-1] = i4_sign ( is ) * abs ( p[i-1] );
  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = - p[i-1];

    if ( 0 <= i1 )
    {
      i0 = i;

      for ( ; ; )
      {
        i2 = p[i1-1];
        p[i1-1] = i0;

        if ( i2 < 0 )
        {
          break;
        }
        i0 = i1;
        i1 = i2;
      }
    }
  }
//
//  Now we can restore the permutation.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + p_min - base;
  }

  return;
}
//****************************************************************************80

void perm_print ( int n, int p[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_PRINT prints a permutation.
//
//  Example:
//
//    Input:
//
//      P = 7 2 4 1 5 3 6
//
//    Printed output:
//
//      "This is the permutation:"
//
//      1 2 3 4 5 6 7
//      7 2 4 1 5 3 6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects permuted.
//
//    Input, int P[N], the permutation, in standard index form.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int ihi;
  int ilo;
  int inc = 20;

  cout << "\n";
  cout << title << "\n";

  for ( ilo = 0; ilo < n; ilo = ilo + inc )
  {
    ihi = ilo + inc;
    if ( n < ihi ) 
    {
      ihi = n;
    }
    cout << "\n";
    for ( i = ilo; i < ihi; i++ )
    {
      cout << setw(4) << i+1;
    }
    cout << "\n";
    for ( i = ilo; i < ihi; i++ )
    {
      cout << setw(4) << p[i];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

int *perm_uniform_new ( int n, int base, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM_NEW selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 October 2008
//
//  Author:
//
//    John Burkardt
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
//    Input, int N, the number of objects to be permuted.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for 
//    a 1-based permutation.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int PERM_UNIFORM_NEW[N], a permutation of 
//    (BASE, BASE+1, ..., BASE+N-1).
//
{
  int i;
  int j;
  int k;
  int *p;

  p = new int[n];
 
  for ( i = 0; i < n; i++ )
  {
    p[i] = i + base;
  }

  for ( i = 0; i < n; i++ )
  {
    j = i4_uniform ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }
 
  return p;
}
//****************************************************************************80

double pounds_to_kilograms ( double lb )

//****************************************************************************80
//
//  Purpose:
//
//    POUNDS_TO_KILOGRAMS converts a measurement in pounds to kilograms.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double LB, the weight in pounds.
//
//    Output, double POUNDS_TO_KILOGRAMS, the corresponding weight in kilograms.
//
{
  double value;

  value = 0.4535924 * lb;

  return value;
}
//****************************************************************************80

int prime ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME returns any of the first PRIME_MAX prime numbers.
//
//  Discussion:
//
//    PRIME_MAX is 1600, and the largest prime stored is 13499.
//
//    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 95-98.
//
//  Parameters:
//
//    Input, int N, the index of the desired prime number.
//    In general, is should be true that 0 <= N <= PRIME_MAX.
//    N = -1 returns PRIME_MAX, the index of the largest prime available.
//    N = 0 is legal, returning PRIME = 1.
//
//    Output, int PRIME, the N-th prime.  If N is out of range, PRIME
//    is returned as -1.
//
{
# define PRIME_MAX 1600

  int npvec[PRIME_MAX] = {
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541,
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657,
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553,
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, 
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, 
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, 
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, 
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, 
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, 
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, 
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, 
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, 
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 };

  if ( n == -1 )
  {
    return PRIME_MAX;
  }
  else if ( n == 0 )
  {
    return 1;
  }
  else if ( n <= PRIME_MAX )
  {
    return npvec[n-1];
  }
  else
  {
    cerr << "\n";
    cerr << "PRIME - Fatal error!\n";
    cerr << "  Unexpected input value of n = " << n << "\n";
    exit ( 1 );
  }

  return 0;
# undef PRIME_MAX
}
//****************************************************************************80

int prime_ge ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_GE returns the smallest prime greater than or equal to N.
//
//  Example:
//
//    N     PRIME_GE
//
//    -10    2
//      1    2
//      2    2
//      3    3
//      4    5
//      5    5
//      6    7
//      7    7
//      8   11
//      9   11
//     10   11
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number to be bounded.
//
//    Output, int PRIME_GE, the smallest prime number that is greater
//    than or equal to N.  However, if N is larger than the
//    largest prime stored, then PRIME_GE is returned as -1.
//
{
  int i_hi;
  int i_lo;
  int i_mid;
  int p;
  int p_hi;
  int p_lo;
  int p_mid;

  if ( n <= 2 )
  {
    p = 2;
  }
  else
  {
    i_lo = 1;
    p_lo = prime(i_lo);
    i_hi = prime(-1);
    p_hi = prime(i_hi);

    if ( p_hi < n )
    {
      p = - p_hi;
    }
    else
    {
      for ( ; ; )
      {
        if ( i_lo + 1 == i_hi )
        {
          p = p_hi;
          break;
        }

        i_mid = ( i_lo + i_hi ) / 2;
        p_mid = prime(i_mid);

        if ( p_mid < n )
        {
          i_lo = i_mid;
          p_lo = p_mid;
        }
        else if ( n <= p_mid )
        {
          i_hi = i_mid;
          p_hi = p_mid;
        }
      }
    }
  }

  return p;
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
    value = - ( int ) ( - x + 0.5 );
  }
  else
  {
    value =   ( int ) (  x + 0.5 );
  }
  return value;
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

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

double r8_log_10 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LOG_10 returns the logarithm base 10 of the absolute value of an R8.
//
//  Discussion:
//
//    value = Log10 ( |X| )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose base 2 logarithm is desired.
//    X should not be 0.
//
//    Output, double R8_LOG_10, the logarithm base 10 of the absolute
//    value of X.  It should be true that |X| = 10**R_LOG_10.
//
{
  double value;

  if ( x == 0.0 )
  {
    value = - r8_huge ( );
  }
  else
  {
    value = log10 ( r8_abs ( x ) );
  }

  return value;
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

double r8_modp ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MODP returns the nonnegative remainder of R8 division.
//
//  Discussion:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360.0) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
//
//    If
//      REM = R8_MODP ( X, Y )
//      RMULT = ( X - REM ) / Y
//    then
//      X = Y * RMULT + REM
//    where REM is always nonnegative.
//
//  Example:
//
//        I         J     MOD  R8_MODP   R8_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
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
//    Input, double X, the number to be divided.
//
//    Input, double Y, the number that divides X.
//
//    Output, double R8_MODP, the nonnegative remainder when X is divided by Y.
//
{
  double value;

  if ( y == 0.0 )
  {
    cerr << "\n";
    cerr << "R8_MODP - Fatal error!\n";
    cerr << "  R8_MODP ( X, Y ) called with Y = " << y << "\n";
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( value < 0.0 )
  {
    value = value + r8_abs ( y );
  }

  return value;
}
//****************************************************************************80

double r8_uniform ( double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM returns a scaled pseudorandom R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double B, C, the minimum and maximum values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8_UNIFORM, the randomly chosen value.
//
{
  double t;

  t = r8_uniform_01 ( seed );

  t = ( 1.0 - t ) * b + t * c;

  return t;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
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
//    11 August 2004
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
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
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
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
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
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
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
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
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
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << ":";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

int r8poly_degree ( int na, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_DEGREE returns the degree of a polynomial.
//
//  Discussion:
//
//    The degree of a polynomial is the index of the highest power
//    of X with a nonzero coefficient.
//
//    The degree of a constant polynomial is 0.  The degree of the
//    zero polynomial is debatable, but this routine returns the
//    degree as 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NA, the dimension of A.
//
//    Input, double A[NA+1], the coefficients of the polynomials.
//
//    Output, int R8POLY_DEGREE, the degree of A.
//
{
  int degree;

  degree = na;

  while ( 0 < degree )
  {
    if ( a[degree] != 0.0 )
    {
      return degree;
    }
    degree = degree - 1;
  }

  return degree;
}
//****************************************************************************80

void r8poly_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_PRINT prints out a polynomial.
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
//    Input, int N, the dimension of A.
//
//    Input, double A[N+1], the polynomial coefficients.
//    A(0) is the constant term and
//    A(N) is the coefficient of X**N.
//
//    Input, string TITLE, a title.
//
{
  int i;
  double mag;
  int n2;
  char plus_minus;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  n2 = r8poly_degree ( n, a );

  if ( n2 <= 0 ) 
  {
    cout << "  p(x) = 0\n";
    return;
  }

  if ( a[n2] < 0.0 )
  {
    plus_minus = '-';
  }
  else
  {
    plus_minus = ' ';
  }

  mag = r8_abs ( a[n2] );

  if ( 2 <= n2 )
  {
    cout << "  p(x) = " << plus_minus 
         << setw(14) << mag << " * x ^ " << n2 << "\n";
  }
  else if ( n2 == 1 )
  {
    cout << "  p(x) = " << plus_minus 
         << setw(14) << mag << " * x\n";
  }
  else if ( n2 == 0 )
  {
    cout << "  p(x) = " << plus_minus 
         << setw(14) << mag << "\n";
  }

  for ( i = n2-1; 0 <= i; i-- )
  {
    if ( a[i] < 0.0 )
    {
      plus_minus = '-';
    }
    else
    {
      plus_minus = '+';
    }

    mag = r8_abs ( a[i] );

    if ( mag != 0.0 )
    {
      if ( 2 <= i )
      {
        cout << "         " << plus_minus 
             << setw(14) << mag << " * x ^ " << i << "\n";
      }
      else if ( i == 1 )
      {
        cout << "         " << plus_minus 
             << setw(14) << mag << " * x\n";
      }
      else if ( i == 0 )
      {
        cout << "         " << plus_minus 
             << setw(14) << mag << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

double *r8vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDICATOR_NEW sets an R8VEC to the indicator vector {1,2,3...}.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, double R8VEC_INDICATOR_NEW[N], the indicator array.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i <= n-1; i++ ) 
  {
    a[i] = ( double ) ( i + 1 );
  }

  return a;
}
//****************************************************************************80

double r8vec_max ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double value;

  value = - r8_huge ( );

  if ( n <= 0 ) 
  {
    return value;
  }

  for ( i = 0; i < n; i++ ) 
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

double r8vec_mean ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MEAN returns the mean of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[N], the vector whose mean is desired.
//
//    Output, double R8VEC_MEAN, the mean, or average, of the vector entries.
//
{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ ) 
  {
    mean = mean + x[i];
  }

  mean = mean / ( double ) n;

  return mean;
}
///****************************************************************************80

double r8vec_min ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], the array to be checked.
//
//    Output, double R8VEC_MIN, the value of the minimum element.
//
{
  int i;
  double value;

  value = r8_huge ( );

  if ( n <= 0 ) 
  {
    return value;
  }

  for ( i = 0; i < n; i++ ) 
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
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
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

double r8vec_variance ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_VARIANCE returns the variance of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[N], the vector whose variance is desired.
//
//    Output, double R8VEC_VARIANCE, the variance of the vector entries.
//
{
  int i;
  double mean;
  double variance;

  mean = r8vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ ) 
  {
    variance = variance + ( x[i] - mean ) * ( x[i] - mean );
  }

  if ( 1 < n ) 
  {
    variance = variance / ( double ) ( n - 1 );
  }
  else
  {
    variance = 0.0;
  }

  return variance;
}
//****************************************************************************80

double *r8vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO_NEW creates and zeroes an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return a;
}
//****************************************************************************80

double radians_to_degrees ( double angle )

//****************************************************************************80
//
//  Purpose:
//
//    RADIANS_TO_DEGREES converts an angle from radians to degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, an angle in radians.
//
//    Output, double RADIANS_TO_DEGREES, the equivalent angle in degrees.
//
{
  double pi = 3.141592653589793;
  double value;

  value = ( angle / pi ) * 180.0;

  return value;
}
//****************************************************************************80

unsigned long rand_initialize ( unsigned long seed )

//****************************************************************************80
//
//  Purpose:
//
//    RAND_INITIALIZE initializes the random number generator.
//
//  Discussion:
//
//    If you don't initialize RAND, the random number generator, 
//    it will behave as though it were seeded with value 1.  
//    This routine will either take a user-specified seed, or
//    (if the user passes a 0) make up a "random" one.  In either
//    case, the seed is passed to SRAND (the appropriate routine 
//    to call when setting the seed for RAND).  The seed is also
//    returned to the user as the value of the function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned long SEED, is either 0, which means that the user
//    wants this routine to come up with a seed, or nonzero, in which
//    case the user has supplied the seed.
//
//    Output, unsigned long RAND_INITIALIZE, is the value of the seed
//    passed to SRAND, which is either the user's input value, or if
//    that was zero, the value selected by this routine.
//
{
  if ( seed != 0 )
  {
    cout << "\n";
    cout << "RAND_INITIALIZE:\n";
    cout << "  Initialize RAND with user SEED = " << seed << "\n";
  }
  else
  {
    seed = get_seed ( );

    cout << "\n";
    cout << "RAND_INITIALIZE:\n";
    cout << "  Initialize RAND with arbitrary SEED = " << seed << "\n";
  }
//
//  Now set the seed.
//
  srand ( seed );

  return seed;
}
//****************************************************************************80

unsigned long random_initialize ( unsigned long seed )

//****************************************************************************80
//
//  Purpose:
//
//    RANDOM_INITIALIZE initializes the RANDOM random number generator.
//
//  Discussion:
//
//    If you don't initialize RANDOM, the random number generator, 
//    it will behave as though it were seeded with value 1.  
//    This routine will either take a user-specified seed, or
//    (if the user passes a 0) make up a "random" one.  In either
//    case, the seed is passed to SRANDOM (the appropriate routine 
//    to call when setting the seed for RANDOM).  The seed is also
//    returned to the user as the value of the function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned long SEED, is either 0, which means that the user
//    wants this routine to come up with a seed, or nonzero, in which
//    case the user has supplied the seed.
//
//    Output, unsigned long RANDOM_INITIALIZE, is the value of the seed
//    passed to SRANDOM, which is either the user's input value, or if
//    that was zero, the value selected by this routine.
//
{
# define DEBUG 0

  if ( seed != 0 )
  {
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE:\n";
      cout << "  Initialize RANDOM with user SEED = " << seed << "\n";
    }
  }
  else
  {
    seed = get_seed ( );
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE:\n";
      cout << "  Initialize RANDOM with arbitrary SEED = " << seed << "\n";
    }
  }
//
//  Now set the seed.
//
  srandom ( seed );

  return seed;
# undef DEBUG
}
//****************************************************************************80

void rat_factor ( int m, int n, int factor_max, int *factor_num, int factor[], 
  int power[], int *mleft, int *nleft )

//****************************************************************************80
//
//  Purpose:
//
//    RAT_FACTOR factors a rational value into a product of prime factors.
//
//  Discussion:
//
//    ( M / N ) = ( MLEFT / NLEFT ) * Product ( I = 1 to FACTOR_NUM )
//      FACTOR(I)**POWER(I).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the top and bottom of a rational value.
//    The ratio of M and N must be positive.
//
//    Input, int FACTOR_MAX, the maximum number of factors for
//    which storage has been allocated.
//
//    Output, int *FACTOR_NUM, the number of prime factors of M/N.
//
//    Output, int FACTOR[FACTOR_MAX], the prime factors of M/N.
//
//    Output, int POWER[FACTOR_MAX].  POWER(I) is the power of
//    the FACTOR(I) in the representation of M/N.
//
//    Output, int *MLEFT, *NLEFT, the top and bottom of the factor of
//    M / N that remains.  If ABS ( MLEFT / NLEFT ) is not 1, then
//    the rational value was not completely factored.
//
{
  int i;
  int p;
  int prime_max;

  *factor_num = 0;

  *mleft = m;
  *nleft = n;
//
//  NLEFT should be nonnegative.
//
  if ( *nleft < 0 )
  {
    *mleft = -(*mleft);
    *nleft = -(*nleft);
  }

  if ( m == 0 || n == 0 )
  {
    return;
  }

  if ( m == n )
  {
    *factor_num = 1;
    factor[0] = 1;
    power[0] = 1;
    return;
  }
//
//  Find out how many primes we stored.
//
  prime_max = prime ( -1 );

  for ( i = 1; i <= prime_max; i++ )
  {
    p = prime ( i );

    if ( ( *nleft % p ) == 0 || ( abs ( *mleft ) % p ) == 0 )
    {
      if ( *factor_num < factor_max )
      {
        *factor_num = *factor_num + 1;
        factor[*factor_num-1] = p;
        power[*factor_num-1] = 0;
//
//  Divide MLEFT by PRIME(I) as often as you can.
//
        if ( ( abs ( *mleft ) % p ) == 0  )
        {
          for ( ; ; )
          {
            power[*factor_num-1] = power[*factor_num-1] + 1;
            *mleft = *mleft / p;

            if ( ( abs ( *mleft ) % p ) != 0 )
            {
              break;
            }
          }
        }
//
//  Divide NLEFT by PRIME(I) as often as you can.
//
        if ( ( *nleft % p ) == 0  )
        {
          for ( ; ; )
          {
            power[*factor_num-1] = power[*factor_num-1] - 1;
            *nleft = *nleft / p;

            if ( ( *nleft % p ) != 0 )
            {
              break;
            }
          }
        }

        if ( power[*factor_num-1] == 0 )
        {
          *factor_num = *factor_num - 1;
        }
      }
    }
  }

  return;
}
//****************************************************************************80

double rickey ( int ab, int bb, int er, double f, int h, int hb, int hp, 
  int r, int so, int tb )

//****************************************************************************80
//
//  Purpose:
//
//    RICKEY evaluates Branch Rickey's baseball index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Schwarz,
//    Looking Beyond the Batting Average,
//    The New York Times, 
//    Sunday, 1 August 2004.
//   
//    Branch Rickey,
//    Goodby to Some Old Baseball Ideas,
//    Life Magazine,
//    2 August 1954.
//
//  Parameters:
//
//    Input, int AB, number of at-bats.
//
//    Input, int BB, base on balls.
//
//    Input, int ER, earned runs.
//
//    Input, double F, the fielding rating.
//
//    Input, int H, number of hits.
//
//    Input, int HB, hit batsmen.
//
//    Input, int HP, hit by pitcher.
//
//    Input, int R, runs.
//
//    Input, int SO, strike outs.
//
//    Input, int TB, total bases.
//
//    Output, double RICKEY, the Branch Rickey index, an estimate for the
//    expected winning percentage of a team with the given statistics.
//    (0.5 has already been subtracted from this value.)
//
{
  double g;;
  double hitting;
  double pitching;

  hitting = 
      ( double ) ( h + bb + hp ) / ( double ) ( ab + bb + hp ) 
    + ( double ) ( 3 * ( tb - h ) ) / ( double ) ( 4 * ab ) 
    + ( double ) ( r ) / ( double ) ( h + bb + hp );

  pitching = 
      ( double ) ( h ) / ( double ) ( ab ) 
    + ( double ) ( bb + hb ) / ( double ) ( ab + bb + hb ) 
    + ( double ) ( er ) / ( double ) ( h + bb + hb ) 
    - ( double ) ( so ) / ( double ) ( 8 * ( ab + bb + hb ) );

  g = hitting - pitching - f;

  return g;
}
//****************************************************************************80

int *roots_to_i4poly ( int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    ROOTS_TO_I4POLY converts polynomial roots to polynomial coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of roots specified.
//
//    Input, int X[N], the roots.
//
//    Output, int ROOTS_TO_I4POLY[N+1], the coefficients of the polynomial.
//
{
  int *c;
  int i;
  int j;

  c = i4vec_zero_new ( n + 1 );
//
//  Initialize C to (0, 0, ..., 0, 1).
//  Essentially, we are setting up a divided difference table.
//
  c[n] = 1;
//
//  Convert to standard polynomial form by shifting the abscissas
//  of the divided difference table to 0.
//
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= n+1-j; i++ )
    {
      c[n-i] = c[n-i] - x[n+1-i-j] * c[n-i+1];
    }
  }
  return c;
}
//****************************************************************************80

double *roots_to_r8poly ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of roots specified.
//
//    Input, double X[N], the roots.
//
//    Output, double ROOTS_TO_R8POLY[N+1], the coefficients of the polynomial.
//
{
  double *c;
  int i;
  int j;

  c = r8vec_zero_new ( n + 1 );
//
//  Initialize C to (0, 0, ..., 0, 1).
//  Essentially, we are setting up a divided difference table.
//
  c[n] = 1.0;
//
//  Convert to standard polynomial form by shifting the abscissas
//  of the divided difference table to 0.
//
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= n+1-j; i++ )
    {
      c[n-i] = c[n-i] - x[n+1-i-j] * c[n-i+1];
    }
  }
  return c;
}
//****************************************************************************80

bool s_eqi ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, S2, two strings.
//
//    Output, bool S_EQI, is true if the strings are equal. 
//
{
  int i;
  int nchar;
  int s1_length;
  int s2_length;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  if ( s1_length < s2_length )
  {
    nchar = s1_length;
  }
  else
  {
    nchar = s2_length;
  }
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ ) 
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) ) 
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < s1_length ) 
  {
    for ( i = nchar; i < s1_length; i++ ) 
    {
      if ( s1[i] != ' ' ) 
      {
        return false;
      }
    } 
  }
  else if ( nchar < s2_length ) 
  {
    for ( i = nchar; i < s2_length; i++ )
    {
      if ( s2[i] != ' ' ) 
      {
        return false;
      }
    } 
  }

  return true;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n ) 
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

int s_to_i4 ( char *s, int *last, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4 reads an integer value from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a string to be examined.
//
//    Output, int *LAST, the last character of S used to make IVAL.
//
//    Output, bool *ERROR is TRUE if an error occurred.
//
//    Output, int *S_TO_I4, the integer value read from the string.
//    If the string is blank, then IVAL will be returned 0.
//
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = false;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  while ( *s ) 
  {
    c = s[i];
    i = i + 1;
//
//  Haven't read anything.
//
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read the sign, expecting digits.
//
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read at least one digit, expecting more.
//
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
//
//  If we read all the characters in the string, see if we're OK.
//
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = charstar_len_trim ( s );
  }
  else
  {
    *error = true;
    *last = 0;
  }

  return ival;
}
//****************************************************************************80

bool s_to_i4vec ( char *s, int n, int ivec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4VEC reads an I4VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, int IVEC[N], the values read from the string.
//
//    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
//
{
  bool error;
  int i;
  int lchar;

  error = false;

  for ( i = 0; i < n; i++ )
  {
    ivec[i] = s_to_i4 ( s, &lchar, &error );

    if ( error )
    {
      cerr << "\n";
      cerr << "S_TO_I4VEC - Fatal error!\n";
      cerr << "  S_TO_I4 returned error while reading item " << i << "\n";
      exit ( 1 );
    }

    s = s + lchar;

  }

  return error;
}
//****************************************************************************80

double s_to_r8 ( char *s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = charstar_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7//
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( char *s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  bool error;
  int i;
  int lchar;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;

  }

  return error;
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
//****************************************************************************80

void tuple_next2 ( int n, int xmin[], int xmax[], int x[], int *rank )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT2 computes the next element of an integer tuple space.
//
//  Discussion:
//
//    The elements X are N vectors.
//
//    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
//
//    The elements are produced one at a time.
//
//    The first element is
//      (XMIN(1), XMIN(2), ..., XMIN(N)),
//    the second is (probably)
//      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
//    and the last element is
//      (XMAX(1), XMAX(2), ..., XMAX(N))
//
//    Intermediate elements are produced in a lexicographic order, with
//    the first index more important than the last, and the ordering of
//    values at a fixed index implicitly defined by the sign of
//    XMAX(I) - XMIN(I).
//
//  Example:
//
//    N = 2,
//    XMIN = (/ 1, 10 /)
//    XMAX = (/ 3,  8 /)
//
//    RANK    X
//    ----  -----
//      1   1 10
//      2   1  9
//      3   1  8
//      4   2 10
//      5   2  9
//      6   2  8
//      7   3 10
//      8   3  9
//      9   3  8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, int XMIN[N], XMAX[N], the "minimum" and "maximum" entry values.
//    These values are minimum and maximum only in the sense of the 
//    lexicographic ordering.  In fact, XMIN(I) may be less than, equal to, 
//    or greater than XMAX(I).
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
//    Input/output, int *RANK, the rank of the item.  On first call,
//    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
//    there are no more items in the sequence.
//
{
  int i;
  int test;

  if ( *rank < 0 )
  {
    cerr << "\n";
    cerr << "TUPLE_NEXT2 - Fatal error!\n";
    cerr << "  Illegal value of RANK = " << *rank << "\n";
    exit ( 1 );
  }

  test = 1;
  for ( i = 0; i < n; i++ )
  {
    test = test * ( 1 + abs ( xmax[i] - xmin[i] ) );
  }

  if ( test < *rank )
  {
    cerr << "\n";
    cerr << "TUPLE_NEXT2 - Fatal error!\n";
    cerr << "  Illegal value of RANK = " << *rank << "\n";
    exit ( 1 );
  }

  if ( *rank == 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = xmin[i];
    }
    *rank = 1;
    return;
  }

  *rank = *rank + 1;
  i = n - 1;

  for ( ; ; )
  {
    if ( x[i] != xmax[i] )
    {
      if ( xmin[i] < xmax[i] )
      {
        x[i] = x[i] + 1;
      }
      else
      {
        x[i] = x[i] - 1;
      }
      break;
    }

    x[i] = xmin[i];

    if ( i == 0 )
    {
      *rank = 0;
      break;
    }

    i = i - 1;

  }

  return;
}
//****************************************************************************80

double *tvec_even ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN computes an evenly spaced set of angles between 0 and 2*PI.
//
//  Discussion:
//
//    The computation realizes that 0 = 2 * PI.
//
//  Example:
//
//    NT = 4
//
//    T = ( 0, PI/2, PI, 3*PI/2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Output, double TVEC[NT], the evenly spaced angles, in radians.
//
{
  int i;
  double pi = 3.141592653589793;
  double *t;

  if ( nt < 1 )
  {
    return NULL;
  }

  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = ( double ) ( 2 * ( i - 1 ) ) * pi / ( double ) ( nt );
  }

  return t;
}
//****************************************************************************80

double *tvec_even2 ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN2 computes evenly spaced angles between 0 and 2*PI.
//
//  Discussion:
//
//    The computation realizes that 0 = 2 * PI.  The values are equally
//    spaced in the circle, do not include 0, and are symmetric about 0.
//
//  Example:
//
//    NT = 4
//
//    T = ( PI/4, 3*PI/4, 5*PI/4, 7*PI/4 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Output, double TVEC[NT], the evenly spaced angles, in radians.
//
{
  int i;
  double pi = 3.141592653589793;
  double *t;

  if ( nt < 1 )
  {
    return NULL;
  }

  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = ( double ) ( 2 * i - 1 ) * pi / ( double ) ( nt );
  }

  return t;
}
//****************************************************************************80

double *tvec_even3 ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN3 computes an evenly spaced set of angles between 0 and 2*PI.
//
//  Discussion:
//
//    The angles begin with 0 and end with 2*PI.
//
//  Example:
//
//    NT = 4
//
//    T = ( 0, 2*PI/3, 4*PI/3 2*PI )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Output, double TVEC[NT], the evenly spaced angles, in radians.
//
{
  int i;
  double pi = 3.141592653589793;
  double *t;

  if ( nt < 1 )
  {
    return NULL;
  }

  t = new double[nt];

  if ( nt == 1 )
  {
    t[0] = pi;
  }
  else
  {
    for ( i = 1; i <= nt; i++ )
    {
      t[i-1] = ( double ) ( 2 * ( i - 1 ) ) * pi / ( double ) ( nt - 1 );
    }
  }

  return t;
}
//****************************************************************************80

double *tvec_even_bracket ( int nt, double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN_BRACKET computes evenly spaced angles between THETA1 and THETA2.
//
//  Example:
//
//    NT = 4
//    THETA1 = 30
//    THETA2 = 90
//
//    T = ( 30, 50, 70, 90 )
//
//  Discussion:
//
//    The interval between THETA1 and THETA2 is divided into NT-1 subintervals.
//
//    The angles returned are the breakpoints of these subintervals,
//    including THETA1 and THETA2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Input, double THETA1, THETA2, the limiting angles.
//
//    Output, double TVEC_EVEN_BRACKET[NT], the evenly spaced angles.
//
{
  int i;
  double *t;

  if ( nt < 1 ) 
  {
    return NULL;
  }

  t = new double[nt];

  if ( nt == 1 )
  {
    t[0] = ( theta1 + theta2 ) / 2.0;
  }
  else
  {
    for ( i = 1; i <= nt; i++ )
    {
      t[i-1] = ( ( double ) ( nt - i     ) * theta1   
               + ( double ) (      i - 1 ) * theta2 ) 
               / ( double ) ( nt     - 1 );
    }
  }

  return t;
}
//****************************************************************************80

double *tvec_even_bracket2 ( int nt, double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN_BRACKET2 computes evenly spaced angles from THETA1 to THETA2.
//
//  Discussion:
//
//    The interval between THETA1 and THETA2 is divided into NT+1 subintervals.
//
//    The angles returned are the internal NT breakpoints of the subintervals.
//
//  Example:
//
//    NT = 5
//    THETA1 = 30
//    THETA2 = 90
//
//    T = ( 40, 50, 60, 70, 80 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Input, double THETA1, THETA2, the limiting angles.
//
//    Output, double TVEC_EVEN_BRACKET2[NT], the evenly spaced angles.
//
{
  int i;
  double *t;

  if ( nt < 1 ) 
  {
    return NULL;
  }

  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = ( ( double ) ( nt + 1 - i ) * theta1   
             + ( double ) (          i ) * theta2 ) 
             / ( double ) ( nt + 1     );
  }

  return t;
}
//****************************************************************************80

double *tvec_even_bracket3 ( int nt, double theta1, double theta2 )

//****************************************************************************80
//
//  Purpose:
//
//    TVEC_EVEN_BRACKET3 computes evenly spaced angles from THETA1 to THETA2.
//
//  Discussion:
//
//    The interval between THETA1 and THETA2 is divided into NT subintervals.
//
//    The angles returned are the midpoints of each subinterval.
//
//  Example:
//
//    NT = 3
//    THETA1 = 30
//    THETA2 = 90
//
//    T = ( 40, 60, 80 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the number of values to compute.
//
//    Input, double THETA1, THETA2, the limiting angles.
//
//    Output, double TVEC_EVEN_BRACKET3[NT], the evenly spaced angles.
//
{
  int i;
  double *t;

  t = new double[nt];

  for ( i = 1; i <= nt; i++ )
  {
    t[i-1] = ( ( double ) ( 2 * nt - 2 * i + 1 ) * theta1   
             + ( double ) (          2 * i - 1 ) * theta2 ) 
             / ( double ) ( 2 * nt             );
  }

  return t;
}
//****************************************************************************80

int upc_check_digit ( int p, int l, int r )

//****************************************************************************80
//
//  Purpose:
//
//    UPC_CHECK_DIGIT returns the check digit of a UPC.
//
//  Discussion:
//
//    UPC stands for Universal Price Code.
//
//    A full UPC is a string of 12 digits, in groups of size 1, 5, 5, and 1,
//    of the form P-LLLLL-RRRRR-C, where:
//
//      P is the one-digit product type code.
//      L is the five-digit manufacturer code.
//      R is the five_digit product code
//      C is the check digit.
//
//  Example:
//
//    0-72890-00011-8
//    0-12345-67890-5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the one-digit product type code.
//
//    Input, int L, the five-digit manufacturer code.
//
//    Input, int R, the five-digit product code.
//
//    Output, int UPC_CHECK_DIGIT, the check digit.
//
{
  int c;
  int *lc;
  int *rc;

  if ( p < 0 || 9 < p )
  {
    cerr << "\n";
    cerr << "UPC_CHECK_DIGIT - Fatal error!\n";
    cerr << "  P < 0 or 9 < P.\n";
    exit ( 1 );
  }

  if ( l < 0 || 99999 < l )
  {
    cerr << "\n";
    cerr << "UPC_CHECK_DIGIT - Fatal error!\n";
    cerr << "  L < 0 or 99999 < L.\n";
    exit ( 1 );
  }

  if ( r < 0 || 99999 < r )
  {
    cerr << "\n";
    cerr << "UPC_CHECK_DIGIT - Fatal error!\n";
    cerr << "  R < 0 or 99999 < R.\n";
    exit ( 1 );
  }

  lc = i4_to_digits_decimal ( l, 5 );
  rc = i4_to_digits_decimal ( r, 5 );

  c = ( p + lc[1] + lc[3] + rc[0] + rc[2] + rc[4] ) * 3 
          + lc[0] + lc[2] + lc[4] + rc[1] + rc[3];

  c = c % 10;

  c = ( 10 - c ) % 10;

  delete [] rc;
  delete [] lc;

  return c;
}
//****************************************************************************80

double versine_pulse ( double t, double ta, double tb, double v1, double amp )

//****************************************************************************80
//
//  Purpose:
//
//    VERSINE_PULSE adds a versine pulse to a constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the current time.
//
//    Input, double TA, the time at which the pulse begins.
//
//    Input, double TB, the time at which the pulse finishes.
//
//    Input, double V1, the constant value.
//
//    Input, double AMP, the amplitude of the pulse.
//
//    Output, double VERSINE_PULSE, the value of the signal at time T.
//
{
  double pi = 3.141592653589793;
  double v;

  v = v1;

  if ( ta <= t && t <= tb )
  {
    v = v + ( 0.5 * amp * ( 1.0 - cos ( 2.0 * pi * ( t - ta ) / ( tb - ta ) ) ) );
  }

  return v;
}

