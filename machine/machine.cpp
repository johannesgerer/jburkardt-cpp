# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "machine.hpp"

//****************************************************************************80

double d1mach ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    D1MACH returns double precision real machine constants.
//
//  Discussion:
//
//    Assuming that the internal representation of a double precision real
//    number is in base B, with T the number of base-B digits in the mantissa,
//    and EMIN the smallest possible exponent and EMAX the largest possible 
//    exponent, then
//
//      D1MACH(1) = B^(EMIN-1), the smallest positive magnitude.
//      D1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
//      D1MACH(3) = B^(-T), the smallest relative spacing.
//      D1MACH(4) = B^(1-T), the largest relative spacing.
//      D1MACH(5) = log10(B).
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
//    Output, double D1MACH, the value of the chosen parameter.
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
  else if ( 5 < i )
  {
    cerr << "\n";
    cerr << "D1MACH - Fatal error!\n";
    cerr << "  The input argument I is out of bounds.\n";
    cerr << "  Legal values satisfy 1 <= I <= 5.\n";
    cerr << "  I = " << i << "\n";
    value = 0.0;
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

int i1mach ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I1MACH returns integer machine constants.
//
//  Discussion:
//
//    Input/output unit numbers.
//
//      I1MACH(1) = the standard input unit.
//      I1MACH(2) = the standard output unit.
//      I1MACH(3) = the standard punch unit.
//      I1MACH(4) = the standard error message unit.
//
//    Words.
//
//      I1MACH(5) = the number of bits per integer storage unit.
//      I1MACH(6) = the number of characters per integer storage unit.
//
//    Integers.
//
//    Assume integers are represented in the S digit base A form:
//
//      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
//
//    where 0 <= X(1:S-1) < A.
//
//      I1MACH(7) = A, the base.
//      I1MACH(8) = S, the number of base A digits.
//      I1MACH(9) = A^S-1, the largest integer.
//
//    Floating point numbers
//
//    Assume floating point numbers are represented in the T digit 
//    base B form:
//
//      Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B^T) )
//
//    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
//
//      I1MACH(10) = B, the base.
//
//    Single precision
//
//      I1MACH(11) = T, the number of base B digits.
//      I1MACH(12) = EMIN, the smallest exponent E.
//      I1MACH(13) = EMAX, the largest exponent E.
//
//    Double precision
//
//      I1MACH(14) = T, the number of base B digits.
//      I1MACH(15) = EMIN, the smallest exponent E.
//      I1MACH(16) = EMAX, the largest exponent E.
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
//    Output, int I1MACH, the value of the chosen parameter.
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
    cerr << "I1MACH - Fatal error!\n";
    cerr << "  The input argument I is out of bounds.\n";
    cerr << "  Legal values satisfy 1 <= I <= 16.\n";
    cerr << "  I = " << i << "\n";
    value = 0;
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

float r1mach ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R1MACH returns single precision real machine constants.
//
//  Discussion:
//
//    Assume that single precision real numbers are stored with a mantissa 
//    of T digits in base B, with an exponent whose value must lie 
//    between EMIN and EMAX.  Then for values of I between 1 and 5, 
//    R1MACH will return the following values:
//
//      R1MACH(1) = B^(EMIN-1), the smallest positive magnitude.
//      R1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
//      R1MACH(3) = B^(-T), the smallest relative spacing.
//      R1MACH(4) = B^(1-T), the largest relative spacing.
//      R1MACH(5) = log10(B)
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
//    Output, float R1MACH, the value of the chosen parameter.
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
    cerr << "R1MACH - Fatal error!\n";
    cerr << "  The input argument I is out of bounds.\n";
    cerr << "  Legal values satisfy 1 <= I <= 5.\n";
    cerr << "  I = " << i << "\n";
    value = 0.0;
    exit ( 1 );
  }
 
  return value;
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
//    24 September 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

