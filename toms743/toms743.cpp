# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "toms743.hpp"

//****************************************************************************80

double bisect ( double xx, int nb, int &ner, int l )

//****************************************************************************80
//
//  Purpose:
//
//    BISECT approximates the W function using bisection.
//
//  Discussion:
//
//    The parameter TOL, which determines the accuracy of the bisection
//    method, is calculated using NBITS (assuming the final bit is lost
//    due to rounding error).
//
//    N0 is the maximum number of iterations used in the bisection
//    method.
//
//    For XX close to 0 for Wp, the exponential approximation is used.
//    The approximation is exact to O(XX^8) so, depending on the value
//    of NBITS, the range of application of this formula varies. Outside
//    this range, the usual bisection method is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
//    Patricia Culligan-Hensley.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
//    Algorithm 743: WAPR - A Fortran routine for calculating real 
//    values of the W-function,
//    ACM Transactions on Mathematical Software,
//    Volume 21, Number 2, June 1995, pages 172-181.
//
//  Parameters:
//
//    Input, double XX, the argument.
//
//    Input, int NB, indicates the branch of the W function.
//    0, the upper branch;
//    nonzero, the lower branch.
//
//    Output, int &NER, the error flag.
//    0, success;
//    1, the routine did not converge.  Perhaps reduce NBITS and try again.
//
//    Input, int L, the offset indicator.
//    1, XX represents the offset of the argument from -exp(-1).
//    not 1, XX is the actual argument.
//
//    Output, double BISECT, the value of W(X), as determined
//
{
  double d;
  double f;
  double fd;
  int i;
  const int n0 = 500;
  static int nbits = 0;
  double r;
  double test;
  double tol;
  double u;
  double value;
  double x;

  value = 0.0;
  ner = 0;

  if ( nbits == 0 )
  {
    nbits = nbits_compute ( );
  }

  if ( l == 1 )
  {
    x = xx - exp ( -1.0 );
  }
  else
  {
    x = xx;
  }

  if ( nb == 0 )
  {
    test = 1.0 / pow ( pow ( 2.0, nbits ), ( 1.0 / 7.0 ) );

    if ( fabs ( x ) < test )
    {
      value = x 
        * exp ( - x 
        * exp ( - x 
        * exp ( - x 
        * exp ( - x 
        * exp ( - x 
        * exp ( - x ))))));

      return value;
    }
    else
    {
      u = crude ( x, nb ) + 1.0E-03;
      tol = fabs ( u ) / pow ( 2.0, nbits );
      d = r8_max ( u - 2.0E-03, -1.0 );

      for ( i = 1; i <= n0; i++ )
      {
        r = 0.5 * ( u - d );
        value = d + r;
//
//  Find root using w*exp(w)-x to avoid ln(0) error.
//
        if ( x < exp ( 1.0 ) )
        {
          f = value * exp ( value ) - x;
          fd = d * exp ( d ) - x;
        }
//
//  Find root using ln(w/x)+w to avoid overflow error.
//
        else
        {
          f = log ( value / x ) + value;
          fd = log ( d / x ) + d;
        }

        if ( f == 0.0 )
        {
          return value;
        }

        if ( fabs ( r ) <= tol )
        {
          return value;
        }

        if ( 0.0 < fd * f )
        {
          d = value;
        }
        else
        {
          u = value;
        }
      }
    }
  }
  else
  {
    d = crude ( x, nb ) - 1.0E-03;
    u = r8_min ( d + 2.0E-03, -1.0 );
    tol = fabs ( u ) / pow ( 2.0, nbits );

    for ( i = 1; i <= n0; i++ )
    {
      r = 0.5 * ( u - d );
      value = d + r;
      f = value * exp ( value ) - x;

      if ( f == 0.0 )
      {
        return value;
      }

      if ( fabs ( r ) <= tol )
      {
        return value;
      }

      fd = d * exp ( d ) - x;

      if ( 0.0 < fd * f )
      {
        d = value;
      }
      else
      {
        u = value;
      }
    }
  }
//
//  The iteration did not converge.
//
  ner = 1;

  return value;
}
//****************************************************************************80

double crude ( double xx, int nb )

//****************************************************************************80
//
//  Purpose:
//
//    CRUDE returns a crude approximation for the W function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
//    Patricia Culligan-Hensley.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
//    Algorithm 743: WAPR - A Fortran routine for calculating real 
//    values of the W-function,
//    ACM Transactions on Mathematical Software,
//    Volume 21, Number 2, June 1995, pages 172-181.
//
//  Parameters:
//
//    Input, double XX, the argument.
//
//    Input, int NB, indicates the desired branch.
//    * 0, the upper branch;
//    * nonzero, the lower branch.
//
//    Output, double CRUDE, the crude approximation to W at XX.
//
{
  double an2;
  static double c13;
  double crude;
  static double em;
  static double em2;
  static double em9;
  double eta;
  static int init = 0;
  double reta;
  static double s2;
  static double s21;
  static double s22;
  static double s23;
  double t;
  double ts;
  double value;
  double zl;

  value = 0.0;
//
//  Various mathematical constants.
//
  if ( init == 0 )
  {
    init = 1;
    em = - exp ( -1.0 );
    em9 = - exp ( -9.0 );
    c13 = 1.0 / 3.0;
    em2 = 2.0 / em;
    s2 = sqrt ( 2.0 );
    s21 = 2.0 * s2 - 3.0;
    s22 = 4.0 - 3.0 * s2;
    s23 = s2 - 2.0;
  }
//
//  Crude Wp.
//
  if ( nb == 0 )
  {
    if ( xx <= 20.0 )
    {
      reta = s2 * sqrt ( 1.0 - xx / em );
      an2 = 4.612634277343749 * sqrt ( sqrt ( reta + 
        1.09556884765625 ) );
      value = reta / ( 1.0 + reta / ( 3.0 
        + ( s21 * an2 + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.0;
    }
    else
    {
      zl = log ( xx );
      value = log ( xx / log ( xx 
        / pow ( zl, exp ( -1.124491989777808 / 
        ( 0.4225028202459761 + zl )) ) ));
    }
  }
//
//  Crude Wm.
//
  else
  {
    if ( xx <= em9 )
    {
      zl = log ( -xx );
      t = -1.0 - zl;
      ts = sqrt ( t );
      value = zl - ( 2.0 * ts ) / ( s2 + ( c13 - t 
        / ( 270.0 + ts * 127.0471381349219 ) ) * ts );
    }
    else
    {
      zl = log ( -xx );
      eta = 2.0 - em2 * xx;
      value = log ( xx / log ( - xx / ( ( 1.0 
        - 0.5043921323068457 * ( zl + 1.0 ) ) 
        * ( sqrt ( eta ) + eta / 3.0 ) + 1.0 ) ) );
     }
  }

  return value;
}
//****************************************************************************80

int nbits_compute ( )

//****************************************************************************80
//
//  Purpose:
//
//    NBITS_COMPUTE computes the mantissa length minus one.
//
//  Discussion:
//
//    NBITS is the number of bits (less 1) in the mantissa of the
//    floating point number number representation of your machine.
//    It is used to determine the level of accuracy to which the W
//    function should be calculated.
//
//    Most machines use a 24-bit matissa for single precision and
//    53-56 bits for double. The IEEE standard is 53
//    bits. The Fujitsu VP2200 uses 56 bits. Long word length
//    machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
//    single precision.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
//    Patricia Culligan-Hensley.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
//    Algorithm 743: WAPR - A Fortran routine for calculating real 
//    values of the W-function,
//    ACM Transactions on Mathematical Software,
//    Volume 21, Number 2, June 1995, pages 172-181.
//
//  Parameters:
//
//    Output, int NBITS_COMPUTE, the mantissa length, in bits, 
//    minus one.
//
{
  double b;
  int i;
  int nbits;
  double v;

  nbits = 0;

  b = 1.0;

  for ( ; ; )
  {
    b = b / 2.0;
    v = b + 1.0;

    if ( v == 1.0 )
    {
      break;
    }
    nbits = nbits + 1;
  }

  return nbits;
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
//    07 May 2006
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
//    07 May 2006
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

double wapr ( double x, int nb, int &nerror, int l )

//****************************************************************************80
//
// WAPR approximates the W function.
//
//  Discussion:
//
//    The call will fail if the input value X is out of range.
//    The range requirement for the upper branch is:
//      -exp(-1) <= X.
//    The range requirement for the lower branch is:
//      -exp(-1) < X < 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
//    Patricia Culligan-Hensley.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
//    Algorithm 743: WAPR - A Fortran routine for calculating real 
//    values of the W-function,
//    ACM Transactions on Mathematical Software,
//    Volume 21, Number 2, June 1995, pages 172-181.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Input, int NB, indicates the desired branch.
//    * 0, the upper branch;
//    * nonzero, the lower branch.
//
//    Output, int &NERROR, the error flag.
//    * 0, successful call.
//    * 1, failure, the input X is out of range.
//
//    Input, int L, indicates the interpretation of X.
//    * 1, X is actually the offset from -(exp-1), so compute W(X-exp(-1)).
//    * not 1, X is the argument; compute W(X);
//
//    Output, double WAPR, the approximate value of W(X).
//
{
  double an2;
  static double an3;
  static double an4;
  static double an5;
  static double an6;
  static double c13;
  static double c23;
  static double d12;
  double delx;
  static double em;
  static double em2;
  static double em9;
  double eta;
  int i;
  static int init = 0;
  int m;
  static int nbits;
  static int niter = 1;
  double reta;
  static double s2;
  static double s21;
  static double s22;
  static double s23;
  double t;
  static double tb;
  static double tb2;
  double temp;
  double temp2;
  double ts;
  double value;
  static double x0;
  static double x1;
  double xx;
  double zl;
  double zn;

  value = 0.0;
  nerror = 0;

  if ( init == 0 )
  {
    init = 1;

    nbits = nbits_compute ( );

    if ( 56 <= nbits )
    {
      niter = 2;
    }
//
//  Various mathematical constants.
//
    em = -exp ( -1.0 );
    em9 = -exp ( -9.0 );
    c13 = 1.0 / 3.0;
    c23 = 2.0 * c13;
    em2 = 2.0 / em;
    d12 = -em2;
    tb = pow ( 0.5, nbits );
    tb2 = sqrt ( tb );
    x0 = pow ( tb, 1.0 / 6.0 ) * 0.5;
    x1 = ( 1.0 - 17.0 * pow ( tb, 2.0 / 7.0 ) ) * em;
    an3 = 8.0 / 3.0;
    an4 = 135.0 / 83.0;
    an5 = 166.0 / 39.0;
    an6 = 3167.0 / 3549.0;
    s2 = sqrt ( 2.0 );
    s21 = 2.0 * s2 - 3.0;
    s22 = 4.0 - 3.0 * s2;
    s23 = s2 - 2.0;
  }

  if ( l == 1 )
  {
    delx = x;

    if ( delx < 0.0 )
    {
      nerror = 1;
      cerr << "\n";
      cerr << "WAPR - Fatal error!\n";
      cerr << "  The offset X is negative.\n";
      cerr << "  It must be nonnegative.\n";
      exit ( 1 );
    }

    xx = x + em;
  }
  else
  {
    if ( x < em )
    {
      nerror = 1;
      return value;
    }
    else if ( x == em )
    {
      value = -1.0;
      return value;
    }
    xx = x;
    delx = xx - em;
  }
//
//  Calculations for Wp.
//
  if ( nb == 0 )
  {
    if ( fabs ( xx ) <= x0 )
    {
      value = xx / ( 1.0 + xx / ( 1.0 + xx 
        / ( 2.0 + xx / ( 0.6 + 0.34 * xx ))));
      return value;
    }
    else if ( xx <= x1 )
    {
      reta = sqrt ( d12 * delx );
      value = reta / ( 1.0 + reta / ( 3.0 + reta / ( reta 
        / ( an4 + reta / ( reta * an6 + an5 ) ) + an3 ) ) ) 
        - 1.0;
      return value;
    }
    else if ( xx <= 20.0 )
    {
      reta = s2 * sqrt ( 1.0 - xx / em );
      an2 = 4.612634277343749 * sqrt ( sqrt ( reta + 
        1.09556884765625 ));
      value = reta / ( 1.0 + reta / ( 3.0 + ( s21 * an2 
        + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.0;
    }
    else
    {
      zl = log ( xx );
      value = log ( xx / log ( xx 
        / pow ( zl, exp ( -1.124491989777808 / 
        ( 0.4225028202459761 + zl )) ) ));
    }
  }
//
//  Calculations for Wm.
//
  else
  {
    if ( 0.0 <= xx )
    {
      nerror = 1;
      return value;
    }
    else if ( xx <= x1 )
    {
      reta = sqrt ( d12 * delx );
      value = reta / ( reta / ( 3.0 + reta / ( reta / ( an4 
        + reta / ( reta * an6 - an5 ) ) - an3 ) ) - 1.0 ) - 1.0;
      return value;
    }
    else if ( xx <= em9 )
    {
      zl = log ( -xx );
      t = -1.0 - zl;
      ts = sqrt ( t );
      value = zl - ( 2.0 * ts ) / ( s2 + ( c13 - t 
        / ( 270.0 + ts * 127.0471381349219 )) * ts );
    }
    else
    {
      zl = log ( -xx );
      eta = 2.0 - em2 * xx;
      value = log ( xx / log ( -xx / ( ( 1.0 
        - 0.5043921323068457 * ( zl + 1.0 ) ) 
        * ( sqrt ( eta ) + eta / 3.0 ) + 1.0 )));
    }

  }

  for ( i = 1; i <= niter; i++ )
  {
    zn = log ( xx / value ) - value;
    temp = 1.0 + value;
    temp2 = temp + c23 * zn;
    temp2 = 2.0 * temp * temp2;
    value = value * ( 1.0 + ( zn / temp ) * ( temp2 - zn ) 
      / ( temp2 - 2.0 * zn ) );
  }

  return value;
}
