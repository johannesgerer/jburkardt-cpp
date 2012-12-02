# include <cstdlib>
# include <iostream>
# include <ctime>

# include "randlc.hpp"

using namespace std;

//****************************************************************************80

double randlc ( double *x )

//****************************************************************************80
//
//  Purpose:
//
//    RANDLC returns a uniform pseudorandom double precision number.
//
//  Discussion:
//
//    This function, when called repeatedly, computes a sequence of
//    pseudorandom double precision numbers in the range (0,1).
//
//    The calculation involves a series of integer seeds X(0), X(1),
//    ..., X(K+1), and a fixed multiplier A, so that
//
//      X(K+1) = A * X(K)  mod 2^46
//
//    and the pseudorandom value is then computed by
//
//      RANDLC = X(K+1) / 2^46.
//
//    This scheme generates 2^44 numbers before repeating.  
//
//    The multiplier A and the seed X must be odd double precision integers
//    in the range (1, 2^46). 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 March 2010
//
//  Author:
//
//    Original C version by David Bailey.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Bailey, Eric Barszcz, John Barton, D Browning, Robert Carter, 
//    Leonardo Dagum, Rod Fatoohi,
//    Samuel Fineberg, Paul Frederickson, Thomas Lasinski, Robert Schreiber, 
//    Horst Simon, V Venkatakrishnan, Sisira Weeratunga,
//    The NAS Parallel Benchmarks,
//    RNR Technical Report RNR-94-007,
//    March 1994.
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 2, Seminumerical Algorithms,
//    Third Edition,
//    Addison Wesley, 1997,
//    ISBN: 0201896842,
//    LC: QA76.6.K64.
//
//  Parameters:
//
//    Input/output, double *X, the seed value.  On input, X should be
//    an odd integer value in the range 1, 2^46.  On output, X has been
//    updated by the linear congruential function.  If the input value
//    is zero, it is replaced by 314159265.  If the input value is negative, 
//    it is replaced by its absolute value.
//
//    Output, double RANDLC, the pseudorandom value corresponding to
//    the input value of X. 
//
{
  static double a = 1220703125.00;
  static double a1;
  static double a2;
  int i;
  int j;
  static int ks = 0;
  static double r23;
  static double r46;
  double t1;
  double t2;
  static double t23;
  double t3;
  double t4;
  static double t46;
  double value;
  double x1;
  double x2;
  double z;
//
//  If this is the first call, compute 
//
//    R23 = 2 ^ -23, 
//    R46 = 2 ^ -46,
//    T23 = 2 ^ 23, 
//    T46 = 2 ^ 46.  
//
//  These are computed in loops, rather than by merely using the power operator, 
//  in order to insure that the results are exact on all systems.  
//
  if ( ks == 0 ) 
  {
    r23 = 1.0;
    r46 = 1.0;
    t23 = 1.0;
    t46 = 1.0;
    
    for ( i = 1; i <= 23; i++ )
    {
      r23 = 0.50 * r23;
      t23 = 2.0 * t23;
    }
    for ( i = 1; i <= 46; i++ )
    {
      r46 = 0.50 * r46;
      t46 = 2.0 * t46;
    }
//
//  Break A into two parts such that A = 2^23 * A1 + A2.
//
    t1 = r23 * a;
    a1 = ( double ) ( int ) t1;
    a2 = a - t23 * a1;

    ks = 1;
  }
//
//  Deal with a 0 input value of X.
//
  if ( *x == 0 )
  {
    *x = 314159265.0;
  }
//
//  Deal somewhat arbitrarily with negative input X.
//
  if ( *x < 0 )
  {
    *x = - ( *x );
  }
//
//  Break X into two parts X1 and X2 such that:
//
//    X = 2^23 * X1 + X2, 
//
//  then compute
//
//    Z = A1 * X2 + A2 * X1  (mod 2^23)
//    X = 2^23 * Z + A2 * X2  (mod 2^46).
//
  t1 = r23 * *x;
  x1 = ( double ) ( int ) t1;
  x2 = *x - t23 * x1;

  t1 = a1 * x2 + a2 * x1;
  t2 = ( double ) ( int ) ( r23 * t1 );
  z = t1 - t23 * t2;

  t3 = t23 * z + a2 * x2;
  t4 = ( double ) ( int ) ( r46 * t3 );
  *x = t3 - t46 * t4;

  value = r46 * ( *x );

  return value;
}
//****************************************************************************80

double randlc_jump ( double x, int k )

//****************************************************************************80
//
//  Purpose:
//
//    RANDLC_JUMP returns the K-th element of a uniform pseudorandom sequence.
//
//  Discussion:
//
//    The sequence uses the linear congruential generator:
//
//      X(K+1) = A * X(K)  mod 2^46
//
//    The K-th element, which can be represented as
//
//      X(K) = A^K * X(0)  mod 2^46
//
//    is computed directly using the binary algorithm for exponentiation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Bailey, Eric Barszcz, John Barton, D Browning, Robert Carter, 
//    Leonardo Dagum, Rod Fatoohi,
//    Samuel Fineberg, Paul Frederickson, Thomas Lasinski, Robert Schreiber, 
//    Horst Simon, V Venkatakrishnan, Sisira Weeratunga,
//    The NAS Parallel Benchmarks,
//    RNR Technical Report RNR-94-007,
//    March 1994.
//
//    Donald Knuth,
//    The Art of Computer Programming,
//    Volume 2, Seminumerical Algorithms,
//    Third Edition,
//    Addison Wesley, 1997,
//    ISBN: 0201896842,
//    LC: QA76.6.K64.
//
//  Parameters:
//
//    Input, double X, the initial seed (with index 0).  
//
//    Input, int K, the index of the desired value.
//
//    Output, double RANDLC_JUMP, the K-th value in the sequence.
//
{
  static double a = 1220703125.0;
  static double a1;
  static double a2;
  double b;
  double b1;
  double b2;
  int i;
  int j;
  int k2;
  static int ks = 0;
  int m;
  static double r23;
  static double r46;
  double t1;
  double t2;
  static double t23;
  double t3;
  double t4;
  static double t46;
  int twom;
  double value;
  double x1;
  double x2;
  double xk;
  double z;
//
//  If this is the first call, compute 
//
//    R23 = 2 ^ -23, 
//    R46 = 2 ^ -46,
//    T23 = 2 ^ 23, 
//    T46 = 2 ^ 46.  
//
//  These are computed in loops, rather than by merely using the power operator, 
//  in order to insure that the results are exact on all systems.  
//
  if ( ks == 0 ) 
  {
    r23 = 1.0;
    r46 = 1.0;
    t23 = 1.0;
    t46 = 1.0;
    
    for ( i = 1; i <= 23; i++ )
    {
      r23 = 0.50 * r23;
      t23 = 2.0 * t23;
    }
    for ( i = 1; i <= 46; i++ )
    {
      r46 = 0.50 * r46;
      t46 = 2.0 * t46;
    }
//
//  Break A into two parts such that A = 2^23 * A1 + A2.
//
    t1 = r23 * a;
    a1 = ( double ) ( int ) t1;
    a2 = a - t23 * a1;

    ks = 1;
  }

  if ( k < 0 )
  {
    cerr << "\n";
    cerr << "RANDLC_JUMP - Fatal error!\n";
    cerr << "  K < 0.\n";
    exit ( 1 );
  }
  else if ( k == 0 )
  {
    xk = x;
  }
//
//  Find M so that K < 2^M.
//
  else
  {
    k2 = k;
    xk = x;

    m = 1;
    twom = 2;
    while ( twom <= k )
    {
      twom = twom * 2;
      m = m + 1;
    }

    b = a;
    b1 = a1;
    b2 = a2;

    for ( i = 1; i <= m; i++ )
    {
      j = k2 / 2;
//
//  Replace X by A * X, if appropriate.
//
      if ( 2 * j != k2 )
      {
        t1 = r23 * xk;
        x1 = ( double ) ( ( int ) ( t1 ) );
        x2 = xk - t23 * x1;

        t1 = b1 * x2 + b2 * x1;
        t2 = ( double ) ( ( int ) ( r23 * t1 ) );
        z = t1 - t23 * t2;

        t3 = t23 * z + b2 * x2;
        t4 = ( double ) ( ( int ) ( r46 * t3 ) );
        xk = t3 - t46 * t4;
      }
//
//  Replace A by A * A mod 2^46.
//
      t1 = r23 * b;
      x1 = ( double ) ( ( int ) ( t1 ) );
      x2 = b - t23 * x1;

      t1 = b1 * x2 + b2 * x1;
      t2 = ( double ) ( ( int ) ( r23 * t1 ) );
      z = t1 - t23 * t2;

      t3 = t23 * z + b2 * x2;
      t4 = ( double ) ( ( int ) ( r46 * t3 ) );
      b = t3 - t46 * t4;
//
//  Update A1, A2.
//
      t1 = r23 * b;
      b1 = ( double ) ( ( int ) ( t1 ) );
      b2 = b - t23 * b1;

      k2 = j;
    }
  }
  value = r46 * xk;

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
