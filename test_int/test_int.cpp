# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_int.hpp"

//****************************************************************************80

char ch_cap ( char ch )

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
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 )
  {
    ch = ch - 32;
  }
  return ch;
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
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EULER_CONSTANT, the value of the Euler-Mascheroni constant.
//
{
  double value;

  value = 0.577215664901532860606512090082402431042;

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

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" I4.
//
{
  return 2147483647;
}
//****************************************************************************80

int i4_log_2 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
//
//  Example:
//
//        I  I4_LOG_10
//    -----  --------
//        0    0
//        1    0
//        2    1
//        3    1
//        4    2
//        5    2
//        7    2
//        8    3
//        9    3
//     1000    9
//     1024   10
//
//  Discussion:
//
//    I4_LOG_2 ( I ) + 1 is the number of binary digits in I.
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
//    Input, int I, the number whose logarithm base 2 is desired.
//
//    Output, int I4_LOG_2, the integer part of the logarithm base 2 of
//    the absolute value of X.
//
{
  int i_abs;
  int two_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    two_pow = 2;

    i_abs = abs ( i );

    while ( two_pow <= i_abs )
    {
      value = value + 1;
      two_pow = two_pow * 2;
    }
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

double *i4_to_halton_number_sequence_new ( int seed, int base, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_HALTON_NUMBER_SEQUENCE: next N elements of a scalar Halton sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, pages 84-90, 1960.
// 
//  Parameters:
//
//    Input, int SEED, the index of the desired element.
//    Only the absolute value of SEED is considered.
//    SEED = 0 is allowed, and returns R = 0.
//
//    Input, int BASE, the Halton base, which should 
//    be a prime number.  This routine only checks that BASE is greater 
//    than 1.
//
//    Input, int N, the number of elements desired.
//
//    Output, double I4_TO_HALTON_NUMBER_SEQUENCE[N], the SEED-th through 
//    (SEED+N-1)-th elements of the Halton sequence for base BASE.
//
{
  double base_inv;
  int digit;
  int i;
  double *r;
  int *seed2;
//
//  Set SEED2 = ( SEED, SEED+1, SEED+2, ..., SEED+N-1 )
//
  seed2 = i4vec_indicator_new ( n );

  for ( i = 0; i < n; i++ )
  {
    seed2[i] = seed2[i] + abs ( seed ) - 1;
  }

  if ( base <= 1 )
  {
    cerr << "\n";
    cerr << "I4_TO_HALTON_NUMBER_SEQUENCE - Fatal error!\n";
    cerr << "  The input base BASE is <= 1!\n";
    cerr << "  BASE = " << base << "\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    r[i] = 0.0;
    base_inv = 1.0 / ( double ) ( base );
    while ( seed2[i] != 0 )
    {
      digit = ( seed2[i] % base );
      r[i] = r[i] + ( double ) ( digit ) * base_inv;
      base_inv = base_inv / ( double ) ( base );
      seed2[i] = seed2[i] / base;
    }
  }

  delete [] seed2;

  return r;
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

double p00_even ( int prob, int int_num )

//****************************************************************************80
//
//  Purpose:
//
//    P00_EVEN uses evenly spaced points to integrate a function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int INT_NUM, the number of sample points.
//
//    Output, double P00_EVEN, the approximate integral.
//
{
  double a;
  double b;
  double *fx;
  double result;
  double *x;

  p00_lim ( prob, a, b );

  x = r8vec_linspace_new ( int_num, a, b );

  fx = p00_fun ( prob, int_num, x );

  result = ( b - a ) * r8vec_sum ( int_num, fx ) / ( double ) ( int_num );

  delete [] fx;
  delete [] x;

  return result;
}
//****************************************************************************80

double p00_exact ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    P00_EXACT returns the exact integral for any problem.
//
//  Discussion:
//
//    This routine provides a "generic" interface to the exact integral
//    routines for the various problems, and allows a problem to be called
//    by number (PROB) rather than by name.
//
//    In some cases, the "exact" value of the integral is in fact
//    merely a respectable approximation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the number of the desired test problem.
//
//    Output, double P00_EXACT, the exact value of the integral.
//
{
  double exact;

  if ( prob == 1 )
  {
    exact = p01_exact ( );
  }
  else if ( prob == 2 )
  {
    exact = p02_exact ( );
  }
  else if ( prob == 3 )
  {
    exact = p03_exact ( );
  }
  else if ( prob == 4 )
  {
    exact = p04_exact ( );
  }
  else if ( prob == 5 )
  {
    exact = p05_exact ( );
  }
  else if ( prob == 6 )
  {
    exact = p06_exact ( );
  }
  else if ( prob == 7 )
  {
    exact = p07_exact ( );
  }
  else if ( prob == 8 )
  {
    exact = p08_exact ( );
  }
  else if ( prob == 9 )
  {
    exact = p09_exact ( );
  }
  else if ( prob == 10 )
  {
    exact = p10_exact ( );
  }
  else if ( prob == 11 )
  {
    exact = p11_exact ( );
  }
  else if ( prob == 12 )
  {
    exact = p12_exact ( );
  }
  else if ( prob == 13 )
  {
    exact = p13_exact ( );
  }
  else if ( prob == 14 )
  {
    exact = p14_exact ( );
  }
  else if ( prob == 15 )
  {
    exact = p15_exact ( );
  }
  else if ( prob == 16 )
  {
    exact = p16_exact ( );
  }
  else if ( prob == 17 )
  {
    exact = p17_exact ( );
  }
  else if ( prob == 18 )
  {
    exact = p18_exact ( );
  }
  else if ( prob == 19 )
  {
    exact = p19_exact ( );
  }
  else if ( prob == 20 )
  {
    exact = p20_exact ( );
  }
  else if ( prob == 21 )
  {
    exact = p21_exact ( );
  }
  else if ( prob == 22 )
  {
    exact = p22_exact ( );
  }
  else if ( prob == 23 )
  {
    exact = p23_exact ( );
  }
  else if ( prob == 24 )
  {
    exact = p24_exact ( );
  }
  else if ( prob == 25 )
  {
    exact = p25_exact ( );
  }
  else if ( prob == 26 )
  {
    exact = p26_exact ( );
  }
  else if ( prob == 27 )
  {
    exact = p27_exact ( );
  }
  else if ( prob == 28 )
  {
    exact = p28_exact ( );
  }
  else if ( prob == 29 )
  {
    exact = p29_exact ( );
  }
  else if ( prob == 30 )
  {
    exact = p30_exact ( );
  }
  else if ( prob == 31 )
  {
    exact = p31_exact ( );
  }
  else if ( prob == 32 )
  {
    exact = p32_exact ( );
  }
  else if ( prob == 33 )
  {
    exact = p33_exact ( );
  }
  else if ( prob == 34 )
  {
    exact = p34_exact ( );
  }
  else if ( prob == 35 )
  {
    exact = p35_exact ( );
  }
  else if ( prob == 36 )
  {
    exact = p36_exact ( );
  }
  else if ( prob == 37 )
  {
    exact = p37_exact ( );
  }
  else if ( prob == 38 )
  {
    exact = p38_exact ( );
  }
  else if ( prob == 39 )
  {
    exact = p39_exact ( );
  }
  else if ( prob == 40 )
  {
    exact = p40_exact ( );
  }
  else if ( prob == 41 )
  {
    exact = p41_exact ( );
  }
  else if ( prob == 42 )
  {
    exact = p42_exact ( );
  }
  else if ( prob == 43 )
  {
    exact = p43_exact ( );
  }
  else if ( prob == 44 )
  {
    exact = p44_exact ( );
  }
  else if ( prob == 45 )
  {
    exact = p45_exact ( );
  }
  else if ( prob == 46 )
  {
    exact = p46_exact ( );
  }
  else if ( prob == 47 )
  {
    exact = p47_exact ( );
  }
  else if ( prob == 48 )
  {
    exact = p48_exact ( );
  }
  else if ( prob == 49 )
  {
    exact = p49_exact ( );
  }
  else if ( prob == 50 )
  {
    exact = p50_exact ( );
  }
  else if ( prob == 51 )
  {
    exact = p51_exact ( );
  }
  else if ( prob == 52 )
  {
    exact = p52_exact ( );
  }
  else if ( prob == 53 )
  {
    exact = p53_exact ( );
  }
  else if ( prob == 54 )
  {
    exact = p54_exact ( );
  }
  else if ( prob == 55 )
  {
    exact = p55_exact ( );
  }
  else if ( prob == 56 )
  {
    exact = p56_exact ( );
  }
  else if ( prob == 57 )
  {
    exact = p57_exact ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_EXACT - Fatal error!\n";
    cerr << "  Illegal problem number = " << prob << "\n";
    exit ( 1 );
  }
  return exact;
}
//****************************************************************************80

double *p00_fun ( int prob, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_FUN evaluates the integrand for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the number of the desired test problem.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  double *fx;

  if ( prob == 1 )
  {
    fx = p01_fun ( n, x );
  }
  else if ( prob == 2 )
  {
    fx = p02_fun ( n, x );
  }
  else if ( prob == 3 )
  {
    fx = p03_fun ( n, x );
  }
  else if ( prob == 4 )
  {
    fx = p04_fun ( n, x );
  }
  else if ( prob == 5 )
  {
    fx = p05_fun ( n, x );
  }
  else if ( prob == 6 )
  {
    fx = p06_fun ( n, x );
  }
  else if ( prob == 7 )
  {
    fx = p07_fun ( n, x );
  }
  else if ( prob == 8 )
  {
    fx = p08_fun ( n, x );
  }
  else if ( prob == 9 )
  {
    fx = p09_fun ( n, x );
  }
  else if ( prob == 10 )
  {
    fx = p10_fun ( n, x );
  }
  else if ( prob == 11 )
  {
    fx = p11_fun ( n, x );
  }
  else if ( prob == 12 )
  {
    fx = p12_fun ( n, x );
  }
  else if ( prob == 13 )
  {
    fx = p13_fun ( n, x );
  }
  else if ( prob == 14 )
  {
    fx = p14_fun ( n, x );
  }
  else if ( prob == 15 )
  {
    fx = p15_fun ( n, x );
  }
  else if ( prob == 16 )
  {
    fx = p16_fun ( n, x );
  }
  else if ( prob == 17 )
  {
    fx = p17_fun ( n, x );
  }
  else if ( prob == 18 )
  {
    fx = p18_fun ( n, x );
  }
  else if ( prob == 19 )
  {
    fx = p19_fun ( n, x );
  }
  else if ( prob == 20 )
  {
    fx = p20_fun ( n, x );
  }
  else if ( prob == 21 )
  {
    fx = p21_fun ( n, x );
  }
  else if ( prob == 22 )
  {
    fx = p22_fun ( n, x );
  }
  else if ( prob == 23 )
  {
    fx = p23_fun ( n, x );
  }
  else if ( prob == 24 )
  {
    fx = p24_fun ( n, x );
  }
  else if ( prob == 25 )
  {
    fx = p25_fun ( n, x );
  }
  else if ( prob == 26 )
  {
    fx = p26_fun ( n, x );
  }
  else if ( prob == 27 )
  {
    fx = p27_fun ( n, x );
  }
  else if ( prob == 28 )
  {
    fx = p28_fun ( n, x );
  }
  else if ( prob == 29 )
  {
    fx = p29_fun ( n, x );
  }
  else if ( prob == 30 )
  {
    fx = p30_fun ( n, x );
  }
  else if ( prob == 31 )
  {
    fx = p31_fun ( n, x );
  }
  else if ( prob == 32 )
  {
    fx = p32_fun ( n, x );
  }
  else if ( prob == 33 )
  {
    fx = p33_fun ( n, x );
  }
  else if ( prob == 34 )
  {
    fx = p34_fun ( n, x );
  }
  else if ( prob == 35 )
  {
    fx = p35_fun ( n, x );
  }
  else if ( prob == 36 )
  {
    fx = p36_fun ( n, x );
  }
  else if ( prob == 37 )
  {
    fx = p37_fun ( n, x );
  }
  else if ( prob == 38 )
  {
    fx = p38_fun ( n, x );
  }
  else if ( prob == 39 )
  {
    fx = p39_fun ( n, x );
  }
  else if ( prob == 40 )
  {
    fx = p40_fun ( n, x );
  }
  else if ( prob == 41 )
  {
    fx = p41_fun ( n, x );
  }
  else if ( prob == 42 )
  {
    fx = p42_fun ( n, x );
  }
  else if ( prob == 43 )
  {
    fx = p43_fun ( n, x );
  }
  else if ( prob == 44 )
  {
    fx = p44_fun ( n, x );
  }
  else if ( prob == 45 )
  {
    fx = p45_fun ( n, x );
  }
  else if ( prob == 46 )
  {
    fx = p46_fun ( n, x );
  }
  else if ( prob == 47 )
  {
    fx = p47_fun ( n, x );
  }
  else if ( prob == 48 )
  {
    fx = p48_fun ( n, x );
  }
  else if ( prob == 49 )
  {
    fx = p49_fun ( n, x );
  }
  else if ( prob == 50 )
  {
    fx = p50_fun ( n, x );
  }
  else if ( prob == 51 )
  {
    fx = p51_fun ( n, x );
  }
  else if ( prob == 52 )
  {
    fx = p52_fun ( n, x );
  }
  else if ( prob == 53 )
  {
    fx = p53_fun ( n, x );
  }
  else if ( prob == 54 )
  {
    fx = p54_fun ( n, x );
  }
  else if ( prob == 55 )
  {
    fx = p55_fun ( n, x );
  }
  else if ( prob == 56 )
  {
    fx = p56_fun ( n, x );
  }
  else if ( prob == 57 )
  {
    fx = p57_fun ( n, x );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_FUN - Fatal error!\n";
    cerr << "  Illegal problem number = " << prob << "\n";
    exit ( 1 );
  }
  return fx;
}
//****************************************************************************80

double p00_gauss_legendre ( int prob, int int_num )

//****************************************************************************80
//
//  Purpose:
//
//    P00_GAUSS_LEGENDRE applies a composite Gauss-Legendre rule.
//
//  Discussion:
//
//    A 4 point rule is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int INT_NUM, the number of subintervals.
//
//    Output, double P00_GAUSS_LEGENDRE, the approximate integral.
//
{
# define GAUSS_NUM 4

  double a;
  double a_sub;
  double b;
  double b_sub;
  double *fx;
  double gauss_abs[GAUSS_NUM] = {
    -0.861136311594052575223946488893, 
    -0.339981043584856264802665759103, 
     0.339981043584856264802665759103, 
     0.861136311594052575223946488893 };
  double gauss_weight[GAUSS_NUM] = {
    0.347854845137453857373063949222, 
    0.652145154862546142626936050778, 
    0.652145154862546142626936050778, 
    0.347854845137453857373063949222 };
  double h;
  int i;
  int int_i;
  double result;
  double x[GAUSS_NUM];

  p00_lim ( prob, a, b );

  h = ( b - a ) / ( double ) ( int_num );

  result = 0.0;

  for ( int_i = 1; int_i <= int_num; int_i++ )
  {
    a_sub = ( ( double ) ( int_num - int_i + 1 ) * a   
            + ( double ) (           int_i - 1 ) * b ) 
            / ( double ) ( int_num             );

    b_sub = ( ( double ) ( int_num - int_i ) * a   
            + ( double ) (           int_i ) * b ) 
            / ( double ) ( int_num         );

    for ( i = 0; i < GAUSS_NUM; i++ )
    {
      x[i] = 0.5 * ( ( 1.0 - gauss_abs[i] ) * a_sub 
                   + ( 1.0 + gauss_abs[i] ) * b_sub );
    }
    fx = p00_fun ( prob, GAUSS_NUM, x );

    result = result 
      + 0.5 * h * r8vec_dot_product ( GAUSS_NUM, gauss_weight, fx );

    delete [] fx;
  }
  return result;
# undef GAUSS_NUM
}
//****************************************************************************80

double p00_halton ( int prob, int int_num )

//****************************************************************************80
//
//  Purpose:
//
//    P00_HALTON applies a Halton sequence rule to integrate a function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, pages 84-90, 1960.
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int INT_NUM, the number of sample points.
//
//    Output, double P00_HALTON, the approximate integral.
//
{
  double a;
  double b;
  int base;
  double *fx;
  int i;
  double result;
  int seed;
  double *x;

  p00_lim ( prob, a, b );

  seed = 1;
  base = 2;
  x = i4_to_halton_number_sequence_new ( seed, base, int_num );

  for ( i = 0; i < int_num; i++ )
  {
    x[i] = a + ( b - a ) * x[i];
  }
  fx = p00_fun ( prob, int_num, x );

  result = ( b - a ) * r8vec_sum ( int_num, fx ) / ( double ) ( int_num );

  delete [] fx;
  delete [] x;

  return result;
}
//****************************************************************************80

void p00_lim ( int prob, double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P00_LIM returns the integration limits for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the number of the desired test problem.
//
//    Output, double &A, &B, the limits of integration.
//
{
  if ( prob == 1 )
  {
    p01_lim ( a, b );
  }
  else if ( prob == 2 )
  {
    p02_lim ( a, b );
  }
  else if ( prob == 3 )
  {
    p03_lim ( a, b );
  }
  else if ( prob == 4 )
  {
    p04_lim ( a, b );
  }
  else if ( prob == 5 )
  {
    p05_lim ( a, b );
  }
  else if ( prob == 6 )
  {
    p06_lim ( a, b );
  }
  else if ( prob == 7 )
  {
    p07_lim ( a, b );
  }
  else if ( prob == 8 )
  {
    p08_lim ( a, b );
  }
  else if ( prob == 9 )
  {
    p09_lim ( a, b );
  }
  else if ( prob == 10 )
  {
    p10_lim ( a, b );
  }
  else if ( prob == 11 )
  {
    p11_lim ( a, b );
  }
  else if ( prob == 12 )
  {
    p12_lim ( a, b );
  }
  else if ( prob == 13 )
  {
    p13_lim ( a, b );
  }
  else if ( prob == 14 )
  {
    p14_lim ( a, b );
  }
  else if ( prob == 15 )
  {
    p15_lim ( a, b );
  }
  else if ( prob == 16 )
  {
    p16_lim ( a, b );
  }
  else if ( prob == 17 )
  {
    p17_lim ( a, b );
  }
  else if ( prob == 18 )
  {
    p18_lim ( a, b );
  }
  else if ( prob == 19 )
  {
    p19_lim ( a, b );
  }
  else if ( prob == 20 )
  {
    p20_lim ( a, b );
  }
  else if ( prob == 21 )
  {
    p21_lim ( a, b );
  }
  else if ( prob == 22 )
  {
    p22_lim ( a, b );
  }
  else if ( prob == 23 )
  {
    p23_lim ( a, b );
  }
  else if ( prob == 24 )
  {
    p24_lim ( a, b );
  }
  else if ( prob == 25 )
  {
    p25_lim ( a, b );
  }
  else if ( prob == 26 )
  {
    p26_lim ( a, b );
  }
  else if ( prob == 27 )
  {
    p27_lim ( a, b );
  }
  else if ( prob == 28 )
  {
    p28_lim ( a, b );
  }
  else if ( prob == 29 )
  {
    p29_lim ( a, b );
  }
  else if ( prob == 30 )
  {
    p30_lim ( a, b );
  }
  else if ( prob == 31 )
  {
    p31_lim ( a, b );
  }
  else if ( prob == 32 )
  {
    p32_lim ( a, b );
  }
  else if ( prob == 33 )
  {
    p33_lim ( a, b );
  }
  else if ( prob == 34 )
  {
    p34_lim ( a, b );
  }
  else if ( prob == 35 )
  {
    p35_lim ( a, b );
  }
  else if ( prob == 36 )
  {
    p36_lim ( a, b );
  }
  else if ( prob == 37 )
  {
    p37_lim ( a, b );
  }
  else if ( prob == 38 )
  {
    p38_lim ( a, b );
  }
  else if ( prob == 39 )
  {
    p39_lim ( a, b );
  }
  else if ( prob == 40 )
  {
    p40_lim ( a, b );
  }
  else if ( prob == 41 )
  {
    p41_lim ( a, b );
  }
  else if ( prob == 42 )
  {
    p42_lim ( a, b );
  }
  else if ( prob == 43 )
  {
    p43_lim ( a, b );
  }
  else if ( prob == 44 )
  {
    p44_lim ( a, b );
  }
  else if ( prob == 45 )
  {
    p45_lim ( a, b );
  }
  else if ( prob == 46 )
  {
    p46_lim ( a, b );
  }
  else if ( prob == 47 )
  {
    p47_lim ( a, b );
  }
  else if ( prob == 48 )
  {
    p48_lim ( a, b );
  }
  else if ( prob == 49 )
  {
    p49_lim ( a, b );
  }
  else if ( prob == 50 )
  {
    p50_lim ( a, b );
  }
  else if ( prob == 51 )
  {
    p51_lim ( a, b );
  }
  else if ( prob == 52 )
  {
    p52_lim ( a, b );
  }
  else if ( prob == 53 )
  {
    p53_lim ( a, b );
  }
  else if ( prob == 54 )
  {
    p54_lim ( a, b );
  }
  else if ( prob == 55 )
  {
    p55_lim ( a, b );
  }
  else if ( prob == 56 )
  {
    p56_lim ( a, b );
  }
  else if ( prob == 57 )
  {
    p57_lim ( a, b );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_LIM - Fatal error!\n";
    cerr << "  Illegal problem number = " << prob << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

double p00_midpoint ( int prob, int int_num )

//****************************************************************************80
//
//  Purpose:
//
//    P00_MIDPOINT applies the composite midpoint rule to integrate a function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int INT_NUM, the number of subintervals.
//
//    Output, double P00_MIDPOINT, the approximate integral.
//
{
  double a;
  double a_sub;
  double b;
  double b_sub;
  double *fx;
  double h;
  int int_i;
  double result;
  double *x;

  p00_lim ( prob, a, b );

  h = b - a;

  x = new double[int_num];

  for ( int_i = 0; int_i < int_num; int_i++ )
  {
    a_sub = ( ( double ) ( int_num - int_i     ) * a   
            + ( double ) (           int_i     ) * b ) 
            / ( double ) ( int_num             );

    b_sub = ( ( double ) ( int_num - int_i - 1 ) * a   
            + ( double ) (           int_i + 1 ) * b ) 
            / ( double ) ( int_num         );

    x[int_i] = 0.5 * ( a_sub + b_sub );
  }

  fx = p00_fun ( prob, int_num, x );

  result = h * r8vec_sum ( int_num, fx ) / ( double ) ( int_num );

  delete [] fx;
  delete [] x;

  return result;
}
//****************************************************************************80

double p00_montecarlo ( int prob, int int_num )

//****************************************************************************80
//
//  Purpose:
//
//    P00_MONTECARLO applies the Monte Carlo rule to integrate a function.
//
//  Discussion:
//
//    This routine originally used an automatic array for X.  However,
//    under the G95 compiler, this was causing bizarre errors.  Replacing
//    the automatic array by an allocatable array made the problems
//    disappear.  Not an entirely satisfactory conclusion//
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int INT_NUM, the number of sample points.
//
//    Output, double P00_MONTECARLO, the approximate integral.
//
{
  double a;
  double b;
  double *fx;
  int i;
  double result;
  int seed;
  double *x;

  p00_lim ( prob, a, b );

  seed = get_seed ( );

  x = r8vec_uniform_01_new ( int_num, &seed );

  for ( i = 0; i < int_num; i++ )
  {
    x[i] = a + ( b - a ) * x[i];
  }
  fx = p00_fun ( prob, int_num, x );

  result = ( b - a ) * r8vec_sum ( int_num, fx ) / ( double ) ( int_num );

  delete [] fx;
  delete [] x;

  return result;
}
//****************************************************************************80

int p00_prob_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P00_PROB_NUM returns the number of test integration problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P00_PROB_NUM, the number of test integration
//    problems.
//
{
  int value;

  value = 57;

  return value;
}
//****************************************************************************80

double p00_simpson ( int prob, int int_num )

//****************************************************************************80
//
//  Purpose:
//
//    P00_SIMPSON applies the composite Simpson rule to integrate a function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int INT_NUM, the number of subintervals.
//
//    Output, double P00_SIMPSON, the approximate integral.
//
{
  double a;
  double b;
  double *fx;
  double h;
  int i;
  double result;
  double sum2;
  double sum4;
  double *x;

  p00_lim ( prob, a, b );

  h = b - a;

  x = r8vec_linspace_new ( 2 * int_num + 1, a, b );
 
  fx = p00_fun ( prob, 2 * int_num + 1, x );

  sum4 = 0.0;
  for ( i = 1; i < 2 * int_num + 1; i = i + 2 )
  {
    sum4 = sum4 + 4.0 * fx[i];
  }
  sum2 = 0.0;
  for ( i = 2; i < 2 * int_num; i = i + 2 )
  {
    sum2 = sum2 + 2.0 * fx[i];
  }
  result = fx[0] + fx[2*int_num] + sum4 + sum2;

  result = h * result / 6.0 / ( double ) int_num;

  delete [] fx;
  delete [] x;

  return result;
}
//****************************************************************************80

double p00_trapezoid ( int prob, int int_num )

//****************************************************************************80
//
//  Purpose:
//
//    P00_TRAPEZOID applies the composite trapezoid rule to integrate a function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int INT_NUM, the number of subintervals.
//
//    Output, double P00_TRAPEZOID, the approximate integral.
//
{
  double a;
  double b;
  double *fx;
  double result;
  double *x;

  p00_lim ( prob, a, b );

  x = r8vec_linspace_new ( int_num + 1, a, b );

  fx = p00_fun ( prob, int_num + 1, x );

  result = ( b - a ) * ( 
           0.5 * fx[0] 
         + r8vec_sum ( int_num - 1, fx + 1 )
         + 0.5 * fx[int_num] ) / ( double ) ( int_num );

  delete [] fx;
  delete [] x;

  return result;
}
//****************************************************************************80

double p01_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_EXACT returns the exact integral for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P01_EXACT, the value of the integral.
//
{
  double exact;

  exact = exp ( 1.0 ) - 1.0;

  return exact;
}
//****************************************************************************80

double *p01_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_FUN evaluates the integrand for problem 1.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    exp ( x )
//
//  Antiderivative:
//
//    exp ( x )
//
//  Exact Integral:
//
//    exp ( 1 ) - 1
//
//  Approximate Integral (25 digits):
//
//    1.718281828459045235360287...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double FX[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = exp ( x[i] );
  }
  return fx;
}
//****************************************************************************80

void p01_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P01_LIM returns the integration limits for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p02_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_EXACT returns the exact integral for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P02_EXACT, the value of the integral.
//
{
  double exact;

  exact = 0.7;

  return exact;
}
//****************************************************************************80

double *p02_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_FUN evaluates the integrand for problem 2.
//
//  Discussion:
//
//    The integrand is discontinuous at X = 0.3.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    if ( x < 0.3 )
//      f(x) = 0
//    else
//      f(x) = 1
//
//  Antiderivative:
//
//    if ( x < 0.3 )
//      g(x) = 0
//    else
//      g(x) = x - 0.3
//
//  Exact Integral:
//
//    0.7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P02_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < 0.3 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = 1.0;
    }
  }
  return fx;
}
//****************************************************************************80

void p02_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P02_LIM returns the integration limits for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p03_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_EXACT returns the exact integral for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P03_EXACT, the value of the integral.
//
{
  double exact;

  exact = 2.0 / 3.0;

  return exact;
}
//****************************************************************************80

double *p03_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_FUN evaluates the integrand for problem 3.
//
//  Discussion:
//
//    The integrand is not differentiable at X = 0.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    sqrt ( x )
//
//  Antiderivative:
//
//    ( 2 / 3 ) * x^(3/2)
//
//  Exact Integral:
//
//    2 / 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P03_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = sqrt ( x[i] );
  }
  return fx;
}
//****************************************************************************80

void p03_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P03_LIM returns the integration limits for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p04_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_EXACT returns the estimated integral for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P04_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.47942822668880166736;

  return exact;
}
//****************************************************************************80

double *p04_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_FUN evaluates the integrand for problem 4.
//
//  Interval:
//
//    -1 <= x <= 1
//
//  Integrand:
//
//    0.92 * cosh ( x ) - cos ( x )
//
//  Antiderivative:
//
//    0.92 * sinh ( x ) - sin ( x )
//
//  Exact Integral:
//
//    1.84 * sinh ( 1 ) - 2 * sin ( 1 )
//
//  Approximate Integral (20 digits):
//
//    0.47942822668880166736...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Clenshaw, Alan Curtis,
//    A Method for Numerical Integration on an Automatic Computer,
//    Numerische Mathematik,
//    Volume 2, Number 1, December 1960, pages 197-205.
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P04_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 0.92 * cosh ( x[i] ) - cos ( x[i] );
  }
  return fx;
}
//****************************************************************************80

void p04_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P04_LIM returns the integration limits for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = -1.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p05_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_EXACT returns the estimated integral for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P05_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 1.5822329637296729331;

  return exact;
}
//****************************************************************************80

double *p05_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_FUN evaluates the integrand for problem 5.
//
//  Interval:
//
//    -1 <= x <= 1
//
//  Integrand:
//
//    1 / ( x^4 + x^2 + 0.9 )
//
//  Approximate Integral (20 digits):
//
//    1.5822329637296729331...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Clenshaw, Alan Curtis,
//    A Method for Numerical Integration on an Automatic Computer,
//    Numerische Mathematik,
//    Volume 2, Number 1, December 1960, pages 197-205.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P05_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0 / ( pow ( x[i], 4 ) + pow ( x[i], 2 ) + 0.9 );
  }
  return fx;
}
//****************************************************************************80

void p05_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P05_LIM returns the integration limits for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = -1.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p06_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_EXACT returns the exact integral for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P06_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 1.460447131787105;

  return exact;
}
//****************************************************************************80

double *p06_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_FUN evaluates the integrand for problem 6.
//
//  Interval:
//
//    -1 <= x <= 1
//
//  Integrand:
//
//    sqrt ( abs ( x + 0.5 ) )
//
//  Exact Integral:
//
//    ( sqrt ( 2 ) + 3 * sqrt ( 6 ) ) / 6 = 1.460447131787105
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Clenshaw, Alan Curtis,
//    A Method for Numerical Integration on an Automatic Computer,
//    Numerische Mathematik,
//    Volume 2, Number 1, December 1960, pages 197-205.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P06_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = sqrt ( abs ( x[i] + 0.5 ) );
  }
  return fx;
}
//****************************************************************************80

void p06_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P06_LIM returns the integration limits for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = -1.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p07_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_EXACT returns the exact integral for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P07_EXACT, the value of the integral.
//
{
  double exact;

  exact = 2.0;

  return exact;
}
//****************************************************************************80

double *p07_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_FUN evaluates the integrand for problem 7.
//
//  Discussion:
//
//    The integrand is singular at x = 0.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    1 / sqrt ( x )
//
//  Antiderivative:
//
//    2 * sqrt ( x )
//
//  Exact Integral:
//
//    2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P07_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < x[i] )
    {
      fx[i] = 1.0 / sqrt ( x[i] );
    }
    else
    {
      fx[i] = 0.0;
    }
  }
  return fx;
}
//****************************************************************************80

void p07_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P07_LIM returns the integration limits for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p08_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_EXACT returns the estimated integral for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P08_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.86697298733991103757;

  return exact;
}
//****************************************************************************80

double *p08_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_FUN evaluates the integrand for problem 8.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    1 / ( 1 + x^4 )
//
//  Antiderivative:
//
//    (1/8) * sqrt ( 2 )
//    * ln ( ( x^2 + sqrt ( 2 ) * x + 1 ) / ( x^2 - sqrt ( 2 ) * x + 1 ) )
//    + (1/4) * sqrt ( 2 ) * arctan ( sqrt ( 2 ) * x / ( 1 - x^2 ) )
//
//  Approximate Integral (20 digits):
//
//    0.86697298733991103757...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P08_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0 / ( 1.0 + pow ( x[i], 4 ) );
  }
  return fx;
}
//****************************************************************************80

void p08_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P08_LIM returns the integration limits for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p09_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_EXACT returns the estimated integral for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P09_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 1.1547005383792515290;

  return exact;
}
//****************************************************************************80

double *p09_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P09_FUN evaluates the integrand for problem 9.
//
//  Discussion:
//
//    The integrand is oscillatory, going through 5 periods in [0,1].
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    2 / ( 2 + sin ( 10 * pi * x ) )
//
//  Antiderivative:
//
//    1 / ( 5 * pi * sqrt ( 3 ) ) *
//    arctan ( ( 1 + 2 * tan ( 5 * pi * x ) ) / sqrt ( 3 ) )
//
//  Exact Integral:
//
//    2 / sqrt ( 3 )
//
//  Approximate Integral (20 digits):
//
//    1.1547005383792515290...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P09_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  static double pi = 3.141592653589793;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 2.0 / ( 2.0 + sin ( 10.0 * pi * x[i] ) );
  }
  return fx;
}
//****************************************************************************80

void p09_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P09_LIM returns the integration limits for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p10_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_EXACT returns the estimated integral for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P10_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.6931471805599453094172321;

  return exact;
}
//****************************************************************************80

double *p10_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P10_FUN evaluates the integrand for problem 10.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    1 / ( 1 + x )
//
//  Antiderivative:
//
//    ln ( 1 + x )
//
//  Exact Integral:
//
//    ln ( 2 )
//
//  Approximate Integral (25 digits):
//
//    0.6931471805599453094172321...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P10_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0 / ( 1.0 + x[i] );
  }
  return fx;
}
//****************************************************************************80

void p10_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P10_LIM returns the integration limits for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p11_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P11_EXACT returns the estimated integral for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P11_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.37988549304172247537;

  return exact;
}
//****************************************************************************80

double *p11_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P11_FUN evaluates the integrand for problem 11.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    1 / ( 1 + exp ( x ) )
//
//  Antiderivative:
//
//    ln ( exp ( x ) / ( 1 + exp ( x ) ) )
//
//  Exact Integral:
//
//    ln ( 2 * e / ( 1 + e ) )
//
//  Approximate Integral (20 digits):
//
//    0.37988549304172247537...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P11_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0 / ( 1.0 + exp ( x[i] ) );
  }
  return fx;
}
//****************************************************************************80

void p11_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P11_LIM returns the integration limits for problem 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p12_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P12_EXACT returns the estimated integral for problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P12_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.77750463411224827642;

  return exact;
}
//****************************************************************************80

double *p12_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P12_FUN evaluates the integrand for problem 12.
//
//  Discussion:
//
//    The integrand has a removable singularity at x = 0.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    x / ( exp ( x ) - 1 )
//
//  Antiderivative:
//
//    The Debye function.
//
//  Approximate Integral (20 digits):
//
//    0.77750463411224827642...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P12_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 1.0;
    }
    else
    {
      fx[i] = x[i] / ( exp ( x[i] ) - 1.0 );
    }
  }
  return fx;
}
//****************************************************************************80

void p12_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P12_LIM returns the integration limits for problem 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p13_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P13_EXACT returns the estimated integral for problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P13_EXACT, the estimated value of the integral.
//
{
  double a;
  double b;
  double exact;
 
  p13_lim ( a, b );

  exact = r8_si ( b ) - r8_si ( a );

  return exact;
}
//****************************************************************************80

double *p13_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P13_FUN evaluates the integrand for problem 13.
//
//  Interval:
//
//    0 <= x <= 10
//
//  Integrand:
//
//    sin ( x ) / x
//
//  Approximate Integral (20 digits):
//
//    1.6583475942188740493...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P13_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 1.0;
    }
    else
    {
      fx[i] = sin ( x[i] ) / x[i];
    }
  }
  return fx;
}
//****************************************************************************80

void p13_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P13_LIM returns the integration limits for problem 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 10.0;

  return;
}
//****************************************************************************80

double p14_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P14_EXACT returns the estimated integral for problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P14_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.500000211166;

  return exact;
}
//****************************************************************************80

double *p14_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P14_FUN evaluates the integrand for problem 14.
//
//  Discussion:
//
//    For X's that aren't actually very big, the function becomes very 
//    small.  Some compilers may product code that fails in these cases.
//    An attempt has been made to return a value of 0 when the computed
//    value of F(X) would be extremely small.
//
//  Interval:
//
//    0 <= x <= 10
//
//  Integrand:
//
//    sqrt ( 50 ) * exp ( - 50 * pi * x * x )
//
//  Exact Integral:
//
//    0.5 * erf ( 50 * sqrt ( 2 * pi ) )
//
//  Approximate Integral:
//
//    0.500000211166...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P14_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  static double pi = 3.141592653589793;
  static double x_max = 0.0;

  if ( x_max == 0.0 )
  {
    x_max = sqrt ( log ( r8_max ( r8_epsilon ( ), 1.0E-10 ) ) 
      / ( - 50.0 * pi ) );
  }

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {

    if ( x_max < r8_abs ( x[i] ) )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = sqrt ( 50.0 ) * exp ( - 50.0 * pi * x[i] * x[i] );
    }
  }
  return fx;
}
//****************************************************************************80

void p14_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P14_LIM returns the integration limits for problem 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 10.0;

  return;
}
//****************************************************************************80

double p15_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P15_EXACT returns the exact integral for problem 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P15_EXACT, the value of the integral.
//
{
  double exact;

  exact = 1.0;

  return exact;
}
//****************************************************************************80

double *p15_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P15_FUN evaluates the integrand for problem 15.
//
//  Interval:
//
//    0 <= x <= 10
//
//  Integrand:
//
//    25 * exp ( - 25 * x )
//
//  Antiderivative:
//
//    - exp ( - 25 * x )
//
//  Exact Integral:
//
//    1 - exp ( - 250 )
//
//  Approximate Integral:
//
//    1.00000000...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P15_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 25.0 * exp ( - 25.0 * x[i] );
  }
  return fx;
}
//****************************************************************************80

void p15_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P15_LIM returns the integration limits for problem 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 10.0;

  return;
}
//****************************************************************************80

double p16_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P16_EXACT returns the exact integral for problem 16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P16_EXACT, the value of the integral.
//
{
  double exact;

  exact = 0.49936338107645674464;

  return exact;
}
//****************************************************************************80

double *p16_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P16_FUN evaluates the integrand for problem 16.
//
//  Interval:
//
//    0 <= x <= 10
//
//  Integrand:
//
//    50.0 / ( pi * ( 2500.0 * x * x + 1.0 ) )
//
//  Antiderivative:
//
//    ( 1 / pi ) * arctan ( 50 * x )
//
//  Approximate Integral (20 digits):
//
//    0.49936338107645674464...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P16_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  static double pi = 3.141592653589793;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 50.0 / pi / ( 2500.0 * x[i] * x[i] + 1.0 );
  }
  return fx;
}
//****************************************************************************80

void p16_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P16_LIM returns the integration limits for problem 16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p17_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P17_EXACT returns the estimated integral for problem 17.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P17_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.5;

  return exact;
}
//****************************************************************************80

double *p17_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P17_FUN evaluates the integrand for problem 17.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    ( sin ( 50 * pi * x ) )^2
//
//  Antiderivative:
//
//    1/2 x - sin ( 100 * pi * x ) / ( 200 * pi )
//
//  Approximate Integral:
//
//    0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P17_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  static double pi = 3.141592653589793;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = pow ( sin ( 50.0 * pi * x[i] ), 2 );
  }
  return fx;
}
//****************************************************************************80

void p17_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P17_LIM returns the integration limits for problem 17.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p18_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P18_EXACT returns the estimated integral for problem 18.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P18_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.17055734950243820437;

  return exact;
}
//****************************************************************************80

double *p18_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P18_FUN evaluates the integrand for problem 18.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    x / ( exp ( x ) + 1 )
//
//  Approximate Integral (20 digits):
//
//    0.17055734950243820437...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Hermann Engels,
//    Numerical Quadrature and Cubature,
//    Academic Press, 1980.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P18_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = x[i] / ( exp ( x[i] ) + 1.0 );
  }
  return fx;
}
//****************************************************************************80

void p18_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P18_LIM returns the integration limits for problem 18.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p19_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P19_EXACT returns the exact integral for problem 19.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P19_EXACT, the value of the integral.
//
{
  double exact;

  exact = - 1.0;

  return exact;
}
//****************************************************************************80

double *p19_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P19_FUN evaluates the integrand for problem 19.
//
//  Discussion:
//
//    The integrand is singular at x = 0.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    ln ( x )
//
//  Antiderivative:
//
//    x * ln ( x ) - x
//
//  Exact Integral:
//
//    -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P19_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] <= 1.0E-15 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = log ( x[i] );
    }
  }
  return fx;
}
//****************************************************************************80

void p19_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P19_LIM returns the integration limits for problem 19.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p20_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P20_EXACT returns the estimated integral for problem 20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P20_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 1.5643964440690497731;

  return exact;
}
//****************************************************************************80

double *p20_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P20_FUN evaluates the integrand for problem 20.
//
//  Interval:
//
//    -1 <= x <= 1
//
//  Integrand:
//
//    1 / ( x^2 + 1.005 )
//
//  Antiderivative:
//
//    ( 1 / sqrt ( 1.005 ) ) * arctan ( x / sqrt ( 1.005 ) )
//
//  Approximate Integral (20 digits):
//
//    1.5643964440690497731...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P20_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0 / ( pow ( x[i], 2 ) + 1.005 );
  }
  return fx;
}
//****************************************************************************80

void p20_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P20_LIM returns the integration limits for problem 20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = -1.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p21_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P21_EXACT returns the estimated integral for problem 21.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P21_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.21080273631018169851;

  return exact;
}
//****************************************************************************80

double *p21_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P21_FUN evaluates the integrand for problem 21.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//       ( sech (   10.0 * ( x - 0.2 ) ) )^2
//     + ( sech (  100.0 * ( x - 0.4 ) ) )^4
//     + ( sech ( 1000.0 * ( x - 0.6 ) ) )^6
//
//  Exact Integral:
//
//    ( 1 + tanh ( 8 ) * tanh ( 2 ) ) / 10.0 + 2 / 150 + 2 / 1875
//
//  Approximate Integral (20 digits):
//
//    0.21080273631018169851...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P21_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 
      pow ( r8_sech (   10.0 * ( x[i] - 0.2 ) ), 2 )
    + pow ( r8_sech (  100.0 * ( x[i] - 0.4 ) ), 4 )
    + pow ( r8_sech ( 1000.0 * ( x[i] - 0.6 ) ), 6 );
  }
  return fx;
}
//****************************************************************************80

void p21_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P21_LIM returns the integration limits for problem 21.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p22_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P22_EXACT returns the estimated integral for problem 22.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P22_EXACT, the estimated value of the integral.
//
{
  double exact;
  static double pi = 3.141592653589793;

  exact = 0.125 * log ( 9.0 ) + pi / sqrt ( 48.0 );

  return exact;
}
//****************************************************************************80

double *p22_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P22_FUN evaluates the integrand for problem 22.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    1 / ( x^4 + x^2 + 1 )
//
//  Exact integral:
//
//    ln ( 9 ) / 8 + pi / sqrt ( 48 )
//
//  Approximate Integral (20 digits):
//
//    0.72810291322558188550...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P22_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0 / ( pow ( x[i], 4 ) + pow ( x[i], 2 ) + 1.0 );
  }
  return fx;
}
//****************************************************************************80

void p22_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P22_LIM returns the integration limits for problem 22.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p23_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P23_EXACT returns the estimated integral for problem 23.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P23_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.62471325642771360429;

  return exact;
}
//****************************************************************************80

double *p23_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P23_FUN evaluates the integrand for problem 23.
//
//  Discussion:
//
//    The integrand has a singularity at X = 0.
//    The integrand is discontinuous at X = 0.
//    The integrand is arbitrarily oscillatory as X decreases to 0.
//    The integrand becomes unbounded as X decreases to 0.
//
//    Integral ( 0 < X < 1 ) ( 1 / X ) sin ( 1 / X ) dX
//    = Integral ( 1 < X < Infinity ) ( 1 / X ) * sin ( X ) dX.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    ( 1 / x ) sin ( 1 / x )
//
//  Approximate Integral (20 digits):
//
//    0.62471325642771360429...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P23_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  
  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = ( 1.0 / x[i] ) * sin ( 1.0 / x[i] );
    }
  }
  return fx;
}
//****************************************************************************80

void p23_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P23_LIM returns the integration limits for problem 23.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p24_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P24_EXACT returns the estimated integral for problem 24.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P24_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = - 0.0067547455;

  return exact;
}
//****************************************************************************80

double *p24_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P24_FUN evaluates the integrand for problem 24.
//
//  Discussion:
//
//    The integrand is continuous, but nowhere differentiable.
//
//  Interval:
//
//    0 <= X <= 0.5
//
//  Integrand:
//
//    ( 1 / pi ) * sum ( 1 <= I < Infinity ) 2^(-I) * cos ( 7^I * pi * X )
//
//  Approximate Integral:
//
//    - 0.0067547455
//
//  Antiderivative:
//
//    ( 1 / pi^2 ) * sum ( 1 <= I < Infinity ) 14^(-I) * sin ( 7^I * pi * X )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Herbert Salzer, Norman Levine,
//    Table of a Weierstrass Continuous Nondifferentiable Function,
//    Mathematics of Computation,
//    Volume 15, pages 120 - 130, 1961.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P24_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  int j;
  static int n_term = 40;
  static double pi = 3.141592653589793;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 0.0;
    for ( j = 1; j <= n_term; j++ )
    {
      fx[i] = fx[i] + cos ( pow ( 7.0, j ) * pi * x[i] ) / pow ( 2.0, j );
    }
    fx[i] = fx[i] / pi;
  }
  return fx;
}
//****************************************************************************80

void p24_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P24_LIM returns the integration limits for problem 24.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 0.5;

  return;
}
//****************************************************************************80

double p25_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P25_EXACT returns the estimated integral for problem 25.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P25_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.3 * log ( 0.3 ) + 0.7 * log ( 0.7 ) - 1.0;

  return exact;
}
//****************************************************************************80

double *p25_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P25_FUN evaluates the integrand for problem 25.
//
//  Interval:
//
//    0 <= X <= 1.
//
//  Integrand:
//
//    ln ( abs ( x - 0.7 ) )
//
//  Exact Integral:
//
//    0.3 * ln ( 0.3 ) + 0.7 * ln ( 0.7 ) - 1
//
//  Approximate Integral (20 digits):
//
//    -1.6108643020548934630
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kendall Atkinson,
//    An Introduction to Numerical Analysis,
//    Prentice Hall, 1984, page 303.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P25_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  
  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.7 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = log ( r8_abs ( x[i] - 0.7 ) );
    }
  }
  return fx;
}
//****************************************************************************80

void p25_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P25_LIM returns the integration limits for problem 25.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p26_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P26_EXACT returns the exact integral for problem 26.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P26_EXACT, the value of the integral.
//
{
  double exact;

  exact = 7.9549265210128452745;

  return exact;
}
//****************************************************************************80

double *p26_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P26_FUN evaluates the integrand for problem 26.
//
//  Interval:
//
//    0 <= x <= 2 pi
//
//  Integrand:
//
//    exp ( cos ( x ) )
//
//  Exact Integral:
//
//    2 * Pi * I0(1)
//
//  Approximate Integral (20 digits):
//
//    7.9549265210128452745...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kendall Atkinson,
//    An Introduction to Numerical Analysis,
//    Prentice Hall, 1984, page 262.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P26_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  
  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = exp ( cos ( x[i] ) );
  }
  return fx;
}
//****************************************************************************80

void p26_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P26_LIM returns the integration limits for problem 26.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  static double pi = 3.141592653589793;

  a = 0.0;
  b = 2.0 * pi;

  return;
}
//****************************************************************************80

double p27_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P27_EXACT returns the exact integral for problem 27.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P27_EXACT, the value of the integral.
//
{
  double exact;

  exact = 5.0 - 6.0 * log ( 2.0 );

  return exact;
}
//****************************************************************************80

double *p27_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P27_FUN evaluates the integrand for problem 27.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    1 / ( X^(1/2) + X^(1/3) )
//
//  Exact Integral:
//
//    5 - 6 * ln ( 2 )
//
//  Approximate Integral (20 digits):
//
//    0.84111691664032814350...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P27_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = 1.0 / ( sqrt ( x[i] ) + pow ( x[i], 1.0 / 3.0 ) );
    }
  }
  return fx;
}
//****************************************************************************80

void p27_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P27_LIM returns the integration limits for problem 27.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p28_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P28_EXACT returns the exact integral for problem 28.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P28_EXACT, the value of the integral.
//
{
  double exact;
  static double pi = 3.141592653589793;

  exact = ( 50.0 / 2501.0 ) * ( 1.0 - exp ( - 2.0 * pi ) );

  return exact;
}
//****************************************************************************80

double *p28_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P28_FUN evaluates the integrand for problem 28.
//
//  Interval:
//
//    0 <= X <= 2 PI
//
//  Integrand:
//
//    exp ( - X ) * sin ( 50 * X )
//
//  Exact Integral:
//
//    50 / ( 2501 ) * ( 1 - exp ( - 2 * PI ) )
//
//  Approximate Integral (20 digits):
//
//    0.019954669277654778312...
//
//  Reference:
//
//    Kendall Atkinson,
//    An Introduction to Numerical Analysis,
//    Prentice Hall, 1984, page 303.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P28_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = exp ( - x[i] ) * sin ( 50.0 * x[i] );
  }
  return fx;
}
//****************************************************************************80

void p28_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P28_LIM returns the integration limits for problem 28.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  static double pi = 3.141592653589793;

  a = 0.0;
  b = 2.0 * pi;

  return;
}
//****************************************************************************80

double p29_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P29_EXACT returns the exact integral for problem 29.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P29_EXACT, the value of the integral.
//
{
  double exact;

  exact = 1.0 - log ( 2.0 );

  return exact;
}
//****************************************************************************80

double *p29_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P29_FUN evaluates the integrand for problem 29.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    F ( X ) = 1 / ( X + 2 )   for 0 < X < E - 2
//            = 0               otherwise
//
//  Exact Integral:
//
//    1 - ln ( 2 )
//
//  Approximate Integral (20 digits):
//
//    0.30685281944005469058...
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P29_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 <= x[i] && x[i] <= exp ( 1.0 ) - 2.0 )
    {
      fx[i] = 1.0 / ( x[i] + 2.0 );
    }
    else
    {
      fx[i] = 0.0;
    }
  }
  return fx;
}
//****************************************************************************80

void p29_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P29_LIM returns the integration limits for problem 29.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p30_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P30_EXACT returns the exact integral for problem 30.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P30_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = -4.5275696251606720278;

  return exact;
}
//****************************************************************************80

double *p30_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P30_FUN evaluates the integrand for problem 30.
//
//  Interval:
//
//    2 <= x <= 7
//
//  Integrand:
//
//          cos (       x )
//    + 5 * cos ( 1.6 * x )
//    - 2 * cos ( 2.0 * x )
//    + 5 * cos ( 4.5 * x )
//    + 7 * cos ( 9.0 * x )
//
//  Antiderivative:
//
//          sin (       x )
//    + 5 * sin ( 1.6 * x ) / 1.6
//    - 2 * sin ( 2.0 * x ) / 2.0
//    + 5 * sin ( 4.5 * x ) / 4.5
//    + 7 * sin ( 9.0 * x ) / 9.0
//
//  Exact Integral:
//
//    -4.5275696251606720278
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dianne OLeary,
//    Scientific Computing with Case Studies,
//    SIAM, 2008,
//    ISBN13: 978-0-898716-66-5,
//    LC: QA401.O44.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P30_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 
              cos (       x[i] ) 
      + 5.0 * cos ( 1.6 * x[i] ) 
      - 2.0 * cos ( 2.0 * x[i] ) 
      + 5.0 * cos ( 4.5 * x[i] ) 
      + 7.0 * cos ( 9.0 * x[i] );
  }
  return fx;
}
//****************************************************************************80

void p30_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P30_LIM returns the integration limits for problem 30.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 2.0;
  b = 7.0;

  return;
}
//****************************************************************************80

double p31_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P31_EXACT returns the exact integral for problem 31.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P31_EXACT, the value of the integral.
//
{
  double exact;

  exact = 2.0 * atan ( 4.0 );

  return exact;
}
//****************************************************************************80

double *p31_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P31_FUN evaluates the integrand for problem 31.
//
//  Discussion:
//
//    A simple Newton-Cotes quadrature rule, in which the order of the
//    rule is increased, but the interval is not subdivided, diverges
//    for this integrand.
//
//    This is Runge's function, a standard example of the perils of
//    using high order polynomial interpolation at equally spaced nodes.
//    Since this is exactly what is really going on in a Newton Cotes
//    rule, it is little wonder that the result is so poor.
//
//  Interval:
//
//    -4 <= x <= 4
//
//  Integrand:
//
//    1 / ( 1 + x^2 )
//
//  Antiderivative:
//
//    arctan ( x )
//
//  Exact Integral:
//
//    2 * arctan ( 4 )
//
//  Approximate Integral (20 digits):
//
//    2.6516353273360649301...
//
//  Reference:
//
//    Kendall Atkinson,
//    An Introduction to Numerical Analysis,
//    Prentice Hall, 1984, page 266.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P31_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0 / ( 1.0 + x[i] * x[i] );
  }
  return fx;
}
//****************************************************************************80

void p31_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P31_LIM returns the integration limits for problem 31.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = - 4.0;
  b =   4.0;

  return;
}
//****************************************************************************80

double p32_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P32_EXACT returns the exact integral for problem 32.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P32_EXACT, the value of the integral.
//
{
  double exact;
  static double pi = 3.141592653589793;

  exact = - 0.5 * ( exp ( pi ) + 1.0 );

  return exact;
}
//****************************************************************************80

double *p32_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P32_FUN evaluates the integrand for problem 32.
//
//  Interval:
//
//    0 <= X <= PI
//
//  Integrand:
//
//    exp ( X ) * cos ( X )
//
//  Antiderivative:
//
//    0.5 * exp ( X ) * ( sin ( X ) + cos ( X ) )
//
//  Exact Integral:
//
//    - 0.5 * ( exp ( PI ) + 1 )
//
//  Approximate Integral (20 digits):
//
//    -12.070346316389634503...
//
//  Reference:
//
//    Kendall Atkinson,
//    An Introduction to Numerical Analysis,
//    Prentice Hall, 1984, page 254, 277, 297.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P32_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = exp ( x[i] ) * cos ( x[i] );
  }
  return fx;
}
//****************************************************************************80

void p32_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P32_LIM returns the integration limits for problem 32.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  static double pi = 3.141592653589793;

  a = 0.0;
  b = pi;

  return;
}
//****************************************************************************80

double p33_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P33_EXACT returns the exact integral for problem 33.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P33_EXACT, the value of the integral.
//
{
  double exact;
  static double pi = 3.141592653589793;

  exact = 0.5 * sqrt ( pi );

  return exact;
}
//****************************************************************************80

double *p33_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P33_FUN evaluates the integrand for problem 33.
//
//  Discussion:
//
//    The integrand is singular at both endpoints of the interval.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    sqrt ( - ln ( X ) )
//
//  Exact Integral:
//
//    sqrt ( pi ) / 2
//
//  Approximate Integral (20 digits):
//
//    0.88622692545275801365...
//
//  Reference:
//
//    Kendall Atkinson,
//    An Introduction to Numerical Analysis,
//    Prentice Hall, 1984, page 307.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P33_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] <= 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = sqrt ( - log ( x[i] ) );
    }
  }
  return fx;
}
//****************************************************************************80

void p33_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P33_LIM returns the integration limits for problem 33.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p34_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P34_EXACT returns the exact integral for problem 34.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P34_EXACT, the value of the integral.
//
{
  double exact;

  exact = 1627879.0 / 1500.0;

  return exact;
}
//****************************************************************************80

double *p34_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P34_FUN evaluates the integrand for problem 34.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    ( 10 * X - 1 ) * ( 10 * X - 1.1 ) * ( 10 * X - 1.2 ) * ( 10 * X - 1.3 )
//
//  Exact Integral:
//
//    1627879 / 1500
//
//  Approximate Integral (20 digits):
//
//    1085.2526666666666666...
//
//  Reference:
//
//    Hermann Engels,
//    Numerical Quadrature and Cubature,
//    Academic Press, 1980.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P34_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = ( 10.0 * x[i] - 1.0 ) * ( 10.0 * x[i] - 1.1 )
      * ( 10.0 * x[i] - 1.2 ) * ( 10.0 * x[i] - 1.3 );
  }
  return fx;
}
//****************************************************************************80

void p34_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P34_LIM returns the integration limits for problem 34.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p35_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P35_EXACT returns the exact integral for problem 35.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P35_EXACT, the value of the integral.
//
{
  double exact;

  exact = 26.0;

  return exact;
}
//****************************************************************************80

double *p35_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P35_FUN evaluates the integrand for problem 35.
//
//  Interval:
//
//    -9 <= X <= 100
//
//  Integrand:
//
//    1 / sqrt ( abs ( X ) )
//
//  Exact Integral:
//
//    26
//
//  Reference:
//
//    Hermann Engels,
//    Numerical Quadrature and Cubature,
//    Academic Press, 1980.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P35_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = 1.0 / sqrt ( r8_abs ( x[i] ) );
    }
  }
  return fx;
}
//****************************************************************************80

void p35_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P35_LIM returns the integration limits for problem 35.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = -9.0;
  b = 100.0;

  return;
}
//****************************************************************************80

double p36_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P36_EXACT returns the exact integral for problem 36.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P36_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;

  alpha = p36_param_get ( );

  exact = 1.0 / pow ( alpha + 1.0, 2 );

  return exact;
}
//****************************************************************************80

double *p36_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P36_FUN evaluates the integrand for problem 36.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P36_PARAM_SET.  It had a default value of -0.9.
//
//    The integrand has an endpoint singularity at X=0.
//
//    Suggested values of ALPHA include -0.9 through 2.6.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    X^alpha * ln ( 1 / X )
//
//  Exact Integral:
//
//    1 / ( alpha + 1 )^2
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 83.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P36_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;

  alpha = p36_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] <= 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = pow ( x[i], alpha ) * log ( 1.0 / x[i] );
    }
  }
  return fx;
}
//****************************************************************************80

void p36_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P36_LIM returns the integration limits for problem 36.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

void p36_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P36_PARAM gets or sets the parameter values for problem 36.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = - 0.9;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P36_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P36_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P36_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p36_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P36_PARAM_GET returns the parameter values for problem 36.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p36_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p36_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P36_PARAM_SET sets the parameter values for problem 36.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p36_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p37_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P37_EXACT returns the exact integral for problem 37.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P37_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;
  static double pi = 3.141592653589793;

  alpha = p37_param_get ( );

  exact = atan ( ( 4.0 - pi ) * pow ( 4.0, alpha - 1.0 ) ) 
        + atan (         pi   * pow ( 4.0, alpha - 1.0 ) );

  return exact;
}
//****************************************************************************80

double *p37_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P37_FUN evaluates the integrand for problem 37.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P37_PARAM_SET.  It had a default value of 5.0.
//
//    The integrand has a peak of height 4^ALPHA at X = PI/4.
//
//    Suggested values of ALPHA include 0 through 20.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    4^(-ALPHA) / ( (X-PI/4)^2 + 16^(-ALPHA) )
//
//  Exact Integral:
//
//    atan ( ( 4 - PI ) * 4^(ALPHA-1) ) + atan ( PI * 4^(ALPHA-1) )
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 83.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P37_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;
  static double pi = 3.141592653589793;

  alpha = p37_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = pow ( 4.0, -alpha ) 
      / ( pow ( x[i] - 0.25 * pi, 2 ) + pow ( 16.0, - alpha ) );
  }
  return fx;
}
//****************************************************************************80

void p37_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P37_LIM returns the integration limits for problem 37.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

void p37_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P37_PARAM gets or sets the parameter values for problem 37.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = 5.0;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P37_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P37_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P37_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p37_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P37_PARAM_GET returns the parameter values for problem 37.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p37_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p37_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P37_PARAM_SET sets the parameter values for problem 37.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p37_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p38_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P38_EXACT returns the exact integral for problem 38.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P38_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;
  static double pi = 3.141592653589793;
  double x;

  alpha = p38_param_get ( );

  x = pow ( 2.0, alpha );

  exact = pi * r8_besj0 ( x );

  return exact;
}
//****************************************************************************80

double *p38_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P38_FUN evaluates the integrand for problem 38.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P38_PARAM_SET.
//
//    The integrand oscillates more strongly as ALPHA is increased.
//
//    The suggested range for ALPHA is 0 to 10.
//
//  Interval:
//
//    0 <= X <= PI
//
//  Integrand:
//
//    cos ( 2^ALPHA * sin ( x ) )
//
//  Exact Integral:
//
//    pi * J0 ( 2^ALPHA )
//
//    where J0 ( x ) is the J Bessel function of order 0.
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 83.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P38_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;

  alpha = p38_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = cos ( pow ( 2.0, alpha ) * sin ( x[i] ) );
  }
  return fx;
}
//****************************************************************************80

void p38_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P38_LIM returns the integration limits for problem 38.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  static double pi = 3.141592653589793;

  a = 0.0;
  b = pi;

  return;
}
//****************************************************************************80

void p38_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P38_PARAM gets or sets the parameter values for problem 38.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = 3.0;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P38_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P38_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P38_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p38_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P38_PARAM_GET returns the parameter values for problem 38.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p38_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p38_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P38_PARAM_SET sets the parameter values for problem 38.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p38_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p39_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P39_EXACT returns the exact integral for problem 39.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P39_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;

  alpha = p39_param_get ( );

  exact = ( pow ( 2.0 / 3.0, alpha + 1.0 ) 
          + pow ( 1.0 / 3.0, alpha + 1.0 ) ) / ( alpha + 1.0 );

  return exact;
}
//****************************************************************************80

double *p39_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P39_FUN evaluates the integrand for problem 39.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P39_PARAM_SET.
//
//    The integrand has a singularity at an internal point ( x = 1/3 )
//    with small binary period.
//
//    The suggested range for ALPHA is -0.8 through 2.1.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    ( abs ( x - 1/3 ) )^alpha
//
//  Exact Integral:
//
//    ( (2/3)^(alpha+1) + (1/3)^(alpha+1) ) / ( alpha + 1 )
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 83.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P39_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;

  alpha = p39_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] - 1.0 / 3.0 == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = pow ( r8_abs ( x[i] - 1.0 / 3.0 ), alpha );
    }
  }
  return fx;
}
//****************************************************************************80

void p39_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P39_LIM returns the integration limits for problem 39.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

void p39_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P39_PARAM gets or sets the parameter values for problem 39.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = - 0.5;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P39_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P39_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P39_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

double p39_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P39_PARAM_GET returns the parameter values for problem 39.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p39_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p39_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P39_PARAM_SET sets the parameter values for problem 39.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p39_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p40_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P40_EXACT returns the exact integral for problem 40.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P40_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;
  static double pi = 3.141592653589793;

  alpha = p40_param_get ( );

  exact = ( pow ( 1.0 - 0.25 * pi, alpha + 1.0 ) 
          + pow (     + 0.25 * pi, alpha + 1.0 ) ) 
          / ( alpha + 1.0 );

  return exact;
}
//****************************************************************************80

double *p40_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P40_FUN evaluates the integrand for problem 40.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P40_PARAM_SET.
//
//    The integrand has a singularity at an internal point ( x = pi/4 ).
//
//    The suggested range for ALPHA is -0.8 through 2.1.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    ( abs ( x - pi/4 ) )^alpha
//
//  Exact Integral:
//
//    ( (1-pi/4)^(alpha+1) + (pi/4)^(alpha+1) ) / ( alpha + 1 )
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 83.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P40_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;
  static double pi = 3.141592653589793;

  alpha = p40_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = pow ( r8_abs ( x[i] - 0.25 * pi ), alpha );
  }
  return fx;
}
//****************************************************************************80

void p40_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P40_LIM returns the integration limits for problem 40.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

void p40_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P40_PARAM gets or sets the parameter values for problem 40.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = - 0.5;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P40_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P40_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P40_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p40_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P40_PARAM_GET returns the parameter values for problem 40.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p40_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p40_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P40_PARAM_SET sets the parameter values for problem 40.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p40_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p41_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P41_EXACT returns the exact integral for problem 41.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P41_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;
  static double pi = 3.141592653589793;

  alpha = p41_param_get ( );

  exact = pi / sqrt ( pow ( 1.0 + pow ( 2.0, - alpha ), 2 ) - 1.0 );

  return exact;
}
//****************************************************************************80

double *p41_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P41_FUN evaluates the integrand for problem 41.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P41_PARAM_SET.
//
//    The integrand has a singularity at both endpoints, whose
//    severity increases with ALPHA.
//
//    The suggested range for ALPHA is 1 through 20.
//
//  Interval:
//
//    -1 <= X <= 1
//
//  Integrand:
//
//    1 / ( sqrt ( 1 - x^2 ) * ( x + 1 + 2^(-alpha) ) )
//
//  Exact Integral:
//
//    pi / sqrt ( ( 1 + 2^(-alpha) ) - 1 )
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 83.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P41_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;

  alpha = p41_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( 1.0 - x[i] * x[i] == 0.0 || x[i] + 1.0 + pow ( 0.5, alpha ) == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = 1.0 / 
        ( sqrt ( 1.0 - x[i] * x[i] ) * ( x[i] + 1.0 + pow ( 0.5, alpha ) ) );
    }
  }
  return fx;
}
//****************************************************************************80

void p41_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P41_LIM returns the integration limits for problem 41.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = -1.0;
  b = 1.0;

  return;
}
//****************************************************************************80

void p41_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P41_PARAM gets or sets the parameter values for problem 41.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = 3.0;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P41_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P41_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P41_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p41_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P41_PARAM_GET returns the parameter values for problem 41.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p41_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p41_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P41_PARAM_SET sets the parameter values for problem 41.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p41_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p42_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P42_EXACT returns the exact integral for problem 42.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P42_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;

  alpha = p42_param_get ( );

  exact = pow ( 2.0, alpha - 2.0 ) * pow ( r8_gamma ( alpha / 2.0 ), 2 )
    / r8_gamma ( alpha );

  return exact;
}
//****************************************************************************80

double *p42_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P42_FUN evaluates the integrand for problem 42.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P42_PARAM_SET.
//
//    The integrand has a singularity at X = 0 if ALPHA < 1.
//
//    The suggested range for ALPHA is 0.1 through 2.
//
//  Interval:
//
//    0 <= X <= pi/2
//
//  Integrand:
//
//    ( sin(x) )^( alpha - 1 )
//
//  Exact Integral:
//
//    2^( alpha - 2 ) * ( Gamma(alpha/2) )^2 / Gamma(alpha)
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 83.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P42_FUN[N], the integrand values.
//
{
  double alpha;
  double base;
  double *fx;
  int i;

  alpha = p42_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    base = sin ( x[i] );

    if ( base == 0.0 )
    {
      if ( 1.0 < alpha )
      {
        fx[i] = 0.0;
      }
      else if ( alpha == 1.0 )
      {
        fx[i] = 1.0;
      }
      else
      {
        fx[i] = 0.0;
      }
    }
    else
    {
      fx[i] = pow ( base, alpha - 1.0 );
    }
  }
  return fx;
}
//****************************************************************************80

void p42_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P42_LIM returns the integration limits for problem 42.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  static double pi = 3.141592653589793;

  a = 0.0;
  b = pi / 2.0;

  return;
}
//****************************************************************************80

void p42_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P42_PARAM gets or sets the parameter values for problem 42.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = 0.3;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P42_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P42_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P42_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p42_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P42_PARAM_GET returns the parameter values for problem 42.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p42_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p42_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P42_PARAM_SET sets the parameter values for problem 42.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p42_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p43_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P43_EXACT returns the exact integral for problem 43.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P43_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;

  alpha = p43_param_get ( );

  exact = r8_gamma ( alpha );

  return exact;
}
//****************************************************************************80

double *p43_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P43_FUN evaluates the integrand for problem 43.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P43_PARAM_SET.
//
//    The suggested parameter range is 0.1 <= ALPHA <= 2.0.
//
//    The integrand has an algebraic endpoint singularity at X = 1
//    times a singular factor.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    ( ln ( 1 / x ) )^( alpha - 1 )
//
//  Exact Integral:
//
//    Gamma(alpha)
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 84.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P43_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;

  alpha = p43_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] <= 0.0 )
    {
      fx[i] = 0.0;
    }
    else if ( x[i] == 0.0 )
    {
      if ( alpha - 1.0 < 0.0 )
      {
        fx[i] = 0.0;
      }
      else if ( alpha - 1.0 == 0.0 )
      {
        fx[i] = 1.0;
      }
      else
      {
        fx[i] = 0.0;
      }
    }
    else if ( x[i] == 1.0 )
    {
      if ( alpha - 1.0 < 0.0 )
      {
        fx[i] = 0.0;
      }
      else if ( alpha - 1.0 == 0.0 )
      {
        fx[i] = 1.0;
      }
      else
      {
        fx[i] = 0.0;
      }
    }
    else
    {
      fx[i] = pow ( log ( 1.0 / x[i] ), alpha - 1.0 );
    }
  }
  return fx;
}
//****************************************************************************80

void p43_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P43_LIM returns the integration limits for problem 43.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

void p43_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P43_PARAM gets or sets the parameter values for problem 43.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = 0.3;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P43_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P43_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P43_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p43_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P43_PARAM_GET returns the parameter values for problem 43.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p43_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p43_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P43_PARAM_SET sets the parameter values for problem 43.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p43_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p44_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P44_EXACT returns the exact integral for problem 44.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P44_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;
  double t;

  alpha = p44_param_get ( );

  t = pow ( 2.0, alpha );

  exact = ( 20.0 * sin ( t ) - t * cos ( t ) 
    + t * exp ( - 20.0 ) ) / ( 400.0 + t * t );

  return exact;
}
//****************************************************************************80

double *p44_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P44_FUN evaluates the integrand for problem 44.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P44_PARAM_SET.
//
//    The suggested parameter range is 0.0 <= ALPHA <= 9.0.
//
//    As ALPHA increases, the integrand becomes more oscillatory.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    exp ( 20 * ( x - 1 ) ) * sin ( 2^alpha * x )
//
//  Exact Integral:
//
//    ( 20 sin ( 2^alpha ) - 2^alpha cos ( 2^alpha )
//    + 2^alpha exp ( -20 ) ) / ( 400 + 4^alpha )
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 84.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P44_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;

  alpha = p44_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = exp ( 20.0 * ( x[i] - 1.0 ) ) * sin ( pow ( 2.0, alpha ) * x[i] );
  }
  return fx;
}
//****************************************************************************80

void p44_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P44_LIM returns the integration limits for problem 44.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

void p44_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P44_PARAM gets or sets the parameter values for problem 44.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = 2.0;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P44_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P44_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P44_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p44_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P44_PARAM_GET returns the parameter values for problem 44.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p44_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p44_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P44_PARAM_SET sets the parameter values for problem 44.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p44_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p45_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P45_EXACT returns the exact integral for problem 45.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P45_EXACT, the value of the integral.
//
{
  double alpha;
  double exact;
  static double pi = 3.141592653589793;
  double t;

  alpha = p45_param_get ( );

  t = pow ( 2.0, alpha - 1.0 );

  exact = pi * cos ( t ) * r8_besj0 ( t );

  return exact;
}
//****************************************************************************80

double *p45_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P45_FUN evaluates the integrand for problem 45.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P45_PARAM_SET.
//
//    The suggested parameter range is 0.0 <= ALPHA <= 8.0.
//
//    The function is singular at 0 and 1.
//
//  Interval:
//
//    0 <= X <= 1
//
//  Integrand:
//
//    cos ( 2^alpha * x ) / sqrt ( x * ( 1 - x ) )
//
//  Exact Integral:
//
//    pi * cos ( 2^(alpha-1) ) * J0 ( 2^(alpha-1) )
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 84.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P45_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;

  alpha = p45_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else if ( x[i] == 1.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = cos ( pow ( 2.0, alpha ) * x[i] ) 
        / sqrt ( x[i] * ( 1.0 - x[i] ) );
    }
  }
  return fx;
}
//****************************************************************************80

void p45_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P45_LIM returns the integration limits for problem 45.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

void p45_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P45_PARAM gets or sets the parameter values for problem 45.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = 2.0;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P45_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P45_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P45_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p45_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P45_PARAM_GET returns the parameter values for problem 45.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p45_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p45_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P45_PARAM_SET sets the parameter values for problem 45.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p45_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p46_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P46_EXACT returns the exact integral for problem 46.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P46_EXACT, the value of the integral.
//
{
  double exact;

  exact = 6.0690909595647754101;

  return exact;
}
//****************************************************************************80

double *p46_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P46_FUN evaluates the integrand for problem 46.
//
//  Discussion:
//
//    The problem has a parameter ALPHA that can be set by calling
//    P46_PARAM_SET.
//
//    The integrand is the radius of an ellipse as a function of angle.
//
//    The integral represents the arc length of the ellipse.
//
//    The suggested parameter range is 0.0 <= ALPHA < 1.0.  ALPHA is
//    the eccentricity of the ellipse.
//
//  Interval:
//
//    0 <= theta <= 2 pi
//
//  Integrand:
//
//    r(theta) = ( 1 - alpha^2 ) / ( 1 - alpha * cos ( theta ) )
//
//  Exact Integral:
//
//    When alpha = sin ( pi / 12 ), then
//
//      6.0690909595647754101
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Crandall,
//    Projects in Scientific Computing,
//    Springer, 2000, page 47.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P46_FUN[N], the integrand values.
//
{
  double alpha;
  double *fx;
  int i;

  alpha = p46_param_get ( );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = ( 1.0 - alpha * alpha ) / ( 1.0 - alpha * cos ( x[i] ) );
  }
  return fx;
}
//****************************************************************************80

void p46_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P46_LIM returns the integration limits for problem 46.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  static double pi = 3.141592653589793;

  a = 0.0;
  b = 2.0 * pi;

  return;
}
//****************************************************************************80

void p46_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P46_PARAM gets or sets the parameter values for problem 46.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.
//    "GET" to get the value.
//    "SET" to set the value.
//
//    Input, string NAME, the name of the parameter.
//    "ALPHA" is the only option.
//
//    Input/output, double &VALUE.
//    If the action is "GET", then VALUE returns the current parameter value.
//    If ACTION is "SET", then the parameter value is set to VALUE.
//
{
  static double alpha = 0.0;
  static double pi = 3.141592653589793;
  static bool set = false;

  if ( !set )
  {
    alpha = sin ( pi / 12.0 );
    set = true;
  }

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      value = alpha;
    }
    else
    {
      cerr << "\n";
      cerr << "P46_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "alpha" ) )
    {
      alpha = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P46_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P46_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double p46_param_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    P46_PARAM_GET returns the parameter values for problem 46.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ALPHA, the current value of the parameter.
//
{
  double alpha;

  p46_param ( "get", "alpha", alpha );

  return alpha;
}
//****************************************************************************80

void p46_param_set ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    P46_PARAM_SET sets the parameter values for problem 46.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 March 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the new value of the parameter.
//
{
  p46_param ( "set", "alpha", alpha );

  return;
}
//****************************************************************************80

double p47_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P47_EXACT returns the exact integral for problem 47.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P47_EXACT, the value of the integral.
//
{
  double exact;

  exact = - 4.0 / 9.0;

  return exact;
}
//****************************************************************************80

double *p47_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P47_FUN evaluates the integrand for problem 47.
//
//  Discussion:
//
//    The function is singular at the left endpoint.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    sqrt ( x ) * ln ( x )
//
//  Exact Integral:
//
//    -4/9 = -0.4444...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 101.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P47_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = sqrt ( x[i] ) * log ( x[i] );
    }
  }
  return fx;
}
//****************************************************************************80

void p47_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P47_LIM returns the integration limits for problem 47.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p48_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P48_EXACT returns the exact integral for problem 48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P48_EXACT, the value of the integral.
//
{
  double exact;

  exact = -4.0;

  return exact;
}
//****************************************************************************80

double *p48_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P48_FUN evaluates the integrand for problem 48.
//
//  Discussion:
//
//    The function is singular at the left endpoint.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    ln ( x ) / sqrt ( x )
//
//  Exact Integral:
//
//    -4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 103.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P48_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = log ( x[i] ) / sqrt ( x[i] );
    }
  }
  return fx;
}
//****************************************************************************80

void p48_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P48_LIM returns the integration limits for problem 48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p49_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P49_EXACT returns the exact integral for problem 49.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P49_EXACT, the value of the integral.
//
{
  double exact;

  exact = 61.0 * log ( 2.0 ) + 77.0 * log ( 7.0 ) / 4.0 - 27.0;

  return exact;
}
//****************************************************************************80

double *p49_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P49_FUN evaluates the integrand for problem 49.
//
//  Discussion:
//
//    The function is singular at two internal points, 1 and sqrt(2).
//
//  Interval:
//
//    0 <= x <= 3
//
//  Integrand:
//
//    x^3 * log ( abs ( ( x^2 - 1 ) * ( x^2 - 2 ) ) )
//
//  Exact Integral:
//
//    61 log ( 2 ) + (77/4) log ( 7 ) - 27 = 52.7408...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 104.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P49_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( ( x[i] * x[i] - 1.0 ) * ( x[i] * x[i] - 2.0 ) == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = pow ( x[i], 3 ) * log ( r8_abs ( ( x[i] * x[i] - 1.0 ) 
        * ( x[i] * x[i] - 2.0 ) ) );
    }
  }
  return fx;
}
//****************************************************************************80

void p49_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P49_LIM returns the integration limits for problem 49.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 3.0;

  return;
}
//****************************************************************************80

double p50_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P50_EXACT returns the exact integral for problem 50.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P50_EXACT, the value of the integral.
//
{
  double exact;
  double pi = 3.141592653589793;
  double t;

  t = 10.0 * pi;

  exact = ( - euler_constant ( ) - log ( t ) + r8_ci ( t ) ) / t;

  return exact;
}
//****************************************************************************80

double *p50_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P50_FUN evaluates the integrand for problem 50.
//
//  Discussion:
//
//    The function has a removable singularity at x = 0.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    log ( x ) * sin ( 10 * pi * x )
//
//  Exact Integral:
//
//    ( - gamma - log ( 10 * pi ) + Ci ( 10 * pi ) ) / 10 * pi = -0.1281316...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 106.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P50_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  double pi = 3.141592653589793;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = log ( x[i] ) * sin ( 10.0 * pi * x[i] );
    }
  }
  return fx;
}
//****************************************************************************80

void p50_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P50_LIM returns the integration limits for problem 50.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p51_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P51_EXACT returns the exact integral for problem 51.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P51_EXACT, the value of the integral.
//
{
  double exact;
  double pi = 3.141592653589793;

  exact = - ( r8_ci ( 1.0 ) * sin ( 1.0 ) + 
    ( 0.5 * pi  - r8_si ( 1.0 ) ) * cos ( 1.0 ) ) / pi;

  return exact;
}
//****************************************************************************80

double *p51_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P51_FUN evaluates the integrand for problem 51.
//
//  Discussion:
//
//    The function has a removable singularity at x = 0.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    ln ( x ) / ( 1 + ( ln(x) )^2 )^2
//
//  Exact Integral:
//
//    - ( ci(1) * sin(1) + ( pi/2 - si(1) ) * cos(1) ) / pi = - 0.1892752...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 108.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P51_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = log ( x[i] ) / pow ( 1.0 + pow ( log ( x[i] ), 2 ), 2 );
    }
  }
  return fx;
}
//****************************************************************************80

void p51_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P51_LIM returns the integration limits for problem 51.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p52_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P52_EXACT returns the exact integral for problem 52.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P52_EXACT, the value of the integral.
//
{
  double exact;

  exact = log ( 125.0 / 631.0 ) / 18.0;

  return exact;
}
//****************************************************************************80

double *p52_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P52_FUN evaluates the integrand for problem 52.
//
//  Discussion:
//
//    The function has a singularity at x = 0.
//
//  Interval:
//
//    -1 <= x <= 5
//
//  Integrand:
//
//    1 / ( x * ( 5 * x^3 + 6 ) )
//
//  Exact Integral:
//
//    ln ( 125 / 631 ) / 18 = -0.08994401...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 109.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P52_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] == 0.0 )
    {
      fx[i] = 0.0;
    }
    else
    {
      fx[i] = 1.0 / ( x[i] * ( 5.0 * pow ( x[i], 3 ) + 6.0 ) );
    }
  }
  return fx;
}
//****************************************************************************80

void p52_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P52_LIM returns the integration limits for problem 52.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = -1.0;
  b = 5.0;

  return;
}
//****************************************************************************80

double p53_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P53_EXACT returns the exact integral for problem 53.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P53_EXACT, the value of the integral.
//
{
  double exact;
  double pi = 3.141592653589793;

  exact = 0.5 * pi - atan ( 1.0 / sqrt ( 2.0 ) ) + log ( 3.0 ) / 2.0;

  return exact;
}
//****************************************************************************80

double *p53_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P53_FUN evaluates the integrand for problem 53.
//
//  Discussion:
//
//    The integrand is singular at x = -1 + sqrt ( 3 ) = 0.732...
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    1 / sqrt ( abs ( x^2 + 2 * x - 2 ) )
//
//  Exact Integral:
//
//    pi / 2 - arctan ( 1 / sqrt ( 2 ) ) + ln ( 3 ) / 2 = 1.504622...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Piessens, Elise de Doncker-Kapenga, 
//    Christian Ueberhuber, David Kahaner,
//    QUADPACK: A Subroutine Package for Automatic Integration,
//    Springer, 1983, page 110.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P53_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0 / sqrt ( r8_abs ( x[i] * x[i] + 2.0 * x[i] - 2.0 ) );
  }
  return fx;
}
//****************************************************************************80

void p53_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P53_LIM returns the integration limits for problem 53.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p54_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P54_EXACT returns the exact integral for problem 54.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P54_EXACT, the value of the integral.
//
{
  double exact;

  exact = 2.0 / sqrt ( 3.0 );

  return exact;
}
//****************************************************************************80

double *p54_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P54_FUN evaluates the integrand for problem 54.
//
//  Discussion:
//
//    The reference claims that this integrand is more closely approximated
//    by the trapezoid rule than by Gauss-Legendre quadrature.
//
//    Points  Trapezoid  Gauss-Legendre
//     4      1.91667    2.53883
//    12      2.1594     2.25809
//
//    However, the stated results hardly give one confindence in
//    the convergence of the trapezoid results, and I am unable to
//    confirm them, because my results for 4 points give good results
//    (about 1.14) for BOTH Trapezoid and Gauss-Legendre//
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    2 / ( 2 + sin ( 10 * PI * X ) )
//
//  Exact Integral:
//
//    2 / sqrt ( 3 ) = 1.1547...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Prem Kythe, Pratap Puri,
//    Computational Methods for Linear Integral Equations,
//    Birkhaeuser, 2002.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P54_FUN[N], the integrand values.
//
{
  double *fx;
  int i;
  double pi = 3.141592653589793;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 2.0 / ( 2.0 + sin ( 10.0 * pi * x[i] ) );
  }
  return fx;
}
//****************************************************************************80

void p54_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P54_LIM returns the integration limits for problem 54.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p55_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P55_EXACT returns the exact integral for problem 55.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P55_EXACT, the value of the integral.
//
{
  double a;
  double b;
  double c;
  double exact;
  double pi = 3.141592653589793;
  double x0;

  p55_lim ( a, b );
  p55_param ( "get", "c", c );
  p55_param ( "get", "x0", x0 );

  exact = sqrt ( pi ) * 
    ( r8_erf ( c * ( b - x0 ) ) - r8_erf ( c * ( a - x0 ) ) ) 
    / ( 2.0 * c );

  return exact;
}
//****************************************************************************80

double *p55_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P55_FUN evaluates the integrand for problem 55.
//
//  Interval:
//
//    a = 0 <= x <= 1 = b
//
//  Integrand:
//
//    exp ( - c^2 * ( x - x0 )^2 )
//
//  Exact Integral:
//
//    sqrt ( pi ) 
//    * ( erf ( c * ( b - x0 ) ) - erf ( c * ( a - x0 ) ) )
//    / ( 2 * c )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P55_FUN[N], the integrand values.
//
{
  double c;
  double *fx;
  int i;
  double x0;

  p55_param ( "get", "c", c );
  p55_param ( "get", "x0", x0 );

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = exp ( - c * c * pow ( x[i] - x0, 2 ) );
  }
  return fx;
}
//****************************************************************************80

void p55_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P55_LIM returns the integration limits for problem 55.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

  return;
}
//****************************************************************************80

void p55_param ( string action, string name, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    P55_PARAM sets or gets real scalar parameters for problem 55.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION,
//    "get" to get a parameter.
//    "set" to set a parameter.
//
//    Input, string NAME, the name of the variable.
//    "c" is the coefficient.
//    "x0" is the base point.
//
//    Input/output, double &VALUE.
//    * If ACTION = "set", then VALUE is an input quantity, and is the
//      new value to be assigned to NAME.
//    * If ACTION = "get", then VALUE is an output quantity, and is the
//      current value of NAME.
//
{
  static double c = 3.0;
  static double x0 = 0.75;

  if ( s_eqi ( action, "get" ) )
  {
    if ( s_eqi ( name, "c" ) )
    {
      value = c;
    }
    else if ( s_eqi ( name, "x0" ) )
    {
      value = x0;
    }
    else
    {
      cerr << "\n";
      cerr << "P55_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else if ( s_eqi ( action, "set" ) )
  {
    if ( s_eqi ( name, "c" ) )
    {
      c = value;
    }
    else if ( s_eqi ( name, "x0" ) )
    {
      x0 = value;
    }
    else
    {
      cerr << "\n";
      cerr << "P55_PARAM - Fatal error!\n";
      cerr << "  Unrecognized name.\n";
      exit ( 1 );
    }
  }
  else
  {
    cerr << "\n";
    cerr << "P55_PARAM - Fatal error!\n";
    cerr << "  Unrecognized action.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

double p56_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P56_EXACT returns the estimated integral for problem 56.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P56_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 1.9922524079504000171;

  return exact;
}
//****************************************************************************80

double *p56_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P56_FUN evaluates the integrand for problem 56.
//
//  Interval:
//
//    -1 <= x <= 1
//
//  Integrand:
//
//    1 / ( x^6 + 0.9 )
//
//  Approximate Integral (20 digits):
//
//    1.9922524079504000171...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software,
//    edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P56_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = 1.0 / ( pow ( x[i], 6 ) + 0.9 );
  }
  return fx;
}
//****************************************************************************80

void p56_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P56_LIM returns the integration limits for problem 56.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = - 1.0;
  b = 1.0;

  return;
}
//****************************************************************************80

double p57_exact ( )

//****************************************************************************80
//
//  Purpose:
//
//    P57_EXACT returns the exact integral for problem 57.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P57_EXACT, the estimated value of the integral.
//
{
  double exact;

  exact = 0.4;

  return exact;
}
//****************************************************************************80

double *p57_fun ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P57_FUN evaluates the integrand for problem 57.
//
//  Interval:
//
//    0 <= x <= 1
//
//  Integrand:
//
//    x^(3/2)
//
//  Antiderivative:
//
//    (2/5) * x^(5/2)
//
//  Exact Integral:
//
//    0.4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner,
//    Comparison of Numerical Quadrature Formulas,
//    in Mathematical Software, edited by John R Rice,
//    Academic Press, 1971.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double P57_FUN[N], the integrand values.
//
{
  double *fx;
  int i;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    fx[i] = sqrt ( pow ( x[i], 3 ) );
  }
  return fx;
}
//****************************************************************************80

void p57_lim ( double &a, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    P57_LIM returns the integration limits for problem 57.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double &A, &B, the limits of integration.
//
{
  a = 0.0;
  b = 1.0;

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
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_add ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ADD adds two R8's.
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
//    Input, double X, Y, the numbers to be added.
//
//    Output, double R8_ADD, the sum of X and Y.
//
{
  double value;

  value = x + y;

  return value;
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

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
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
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  double one;
  double temp;
  double test;
  double value;

  one = ( double ) ( 1 );

  value = one;
  temp = value / 2.0;
  test = r8_add ( one, temp );

  while ( one < test )
  {
    value = temp;
    temp = value / 2.0;
    test = r8_add ( one, temp );
  }
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

  if ( i < 1 )
  {
    cout << "\n";
    cout << "R8_MACH - Fatal error!\n";
    cout << "  The input argument I is out of bounds.\n";
    cout << "  Legal values satisfy 1 <= I <= 5.\n";
    cout << "  I = " << i << "\n";
    value = 0.0;
    exit ( 1 );
  }
  else if ( i == 1 )
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
    cout << "\n";
    cout << "R8_MACH - Fatal error!\n";
    cout << "  The input argument I is out of bounds.\n";
    cout << "  Legal values satisfy 1 <= I <= 5.\n";
    cout << "  I = " << i << "\n";
    value = 0.0;
    exit ( 1 );
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

double r8_sech ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SECH evaluates the hyperbolic secant, while avoiding COSH overflow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_SECH, the value of the function.
//
{
  double log_huge = 80.0;
  double value;

  if ( log_huge < r8_abs ( x ) )
  {
    value = 0.0;
  }
  else
  {
    value = 1.0 / cosh ( x );
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

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
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
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
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
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
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

