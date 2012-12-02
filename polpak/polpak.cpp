# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

# include "polpak.hpp"

//****************************************************************************80

double agm ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    AGM computes the arithmetic-geometric mean of A and B.
//
//  Discussion:
//
//    The AGM is defined for nonnegative A and B.
//
//    The AGM of numbers A and B is defined by setting
//
//      A(0) = A,
//      B(0) = B
//
//      A(N+1) = ( A(N) + B(N) ) / 2
//      B(N+1) = sqrt ( A(N) * B(N) )
//
//    The two sequences both converge to AGM(A,B).
//
//    In Mathematica, the AGM can be evaluated by
//
//      ArithmeticGeometricMean [ a, b ]
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
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input, double A, B, the arguments whose AGM is to be computed.
//
//    Output, double AGM, the arithmetic-geometric mean of A and B.
//
{
  double c;
  double d;
  int it;
  int it_max = 1000;
  double tol;
  double value;

  if ( a < 0.0 )
  {
    exit ( 1 );
  }

  if ( b < 0.0 )
  {
    exit ( 1 );
  }

  if ( a == 0.0 || b == 0.0 )
  {
    value = 0.0;
    return value;
  }

  it = 0;
  tol = 100.0 * r8_epsilon ( );

  for ( ; ; )
  {
    it = it + 1;

    c = ( a + b ) / 2.0;
    d = sqrt ( a * b );

    if ( r8_abs ( c - d ) <= tol * ( c + d ) )
    {
      break;
    }

    if ( it_max < it )
    {
      break;
    }

    a = c;
    b = d;
  }
  value = c;

  return value;
}
//****************************************************************************80

void agm_values ( int *n_data, double *a, double *b, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    AGM_VALUES returns some values of the AGM.
//
//  Discussion:
//
//    The AGM is defined for nonnegative A and B.
//
//    The AGM of numbers A and B is defined by setting
//
//      A(0) = A,
//      B(0) = B
//
//      A(N+1) = ( A(N) + B(N) ) / 2
//      B(N+1) = sqrt ( A(N) * B(N) )
//
//    The two sequences both converge to AGM(A,B).
//
//    In Mathematica, the AGM can be evaluated by
//
//      ArithmeticGeometricMean [ a, b ]
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
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, *B, the argument ofs the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 15

  double a_vec[N_MAX] = {
        22.0, 
        83.0, 
        42.0, 
        26.0, 
         4.0, 
         6.0, 
        40.0,
        80.0,
        90.0,
         9.0,
        53.0,
         1.0,
         1.0,
         1.0,
         1.5 }; 
  double b_vec[N_MAX] = {
        96.0,
        56.0,
         7.0,
        11.0,
        63.0,
        45.0,
        75.0,
         0.0,
        35.0,
         1.0,
        53.0,
         2.0,
         4.0,
         8.0,
         8.0 };    
  double fx_vec[N_MAX] = {
        52.274641198704240049,
        68.836530059858524345,
        20.659301196734009322,
        17.696854873743648823,
        23.867049721753300163,
        20.717015982805991662,
        56.127842255616681863,
         0.000000000000000000,
        59.269565081229636528,
        3.9362355036495554780,
        53.000000000000000000,
        1.4567910310469068692,
        2.2430285802876025701,
        3.6157561775973627487,
        4.0816924080221632670 };
 
  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0;
    *b = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double agud ( double gamma )

//****************************************************************************80
//
//  Purpose:
//
//    AGUD evaluates the inverse Gudermannian function.
//
//  Discussion:
//
//    The Gudermannian function relates the hyperbolic and trigonomentric
//    functions.  For any argument X, there is a corresponding value
//    GAMMA so that
//
//      SINH(X) = TAN(GAMMA).
//
//    This value GAMMA(X) is called the Gudermannian of X.  The inverse
//    Gudermannian function is given as input a value GAMMA and computes
//    the corresponding value X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double GAMMA, the value of the Gudermannian.
//
//    Output, double AGUD, the argument of the Gudermannian.
//
{
  double value;

  value = log ( tan ( 0.25 * r8_pi ( ) + 0.5 * gamma ) );

  return value;
}
//****************************************************************************80

int align_enum ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    ALIGN_ENUM counts the alignments of two sequences of M and N elements.
//
//  Discussion:
//
//    We assume that we have sequences A and B of M and N characters each.
//    An alignment of the two sequences is a rule matching corresponding
//    elements of one sequence to another.  Some elements of either sequence
//    can be matched to a null element.  If A(I1) and A(I2) are matched
//    to B(J1) and B(J2), and I1 < I2, then it must be the case that J1 < J2.
//
//    The 5 alignments of a sequence of 1 to a sequence of 2 are:
//
//          _1_   _2_   __3__   __4__   __5__
//
//      A:  1 -   - 1   - 1 -   - - 1   1 - -
//      B:  1 2   1 2   1 - 2   1 2 -   - 1 2
//
//    The formula is:
//
//      F(0,0) = 1
//      F(1,0) = 1
//      F(0,1) = 1
//      F(M,N) = F(M-1,N) + F(M-1,N-1) + F(M,N-1)
//
//    To compute F(M,N), it is not necessary to keep an M+1 by N+1
//    array in memory.  A vector of length N will do.
//
//    F(N,N) is approximately ( 1 + sqrt(2) )^(2*N+1) / sqrt ( N )
//
//  Example:
//
//    The initial portion of the table is:
//
//
//  M/N   0    1    2    3    4       5       6       7       8       9      10
//
//   0    1    1    1    1    1       1       1       1       1       1       1
//   1    1    3    5    7    9      11      13      15      17      19      21
//   2    1    5   13   25   41      61      85     113     145     181     221
//   3    1    7   25   63  129     231     377     575     833    1159    1561
//   4    1    9   41  129  321     681    1289    2241    3649    5641    8361
//   5    1   11   61  231  681    1683    3653    7183   13073   22363   36365
//   6    1   13   85  377 1289    3653    8989   19825   40081   75517  134245
//   7    1   15  113  575 2241    7183   19825   48639  108545  224143  433905
//   8    1   17  145  833 3649   13073   40081  108545  265729  598417 1256465
//   9    1   19  181 1159 5641   22363   75517  224143  598417 1462563 3317445
//  10    1   21  221 1561 8361   36365  134245  433905 1256465 3317445 8097453
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Michael Waterman,
//    Introduction to Computational Biology,
//    Chapman and Hall, 1995, pages 186-190.
//
//  Parameters:
//
//    Input, int M, N, the number of elements of the two sequences.
//
//    Output, int ALIGN_ENUM, the number of possible alignments of the
//    sequences.
//
{
  int *fi;
  int fim1j;
  int fim1jm1;
  int i;
  int j;
  int value;

  if ( m < 0 )
  {
    return 0;
  }
  else if ( n < 0 )
  {
    return 0;
  }
  else if ( m == 0 )
  {
    return 1;
  }
  else if ( n == 0 )
  {
    return 1;
  }

  fi = new int[n+1];
  for ( i = 0; i <= n; i++ )
  {
    fi[i] = 1;
  }

  for ( i = 1; i <= m; i++ )
  {
    fim1jm1 = 1;

    for ( j = 1; j <= n; j++ )
    {

      fim1j = fi[j];

      fi[j] = fi[j] + fi[j-1] + fim1jm1;

      fim1jm1 = fim1j;

    }
  }

  value = fi[n];
  delete [] fi;

  return value;
}
//****************************************************************************80

double arc_cosine ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    ARC_COSINE computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//
//    In particular, you may get the value NaN returned.
//
//    This routine truncates arguments outside the range, avoiding the problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument.
//
//    Output, double ARC_COSINE, an angle whose cosine is C.
//
{
  double value;

  c = r8_max ( c, -1.0 );
  c = r8_min ( c, +1.0 );

  value = acos ( c );

  return value;
}
//****************************************************************************80

double arc_sine ( double s )

//****************************************************************************80
//
//  Purpose:
//
//    ARC_SINE computes the arc sine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ASIN routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//
//    In particular, you may get the value NaN returned.
//
//    This routine truncates arguments outside the range, avoiding the problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double S, the argument.
//
//    Output, double ARC_SINE, an angle whose sine is S.
//
{
  double value;

  s = r8_max ( s, -1.0 );
  s = r8_min ( s, +1.0 );

  value = asin ( s );

  return value;
}
//****************************************************************************80

double asinh2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    ASINH2 returns the inverse hyperbolic sine of a number.
//
//  Discussion:
//
//    The assertion that:
//
//      Y = ASINH2(X)
//
//    implies that
//
//      X = SINH(Y) = 0.5 * ( EXP(Y) - EXP(-Y) ).
//
//    Since a library function ASINH may be available on some systems,
//    this routine is named ASINH2 to avoid name conflicts.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose inverse hyperbolic sine is desired.
//
//    Output, double ASINH2, the inverse hyperbolic sine of X.
//
{
  double value;

  value = log ( x + sqrt ( x * x + 1.0 ) );

  return value;
}
//****************************************************************************80

double atan4 ( double y, double x )

//****************************************************************************80
//
//  Purpose:
//
//    ATAN4 computes the inverse tangent of the ratio Y / X.
//
//  Discussion:
//
//    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
//    the built in functions ATAN and ATAN2 already do.
//
//    However:
//
//    * ATAN4 always returns a positive angle, between 0 and 2 PI,
//      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
//      and [-PI,+PI] respectively;
//
//    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
//     function by contrast always returns an angle in the first or fourth
//     quadrants.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double Y, X, two quantities which represent the tangent of
//    an angle.  If Y is not zero, then the tangent is (Y/X).
//
//    Output, double ATAN4, an angle between 0 and 2 * PI, whose tangent is
//    (Y/X), and which lies in the appropriate quadrant so that the signs
//    of its cosine and sine match those of X and Y.
//
{
  double abs_x;
  double abs_y;
  double theta;
  double theta_0;
//
//  Special cases:
//
  if ( x == 0.0 )
  {
    if ( 0.0 < y )
    {
      theta = r8_pi ( ) / 2.0;
    }
    else if ( y < 0.0 )
    {
      theta = 3.0 * r8_pi ( ) / 2.0;
    }
    else if ( y == 0.0 )
    {
      theta = 0.0;
    }
  }
  else if ( y == 0.0 )
  {
    if ( 0.0 < x )
    {
      theta = 0.0;
    }
    else if ( x < 0.0 )
    {
      theta = r8_pi ( );
    }
  }
//
//  We assume that ATAN2 is correct when both arguments are positive.
//
  else
  {
    abs_y = fabs ( y );
    abs_x = fabs ( x );

    theta_0 = atan2 ( abs_y, abs_x );

    if ( 0.0 < x && 0.0 < y )
    {
      theta = theta_0;
    }
    else if ( x < 0.0 && 0.0 < y )
    {
      theta = r8_pi ( ) - theta_0;
    }
    else if ( x < 0.0 && y < 0.0 )
    {
      theta = r8_pi ( ) + theta_0;
    }
    else if ( 0.0 < x && y < 0.0 )
    {
      theta = 2.0 * r8_pi ( ) - theta_0;
    }

  }

  return theta;
}
//****************************************************************************80

double atanh2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    ATANH2 returns the inverse hyperbolic tangent of a number.
//
//  Definition:
//
//    Y = ATANH2(X) implies that
//    X = TANH(Y) = ( EXP(Y) - EXP(-Y) ) / ( EXP(Y) + EXP(-Y) )
//
//  Discussion:
//
//    Since a library function ATANH may be available on some systems,
//    this routine is named ATANH2 to avoid name conflicts.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose inverse hyperbolic tangent is desired.
//    The absolute value of X should be less than or equal to 1.
//
//    Output, double ATANH2, the inverse hyperbolic tangent of X.
//
{
  double value;

  if ( x <= -1.0 )
  {
    value = -HUGE_VAL;
  }
  else if ( 1.0 <= x )
  {
    value = HUGE_VAL;
  }
  else
  {
    value = 0.5 * log ( ( 1.0 + x ) / ( 1.0 - x ) );
  }
  return value;
}
//****************************************************************************80

void bell ( int n, int b[] )

//****************************************************************************80
//
//  Purpose:
//
//    BELL returns the Bell numbers from 0 to N.
//
//  Discussion:
//
//    The Bell number B(N) is the number of restricted growth functions
//    on N.
//
//    Note that the Stirling numbers of the second kind, S^m_n, count the
//    number of partitions of N objects into M classes, and so it is
//    true that
//
//      B(N) = S^1_N + S^2_N + ... + S^N_N.
//
//  Definition:
//
//    The Bell number B(N) is defined as the number of partitions (of
//    any size) of a set of N distinguishable objects.
//
//    A partition of a set is a division of the objects of the set into
//    subsets.
//
//  Examples:
//
//    There are 15 partitions of a set of 4 objects:
//
//      (1234), (123)(4), (124)(3), (12)(34), (12)(3)(4),
//      (134)(2), (13)(24), (13)(2)(4), (14)(23), (1)(234),
//      (1)(23)(4), (14)(2)(3), (1)(24)(3), (1)(2)(34), (1)(2)(3)(4)
//
//    and so B(4) = 15.
//
//  First values:
//
//     N         B(N)
//     0           1
//     1           1
//     2           2
//     3           5
//     4          15
//     5          52
//     6         203
//     7         877
//     8        4140
//     9       21147
//    10      115975
//
//  Recursion:
//
//    B(I) = sum ( 1 <= J <= I ) Binomial ( I-1, J-1 ) * B(I-J)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of Bell numbers desired.
//
//    Output, int B[N+1], the Bell numbers from 0 to N.
//
{
  int i;
  int j;

  b[0] = 1;

  for ( i = 1; i <= n; i++ )
  {
    b[i] = 0;
    for ( j = 1; j <= i; j++ )
    {
      b[i] = b[i] + b[i-j] * i4_choose ( i-1, j-1 );
    }
  }

  return;
}
//****************************************************************************80

void bell_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    BELL_VALUES returns some values of the Bell numbers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the order of the Bell number.
//
//    Output, int *C, the value of the Bell number.
//
{
# define N_MAX 11

  int c_vec[N_MAX] = { 1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 };
  int n_vec[N_MAX] = { 0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10};

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double benford ( int ival )

//****************************************************************************80
//
//  Purpose:
//
//    BENFORD returns the Benford probability of one or more significant digits.
//
//  Discussion:
//
//    Benford's law is an empirical formula explaining the observed
//    distribution of initial digits in lists culled from newspapers,
//    tax forms, stock market prices, and so on.  It predicts the observed
//    high frequency of the initial digit 1, for instance.
//
//    Note that the probabilities of digits 1 through 9 are guaranteed
//    to add up to 1, since
//      LOG10 ( 2/1 ) + LOG10 ( 3/2) + LOG10 ( 4/3 ) + ... + LOG10 ( 10/9 )
//      = LOG10 ( 2/1 * 3/2 * 4/3 * ... * 10/9 ) = LOG10 ( 10 ) = 1.
//
//    The formula is:
//
//      Prob ( First significant digits are IVAL ) =
//        LOG10 ( ( IVAL + 1 ) / IVAL ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Frank Benford,
//    The Law of Anomalous Numbers,
//    Proceedings of the American Philosophical Society,
//    Volume 78, pages 551-572, 1938.
//
//    T P Hill,
//    The First Digit Phenomenon,
//    American Scientist,
//    Volume 86, July/August 1998, pages 358 - 363.
//
//    R Raimi,
//    The Peculiar Distribution of First Digits,
//    Scientific American,
//    December 1969, pages 109-119.
//
//  Parameters:
//
//    Input, int IVAL, the string of significant digits to be checked.
//    If IVAL is 1, then we are asking for the Benford probability that
//    a value will have first digit 1.  If IVAL is 123, we are asking for
//    the probability that the first three digits will be 123, and so on.
//    Note that IVAL must not be 0 or negative.
//
//    Output, double BENFORD, the Benford probability that an item taken
//    from a real world distribution will have the initial digits IVAL.
//
{
  double value;

  if ( ival <= 0 )
  {
    cerr << "\n";
    cerr << "BENFORD - Fatal error!\n";
    cerr << "  The input argument must be positive.\n";
    cerr << "  Your value was " << ival << "\n";
    exit ( 1 );
  }

  value = log10 ( ( double ) ( ival + 1 ) / ( double ) ival  );

  return value;
}
//****************************************************************************80

void bernoulli_number ( int n, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_NUMBER computes the value of the Bernoulli numbers B(0) through B(N).
//
//  Discussion:
//
//    The Bernoulli numbers are rational.
//
//    If we define the sum of the M-th powers of the first N integers as:
//
//      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
//
//    and let C(I,J) be the combinatorial coefficient:
//
//      C(I,J) = I! / ( ( I - J )! * J! )
//
//    then the Bernoulli numbers B(J) satisfy:
//
//      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)^(M+1-J)
//
//  First values:
//
//   B0  1                   =         1.00000000000
//   B1 -1/2                 =        -0.50000000000
//   B2  1/6                 =         1.66666666666
//   B3  0                   =         0
//   B4 -1/30                =        -0.03333333333
//   B5  0                   =         0
//   B6  1/42                =         0.02380952380
//   B7  0                   =         0
//   B8 -1/30                =        -0.03333333333
//   B9  0                   =         0
//  B10  5/66                =         0.07575757575
//  B11  0                   =         0
//  B12 -691/2730            =        -0.25311355311
//  B13  0                   =         0
//  B14  7/6                 =         1.16666666666
//  B15  0                   =         0
//  B16 -3617/510            =        -7.09215686274
//  B17  0                   =         0
//  B18  43867/798           =        54.97117794486
//  B19  0                   =         0
//  B20 -174611/330          =      -529.12424242424
//  B21  0                   =         0
//  B22  854,513/138         =      6192.123
//  B23  0                   =         0
//  B24 -236364091/2730      =    -86580.257
//  B25  0                   =         0
//  B26  8553103/6           =   1425517.16666
//  B27  0                   =         0
//  B28 -23749461029/870     = -27298231.0678
//  B29  0                   =         0
//  B30  8615841276005/14322 = 601580873.901
//
//  Recursion:
//
//    With C(N+1,K) denoting the standard binomial coefficient,
//
//    B(0) = 1.0
//    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
//
//  Warning:
//
//    This recursion, which is used in this routine, rapidly results
//    in significant errors.
//
//  Special Values:
//
//    Except for B(1), all Bernoulli numbers of odd index are 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the highest Bernoulli number to compute.
//
//    Output, double B[N+1], B(I) contains the I-th Bernoulli number.
//
{
  double b_sum;
  int *combo;
  int i;
  int j;
  bool next;

  if ( n < 0 )
  {
    return;
  }

  b[0] = 1.0;

  if ( n < 1 )
  {
    return;
  }

  b[1] = -0.5;

  next = false;

  combo = new int[n+2];

  for ( i = 2; i <= n; i++ )
  {
    comb_row ( next, i+1, combo );
    next = true;

    if ( ( i % 2 ) == 1 )
    {
      b[i] = 0.0;
    }
    else
    {
      b_sum = 0.0;
      for ( j = 0; j <= i-1; j++ )
      {
        b_sum = b_sum + b[j] * ( double ) ( combo[j] );
      }

      b[i] = - b_sum / ( double ) ( combo[i] );

    }

  }

  delete [] combo;

  return;
}
//****************************************************************************80

void bernoulli_number2 ( int n, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_NUMBER2 evaluates the Bernoulli numbers.
//
//  Discussion:
//
//    The Bernoulli numbers are rational.
//
//    If we define the sum of the M-th powers of the first N integers as:
//
//      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
//
//    and let C(I,J) be the combinatorial coefficient:
//
//      C(I,J) = I! / ( ( I - J )! * J! )
//
//    then the Bernoulli numbers B(J) satisfy:
//
//      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)^(M+1-J)
//
//    Note that the Bernoulli numbers grow rapidly.  Bernoulli number
//    62 is probably the last that can be computed on the VAX without
//    overflow.
//
//    A different method than that used in BERN is employed.
//
//  First values:
//
//   B0  1                   =         1.00000000000
//   B1 -1/2                 =        -0.50000000000
//   B2  1/6                 =         1.66666666666
//   B3  0                   =         0
//   B4 -1/30                =        -0.03333333333
//   B5  0                   =         0
//   B6  1/42                =         0.02380952380
//   B7  0                   =         0
//   B8 -1/30                =        -0.03333333333
//   B9  0                   =         0
//  B10  5/66                =         0.07575757575
//  B11  0                   =         0
//  B12 -691/2730            =        -0.25311355311
//  B13  0                   =         0
//  B14  7/6                 =         1.16666666666
//  B15  0                   =         0
//  B16 -3617/510            =        -7.09215686274
//  B17  0                   =         0
//  B18  43867/798           =        54.97117794486
//  B19  0                   =         0
//  B20 -174611/330          =      -529.12424242424
//  B21  0                   =         0
//  B22  854,513/138         =      6192.123
//  B23  0                   =         0
//  B24 -236364091/2730      =    -86580.257
//  B25  0                   =         0
//  B26  8553103/6           =   1425517.16666
//  B27  0                   =         0
//  B28 -23749461029/870     = -27298231.0678
//  B29  0                   =         0
//  B30  8615841276005/14322 = 601580873.901
//
//  Recursion:
//
//    With C(N+1,K) denoting the standard binomial coefficient,
//
//    B(0) = 1.0
//    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
//
//  Special Values:
//
//    Except for B(1), all Bernoulli numbers of odd index are 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the highest order Bernoulli number to compute.
//
//    Output, double B[N+1], the requested Bernoulli numbers.
//
{
  double altpi;
  int i;
  int k;
  int kmax = 400;
  double sgn;
  double sum2;
  double t;
  double term;
  double tol = 1.0E-06;

  if ( n < 0 )
  {
    return;
  }

  b[0] = 1.0;

  if ( n < 1 )
  {
    return;
  }

  b[1] = -0.5;

  if ( n < 2 )
  {
    return;
  }

  altpi = log ( 2.0 * r8_pi ( ) );
//
//  Initial estimates for B(I), I = 2 to N
//
  b[2] = log ( 2.0 );

  for ( i = 3; i <= n; i++ )
  {
    if ( ( i % 2 ) == 1 )
    {
      b[i] = 0.0;
    }
    else
    {
      b[i] = log ( ( double ) ( i * ( i - 1 ) ) ) + b[i-2];
    }
  }

  b[2] = 1.0 / 6.0;

  if ( n <= 3 )
  {
    return;
  }

  b[4] = -1.0 / 30.0;

  sgn = -1.0;

  for ( i = 6; i <= n; i = i + 2 )
  {
    sgn = -sgn;
    t = 2.0 * sgn * exp ( b[i] - ( double ) ( i ) * altpi );

    sum2 = 1.0;

    for ( k = 2; k <= kmax; k++ )
    {
      term = pow ( ( double ) k, ( double ) -i );

      sum2 = sum2 + term;

      if ( term <= tol * sum2 )
      {
        break;
      }

    }

    b[i] = t * sum2;

  }

  return;
}
//****************************************************************************80

double bernoulli_number3 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_NUMBER3 computes the value of the Bernoulli number B(N).
//
//  Discussion:
//
//    The Bernoulli numbers are rational.
//
//    If we define the sum of the M-th powers of the first N integers as:
//
//      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
//
//    and let C(I,J) be the combinatorial coefficient:
//
//      C(I,J) = I! / ( ( I - J )! * J! )
//
//    then the Bernoulli numbers B(J) satisfy:
//
//      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) 
//        C(M+1,J) B(J) * (N+1)^(M+1-J)
//
//  First values:
//
//     B0  1                   =         1.00000000000
//     B1 -1/2                 =        -0.50000000000
//     B2  1/6                 =         1.66666666666
//     B3  0                   =         0
//     B4 -1/30                =        -0.03333333333
//     B5  0                   =         0
//     B6  1/42                =         0.02380952380
//     B7  0                   =         0
//     B8 -1/30                =        -0.03333333333
//     B9  0                   =         0
//    B10  5/66                =         0.07575757575
//    B11  0                   =         0
//    B12 -691/2730            =        -0.25311355311
//    B13  0                   =         0
//    B14  7/6                 =         1.16666666666
//    B15  0                   =         0
//    B16 -3617/510            =        -7.09215686274
//    B17  0                   =         0
//    B18  43867/798           =        54.97117794486
//    B19  0                   =         0
//    B20 -174611/330          =      -529.12424242424
//    B21  0                   =         0
//    B22  854513/138          =      6192.123
//    B23  0                   =         0
//    B24 -236364091/2730      =    -86580.257
//    B25  0                   =         0
//    B26  8553103/6           =   1425517.16666
//    B27  0                   =         0
//    B28 -23749461029/870     = -27298231.0678
//    B29  0                   =         0
//    B30  8615841276005/14322 = 601580873.901
//
//  Recursion:
//
//    With C(N+1,K) denoting the standard binomial coefficient,
//
//    B(0) = 1.0
//    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
//
//  Special Values:
//
//    Except for B(1), all Bernoulli numbers of odd index are 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the Bernoulli number to compute.
//
//    Output, double BERNOULLI_NUMBER3, the desired Bernoulli number.
//
{
  int i;
  int itmax = 1000;
  double sum2;
  double term;
  double tol = 5.0E-07;
  double value;

  if ( n < 0 )
  {
    value = 0.0;
  }
  else if ( n == 0 )
  {
    value = 1.0;
  }
  else if ( n == 1 )
  {
    value = -0.5;
  }
  else if ( n == 2 )
  {
    value = 1.0 / 6.0;
  }
  else if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    sum2 = 0.0;

    for ( i = 1; i <= itmax; i++ )
    {
      term = 1.0 / pow ( ( double ) i, n );
      sum2 = sum2 + term;

      if ( fabs ( term ) < tol || fabs ( term ) < tol * fabs ( sum2 ) )
      {
        break;
      }

    }

    value = 2.0 * sum2 * r8_factorial ( n ) 
      / pow ( ( 2.0 * r8_pi ( ) ), n );

    if ( ( n % 4 ) == 0 )
    {
      value = -value;
    }

  }

  return value;
}
//****************************************************************************80

void bernoulli_number_values ( int *n_data, int *n, double *c )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_NUMBER_VALUES returns some values of the Bernoulli numbers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the order of the Bernoulli number.
//
//    Output, double *C, the value of the Bernoulli number.
//
{
# define N_MAX 10

  double c_vec[N_MAX] = {
    1.0000000000, -0.5000000000,  0.1666666667, 
    0.0000000000, -0.0333333333, -0.02380952381, 
   -0.0333333333,  0.0757575757, -529.1242424, 
    601580873.9 };

  int n_vec[N_MAX] = {
     0,  1,  2,  3,  4, 6,  8, 10, 20, 30 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double bernoulli_poly ( int n, double x )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_POLY evaluates the Bernoulli polynomial of order N at X.
//
//  Discussion:
//
//    Thanks to Bart Vandewoestyne for pointing out an error in the previous
//    documentation, 31 January 2008.
//
//    Special values of the Bernoulli polynomial include:
//
//      B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
//
//      B'(N,X) = N * B(N-1,X)
//
//      B(N,X+1) - B(N,X) = N * X^(N-1)
//      B(N,X) = (-1)^N * B(N,1-X)
//
//    A formula for the Bernoulli polynomial in terms of the Bernoulli
//    numbers is:
//
//      B(N,X) = sum ( 0 <= K <= N ) B(K) * C(N,K) * X^(N-K)
//
//    The first few polynomials include:
//
//      B(0,X) = 1
//      B(1,X) = X    - 1/2
//      B(2,X) = X^2 -   X      +  1/6
//      B(3,X) = X^3 - 3/2*X^2 +  1/2*X
//      B(4,X) = X^4 - 2*X^3   +      X^2 - 1/30
//      B(5,X) = X^5 - 5/2*X^4 +  5/3*X^3 - 1/6*X
//      B(6,X) = X^6 - 3*X^5   +  5/2*X^4 - 1/2*X^2 + 1/42
//      B(7,X) = X^7 - 7/2*X^6 +  7/2*X^5 - 7/6*X^3 + 1/6*X
//      B(8,X) = X^8 - 4*X^7   + 14/3*X^6 - 7/3*X^4 + 2/3*X^2 - 1/30
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the Bernoulli polynomial to
//    be evaluated.  N must be 0 or greater.
//
//    Input, double X, the value of X at which the polynomial is to
//    be evaluated.
//
//    Output, double BERNOULLI_POLY, the value of B(N,X).
//
{
  int i;
  int *iwork;
  bool next;
  double value;
  double *work;

  work = new double[n+1];
  bernoulli_number ( n, work );

  iwork = new int[n+1];
  next = false;
  comb_row ( next, n, iwork );

  value = 1.0;
  for ( i = 1; i <= n; i++ )
  {
    value = value * x + work[i] * ( double ) iwork[i];
  }

  delete [] iwork;
  delete [] work;

  return value;
}
//****************************************************************************80

double bernoulli_poly2 ( int n, double x )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_POLY2 evaluates the N-th Bernoulli polynomial at X.
//
//  Discussion:
//
//    Thanks to Bart Vandewoestyne for pointing out an error in the previous
//    documentation, 31 January 2008.
//
//    Special values of the Bernoulli polynomial include:
//
//      B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
//
//      B'(N,X) = N * B(N-1,X)
//
//      B(N,X+1) - B(N,X) = N * X^(N-1)
//      B(N,X) = (-1)^N * B(N,1-X)
//
//    A formula for the Bernoulli polynomial in terms of the Bernoulli
//    numbers is:
//
//      B(N,X) = sum ( 0 <= K <= N ) B(K) * C(N,K) * X^(N-K)
//
//    The first few polynomials include:
//
//      B(0,X) = 1
//      B(1,X) = X    - 1/2
//      B(2,X) = X^2 -   X      +  1/6
//      B(3,X) = X^3 - 3/2*X^2 +  1/2*X
//      B(4,X) = X^4 - 2*X^3   +      X^2 - 1/30
//      B(5,X) = X^5 - 5/2*X^4 +  5/3*X^3 - 1/6*X
//      B(6,X) = X^6 - 3*X^5   +  5/2*X^4 - 1/2*X^2 + 1/42
//      B(7,X) = X^7 - 7/2*X^6 +  7/2*X^5 - 7/6*X^3 + 1/6*X
//      B(8,X) = X^8 - 4*X^7   + 14/3*X^6 - 7/3*X^4 + 2/3*X^2 - 1/30
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the Bernoulli polynomial to
//    be evaluated.  N must be 0 or greater.
//
//    Input, double X, the value at which the polynomial is to
//    be evaluated.
//
//    Output, double BERNOULLI_POLY2, the value of B(N,X).
//
{
  double fact;
  int i;
  double value;

  fact = 1.0;

  value = bernoulli_number3 ( 0 );

  for ( i = 1; i <= n; i++ )
  {
    fact = fact * ( double ) ( n + 1 - i ) / ( double ) i;
    value = value * x + fact * bernoulli_number3 ( i );
  }

  return value;
}
//****************************************************************************80

void bernstein_poly ( int n, double x, double bern[] )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY evaluates the Bernstein polynomials at a point X.
//
//  Discussion:
//
//    The Bernstein polynomials are assumed to be based on [0,1].
//
//    The formula is:
//
//      B(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
//
//  First values:
//
//    B(0,0,X) = 1
//
//    B(1,0,X) =      1-X
//    B(1,1,X) =               X
//
//    B(2,0,X) =     (1-X)^2
//    B(2,1,X) = 2 * (1-X)   * X
//    B(2,2,X) =               X^2
//
//    B(3,0,X) =     (1-X)^3
//    B(3,1,X) = 3 * (1-X)^2 * X
//    B(3,2,X) = 3 * (1-X)   * X^2
//    B(3,3,X) =               X^3
//
//    B(4,0,X) =     (1-X)^4
//    B(4,1,X) = 4 * (1-X)^3 * X
//    B(4,2,X) = 6 * (1-X)^2 * X^2
//    B(4,3,X) = 4 * (1-X)   * X^3
//    B(4,4,X) =               X^4
//
//  Special values:
//
//    B(N,I,X) has a unique maximum value at X = I/N.
//
//    B(N,I,X) has an I-fold zero at 0 and and N-I fold zero at 1.
//
//    B(N,I,1/2) = C(N,K) / 2^N
//
//    For a fixed X and N, the polynomials add up to 1:
//
//      Sum ( 0 <= I <= N ) B(N,I,X) = 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein polynomials to be
//    used.  For any N, there is a set of N+1 Bernstein polynomials,
//    each of degree N, which form a basis for polynomials on [0,1].
//
//    Input, double X, the point at which the polynomials are to be evaluated.
//
//    Output, double BERN[N+1], the values of the Bernstein polynomials 
//    of orders 0 through N at X.
//
{
  int i;
  int j;

  if ( n == 0 )
  {
    bern[0] = 1.0;
  }
  else if ( 0 < n )
  {
    bern[0] = 1.0 - x;
    bern[1] = x;

    for ( i = 2; i <= n; i++ )
    {
      bern[i] = x * bern[i-1];
      for ( j = i-1; 1 <= j; j-- )
      {
        bern[j] = x * bern[j-1] + ( 1.0 - x ) * bern[j];
      }
      bern[0] = ( 1.0 - x ) * bern[0];
    }

  }

  return;
}
//****************************************************************************80

void bernstein_poly_values ( int *n_data, int *n, int *k, double *x, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_VALUES returns some values of the Bernstein polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the degree of the polynomial.
//
//    Output, int *K, the index of the polynomial.
//
//    Output, double *X, the argument of the polynomial.
//
//    Output, double *B, the value of the polynomial B(N,K,X).
//
{
# define N_MAX 15
 
  double b_vec[N_MAX] = {
    1.0, 
    0.75,       0.25, 
    0.5625,     0.3750,   0.0625, 
    0.421875,   0.421875, 0.140625,  0.015625, 
    0.31640625, 0.421875, 0.2109375, 0.046875, 0.00390625 };
  int k_vec[N_MAX] = {
    0, 
    0, 1, 
    0, 1, 2, 
    0, 1, 2, 3, 
    0, 1, 2, 3, 4 };
  int n_vec[N_MAX] = {
    0, 
    1, 1, 
    2, 2, 2, 
    3, 3, 3, 3, 
    4, 4, 4, 4, 4 };
  double x_vec[N_MAX] = {
    0.25, 
    0.25, 0.25, 
    0.25, 0.25, 0.25, 
    0.25, 0.25, 0.25, 0.25, 
    0.25, 0.25, 0.25, 0.25, 0.25 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *k = 0;
    *x = 0.0;
    *b = 0.0;
  }
  else
  {
    *n = n_vec[*n_data];
    *k = k_vec[*n_data];
    *x = x_vec[*n_data];
    *b = b_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double beta ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    BETA returns the value of the Beta function.
//
//  Discussion:
//
//    The formula is:
//
//      BETA(X,Y) = ( GAMMA(X) * GAMMA(Y) ) / GAMMA(X+Y)
//
//    Both X and Y must be greater than 0.
//
//    BETA(X,Y) = BETA(Y,X).
//    BETA(X,Y) = Integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
//    BETA(X,Y) = GAMMA(X) * GAMMA(Y) / GAMMA(X+Y)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the two parameters that define the Beta function.
//    X and Y must be greater than 0.
//
//    Output, double BETA, the value of the Beta function.
//
{
  double value;

  if ( x <= 0.0 || y <= 0.0 )
  {
    cerr << "\n";
    cerr << "BETA - Fatal error!\n";
    cerr << "  Both X and Y must be greater than 0.\n";
    exit ( 1 );
  }

  value = exp ( lgamma ( x ) 
              + lgamma ( y ) 
              - lgamma ( x + y ) );

  return value;
}
//****************************************************************************80

void beta_values ( int *n_data, double *x, double *y, double *fxy )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_VALUES returns some values of the Beta function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, double *X, *Y, the arguments of the function.
//
//    Output, double *FXY, the value of the function.
//
{
# define N_MAX 17

  double b_vec[N_MAX] = {
    5.000000, 2.500000, 1.666667, 1.250000, 
    5.000000, 2.500000, 1.000000, 1.666667E-01, 
    0.333333E-01, 7.142857E-03, 1.587302E-03, 0.238095E-01, 
    5.952381E-03, 1.984127E-03, 7.936508E-04, 3.607504E-04, 
    8.325008E-05 };

  double x_vec[N_MAX] = {
    0.2, 0.4, 0.6, 0.8, 
    1.0, 1.0, 1.0, 2.0, 
    3.0, 4.0, 5.0, 6.0, 
    6.0, 6.0, 6.0, 6.0, 
    7.0 };

  double y_vec[N_MAX] = {
    1.0, 1.0, 1.0, 1.0, 
    0.2, 0.4, 1.0, 2.0, 
    3.0, 4.0, 5.0, 2.0, 
    3.0, 4.0, 5.0, 6.0, 
    7.0 };

  if ( n_data < 0 )
  {
    *n_data = 0;
  }
  
  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *y = 0.0;
    *fxy = 0.0;
  }
  else
  {
    *x = x_vec[*n_data];
    *y = y_vec[*n_data];
    *fxy = b_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bpab ( int n, double x, double a, double b, double bern[] )

//****************************************************************************80
//
//  Purpose:
//
//    BPAB evaluates at X the Bernstein polynomials based in [A,B].
//
//  Discussion:
//
//    The formula is:
//
//      BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
//
//  First values:
//
//    B(0,0,X) =   1
//
//    B(1,0,X) = (      B-X              ) / (B-A)
//    B(1,1,X) = (                X-A    ) / (B-A)
//
//    B(2,0,X) = (     (B-X)^2           ) / (B-A)^2
//    B(2,1,X) = ( 2 * (B-X)   * (X-A)   ) / (B-A)^2
//    B(2,2,X) = (               (X-A)^2 ) / (B-A)^2
//
//    B(3,0,X) = (     (B-X)^3           ) / (B-A)^3
//    B(3,1,X) = ( 3 * (B-X)^2 * (X-A)   ) / (B-A)^3
//    B(3,2,X) = ( 3 * (B-X)   * (X-A)^2 ) / (B-A)^3
//    B(3,3,X) = (               (X-A)^3 ) / (B-A)^3
//
//    B(4,0,X) = (     (B-X)^4           ) / (B-A)^4
//    B(4,1,X) = ( 4 * (B-X)^3 * (X-A)   ) / (B-A)^4
//    B(4,2,X) = ( 6 * (B-X)^2 * (X-A)^2 ) / (B-A)^4
//    B(4,3,X) = ( 4 * (B-X)   * (X-A)^3 ) / (B-A)^4
//    B(4,4,X) = (               (X-A)^4 ) / (B-A)^4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein polynomials to be used.
//    For any N, there is a set of N+1 Bernstein polynomials, each of
//    degree N, which form a basis for polynomials on [A,B].
//
//    Input, double X, the point at which the polynomials are to be evaluated.
//
//    Input, double A, B, the endpoints of the interval on which the
//    polynomials are to be based.  A and B should not be equal.
//
//    Output, double BERN[N+1], the values of the N+1 Bernstein polynomials at X.
//
{
  int i;
  int j;

  if ( b == a )
  {
    cerr << "\n";
    cerr << "BPAB - Fatal error!\n";
    cerr << "  A = B = " << a << "\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    bern[0] = 1.0;
  }
  else if ( 0 < n );
  {
    bern[0] = ( b - x ) / ( b - a );
    bern[1] = ( x - a ) / ( b - a );

    for ( i = 2; i <= n; i++ )
    {
      bern[i] = ( x - a ) * bern[i-1] / ( b - a );
      for ( j = i-1; 1 <= j; j-- )
      {
        bern[j] = ( ( b - x     ) * bern[j] 
                  + (     x - a ) * bern[j-1] ) 
                  / ( b     - a );
      }
      bern[0] = ( b - x ) * bern[0] / ( b - a );
    }

  }

  return;
}
//****************************************************************************80

double *cardan ( int n, double x, double s )

//****************************************************************************80
//
//  Purpose:
//
//    CARDAN evaluates the Cardan polynomials.
//
//  First terms:
//
//     N  C(N,S,X)
//
//     0  2
//     1  X
//     2  X^2  -  2 S
//     3  X^3  -  3 S X
//     4  X^4  -  4 S X^2 +  2 S^2
//     5  X^5  -  5 S X^3 +  5 S^2 X
//     6  X^6  -  6 S X^4 +  9 S^2 X^2 -  2 S^3
//     7  X^7  -  7 S X^5 + 14 S^2 X^3 -  7 S^3 X
//     8  X^8  -  8 S X^6 + 20 S^2 X^4 - 16 S^3 X^2 +  2 S^4
//     9  X^9  -  9 S X^7 + 27 S^2 X^5 - 30 S^3 X^3 +  9 S^4 X
//    10  X^10 - 10 S X^8 + 35 S^2 X^6 - 50 S^3 X^4 + 25 S^4 X^2 -  2 S^5
//    11  X^11 - 11 S X^9 + 44 S^2 X^7 - 77 S^3 X^5 + 55 S^4 X^3 - 11 S^5 X
//
//  Recursion:
//
//    Writing the N-th polynomial in terms of its coefficients:
//
//      C(N,S,X) = sum ( 0 <= I <= N ) D(N,I) * S^(N-I)/2 * X^I
//
//    then
//
//    D(0,0) = 1
//
//    D(1,1) = 1
//    D(1,0) = 0
//
//    D(N,N) = 1
//    D(N,K) = D(N-1,K-1) - D(N-2,K)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Thomas Osler,
//    Cardan Polynomials and the Reduction of Radicals,
//    Mathematics Magazine,
//    Volume 74, Number 1, February 2001, pages 26-32.
//
//  Parameters:
//
//    Input, int N, the highest polynomial to compute.
//
//    Input, double X, the point at which the polynomials are to be computed.
//
//    Input, double S, the value of the parameter, which must be positive.
//
//    Output, double CARDAN[N+1], the values of the Cardan polynomials at X.
//
{
  double fact;
  int i;
  double s2;
  double *v;
  double x2[1];

  s2 = sqrt ( s );
  x2[0] = 0.5 * x / s2;

  v = cheby_t_poly ( 1, n, x2 );

  fact = 1.0;

  for ( i = 0; i <= n; i++ )
  {
    v[i] = 2.0 * fact * v[i];
    fact = fact * s2;
  }

  return v;
}
//****************************************************************************80

void cardan_poly_coef ( int n, double s, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    CARDAN_POLY_COEF computes the coefficients of the N-th Cardan polynomial.
//
//  First terms:
//
//    2
//    0       1
//   -2 S     0       1
//    0      -3 S     0       1
//    2 S^2  0      -4 S     0       1
//    0       5 S^2  0      -5 S     0       1
//   -2 S**3  0       9 S^2  0      -6 S     0       1
//    0       7 S**3  0      14 S^2  0      -7 S     0       1
//    2 S**4  0     -16 S**3  0      20 S^2  0      -8 S     0        1
//    0       9 S**4  0     -30 S**3  0      27 S^2  0      -9 S      0     1
//   -2 S**5  0      25 S**4  0     -50 S**3  0      35 S^2  0      -10 S   0   1
//    0     -11 S**5  0      55 S**4  0     -77 S**3  0     +44 S^2   0   -11 S 0 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Thomas Osler,
//    Cardan Polynomials and the Reduction of Radicals,
//    Mathematics Magazine,
//    Volume 74, Number 1, February 2001, pages 26-32.
//
//  Parameters:
//
//    Input, int N, the order of the polynomial
//
//    Input, double S, the value of the parameter, which must be positive.
//
//    Output, double C[N+1], the coefficients.  C(0) is the constant term,
//    and C(N) is the coefficient of X^N.
//
{
  double *cm1;
  double *cm2;
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  c[0] = 2.0;
  for ( i = 1; i <= n; i++ )
  {
    c[i] = 0.0;
  }

  if ( n == 0 )
  {
    return;
  }

  cm1 = new double[n+1];
  cm2 = new double[n+1];

  for ( i = 0; i <= n; i++ )
  {
    cm1[i] = c[i];
  }

  c[0] = 0.0;
  c[1] = 1.0;
  for ( i = 2; i <= n; i++ )
  {
    c[i] = 0.0;
  }

  for ( i = 2; i <= n; i++ )
  {

    for ( j = 0; j <= i-2; j++ )
    {
      cm2[j] = cm1[j];
    }

    for ( j = 0; j <= i-1; j++ )
    {
      cm1[j] = c[j];
    }

    c[0] = 0.0;
    for ( j = 1; j <= i; j++ )
    {
      c[j] = cm1[j-1];
    }

    for ( j = 0; j <= i-2; j++ )
    {
      c[j] = c[j] - s * cm2[j];
    }
  }

  delete [] cm1;
  delete [] cm2;

  return;
}
//****************************************************************************80

void catalan ( int n, int c[] )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN computes the Catalan numbers, from C(0) to C(N).
//
//  Discussion:
//
//    The Catalan number C(N) counts:
//
//    1) the number of binary trees on N vertices;
//    2) the number of ordered trees on N+1 vertices;
//    3) the number of full binary trees on 2N+1 vertices;
//    4) the number of well formed sequences of 2N parentheses;
//    5) the number of ways 2N ballots can be counted, in order,
//       with N positive and N negative, so that the running sum
//       is never negative;
//    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
//    7) the number of monotone functions from [1..N} to [1..N} which
//       satisfy f(i) <= i for all i;
//    8) the number of ways to triangulate a polygon with N+2 vertices.
//
//    The formula is:
//
//      C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
//           = 1 / (N+1) * COMB ( 2N, N )
//           = 1 / (2N+1) * COMB ( 2N+1, N+1).
//
//  First values:
//
//     C(0)     1
//     C(1)     1
//     C(2)     2
//     C(3)     5
//     C(4)    14
//     C(5)    42
//     C(6)   132
//     C(7)   429
//     C(8)  1430
//     C(9)  4862
//    C(10) 16796
//
//  Recursion:
//
//    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
//    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
//
//  Example:
//
//    N = 3
//
//    ()()()
//    ()(())
//    (()())
//    (())()
//    ((()))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dennis Stanton and Dennis White,
//    Constructive Combinatorics,
//    Springer Verlag, New York, 1986.
//
//  Parameters:
//
//    Input, int N, the number of Catalan numbers desired.
//
//    Output, int C[N+1], the Catalan numbers from C(0) to C(N).
//
{
  int i;

  if ( n < 0 )
  {
    return;
  }

  c[0] = 1;
//
//  The extra parentheses ensure that the integer division is
//  done AFTER the integer multiplication.
//
  for ( i = 1; i <= n; i++ )
  {
    c[i] = ( c[i-1] * 2 * ( 2 * i - 1 ) ) / ( i + 1 );
  }

  return;
}
//****************************************************************************80

void catalan_row_next ( bool next, int n, int irow[] )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN_ROW_NEXT computes row N of Catalan's triangle.
//
//  Example:
//
//    I\J 0   1   2   3   4   5   6
//
//    0   1
//    1   1   1
//    2   1   2   2
//    3   1   3   5   5
//    4   1   4   9  14  14
//    5   1   5  14  28  42  42
//    6   1   6  20  48  90 132 132
//
//  Recursion:
//
//    C(0,0) = 1
//    C(I,0) = 1
//    C(I,J) = 0 for I < J
//    C(I,J) = C(I,J-1) + C(I-1,J)
//    C(I,I) is the I-th Catalan number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, bool NEXT, indicates whether this is a call for
//    the 'next' row of the triangle.
//    NEXT = FALSE, this is a startup call.  Row N is desired, but
//    presumably this is a first call, or row N-1 was not computed
//    on the previous call.
//    NEXT = TRUE, this is not the first call, and row N-1 was computed
//    on the previous call.  In this case, much work can be saved
//    by using the information from the previous values of IROW
//    to build the next values.
//
//    Input, int N, the index of the row of the triangle desired.
//
//    Input/output, int IROW(0:N), the row of coefficients.
//    If NEXT = FALSE, then IROW is not required to be set on input.
//    If NEXT = TRUE, then IROW must be set on input to the value of
//    row N-1.
//
{
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  if ( !next )
  {
    irow[0] = 1;
    for ( i = 1; i <= n; i++ )
    {
      irow[i] = 0;
    }

    for ( i = 1; i <= n; i++ )
    {
      irow[0] = 1;

      for ( j = 1; j <= i-1; j++ )
      {
        irow[j] = irow[j] + irow[j-1];
      }

      irow[i] = irow[i-1];

    }
  }
  else
  {
    irow[0] = 1;

    for ( j = 1; j <= n-1; j++ )
    {
      irow[j] = irow[j] + irow[j-1];
    }

    if ( 1 <= n )
   {
      irow[n] = irow[n-1];
    }

  }

  return;
}
//****************************************************************************80

void catalan_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN_VALUES returns some values of the Catalan numbers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the order of the Catalan number.
//
//    Output, int *C, the value of the Catalan number.
//
{
# define N_MAX 11

  int c_vec[N_MAX] = { 1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 };
  int n_vec[N_MAX] = { 0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void charlier ( int n, double a, double x, double value[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHARLIER evaluates Charlier polynomials at a point.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J Simoes Pereira,
//    Algorithm 234: Poisson-Charliers Polynomials,
//    Communications of the ACM,
//    Volume 7, Number 7, page 420, July 1964.
//
//    Walter Gautschi,
//    Orthogonal Polynomials: Computation and Approximation,
//    Oxford, 2004,
//    ISBN: 0-19-850672-4,
//    LC: QA404.5 G3555.
//
//    Gabor Szego,
//    Orthogonal Polynomials,
//    American Mathematical Society, 1975,
//    ISBN: 0821810235,
//    LC: QA3.A5.v23.
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input, int N, the maximum order of the polynomial.  
//    N must be at least 0.
//
//    Input, double A, the parameter.  A must not be 0.
//
//    Input, double X, the evaluation point.
//
//    Output, double VALUE[0:N], the value of the polynomials at X.
//
{
  int i;

  if ( a == 0.0 )
  {
    cerr << "\n";
    cerr << "CHARLIER - Fatal error!\n";
    cerr << "  Parameter A cannot be zero.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "CHARLIER - Fatal error!\n";
    cerr << "  Parameter N must be nonnegative.\n";
    exit ( 1 );
  }

  value[0] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  value[1] = - x / a;

  if ( n == 1 )
  {
    return;
  }

  for ( i = 1; i < n; i++ )
  {
    value[i+1] = ( ( i + a - x ) * value[i] - i * value[i-1] ) / a;
  }

  return;
}
//****************************************************************************80

double *cheby_t_poly ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_POLY evaluates Chebyshev polynomials T(n,x).
//
//  Discussion:
//
//    Chebyshev polynomials are useful as a basis for representing the
//    approximation of functions since they are well conditioned, in the sense
//    that in the interval [-1,1] they each have maximum absolute value 1.
//    Hence an error in the value of a coefficient of the approximation, of
//    size epsilon, is exactly reflected in an error of size epsilon between
//    the computed approximation and the theoretical approximation.
//
//    Typical usage is as follows, where we assume for the moment
//    that the interval of approximation is [-1,1].  The value
//    of N is chosen, the highest polynomial to be used in the
//    approximation.  Then the function to be approximated is
//    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
//    Chebyshev polynomial.  Let these values be denoted by F(XJ).
//
//    The coefficients of the approximation are now defined by
//
//      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I,XJ)
//
//    except that C(0) is given a value which is half that assigned
//    to it by the above formula,
//
//    and the representation is
//
//    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
//
//    Now note that, again because of the fact that the Chebyshev polynomials
//    have maximum absolute value 1, if the higher order terms of the
//    coefficients C are small, then we have the option of truncating
//    the approximation by dropping these terms, and we will have an
//    exact value for maximum perturbation to the approximation that
//    this will cause.
//
//    It should be noted that typically the error in approximation
//    is dominated by the first neglected basis function (some multiple of
//    T(N+1,X) in the example above).  If this term were the exact error,
//    then we would have found the minimax polynomial, the approximating
//    polynomial of smallest maximum deviation from the original function.
//    The minimax polynomial is hard to compute, and another important
//    feature of the Chebyshev approximation is that it tends to behave
//    like the minimax polynomial while being easy to compute.
//
//    To evaluate a sum like
//
//      sum ( 0 <= J <= N ) C(J) T(J,X),
//
//    Clenshaw's recurrence formula is recommended instead of computing the
//    polynomial values, forming the products and summing.
//
//    Assuming that the coefficients C(J) have been computed
//    for J = 0 to N, then the coefficients of the representation of the
//    indefinite integral of the function may be computed by
//
//      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1,
//
//    with
//
//      C(N+1)=0
//      B(0) arbitrary.
//
//    Also, the coefficients of the representation of the derivative of the
//    function may be computed by:
//
//      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0,
//
//    with
//
//      D(N+1) = D(N)=0.
//
//    Some of the above may have to adjusted because of the irregularity of C(0).
//
//    The formula is:
//
//      T(N,X) = COS(N*ARCCOS(X))
//
//  Differential equation:
//
//    (1-X*X) Y'' - X Y' + N N Y = 0
//
//  First terms:
//
//    T(0,X) =  1
//    T(1,X) =  1 X
//    T(2,X) =  2 X^2 -   1
//    T(3,X) =  4 X^3 -   3 X
//    T(4,X) =  8 X^4 -   8 X^2 +  1
//    T(5,X) = 16 X^5 -  20 X^3 +  5 X
//    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
//    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
//
//  Inequality:
//
//    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
//
//  Orthogonality:
//
//    For integration over [-1,1] with weight
//
//      W(X) = 1 / sqrt(1-X*X),
//
//    if we write the inner product of T(I,X) and T(J,X) as
//
//      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
//
//    then the result is:
//
//      0 if I /= J
//      PI/2 if I == J /= 0
//      PI if I == J == 0
//
//    A discrete orthogonality relation is also satisfied at each of
//    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
//                              = 0 if I /= J
//                              = N/2 if I == J /= 0
//                              = N if I == J == 0
//
//  Recursion:
//
//    T(0,X) = 1,
//    T(1,X) = X,
//    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
//
//    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
//
//  Special values:
//
//    T(N,1) = 1
//    T(N,-1) = (-1)^N
//    T(2N,0) = (-1)^N
//    T(2N+1,0) = 0
//    T(N,X) = (-1)^N * T(N,-X)
//
//  Zeroes:
//
//    M-th zero of T(N,X) is cos((2*M-1)*PI/(2*N)), M = 1 to N
//
//  Extrema:
//
//    M-th extremum of T(N,X) is cos(PI*M/N), M = 0 to N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest polynomial to compute.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double CHEBY_T_POLY[M*(N+1)], the values of the Chebyshev polynomials.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i] = 1.0;
  }
  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = 2.0 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
    }
  }

  return v;
}
//****************************************************************************80

double *cheby_t_poly_coef ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_POLY_COEF evaluates coefficients of Chebyshev polynomials T(n,x).
//
//  First terms:
//
//    N/K     0     1      2      3       4     5      6    7      8    9   10
//
//     0      1
//     1      0     1
//     2     -1     0      2
//     3      0    -3      0      4
//     4      1     0     -8      0       8
//     5      0     5      0    -20       0    16
//     6     -1     0     18      0     -48     0     32
//     7      0    -7      0     56       0  -112      0    64
//
//  Recursion:
//
//    T(0,X) = 1,
//    T(1,X) = X,
//    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
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
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Output, double CHEBY_T_POLY_COEF[(N+1)*(N+1)], the coefficients of 
//    the Chebyshev T polynomials.
//
{
  double *c;
  int i;
  int j;

  if ( n < 0 )
  {
    return NULL;
  }

  c = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n == 0 )
  {
    return c;
  }

  c[1+1*(n+1)] = 1.0;

  for ( i = 2; i <= n; i++ )
  {
    c[i+0*(n+1)]     =                    - c[i-2+0*(n+1)];
    for ( j = 1; j <= i - 2; j++ )
    {
      c[i+j*(n+1)] = 2.0 * c[i-1+(j-1)*(n+1)] - c[i-2+j*(n+1)];
    }
    c[i+(i-1)*(n+1)] = 2.0 * c[i-1+(i-2)*(n+1)];
    c[i+ i   *(n+1)] = 2.0 * c[i-1+(i-1)*(n+1)];
  }

  return c;
}
//****************************************************************************80

void cheby_t_poly_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_POLY_VALUES returns values of Chebyshev polynomials T(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ChebyshevT[n,x]
//
//    Chebyshev polynomials are useful as a basis for representing the
//    approximation of functions since they are well conditioned, in the sense
//    that in the interval [-1,1] they each have maximum absolute value 1.
//    Hence an error in the value of a coefficient of the approximation, of
//    size epsilon, is exactly reflected in an error of size epsilon between
//    the computed approximation and the theoretical approximation.
//
//    Typical usage is as follows, where we assume for the moment
//    that the interval of approximation is [-1,1].  The value
//    of N is chosen, the highest polynomial to be used in the
//    approximation.  Then the function to be approximated is
//    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
//    Chebyshev polynomial.  Let these values be denoted by F(XJ).
//
//    The coefficients of the approximation are now defined by
//
//      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I)(XJ)
//
//    except that C(0) is given a value which is half that assigned
//    to it by the above formula,
//
//    and the representation is
//
//    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J)(X)
//
//    Now note that, again because of the fact that the Chebyshev polynomials
//    have maximum absolute value 1, if the higher order terms of the
//    coefficients C are small, then we have the option of truncating
//    the approximation by dropping these terms, and we will have an
//    exact value for maximum perturbation to the approximation that
//    this will cause.
//
//    It should be noted that typically the error in approximation
//    is dominated by the first neglected basis function (some multiple of
//    T(N+1)(X) in the example above).  If this term were the exact error,
//    then we would have found the minimax polynomial, the approximating
//    polynomial of smallest maximum deviation from the original function.
//    The minimax polynomial is hard to compute, and another important
//    feature of the Chebyshev approximation is that it tends to behave
//    like the minimax polynomial while being easy to compute.
//
//    To evaluate a sum like
//
//      sum ( 0 <= J <= N ) C(J) T(J)(X),
//
//    Clenshaw's recurrence formula is recommended instead of computing the
//    polynomial values, forming the products and summing.
//
//    Assuming that the coefficients C(J) have been computed
//    for J = 0 to N, then the coefficients of the representation of the
//    indefinite integral of the function may be computed by
//
//      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1,
//
//    with
//
//      C(N+1)=0
//      B(0) arbitrary.
//
//    Also, the coefficients of the representation of the derivative of the
//    function may be computed by:
//
//      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0,
//
//    with
//
//      D(N+1) = D(N) = 0.
//
//    Some of the above may have to adjusted because of the irregularity of C(0).
//
//  Differential equation:
//
//    (1-X*X) * Y'' - X * Y' + N * N * Y = 0;
//
//  Formula:
//
//    T(N)(X) = cos ( N * arccos(X) )
//
//  First terms:
//
//    T(0)(X) =  1
//    T(1)(X) =  1 X
//    T(2)(X) =  2 X^2 -   1
//    T(3)(X) =  4 X^3 -   3 X
//    T(4)(X) =  8 X^4 -   8 X^2 +  1
//    T(5)(X) = 16 X^5 -  20 X^3 +  5 X
//    T(6)(X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
//    T(7)(X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
//
//  Inequality:
//
//    abs ( T(N)(X) ) <= 1 for -1 <= X <= 1
//
//  Orthogonality:
//
//    For integration over [-1,1] with weight
//
//      W(X) = 1 / sqrt(1-X*X),
//
//    if we write the inner product of T(I)(X) and T(J)(X) as
//
//      < T(I)(X), T(J)(X) > = integral ( -1 <= X <= 1 ) W(X) T(I)(X) T(J)(X) dX
//
//    then the result is:
//
//      0 if I /= J
//      PI/2 if I == J /= 0;
//      PI if I == J == 0;
//
//    A discrete orthogonality relation is also satisfied at each of
//    the N zeroes of T(N)(X):  sum ( 1 <= K <= N ) T(I)(X) * T(J)(X)
//                              = 0 if I /= J
//                              = N/2 if I == J /= 0;
//                              = N if I == J == 0;
//
//  Recursion:
//
//    T(0)(X) = 1,
//    T(1)(X) = X,
//    T(N)(X) = 2 * X * T(N-1)(X) - T(N-2)(X)
//
//    T'(N)(X) = N * ( -X * T(N)(X) + T(N-1)(X) ) / ( 1 - X^2 )
//
//  Special values:
//
//    T(N)(1) = 1
//    T(N)(-1) = (-1)^N
//    T(2N)(0) = (-1)^N
//    T(2N+1)(0) = 0;
//    T(N)(X) = (-1)^N * T(N)(-X)
//
//  Zeroes:
//
//    M-th zero of T(N)(X) is cos((2*M-1)*PI/(2*N)), M = 1 to N
//
//  Extrema:
//
//    M-th extremum of T(N)(X) is cos(PI*M/N), M = 0 to N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2004
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
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.8000000000000000E+00,
      0.2800000000000000E+00,
     -0.3520000000000000E+00,
     -0.8432000000000000E+00,
     -0.9971200000000000E+00,
     -0.7521920000000000E+00,
     -0.2063872000000000E+00,
      0.4219724800000000E+00,
      0.8815431680000000E+00,
      0.9884965888000000E+00,
      0.7000513740800000E+00,
      0.1315856097280000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *cheby_t_poly_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_POLY_ZERO returns zeroes of Chebyshev polynomials T(n,x).
//
//  Discussion:
//
//    The I-th zero of T(N,X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
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
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double CHEBY_T_POLY_ZERO[N], the zeroes of T(N,X).
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;
  double *z;

  z = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    angle = ( double) ( 2 * i - 1 ) * pi / ( double ) ( 2 * n );
    z[i-1] = cos ( angle );
  }

  return z;
}
//****************************************************************************80

void cheby_u_poly ( int n, double x, double cx[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_POLY evaluates Chebyshev polynomials U(n,x).
//
//  Differential equation:
//
//    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
//
//  First terms:
//
//    U(0,X) =   1
//    U(1,X) =   2 X
//    U(2,X) =   4 X^2 -   1
//    U(3,X) =   8 X^3 -   4 X
//    U(4,X) =  16 X^4 -  12 X^2 +  1
//    U(5,X) =  32 X^5 -  32 X^3 +  6 X
//    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
//    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
//
//  Recursion:
//
//    U(0,X) = 1,
//    U(1,X) = 2 * X,
//    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
//
//  Norm:
//
//    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = PI/2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the highest polynomial to compute.
//
//    Input, double X, the point at which the polynomials are to be computed.
//
//    Output, double CX[N+1], the values of the N+1 Chebyshev polynomials.
//
{
  int i;

  if ( n < 0 )
  {
    return;
  }

  cx[0] = 1.0;

  if ( n < 1 )
  {
    return;
  }

  cx[1] = 2.0 * x;

  for ( i = 2; i <= n; i++ )
  {
    cx[i] = 2.0 * x * cx[i-1] - cx[i-2];
  }

  return;
}
//****************************************************************************80

void cheby_u_poly_coef ( int n, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_POLY_COEF evaluates coefficients of Chebyshev polynomials U(n,x).
//
//  First terms:
//
//    N/K     0     1      2      3       4     5      6    7      8    9   10
//
//     0      1
//     1      0     2
//     2     -1     0      4
//     3      0    -4      0      8
//     4      1     0    -12      0      16
//     5      0     6      0    -32       0    32
//     6     -1     0     24      0     -80     0     64
//     7      0    -8      0     80       0  -192      0   128
//
//  Recursion:
//
//    U(0,X) = 1,
//    U(1,X) = 2*X,
//    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 February 2003
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
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Output, double C[(N+1)*((N+1)], the coefficients of the Chebyshev U
//    polynomials.
//
{
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  c[1+1*(n+1)] = 2.0;

  for ( i = 2; i <= n; i++ )
  {
    c[i+0*(n+1)]     =                    - c[i-2+0*(n+1)];
    for ( j = 1; j <= i-2; j++ )
    {
      c[i+j*(n+1)] = 2.0 * c[i-1+(j-1)*(n+1)] - c[i-2+j*(n+1)];
    }

    c[i+(i-1)*(n+1)] = 2.0 * c[i-1+(i-2)*(n+1)];
    c[i+ i   *(n+1)] = 2.0 * c[i-1+(i-1)*(n+1)];
  }

  return;
}
//****************************************************************************80

void cheby_u_poly_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_POLY_VALUES returns values of Chebyshev polynomials U(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ChebyshevU[n,x]
//
//    The Chebyshev U polynomial is a solution to the differential equation:
//
//    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0;
//
//  First terms:
//
//    U(0)(X) =   1
//    U(1)(X) =   2 X
//    U(2)(X) =   4 X^2 -   1
//    U(3)(X) =   8 X^3 -   4 X
//    U(4)(X) =  16 X^4 -  12 X^2 +  1
//    U(5)(X) =  32 X^5 -  32 X^3 +  6 X
//    U(6)(X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
//    U(7)(X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
//
//  Recursion:
//
//    U(0)(X) = 1,
//    U(1)(X) = 2 * X,
//    U(N)(X) = 2 * X * U(N-1)(X) - U(N-2)(X)
//
//  Norm:
//
//    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N)(X)^2 dX = PI/2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2012
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
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.1600000000000000E+01,
      0.1560000000000000E+01,
      0.8960000000000000E+00,
     -0.1264000000000000E+00,
     -0.1098240000000000E+01,
     -0.1630784000000000E+01,
     -0.1511014400000000E+01,
     -0.7868390400000000E+00,
      0.2520719360000000E+00,
      0.1190154137600000E+01,
      0.1652174684160000E+01,
      0.1453325357056000E+01 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *cheby_u_poly_zero ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_POLY_ZERO returns zeroes of Chebyshev polynomials U(n,x).
//
//  Discussion:
//
//    The I-th zero of U(N,X) is cos((I-1)*PI/(N-1)), I = 1 to N
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
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double CHEBY_U_POLY_ZERO[N], the zeroes of U(N,X).
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;
  double *z;

  z = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    angle = ( double) ( i ) * pi / ( double ) ( n + 1 );
    z[i-1] = cos ( angle );
  }

  return z;
}
//****************************************************************************80

void chebyshev_discrete ( int n, int m, double x, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials at a point.
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
//  Reference:
//
//    Walter Gautschi,
//    Orthogonal Polynomials: Computation and Approximation,
//    Oxford, 2004,
//    ISBN: 0-19-850672-4,
//    LC: QA404.5 G3555.
//
//  Parameters:
//
//    Input, int N, the highest order of the polynomials to 
//    be evaluated.  0 <= N <= M.
//
//    Input, int M, the maximum order of the polynomials.
//    0 <= M.
//
//    Input, double X, the evaluation point.
//
//    Output, double V[N+1], the value of the polynomials at X.
//
{
  int i;

  if ( m < 0 )
  {
    cerr << "\n";
    cerr << "CHEBYSHEV_DISCRETE - Fatal error!\n";
    cerr << "  Parameter M must be nonnegative.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "CHEBYSHEV_DISCRETE - Fatal error!\n";
    cerr << "  Parameter N must be nonnegative.\n";
    exit ( 1 );
  }

  if ( m < n )
  {
    cerr << "\n";
    cerr << "CHEBYSHEV_DISCRETE - Fatal error!\n";
    cerr << "  Parameter N must be no greater than M.\n";
    exit ( 1 );
  }

  v[0] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  v[1] = 2.0 * x + ( double ) ( 1 - m );

  if ( n == 1 )
  {
    return;
  }

  for ( i = 1; i < n; i++ )
  {
    v[i+1] = ( 
      ( double ) ( 2 * i + 1 ) 
      * ( 2.0 * x + ( double ) ( 1 - m ) ) * v[i]
      - ( double ) ( i * ( m + i ) * ( m - i ) ) * v[i-1]
    ) / ( double ) ( i + 1 );
  }

  return;
}
//****************************************************************************80

int collatz_count ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    COLLATZ_COUNT counts the number of terms in a Collatz sequence.
//
//  Discussion:
//
//    The rules for generation of the Collatz sequence are recursive.
//    If T is the current entry of the sequence, (T is
//    assumed to be a positive integer), then the next
//    entry, U is determined as follows:
//
//      if T is 1 (or less)
//        terminate the sequence;
//      else if T is even
//        U = T/2.
//      else (if T is odd and not 1)
//        U = 3*T+1;
//
//     N  Sequence                                                Length
//
//     1                                                               1
//     2   1                                                           2
//     3  10,  5, 16,  8,  4,  2,  1                                   8
//     4   2   1                                                       3
//     5  16,  8,  4,  2,  1                                           6
//     6   3, 10,  5, 16,  8,  4,  2,  1                               9
//     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
//     8   4,  2,  1                                                   4
//     9  28, 14,  7, ...                                             20
//    10   5, 16,  8,  4,  2,  1                                       7
//    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
//    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
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
//  Reference:
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input, int N, the first element of the sequence.
//
//    Output, int COLLATZ_COUNT, the number of elements in
//    the Collatz sequence that begins with N.
//
{
  int count;

  count = 1;

  for ( ; ; )
  {
    if ( n <= 1 )
    {
      break;
    }
    else if ( ( n % 2 ) == 0 )
    {
      n = n / 2;
    }
    else
    {
      n = 3 * n + 1;
    }
    count = count + 1;
  }

  return count;
}
//****************************************************************************80

void collatz_count_max ( int n, int *i_max, int *j_max )

//****************************************************************************80
//
//  Purpose:
//
//    COLLATZ_COUNT_MAX seeks the maximum Collatz count for 1 through N.
//
//  Discussion:
//
//    For each integer I, we compute a sequence of values that 
//    terminate when we reach 1.  The number of steps required to
//    reach 1 is the "rank" of I, and we are searching the numbers
//    from 1 to N for the number with maximum rank.
//
//    For a given I, the sequence is produced by:
//
//    1) J = 1, X(J) = I;
//    2) If X(J) = 1, stop.
//    3) J = J + 1; 
//       if X(J-1) was even, X(J) = X(J-1)/2;
//       else                X(J) = 3 * X(J-1) + 1;
//    4) Go to 3
//
//  Example:
//
//            N      I_MAX J_MAX
//
//           10          9    20
//          100         97   119
//        1,000        871   179
//       10,000      6,171   262
//      100,000     77,031   351
//    1,000,000    837,799   525
//   10,000,000  8,400,511   686
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
//  Parameters:
//
//    Input, int N, the maximum integer to check.
//
//    Output, int *I_MAX, *J_MAX, an integer I with the maximum rank,
//    and the value of the maximum rank.
//
{
  int i;
  int j;
  int x;

  *i_max = -1;
  *j_max = -1;

  for ( i = 1; i <= n; i++ )
  {
    j = 1;
    x = i;

    while ( x != 1 )
    {
      j = j + 1;
      if ( ( x % 2 ) == 0 )
      {
        x = x / 2;
      }
      else
      {
        x = 3 * x + 1;
      }
    }

    if ( *j_max < j )
    {
      *i_max = i;
      *j_max = j;
    }
  }

  return;
}
//****************************************************************************80

void collatz_count_values ( int *n_data, int *n, int *count )

//****************************************************************************80
//
//  Purpose:
//
//    COLLATZ_COUNT_VALUES returns some values of the Collatz count function.
//
//  Discussion:
//
//    The rules for generation of the Collatz sequence are recursive.
//    If T is the current entry of the sequence, (T is
//    assumed to be a positive integer), then the next
//    entry, U is determined as follows:
//   
//      if T is 1 (or less)
//        terminate the sequence;
//      else if T is even
//        U = T/2.
//      else (if T is odd and not 1)
//        U = 3*T+1;
//
//    The Collatz count is the length of the Collatz sequence for a given
//    starting value.  By convention, we include the initial value in the
//    count, so the minimum value of the count is 1.
//
//     N  Sequence                                                 Count
//
//     1                                                               1
//     2   1                                                           2
//     3  10,  5, 16,  8,  4,  2,  1                                   8
//     4   2   1                                                       3
//     5  16,  8,  4,  2,  1                                           6
//     6   3, 10,  5, 16,  8,  4,  2,  1                               9
//     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
//     8   4,  2,  1                                                   4
//     9  28, 14,  7, ...                                             20
//    10   5, 16,  8,  4,  2,  1                                       7
//    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
//    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *N, the initial value of a Collatz sequence.
//
//    Output, int *COUNT, the length of the Collatz sequence starting
//    with N.
//
{
# define N_MAX 20

  int count_vec[N_MAX] = {
      1,   2,   8,   3,   6,   9,   17,   4,  20,   7, 
    112,  25,  26,  27,  17,  28,  111,  18,  83,  29 };
  int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, 
     27,  50, 100, 200, 300, 400, 500, 600, 700, 800 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *count = 0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *count = count_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void comb_row ( bool next, int n, int row[] )

//****************************************************************************80
//
//  Purpose:
//
//    COMB_ROW computes row N of Pascal's triangle.
//
//  Discussion:
//
//    Row N contains the N+1 combinatorial coefficients
//
//      C(N,0), C(N,1), C(N,2), ... C(N,N)
//
//  Discussion:
//
//    The sum of the elements of row N is equal to 2^N.
//
//    The formula is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  First terms:
//
//     N K:0  1   2   3   4   5   6   7  8  9 10
//
//     0   1
//     1   1  1
//     2   1  2   1
//     3   1  3   3   1
//     4   1  4   6   4   1
//     5   1  5  10  10   5   1
//     6   1  6  15  20  15   6   1
//     7   1  7  21  35  35  21   7   1
//     8   1  8  28  56  70  56  28   8  1
//     9   1  9  36  84 126 126  84  36  9  1
//    10   1 10  45 120 210 252 210 120 45 10  1
//
//  Recursion:
//
//    C(N,K) = C(N-1,K-1)+C(N-1,K)
//
//  Special values:
//
//    C(N,0) = C(N,N) = 1
//    C(N,1) = C(N,N-1) = N
//    C(N,N-2) = sum ( 1 <= I <= N ) N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, bool NEXT, indicates whether this is a call for
//    the 'next' row of the triangle.
//    NEXT = FALSE means this is a startup call.  Row N is desired, but
//    presumably this is a first call, or row N-1 was not computed
//    on the previous call.
//    NEXT = TRUE means this is not the first call, and row N-1 was computed
//    on the previous call.  In this case, much work can be saved
//    by using the information from the previous values of ROW
//    to build the next values.
//
//    Input, int N, the row of the triangle desired.  The triangle
//    begins with row 0.
//
//    Output, int ROW[N+1], the row of coefficients.
//    ROW(I) = C(N,I-1).
//
{
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  if ( next )
  {
    for ( i = n-1; 1 <= i; i-- )
    {
      row[i] = row[i] + row[i-1];
    }

    row[n] = 1;
  }
  else
  {
    row[0] = 1;
    for ( i = 1; i <= n; i++ )
    {
      row[i] = 0;
    }

    for ( j = 1; j <= n; j++ )
    {
      for ( i = j; 1 <= i; i-- )
      {
        row[i] = row[i] + row[i-1];
      }
    }

  }

  return;
}
//****************************************************************************80

int commul ( int iarray[], int n, int nfact )

//****************************************************************************80
//
//  Purpose:
//
//    COMMUL computes a multinomial combinatorial coefficient.
//
//  Discussion:
//
//    The multinomial coefficient is a generalization of the binomial
//    coefficient.  It may be interpreted as the number of combinations of
//    N objects, where IARRAY(1) objects are indistinguishable of type 1,
//    ... and IARRAY(K) are indistinguishable of type NFACT.
//
//    The formula is:
//
//      COMMUL = N! / ( IARRAY(1)! IARRAY(2)! ... IARRAY(NFACT)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IARRAY(NFACT).
//    IARRAY contains the NFACT values used in the denominator.
//    Note that the sum of these entries should be N,
//    and that all entries should be nonnegative.
//
//    Input, inte N, determines the numerator.
//
//    Input, int NFACT, the number of factors in the numerator.
//
//    Output, int COMMUL, the value of the multinomial coefficient.
//
{
  double arg;
  double fack;
  double facn;
  int i;
  int isum;

  for ( i = 0; i < nfact; i++ )
  {
    if ( iarray[i] < 0 )
    {
      cerr << "\n";
      cerr << "COMMUL - Fatal error\n";
      cerr << "  Entry " << i << " of IARRAY = " << iarray[i] << "\n";
      cerr << "  But this value must be nonnegative.\n";
      exit ( 1 );
    }
  }

  isum = 0;
  for ( i = 0; i < nfact; i++ )
  {
    isum = isum + iarray[i];
  }

  if ( isum != n )
  {
    cerr << "\n";
    cerr << "COMMUL - Fatal error!\n";
    cerr << "  The sum of the IARRAY entries is " << isum << "\n";
    cerr << "  But it must equal N = " << n << "\n";
    exit ( 1 );
  }

  arg = ( double ) ( n + 1 );
  facn = lgamma ( arg );

  for ( i = 0; i < nfact; i++ )
  {
    arg = ( double ) ( iarray[i] + 1 );
    fack = lgamma ( arg );
    facn = facn - fack;
  }

  return ( r8_nint ( exp ( facn ) ) );
}
//****************************************************************************80

double cos_deg ( double angle )

//****************************************************************************80
//
//  Purpose:
//
//    COS_DEG returns the cosine of an angle given in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, the angle, in degrees.
//
//    Output, double COS_DEG, the cosine of the angle.
//
{
  double degrees_to_radians = 3.141592653589793 / 180.0;
  double value;

  value = cos ( degrees_to_radians * angle );

  return value;
}
//****************************************************************************80

double cos_power_int ( double a, double b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    COS_POWER_INT evaluates the cosine power integral.
//
//  Discussion:
//
//    The function is defined by
//
//      COS_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( cos ( t ))^n dt
//
//    The algorithm uses the following fact:
//
//      Integral cos^n ( t ) = -(1/n) * (
//        cos^(n-1)(t) * sin(t) + ( n-1 ) * Integral cos^(n-2) ( t ) dt )
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
//  Parameters
//
//    Input, double A, B, the limits of integration.
//
//    Input, integer N, the power of the sine function.
//
//    Output, double COS_POWER_INT, the value of the integral.
//
{
  double ca;
  double cb;
  int m;
  int mlo;
  double sa;
  double sb;
  double value;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "COS_POWER_INT - Fatal error!\n";
    cerr << "  Power N < 0.\n";
    exit ( 1 );
  }

  sa = sin ( a );
  sb = sin ( b );
  ca = cos ( a );
  cb = cos ( b );

  if ( ( n % 2 ) == 0 )
  {
    value = b - a;
    mlo = 2;
  }
  else
  {
    value = sb - sa;
    mlo = 3;
  }

  for ( m = mlo; m <= n; m = m + 2 )
  {
    value = ( ( double ) ( m - 1 ) * value 
      - pow ( ca, (m-1) ) * sa + pow ( cb, (m-1) ) * sb ) 
      / ( double ) ( m );
  }

  return value;
}
//****************************************************************************80

void cos_power_int_values ( int &n_data, double &a, double &b, int &n,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    COS_POWER_INT_VALUES returns some values of the sine power integral.
//
//  Discussion:
//
//    The function has the form
//
//      COS_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( cos(T) )^N dt
//
//    In Mathematica, the function can be evaluated by:
//
//      Integrate [ ( Cos[x] )^n, { x, a, b } ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the limits of integration.
//
//    Output, int &N, the power.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00 };

  static double b_vec[N_MAX] = {
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793 };

  static double fx_vec[N_MAX] = {
     3.141592653589793, 
     0.0, 
     1.570796326794897, 
     0.0, 
     1.178097245096172, 
     0.0, 
     0.9817477042468104, 
     0.0, 
     0.8590292412159591, 
     0.0, 
     0.7731263170943632 };

  static int n_vec[N_MAX] = {
     0,
     1,
     2,
     3,
     4,
     5,
     6,
     7,
     8,
     9,
    10 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    n = 0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    n = n_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double e_constant ( )

//****************************************************************************80
//
//  Purpose:
//
//    E_CONSTANT returns the value of the base of the natural logarithm system.
//
//  Definition:
//
//    E = Limit ( N -> +oo ) ( 1 + 1 / N )^N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double E_CONSTANT, the base of the natural logarithm system.
//
{
  return 2.718281828459045235360287;
}
//****************************************************************************80

void erf_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    ERF_VALUES returns some values of the ERF or "error" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  double bvec[N_MAX] = {
    0.0000000000, 0.1124629160, 0.2227025892, 0.3286267595, 
    0.4283923550, 0.5204998778, 0.6038560908, 0.6778011938, 
    0.7421009647, 0.7969082124, 0.8427007929, 0.8802050696, 
    0.9103139782, 0.9340079449, 0.9522851198, 0.9661051465, 
    0.9763483833, 0.9837904586, 0.9890905016, 0.9927904292, 
    0.9953222650 };
  double xvec[N_MAX] = {
    0.0, 0.1, 0.2, 0.3, 
    0.4, 0.5, 0.6, 0.7, 
    0.8, 0.9, 1.0, 1.1, 
    1.2, 1.3, 1.4, 1.5, 
    1.6, 1.7, 1.8, 1.9, 
    2.0 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = xvec[*n_data];
    *fx = bvec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double error_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    ERROR_F evaluates the error function ERF(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 May 2007
//
//  Author:
//
//    Original FORTRAN77 version by William Cody.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    "Rational Chebyshev approximations for the error function",
//    Mathematics of Computation, 
//    1969, pages 631-638.
//
//  Parameters:
//
//    Input, double X, the argument of the error function.
//
//    Output, double ERROR_F, the value of ERF(X).
//
{
  double a[5] = {
    3.16112374387056560, 
    1.13864154151050156E+02, 
    3.77485237685302021E+02, 
    3.20937758913846947E+03, 
    1.85777706184603153E-01 };
  double b[4] = {
    2.36012909523441209E+01, 
    2.44024637934444173E+02, 
    1.28261652607737228E+03, 
    2.84423683343917062E+03 };
  double c[9] = {
    5.64188496988670089E-01, 
    8.88314979438837594, 
    6.61191906371416295E+01, 
    2.98635138197400131E+02, 
    8.81952221241769090E+02, 
    1.71204761263407058E+03, 
    2.05107837782607147E+03, 
    1.23033935479799725E+03, 
    2.15311535474403846E-08 };
  double d[8] = {
    1.57449261107098347E+01, 
    1.17693950891312499E+02, 
    5.37181101862009858E+02, 
    1.62138957456669019E+03, 
    3.29079923573345963E+03, 
    4.36261909014324716E+03, 
    3.43936767414372164E+03, 
    1.23033935480374942E+03 };
  double del;
  double erfx;
  int i;
  double p[6] = {
    3.05326634961232344E-01, 
    3.60344899949804439E-01, 
    1.25781726111229246E-01, 
    1.60837851487422766E-02, 
    6.58749161529837803E-04, 
    1.63153871373020978E-02 };
  double q[5] = {
    2.56852019228982242, 
    1.87295284992346047, 
    5.27905102951428412E-01, 
    6.05183413124413191E-02, 
    2.33520497626869185E-03 };
  double sqrpi = 0.56418958354775628695;
  double thresh = 0.46875;
  double pxabs;
  double xbig = 26.543;
  double pxden;
  double xabs;
  double xden;
  double xnum;
  double xsmall = 1.11E-16;
  double xsq;

  xabs = fabs ( x );
//
//  Evaluate ERF(X) for |X| <= 0.46875.
//
  if ( xabs <= thresh )
  {
    if ( xsmall < xabs )
    {
      xsq = xabs * xabs;
    }
    else
    {
      xsq = 0.0;
    }

    xnum = a[4] * xsq;
    xden = xsq;
    for ( i = 0; i < 3; i++ )
    {
      xnum = ( xnum + a[i] ) * xsq;
      xden = ( xden + b[i] ) * xsq;
    }

    erfx = x * ( xnum + a[3] ) / ( xden + b[3] );
  }
//
//  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
//
  else if ( xabs <= 4.0 )
  {
    xnum = c[8] * xabs;
    xden = xabs;
    for ( i = 0; i < 7; i++ )
    {
      xnum = ( xnum + c[i] ) * xabs;
      xden = ( xden + d[i] ) * xabs;
    }

    erfx = ( xnum + c[7] ) / ( xden + d[7] );
    xsq = ( double ) ( ( int ) ( ( xabs * 16.0 ) / 16.0 ) );
    del = ( xabs - xsq ) * ( xabs + xsq );
    erfx = exp ( - xsq * xsq ) * exp ( - del ) * erfx;

    erfx = ( 0.5 - erfx ) + 0.5;

    if ( x < 0.0 )
    {
      erfx = - erfx;
    }
  }
//
//  Evaluate ERFC(X) for 4.0 < |X|.
//
  else
  {
    if ( xbig <= xabs )
    {
      if ( 0.0 < x )
      {
        erfx = 1.0;
      }
      else
      {
        erfx = - 1.0;
      }
    }
    else
    {
      xsq = 1.0 / ( xabs * xabs );

      xnum = p[5] * xsq;
      xden = xsq;
      for ( i = 0; i < 4; i++ )
      {
        xnum = ( xnum + p[i] ) * xsq;
        xden = ( xden + q[i] ) * xsq;
      }

      erfx = xsq * ( xnum + p[4] ) / ( xden + q[4] );
      erfx = ( sqrpi - erfx ) / xabs;
      xsq = ( double ) ( ( int ) ( ( xabs * 16.0 ) / 16.0 ) );
      del = ( xabs - xsq ) * ( xabs + xsq );
      erfx = exp ( - xsq * xsq ) * exp ( - del ) * erfx;

      erfx = ( 0.5 - erfx ) + 0.5;
      if ( x < 0.0 )
      {
        erfx = - erfx;
      }

    }

  }

  return erfx;
}
//****************************************************************************80

double error_f_inverse ( double y )

//****************************************************************************80
//
//  Purpose:
//
//    ERROR_F_INVERSE inverts the error function ERF.
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
//  Parameters:
//
//    Input, double Y, the value of the error function.
//
//    Output, double ERROR_F_INVERSE, the value X such that
//    ERROR_F(X) = Y.
//
{
  double value;
  double x;
  double z;

  z = ( y + 1.0 ) / 2.0;

  x = normal_01_cdf_inv ( z );

  value = x / sqrt ( 2.0 );

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
//      Gamma = limit ( M -> +oo )
//        ( sum ( 1 <= N <= M ) 1 / N ) - log ( M )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
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
  return ( 0.577215664901532860606512090082402431042 );
}
//****************************************************************************80

void euler_number ( int n, int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_NUMBER computes the Euler numbers.
//
//  Discussion:
//
//    The Euler numbers can be evaluated in Mathematica with the call
//
//      EulerE[n]
//
//    These numbers rapidly get too big to store in an ordinary integer!
//
//    The terms of odd index are 0.
//
//    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
//
//  First terms:
//
//    E0  = 1
//    E1  = 0
//    E2  = -1
//    E3  = 0
//    E4  = 5
//    E5  = 0
//    E6  = -61
//    E7  = 0
//    E8  = 1385
//    E9  = 0
//    E10 = -50521
//    E11 = 0
//    E12 = 2702765
//    E13 = 0
//    E14 = -199360981
//    E15 = 0
//    E16 = 19391512145
//    E17 = 0
//    E18 = -2404879675441
//    E19 = 0
//    E20 = 370371188237525
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input, int N, the index of the last Euler number to compute.
//
//    Output, int E[N+1], the Euler numbers from index 0 to N.
//
{
  int cnk;
  int i;
  int j;
  int sgn;

  if ( n < 0 )
  {
    return;
  }

  e[0] = 1;

  if ( n == 0 )
  {
    return;
  }

  e[1] = 0;

  if ( n == 1 )
  {
    return;
  }

  e[2] = -1;

  for ( i = 3; i <= n; i++ )
  {
    e[i] = 0;

    if ( ( i % 2 ) == 0 )
    {
      for ( j = 2; j <= i; j = j + 2 )
      {
        e[i] = e[i] - i4_choose ( i, j ) * e[i-j];
      }
    }
  }
  return;
}
//****************************************************************************80

double euler_number2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_NUMBER2 computes the Euler numbers.
//
//  Discussion:
//
//    The Euler numbers can be evaluated in Mathematica with the call
//
//      EulerE[n]
//
//  First terms:
//
//    E0  = 1
//    E1  = 0
//    E2  = -1
//    E3  = 0
//    E4  = 5
//    E5  = 0
//    E6  = -61
//    E7  = 0
//    E8  = 1385
//    E9  = 0
//    E10 = -50521
//    E11 = 0
//    E12 = 2702765
//    E13 = 0
//    E14 = -199360981
//    E15 = 0
//    E16 = 19391512145
//    E17 = 0
//    E18 = -2404879675441
//    E19 = 0
//    E20 = 370371188237525
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input, int N, the index of the Euler number to compute.
//
//    Output, double EULER_NUMBER2, the value of E(N).
//
{
  double e[7] = { 1.0, -1.0, 5.0, -61.0, 1385.0, 
       -50521.0, 2702765.0 };
  int i;
  int itmax = 1000;
  double sum1;
  double term;
  double value;

  if ( n < 0 )
  {
    return 0.0;
  }

  if ( n == 0 )
  {
    return e[0];
  }

  if ( ( n % 2 ) == 1 )
  {
    return 0.0;
  }

  if ( n <= 12 )
  {
    return e[n/2];
  }

  sum1 = 0.0;

  for ( i = 1; i <= itmax; i++ )
  {
    term = 1.0 / pow ( ( double ) ( 2 * i - 1 ), n + 1 );

    if ( ( i % 2 ) == 1 )
    {
      sum1 = sum1 + term;
    }
    else
    {
      sum1 = sum1 - term;
    }

    if ( fabs ( term ) < 1.0E-10 )
    {
      break;
    }
    else if ( fabs ( term ) < 1.0E-08 * fabs ( sum1 ) )
    {
      break;
    }

  }

  value = pow ( 2.0, n + 2 ) * sum1 * r8_factorial ( n ) 
    / pow ( r8_pi ( ), n + 1 );

  if ( ( n % 4 ) != 0 )
  {
    value = -value;
  }

  return value;
}
//****************************************************************************80

void euler_number_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_NUMBER_VALUES returns some values of the Euler numbers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the order of the Euler number.
//
//    Output, int *C, the value of the Euler number.
//
{
# define N_MAX 8
 
  int c_vec[N_MAX] = {
    1, 0, -1, 5, 61, 1385, -50521, 2702765 };

  int n_vec[N_MAX] = {
     0, 1, 2, 4, 6, 8, 10, 12 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double euler_poly ( int n, double x )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_POLY evaluates the N-th Euler polynomial at X.
//
//  First values:
//
//    E(0,X) = 1
//    E(1,X) = X - 1/2
//    E(2,X) = X^2 - X
//    E(3,X) = X^3 - 3/2 X^2 + 1/4
//    E(4,X) = X^4 - 2*X^3 + X
//    E(5,X) = X^5 - 5/2 X^4 + 5/2 X^2 - 1/2
//    E(6,X) = X^6 - 3 X^5 + 5 X^3 - 3 X
//    E(7,X) = X^7 - 7/2 X^6 + 35/4 X^4 - 21/2 X^2 + 17/8
//    E(8,X) = X^8 - 4 X^7 + 14 X^5 - 28 X^3 + 17 X
//
//  Special values:
//
//    E'(N,X) = N * E(N-1,X)
//
//    E(N,1/2) = E(N) / 2^N, where E(N) is the N-th Euler number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the Euler polynomial to
//    be evaluated.  N must be 0 or greater.
//
//    Input, double X, the value at which the polynomial is to
//    be evaluated.
//
//    Output, double EULER_POLY, the value of E(N,X).
//
{
  double bx1;
  double bx2;
  double value;

  bx1 = bernoulli_poly2 ( n+1, x );
  bx2 = bernoulli_poly2 ( n+1, 0.5 * x );

  value = 2.0 * ( bx1 - bx2 * pow ( ( double ) 2, n+1 ) ) 
    / ( double ) ( n + 1 );

  return value;
}
//****************************************************************************80

void eulerian ( int n, int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    EULERIAN computes the Eulerian number E(N,K).
//
//  Definition:
//
//    A run in a permutation is a sequence of consecutive ascending values.
//
//    E(N,K) is the number of permutations of N objects which contain
//    exactly K runs.
//
//  Examples:
//
//     N = 7
//
//     1     0     0     0     0     0     0
//     1     1     0     0     0     0     0
//     1     4     1     0     0     0     0
//     1    11    11     1     0     0     0
//     1    26    66    26     1     0     0
//     1    57   302   302    57     1     0
//     1   120  1191  2416  1191   120     1
//
//  Recursion:
//
//    E(N,K) = K * E(N-1,K) + (N-K+1) * E(N-1,K-1).
//
//  Properties:
//
//    E(N,1) = E(N,N) = 1.
//    E(N,K) = 0 if K <= 0 or N < K.
//    sum ( 1 <= K <= N ) E(N,K) = N!.
//    X^N = sum ( 0 <= K <= N ) COMB(X+K-1, N ) E(N,K)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dennis Stanton and Dennis White,
//    Constructive Combinatorics,
//    Springer Verlag, 1986
//
//  Parameters:
//
//    Input, int N, the number of rows desired.
//
//    Output, int E[N*N], the first N rows of Eulerian numbers.
//
{
  int i;
  int j;

  if ( n < 1 )
  {
    return;
  }
//
//  Construct rows 1, 2, ..., N of the Eulerian triangle.
//
  e[1-1+(1-1)*n] = 1;
  for ( j = 2; j <= n; j++ )
  {
    e[1-1+(j-1)*n] = 0;
  }

  for ( i = 2; i <= n; i++ )
  {
    e[i-1+(1-1)*n] = 1;
    for ( j = 2; j <= n; j++ )
    {
      e[i-1+(j-1)*n] = j * e[i-2+(j-1)*n] + ( i - j + 1 ) * e[i-2+(j-2)*n];
    }
  }

  return;
}
//****************************************************************************80

int f_hofstadter ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    F_HOFSTADTER computes the Hofstadter F sequence.
//
//  Discussion:
//
//    F(N) = 0                if N = 0
//         = N - F ( N - 1 ), otherwise.
//
//    F(N) is defined for all nonnegative integers, and turns out
//    to be equal to int ( ( N + 1 ) / 2 ).
//
//  Table:
//
//     N  F(N)
//    --  ----
//
//     0     0
//     1     1
//     2     1
//     3     2
//     4     2
//     5     3
//     6     3
//     7     4
//     8     4
//     9     5
//    10     5
//    11     6
//    12     6
//    13     7
//    14     7
//    15     8
//    16     8
//    17     9
//    18     9
//    19    10
//    20    10
//    21    11
//    22    11
//    23    12
//    24    12
//    25    13
//    26    13
//    27    14
//    28    14
//    29    15
//    30    15
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Douglas Hofstadter,
//    Goedel, Escher, Bach,
//    Basic Books, 1979.
//
//  Parameters:
//
//    Input, int N, the argument of the function.
//
//    Output, int F_HOFSTADTER, the value of the function.
//
{
  if ( n <= 0 )
  {
    return 0;
  }
  else
  {
    return ( n - f_hofstadter ( n-1 ) );
  }
}
//****************************************************************************80

int fibonacci_direct ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_DIRECT computes the N-th Fibonacci number directly.
//
//  Discussion:
//
//    The formula is:
//
//      F(N) = ( PHIP**N - PHIM**N ) / sqrt(5)
//
//    where
//
//      PHIP = ( 1 + sqrt(5) ) / 2,
//      PHIM = ( 1 - sqrt(5) ) / 2.
//
//  Example:
//
//     N   F
//    --  --
//     0   0
//     1   1
//     2   1
//     3   2
//     4   3
//     5   5
//     6   8
//     7  13
//     8  21
//     9  34
//    10  55
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the index of the Fibonacci number to compute.
//    N should be nonnegative.
//
//    Output, int FIBONACCI_DIRECT, the value of the N-th Fibonacci number.
//
{
  int f;
  double sqrt5 = 2.236068;
  double phim = ( 1.0 - sqrt5 ) / 2.0;
  double phip = ( 1.0 + sqrt5 ) / 2.0;

  if ( n < 0 )
  {
    f = 0;
  }
  else
  {
    f = r8_nint ( ( pow ( phip, n ) - pow ( phim, n ) ) / sqrt ( 5.0 ) );
  }

  return f;
}
//****************************************************************************80

void fibonacci_floor ( int n, int *f, int *i )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_FLOOR returns the largest Fibonacci number less or equal to N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the positive integer whose Fibonacci "floor" is desired.
//
//    Output, int *F, the largest Fibonacci number less than or equal to N.
//
//    Output, int *I, the index of the F.
//
{
  if ( n <= 0 )
  {
    *i = 0;
    *f = 0;
  }
  else
  {
    *i = ( int ) ( 
        log ( 0.5 * ( double ) ( 2 * n + 1 ) * sqrt ( 5.0 ) )
      / log ( 0.5 * ( 1.0 + sqrt ( 5.0 ) ) ) );

    *f = fibonacci_direct ( *i );

    if ( n < *f )
    {
      *i = *i - 1;
      *f = fibonacci_direct ( *i );
    }

  }

  return;
}
//****************************************************************************80

void fibonacci_recursive ( int n, int f[] )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_RECURSIVE computes the first N Fibonacci numbers.
//
//  Discussion:
//
//    The 'golden ratio' PHI = (1+sqrt(5))/2 satisfies the equation
//
//      X*X-X-1=0
//
//    which is often written as:
//
//       X        1
//      --- =  ------
//       1      X - 1
//
//    expressing the fact that a rectangle, whose sides are in proportion X:1,
//    is similar to the rotated rectangle after a square of side 1 is removed.
//
//      <----X---->
//
//      +-----*---*
//      |     |   |  1
//      |     |   |
//      +-----*---+
//      <--1-><X-1>
//
//    The formula is:
//
//      PHIP = ( 1 + sqrt(5) ) / 2
//      PHIM = ( 1 - sqrt(5) ) / 2
//      F(N) = ( PHIP^N + PHIM^N ) / sqrt(5)
//
//    Moreover, F(N) can be computed by computing PHIP**N / sqrt(5) and rounding
//    to the nearest whole number.
//
//  First terms:
//
//      1
//      1
//      2
//      3
//      5
//      8
//     13
//     21
//     34
//     55
//     89
//    144
//
//    The 40th number is                  102,334,155.
//    The 50th number is               12,586,269,025.
//    The 100th number is 354,224,848,179,261,915,075.
//
//  Recursion:
//
//    F(1) = 1
//    F(2) = 1
//
//    F(N) = F(N-1) + F(N-2)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the highest Fibonacci number to compute.
//
//    Output, int F[N], the first N Fibonacci numbers.
//
{
  int i;

  if ( n <= 0 )
  {
    return;
  }

  f[0] = 1;

  if ( n <= 1 )
  {
    return;
  }

  f[1] = 1;

  for ( i = 2; i < n; i++ )
  {
    f[i] = f[i-1] + f[i-2];
  }

  return;
}
//****************************************************************************80

int g_hofstadter ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    G_HOFSTADTER computes the Hofstadter G sequence.
//
//  Discussion:
//
//    G(N) = 0                      if N = 0
//         = N - G ( G ( N - 1 ) ), otherwise.
//
//    G(N) is defined for all nonnegative integers.
//
//    The value of G(N) turns out to be related to the Zeckendorf
//    representation of N as a sum of non-consecutive Fibonacci numbers.
//    To compute G(N), determine the Zeckendorf representation:
//
//      N = sum ( 1 <= I <= M ) F(I)
//
//    and reduce the index of each Fibonacci number by 1:
//
//      G(N) = sum ( 1 <= I <= M ) F(I-1)
//
//    However, this is NOT how the computation is done in this routine.
//    Instead, a straightforward recursive function call is defined
//    to correspond to the definition of the mathematical function.
//
//  Table:
//
//     N  G(N)  Zeckendorf   Decremented
//    --  ----  ----------   -----------
//
//     1   1    1            1
//     2   1    2            1
//     3   2    3            2
//     4   3    3 + 1        2 + 1
//     5   3    5            3
//     6   4    5 + 1        3 + 1
//     7   4    5 + 2        3 + 1
//     8   5    8            5
//     9   6    8 + 1        5 + 1
//    10   6    8 + 2        5 + 1
//    11   7    8 + 3        5 + 2
//    12   8    8 + 3 + 1    5 + 2 + 1
//    13   8    13           8
//    14   9    13 + 1       8 + 1
//    15   9    13 + 2       8 + 1
//    16  10    13 + 3       8 + 2
//    17  11    13 + 3 + 1   8 + 2 + 1
//    18  11    13 + 5       8 + 3
//    19  12    13 + 5 + 1   8 + 3 + 1
//    20  12    13 + 5 + 2   8 + 3 + 1
//    21  13    21           13
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Douglas Hofstadter,
//    Goedel, Escher, Bach,
//    Basic Books, 1979.
//
//  Parameters:
//
//    Input, int N, the argument of the function.
//
//    Output, int G_HOFSTADTER, the value of the function.
//
{
  if ( n <= 0 )
  {
    return 0;
  }
  else
  {
    return ( n - g_hofstadter ( g_hofstadter ( n-1 ) ) );
  }

}
//****************************************************************************80

void gamma_log_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LOG_VALUES returns some values of the Log Gamma function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA isincremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 18

  double bvec[N_MAX] = {
     1.524064183,    0.7966780066,   0.3982337117,  
     0.1520599127,   0.000000000,   -0.04987246543, 
    -0.08537410945, -0.1081747934,  -0.1196128950,  
    -0.1207822040,  -0.1125917658,  -0.09580771625, 
    -0.07108385116, -0.03898428380,  0.000000000,   
    12.80182743,    39.33988571,    71.25704193 };
  double xvec[N_MAX] = {
    0.2,  0.4,  0.6,  0.8, 
    1.0,  1.1,  1.2,  1.3, 
    1.4,  1.5,  1.6,  1.7, 
    1.8,  1.9,  2.0, 10.0, 
   20.0, 30.0 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = xvec[*n_data];
    *fx = bvec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_VALUES returns some values of the Gamma function.
//
//  Discussion:
//
//    The Gamma function is defined as:
//
//      Gamma(Z) = Integral ( 0 <= T < +oo) T^(Z-1) exp(-T) dT
//
//    It satisfies the recursion:
//
//      Gamma(X+1) = X * Gamma(X)
//
//    Gamma is undefined for nonpositive integral X.
//    Gamma(0.5) = sqrt(PI)
//    For N a positive integer, Gamma(N+1) = Factorial ( N ).
//
//    In Mathematica, the function can be evaluated by:
//
//      Gamma[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 May 2007
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
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 25

  double fx_vec[N_MAX] = { 
     -0.3544907701811032E+01,  
     -0.1005871979644108E+03,  
      0.9943258511915060E+02,  
      0.9513507698668732E+01,  
      0.4590843711998803E+01,  
      0.2218159543757688E+01,  
      0.1772453850905516E+01,  
      0.1489192248812817E+01,  
      0.1164229713725303E+01,  
      0.1000000000000000E+01,  
      0.9513507698668732E+00,  
      0.9181687423997606E+00,  
      0.8974706963062772E+00,  
      0.8872638175030753E+00,  
      0.8862269254527580E+00,  
      0.8935153492876903E+00,  
      0.9086387328532904E+00,  
      0.9313837709802427E+00,  
      0.9617658319073874E+00,  
      0.1000000000000000E+01,  
      0.2000000000000000E+01,  
      0.6000000000000000E+01,  
      0.3628800000000000E+06,  
      0.1216451004088320E+18,  
      0.8841761993739702E+31 };

  double x_vec[N_MAX] = { 
     -0.50E+00,  
     -0.01E+00,  
      0.01E+00,  
      0.10E+00,  
      0.20E+00,  
      0.40E+00,  
      0.50E+00,  
      0.60E+00,  
      0.80E+00,  
      1.00E+00,  
      1.10E+00,  
      1.20E+00,  
      1.30E+00,  
      1.40E+00,  
      1.50E+00,  
      1.60E+00,  
      1.70E+00,  
      1.80E+00,  
      1.90E+00,  
      2.00E+00,  
      3.00E+00,  
      4.00E+00,  
     10.00E+00,  
     20.00E+00,  
     30.00E+00 }; 

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gegenbauer_poly ( int n, double alpha, double x, double cx[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_POLY computes the Gegenbauer polynomials C(0:N,ALPHA,X).
//
//  Discussion:
//
//    The Gegenbauer polynomial can be evaluated in Mathematica with
//    the command
//
//      GegenbauerC[n,m,x]
//
//  Differential equation:
//
//    (1-X*X) Y'' - (2 ALPHA + 1) X Y' + N (N + 2 ALPHA) Y = 0
//
//  Recursion:
//
//    C(0,ALPHA,X) = 1,
//    C(1,ALPHA,X) = 2*ALPHA*X
//    C(N,ALPHA,X) = ( 2*(N-1+ALPHA)*X * C(N-1,ALPHA,X) - (N-2+2*ALPHA) * C(N-2,ALPHA,X) )/N
//
//  Restrictions:
//
//    ALPHA must be greater than -0.5.
//
//  Special values:
//
//    If ALPHA = 1, the Gegenbauer polynomials reduce to the Chebyshev
//    polynomials of the second kind.
//
//  Norm:
//
//    Integral ( -1 <= X <= 1 )
//      ( 1 - X^2 )^( ALPHA - 0.5 ) * C(N,ALPHA,X)^2 dX
//
//    = PI * 2^( 1 - 2 * ALPHA ) * Gamma ( N + 2 * ALPHA )
//      / ( N! * ( N + ALPHA ) * ( Gamma ( ALPHA ) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, double ALPHA, a parameter which is part of the definition of
//    the Gegenbauer polynomials.  It must be greater than -0.5.
//
//    Input, double X, the point at which the polynomials are to be evaluated.
//
//    Output, double CX[N+1], the values of the first N+1 Gegenbauer
//    polynomials at the point X.
//
{
  int i;

  if ( alpha <= -0.5 )
  {
    cerr << "\n";
    cerr << "GEGENBAUER_POLY - Fatal error!\n";
    cerr << "  Illegal value of ALPHA = " <<alpha << "\n";
    cerr << "  but ALPHA must be greater than -0.5.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    return;
  }

  cx[0] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  cx[1] = 2.0 * alpha * x;

  for ( i = 2; i <= n; i++ )
  {
    cx[i] = ( ( ( double ) ( 2 * i - 2 ) + 2.0 * alpha ) * x * cx[i-1] 
            + ( ( double ) (   - i + 2 ) - 2.0 * alpha )     * cx[i-2] ) 
            /   ( double )       i;
  }

  return;
}
//****************************************************************************80

void gegenbauer_poly_values ( int *n_data, int *n, double *a, double *x, 
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_POLY_VALUES returns some values of the Gegenbauer polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 August 2004
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
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *N, the order parameter of the function.
//
//    Output, double A, the real parameter of the function.
//
//    Output, double X, the argument of the function.
//
//    Output, double FX, the value of the function.
//
{
# define N_MAX 38

  double a_vec[N_MAX] = {
     0.5,  0.5,  0.5, 
     0.5,  0.5,  0.5, 
     0.5,  0.5,  0.5, 
     0.5,  0.5,  0.0, 
     1.0,  2.0,  3.0, 
     4.0,  5.0,  6.0, 
     7.0,  8.0,  9.0, 
    10.0,  3.0,  3.0, 
     3.0,  3.0,  3.0, 
     3.0,  3.0,  3.0, 
     3.0,  3.0,  3.0, 
     3.0,  3.0,  3.0, 
     3.0,  3.0 };
  double fx_vec[N_MAX] = {
    1.0000000000,   0.2000000000,  -0.4400000000, 
   -0.2800000000,   0.2320000000,   0.3075200000, 
   -0.0805760000,  -0.2935168000,  -0.0395648000, 
    0.2459712000,   0.1290720256,   0.0000000000, 
   -0.3600000000,  -0.0800000000,   0.8400000000, 
    2.4000000000,   4.6000000000,   7.4400000000, 
   10.9200000000,  15.0400000000,  19.8000000000, 
   25.2000000000,  -9.0000000000,  -0.1612800000, 
   -6.6729600000,  -8.3750400000,  -5.5267200000, 
    0.0000000000,   5.5267200000,   8.3750400000, 
    6.6729600000,   0.1612800000,  -9.0000000000, 
  -15.4252800000,  -9.6969600000,  22.4409600000, 
  100.8892800000, 252.0000000000 };
  int n_vec[N_MAX] = {
     0,  1,  2, 
     3,  4,  5, 
     6,  7,  8, 
     9, 10,  2, 
     2,  2,  2, 
     2,  2,  2, 
     2,  2,  2, 
     2,  5,  5, 
     5,  5,  5, 
     5,  5,  5, 
     5,  5,  5, 
     5,  5,  5, 
     5,  5 };
  double x_vec[N_MAX] = {
    0.20,  0.20,  0.20, 
    0.20,  0.20,  0.20, 
    0.20,  0.20,  0.20, 
    0.20,  0.20,  0.40, 
    0.40,  0.40,  0.40, 
    0.40,  0.40,  0.40, 
    0.40,  0.40,  0.40, 
    0.40, -0.50, -0.40, 
   -0.30, -0.20, -0.10, 
    0.00,  0.10,  0.20, 
    0.30,  0.40,  0.50, 
    0.60,  0.70,  0.80, 
    0.90,  1.00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *a = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//******************************************************************************

void gen_hermite_poly ( int n, double x, double mu, double p[] )

//******************************************************************************
//
//  Purpose:
//
//    GEN_HERMITE_POLY evaluates the generalized Hermite polynomials at X.
//
//  Discussion:
//
//    The generalized Hermite polynomials are orthogonal under the weight
//    function:
//
//      w(x) = |x|^(2*MU) * exp ( - x^2 )
//
//    over the interval (-oo,+oo).
//
//    When MU = 0, the generalized Hermite polynomial reduces to the standard
//    Hermite polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Theodore Chihara,
//    An Introduction to Orthogonal Polynomials,
//    Gordon and Breach, 1978,
//    ISBN: 0677041500,
//    LC: QA404.5 C44.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//
//    Input, double X, the point at which the polynomials are 
//    to be evaluated.
//
//    Input, double MU, the parameter.
//    - 1 / 2 < MU.
//
//    Output, double P[N+1], the values of the first N+1
//    polynomials at the point X.
//
{
  int i;
  double theta;

  if ( n < 0 )
  {
    return;
  }

  p[0] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  p[1] = 2.0 * x;
 
  for ( i = 1; i < n; i++ )
  {
    if ( ( i % 2 ) == 0 )
    {
      theta = 0.0;
    }
    else
    {
      theta = 2.0 * mu;
    }

    p[i+1] = 2.0 * x * p[i] - 2.0 * ( ( double ) ( i ) + theta ) * p[i-1];
  }
 
  return;
}
//****************************************************************************80

void gen_laguerre_poly ( int n, double alpha, double x, double cx[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_POLY evaluates generalized Laguerre polynomials.
//
//  Differential equation:
//
//    X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0
//
//  Recursion:
//
//    L(0,ALPHA,X) = 1
//    L(1,ALPHA,X) = 1+ALPHA-X
//
//    L(N,ALPHA,X) = ( (2*N-1+ALPHA-X) * L(N-1,ALPHA,X)
//                   - (N-1+ALPHA) * L(N-2,ALPHA,X) ) / N
//
//  Restrictions:
//
//    -1 < ALPHA
//
//  Special values:
//
//    For ALPHA = 0, the generalized Laguerre polynomial L(N,ALPHA,X)
//    is equal to the Laguerre polynomial L(N,X).
//
//    For ALPHA integral, the generalized Laguerre polynomial
//    L(N,ALPHA,X) equals the associated Laguerre polynomial L(N,ALPHA,X).
//
//  Norm:
//
//    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,ALPHA,X)^2 dX
//    = Gamma ( N + ALPHA + 1 ) / N!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2010
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
//  Parameters:
//
//    Input, int N, the highest order function to compute.
//
//    Input, double ALPHA, the parameter.  -1 < ALPHA is required.
//
//    Input, double X, the point at which the functions are to be
//    evaluated.
//
//    Output, double CX[N+1], the polynomials of
//    degrees 0 through N evaluated at the point X.
//
{
  int i;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "GEN_LAGUERRE_POLY - Fatal error!\n";
    cerr << "  The input value of ALPHA is " << alpha << "\n";
    cerr << "  but ALPHA must be greater than -1.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    return;
  }

  cx[0] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  cx[1] = 1.0 + alpha - x;

  for ( i = 2; i <= n; i++ )
  {
    cx[i] = ( ( ( double ) ( 2 * i - 1 ) + alpha - x ) * cx[i-1] 
            + ( ( double ) (   - i + 1 ) - alpha     ) * cx[i-2] ) 
              / ( double )       i;
  }

  return;
}
//****************************************************************************80

double gud ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    GUD evaluates the Gudermannian function.
//
//  Definition:
//
//    The Gudermannian function relates the hyperbolic and trigonometric
//    functions.  For any argument X, there is a corresponding value
//    GAMMA so that
//
//      sinh(x) = tan(gamma).
//
//    The value GAMMA is called the Gudermannian of X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the Gudermannian.
//
//    Output, double GUD, the value of the Gudermannian.
//
{
  double value;

  value = 2.0 * atan ( tanh ( 0.5 * x ) );

  return value;
}
//****************************************************************************80

void gud_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GUD_VALUES returns some values of the Gudermannian function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 13

  double f_vec[N_MAX] = {
    -1.301760336,  -0.8657694832, 0.0000000000, 
     0.09983374879, 0.1986798470, 0.4803810791, 
     0.8657694832,  1.131728345,  1.301760336,  
     1.406993569,   1.471304341,  1.510419908,  
     1.534169144 };
  double x_vec[N_MAX] = {
    -2.0, -1.0,  0.0, 
     0.1,  0.2,  0.5, 
     1.0,  1.5,  2.0, 
     2.5,  3.0,  3.5, 
     4.0};

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data];
    *fx = f_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int h_hofstadter ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    H_HOFSTADTER computes the Hofstadter H sequence.
//
//  Discussion:
//
//    H(N) = 0                          if N = 0
//         = N - H ( H ( H ( N - 1 ) ), otherwise.
//
//    H(N) is defined for all nonnegative integers.
//
//  Table:
//
//     N  H(N)
//    --  ----
//
//     0     0
//     1     1
//     2     1
//     3     2
//     4     3
//     5     4
//     6     4
//     7     5
//     8     5
//     9     6
//    10     7
//    11     7
//    12     8
//    13     9
//    14    10
//    15    10
//    16    11
//    17    12
//    18    13
//    19    13
//    20    14
//    21    14
//    22    15
//    23    16
//    24    17
//    25    17
//    26    18
//    27    18
//    28    19
//    29    20
//    30    20
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Douglas Hofstadter,
//    Goedel, Escher, Bach,
//    Basic Books, 1979.
//
//  Parameters:
//
//    Input, int N, the argument of the function.
//
//    Output, int H_HOFSTADTER, the value of the function.
//
{
  if ( n <= 0 )
  {
    return 0;
  }
  else
  {
    return ( n - h_hofstadter ( h_hofstadter ( h_hofstadter ( n-1 ) ) ) );
  }

}
//****************************************************************************80

int hail ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HAIL computes the hail function.
//
//  Discussion:
//
//    Starting with a positive integer N, we divide it by 2 if it is
//    even, or triple and add 1 if odd, and repeat this process until
//    reaching the value 1.  The number of times the process is carried
//    out is the value of the hail function for the given starting value.
//
//    Actually, HAIL is not well defined, since it is not known if
//    the above process actually terminates, let alone terminates at 1,
//    for every starting value N.
//
//  Example:
//
//     N  Sequence                                                  Hail
//
//     1                                                               0
//     2   1                                                           1
//     3  10,  5, 16,  8,  4,  2,  1                                   7
//     4   2   1                                                       2
//     5  16,  8,  4,  2,  1                                           5
//     6   3, 10,  5, 16,  8,  4,  2,  1                               8
//     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   16
//     8   4,  2,  1                                                   3
//     9  28, 14,  7, ...                                             19
//    10   5, 16,  8,  4,  2,  1                                       6
//    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          14
//    12   6,  3, 10,  5, 16,  8,  4,  2,  1                           9
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the starting value for the hail sequence.
//
//    Output, int HAIL, the number of steps before the hail sequence
//    reached 1.
//
{
  int k;
  int m;

  k = 0;
  m = n;

  if ( 0 < n )
  {
    while ( m != 1 )
    {
      k = k + 1;
      if ( ( m % 2 ) == 0 )
      {
        m = m / 2;
      }
      else
      {
        m = 3 * m + 1;
      }
    }

  }

  return k;
}
//****************************************************************************80

void hermite_poly ( int n, double x, double cx[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLY evaluates the Hermite polynomials at X.
//
//  Differential equation:
//
//    Y'' - 2 X Y' + 2 N Y = 0
//
//  First terms:
//
//      1
//      2 X
//      4 X^2     -  2
//      8 X^3     - 12 X
//     16 X^4     - 48 X^2     + 12
//     32 X^5    - 160 X^3    + 120 X
//     64 X^6    - 480 X^4    + 720 X^2    - 120
//    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
//    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
//    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
//   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
//
//  Recursion:
//
//    H(0,X) = 1,
//    H(1,X) = 2*X,
//    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
//
//  Norm:
//
//    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
//    = sqrt ( PI ) * 2^N * N!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
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
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, double X, the point at which the polynomials are to be evaluated.
//
//    Output, double CX[N+1], the values of the first N+1 Hermite
//    polynomials at the point X.
//
{
  int i;

  if ( n < 0 )
  {
    return;
  }

  cx[0] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  cx[1] = 2.0 * x;

  for ( i = 2; i <= n; i++ )
  {
    cx[i] = 2.0 * x * cx[i-1] - 2.0 * ( double ) ( i - 1 ) * cx[i-2];
  }

  return;
}
//****************************************************************************80

void hermite_poly_coef ( int n, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLY_COEF evaluates the coefficients of the physicist's Hermite polynomial H(n,x).
//
//  First terms:
//
//    N/K     0     1      2      3       4     5      6    7      8    9   10
//
//     0      1
//     1      0     2
//     2     -2     0      4
//     3      0   -12      0      8
//     4     12     0    -48      0      16
//     5      0   120      0   -160       0    32
//     6   -120     0    720      0    -480     0     64
//     7      0 -1680      0   3360       0 -1344      0   128
//     8   1680     0 -13440      0   13440     0  -3584     0    256
//     9      0 30240      0 -80640       0 48384      0 -9216      0 512
//    10 -30240     0 302400      0 -403200     0 161280     0 -23040   0 1024
//
//  Recursion:
//
//    H(0,X) = 1,
//    H(1,X) = 2*X,
//    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2003
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
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Output, double C[(N+1)*(N+1)], the coefficients of the Hermite
//    polynomials.
//
{
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  c[1+1*(n+1)] = 2.0;

  for ( i = 2; i <= n; i++ )
  {
    c[i+0*(n+1)] = - 2.0 * ( double ) ( i - 1 ) * c[i-2+0*(n+1)];
    for ( j = 1; j <= i-2; j++ )
    {
      c[i+j*(n+1)] =   2.0 * c[i-1+(j-1)*(n+1)] 
                     - 2.0 * ( double ) ( i - 1 ) * c[i-2+j*(n+1)];
    }
    c[i+(i-1)*(n+1)] =  2.0 * c[i-1+(i-2)*(n+1)];
    c[i+ i   *(n+1)] =  2.0 * c[i-1+(i-1)*(n+1)];
  }

  return;
}
//****************************************************************************80

void hermite_poly_values ( int *n_data, int *n, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLY_VALUES returns some values of the Hermite polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the order of the polynomial.
//
//    Output, double *X, the point where the polynomial is evaluated.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 17

  double fx_vec[N_MAX] = {
     1.0,            10.0,           98.0, 
     940.0,          8812.0,         80600.0, 
     717880.0,       6211600.0,      520656800.0, 
     421271200,      3275529760.0,   24329873600.0, 
     171237081280.0, 41.0,          -8.0, 
     3816.0,         3041200.0 };
  int n_vec[N_MAX] = {
     0,  1,  2, 
     3,  4,  5, 
     6,  7,  8, 
     9, 10, 11, 
    12,  5,  5, 
     5,  5 };
  double x_vec[N_MAX] = {
    5.0,  5.0,  5.0, 
    5.0,  5.0,  5.0, 
    5.0,  5.0,  5.0, 
    5.0,  5.0,  5.0, 
    5.0,  0.5,  1.0, 
    3.0, 10.0 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data];
    *x = x_vec[*n_data];
    *fx = fx_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hyper_2f1_values ( int *n_data, double *a, double *b, double *c, 
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric function 2F1.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      fx = Hypergeometric2F1 [ a, b, c, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2007
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
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Shanjie Zhang, Jianming Jin,
//    Computation of Special Functions,
//    Wiley, 1996,
//    ISBN: 0-471-11963-6,
//    LC: QA351.C45
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, *B, *C, *X, the parameters of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 24

  double a_vec[N_MAX] = { 
   -2.5, 
   -0.5, 
    0.5, 
    2.5, 
   -2.5, 
   -0.5, 
    0.5, 
    2.5, 
   -2.5, 
   -0.5, 
    0.5, 
    2.5, 
    3.3, 
    1.1, 
    1.1, 
    3.3, 
    3.3, 
    1.1, 
    1.1, 
    3.3, 
    3.3, 
    1.1, 
    1.1, 
    3.3 };
  double b_vec[N_MAX] = { 
    3.3, 
    1.1, 
    1.1, 
    3.3, 
    3.3, 
    1.1, 
    1.1, 
    3.3, 
    3.3, 
    1.1, 
    1.1, 
    3.3, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7 };
  double c_vec[N_MAX] = { 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
    6.7, 
   -5.5, 
   -0.5, 
    0.5, 
    4.5, 
   -5.5, 
   -0.5, 
    0.5, 
    4.5, 
   -5.5, 
   -0.5, 
    0.5, 
    4.5 };
  double fx_vec[N_MAX] = { 
    0.72356129348997784913, 
    0.97911109345277961340, 
    1.0216578140088564160, 
    1.4051563200112126405, 
    0.46961431639821611095, 
    0.95296194977446325454, 
    1.0512814213947987916, 
    2.3999062904777858999, 
    0.29106095928414718320, 
    0.92536967910373175753, 
    1.0865504094806997287, 
    5.7381565526189046578, 
    15090.669748704606754, 
   -104.31170067364349677, 
    21.175050707768812938, 
    4.1946915819031922850, 
    1.0170777974048815592E+10, 
   -24708.635322489155868, 
    1372.2304548384989560, 
    58.092728706394652211, 
    5.8682087615124176162E+18, 
   -4.4635010147295996680E+08, 
    5.3835057561295731310E+06, 
    20396.913776019659426 };
  double x_vec[N_MAX] = { 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.55, 
    0.55, 
    0.55, 
    0.55, 
    0.85, 
    0.85, 
    0.85, 
    0.85, 
    0.25, 
    0.25, 
    0.25, 
    0.25, 
    0.55, 
    0.55, 
    0.55, 
    0.55, 
    0.85, 
    0.85, 
    0.85, 
    0.85 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }
 
  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0;
    *b = 0.0;
    *c = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *c = c_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int i4_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE computes the binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
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
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, are the values of N and K.
//
//    Output, int I4_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0;
  }
  else if ( mn == 0 )
  {
    value = 1;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
    }
  }

  return value;
}
//****************************************************************************80

void i4_factor ( int n, int maxfactor, int *nfactor, int factor[], 
  int power[], int *nleft )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTOR factors an integer into prime factors.
//
//  Discussion:
//
//    N = NLEFT * Product ( 1 <= I <= NFACTOR ) FACTOR(I)^POWER(I).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer to be factored.  N may be positive,
//    negative, or 0.
//
//    Input, int MAXFACTOR, the maximum number of prime factors for
//    which storage has been allocated.
//
//    Output, int *NFACTOR, the number of prime factors of N discovered
//    by the routine.
//
//    Output, int FACTOR[MAXFACTOR], the prime factors of N.
//
//    Output, int POWER[MAXFACTOR].  POWER(I) is the power of
//    the FACTOR(I) in the representation of N.
//
//    Output, int *NLEFT, the factor of N that the routine could not
//    divide out.  If NLEFT is 1, then N has been completely factored.
//    Otherwise, NLEFT represents factors of N involving large primes.
//
{
  int i;
  int maxprime;
  int p;

  *nfactor = 0;

  for ( i = 0; i < maxfactor; i++ )
  {
    factor[i] = 0;
  }

  for ( i = 0; i < maxfactor; i++ )
  {
    power[i] = 0;
  }

  *nleft = n;

  if ( n == 0 )
  {
    return;
  }

  if ( abs ( n ) == 1 )
  {
    *nfactor = 1;
    factor[0] = 1;
    power[0] = 1;
    return;
  }
//
//  Find out how many primes we stored.
//
  maxprime = prime ( -1 );
//
//  Try dividing the remainder by each prime.
//
  for ( i = 1; i <= maxprime; i++ )
  {
    p = prime ( i );

    if ( abs ( *nleft ) % p == 0 )
    {
      if ( *nfactor < maxfactor )
      {
        *nfactor = *nfactor + 1;
        factor[*nfactor-1] = p;

        for ( ; ; )
        {
          power[*nfactor-1] = power[*nfactor-1] + 1;
          *nleft = *nleft / p;

          if ( abs ( *nleft ) % p != 0 )
          {
            break;
          }

        }

        if ( abs ( *nleft ) == 1 )
        {
          break;
        }
      }
    }
  }

  return;
}
//****************************************************************************80

int i4_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//    0 <= N <= 13 is required.
//
//    Output, int I4_FACTORIAL, the factorial of N.
//
{
  int i;
  int value;

  value = 1;

  if ( 13 < n ) 
  {
    value = - 1;
    cerr << "I4_FACTORIAL - Fatal error!\n";
    cerr << "  I4_FACTORIAL(N) cannot be computed as an integer\n";
    cerr << "  for 13 < N.\n";
    cerr << "  Input value N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 1; i <= n; i++ )
  {
    value = value * i;
  }

  return value;
}
//****************************************************************************80

void i4_factorial_values ( int *n_data, int *n, int *fn )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL_VALUES returns values of the factorial function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to the index of the test data.  On each subsequent call, N_DATA is
//    incremented and that test data is returned.  When there is no more
//    test data, N_DATA is set to 0.
//
//    Output, int *N, the argument of the function.
//
//    Output, int *FN, the value of the function.
//
{
# define N_MAX 13

  int fnvec[N_MAX] = {
    1, 1, 2, 6, 
    24, 120, 720, 5040, 
    40320, 362880, 3628800, 39916800, 
    479001600 };
  int nvec[N_MAX] = {
     0,  1,  2,  3, 
     4,  5,  6,  7, 
     8,  9, 10, 11, 
    12 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *fn = 0;
  }
  else
  {
    *n = nvec[*n_data];
    *fn = fnvec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int i4_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL2 computes the double factorial function.
//
//  Discussion:
//
//    The formula is:
//
//      FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                      = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N    FACTORIAL2(N)
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the double factorial function.
//    If N is less than 1, I4_FACTORIAL2 is returned as 1.
//
//    Output, int I4_FACTORIAL2, the value of the double factorial function.
//
{
  int value;

  if ( n < 1 )
  {
    return 1;
  }

  value = 1;

  while ( 1 < n )
  {
    value = value * n;
    n = n - 2;
  }

  return value;
}
//****************************************************************************80

void i4_factorial2_values ( int *n_data, int *n, int *fn )

//****************************************************************************80
//
//  Purpose: 
//
//    I4_FACTORIAL2_VALUES returns values of the double factorial function.
//
//  Discussion:
//
//    The formula is:
//
//      FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                      = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N    FACTORIAL2(N)
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 August 2004
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
//    CRC Press, 1996, page 16.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *N, the argument of the function.
//
//    Output, int *FN, the value of the function.
//
{
# define N_MAX 16

  int fn_vec[N_MAX] = { 
        1, 
        1,     2,      3,      8,      15, 
       48,   105,    384,    945,    3840, 
    10395, 46080, 135135, 645120, 2027025 };
  int n_vec[N_MAX] = { 
     0, 
     1,  2,  3,  4,  5, 
     6,  7,  8,  9, 10, 
    11, 12, 13, 14, 15 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *fn = 0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *fn = fn_vec[*n_data-1];
  }
  return;
# undef N_MAX
}
//****************************************************************************80

int i4_gcd ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_GCD finds the greatest common divisor of two I4's.
//
//  Discussion:
//
//    Only the absolute values of I and J are
//    considered, so that the result is always nonnegative.
//
//    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
//
//    If I and J have no common factor, I4_GCD is returned as 1.
//
//    Otherwise, using the Euclidean algorithm, I4_GCD is the
//    largest common factor of I and J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, two numbers whose greatest common divisor
//    is desired.
//
//    Output, int I4_GCD, the greatest common divisor of I and J.
//
{
  int ip;
  int iq;
  int ir;
//
//  Return immediately if either I or J is zero.
//
  if ( i == 0 )
  {
    return i4_max ( 1, abs ( j ) );
  }
  else if ( j == 0 )
  {
    return i4_max ( 1, abs ( i ) );
  }
//
//  Set IP to the larger of I and J, IQ to the smaller.
//  This way, we can alter IP and IQ as we go.
//
  ip = i4_max ( abs ( i ), abs ( j ) );
  iq = i4_min ( abs ( i ), abs ( j ) );
//
//  Carry out the Euclidean algorithm.
//
  for ( ; ; )
  {
    ir = ip % iq;

    if ( ir == 0 )
    {
      break;
    }

    ip = iq;
    iq = ir;
  }

  return iq;
}
//****************************************************************************80

bool i4_is_prime ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_IS_PRIME reports whether an I4 is prime.
//
//  Method:
//
//    A simple, unoptimized sieve of Erasthosthenes is used to
//    check whether N can be divided by any integer between 2
//    and SQRT(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer to be tested.
//
//    Output, bool I4_IS_PRIME, is TRUE if N is prime, and FALSE
//    otherwise.  Note that negative numbers and 0 are not
//    considered prime.
//
{
  int i;
  int nhi;

  if ( n <= 0 )
  {
    return false;
  }

  if ( n <= 3 )
  {
    return true;
  }

  nhi = ( int ) ( sqrt ( ( double ) n ) );

  for ( i = 2; i <= nhi; i++ )
  {
    if ( ( n % i ) == 0 )
    {
      return false;
    }
  }

  return true;
}
//****************************************************************************80

bool i4_is_triangular ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_IS_TRIANGULAR determines whether an I4 is triangular.
//
//  Discussion:
//
//    The N-th triangular number is equal to the sum of the first
//    N integers.
//
//  First Values:
//
//    Index  Value
//     0      0
//     1      1
//     2      3
//     3      6
//     4     10
//     5     15
//     6     21
//     7     28
//     8     36
//     9     45
//    10     55
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer to be checked.
//
//    Output, bool I4_IS_TRIANGULAR, is TRUE if I is triangular.
//
{
  int j;
  int k;

  if ( i < 0 )
  {
    return false;
  }
  else if ( i == 0 )
  {
    return true;
  }
  else
  {
    i4_to_triangle ( i, &j, &k );

    if ( j == k )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

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
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int i1 and i2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of i1 and i2.
//
{
  if ( i1 > i2 ) 
  {
    return i1;
  }
  else 
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int i1 and i2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
//
{
  if ( i1 < i2 ) 
  {
    return i1;
  }
  else 
  {
    return i2;
  }

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
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//    The formula is:
//
//      If 
//        NREM = I4_MODP ( I, J ) 
//        NMULT = ( I - NREM ) / J
//      then
//        I = J * NMULT + NREM
//      where NREM is always nonnegative.
//
//  Example:
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
//    12 May 2003
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

int i4_partition_distinct_count ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_DISTINCT_COUNT returns any value of Q(N).
//
//  Discussion:
//
//    A partition of an integer N is a representation of the integer
//    as the sum of nonzero positive integers.  The order of the summands
//    does not matter.  The number of partitions of N is symbolized
//    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
//    following partitions:
//
//    5 = 5
//      = 4 + 1
//      = 3 + 2
//      = 3 + 1 + 1
//      = 2 + 2 + 1
//      = 2 + 1 + 1 + 1
//      = 1 + 1 + 1 + 1 + 1
//
//    However, if we require that each member of the partition
//    be distinct, we are computing something symbolized by Q(N).
//    The number 5 has Q(N) = 3, because it has the following partitions
//    into distinct parts:
//
//    5 = 5
//      = 4 + 1
//      = 3 + 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 February 2003
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
//  Parameters:
//
//    Input, int N, the integer to be partitioned.
//
//    Output, int I4_PARTITION_DISTINCT_COUNT, the number of partitions 
//    of the integer into distinct parts.
//
{
  int *c;
  int i;
  int k;
  int k2;
  int k_sign;
  int value;

  c = new int[n+1];

  c[0] = 1;

  for ( i = 1; i <= n; i++ )
  {
    if ( i4_is_triangular ( i ) )
    {
      c[i] = 1;
    }
    else
    {
      c[i] = 0;
    }

    k = 0;
    k_sign = -1;

    for ( ; ; )
    {
      k = k + 1;
      k_sign = -k_sign;
      k2 = k * ( 3 * k + 1 );

      if ( i < k2 )
      {
        break;
      }

      c[i] = c[i] + k_sign * c[i-k2];

    }

    k = 0;
    k_sign = -1;

    for ( ; ; )
    {
      k = k + 1;
      k_sign = -k_sign;
      k2 = k * ( 3 * k - 1 );

      if ( i < k2 )
      {
        break;
      }

      c[i] = c[i] + k_sign * c[i-k2];

    }

  }

  value = c[n];

  delete [] c;

  return value;
}
//****************************************************************************80

int i4_pochhammer ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POCHHAMMER returns the value of ( I * (I+1) * ... * (J-1) * J ).
//
//  Discussion:
//
//    Pochhammer's symbol (A)_N is the value
//
//      (A)_N = Gamma ( A + N ) / Gamma ( A )
//
//    or, for integer arguments,
//
//      (I)_N = I * ( I + 1 ) * ... * ( I + N - 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, inte I, J, values that define the product.
//
//    Output, int I4_POCHHAMMER, the value of the product.
//
{
  int k;
  int value;

  value = 1;

  for ( k = i; k <= j; k++ )
  {
    value = value * k;
  }

  return value;
}
//****************************************************************************80

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an integer.
//
//  Discussion:
//
//    The sign of 0 and all positive integers is taken to be +1.
//    The sign of all negative integers is -1.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
  if ( i < 0 ) 
  {
    return (-1);
  }
  else
  {
    return 1;
  }

}
//****************************************************************************80

void i4_swap ( int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two integer values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
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

void i4_to_triangle ( int k, int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_TRIANGLE converts an integer to triangular coordinates.
//
//  Discussion:
//
//    Triangular coordinates are handy when storing a naturally triangular
//    array (such as the lower half of a matrix) in a linear array.
//
//    Thus, for example, we might consider storing
//
//    (0,0)
//    (1,0) (1,1)
//    (2,0) (2,1) (2,2)
//    (3,0) (3,1) (3,2) (3,3)
//
//    as the linear array
//
//    (0,0) (1,0) (1,1) (2,0) (2,1) (2,2) (3,0) (3,1) (3,2) (3,3)
//
//    Here, the quantities in parenthesis represent the natural row and
//    column indices of a single number when stored in a rectangular array.
//
//    In this routine, we are given the location K of an item in the
//    linear array, and wish to determine the row I and column J
//    of the item when stored in the triangular array.
//
//  Example:
//
//    K  I  J
//
//    0  0  0
//    1  1  0
//    2  1  1
//    3  2  0
//    4  2  1
//    5  2  2
//    6  3  0
//    7  3  1
//    8  3  2
//    9  3  3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int K, the linear index of the (I,J) element, which
//    must be nonnegative.
//
//    Output, int *I, *J, the row and column indices.
//
{
  int c;
  int r;

  if ( k < 0 )
  {
    cerr << "\n";
    cerr << "I4_TO_TRIANGLE - Fatal error!\n";
    cerr << "  K < 0.\n";
    cerr << "  K = " << k << "\n";
    exit ( 1 );
  }
  else if ( k == 0 )
  {
    *i = 0;
    *j = 0;
    return;
  }
//
//   ( N - 1 )^2 + ( N - 1 ) < 2 * K <= N^2 + N
//
  r = ( int ) ( sqrt ( ( double ) ( 2 * ( k + 1 ) ) ) );

  if ( r * r + r < 2 * ( k + 1 ) )
  {
    r = r + 1;
  }

  r = r - 1;

  c = k - ( r * ( r + 1 ) ) / 2;

  *i = r;
  *j = c;

  return;
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

void i4mat_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT prints an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
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
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT_SOME prints some of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
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
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
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
    cout << "  Col:";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << "  " << setw(6) << j - 1;
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
//  Print out (up to INCX) entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ":";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << "  " << setw(6) << a[i-1+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *jacobi_poly ( int n, double alpha, double beta, double x )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_POLY evaluates the Jacobi polynomials at X.
//
//  Differential equation:
//
//    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
//
//  Recursion:
//
//    P(0,ALPHA,BETA,X) = 1,
//
//    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
//
//    P(N,ALPHA,BETA,X)  =
//      (
//        (2*N+ALPHA+BETA-1)
//        * ((ALPHA^2-BETA^2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X)
//        * P(N-1,ALPHA,BETA,X)
//        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
//      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
//
//  Restrictions:
//
//    -1 < ALPHA
//    -1 < BETA
//
//  Norm:
//
//    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA
//      * P(N,ALPHA,BETA,X)^2 dX
//    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
//      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
//
//  Special values:
//
//    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
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
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.  Note
//    that polynomials 0 through N will be computed.
//
//    Input, double ALPHA, one of the parameters defining the Jacobi
//    polynomials, ALPHA must be greater than -1.
//
//    Input, double BETA, the second parameter defining the Jacobi
//    polynomials, BETA must be greater than -1.
//
//    Input, double X, the point at which the polynomials are to be evaluated.
//
//    Output, double JACOBI_POLY[N+1], the values of the first N+1 Jacobi
//    polynomials at the point X.
//
{
  double c1;
  double c2;
  double c3;
  double c4;
  double *cx;
  int i;

  if ( alpha <= -1.0 )
  {
    cerr << "\n";
    cerr << "JACOBI_POLY - Fatal error!\n";
    cerr << "  Illegal input value of ALPHA = " << alpha << "\n";
    cerr << "  But ALPHA must be greater than -1.\n";
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    cerr << "\n";
    cerr << "JACOBI_POLY - Fatal error!\n";
    cerr << "  Illegal input value of BETA = " << beta << "\n";
    cerr << "  But BETA must be greater than -1.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    return NULL;
  }

  cx = new double[n+1];

  cx[0] = 1.0;

  if ( n == 0 )
  {
    return cx;
  }

  cx[1] = ( 1.0 + 0.5 * ( alpha + beta ) ) * x 
    + 0.5 * ( alpha - beta );

  for ( i = 2; i <= n; i++ )
  {
    c1 = 2.0 * ( double ) ( i ) * ( ( double ) ( i ) + alpha + beta ) 
      * ( ( double ) ( 2 * i - 2 ) + alpha + beta );

    c2 = ( ( double ) ( 2 * i - 1 ) + alpha + beta ) 
      * ( ( double ) ( 2 * i ) + alpha + beta ) 
      * ( ( double ) ( 2 * i - 2 ) + alpha + beta );

    c3 = ( ( double ) ( 2 * i - 1 ) + alpha + beta ) 
      * ( alpha + beta ) * ( alpha - beta );

    c4 = - ( double ) ( 2 ) * ( ( double ) ( i - 1 ) + alpha ) 
      * ( ( double ) ( i - 1 ) + beta )  
      * ( ( double ) ( 2 * i ) + alpha + beta );

    cx[i] = ( ( c3 + c2 * x ) * cx[i-1] + c4 * cx[i-2] ) / c1;
  }

  return cx;
}
//****************************************************************************80

void jacobi_poly_values ( int &n_data, int &n, double &a, double &b, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_POLY_VALUES returns some values of the Jacobi polynomial.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      JacobiP[ n, a, b, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
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
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the degree of the polynomial.
//
//    Output, double &A, &B, parameters of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 26

  static int a_vec[N_MAX] = {
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 1.0, 2.0,
     3.0, 4.0, 5.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0 };

  static int b_vec[N_MAX] = {
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 2.0,
    3.0, 4.0, 5.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0 };

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.2500000000000000E+00,
     -0.3750000000000000E+00,
     -0.4843750000000000E+00,
     -0.1328125000000000E+00,
      0.2753906250000000E+00,
     -0.1640625000000000E+00,
     -0.1174804687500000E+01,
     -0.2361328125000000E+01,
     -0.2616210937500000E+01,
      0.1171875000000000E+00,
      0.4218750000000000E+00,
      0.5048828125000000E+00,
      0.5097656250000000E+00,
      0.4306640625000000E+00,
     -0.6000000000000000E+01,
      0.3862000000000000E-01,
      0.8118400000000000E+00,
      0.3666000000000000E-01,
     -0.4851200000000000E+00,
     -0.3125000000000000E+00,
      0.1891200000000000E+00,
      0.4023400000000000E+00,
      0.1216000000000000E-01,
     -0.4396200000000000E+00,
      0.1000000000000000E+01 };

  static int n_vec[N_MAX] = {
     0, 1, 2, 3,
     4, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5 };

  static double x_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
     -1.0E+00,
     -0.8E+00,
     -0.6E+00,
     -0.4E+00,
     -0.2E+00,
      0.0E+00,
      0.2E+00,
      0.4E+00,
      0.6E+00,
      0.8E+00,
      1.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int jacobi_symbol ( int q, int p )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_SYMBOL evaluates the Jacobi symbol (Q/P).
//
//  Discussion:
//
//    If P is prime, then
//
//      Jacobi Symbol (Q/P) = Legendre Symbol (Q/P)
//
//    Else
//
//      let P have the prime factorization
//
//        P = Product ( 1 <= I <= N ) P(I)^E(I)
//
//      Jacobi Symbol (Q/P) =
//
//        Product ( 1 <= I <= N ) Legendre Symbol (Q/P(I))^E(I)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2000
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 86-87.
//
//  Parameters:
//
//    Input, int Q, an integer whose Jacobi symbol with
//    respect to P is desired.
//
//    Input, int P, the number with respect to which the Jacobi
//    symbol of Q is desired.  P should be 2 or greater.
//
//    Output, int JACOBI_SYMBOL, the Jacobi symbol (Q/P).
//    Ordinarily, L will be -1, 0 or 1.
//    If JACOBI_SYMBOL is -10, an error occurred.
//
{
# define FACTOR_MAX 20

  int factor[FACTOR_MAX];
  int i;
  int l;
  int nfactor;
  int nleft;
  int power[FACTOR_MAX];
  int value;
//
//  P must be greater than 1.
//
  if ( p <= 1 )
  {
    cerr << "\n";
    cerr << "JACOBI_SYMBOL - Fatal error!\n";
    cerr << "  P must be greater than 1.\n";
    exit ( 1 );
  }
//
//  Decompose P into factors of prime powers.
//
  i4_factor ( p, FACTOR_MAX, &nfactor, factor, power, &nleft );

  if ( nleft != 1 )
  {
    cerr << "\n";
    cerr << "JACOBI_SYMBOL - Fatal error!\n";
    cerr << "  Not enough factorization space.\n";
    exit ( 1 );
  }
//
//  Force Q to be nonnegative.
//
  while ( q < 0 )
  {
    q = q + p;
  }
//
//  For each prime factor, compute the Legendre symbol, and
//  multiply the Jacobi symbol by the appropriate factor.
//
  value = 1;

  for ( i = 0; i < nfactor; i++ )
  {
    l = legendre_symbol ( q, factor[i] );

    if ( l < -1 )
    {
      cerr << "\n";
      cerr << "JACOBI_SYMBOL - Fatal error!\n";
      cerr << "  Error during Legendre symbol calculation.\n";
      exit ( 1 );
    }
    value = value * ( int ) pow ( ( double ) l, power[i] );
  }

  return value;
# undef FACTOR_MAX
}
//****************************************************************************80

void krawtchouk ( int n, double p, double x, int m, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    KRAWTCHOUK evaluates the Krawtchouk polynomials at X.
//
//  Discussion:
//
//    The polynomial has a parameter P, which must be striclty between
//    0 and 1, and a parameter M which must be a nonnegative integer.
//
//    The Krawtchouk polynomial of order N, with parameters P and M,
//    evaluated at X, may be written K(N,P,X,M).
//
//    The first two terms are:
//
//      K(0,P,X,M) = 1
//      K(1,P,X,M) = X - P * M
//
//    and the recursion, for fixed P and M is
//
//                             ( N + 1 ) * K(N+1,P,X,M) =
//        ( X - ( N + P * ( M - 2 * N))) * K(N,  P,X,M)
//       - ( M - N + 1 ) * P * ( 1 - P ) * K(N-1,P,X,M)
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
//  Reference:
//
//    Walter Gautschi,
//    Orthogonal Polynomials: Computation and Approximation,
//    Oxford, 2004,
//    ISBN: 0-19-850672-4,
//    LC: QA404.5 G3555.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to evaluate.
//    0 <= N.
//
//    Input, double P, the parameter.  0 < P < 1.
//
//    Input, double X, the evaluation parameter.
//
//    Input, int M, the parameter.  0 <= M.
//
//    Output, double V[N+1], the values of the Krawtchouk polynomials
//    of orders 0 through N at X.
//
{
  int i;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "KRAWTCHOUK - Fatal error!\n";
    cerr << "  0 <= N is required.\n";
    exit ( 1 );
  }

  if ( p <= 0.0 || 1.0 <= p )
  {
    cerr << "\n";
    cerr << "KRAWTCHOUK - Fatal error!\n";
    cerr << "  0 < P < 1 is required.\n";
    exit ( 1 );
  }

  if ( m < 0 )
  {
    cerr << "\n";
    cerr << "KRAWTCHOUK - Fatal error!\n";
    cerr << "  0 <= M is required.\n";
    exit ( 1 );
  }

  v[0] = 1.0;

  if ( 1 <= n )
  {
    v[1] = x - p * ( double ) ( m );
  }

  for ( i = 1; i < n; i++ )
  {
    v[i+1] = ( 
      ( x - ( ( double ) ( i ) + p * ( double ) ( m - 2 * i ) ) ) * v[i] 
      - ( double ) ( m - i + 1 ) * p * ( 1.0 - p ) * v[i-1]
      ) / ( double ) ( i + 1 );
  }

  return;
}
//****************************************************************************80

void laguerre_associated ( int n, int m, double x, double cx[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_ASSOCIATED evaluates the associated Laguerre polynomials L(N,M,X) at X.
//
//  Differential equation:
//
//    X Y'' + (M+1-X) Y' + (N-M) Y = 0
//
//  First terms:
//
//    M = 0
//
//    L(0,0,X) =   1
//    L(1,0,X) =  -X    +  1
//    L(2,0,X) =   X^2 -  4 X     +  2
//    L(3,0,X) =  -X^3 +  9 X^2 -  18 X    +    6
//    L(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +      24
//    L(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120
//    L(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
//
//    M = 1
//
//    L(0,1,X) =    0
//    L(1,1,X) =   -1,
//    L(2,1,X) =    2 X - 4,
//    L(3,1,X) =   -3 X^2 + 18 X - 18,
//    L(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
//
//    M = 2
//
//    L(0,2,X) =    0
//    L(1,2,X) =    0,
//    L(2,2,X) =    2,
//    L(3,2,X) =   -6 X + 18,
//    L(4,2,X) =   12 X^2 - 96 X + 144
//
//    M = 3
//
//    L(0,3,X) =    0
//    L(1,3,X) =    0,
//    L(2,3,X) =    0,
//    L(3,3,X) =   -6,
//    L(4,3,X) =   24 X - 96
//
//    M = 4
//
//    L(0,4,X) =    0
//    L(1,4,X) =    0
//    L(2,4,X) =    0
//    L(3,4,X) =    0
//    L(4,4,X) =   24
//
//  Recursion:
//
//    if N = 0:
//
//      L(N,M,X)   = 0
//
//    if N = 1:
//
//      L(N,M,X)   = (M+1-X)
//
//    if N => 2:
//
//      L(N,M,X)   = ( (M+2*N-1-X) * L(N-1,M,X)
//                  +   (1-M-N)     * L(N-2,M,X) ) / N
//
//  Special values:
//
//    For M = 0, the associated Laguerre polynomials L(N,M,X) are equal
//    to the Laguerre polynomials L(N,X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2003
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
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, int M, the parameter.  M must be nonnegative.
//
//    Input, double X, the point at which the polynomials are to be evaluated.
//
//    Output, double CX[N+1], the associated Laguerre polynomials of
//    degrees 0 through N evaluated at the point X.
//
{
  int i;
  int ifact;

  if ( m < 0 )
  {
    cerr << "\n";
    cerr << "LAGUERRE_ASSOCIATED - Fatal error!\n";
    cerr << "  Input value of M = " << m << "\n";
    cerr << "  but M must be nonnegative.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    return;
  }

  cx[0] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  cx[1] = ( double ) ( m + 1 ) - x;

  for ( i = 2; i <= n; i++ )
  {
    cx[i] = ( ( ( double ) ( 2 * i + m - 1 ) - x ) * cx[i-1]
              + ( double ) (   - i - m + 1 )       * cx[i-2] ) 
              / ( double )       i;
  }

  return;
}
//****************************************************************************80

void laguerre_poly ( int n, double x, double cx[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLY evaluates the Laguerre polynomials at X.
//
//  Differential equation:
//
//    X * Y'' + (1-X) * Y' + N * Y = 0
//
//  First terms:
//
//      1
//     -X    +  1
//   (  X^2 -  4 X     +  2 ) / 2
//   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
//   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
//   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
//   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
//   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
//      + 52920 X^2 - 35280 X + 5040 ) / 5040
//
//  Recursion:
//
//    L(0,X) = 1,
//    L(1,X) = 1-X,
//    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
//
//  Orthogonality:
//
//    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,X) * L(M,X) dX
//    = 0 if N /= M
//    = 1 if N == M
//
//  Special values:
//
//    L(N,0) = 1.
//
//  Relations:
//
//    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2003
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
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, double X, the point at which the polynomials are to be evaluated.
//
//    Output, double CX[N+1], the Laguerre polynomials of degree 0 through
//    N evaluated at the point X.
//
{
  int i;

  if ( n < 0 )
  {
    return;
  }

  cx[0] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  cx[1] = 1.0 - x;

  for ( i = 2; i <= n; i++ )
  {
    cx[i] = ( ( ( double ) ( 2 * i - 1 ) - x ) * cx[i-1] 
              + ( double ) (   - i + 1 )       * cx[i-2] ) 
              / ( double ) (     i     );

  }

  return;
}
//****************************************************************************80

void laguerre_poly_coef ( int n, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLY_COEF evaluates the Laguerre polynomial coefficients.
//
//  First terms:
//
//    0: 1
//    1: 1  -1
//    2: 1  -2  1/2
//    3: 1  -3  3/2  1/6
//    4: 1  -4  4   -2/3  1/24
//    5: 1  -5  5   -5/3  5/24  -1/120
//
//  Recursion:
//
//    L(0) = ( 1,  0, 0, ..., 0 )
//    L(1) = ( 1, -1, 0, ..., 0 )
//    L(N) = (2*N-1-X) * L(N-1) - (N-1) * L(N-2) / N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2003
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
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Output, double C[(N+1)*(N+1)], the coefficients of the Laguerre 
//    polynomials of degree 0 through N.  Each polynomial is stored as a row.
//
{
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  for ( i = 0; i <= n; i++ )
  {
    c[i+0*(n+1)] = 1.0;
  }

  if ( n == 0 )
  {
    return;
  }

  c[1+1*(n+1)] = -1.0;

  for ( i = 2; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c[i+j*(n+1)] = ( 
          ( double ) ( 2 * i - 1 ) * c[i-1+j*(n+1)] 
        + ( double ) (   - i + 1 ) * c[i-2+j*(n+1)] 
        -                          c[i-1+(j-1)*(n+1)] ) 
        / ( double )       i;
    }
 
  }

  return;
}
//****************************************************************************80

void laguerre_polynomial_values ( int *n_data, int *n, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_VALUES returns some values of the Laguerre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the order of the polynomial.
//
//    Output, double *X, the point where the polynomial is evaluated.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 17

  double fx_vec[N_MAX] = {
     1.0000000000,  0.0000000000, -0.5000000000, 
    -0.6666666667, -0.6250000000, -0.4666666667, 
    -0.2569444444, -0.0404761905,  0.1539930556, 
     0.3097442681,  0.4189459325,  0.4801341791, 
     0.4962122235, -0.4455729167,  0.8500000000, 
    -3.1666666667, 34.3333333333 };
  int n_vec[N_MAX] = {
     0,  1,  2, 
     3,  4,  5, 
     6,  7,  8, 
     9, 10, 11, 
    12,  5,  5, 
     5,  5 };
  double x_vec[N_MAX] = {
    1.0,  1.0,  1.0, 
    1.0,  1.0,  1.0, 
    1.0,  1.0,  1.0, 
    1.0,  1.0,  1.0, 
    1.0,  0.5,  3.0, 
    5.0, 10.0 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data];
    *x = x_vec[*n_data];
    *fx = fx_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_associated ( int n, int m, double x, double cx[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED evaluates the associated Legendre functions.
//
//  Differential equation:
//
//    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
//
//  First terms:
//
//    M = 0  ( = Legendre polynomials of first kind P(N,X) )
//
//    P00 =    1
//    P10 =    1 X
//    P20 = (  3 X^2 -   1)/2
//    P30 = (  5 X^3 -   3 X)/2
//    P40 = ( 35 X^4 -  30 X^2 +   3)/8
//    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
//    P60 = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
//    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
//
//    M = 1
//
//    P01 =   0
//    P11 =   1 * SQRT(1-X*X)
//    P21 =   3 * SQRT(1-X*X) * X
//    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
//    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
//
//    M = 2
//
//    P02 =   0
//    P12 =   0
//    P22 =   3 * (1-X*X)
//    P32 =  15 * (1-X*X) * X
//    P42 = 7.5 * (1-X*X) * (7*X*X-1)
//
//    M = 3
//
//    P03 =   0
//    P13 =   0
//    P23 =   0
//    P33 =  15 * (1-X*X)^1.5
//    P43 = 105 * (1-X*X)^1.5 * X
//
//    M = 4
//
//    P04 =   0
//    P14 =   0
//    P24 =   0
//    P34 =   0
//    P44 = 105 * (1-X*X)^2
//
//  Recursion:
//
//    if N < M:
//      P(N,M) = 0
//    if N = M:
//      P(N,M) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
//      all the odd integers less than or equal to N.
//    if N = M+1:
//      P(N,M) = X*(2*M+1)*P(M,M)
//    if M+1 < N:
//      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
//
//  Special values:
//
//    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
//    function of the first kind equals the Legendre polynomial of the
//    first kind.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2005
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
//  Parameters:
//
//    Input, int N, the maximum first index of the Legendre
//    function, which must be at least 0.
//
//    Input, int M, the second index of the Legendre function,
//    which must be at least 0, and no greater than N.
//
//    Input, double X, the point at which the function is to be
//    evaluated.  X must satisfy -1 <= X <= 1.
//
//    Output, double CX[N+1], the values of the first N+1 function.
//
{
  double factor;
  int i;
  double somx2;

  if ( m < 0 )
  {
    cerr << "\n";
    cerr << "LEGENDRE_ASSOCIATED - Fatal error!\n";
    cerr << "  Input value of M is " << m << "\n";
    cerr << "  but M must be nonnegative.\n";
    exit ( 1 );
  }

  if ( n < m )
  {
    cerr << "\n";
    cerr << "LEGENDRE_ASSOCIATED - Fatal error!\n";
    cerr << "  Input value of M = " << m << "\n";
    cerr << "  Input value of N = " << n << "\n";
    cerr << "  but M must be less than or equal to N.\n";
    exit ( 1 );
  }

  if ( x < -1.0 )
  {
    cerr << "\n";
    cerr << "LEGENDRE_ASSOCIATED - Fatal error!\n";
    cerr << "  Input value of X = " << x << "\n";
    cerr << "  but X must be no less than -1.\n";
    exit ( 1 );
  }

  if ( 1.0 < x )
  {
    cerr << "\n";
    cerr << "LEGENDRE_ASSOCIATED - Fatal error!\n";
    cerr << "  Input value of X = " << x << "\n";
    cerr << "  but X must be no more than 1.\n";
    exit ( 1 );
  }

  for ( i = 0; i <= m-1; i++ )
  {
    cx[i] = 0.0;
  }
  cx[m] = 1.0;

  somx2 = sqrt ( 1.0 - x * x );

  factor = 1.0;
  for ( i = 1; i <= m; i++ )
  {
    cx[m] = -cx[m] * factor * somx2;
    factor = factor + 2.0;
  }

  if ( m == n )
  {
    return;
  }

  cx[m+1] = x * ( double ) ( 2 * m + 1 ) * cx[m];

  for ( i = m+2; i <= n; i++ )
  {
    cx[i] = ( ( double ) ( 2 * i     - 1 ) * x * cx[i-1] 
            + ( double ) (   - i - m + 1 )     * cx[i-2] ) 
            / ( double ) (     i - m     );
  }

  return;
}
//****************************************************************************80

void legendre_associated_normalized ( int n, int m, double x, double cx[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_NORMALIZED: normalized associated Legendre functions.
//
//  Discussion:
//
//    The unnormalized associated Legendre functions P_N^M(X) have
//    the property that
//
//      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX
//      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
//
//    By dividing the function by the square root of this term,
//    the normalized associated Legendre functions have norm 1.
//
//    However, we plan to use these functions to build spherical
//    harmonics, so we use a slightly different normalization factor of
//
//      sqrt ( ( ( 2 * N + 1 ) * ( N - M )! ) / ( 4 * pi * ( N + M )! ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2005
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
//  Parameters:
//
//    Input, int N, the maximum first index of the Legendre
//    function, which must be at least 0.
//
//    Input, int M, the second index of the Legendre function,
//    which must be at least 0, and no greater than N.
//
//    Input, double X, the point at which the function is to be
//    evaluated.  X must satisfy -1 <= X <= 1.
//
//    Output, double CX[N+1], the values of the first N+1 function.
//
{
  double factor;
  int i;
  int mm;
  double pi = 3.141592653589793;
  double somx2;

  if ( m < 0 )
  {
    cerr << "\n";
    cerr << "LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!\n";
    cerr << "  Input value of M is " << m << "\n";
    cerr << "  but M must be nonnegative.\n";
    exit ( 1 );
  }

  if ( n < m )
  {
    cerr << "\n";
    cerr << "LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!\n";
    cerr << "  Input value of M = " << m << "\n";
    cerr << "  Input value of N = " << n << "\n";
    cerr << "  but M must be less than or equal to N.\n";
    exit ( 1 );
  }

  if ( x < -1.0 )
  {
    cerr << "\n";
    cerr << "LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!\n";
    cerr << "  Input value of X = " << x << "\n";
    cerr << "  but X must be no less than -1.\n";
    exit ( 1 );
  }

  if ( 1.0 < x )
  {
    cerr << "\n";
    cerr << "LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!\n";
    cerr << "  Input value of X = " << x << "\n";
    cerr << "  but X must be no more than 1.\n";
    exit ( 1 );
  }

  for ( i = 0; i <= m-1; i++ )
  {
    cx[i] = 0.0;
  }
  cx[m] = 1.0;

  somx2 = sqrt ( 1.0 - x * x );

  factor = 1.0;
  for ( i = 1; i <= m; i++ )
  {
    cx[m] = -cx[m] * factor * somx2;
    factor = factor + 2.0;
  }

  if ( m+1 <= n )
  {
    cx[m+1] = x * ( double ) ( 2 * m + 1 ) * cx[m];
  }

  for ( i = m+2; i <= n; i++ )
  {
    cx[i] = ( ( double ) ( 2 * i     - 1 ) * x * cx[i-1] 
            + ( double ) (   - i - m + 1 )     * cx[i-2] ) 
            / ( double ) (     i - m     );
  }
//
//  Normalization.
//
  for ( mm = m; mm <= n; mm++ )
  {
    factor = sqrt ( ( ( double ) ( 2 * mm + 1 ) * r8_factorial ( mm - m ) ) 
      / ( 4.0 * pi * r8_factorial ( mm + m ) ) );
    cx[mm] = cx[mm] * factor;
  }

  return;
}
//****************************************************************************80

void legendre_associated_normalized_values ( int *n_data, int *n, int *m, 
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_NORMALIZED_VALUES: normalied associated Legendre.
//
//  Discussion:
//
//    The function considered is the associated Legendre polynomial P^M_N(X).
//
//    In Mathematica, the function can be evaluated by:
//
//      LegendreP [ n, m, x ]
//
//    The function is normalized by dividing by 
//
//      sqrt ( 4 * pi * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
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
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *N, integer M, double X, 
//    the arguments of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = { 
     0.2820947917738781, 
     0.2443012559514600, 
    -0.2992067103010745, 
    -0.07884789131313000, 
    -0.3345232717786446, 
     0.2897056515173922, 
    -0.3265292910163510, 
    -0.06997056236064664, 
     0.3832445536624809, 
    -0.2709948227475519, 
    -0.2446290772414100, 
     0.2560660384200185, 
     0.1881693403754876, 
    -0.4064922341213279, 
     0.2489246395003027, 
     0.08405804426339821, 
     0.3293793022891428, 
    -0.1588847984307093, 
    -0.2808712959945307, 
     0.4127948151484925, 
    -0.2260970318780046 };

  static int m_vec[N_MAX] = { 
    0, 0, 1, 0, 
    1, 2, 0, 1, 
    2, 3, 0, 1, 
    2, 3, 4, 0, 
    1, 2, 3, 4, 
    5 };

  static int n_vec[N_MAX] = { 
    0,  1,  1,  2, 
    2,  2,  3,  3, 
    3,  3,  4,  4, 
    4,  4,  4,  5, 
    5,  5,  5,  5, 
    5 };

  static double x_vec[N_MAX] = { 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50, 
    0.50 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *m = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *m = m_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_associated_values ( int *n_data, int *n, int *m, double *x, 
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_VALUES returns values of associated Legendre functions.
//
//  Discussion:
//
//    The function considered is the associated Legendre function P^M_N(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, int *M, double *X, the arguments of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 19

  int n_vec[N_MAX] = {
     1,  1,  1,  1, 
     1,  2,  2,  2, 
     3,  3,  3,  3, 
     4,  5,  6,  7, 
     8,  9, 10 };
  int m_vec[N_MAX] = {
     0,  0,  0,  0, 
     1,  0,  1,  2, 
     0,  1,  2,  3, 
     2,  2,  3,  3, 
     4,  4,  5 };
  double fx_vec[N_MAX] = {
     0.000000,  0.500000,  0.707107,  1.000000, 
    -0.866025, -0.125000, -1.29904,   2.25000, 
    -0.437500, -0.324759,  5.62500,  -9.74278, 
     4.21875,  -4.92187,   12.7874,   116.685, 
    -1050.67,  -2078.49,   30086.2 };
  double x_vec[N_MAX] = {
    0.0,       0.5,       0.7071067, 1.0, 
    0.5,       0.5,       0.5,       0.5, 
    0.5,       0.5,       0.5,       0.5, 
    0.5,       0.5,       0.5,       0.5, 
    0.5,       0.5,       0.5 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *m = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data];
    *m = m_vec[*n_data];
    *x = x_vec[*n_data];
    *fx = fx_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_function_q ( int n, double x, double cx[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_FUNCTION_Q evaluates the Legendre Q functions.
//
//  Differential equation:
//
//    (1-X*X) Y'' - 2 X Y' + N (N+1) = 0
//
//  First terms:
//
//    Q(0,X) = 0.5 * log((1+X)/(1-X))
//    Q(1,X) = Q(0,X)*X - 1
//    Q(2,X) = Q(0,X)*(3*X*X-1)/4 - 1.5*X
//    Q(3,X) = Q(0,X)*(5*X*X*X-3*X)/4 - 2.5*X^2 + 2/3
//    Q(4,X) = Q(0,X)*(35*X^4-30*X^2+3)/16 - 35/8 * X^3 + 55/24 * X
//    Q(5,X) = Q(0,X)*(63*X^5-70*X^3+15*X)/16 - 63/8*X^4 + 49/8*X^2 - 8/15
//
//  Recursion:
//
//    Q(0) = 0.5 * log ( (1+X) / (1-X) )
//    Q(1) = 0.5 * X * log ( (1+X) / (1-X) ) - 1.0
//
//    Q(N) = ( (2*N-1) * X * Q(N-1) - (N-1) * Q(N-2) ) / N
//
//  Restrictions:
//
//    -1 < X < 1
//
//  Special values:
//
//    Note that the Legendre function Q(N,X) is equal to the
//    associated Legendre function of the second kind,
//    Q(N,M,X) with M = 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2003
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
//  Parameters:
//
//    Input, int N, the highest order function to evaluate.
//
//    Input, double X, the point at which the functions are to be
//    evaluated.  X must satisfy -1 < X < 1.
//
//    Output, double CX[N+1], the values of the first N+1 Legendre
//    functions at the point X.
//
{
  int i;
//
//  Check the value of X.
//
  if ( x <= -1.0 || 1.0 <= x )
  {
    cerr << "\n";
    cerr << "LEGENDRE_FUNCTION_Q - Fatal error!\n";
    cerr << "  Illegal input value of X = " << x << "\n";
    cerr << "  But X must be between -1 and 1.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    return;
  }

  cx[0] = 0.5 * log ( ( 1.0 + x ) / ( 1.0 - x ) );

  if ( n == 0 )
  {
    return;
  }

  cx[1] = x * cx[0] - 1.0;

  for ( i = 2; i <= n; i++ )
  {
    cx[i] = ( ( double ) ( 2 * i - 1 ) * x * cx[i-1] 
            + ( double ) (   - i + 1 ) *     cx[i-2] ) 
            / ( double )       i;
  }

  return;
}
//****************************************************************************80

void legendre_function_q_values ( int *n_data, int *n, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_FUNCTION_Q_VALUES returns values of the Legendre Q function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the order of the function.
//
//    Output, double *X, the point where the function is evaluated.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double fx_vec[N_MAX] = {
     0.00000000, -1.00000000,  0.00000000, 
     0.66666667, -0.40634921,  0.00000000, 
     0.54930614, -0.72534693, -0.81866327, 
    -0.19865477, -0.11616303,  0.29165814 };
  int n_vec[N_MAX] = {
     0,  1,  2, 
     3,  9, 10, 
     0,  1,  2, 
     3,  9, 10 };
  double x_vec[N_MAX] = {
    0.0,  0.0,  0.0, 
    0.0,  0.0,  0.0, 
    0.5,  0.5,  0.5, 
    0.5,  0.5,  0.5  };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data];
    *x = x_vec[*n_data];
    *fx = fx_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_poly ( int n, double x, double cx[], double cpx[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLY evaluates the Legendre polynomials.
//
//  Discussion:
//
//    P(N,1) = 1.
//    P(N,-1) = (-1)^N.
//    | P(N,X) | <= 1 in [-1,1].
//
//    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
//    function of the first kind and order N equals the Legendre polynomial
//    of the first kind and order N.
//
//    The N zeroes of P(N,X) are the abscissas used for Gauss-Legendre
//    quadrature of the integral of a function F(X) with weight function 1
//    over the interval [-1,1].
//
//    The Legendre polynomials are orthogonal under the inner product defined
//    as integration from -1 to 1:
//
//      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
//        = 0 if I =/= J
//        = 2 / ( 2*I+1 ) if I = J.
//
//    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
//
//    A function F(X) defined on [-1,1] may be approximated by the series
//      C0*P(0,X) + C1*P(1,X) + ... + CN*P(N,X)
//    where
//      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,X) dx.
//
//    The formula is:
//
//      P(N,X) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
//
//  Differential equation:
//
//    (1-X*X) * P(N,X)'' - 2 * X * P(N,X)' + N * (N+1) = 0
//
//  First terms:
//
//    P( 0,X) =       1
//    P( 1,X) =       1 X
//    P( 2,X) =  (    3 X^2 -       1)/2
//    P( 3,X) =  (    5 X^3 -     3 X)/2
//    P( 4,X) =  (   35 X^4 -    30 X^2 +     3)/8
//    P( 5,X) =  (   63 X^5 -    70 X^3 +    15 X)/8
//    P( 6,X) =  (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
//    P( 7,X) =  (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
//    P( 8,X) =  ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
//    P( 9,X) =  (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
//    P(10,X) =  (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2
//                 -63 ) /256
//
//  Recursion:
//
//    P(0,X) = 1
//    P(1,X) = X
//    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
//
//    P'(0,X) = 0
//    P'(1,X) = 1
//    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
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
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Input, double X, the point at which the polynomials are to be evaluated.
//
//    Output, double CX[N+1], the values of the Legendre polynomials
//    of order 0 through N at the point X.
//
//    Output, double CPX[N+1], the values of the derivatives of the
//    Legendre polynomials of order 0 through N at the point X.
//
{
  int i;

  if ( n < 0 )
  {
    return;
  }

  cx[0] = 1.0;
  cpx[0] = 0.0;

  if ( n < 1 )
  {
    return;
  }

  cx[1] = x;
  cpx[1] = 1.0;

  for ( i = 2; i <= n; i++ )
  {
    cx[i] = ( ( double ) ( 2 * i - 1 ) * x * cx[i-1] 
            + ( double ) (   - i + 1 )     * cx[i-2] ) 
            / ( double ) (     i     );

    cpx[i] = ( ( double ) ( 2 * i - 1 ) * ( cx[i-1] + x * cpx[i-1] ) 
             + ( double ) (   - i + 1 ) * cpx[i-2] ) 
             / ( double ) (     i     );

  }

  return;
}
//****************************************************************************80

void legendre_poly_coef ( int n, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLY_COEF evaluates the Legendre polynomial coefficients.
//
//  First terms:
//
//     1
//     0     1
//    -1/2   0      3/2
//     0    -3/2    0     5/2
//     3/8   0    -30/8   0     35/8
//     0    15/8    0   -70/8    0     63/8
//    -5/16  0    105/16  0   -315/16   0    231/16
//     0   -35/16   0   315/16   0   -693/16   0    429/16
//
//     1.00000
//     0.00000  1.00000
//    -0.50000  0.00000  1.50000
//     0.00000 -1.50000  0.00000  2.5000
//     0.37500  0.00000 -3.75000  0.00000  4.37500
//     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
//    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
//     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2003
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
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Output, double C[(N+1)*(N+1)], the coefficients of the Legendre polynomials
//    of degree 0 through N.  Each polynomial is stored as a row.
//
{
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n <= 0 )
  {
    return;
  }

  c[1+1*(n+1)] = 1.0;

  for ( i = 2; i <= n; i++ )
  {
    for ( j = 0; j <= i-2; j++ )
    {
      c[i+j*(n+1)] =          
          ( double ) (   - i + 1 ) * c[i-2+j*(n+1)] / ( double ) i;
    }
    for ( j = 1; j <= i; j++ )
    {
      c[i+j*(n+1)] = c[i+j*(n+1)] 
        + ( double ) ( i + i - 1 ) * c[i-1+(j-1)*(n+1)] / ( double ) i;
    }
  }

  return;
}
//****************************************************************************80

void legendre_poly_values ( int *n_data, int *n, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLY_VALUES returns values of the Legendre polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2003
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
//    LC: QA47.A34..
//
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the order of the function.
//
//    Output, double *X, the point where the function is evaluated.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 13

  double fx_vec[N_MAX] = {
     1.0000000,   0.25000000, -0.4062500, 
    -0.3359375,   0.17682442,  0.2212002, 
     0.0000000,  -0.14750000, -0.2800000, 
    -0.3825000,  -0.44000000, -0.4375000, 
     1.0000000 };

  int n_vec[N_MAX] = {
     0,  1,  2, 
     3,  9, 10, 
     3,  3,  3, 
     3,  3,  3, 
     3 };

  double x_vec[N_MAX] = {
    0.25,  0.25,  0.25, 
    0.25,  0.25,  0.25, 
    0.00,  0.10,  0.20, 
    0.30,  0.40,  0.50, 
    1.0 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n = n_vec[*n_data];
    *x = x_vec[*n_data];
    *fx = fx_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int legendre_symbol ( int q, int p )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SYMBOL evaluates the Legendre symbol (Q/P).
//
//  Definition:
//
//    Let P be an odd prime.  Q is a QUADRATIC RESIDUE modulo P
//    if there is an integer R such that R^2 = Q ( mod P ).
//    The Legendre symbol ( Q / P ) is defined to be:
//
//      + 1 if Q ( mod P ) /= 0 and Q is a quadratic residue modulo P,
//      - 1 if Q ( mod P ) /= 0 and Q is not a quadratic residue modulo P,
//        0 if Q ( mod P ) == 0.
//
//    We can also define ( Q / P ) for P = 2 by:
//
//      + 1 if Q ( mod P ) /= 0
//        0 if Q ( mod P ) == 0
//
//  Example:
//
//    (0/7) =   0
//    (1/7) = + 1  ( 1^2 = 1 mod 7 )
//    (2/7) = + 1  ( 3^2 = 2 mod 7 )
//    (3/7) = - 1
//    (4/7) = + 1  ( 2^2 = 4 mod 7 )
//    (5/7) = - 1
//    (6/7) = - 1
//
//  Discussion:
//
//    For any prime P, exactly half of the integers from 1 to P-1
//    are quadratic residues.
//
//    ( 0 / P ) = 0.
//
//    ( Q / P ) = ( mod ( Q, P ) / P ).
//
//    ( Q / P ) = ( Q1 / P ) * ( Q2 / P ) if Q = Q1 * Q2.
//
//    If Q is prime, and P is prime and greater than 2, then:
//
//      if ( Q == 1 ) then
//
//        ( Q / P ) = 1
//
//      else if ( Q == 2 ) then
//
//        ( Q / P ) = + 1 if mod ( P, 8 ) = 1 or mod ( P, 8 ) = 7,
//        ( Q / P ) = - 1 if mod ( P, 8 ) = 3 or mod ( P, 8 ) = 5.
//
//      else
//
//        ( Q / P ) = - ( P / Q ) if Q = 3 ( mod 4 ) and P = 3 ( mod 4 ),
//                  =   ( P / Q ) otherwise.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 March 2001
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Pinter,
//    A Book of Abstract Algebra,
//    McGraw Hill, 1982, pages 236-237.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 86-87.
//
//  Parameters:
//
//    Input, int Q, an integer whose Legendre symbol with
//    respect to P is desired.
//
//    Input, int P, a prime number, greater than 1, with respect
//    to which the Legendre symbol of Q is desired.
//
//    Output, int LEGENDRE_SYMBOL, the Legendre symbol (Q/P).
//    Ordinarily, this will be -1, 0 or 1.
//    L = -2, P is less than or equal to 1.
//    L = -3, P is not prime.
//    L = -4, the internal stack of factors overflowed.
//    L = -5, not enough factorization space.
//
{
# define FACTOR_MAX 20
# define STACK_MAX 50

  int factor[FACTOR_MAX];
  int i;
  int l;
  int nfactor;
  int nleft;
  int nmore;
  int nstack;
  int power[FACTOR_MAX];
  int pstack[STACK_MAX];
  int qstack[STACK_MAX];
//
//  P must be greater than 1.
//
  if ( p <= 1 )
  {
    cerr << "\n";
    cerr << "LEGENDRE_SYMBOL - Fatal error!\n";
    cerr << "  P must be greater than 1.\n";
    exit ( 1 );
  }
//
//  P must be prime.
//
  if ( !( i4_is_prime ( p ) ) )
  {
    cerr << "\n";
    cerr << "LEGENDRE_SYMBOL - Fatal error!\n";
    cerr << "  P is not prime.\n";
    exit ( 1 );
  }
//
//  ( k*P / P ) = 0.
//
  if ( ( q % p ) == 0 )
  {
    return 0;
  }
//
//  For the special case P = 2, (Q/P) = 1 for all odd numbers.
//
  if ( p == 2 )
  {
    return 1;
  }
//
//  Force Q to be nonnegative.
//
  while ( q < 0 )
  {
    q = q + p;
  }

  nstack = 0;
  l = 1;

  for ( ; ; )
  {
    q = q % p;
//
//  Decompose Q into factors of prime powers.
//
    i4_factor ( q, FACTOR_MAX, &nfactor, factor, power, &nleft );

    if ( nleft != 1 )
    {
      cerr << "\n";
      cerr << "LEGENDRE_SYMBOL - Fatal error!\n";
      cerr << "  Not enough factorization space.\n";
      exit ( 1 );
    }
//
//  Each factor which is an odd power is added to the stack.
//
    nmore = 0;

    for ( i = 0; i < nfactor; i++ )
    {
      if ( ( power[i] % 2 ) == 1 )
      {
        nmore = nmore + 1;

        if ( STACK_MAX <= nstack )
        {
          cerr << "\n";
          cerr << "LEGENDRE_SYMBOL - Fatal error!\n";
          cerr << "  Stack overflow!\n";
          exit ( 1 );
        }

        pstack[nstack] = p;
        qstack[nstack] = factor[i];
        nstack = nstack + 1;
      }

    }

    if ( nmore != 0 )
    {
      nstack = nstack - 1;
      q = qstack[nstack];
//
//  Check for a Q of 1 or 2.
//
      if ( q == 1 )
      {
        l = + 1 * l;
      }
      else if ( q == 2 && ( ( p % 8 ) == 1 || ( p % 8 ) == 7 ) )
      {
        l = + 1 * l;
      }
      else if ( q == 2 && ( ( p % 8 ) == 3 || ( p % 8 ) == 5 ) )
      {
        l = - 1 * l;
      }
      else
      {
        if ( ( p % 4 ) == 3 && ( q % 4 ) == 3 )
        {
          l = - 1 * l;
        }

        i4_swap ( &p, &q );

        continue;

      }

    }
//
//  If the stack is empty, we're done.
//
    if ( nstack == 0 )
    {
      break;
    }
//
//  Otherwise, get the last P and Q from the stack, and process them.
//
    nstack = nstack - 1;
    p = pstack[nstack];
    q = qstack[nstack];
  }

  return l;

# undef FACTOR_MAX
# undef STACK_MAX
}
//****************************************************************************80

double lerch ( double z, int s, double a )

//****************************************************************************80
//
//  Purpose:
//
//    LERCH estimates the Lerch transcendent function.
//
//  Definition:
//
//    The Lerch transcendent function is defined as:
//
//      LERCH ( Z, S, A ) = Sum ( 0 <= K < +oo ) Z**K / ( A + K )^S
//
//    excluding any term with ( A + K ) = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      LerchPhi[z,s,a]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Thanks:
//
//    Oscar van Vlijmen
//
//  Parameters:
//
//    Input, double Z, int S, double A, the parameters of the function.
//
//    Output, double LERCH, an approximation to the Lerch
//    transcendent function.
//
{
  double eps = 1.0E-10;
  int k;
  double total;
  double term;
  double z_k;

  total = 0.0;

  if ( z <= 0.0 )
  {
    return total;
  }

  k = 0;
  z_k = 1.0;

  for ( ; ; )
  {
    if ( a + ( double ) ( k ) != 0.0 )
    {
      term = z_k / pow ( a + ( double ) k, s );

      total = total + term;

      if ( fabs ( term ) <= eps * ( 1.0 + fabs ( total ) ) )
      {
        break;
      }
    }

    k = k + 1;
    z_k = z_k * z;

  }

  return total;
}
//****************************************************************************80

void lerch_values ( int *n_data, double *z, int *s, double *a, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LERCH_VALUES returns some values of the Lerch transcendent function.
//
//  Discussion:
//
//    The Lerch function is defined as
//
//      Phi(z,s,a) = Sum ( 0 <= k < +oo ) z^k / ( a + k )^s
//
//    omitting any terms with ( a + k ) = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      LerchPhi[z,s,a]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *Z, the parameters of the function.
//
//    Output, int *S, the parameters of the function.
//
//    Output, double *A, the parameters of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double a_vec[N_MAX] = {
    0.0, 
    0.0, 
    0.0, 
    1.0, 
    1.0, 
    1.0, 
    2.0, 
    2.0, 
    2.0, 
    3.0, 
    3.0, 
    3.0 };

  double fx_vec[N_MAX] = {
    0.1644934066848226E+01, 
    0.1202056903159594E+01, 
    0.1000994575127818E+01, 
    0.1164481052930025E+01, 
    0.1074426387216080E+01, 
    0.1000492641212014E+01, 
    0.2959190697935714, 
    0.1394507503935608, 
    0.9823175058446061E-03, 
    0.1177910993911311, 
    0.3868447922298962E-01, 
    0.1703149614186634E-04 };

  int s_vec[N_MAX] = {
     2, 3, 10, 
     2, 3, 10, 
     2, 3, 10, 
     2, 3, 10 };

  double z_vec[N_MAX] = {
    0.1000000000000000E+01, 
    0.1000000000000000E+01, 
    0.1000000000000000E+01, 
    0.5000000000000000, 
    0.5000000000000000, 
    0.5000000000000000, 
    0.3333333333333333, 
    0.3333333333333333, 
    0.3333333333333333, 
    0.1000000000000000, 
    0.1000000000000000, 
    0.1000000000000000 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }
 
  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *z = 0.0;
    *s = 0;
    *a = 0.0;
    *fx = 0.0;
  }
  else
  {
    *z = z_vec[*n_data-1];
    *s = s_vec[*n_data-1];
    *a = a_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void lock ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    LOCK returns the number of codes for a lock with N buttons.
//
//  Discussion:
//
//    A button lock has N numbered buttons.  To open the lock, groups
//    of buttons must be pressed in the correct order.  Each button
//    may be pushed no more than once.  Thus, a code for the lock is
//    an ordered list of the groups of buttons to be pushed.
//
//    For this discussion, we will assume that EVERY button is pushed
//    at some time, as part of the code.  To count the total number
//    of codes, including those which don't use all the buttons, then
//    the number is 2 * A(N), or 2 * A(N) - 1 if we don't consider the
//    empty code to be valid.
//
//  Examples:
//
//    If there are 3 buttons, then there are 13 possible "full button" codes:
//
//      (123)
//      (12) (3)
//      (13) (2)
//      (23) (1)
//      (1) (23)
//      (2) (13)
//      (3) (12)
//      (1) (2) (3)
//      (1) (3) (2)
//      (2) (1) (3)
//      (2) (3) (1)
//      (3) (1) (2)
//      (3) (2) (1)
//
//
//    and, if we don't need to push all the buttons, every "full button" code above
//    yields a distinct "partial button" code by dropping the last set of buttons:
//
//      ()
//      (12)
//      (13)
//      (23)
//      (1)
//      (2)
//      (3)
//      (1) (2)
//      (1) (3)
//      (2) (1)
//      (2) (3)
//      (3) (1)
//      (3) (2)
//
//  First values:
//
//     N         A(N)
//     0           1
//     1           1
//     2           3
//     3          13
//     4          75
//     5         541
//     6        4683
//     7       47293
//     8      545835
//     9     7087261
//    10   102247563
//
//  Recursion:
//
//    A(I) = sum ( 0 <= J < I ) Binomial ( I, N-J ) * A(J)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Velleman, Gregory Call,
//    Permutations and Combination Locks,
//    Mathematics Magazine,
//    Volume 68, Number 4, October 1995, pages 243-253.
//
//  Parameters:
//
//    Input, int N, the maximum number of lock buttons.
//
//    Output, int A[N+1], the number of lock codes with 0 to N buttons.
//
{
  int combo;
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  a[0] = 1;

  for ( i = 1; i <= n; i++ )
  {
    a[i] = 0;
    for ( j = 0; j <= i-1; j++ )
    {
      a[i] = a[i] + i4_choose ( i, i-j ) * a[j];
    }
  }

  return;
}
//****************************************************************************80

void meixner ( int n, double beta, double c, double x, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    MEIXNER evaluates Meixner polynomials at a point.
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
//  Reference:
//
//    Walter Gautschi,
//    Orthogonal Polynomials: Computation and Approximation,
//    Oxford, 2004,
//    ISBN: 0-19-850672-4,
//    LC: QA404.5 G3555.
//
//  Parameters:
//
//    Input, int N, the maximum order of the polynomial.  
//    N must be at least 0.
//
//    Input, double BETA, the Beta parameter.  0 < BETA.
//
//    Input, double C, the C parameter.  0 < C < 1.
//
//    Input, double X, the evaluation point.
//
//    Output, double V[N+1], the value of the polynomials at X.
//
{
  int i;

  if ( beta <= 0.0 )
  {
    cerr << "\n";
    cerr << "MEIXNER - Fatal error!\n";
    cerr << "  Parameter BETA must be positive.\n";
    exit ( 1 );
  }

  if ( c <= 0.0 || 1.0 <= c )
  {
    cerr << "\n";
    cerr << "MEIXNER - Fatal error!\n";
    cerr << "  Parameter C must be strictly between 0 and 1.\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "MEIXNER - Fatal error!\n";
    cerr << "  Parameter N must be nonnegative.\n";
    exit ( 1 );
  }

  v[0] = 1.0;

  if ( n == 0 )
  {
    return;
  }

  v[1] = ( c - 1.0 ) * x / beta / c + 1.0;

  if ( n == 1 )
  {
    return;
  }

  for ( i = 1; i < n; i++ )
  {
    v[i+1] = ( 
      ( ( c - 1.0 ) * x + ( 1.0 + c ) 
      * ( double ) ( i ) + beta * c ) * v[i]
      - ( double ) ( i ) * v[i-1] 
      ) / ( ( double ) ( i ) + beta );
  }

  return;
}
//****************************************************************************80

int mertens ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    MERTENS evaluates the Mertens function.
//
//  Discussion:
//
//    The Mertens function M(N) is the sum from 1 to N of the Moebius
//    function MU.  That is,
//
//    M(N) = sum ( 1 <= I <= N ) MU(I)
//
//        N   M(N)
//        --  ----
//         1     1
//         2     0
//         3    -1
//         4    -1
//         5    -2
//         6    -1
//         7    -2
//         8    -2
//         9    -2
//        10    -1
//        11    -2
//        12    -2
//       100     1
//      1000     2
//     10000   -23
//    100000   -48
//
//    The determinant of the Redheffer matrix of order N is equal
//    to the Mertens function M(N).
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    M Deleglise, J Rivat,
//    Computing the Summation of the Moebius Function,
//    Experimental Mathematics,
//    Volume 5, 1996, pages 291-295.
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45
//
//  Parameters:
//
//    Input, int N, the argument.
//
//    Output, int MERTENS, the value.
//
{
  int i;
  int value;

  value = 0;

  for ( i = 1; i <= n; i++ )
  {
    value = value + moebius ( i );
  }
  return value;
}
//****************************************************************************80

void mertens_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    MERTENS_VALUES returns some values of the Mertens function.
//
//  Discussion:
//
//    The Mertens function M(N) is the sum from 1 to N of the Moebius
//    function MU.
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
//  Reference:
//
//    M Deleglise, J Rivat,
//    Computing the Summation of the Moebius Function,
//    Experimental Mathematics,
//    Volume 5, 1996, pages 291-295.
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When 
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the argument of the Mertens function.
//
//    Output, int *C, the value of the Mertens function.
//
{
# define N_MAX 15

  int c_vec[N_MAX] = {
      1,   0,  -1,   -1,  -2,  -1,  -2,  -2,   -2,  -1, 
     -2,  -2,   1,    2, -23 };
  int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, 
     11,  12,  100, 1000, 10000 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int moebius ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    MOEBIUS returns the value of MU(N), the Moebius function of N.
//
//  Discussion:
//
//    MU(N) is defined as follows:
//
//      MU(N) = 1 if N = 1;
//              0 if N is divisible by the square of a prime;
//              (-1)^K, if N is the product of K distinct primes.
//
//    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
//    if N is a square, cube, etc.
//
//    The Moebius function is related to Euler's totient function:
//
//      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
//
//  First values:
//
//     N  MU(N)
//
//     1    1
//     2   -1
//     3   -1
//     4    0
//     5   -1
//     6    1
//     7   -1
//     8    0
//     9    0
//    10    1
//    11   -1
//    12    0
//    13   -1
//    14    1
//    15    1
//    16    0
//    17   -1
//    18    0
//    19   -1
//    20    0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 March 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the value to be analyzed.
//
//    Output, int MOEBIUS, the value of MU(N).
//    If N is less than or equal to 0, or there was not enough internal 
//    space for factoring, MOEBIUS is returned as -1.
//
{
# define FACTOR_MAX 20

  int factor[FACTOR_MAX];
  int i;
  int nfactor;
  int nleft;
  int power[FACTOR_MAX];
  int value;

  if ( n <= 0 )
  {
    return (-1);
  }

  if ( n == 1 )
  {
    return 1;
  }
//
//  Factor N.
//
  i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );

  if ( nleft != 1 )
  {
    cerr << "\n";
    cerr << "MOEBIUS - Fatal error!\n";
    cerr <<  "  Not enough factorization space.\n";
    exit ( 1 );
  }

  value = 1;

  for ( i = 0; i < nfactor; i++ )
  {
    value = -value;

    if ( 1 < power[i] )
    {
      return 0;
    }

  }

  return value;
# undef FACTOR_MAX
}
//****************************************************************************80

void moebius_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    MOEBIUS_VALUES returns some values of the Moebius function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the argument of the Moebius function.
//
//    Output, int *C, the value of the Moebius function.
//
{
# define N_MAX 20

  int c_vec[N_MAX] = {
      1,  -1,  -1,   0,  -1,   1,  -1,   0,   0,   1, 
     -1,   0,  -1,   1,   1,   0,  -1,   0,  -1,   0 };

  int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, 
     11,  12,  13,  14,  15,  16,  17,  18,  19,  20 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void motzkin ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    MOTZKIN returns the Motzkin numbers up to order N.
//
//  Discussion:
//
//    The Motzkin number A(N) counts the number of distinct paths
//    from (0,0) to (0,N) in which the only steps used are
//    (1,1), (1,-1) and (1,0), and the path is never allowed to
//    go below the X axis.
//
//  First values:
//
//     N  A(N)
//
//     0    1
//     1    1
//     2    2
//     3    4
//     4    9
//     5   21
//     6   51
//     7  127
//     8  323
//     9  835
//    10 2188
//
//  Recursion:
//
//    A(N) = A(N-1) + sum ( 0 <= K <= N-2 ) A(K) * A(N-2-K)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input, int N, the highest order Motzkin number to compute.
//
//    Output, int A[N+1], the Motzkin numbers.
//
{
  int i;
  int j;

  if ( n < 0 )
  {
    return;
  }

  a[0] = 1;

  for ( i = 1; i <= n; i++ )
  {
    a[i] = a[i-1];
    for ( j = 0; j <= i-2; j++ )
    {
      a[i] = a[i] + a[j] * a[i-2-j];
    }
  }

  return;
}
//****************************************************************************80

double normal_01_cdf_inv ( double p )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF_INV inverts the standard normal CDF.
//
//  Discussion:
//
//    The result is accurate to about 1 part in 10**16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2004
//
//  Author:
//
//    Original FORTRAN77 version by Michael Wichura.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Michael Wichura,
//    The Percentage Points of the Normal Distribution,
//    Algorithm AS 241,
//    Applied Statistics,
//    Volume 37, Number 3, pages 477-484, 1988.
//
//  Parameters:
//
//    Input, double P, the value of the cumulative probability 
//    densitity function.  0 < P < 1.  If P is outside this range, an
//    "infinite" value is returned.
//
//    Output, double NORMAL_01_CDF_INVERSE, the normal deviate value 
//    with the property that the probability of a standard normal deviate being 
//    less than or equal to this value is P.
//
{
  double a[8] = { 
    3.3871328727963666080,     1.3314166789178437745e+2,
    1.9715909503065514427e+3,  1.3731693765509461125e+4,
    4.5921953931549871457e+4,  6.7265770927008700853e+4,
    3.3430575583588128105e+4,  2.5090809287301226727e+3 };
  double b[8] = {
    1.0,                       4.2313330701600911252e+1,
    6.8718700749205790830e+2,  5.3941960214247511077e+3,
    2.1213794301586595867e+4,  3.9307895800092710610e+4,
    2.8729085735721942674e+4,  5.2264952788528545610e+3 };
  double c[8] = {
    1.42343711074968357734,     4.63033784615654529590,
    5.76949722146069140550,     3.64784832476320460504,
    1.27045825245236838258,     2.41780725177450611770e-1,
    2.27238449892691845833e-2,  7.74545014278341407640e-4 };
  double const1 = 0.180625;
  double const2 = 1.6;
  double d[8] = {
    1.0,                        2.05319162663775882187,
    1.67638483018380384940,     6.89767334985100004550e-1,
    1.48103976427480074590e-1,  1.51986665636164571966e-2,
    5.47593808499534494600e-4,  1.05075007164441684324e-9 };
  double e[8] = {
    6.65790464350110377720,     5.46378491116411436990,
    1.78482653991729133580,     2.96560571828504891230e-1,
    2.65321895265761230930e-2,  1.24266094738807843860e-3,
    2.71155556874348757815e-5,  2.01033439929228813265e-7 };
  double f[8] = {
    1.0,                        5.99832206555887937690e-1,
    1.36929880922735805310e-1,  1.48753612908506148525e-2,
    7.86869131145613259100e-4,  1.84631831751005468180e-5,
    1.42151175831644588870e-7,  2.04426310338993978564e-15 };
  double q;
  double r;
  double split1 = 0.425;
  double split2 = 5.0;
  double value;

  if ( p <= 0.0 )
  {
    value = -r8_huge ( );
    return value;
  }

  if ( 1.0 <= p )
  {
    value = r8_huge ( );
    return value;
  }

  q = p - 0.5;

  if ( r8_abs ( q ) <= split1 )
  {
    r = const1 - q * q;
    value = q * r8poly_value ( 8, a, r ) / r8poly_value ( 8, b, r );
  }
  else
  {
    if ( q < 0.0 )
    {
      r = p;
    }
    else
    {
      r = 1.0 - p;
    }

    if ( r <= 0.0 )
    {
      value = r8_huge ( );
    }
    else
    {
      r = sqrt ( - log ( r ) );

      if ( r <= split2 )
      {
        r = r - const2;
        value = r8poly_value ( 8, c, r ) / r8poly_value ( 8, d, r ); 
       }
       else
       {
         r = r - split2;
         value = r8poly_value ( 8, e, r ) / r8poly_value ( 8, f, r );
      }
    }

    if ( q < 0.0 )
    {
      value = - value;
    }

  }

  return value;
}
//****************************************************************************80

int omega ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    OMEGA returns OMEGA(N), the number of distinct prime divisors of N.
//
//  Discussion:
//
//    If N = 1, then
//
//      OMEGA(N) = 1
//
//    else if the prime factorization of N is
//
//      N = P1**E1 * P2^E2 * ... * PM^EM,
//
//    then
//
//      OMEGA(N) = M
//
//  First values:
//
//     N   OMEGA(N)
//
//     1    1
//     2    1
//     3    1
//     4    1
//     5    1
//     6    2
//     7    1
//     8    1
//     9    1
//    10    2
//    11    1
//    12    2
//    13    1
//    14    2
//    15    2
//    16    1
//    17    1
//    18    2
//    19    1
//    20    2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 November 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the value to be analyzed.  N must be 1 or
//    greater.
//
//    Output, int OMEGA, the value of OMEGA(N).  But if N is 0 or
//    less, OMEGA is returned as 0, a nonsense value.  If there is
//    not enough room for factoring, OMEGA is returned as -1.
//
{
# define FACTOR_MAX 20

  int factor[FACTOR_MAX];
  int nfactor;
  int nleft;
  int power[FACTOR_MAX];

  if ( n <= 0 )
  {
    return 0;
  }

  if ( n == 1 )
  {
    return 1;
  }
//
//  Factor N.
//
  i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );

  if ( nleft != 1 )
  {
    cerr << "\n";
    cerr << "OMEGA - Fatal error!\n";
    cerr << "  Not enough factorization space.\n";
    exit ( 1 );
  }

  return nfactor;

# undef FACTOR_MAX
}
//****************************************************************************80

void omega_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    OMEGA_VALUES returns some values of the OMEGA function.
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
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the argument of the OMEGA function.
//
//    Output, int *C, the value of the OMEGA function.
//
{
# define N_MAX 23

  int c_vec[N_MAX] = {
      1,   1,   1,   1,   1, 
      2,   1,   1,   1,   2, 
      3,   1,   4,   4,   3, 
      1,   5,   2,   2,   1, 
      6,   7,   8 };
  int n_vec[N_MAX] = {
           1,        2,        3,        4,        5, 
           6,        7,        8,        9,       10, 
          30,      101,      210,     1320,     1764, 
        2003,     2310,     2827,     8717,    12553, 
       30030,   510510,  9699690 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void partition_count_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    PARTITION_COUNT_VALUES returns values of the integer partition count.
//
//  Discussion:
//
//    A partition of an integer N is a representation of the integer
//    as the sum of nonzero positive integers.  The order of the summands
//    does not matter.  The number of partitions of N is symbolized
//    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
//    following partitions:
//
//    5 = 5
//      = 4 + 1
//      = 3 + 2
//      = 3 + 1 + 1
//      = 2 + 2 + 1
//      = 2 + 1 + 1 + 1
//      = 1 + 1 + 1 + 1 + 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the integer.
//
//    Output, int *C, the number of partitions of the integer.
//
{
# define N_MAX 21

  int c_vec[N_MAX] = {
      1, 
      1,   2,   3,   5,   7,  11,  15,  22,  30,  42, 
     56,  77, 101, 135, 176, 231, 297, 385, 490, 627 };
  int n_vec[N_MAX] = {
     0,  
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void partition_distinct_count_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    PARTITION_DISTINCT_COUNT_VALUES returns some values of Q(N).
//
//  Discussion:
//
//    A partition of an integer N is a representation of the integer
//    as the sum of nonzero positive integers.  The order of the summands
//    does not matter.  The number of partitions of N is symbolized
//    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
//    following partitions:
//
//    5 = 5
//      = 4 + 1
//      = 3 + 2
//      = 3 + 1 + 1
//      = 2 + 2 + 1
//      = 2 + 1 + 1 + 1
//      = 1 + 1 + 1 + 1 + 1
//
//    However, if we require that each member of the partition
//    be distinct, so that no nonzero summand occurs more than once,
//    we are computing something symbolized by Q(N).
//    The number 5 has Q(N) = 3, because it has the following partitions
//    into distinct parts:
//
//    5 = 5
//      = 4 + 1
//      = 3 + 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the integer.
//
//    Output, int *C, the number of partitions of the integer
//    into distinct parts.
//
{
# define N_MAX 21

  int c_vec[N_MAX] = {
      1, 
      1,   1,   2,   2,   3,   4,   5,   6,   8,  10, 
     12,  15,  18,  22,  27,  32,  38,  46,  54,  64 };
  int n_vec[N_MAX] = {
     0,  
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int pentagon_num ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PENTAGON_NUM computes the N-th pentagonal number.
//
//  Discussion:
//
//    The pentagonal number P(N) counts the number of dots in a figure of
//    N nested pentagons.  The pentagonal numbers are defined for both
//    positive and negative N.
//
//    The formula is:
//
//      P(N) = ( N * ( 3 * N - 1 ) ) / 2
//
//  First values:
//
//    N   P
//
//   -5   40
//   -4   26
//   -3   15
//   -2    7
//   -1    2
//    0    0
//    1    1
//    2    5
//    3   12
//    4   22
//    5   35
//    6   51
//    7   70
//    8   92
//    9  117
//   10  145
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the index of the pentagonal number desired.
//
//    Output, int PENTAGON_NUM, the value of the N-th pentagonal number.
//
{
  int p;

  p = ( n * ( 3 * n - 1 ) ) / 2;

  return p;
}
//****************************************************************************80

int phi ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PHI computes the number of relatively prime predecessors of an integer.
//
//  Discussion:
//
//    PHI(N) is the number of integers between 1 and N which are
//    relatively prime to N.  I and J are relatively prime if they
//    have no common factors.  The function PHI(N) is known as
//    "Euler's totient function".
//
//    By convention, 1 and N are relatively prime.
//
//    The formula is:
//
//      PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
// 
//      PHI(P**K) = P^(K-1) * ( P - 1 ) if P is prime.
//
//      PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
//
//      N = Sum ( D divides N ) PHI(D).
//
//  First values:
//
//     N  PHI(N)
//
//     1    1
//     2    1
//     3    2
//     4    2
//     5    4
//     6    2
//     7    6
//     8    4
//     9    6
//    10    4
//    11   10
//    12    4
//    13   12
//    14    6
//    15    8
//    16    8
//    17   16
//    18    6
//    19   18
//    20    8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the value to be analyzed.
//
//    Output, int PHI, the value of PHI(N).  If N is less than
//    or equal to 0, PHI will be returned as 0.  If there is not enough
//    room for full factoring of N, PHI will be returned as -1.
//
{
# define FACTOR_MAX 20

  int factor[FACTOR_MAX];
  int i;
  int nfactor;
  int nleft;
  int power[FACTOR_MAX];
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  if ( n == 1 )
  {
    return 1;
  }
//
//  Factor N.
//
  i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );

  if ( nleft != 1 )
  {
    cerr << "\n";
    cerr << "PHI - Fatal error!\n";
    cerr << "  Not enough factorization space.\n";
    exit ( 1 );
  }

  value = 1;
  for ( i = 0; i < nfactor; i++ )
  {
    value = value * ( int ) pow ( ( double ) factor[i], power[i]-1 ) 
      * ( factor[i] - 1 );
  }

  return value;
# undef FACTOR_MAX
}
//****************************************************************************80

void phi_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    PHI_VALUES returns some values of the PHI function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the argument of the PHI function.
//
//    Output, int *C, the value of the PHI function.
//
{
# define N_MAX 20

  int c_vec[N_MAX] = {
      1,   1,   2,   2,   4,   2,   6,   4,   6,   4, 
      8,   8,  16,  20,  16,  40, 148, 200, 200, 648 };
  int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, 
     20,  30,  40,  50,  60, 100, 149, 500, 750, 999 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int poly_bernoulli ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    POLY_BERNOULLI evaluates the poly-Bernolli numbers with negative index.
//
//  Discussion:
//
//    The poly-Bernoulli numbers B_n^k were defined by M Kaneko
//    formally as the coefficients of X^n/n! in a particular power
//    series.  He also showed that, when the super-index is negative,
//    we have
//
//      B_n^(-k) = Sum ( 0 <= j <= min ( n, k ) )
//        (j!)^2 * S(n+1,j+1) * S(k+1,j+1)
//
//    where S(n,k) is the Stirling number of the second kind, the number of
//    ways to partition a set of size n into k nonempty subset.
//
//    B_n^(-k) is also the number of "lonesum matrices", that is, 0-1
//    matrices of n rows and k columns which are uniquely reconstructable
//    from their row and column sums.
//
//    The poly-Bernoulli numbers get large very quickly.
//
//  Table:
//
//    \ K 0  1    2     3      4       5        6
//    N
//    0   1  1    1     1      1       1        1
//    1   1  2    4     8     16      32       64
//    2   1  4   14    46    146     454     1394
//    3   1  8   46   230   1066    4718    20266
//    4   1 16  146  1066   6902   41506   237686
//    5   1 32  454  4718  41506  329462  2441314
//    6   1 64 1394 20266 237686 2441314 22934774
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
//  Reference:
//
//    Chad Brewbaker,
//    Lonesum (0,1) Matrices and Poly-Bernoulli Numbers of Negative Index,
//    MS Thesis,
//    Iowa State University, 2005.
//
//    M Kaneko,
//    Poly-Bernoulli Numbers,
//    Journal Theorie des Nombres Bordeaux,
//    Volume 9, 1997, pages 221-228.
//
//  Parameters:
//
//    Input, int N, K, the indices.  N and K should be nonnegative.
//
//    Output, int POLY_BERNOULLI, the value of B_N^(-K).
//
{
  int b;
  int j;
  int jfact;
  int jhi;
  int m;
  int *s;

  if ( n < 0 )
  {
    b = 0;
    return b;
  }
  else if ( n == 0 )
  {
    b = 1;
    return b;
  }

  if ( k <= 0 )
  {
    b = 0;
    return b;
  }
  else if ( k == 0 )
  {
    b = 1;
    return b;
  }

  jhi = i4_min ( n, k );
  m = i4_max ( n, k ) + 1;

  s = stirling2 ( m, m );

  jfact = 1;
  b = 0;

  for ( j = 0; j <= jhi; j++ )
  {
    b = b + jfact * jfact * s[n+j*m] * s[k+j*m];

    jfact = jfact * ( j + 1 );
  }

  delete [] s;

  return b;
}
//****************************************************************************80

int poly_coef_count ( int dim, int degree )

//****************************************************************************80
//
//  Purpose:
//
//    POLY_COEF_COUNT: polynomial coefficient count given dimension and degree.
//
//  Discussion:
//
//    To count all monomials of degree 5 or less in dimension 3,
//    we can count all monomials of degree 5 in dimension 4.
//
//    To count all monomials of degree 5 in dimension 4, we imagine
//    that each of the variables X, Y, Z and W is a "box" and that
//    we need to drop 5 pebbles into these boxes.  Every distinct
//    way of doing this represents a degree 5 monomial in dimension 4.
//    Ignoring W gives us monomials up to degree five in dimension 3.
//
//    To count them, we draw 3 lines as separators to indicate the
//    4 boxes, and then imagine all disctinct sequences involving
//    the three lines and the 5 pebbles.  Indicate the lines by 1's
//    and the pebbles by 0's and we're asking for the number of
//    permutations of 3 1's and 5 0's, which is 8! / (3! 5!)
//
//    In other words, 56 = 8! / (3! 5!) is:
//    * the number of monomials of degree exactly 5 in dimension 4, 
//    * the number of monomials of degree 5 or less in dimension 3, 
//    * the number of polynomial coefficients of a polynomial of 
//      degree 5 in (X,Y,Z).
//
//    In general, the formula for the number of monomials of degree DEG
//    or less in dimension DIM is
//
//      (DEG+DIM)! / (DEG! * DIM!)
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
//  Parameters:
//
//    Input, int DIM, the dimension of the polynomial.
//    0 <= DIM.
//
//    Input, int DEGREE, the degree of the polynomnial
//    0 <= DEGREE
//
//    Output, int POLY_COEF_COUNT, the number of coefficients 
//    in the general polynomial of dimension DIM and degree DEGREE.
//
{
  int value;

  if ( dim < 0 )
  {
    value = -1;
  }
  else if ( degree < 0 )
  {
    value = -1;
  }
  else
  {
    value = i4_choose ( degree + dim, degree );
  }

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

void psi_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    PSI_VALUES returns some values of the Psi or Digamma function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      PolyGamma[x]
//
//    or
//
//      Polygamma[0,x]
//
//    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
//
//    PSI(1) = -Euler's constant.
//
//    PSI(X+1) = PSI(X) + 1 / X.
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
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 11

  double fx_vec[N_MAX] = { 
     -0.5772156649015329E+00,  
     -0.4237549404110768E+00,  
     -0.2890398965921883E+00,  
     -0.1691908888667997E+00,  
     -0.6138454458511615E-01,  
      0.3648997397857652E-01,  
      0.1260474527734763E+00,  
      0.2085478748734940E+00,  
      0.2849914332938615E+00,  
      0.3561841611640597E+00,  
      0.4227843350984671E+00 };

  double x_vec[N_MAX] = { 
     1.0E+00,  
     1.1E+00,  
     1.2E+00,  
     1.3E+00,  
     1.4E+00,  
     1.5E+00,  
     1.6E+00,  
     1.7E+00,  
     1.8E+00,  
     1.9E+00,  
     2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int pyramid_num ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_NUM returns the N-th pyramidal number.
//
//  Discussion:
//
//    The N-th pyramidal number P(N) is formed by the sum of the first
//    N triangular numbers T(J):
//
//      T(J) = sum ( 1 <= J <= N ) J
//
//      P(N) = sum ( 1 <= I <= N ) T(I)
//
//    By convention, T(0) = 0.
//
//    The formula is:
//
//      P(N) = ( (N+1)^3 - (N+1) ) / 6
//
//  First Values:
//
//      0
//      1
//      4
//     10
//     20
//     35
//     56
//     84
//    120
//    165
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the index of the desired number, which must be
//    at least 0.
//
//    Output, int PYRAMID_NUM, the N-th pyramidal number.
//
{
  int value;

  value = ( ( n + 1 ) * ( n + 1 ) * ( n + 1 ) - ( n + 1 ) ) / 6;

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
//    02 April 2005
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
    value = -x;
  }
  return value;
}
//****************************************************************************80

double r8_acosh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ACOSH returns the inverse hyperbolic cosine of a number.
//
//  Discussion:
//
//    Applying the inverse function
//
//      Y = R8_ACOSH(X)
//
//    implies that
//
//      X = COSH(Y) = 0.5 * ( EXP(Y) + EXP(-Y) ).
//
//    For every X greater than or equal to 1, there are two possible
//    choices Y such that X = COSH(Y), differing only in sign.  It
//    is usual to resolve this choice by taking the value of ACOSH(X)
//    to be nonnegative.
//
//  Method:
//
//    One formula is:
//
//      R8_ACOSH = LOG ( X + SQRT ( X^2 - 1.0 ) )
//
//    but this formula suffers from roundoff and overflow problems.
//    The formula used here was recommended by W Kahan, as discussed
//    by Moler.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Cleve Moler,
//    Trigonometry is a Complex Subject,
//    MATLAB News and Notes,
//    Summer 1998.
//
//  Parameters:
//
//    Input, double X, the number whose inverse hyperbolic cosine is desired.
//    X should be greater than or equal to 1.
//
//    Output, double R8_ACOSH, the inverse hyperbolic cosine of X.  The
//    principal value (that is, the positive value of the two ) is returned.
//
{
  double value;

  if ( x < 1.0 )
  {
    cerr << "\n";
    cerr << "R8_ACOSH - Fatal error!\n";
    cerr << "  Argument X must satisfy 1 <= X.\n";
    cerr << "  The input X = " << x << "\n";
    exit ( 1 );
  }

  value = 2.0 * log ( 
    sqrt ( 0.5 * ( x + 1.0 ) ) + sqrt ( 0.5 * ( x - 1.0 ) ) );

  return value;
}
//****************************************************************************80

double r8_asinh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ASINH returns the inverse hyperbolic sine of a number.
//
//  Discussion:
//
//    The assertion that:
//
//      Y = R8_ASINH ( X )
//
//    implies that
//
//      X = SINH(Y) = 0.5 * ( EXP(Y) - EXP(-Y) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose inverse hyperbolic 
//    sine is desired.
//
//    Output, double R8_ASINH, the inverse hyperbolic sine of X.
//
{
  double value;

  value = log ( x + sqrt ( x * x + 1.0 ) );

  return value;
}
//****************************************************************************80

double r8_atanh ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ATANH returns the inverse hyperbolic tangent of a number.
//
//  Discussion:
//
//    Y = R8_ATANH ( X )
//
//    implies that
//
//    X = TANH(Y) = ( EXP(Y) - EXP(-Y) ) / ( EXP(Y) + EXP(-Y) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose inverse hyperbolic 
//    tangent is desired.  The absolute value of X should be less than 
//    or equal to 1.
//
//    Output, double R8_ATANH, the inverse hyperbolic tangent of X.
//
{
  double value;

  if ( 1.0 <= r8_abs ( x ) )
  {
    cerr << "\n";
    cerr << "R8_ATANH - Fatal error!\n";
    cerr << "  ABS(X) must be < 1.\n";
    cerr << "  Your input is X = " << x << "\n";
    exit ( 1 );
  }

  value = 0.5 * log ( ( 1.0 + x ) / ( 1.0 - x ) );

  return value;
}
//****************************************************************************80

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the combinatorial coefficient C(N,K).
//
//    Real arithmetic is used, and C(N,K) is computed directly, via
//    Gamma functions, rather than recursively.
//
//    C(N,K) is the number of distinct combinations of K objects
//    chosen from a set of N distinct objects.  A combination is
//    like a set, in that order does not matter.
//
//    The number of combinations of 2 things chosen from 5 is 10.
//
//    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
//
//    The actual combinations may be represented as:
//
//      (1,2), (1,3), (1,4), (1,5), (2,3),
//      (2,4), (2,5), (3,4), (3,5), (4,5).
//
//    The formula for C(N,K) may be written:
//
//      C(N,K) = N! / ( (N-K)! * K! )
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
//  Parameters:
//
//    Input, int N, the value of N.
//
//    Input, int K, the value of K.
//
//    Output, double R8_CHOOSE, the value of C(N,K)
//
{
  double arg;
  double fack;
  double facn;
  double facnmk;
  double value;

  if ( n < 0 )
  {
    value = 0.0;
  }
  else if ( k == 0 )
  {
    value = 1.0;
  }
  else if ( k == 1 )
  {
    value = ( double ) n;
  }
  else if ( 1 < k && k < n-1 )
  {
    arg = ( double ) ( n + 1 );
    facn = lgamma ( arg );

    arg = ( double ) ( k + 1 );
    fack = lgamma ( arg );

    arg = ( double ) ( n - k + 1 );
    facnmk = lgamma ( arg );

    value = ( int ) ( 0.5 + exp ( facn - fack - facnmk ) );
  }
  else if ( k == n-1 )
  {
    value = ( double ) n;
  }
  else if ( k == n )
  {
    value = 1.0;
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

double r8_cot ( double angle )

//****************************************************************************80
//
//  Purpose:
//
//    R8_COT returns the cotangent of an angle.
//
//  Discussion:
//
//    R8_COT ( THETA ) = COS ( THETA ) / SIN ( THETA )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, the angle, in radians.
//
//    Output, double R8_COT, the cotangent of the angle.
//
{
  double value;

  value = cos ( angle ) / sin ( angle );

  return value;
}
//****************************************************************************80

double r8_cot_deg ( double angle )

//****************************************************************************80
//
//  Purpose:
//
//    R8_COT_DEG returns the cotangent of an angle given in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, the angle, in degrees.
//
//    Output, double R8_COT_DEG, the cotangent of the angle.
//
{
  double value;

  angle = r8_pi ( ) * angle / 180.0;

  value = cos ( angle ) / sin ( angle );

  return value;
}
//****************************************************************************80

double r8_csc ( double theta )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CSC returns the cosecant of X.
//
//  Discussion:
//
//    R8_CSC ( THETA ) = 1.0 / SIN ( THETA )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double THETA, the angle, in radians, whose cosecant is desired.
//    It must be the case that SIN ( THETA ) is not zero.
//
//    Output, double R8_CSC, the cosecant of THETA.
//
{
  double value;

  value = 1.0 / sin ( theta );

  return value;
}
//****************************************************************************80

double r8_csc_deg ( double angle )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CSC_DEG returns the cosecant of an angle given in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, the angle, in degrees.
//
//    Output, double R8_CSC_DEG, the cosecant of the angle.
//
{
  double degrees_to_radians = 3.141592653589793 / 180.0;
  double value;

  value = 1.0 / cos ( degrees_to_radians * angle );

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
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 roundoff unit.
//
{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL returns the factorial function as an R8.
//
//  Discussion:
//
//    Factorial ( N ) = Product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    0 <= N.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  double fact;
  int i;
//
//  Check.
//
  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "R8_FACTORIAL - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  fact = 1.0;

  for ( i = 2; i <= n; i++ )
  {
    fact = fact * ( double ) i;
  }

  return fact;
}
//****************************************************************************80

double r8_factorial_log ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_LOG computes the natural logarithm of the factorial function.
//
//  Discussion:
//
//    LOG ( FACTORIAL ( N ) )
//      = LOG ( product ( 1 <= I <= N ) I )
//      = sum ( ( 1 <= I <= N ) LOG ( I ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, R8_FACTORIAL_LOG is returned as 0.
//
//    Output, double R8_FACTORIAL_LOG, the logarithm of the factorial of N.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value + log ( ( double ) i );
  }

  return value;
}
//****************************************************************************80

void r8_factorial_log_values ( int *n_data, int *n, double *fn )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_LOG_VALUES returns values of log(factorial(n)).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to the index of the test data.  On each subsequent call, N_DATA is
//    incremented and that test data is returned.  When there is no more
//    test data, N_DATA is set to 0.
//
//    Output, int *N, the argument of the function.
//
//    Output, double *FN, the value of the function.
//
{
# define N_MAX 27

  double fnvec[N_MAX] = {
      0.0,         0.0,       0.6931472,  1.791757, 
      3.178051,    4.787489,  6.579246,   8.525160, 
     10.60460,    12.80182,  15.10441,   17.50232,   
     19.98722,    22.55216,  25.19123,   27.89927, 
     30.67186,    33.50508,  36.39544,   39.33987, 
     42.33561,    58.00362, 148.4778,   363.7394, 
    605.0201,   2611.331,   5912.128 };
  int nvec[N_MAX] = {
     0,   1,    2,   3, 
     4,   5,    6,   7, 
     8,   9,   10,  11, 
    12,  13,   14,  15, 
    16,  17,   18,  19, 
    20,  25,   50, 100, 
   150, 500, 1000 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *fn = 0.0;
  }
  else
  {
    *n = nvec[*n_data];
    *fn = fnvec[*n_data];
    *n_data = *n_data + 1;
  }

  return;

# undef N_MAX
}
//****************************************************************************80

void r8_factorial_values ( int *n_data, int *n, double *fn )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_VALUES returns values of the real factorial function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to the index of the test data.  On each subsequent call, N_DATA is
//    incremented and that test data is returned.  When there is no more
//    test data, N_DATA is set to 0.
//
//    Output, int *N, the argument of the function.
//
//    Output, double *FN, the value of the function.
//
{
# define N_MAX 23

  double fnvec[N_MAX] = {
    1.0, 1.0, 2.0, 6.0, 
    24.0, 120.0, 720.0, 5040.0, 
    40320.0, 362880.0, 3628800.0, 39916800.0, 
    479001600.0, 6227020800.0, 87178291200.0, 1307674368000.0, 
    2.0922789888E+13, 3.5568742810E+14, 6.4023737057E+15, 1.2164510041E+17, 
    2.4329020082E+18, 1.5511210043E+25, 2.6525285981E+32 };
  int nvec[N_MAX] = {
     0,  1,  2,  3, 
     4,  5,  6,  7, 
     8,  9, 10, 11, 
    12, 13, 14, 15, 
    16, 17, 18, 19, 
    20, 25, 30 };
//
  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *fn = 0.0;
  }
  else
  {
    *n = nvec[*n_data];
    *fn = fnvec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double r8_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL2 computes the double factorial function.
//
//  Discussion:
//
//    The formula is:
//
//      FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                      = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N    Factorial2(N)
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the double factorial 
//    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
//
//    Output, double R8_FACTORIAL2, the value of the double factorial.
//
{
  double r8_n;
  double value;

  if ( n < 1 )
  {
    value = 1.0;
    return value;
  }

  r8_n = ( double ) ( n );
  value = 1.0;

  while ( 1.0 < r8_n )
  {
    value = value * r8_n;
    r8_n = r8_n - 2.0;
  }

  return value;
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    The C MATH library includes a function GAMMA ( X ) which should be
//    invoked instead of this function.
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = false;
  fact = 1.0;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
//****************************************************************************80

double r8_gamma_log ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
//
//  Discussion:
//
//    The C MATH library includes a function LGAMMA ( X ) which should be
//    invoked instead of this function.
//
//    The program uses rational functions that theoretically approximate
//    LOG(GAMMA(X)) to at least 18 significant decimal digits.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody and Kenneth Hillstrom,
//    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
//    Mathematics of Computation,
//    Volume 21, 1967, pages 198-203.
//
//    Kenneth Hillstrom,
//    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
//    May 1969.
//
//    Hart, Et. Al.,
//    Computer Approximations,
//    Wiley and sons, New York, 1968.
//
//  Parameters:
//
//    Input, double X, the argument of the Gamma function.  X must be positive.
//
//    Output, double R8_GAMMA_LOG, the logarithm of the Gamma function of X.
//    If X <= 0.0, or if overflow would occur, the program returns the
//    value XINF, the largest representable floating point number.
//
//  Local Parameters:
//
//  BETA   - radix for the floating-point representation.
//
//  MAXEXP - the smallest positive power of BETA that overflows.
//
//  XBIG   - largest argument for which LN(GAMMA(X)) is representable
//           in the machine, i.e., the solution to the equation
//             LN(GAMMA(XBIG)) = BETA**MAXEXP.
//
//  FRTBIG - Rough estimate of the fourth root of XBIG
//
//
//  Approximate values for some important machines are:
//
//                            BETA      MAXEXP         XBIG
//
//  CRAY-1        (S.P.)        2        8191       9.62E+2461
//  Cyber 180/855
//    under NOS   (S.P.)        2        1070       1.72E+319
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)        2         128       4.08E+36
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)        2        1024       2.55D+305
//  IBM 3033      (D.P.)       16          63       4.29D+73
//  VAX D-Format  (D.P.)        2         127       2.05D+36
//  VAX G-Format  (D.P.)        2        1023       1.28D+305
//
//
//                           FRTBIG
//
//  CRAY-1        (S.P.)   3.13E+615
//  Cyber 180/855
//    under NOS   (S.P.)   6.44E+79
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)   1.42E+9
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)   2.25D+76
//  IBM 3033      (D.P.)   2.56D+18
//  VAX D-Format  (D.P.)   1.20D+9
//  VAX G-Format  (D.P.)   1.89D+76
//
{
  double c[7] = {
    -1.910444077728E-03, 
     8.4171387781295E-04, 
    -5.952379913043012E-04, 
     7.93650793500350248E-04, 
    -2.777777777777681622553E-03, 
     8.333333333333333331554247E-02, 
     5.7083835261E-03 };
  double corr;
  double d1 = - 5.772156649015328605195174E-01;
  double d2 =   4.227843350984671393993777E-01;
  double d4 =   1.791759469228055000094023;
  double frtbig = 1.42E+09;
  int i;
  double p1[8] = {
    4.945235359296727046734888, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = {
    4.974607845568932035012064, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+010,
    1.702665737765398868392998E+011,
    4.926125793377430887588120E+011, 
    5.606251856223951465078242E+011 };
  double pnt68 = 0.6796875;
  double q1[8] = {
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = {
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = {
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+010, 
    1.016803586272438228077304E+011, 
    3.417476345507377132798597E+011, 
    4.463158187419713286462081E+011 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double xbig = 4.08E+36;
  double xden;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double xsq;
//
//  Return immediately if the argument is out of range.
//
  if ( x <= 0.0 || xbig < x )
  {
    return r8_huge ( );
  }

  if ( x <= r8_epsilon ( ) )
  {
    res = - log ( x );
  }
  else if ( x <= 1.5 )
  {
    if ( x < pnt68 )
    {
      corr = - log ( x );
      xm1 = x;
    }
    else
    {
      corr = 0.0;
      xm1 = ( x - 0.5 ) - 0.5;
    }

    if ( x <= 0.5 || pnt68 <= x )
    {
      xden = 1.0;
      xnum = 0.0;

      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm1 + p1[i];
        xden = xden * xm1 + q1[i];
      }

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
    }
    else
    {
      xm2 = ( x - 0.5 ) - 0.5;
      xden = 1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm2 + p2[i];
        xden = xden * xm2 + q2[i];
      }

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );

    }
  }
  else if ( x <= 4.0 )
  {
    xm2 = x - 2.0;
    xden = 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm2 + p2[i];
      xden = xden * xm2 + q2[i];
    }

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
  }
  else if ( x <= 12.0 )
  {
    xm4 = x - 4.0;
    xden = - 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm4 + p4[i];
      xden = xden * xm4 + q4[i];
    }

    res = d4 + xm4 * ( xnum / xden );
  }
  else
  {
    res = 0.0;

    if ( x <= frtbig )
    {

      res = c[6];
      xsq = x * x;

      for ( i = 0; i < 6; i++ )
      {
        res = res / xsq + c[i];
      }

    }

    res = res / x;
    corr = log ( x );
    res = res + sqrtpi - 0.5 * corr;
    res = res + x * ( corr - 1.0 );

  }

  return res;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" real value.
//
{
  return 1.0E+30;
}
//****************************************************************************80

double r8_hyper_2f1 ( double a, double b, double c, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
//
//  Discussion:
//
//    A minor bug was corrected.  The HW variable, used in several places as
//    the "old" value of a quantity being iteratively improved, was not
//    being initialized.  JVB, 11 February 2008.
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
//    Original FORTRAN77 original by Shanjie Zhang, Jianming Jin.
//    C++ version by John Burkardt.
//
//    The original FORTRAN77 version of this routine is copyrighted by
//    Shanjie Zhang and Jianming Jin.  However, they give permission to
//    incorporate this routine into a user program provided that the copyright
//    is acknowledged.
//
//  Reference:
//
//    Shanjie Zhang, Jianming Jin,
//    Computation of Special Functions,
//    Wiley, 1996,
//    ISBN: 0-471-11963-6,
//    LC: QA351.C45
//
//  Parameters:
//
//    Input, double A, B, C, X, the arguments of the function.
//    C must not be equal to a nonpositive integer.
//    X < 1.
//
//    Output, double R8_HYPER_2F1, the value of the function.
//
{
  double a0;
  double aa;
  double bb;
  double c0;
  double c1;
  double el = 0.5772156649015329;
  double eps;
  double f0;
  double f1;
  double g0;
  double g1;
  double g2;
  double g3;
  double ga;
  double gabc;
  double gam;
  double gb;
  double gbm;
  double gc;
  double gca;
  double gcab;
  double gcb;
  double gm;
  double hf;
  double hw;
  int j;
  int k;
  bool l0;
  bool l1;
  bool l2;
  bool l3;
  bool l4;
  bool l5;
  int m;
  int nm;
  double pa;
  double pb;
  double pi = 3.141592653589793;
  double r;
  double r0;
  double r1;
  double rm;
  double rp;
  double sm;
  double sp;
  double sp0;
  double x1;

  l0 = ( c == ( int ) ( c ) ) && ( c < 0.0 );
  l1 = ( 1.0 - x < 1.0E-15 ) && ( c - a - b <= 0.0 );
  l2 = ( a == ( int ) ( a ) ) && ( a < 0.0 );
  l3 = ( b == ( int ) ( b ) ) && ( b < 0.0 );
  l4 = ( c - a == ( int ) ( c - a ) ) && ( c - a <= 0.0 );
  l5 = ( c - b == ( int ) ( c - b ) ) && ( c - b <= 0.0 );

  if ( l0 )
  {
    cerr << "\n";
    cerr << "R8_HYPER_2F1 - Fatal error!\n";
    cerr << "  The hypergeometric series is divergent.\n";
    cerr << "  C is integral and negative.\n";
    cerr << "  C = " << c << "\n";
    exit ( 1 );
  }

  if ( l1 )
  {
    cerr << "\n";
    cerr << "R8_HYPER_2F1 - Fatal error!\n";
    cerr << "  The hypergeometric series is divergent.\n";
    cerr << "  1 - X < 0, C - A - B <= 0\n";
    cerr << "  A = " << a << "\n";
    cerr << "  B = " << b << "\n";
    cerr << "  C = " << c << "\n";
    cerr << "  X = " << x << "\n";
    exit ( 1 );
  }

  if ( 0.95 < x )
  {
    eps = 1.0E-08;
  }
  else
  {
    eps = 1.0E-15;
  }

  if ( x == 0.0 || a == 0.0 || b == 0.0 )
  {
    hf = 1.0;
    return hf;
  }
  else if ( 1.0 - x == eps && 0.0 < c - a - b )
  {
    gc = gamma ( c );
    gcab = gamma ( c - a - b );
    gca = gamma ( c - a );
    gcb = gamma ( c - b );
    hf = gc * gcab / ( gca * gcb );
    return hf;
  }
  else if ( 1.0 + x <= eps && r8_abs ( c - a + b - 1.0 ) <= eps )
  {
    g0 = sqrt ( pi ) * pow ( 2.0, - a );
    g1 = gamma ( c );
    g2 = gamma ( 1.0 + a / 2.0 - b );
    g3 = gamma ( 0.5 + 0.5 * a );
    hf = g0 * g1 / ( g2 * g3 );
    return hf;
  }
  else if ( l2 || l3 )
  {
    if ( l2 )
    {
      nm = ( int ) ( r8_abs ( a ) );
    }

    if ( l3 )
    {
      nm = ( int ) ( r8_abs ( b ) );
    }

    hf = 1.0;
    r = 1.0;

    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }

    return hf;
  }
  else if ( l4 || l5 )
  {
    if ( l4 )
    {
      nm = ( int ) ( r8_abs ( c - a ) );
    }

    if ( l5 )
    {
      nm = ( int ) ( r8_abs ( c - b ) );
    }

    hf = 1.0;
    r  = 1.0;
    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }
    hf = pow ( 1.0 - x, c - a - b ) * hf;
    return hf;
  }

  aa = a;
  bb = b;
  x1 = x;

  if ( x < 0.0 )
  {
    x = x / ( x - 1.0 );
    if ( a < c && b < a && 0.0 < b )
    {
      a = bb;
      b = aa;
    }
    b = c - b;
  }

  if ( 0.75 <= x )
  {
    gm = 0.0;

    if ( r8_abs ( c - a - b - ( int ) ( c - a - b ) ) < 1.0E-15 )
    {
      m = int ( c - a - b );
      ga = gamma ( a );
      gb = gamma ( b );
      gc = gamma ( c );
      gam = gamma ( a + m );
      gbm = gamma ( b + m );

      pa = r8_psi ( a );
      pb = r8_psi ( b );

      if ( m != 0 )
      {
        gm = 1.0;
      }

      for ( j = 1; j <= abs ( m ) - 1; j++ )
      {
        gm = gm * j;
      }

      rm = 1.0;
      for ( j = 1; j <= abs ( m ); j++ )
      {
        rm = rm * j;
      }

      f0 = 1.0;
      r0 = 1.0;;
      r1 = 1.0;
      sp0 = 0.0;;
      sp = 0.0;

      if ( 0 <= m )
      {
        c0 = gm * gc / ( gam * gbm );
        c1 = - gc * pow ( x - 1.0, m ) / ( ga * gb * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( a + k - 1.0 ) + 1.0 / ( b + k - 1.0 ) 
          - 1.0 / ( double ) ( k );
        }

        f1 = pa + pb + sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + ( 1.0 - a ) 
              / ( ( j + k ) * ( a + j + k - 1.0 ) ) 
              + 1.0 / ( b + j + k - 1.0 );
          }

          rp = pa + pb + 2.0 * el + sp + sm + log ( 1.0 - x );

          r1 = r1 * ( a + m + k - 1.0 ) * ( b + m + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }
          hw = f1;
        }
        hf = f0 * c0 + f1 * c1;
      }
      else if ( m < 0 )
      {
        m = - m;
        c0 = gm * gc / ( ga * gb * pow ( 1.0 - x, m ) );
        c1 = - pow ( - 1.0, m ) * gc / ( gam * gbm * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a - m + k - 1.0 ) * ( b - m + k - 1.0 ) 
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( double ) ( k );
        }

        f1 = pa + pb - sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) 
            / ( k * ( a + k - 1.0 ) ) 
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + 1.0 / ( double ) ( j + k );
          }

          rp = pa + pb + 2.0 * el + sp - sm + log ( 1.0 - x );

          r1 = r1 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }

          hw = f1;
        }

        hf = f0 * c0 + f1 * c1;
      }
    }
    else
    {
      ga = gamma ( a );
      gb = gamma ( b );
      gc = gamma ( c );
      gca = gamma ( c - a );
      gcb = gamma ( c - b );
      gcab = gamma ( c - a - b );
      gabc = gamma ( a + b - c );
      c0 = gc * gcab / ( gca * gcb );
      c1 = gc * gabc / ( ga * gb ) * pow ( 1.0 - x, c - a - b );
      hf = 0.0;
      hw = hf;
      r0 = c0;
      r1 = c1;

      for ( k = 1; k <= 250; k++ )
      {
        r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
          / ( k * ( a + b - c + k ) ) * ( 1.0 - x );

        r1 = r1 * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 ) 
          / ( k * ( c - a - b + k ) ) * ( 1.0 - x );

        hf = hf + r0 + r1;

        if ( r8_abs ( hf - hw ) < r8_abs ( hf ) * eps )
        {
          break;
        }
        hw = hf;
      }
      hf = hf + c0 + c1;
    }
  }
  else
  {
    a0 = 1.0;

    if ( a < c && c < 2.0 * a && b < c && c < 2.0 * b )
    {
      a0 = pow ( 1.0 - x, c - a - b );
      a = c - a;
      b = c - b;
    }

    hf = 1.0;
    hw = hf;
    r = 1.0;

    for ( k = 1; k <= 250; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 ) 
        / ( k * ( c + k - 1.0 ) ) * x;

      hf = hf + r;

      if ( r8_abs ( hf - hw ) <= r8_abs ( hf ) * eps )
      {
        break;
      }

      hw = hf;
    }
    hf = a0 * hf;
  }

  if ( x1 < 0.0 )
  {
    x = x1;
    c0 = 1.0 / pow ( 1.0 - x, aa );
    hf = c0 * hf;
  }

  a = aa;
  b = bb;

  if ( 120 < k )
  {
    cout << "\n";
    cout << "R8_HYPER_2F1 - Warning!\n";
    cout << "  A large number of iterations were needed.\n";
    cout << "  The accuracy of the results should be checked.\n";
  }

  return hf;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
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
  if ( y < x )
  {
    return x;
  } 
  else
  {
    return y;
  }
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2003
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
  if ( y < x )
  {
    return y;
  } 
  else
  {
    return x;
  }
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the integer that is nearest to a real value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the real number.
//
//    Output, int R8_NINT, the nearest integer.
//
{
  double d;
  int i;

  i = int ( x );
  d = x - i;

  if ( fabs ( d ) <= 0.5 )
  {
    return i;
  }
  else if ( x < i ) 
  {
    return (i-1);
  }
  else
  {
    return (i+1);
  }
}
//****************************************************************************80

double r8_pi ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PI returns the value of PI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_PI, the value of PI.
//
{
  double value = 3.141592653589793;

  return value;
}
//****************************************************************************80

double r8_psi ( double xx )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PSI evaluates the function Psi(X).
//
//  Discussion:
//
//    This routine evaluates the logarithmic derivative of the
//    Gamma function,
//
//      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
//             = d/dX LN ( GAMMA(X) )
//
//    for real X, where either
//
//      - XMAX1 < X < - XMIN, and X is not a negative integer,
//
//    or
//
//      XMIN < X.
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
//    Original FORTRAN77 version by William Cody
//    C++ version by John Burkardt
//
//  Reference:
//
//    William Cody, Anthony Strecok, Henry Thacher,
//    Chebyshev Approximations for the Psi Function,
//    Mathematics of Computation,
//    Volume 27, Number 121, January 1973, pages 123-127.
//
//  Parameters:
//
//    Input, double XX, the argument of the function.
//
//    Output, double, the value of the function.
//
{
  double aug;
  double den;
  int i;
  int n;
  int nq;
  double p1[9] = { 
   4.5104681245762934160E-03, 
   5.4932855833000385356, 
   3.7646693175929276856E+02, 
   7.9525490849151998065E+03, 
   7.1451595818951933210E+04, 
   3.0655976301987365674E+05, 
   6.3606997788964458797E+05, 
   5.8041312783537569993E+05, 
   1.6585695029761022321E+05 };
  double p2[7] = { 
  -2.7103228277757834192, 
  -1.5166271776896121383E+01, 
  -1.9784554148719218667E+01, 
  -8.8100958828312219821, 
  -1.4479614616899842986, 
  -7.3689600332394549911E-02, 
  -6.5135387732718171306E-21 };
  double piov4 = 0.78539816339744830962;
  double q1[8] = { 
   9.6141654774222358525E+01, 
   2.6287715790581193330E+03, 
   2.9862497022250277920E+04, 
   1.6206566091533671639E+05, 
   4.3487880712768329037E+05, 
   5.4256384537269993733E+05, 
   2.4242185002017985252E+05, 
   6.4155223783576225996E-08 };
  double q2[6] = { 
   4.4992760373789365846E+01, 
   2.0240955312679931159E+02, 
   2.4736979003315290057E+02, 
   1.0742543875702278326E+02, 
   1.7463965060678569906E+01, 
   8.8427520398873480342E-01 };
  double sgn;
  double upper;
  double value;
  double w;
  double x;
  double x01 = 187.0;
  double x01d = 128.0;
  double x02 = 6.9464496836234126266E-04;
  double xinf = 1.70E+38;
  double xlarge = 2.04E+15;
  double xmax1 = 3.60E+16;
  double xmin1 = 5.89E-39;
  double xsmall = 2.05E-09;
  double z;

  x = xx;
  w = r8_abs ( x );
  aug = 0.0;
//
//  Check for valid arguments, then branch to appropriate algorithm.
//
  if ( xmax1 <= - x || w < xmin1 )
  {
    if ( 0.0 < x )
    {
      value = - xinf;
    }
    else
    {
      value = xinf;
    }
    return value;
  }

  if ( x < 0.5 )
  {
//
//  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
//  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
//
    if ( w <= xsmall )
    {
      aug = - 1.0 / x;
    }
//
//  Argument reduction for cotangent.
//
    else
    {
      if ( x < 0.0 )
      {
        sgn = piov4;
      }
      else
      {
        sgn = - piov4;
      }

      w = w - ( double ) ( ( int ) ( w ) );
      nq = int ( w * 4.0 );
      w = 4.0 * ( w - ( double ) ( nq ) * 0.25 );
//
//  W is now related to the fractional part of 4.0 * X.
//  Adjust argument to correspond to values in the first
//  quadrant and determine the sign.
//
      n = nq / 2;

      if ( n + n != nq )
      {
        w = 1.0 - w;
      }

      z = piov4 * w;

      if ( ( n % 2 ) != 0 )
      {
        sgn = - sgn;
      }
//
//  Determine the final value for  -pi * cotan(pi*x).
//
      n = ( nq + 1 ) / 2;
      if ( ( n % 2 ) == 0 )
      {
//
//  Check for singularity.
//
        if ( z == 0.0 )
        {
          if ( 0.0 < x )
          {
            value = -xinf;
          }
          else
          {
            value = xinf;
          }
          return value;
        }
        aug = sgn * ( 4.0 / tan ( z ) );
      }
      else
      {
        aug = sgn * ( 4.0 * tan ( z ) );
      }
    }
    x = 1.0 - x;
  }
//
//  0.5 <= X <= 3.0.
//
  if ( x <= 3.0 )
  {
    den = x;
    upper = p1[0] * x;
    for ( i = 1; i <= 7; i++ )
    {
      den = ( den + q1[i-1] ) * x;
      upper = ( upper + p1[i]) * x;
    }
    den = ( upper + p1[8] ) / ( den + q1[7] );
    x = ( x - x01 / x01d ) - x02;
    value = den * x + aug;
    return value;
  }
//
//  3.0 < X.
//
  if ( x < xlarge )
  {
    w = 1.0 / ( x * x );
    den = w;
    upper = p2[0] * w;
    for ( i = 1; i <= 5; i++ )
    {
      den = ( den + q2[i-1] ) * w;
      upper = ( upper + p2[i] ) * w;
    }
    aug = ( upper + p2[6] ) / ( den + q2[5] ) - 0.5 / x + aug;
  }

  value = aug + log ( x );

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 is a portable pseudorandom number generator.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
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
//    07 March 2003
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
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and the value will be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

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
  r = ( ( double ( *seed ) ) * 4.656612875E-10 );

  return r;
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
//    10 May 2003
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
//    Output, int R8POLY_DEGREE, the degree of the polynomial.
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
//****************************************************************************80*

void r8poly_print ( int n, double a[], string title )

//****************************************************************************80*
//
//  Purpose:
//
//    R8POLY_PRINT prints out a polynomial.
//
//  Discussion:
//
//    The power sum form is:
//
//      p(x) = a(0) + a(1)*x + ... + a(n-1)*x^(n-1) + a(n)*x^(n)
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
//    Input, int N, the dimension of A.
//
//    Input, double A[N+1], the polynomial coefficients.
//    A(0) is the constant term and
//    A(N) is the coefficient of X^N.
//
//    Input, string TITLE, an optional title.
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

  if ( a[n2] < 0.0 )
  {
    plus_minus = '-';
  }
  else
  {
    plus_minus = ' ';
  }

  mag = fabs ( a[n2] );

  if ( 2 <= n2 )
  {
    cout << "  p(x) = " << plus_minus << mag << " * x^" << n2 << "\n";
  }
  else if ( n2 == 1 )
  {
    cout << "  p(x) = " << plus_minus << mag << " * x" << "\n";
  }
  else if ( n2 == 0 )
  {
    cout << "  p(x) = " << plus_minus << mag << "\n";
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

    mag = fabs ( a[i] );

    if ( mag != 0.0 )
    {
      if ( 2 <= i )
      {
        cout << "         " << plus_minus << mag << " * x^" << i << "\n";
      }
      else if ( i == 1 )
      {
        cout << "         " << plus_minus << mag << " * x" << "\n";
      }
      else if ( i == 0 )
      {
        cout << "         " << plus_minus << mag << "\n";
      }
    }

  }

  return;
}
//****************************************************************************80

double r8poly_value ( int n, double a[], double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_VALUE evaluates a double precision polynomial.
//
//  Discussion:
//
//    For sanity's sake, the value of N indicates the NUMBER of 
//    coefficients, or more precisely, the ORDER of the polynomial,
//    rather than the DEGREE of the polynomial.  The two quantities
//    differ by 1, but cause a great deal of confusion.
//
//    Given N and A, the form of the polynomial is:
//
//      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double A[N], the coefficients of the polynomial.
//    A[0] is the constant term.
//
//    Input, double X, the point at which the polynomial is to be evaluated.
//
//    Output, double R8POLY_VALUE, the value of the polynomial at X.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = n-1; 0 <= i; i-- )
  {
    value = value * x + a[i];
  }

  return value;
}
//****************************************************************************80

void r8vec_zero ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO zeroes a real vector.
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
//    Input, int N, the number of entries in the vector.
//
//    Output, double A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
//****************************************************************************80

int s_len_trim ( char* s )

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
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char* t;

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

double sec_deg ( double angle )

//****************************************************************************80
//
//  Purpose:
//
//    SEC_DEG returns the secant of an angle given in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, the angle, in degrees.
//
//    Output, double SEC_DEG, the secant of the angle.
//
{
  double degrees_to_radians = 3.141592653589793 / 180.0;
  double value;

  value = 1.0 / sin ( degrees_to_radians * angle );

  return value;
}
//****************************************************************************80

int sigma ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIGMA returns the value of SIGMA(N), the divisor sum.
//
//  Discussion:
//
//    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
//
//    The formula is:
//
//      SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
//
//      SIGMA(P**K) = ( P^(K+1) - 1 ) / ( P - 1 ) if P is prime.
//
//  First values:
//
//     N  SIGMA(N)
//
//     1    1
//     2    3
//     3    4
//     4    7
//     5    6
//     6   12
//     7    8
//     8   15
//     9   13
//    10   18
//    11   12
//    12   28
//    13   14
//    14   24
//    15   24
//    16   31
//    17   18
//    18   39
//    19   20
//    20   42
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the value to be analyzed.
//
//    Output, int SIGMA, the value of SIGMA(N).  If N is less than
//    or equal to 0, SIGMA will be returned as 0.  If there is not
//    enough room for factoring N, SIGMA is returned as -1.
//
{
# define FACTOR_MAX 20

  int factor[FACTOR_MAX];
  int i;
  int nfactor;
  int nleft;
  int power[FACTOR_MAX];
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  if ( n == 1 )
  {
    return 1;
  }
//
//  Factor N.
//
  i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );

  if ( nleft != 1 )
  {
    cerr << "\n";
    cerr << "SIGMA - Fatal error!\n";
    cerr << "  Not enough factorization space.\n";
    exit ( 1 );
  }

  value = 1;
  for ( i = 0; i < nfactor; i++ )
  {
    value = ( value * 
      ( ( int ) pow ( ( double ) factor[i], power[i] + 1 ) - 1 ) ) 
      / ( factor[i] - 1 );
  }

  return value;
# undef FACTOR_MAX
}
//****************************************************************************80

void sigma_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    SIGMA_VALUES returns some values of the Sigma function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the argument of the Sigma function.
//
//    Output, int *C, the value of the Sigma function.
//
{
# define N_MAX 20

  int c_vec[N_MAX] = {
     1,    3,    4,    7,    6,   12,    8,   15,   13,   18, 
    72,  128,  255,  176,  576, 1170,  618,  984, 2232, 2340 };

  int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, 
     30, 127, 128, 129, 210, 360, 617, 815, 816,1000 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double sin_deg ( double angle )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_DEG returns the sine of an angle given in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, the angle, in degrees.
//
//    Output, double SIN_DEG, the sine of the angle.
//
{
  double degrees_to_radians = 3.141592653589793 / 180.0;
  double value;

  value = sin ( degrees_to_radians * angle );

  return value;
}
//****************************************************************************80

double sin_power_int ( double a, double b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_POWER_INT evaluates the sine power integral.
//
//  Discussion:
//
//    The function is defined by
//
//      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
//
//    The algorithm uses the following fact:
//
//      Integral sin^n ( t ) = (1/n) * (
//        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, double A, B, the limits of integration.
//
//    Input, integer N, the power of the sine function.
//
//    Output, double SIN_POWER_INT, the value of the integral.
//
{
  double ca;
  double cb;
  int m;
  int mlo;
  double sa;
  double sb;
  double value;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "SIN_POWER_INT - Fatal error!\n";
    cerr << "  Power N < 0.\n";
    exit ( 1 );
  }

  sa = sin ( a );
  sb = sin ( b );
  ca = cos ( a );
  cb = cos ( b );

  if ( ( n % 2 ) == 0 )
  {
    value = b - a;
    mlo = 2;
  }
  else
  {
    value = ca - cb;
    mlo = 3;
  }

  for ( m = mlo; m <= n; m = m + 2 )
  {
    value = ( ( double ) ( m - 1 ) * value 
      + pow ( sa, (m-1) ) * ca - pow ( sb, (m-1) ) * cb ) 
      / ( double ) ( m );
  }

  return value;
}
//****************************************************************************80

void sin_power_int_values ( int *n_data, double *a, double *b, int *n, 
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_POWER_INT_VALUES returns some values of the sine power integral.
//
//  Discussion:
//
//    The function has the form
//
//      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin(T) )^N dt
//
//    In Mathematica, the function can be evaluated by:
//
//      Integrate [ ( Sin[x] )^n, { x, a, b } ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, double *B, the limits of integration.
//
//    Output, int *N, the power.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 10

  double a_vec[N_MAX] = {
    0.10E+02, 
    0.00, 
    0.00, 
    0.00, 
    0.00,
    0.00, 
    0.00, 
    0.10E+01, 
    0.00, 
    0.00 };
  double b_vec[N_MAX] = {
    0.20E+02,
    0.10E+01,
    0.10E+01,
    0.10E+01,
    0.10E+01,
    0.10E+01,
    0.20E+01,
    0.20E+01,
    0.10E+01,
    0.10E+01 };
  double fx_vec[N_MAX] = {
    0.10000000000000000000E+02, 
    0.45969769413186028260, 
    0.27267564329357957615, 
    0.17894056254885809051, 
    0.12402556531520681830, 
    0.88974396451575946519E-01, 
    0.90393123848149944133, 
    0.81495684202992349481, 
    0.21887522421729849008E-01, 
    0.17023439374069324596E-01 };
  int n_vec[N_MAX] = {
     0, 
     1, 
     2, 
     3, 
     4, 
     5, 
     5, 
     5, 
    10, 
    11 };
//
  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0;
    *b = 0.0;
    *n = 0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *n = n_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int slice ( int dim_num, int slice_num )

//****************************************************************************80
//
//  Purpose:
//
//    SLICE: maximum number of pieces created by a given number of slices.
//
//  Discussion:
//
//    If we imagine slicing a pizza, each slice produce more pieces.  
//    The position of the slice affects the number of pieces created, but there
//    is a maximum.  
//
//    This function determines the maximum number of pieces created by a given
//    number of slices, applied to a space of a given dimension.
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
//  Reference:
//
//    Robert Banks,
//    Slicing Pizzas, Racing Turtles, and Further Adventures in 
//    Applied Mathematics,
//    Princeton, 1999,
//    ISBN13: 9780691059471,
//    LC: QA93.B358.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int SLICE_NUM, the number of slices.
//
//    Input, int SLICE, the maximum number of pieces that can
//    be created by the given number of slices applied in the given dimension.
//
{
  int j;
  int piece_num;

  piece_num = 0;
  for ( j = 0; j <= i4_min ( dim_num, slice_num ); j++ )
  {
    piece_num = piece_num + i4_choose ( slice_num, j );
  }

  return piece_num;
}
//****************************************************************************80

void spherical_harmonic ( int l, int m, double theta, double phi, 
  double c[], double s[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERICAL_HARMONIC evaluates spherical harmonic functions.
//
//  Discussion:
//
//    The spherical harmonic function Y(L,M,THETA,PHI,X) is the
//    angular part of the solution to Laplace's equation in spherical
//    coordinates.
//
//    Y(L,M,THETA,PHI,X) is related to the associated Legendre
//    function as follows:
//
//      Y(L,M,THETA,PHI,X) = FACTOR * P(L,M,cos(THETA)) * exp ( i * M * PHI )
//
//    Here, FACTOR is a normalization factor:
//
//      FACTOR = sqrt ( ( 2 * L + 1 ) * ( L - M )! / ( 4 * PI * ( L + M )! ) )
//
//    In Mathematica, a spherical harmonic function can be evaluated by
//
//      SphericalHarmonicY [ l, m, theta, phi ]
//
//    Note that notational tradition in physics requires that THETA
//    and PHI represent the reverse of what they would normally mean
//    in mathematical notation; that is, THETA goes up and down, and
//    PHI goes around.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2005
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
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input, int L, the first index of the spherical harmonic function.
//    Normally, 0 <= L.
//
//    Input, int M, the second index of the spherical harmonic function.
//    Normally, -L <= M <= L.
//
//    Input, double THETA, the polar angle, for which
//    0 <= THETA <= PI.
//
//    Input, double PHI, the longitudinal angle, for which
//    0 <= PHI <= 2*PI.
//
//    Output, double C[L+1], S[L+1], the real and imaginary
//    parts of the functions Y(L,0:L,THETA,PHI).
//
{
  double angle;
  int i;
  int m_abs;
  double *plm;

  m_abs = abs ( m );

  plm = new double[l+1];

  legendre_associated_normalized ( l, m_abs, cos ( theta ), plm );

  angle = ( double ) ( m ) * phi;

  if ( 0 <= m )
  {
    for ( i = 0; i <= l; i++ )
    {
      c[i] = plm[i] * cos ( angle );
      s[i] = plm[i] * sin ( angle );
    }
  }
  else
  {
    for ( i = 0; i <= l; i++ )
    {
      c[i] = -plm[i] * cos ( angle );
      s[i] = -plm[i] * sin ( angle );
    }
  }

  delete [] plm;

  return;
}
//****************************************************************************80

void spherical_harmonic_values ( int *n_data, int *l, int *m, double *theta,
  double *phi, double *yr, double *yi ) 

//****************************************************************************80
//
//  Purpose:
//
//    SPHERICAL_HARMONIC_VALUES returns values of spherical harmonic functions.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by
//
//      SphericalHarmonicY [ l, m, theta, phi ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 March 2005
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
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *L, int *M, double THETA, PHI, the arguments
//    of the function.
//
//    Output, double *YR, *YI, the real and imaginary parts of
//    the function.
//
{
# define N_MAX 20

  int l_vec[N_MAX] = {
     0,  1,  2,  
     3,  4,  5,  
     5,  5,  5,  
     5,  4,  4,  
     4,  4,  4,  
     3,  3,  3,  
     3,  3 };
  int m_vec[N_MAX] = {
     0,  0,  1,  
     2,  3,  5,  
     4,  3,  2,  
     1,  2,  2,  
     2,  2,  2,  
    -1, -1, -1,  
    -1, -1 };
  double phi_vec[N_MAX] = {
    0.1047197551196598E+01, 0.1047197551196598E+01, 0.1047197551196598E+01, 
    0.1047197551196598E+01, 0.1047197551196598E+01, 0.6283185307179586, 
    0.6283185307179586, 0.6283185307179586, 0.6283185307179586, 
    0.6283185307179586, 0.7853981633974483, 0.7853981633974483, 
    0.7853981633974483, 0.7853981633974483, 0.7853981633974483, 
    0.4487989505128276, 0.8975979010256552, 0.1346396851538483E+01, 
    0.1795195802051310E+01, 0.2243994752564138E+01 };
  double theta_vec[N_MAX] = {
    0.5235987755982989, 0.5235987755982989, 0.5235987755982989, 
    0.5235987755982989, 0.5235987755982989, 0.2617993877991494, 
    0.2617993877991494, 0.2617993877991494, 0.2617993877991494, 
    0.2617993877991494, 0.6283185307179586, 0.1884955592153876E+01, 
    0.3141592653589793E+01, 0.4398229715025711E+01, 0.5654866776461628E+01, 
    0.3926990816987242, 0.3926990816987242, 0.3926990816987242, 
    0.3926990816987242, 0.3926990816987242 };
  double yi_vec[N_MAX] = {
    0.0000000000000000,  0.0000000000000000, -0.2897056515173922, 
    0.1916222768312404,  0.0000000000000000,  0.0000000000000000, 
    0.3739289485283311E-02, -0.4219517552320796E-01,  0.1876264225575173, 
   -0.3029973424491321,  0.4139385503112256, -0.1003229830187463, 
    0.0000000000000000, -0.1003229830187463,  0.4139385503112256, 
   -0.1753512375142586, -0.3159720118970196, -0.3940106541811563, 
   -0.3940106541811563, -0.3159720118970196 };
  double yr_vec[N_MAX] = {
   0.2820947917738781,  0.4231421876608172, -0.1672616358893223, 
  -0.1106331731112457,  0.1354974113737760,  0.5390423109043568E-03, 
  -0.5146690442951909E-02,  0.1371004361349490E-01,  0.6096352022265540E-01, 
  -0.4170400640977983,  0.0000000000000000,  0.0000000000000000, 
   0.0000000000000000,  0.0000000000000000,  0.0000000000000000, 
   0.3641205966137958,  0.2519792711195075,  0.8993036065704300E-01, 
  -0.8993036065704300E-01, -0.2519792711195075 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *l = 0;
    *m = 0;
    *theta = 0.0;
    *phi = 0.0;
    *yr = 0.0;
    *yi = 0.0;
  }
  else
  {
    *l = l_vec[*n_data-1];
    *m = m_vec[*n_data-1];
    *theta = theta_vec[*n_data-1];
    *phi = phi_vec[*n_data-1];
    *yr = yr_vec[*n_data-1];
    *yi = yi_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int *stirling1 ( int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    STIRLING1 computes the Stirling numbers of the first kind.
//
//  Discussion:
//
//    The absolute value of the Stirling number S1(N,M) gives the number
//    of permutations on N objects having exactly M cycles, while the
//    sign of the Stirling number records the sign (odd or even) of
//    the permutations.  For example, there are six permutations on 3 objects:
//
//      A B C   3 cycles (A) (B) (C)
//      A C B   2 cycles (A) (BC)
//      B A C   2 cycles (AB) (C)
//      B C A   1 cycle  (ABC)
//      C A B   1 cycle  (ABC)
//      C B A   2 cycles (AC) (B)
//
//    There are
//
//      2 permutations with 1 cycle, and S1(3,1) = 2
//      3 permutations with 2 cycles, and S1(3,2) = -3,
//      1 permutation with 3 cycles, and S1(3,3) = 1.
//
//    Since there are N! permutations of N objects, the sum of the absolute
//    values of the Stirling numbers in a given row,
//
//      sum ( 1 <= I <= N ) abs ( S1(N,I) ) = N!
//
//  First terms:
//
//    N/M:  1     2      3     4     5    6    7    8
//
//    1     1     0      0     0     0    0    0    0
//    2    -1     1      0     0     0    0    0    0
//    3     2    -3      1     0     0    0    0    0
//    4    -6    11     -6     1     0    0    0    0
//    5    24   -50     35   -10     1    0    0    0
//    6  -120   274   -225    85   -15    1    0    0
//    7   720 -1764   1624  -735   175  -21    1    0
//    8 -5040 13068 -13132  6769 -1960  322  -28    1
//
//  Recursion:
//
//    S1(N,1) = (-1)^(N-1) * (N-1)! for all N.
//    S1(I,I) = 1 for all I.
//    S1(I,J) = 0 if I < J.
//
//    S1(N,M) = S1(N-1,M-1) - (N-1) * S1(N-1,M)
//
//  Properties:
//
//    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
//
//    X_N = sum ( 0 <= K <= N ) S1(N,K) X^K
//    where X_N is the falling factorial function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows of the table.
//
//    Input, int M, the number of columns of the table.
//
//    Output, int STIRLING1[N*M], the Stirling numbers of the first kind.
//
{
  int i;
  int j;
  int *s1;

  if ( n <= 0 )
  {
    return NULL;
  }

  if ( m <= 0 )
  {
    return NULL;
  }

  s1 = new int[n*m];

  s1[0+0*n] = 1;
  for ( j = 2; j <= m; j++ )
  {
    s1[0+(j-1)*n] = 0;
  }

  for ( i = 2; i <= n; i++ )
  {
    s1[i-1+0*n] = - ( i - 1 ) * s1[i-2+0*n];
    for ( j = 2; j <= m; j++ )
    {
      s1[i-1+(j-1)*n] = s1[i-2+(j-2)*n] - ( i - 1 ) * s1[i-2+(j-1)*n];
    }

  }

  return s1;
}
//****************************************************************************80

int *stirling2 ( int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    STIRLING2 computes the Stirling numbers of the second kind.
//
//  Discussion:
//
//    S2(N,M) represents the number of distinct partitions of N elements
//    into M nonempty sets.  For a fixed N, the sum of the Stirling
//    numbers S2(N,M) is represented by B(N), called "Bell's number",
//    and represents the number of distinct partitions of N elements.
//
//    For example, with 4 objects, there are:
//
//    1 partition into 1 set:
//
//      (A,B,C,D)
//
//    7 partitions into 2 sets:
//
//      (A,B,C) (D)
//      (A,B,D) (C)
//      (A,C,D) (B)
//      (A) (B,C,D)
//      (A,B) (C,D)
//      (A,C) (B,D)
//      (A,D) (B,C)
//
//    6 partitions into 3 sets:
//
//      (A,B) (C) (D)
//      (A) (B,C) (D)
//      (A) (B) (C,D)
//      (A,C) (B) (D)
//      (A,D) (B) (C)
//      (A) (B,D) (C)
//
//    1 partition into 4 sets:
//
//      (A) (B) (C) (D)
//
//    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
//
//
//  First terms:
//
//    N/M: 1    2    3    4    5    6    7    8
//
//    1    1    0    0    0    0    0    0    0
//    2    1    1    0    0    0    0    0    0
//    3    1    3    1    0    0    0    0    0
//    4    1    7    6    1    0    0    0    0
//    5    1   15   25   10    1    0    0    0
//    6    1   31   90   65   15    1    0    0
//    7    1   63  301  350  140   21    1    0
//    8    1  127  966 1701 1050  266   28    1
//
//  Recursion:
//
//    S2(N,1) = 1 for all N.
//    S2(I,I) = 1 for all I.
//    S2(I,J) = 0 if I < J.
//
//    S2(N,M) = M * S2(N-1,M) + S2(N-1,M-1)
//
//  Properties:
//
//    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
//
//    X^N = sum ( 0 <= K <= N ) S2(N,K) X_K
//    where X_K is the falling factorial function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows of the table.
//
//    Input, int M, the number of columns of the table.
//
//    Output, int STIRLING2[N*M], the Stirling numbers of the second kind.
//
{
  int i;
  int j;
  int *s2;

  if ( n <= 0 )
  {
    return NULL;
  }

  if ( m <= 0 )
  {
    return NULL;
  }

  s2 = new int[n*m];

  s2[0+(0)*n] = 1;
  for ( j = 2; j <= m; j++ )
  {
    s2[0+(j-1)*n] = 0;
  }

  for ( i = 2; i <= n; i++ )
  {
    s2[i-1+(0)*n] = 1;

    for ( j = 2; j <= m; j++ )
    {
      s2[i-1+(j-1)*n] = j * s2[i-2+(j-1)*n] + s2[i-2+(j-2)*n];
    }

  }

  return s2;
}
//****************************************************************************80

double tan_deg ( double angle )

//****************************************************************************80
//
//  Purpose:
//
//    TAN_DEG returns the tangent of an angle given in degrees.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, the angle, in degrees.
//
//    Output, double TAN_DEG, the tangent of the angle.
//
{
  return ( tan ( r8_pi ( ) * angle / 180.0 ) );
}
//****************************************************************************80

int tau ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TAU returns the value of TAU(N), the number of distinct divisors of N.
//
//  Discussion:
//
//    TAU(N) is the number of divisors of N, including 1 and N.
//
//    The formula is:
//
//      If the prime factorization of N is
//
//        N = P1^E1 * P2^E2 * ... * PM^EM,
//
//      then
//
//        TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
//
//  First values:
//
//     N   TAU(N)
//
//     1    1
//     2    2
//     3    2
//     4    3
//     5    2
//     6    4
//     7    2
//     8    4
//     9    3
//    10    4
//    11    2
//    12    6
//    13    2
//    14    4
//    15    4
//    16    5
//    17    2
//    18    6
//    19    2
//    20    6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the value to be analyzed.  N must be 1 or
//    greater.
//
//    Output, int TAU, the value of TAU(N).  But if N is 0 or
//    less, TAU is returned as 0, a nonsense value.  If there is
//    not enough room for factoring, TAU is returned as -1.
//
{
# define FACTOR_MAX 20

  int factor[FACTOR_MAX];
  int i;
  int nfactor;
  int nleft;
  int power[FACTOR_MAX];
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  if ( n == 1 )
  {
    return 1;
  }
//
//  Factor N.
//
  i4_factor ( n, FACTOR_MAX, &nfactor, factor, power, &nleft );

  if ( nleft != 1 )
  {
    cerr << "\n";
    cerr << "TAU - Fatal error!\n";
    cerr << "  Not enough factorization space.\n";
    exit ( 1 );
  }

  value = 1;
  for ( i = 0; i < nfactor; i++ )
  {
    value = value * ( power[i] + 1 );
  }

  return value;
}
//****************************************************************************80

void tau_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    TAU_VALUES returns some values of the Tau function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2003
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
//  Parameters:
//
//    Input/output, int *N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int *N, the argument of the Tau function.
//
//    Output, int *C, the value of the Tau function.
//
{
# define N_MAX 20

  int c_vec[N_MAX] = {
    1,  2,  2,  3,  2,  4,  2,  4,  3,  4, 
    2, 12, 12,  4, 18, 24,  2,  8, 14, 28 };
  int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10, 
     23,  72, 126, 226, 300, 480, 521, 610, 832, 960 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  if ( N_MAX <= *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data];
    *c = c_vec[*n_data];
    *n_data = *n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int tetrahedron_num ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_NUM returns the N-th tetrahedral number.
//
//  Discussion:
//
//    The N-th tetrahedral number T3(N) is formed by the sum of the first
//    N triangular numbers:
//
//      T3(N) = sum ( 1 <= I <= N ) T2(I)
//            = sum ( 1 <= I <= N ) sum ( 1 <= J < I ) J
//
//    By convention, T3(0) = 0.
//
//    The formula is:
//
//      T3(N) = ( N * ( N + 1 ) * ( N + 2 ) ) / 6
//
//  First Values:
//
//     0
//     1
//     4
//    10
//    20
//    35
//    56
//    84
//   120
//   165
//   220
//   275
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the index of the desired number, which must be
//    at least 0.
//
//    Output, int TETRAHEDRON_NUM, the N-th tetrahedron number.
//
{
  int value;

  value = ( n * ( n + 1 ) * ( n + 2 ) ) / 6;

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

int triangle_num ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NUM returns the N-th triangular number.
//
//  Discussion:
//
//    The N-th triangular number T(N) is formed by the sum of the first
//    N integers:
//
//      T(N) = sum ( 1 <= I <= N ) I
//
//    By convention, T(0) = 0.
//
//    The formula is:
//
//      T(N) = ( N * ( N + 1 ) ) / 2
//
//  First Values:
//
//     0
//     1
//     3
//     6
//    10
//    15
//    21
//    28
//    36
//    45
//    55
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the index of the desired number, which must be
//    at least 0.
//
//    Output, int TRIANGLE_NUM, the N-th triangular number.
//
{
  int value;

  value = ( n * ( n + 1 ) ) / 2;

  return value;
}
//****************************************************************************80

int triangle_to_i4 ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_TO_I4 converts a triangular coordinate to an integer.
//
//  Discussion:
//
//    Triangular coordinates are handy when storing a naturally triangular
//    array (such as the lower half of a matrix) in a linear array.
//
//    Thus, for example, we might consider storing
//
//    (0,0)
//    (1,0) (1,1)
//    (2,0) (2,1) (2,2)
//    (3,0) (3,1) (3,2) (3,3)
//
//    as the linear array
//
//    (0,0) (1,0) (1,1) (2,0) (2,1) (2,2) (3,0) (3,1) (3,2) (3,3)
//
//    Here, the quantities in parenthesis represent the natural row and
//    column indices of a single number when stored in a rectangular array.
//
//    Thus, our goal is, given the row I and column J of the data,
//    to produce the value K which indicates its position in the linear
//    array.
//
//    The triangular numbers are the indices associated with the
//    diagonal elements of the original array, T(0,0), T(1,1), T(2,2), 
//    T(3,3) and so on.
//
//    The formula is:
//
//      K = J + ( ( I * ( I + 1 ) ) / 2
//
//  Example:
//
//    I  J  K
//
//    0  0  0
//    1  0  1
//    1  1  2
//    2  0  3
//    2  1  4
//    2  2  5
//    3  0  6
//    3  1  7
//    3  2  8
//    3  3  9
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the row and column indices.  I and J must
//    be nonnegative, and J must not be greater than I.
//
//    Output, int TRIANGLE_TO_I4, the linear index of the (I,J) element.
//
{
  int value;

  if ( i < 0 )
  {
    cerr << "\n";
    cerr << "TRIANGLE_TO_I4 - Fatal error!\n";
    cerr << "  I < 0.\n";
    cerr << "  I = " << i << "\n";
    exit ( 1 );
  }
  else if ( j < 0 )
  {
    cerr << "\n";
    cerr << "TRIANGLE_TO_I4 - Fatal error!\n";
    cerr << "  J < 0.\n";
    cerr << "  J = " << j << "\n";
    exit ( 1 );
  }
  else if ( i < j )
  {
    cerr << "\n";
    cerr << "TRIANGLE_TO_I4 - Fatal error!\n";
    cerr << "  I < J.\n";
    cerr << "  I = " << i << "\n";
    cerr << "  J = " << j << "\n";
    exit ( 1 );
  }

  value = j + (  i * ( i + 1 ) ) / 2;

  return value;
}
//****************************************************************************80

int v_hofstadter ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    V_HOFSTADTER computes the Hofstadter V sequence.
//
//  Discussion:
//
//    V(N) = 0                          if N = 0
//         = 1                          if 1 <= N <= 4
//         = V(N-V(N-1)) + V(N-V(N-4)), otherwise.
//
//    V(N) is defined for all nonnegative integers.
//
//  Table:
//
//     N  V(N)
//    --  ----
//
//     0     0
//     1     1
//     2     1
//     3     1
//     4     1
//     5     2
//     6     3
//     7     4
//     8     5
//     9     5
//    10     6
//    11     6
//    12     7
//    13     8
//    14     8
//    15     9
//    16     9
//    17    10
//    18    11
//    19    11
//    20    11
//    21    12
//    22    12
//    23    13
//    24    14
//    25    14
//    26    15
//    27    15
//    28    16
//    29    17
//    30    17
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the function.
//
//    Output, int V_HOFSTADTER, the value of the function.
//
{
  if ( n <= 0 )
  {
    return 0;
  }
  else if ( n <= 4 )
  {
    return 1;
  }
  else
  {
    return (  v_hofstadter ( n - v_hofstadter(n-1) ) 
           + v_hofstadter ( n - v_hofstadter(n-4) ) );
  }

}
//****************************************************************************80

void vibonacci ( int n, int *seed, int v[] )

//****************************************************************************80
//
//  Purpose:
//
//    VIBONACCI computes the first N Vibonacci numbers.
//
//  Discussion:
//
//    The "Vibonacci numbers" are a generalization of the Fibonacci numbers:
//      V(N+1) = +/- V(N) +/- V(N-1)
//    where the signs are chosen randomly.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Brian Hayes,
//    The Vibonacci Numbers,
//    American Scientist,
//    July-August 1999, Volume 87, Number 4.
//
//    Divakar Viswanath,
//    Random Fibonacci sequences and the number 1.13198824,
//    Mathematics of Computation, 1998.
//
//  Parameters:
//
//    Input, int N, the highest number to compute.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int V(N), the first N Vibonacci numbers.  By convention,
//    V(1) and V(2) are taken to be 1.
//
{
  int i;
  int j;
  int s1;
  int s2;

  if ( n <= 0 )
  {
    return;
  }

  v[0] = 1;

  if ( n <= 1 )
  {
    return;
  }

  v[1] = 1;

  for ( i = 2; i < n; i++ )
  {
    j = i4_uniform ( 0, 1, seed );

    if ( j == 0 )
    {
      s1 = -1;
    }
    else
    {
      s1 = +1;
    }

    j = i4_uniform ( 0, 1, seed );

    if ( j == 0 )
    {
      s2 = -1;
    }
    else
    {
      s2 = +1;
    }

    v[i] = s1 * v[i-1] + s2 * v[i-2];

  }

  return;
}
//****************************************************************************80

void zeckendorf ( int n, int m_max, int *m, int i_list[], int f_list[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZECKENDORF produces the Zeckendorf decomposition of a positive integer.
//
//  Discussion:
//
//    Zeckendorf proved that every positive integer can be represented
//    uniquely as the sum of non-consecutive Fibonacci numbers.
//
//    N = sum ( 1 <= I <= M ) F_LIST(I)
//
//  Example:
//
//     N    Decomposition
//
//    50    34 + 13 + 3
//    51    34 + 13 + 3 + 1
//    52    34 + 13 + 5
//    53    34 + 13 + 5 + 1
//    54    34 + 13 + 5 + 2
//    55    55
//    56    55 + 1
//    57    55 + 2
//    58    55 + 3
//    59    55 + 3 + 1
//    60    55 + 5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the positive integer to be decomposed.
//
//    Input, int M_MAX, the maximum number of parts in the decomposition.
//    Set M_MAX = 100 to be safe.
//
//    Output, int M, the number of parts in the decomposition.
//
//    Output, int I_LIST[M], the index of the Fibonacci numbers
//    in the decomposition.
//
//    Output, int F_LIST[M], the value of the Fibonacci numbers
//    in the decomposition.
//
{
  int f;
  int i;
  int j;

  *m = 0;
//
//  Extract a sequence of Fibonacci numbers.
//
  while ( 0 < n && *m < m_max )
  {
    fibonacci_floor ( n, &f, &i );

    i_list[*m] = i;
    *m = *m + 1;
    n = n - f;
  }
//
//  Replace any pair of consecutive indices ( I, I-1 ) by I+1.
//
  for ( i = *m; 2 <= i; i-- )
  {
    if ( i_list[i-2] == i_list[i-1] + 1 )
    {
      i_list[i-2] = i_list[i-2] + 1;
      for ( j = i; j <= *m-1; j++ )
      {
        i_list[j-1] = i_list[j];
      }
      i_list[*m-1] = 0;
      *m = *m - 1;
    }

  }
//
//  Fill in the actual values of the Fibonacci numbers.
//
  for ( i = 0; i < *m; i++ )
  {
    f_list[i] = fibonacci_direct ( i_list[i] );
  }

  return;
}
//****************************************************************************80

double zernike_poly ( int m, int n, double rho )

//****************************************************************************80
//
//  Purpose:
//
//    ZERNIKE_POLY evaluates a Zernike polynomial at RHO.
//
//  Discussion:
//
//    This routine uses the facts that:
//
//    *) R^M_N = 0 if M < 0, or N < 0, or N < M.
//    *) R^M_M = RHO^M
//    *) R^M_N = 0 if mod ( N - M ) = 1.
//
//    and the recursion:
//
//    R^M_(N+2) = A * [ ( B * RHO^2 - C ) * R^M_N - D * R^M_(N-2) ]
//
//    where
//
//    A = ( N + 2 ) / ( ( N + 2 )^2 - M^2 )
//    B = 4 * ( N + 1 )
//    C = ( N + M )^2 / N + ( N - M + 2 )^2 / ( N + 2 )
//    D = ( N^2 - M^2 ) / N
//
//    I wish I could clean up the recursion in the code, but for
//    now, I have to treat the case M = 0 specially.
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
//  Reference:
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input, int M, the upper index.
//
//    Input, int N, the lower index.
//
//    Input, double RHO, the radial coordinate.
//
//    Output, double ZERNIKE_POLY, the value of the Zernike
//    polynomial R^M_N at the point RHO.
//
{
  double a;
  double b;
  double c;
  double d;
  int nn;
  double z;
  double zm2;
  double zp2;
//
//  Do checks.
//
  if ( m < 0 )
  {
    z = 0.0;
    return z;
  }

  if ( n < 0 )
  {
    z = 0.0;
    return z;
  }

  if ( n < m )
  {
    z = 0.0;
    return z;
  }

  if ( ( n - m ) % 2 == 1 )
  {
    z = 0.0;
    return z;
  }

  zm2 = 0.0;
  z = pow ( rho, m );

  if ( m == 0 )
  {
    if ( n == 0 )
    {
      return z;
    }

    zm2 = z;
    z = 2.0 * rho * rho - 1.0;

    for ( nn = m+2; nn <= n-2; nn = nn + 2 )
    {
      a = ( double ) ( nn + 2 ) 
        / ( double ) ( ( nn + 2 ) * ( nn + 2 ) - m * m );

      b = ( double ) ( 4 * ( nn + 1 ) );

      c = ( double ) ( ( nn + m ) * ( nn + m ) ) / ( double ) ( nn ) 
        + ( double ) ( ( nn - m + 2 ) * ( nn - m + 2 ) ) 
        / ( double ) ( nn + 2 );

      d = ( double ) ( nn * nn - m * m ) / ( double ) ( nn );

      zp2 = a * ( ( b * rho * rho - c ) * z - d * zm2 );
      zm2 = z;
      z = zp2;
    }
  }
  else
  {
    for ( nn = m; nn <= n-2; nn = nn + 2 )
    {
      a = ( double ) ( nn + 2 ) 
        / ( double ) ( ( nn + 2 ) * ( nn + 2 ) - m * m );

      b = ( double ) ( 4 * ( nn + 1 ) );

      c = ( double ) ( ( nn + m ) * ( nn + m ) ) / ( double ) ( nn ) 
        + ( double ) ( ( nn - m + 2 ) * ( nn - m + 2 ) ) 
        / ( double ) ( nn + 2 );

      d = ( double ) ( nn * nn - m * m ) / ( double ) ( nn );

      zp2 = a * ( ( b * rho * rho - c ) * z - d * zm2 );
      zm2 = z;
      z = zp2;
    }
  }

  return z;
}
//****************************************************************************80

double *zernike_poly_coef ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    ZERNIKE_POLY_COEF: coefficients of a Zernike polynomial.
//
//  Discussion:
//
//    With our coefficients stored in C(0:N), the
//    radial function R^M_N(RHO) is given by
//
//      R^M_N(RHO) = C(0)
//                 + C(1) * RHO
//                 + C(2) * RHO^2
//                 + ...
//                 + C(N) * RHO^N
//
//    and the odd and even Zernike polynomials are
//
//      Z^M_N(RHO,PHI,odd)  = R^M_N(RHO) * sin(PHI)
//      Z^M_N(RHO,PHI,even) = R^M_N(RHO) * cos(PHI)
//
//    The first few "interesting" values of R are:
//
//    R^0_0 = 1
//
//    R^1_1 = RHO
//
//    R^0_2 = 2 * RHO^2 - 1
//    R^2_2 =     RHO^2
//
//    R^1_3 = 3 * RHO^3 - 2 * RHO
//    R^3_3 =     RHO^3
//
//    R^0_4 = 6 * RHO^4 - 6 * RHO^2 + 1
//    R^2_4 = 4 * RHO^4 - 3 * RHO^2
//    R^4_4 =     RHO^4
//
//    R^1_5 = 10 * RHO^5 - 12 * RHO^3 + 3 * RHO
//    R^3_5 =  5 * RHO^5 -  4 * RHO^3
//    R^5_5 =      RHO^5
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
//  Reference:
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input, int M, N, the parameters of the polynomial.
//    Normally, 0 <= M <= N and 0 <= N.
//
//    Output, double ZERNIKE_POLY_COEF[N+1], the coefficients of the polynomial.
//
{
  double *c;
  int l;
  int nm_minus;
  int nm_plus;

  c = new double[n+1];

  r8vec_zero ( n+1, c );

  if ( n < 0 )
  {
    return c;
  }

  if ( m < 0 )
  {
    return c;
  }

  if ( n < m )
  {
    return c;
  }

  if ( ( n - m ) % 2 == 1 )
  {
    return c;
  }

  nm_plus = ( m + n ) / 2;
  nm_minus = ( n - m ) / 2;

  c[n] = r8_choose ( n, nm_plus );

  for ( l = 0; l <= nm_minus - 1; l++ )
  {
    c[n-2*l-2] = - ( double ) ( ( nm_plus - l ) * ( nm_minus - l ) ) 
      * c[n-2*l] / ( double ) ( ( n - l ) * ( l + 1 ) );

  }

  return c;
}
//****************************************************************************80

double zeta ( double p )

//****************************************************************************80
//
//  Purpose:
//
//    ZETA estimates the Riemann Zeta function.
//
//  Definition:
//
//    For 1 < P, the Riemann Zeta function is defined as:
//
//      ZETA ( P ) = Sum ( 1 <= N < +oo ) 1 / N^P
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, double P, the power to which the integers are raised.
//    P must be greater than 1.  
//
//    Output, double ZETA, an approximation to the Riemann
//    Zeta function.
//
{
  int n;
  double value;
  double value_old;

  if ( p <= 1.0 )
  {
    cerr << "\n";
    cerr << "ZETA - Fatal error!\n";
    cerr << "  Exponent P <= 1.0.\n";
    exit ( 1 );
  }

  value = 0.0;
  n = 0;

  for ( ; ; )
  {
    n = n + 1;
    value_old = value;
    value = value + 1.0 / pow ( ( double ) n, p );

    if ( value <= value_old || 1000 <= n )
    {
      break;
    }

  }

  return value;
}
//****************************************************************************80

void zeta_values ( int *n_data, int *n, double *zeta )

//****************************************************************************80
//
//  Purpose:
//
//    ZETA_VALUES returns some values of the Riemann Zeta function.
//
//  Discussion:
//
//    ZETA(N) = sum ( 1 <= I < +oo ) 1 / I**N
//
//    In Mathematica, the function can be evaluated by:
//
//      Zeta[n]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 August 2004
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
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *N, the argument of the Zeta function.
//
//    Output, double *ZETA, the value of the Zeta function.
//
{
# define N_MAX 15

  int n_vec[N_MAX] = { 
     2, 
     3, 
     4, 
     5, 
     6, 
     7, 
     8, 
     9, 
    10, 
    11, 
    12, 
    16, 
    20, 
    30, 
    40 };

  double zeta_vec[N_MAX] = { 
     0.164493406684822643647E+01,  
     0.120205690315959428540E+01,  
     0.108232323371113819152E+01,  
     0.103692775514336992633E+01,  
     0.101734306198444913971E+01,  
     0.100834927738192282684E+01,  
     0.100407735619794433939E+01,  
     0.100200839292608221442E+01,  
     0.100099457512781808534E+01,  
     0.100049418860411946456E+01,  
     0.100024608655330804830E+01,  
     0.100001528225940865187E+01,  
     0.100000095396203387280E+01,  
     0.100000000093132743242E+01,  
     0.100000000000090949478E+01 }; 

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *zeta = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *zeta = zeta_vec[*n_data-1];
  }

  return;
# undef N_MAX
}

