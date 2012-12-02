# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "prob.hpp"

//****************************************************************************80

double angle_cdf ( double x, int n )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_CDF evaluates the Angle CDF.
//
//  Modified:
//
//    19 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation and Sensitivity of Queueing Networks,
//    Wiley, 1986.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, int N, the spatial dimension.
//    N must be at least 2.
//
//    Output, real CDF, the value of the CDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;
  double zero = 0.0;

  if ( n < 2 )
  {
    cout << "\n";
    cout << "ANGLE_CDF - Fatal error!\n";
    cout << "  N must be at least 2.\n";
    cout << "  The input value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( x < 0.0 )
  {
    cdf = 0.0;
  }
  else if ( pi < x )
  {
    cdf = 1.0;
  }
  else if ( n == 2 )
  {
    cdf = x / pi;
  }
  else
  {
    cdf = sin_power_int ( zero, x, n-2 )
      * r8_gamma ( ( double ) ( n ) / 2.0 )
      / ( sqrt ( pi ) * r8_gamma ( ( double ) ( n - 1 ) / 2.0 ) );
  }

  return cdf;
}
//****************************************************************************80

double angle_mean ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_MEAN returns the mean of the Angle PDF.
//
//  Modified:
//
//    02 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//    N must be at least 2.
//
//    Output, double ANGLE_MEAN, the mean of the PDF.
//
{
  double mean;
  const double pi = 3.14159265358979323;

  mean = pi / 2.0;

  return mean;
}
//****************************************************************************80

double angle_pdf ( double x, int n )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_PDF evaluates the Angle PDF.
//
//  Discussion:
//
//    X is an angle between 0 and PI, corresponding to the angle
//    made in an N-dimensional space, between a fixed line passing
//    through the origin, and an arbitrary line that also passes
//    through the origin, which is specified by a choosing any point
//    on the N-dimensional sphere with uniform probability.
//
//  Formula:
//
//    PDF(X) = ( sin ( X ) )**(N-2) * Gamma ( N / 2 )
//             / ( sqrt ( PI ) * Gamma ( ( N - 1 ) / 2 ) )
//
//    PDF(X) = 1 / PI if N = 2.
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
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation and Sensitivity of Queueing Networks,
//    Wiley, 1986.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, int N, the spatial dimension.
//    N must be at least 2.
//
//    Output, real PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  if ( n < 2 )
  {
    cout << "\n";
    cout << "ANGLE_PDF - Fatal error!\n";
    cout << "  N must be at least 2.\n";
    cout << "  The input value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( x < 0.0 || pi < x )
  {
    pdf = 0.0;
  }
  else if ( n == 2 )
  {
    pdf = 1.0 / pi;
  }
  else
  {
    pdf = pow ( ( sin ( x ) ), ( n - 2 ) )
      * r8_gamma ( ( double ) ( n ) / 2.0 )
      / ( sqrt ( pi ) * r8_gamma ( ( double ) ( n - 1 ) / 2.0 ) );
  }

  return pdf;
}
//****************************************************************************80

double anglit_cdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLIT_CDF evaluates the Anglit CDF.
//
//  Modified:
//
//    29 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Output, double ANGLIT_CDF, the value of the CDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;

  if ( x <  - 0.25 * pi )
  {
    cdf = 0.0;
  }
  else if ( x < 0.25 * pi )
  {
    cdf = 0.5 - 0.5 * cos ( 2.0 * x + pi / 2.0 );
  }
  else
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double anglit_cdf_inv ( double cdf )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLIT_CDF_INV inverts the Anglit CDF.
//
//  Modified:
//
//    29 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Output, double ANGLIT_CDF_INV, the corresponding argument.
//
{
  const double pi = 3.14159265358979323;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "ANGLIT_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = 0.5 * ( acos ( 1.0 - 2.0 * cdf ) - pi / 2.0 );

  return x;
}
//****************************************************************************80

double anglit_mean ( )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLIT_MEAN returns the mean of the Anglit PDF.
//
//  Modified:
//
//    28 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ANGLIT_MEAN, the mean of the PDF.
//
{
  return 0.0;
}
//****************************************************************************80

double anglit_pdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLIT_PDF evaluates the Anglit PDF.
//
//  Formula:
//
//    PDF(X) = SIN ( 2 * X + PI / 2 ) for -PI/4 <= X <= PI/4
//
//  Modified:
//
//    28 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Output, double ANGLIT_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  if ( x <= - 0.25 * pi || 0.25 * pi <= x )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = sin ( 2.0 * x + 0.25 * pi );
  }

  return pdf;
}
//****************************************************************************80

double anglit_sample ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLIT_SAMPLE samples the Anglit PDF.
//
//  Modified:
//
//    28 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double *X, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = anglit_cdf_inv ( cdf );

  return x;
}
//****************************************************************************80

double anglit_variance ( )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLIT_VARIANCE returns the variance of the Anglit PDF.
//
//  Discussion:
//
//    Variance =
//      Integral ( -PI/4 <= X <= PI/4 ) X^2 * SIN ( 2 * X + PI / 2 )
//
//    Antiderivative =
//      0.5 * X * SIN ( 2 * X + PI / 2 )
//      + ( 0.25 - 0.5 * X^2 ) * COS ( 2 * X + PI / 2 )
//
//  Modified:
//
//    29 December 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double ANGLIT_VARIANCE, the variance of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double variance;

  variance = 0.0625 * pi * pi - 0.5;

  return variance;
}
//****************************************************************************80

double arcsin_cdf ( double x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_CDF evaluates the Arcsin CDF.
//
//  Modified:
//
//    20 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, the parameter of the CDF.
//    A must be positive.
//
//    Output, double ARCSIN_CDF, the value of the CDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;

  if ( x <= -a )
  {
    cdf = 0.0;
  }
  else if ( x < a )
  {
    cdf = 0.5 + asin ( x / a ) / pi;
  }
  else
  {
    cdf = 1.0;
  }

  return cdf;
}
//*****************************************************************************

double arcsin_cdf_inv ( double cdf, double a )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_CDF_INV inverts the Arcsin CDF.
//
//  Modified:
//
//    20 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, the parameter of the CDF.
//    A must be positive.
//
//    Output, double ARCSIN_CDF_INV, the corresponding argument.
//
{
  const double pi = 3.14159265358979323;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "ARCSIN_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a * sin ( pi * ( cdf - 0.5 ) );

  return x;
}
//****************************************************************************80

bool arcsin_check ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_CHECK checks the parameter of the Arcsin CDF.
//
//  Modified:
//
//    27 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A.
//
//    Output, bool ARCSIN_CHECK, is TRUE if the data is acceptable.
{
  if ( a <= 0.0 )
  {
    return false;
  }
  else
  {
    return true;
  }
}
//*****************************************************************************

double arcsin_mean ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_MEAN returns the mean of the Arcsin PDF.
//
//
//  Modified:
//
//    20 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the CDF.
//    A must be positive.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean = 0.0;

  return mean;
}
//*****************************************************************************

double arcsin_pdf ( double x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_PDF evaluates the Arcsin PDF.
//
//  Discussion:
//
//    The LOGISTIC EQUATION has the form:
//
//      X(N+1) = 4.0 * LAMBDA * ( 1.0 - X(N) ).
//
//    where 0 < LAMBDA <= 1.  This nonlinear difference equation maps
//    the unit interval into itself, and is a simple example of a system
//    exhibiting chaotic behavior.  Ulam and von Neumann studied the
//    logistic equation with LAMBDA = 1, and showed that iterates of the
//    function generated a sequence of pseudorandom numbers with
//    the Arcsin probability density function.
//
//    The derived sequence
//
//      Y(N) = ( 2 / PI ) * Arcsin ( SQRT ( X(N) ) )
//
//    is a pseudorandom sequence with the uniform probability density
//    function on [0,1].  For certain starting values, such as X(0) = 0, 0.75,
//    or 1.0, the sequence degenerates into a constant sequence, and for
//    values very near these, the sequence takes a while before becoming
//    chaotic.
//
//  Formula:
//
//    PDF(X) = 1 / ( PI * Sqrt ( A * A - X * X ) )
//
//  Reference:
//
//    Daniel Zwillinger, Stephen Kokoska,
//    CRC Standard Probability and Statistics Tables and Formulae,
//    Chapman and Hall/CRC, 2000, pages 114-115.
//
//  Modified:
//
//    20 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 < X < 1.0.
//
//    Input, double A, the parameter of the CDF.
//    A must be positive.
//
//    Output, double ARCSIN_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  if ( a <= 0 )
  {
    cout << "\n";
    cout << "ARCSIN_PDF - Fatal error!\n";
    cout << "  Parameter A <= 0.\n";
    exit ( 1 );
  }

  if ( x <= -a || a <= x )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = 1.0 / ( pi * sqrt ( a * a - x * x ) );
  }

  return pdf;
}
//*****************************************************************************

double arcsin_sample ( double a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_SAMPLE samples the Arcsin PDF.
//
//  Modified:
//
//    20 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the CDF.
//    A must be positive.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double ARCSIN_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = arcsin_cdf_inv ( cdf, a );

  return x;
}
//*****************************************************************************

double arcsin_variance ( double a )

//*****************************************************************************
//
//  Purpose:
//
//    ARCSIN_VARIANCE returns the variance of the Arcsin PDF.
//
//  Modified:
//
//    20 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the CDF.
//    A must be positive.
//
//    Output, double ARCSIN_VARIANCE, the variance of the PDF.
//
{
  double variance = a * a / 2.0;

  return variance;
}
//****************************************************************************80

double benford_pdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    BENFORD_PDF returns the Benford probability of one or more significant digits.
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
//  Formula:
//
//    PDF(X) = LOG10 ( ( X + 1 ) / X ).
//
//  Reference:
//
//    Frank Benford,
//    The Law of Anomalous Numbers,
//    Proceedings of the American Philosophical Society,
//    Volume 78, pages 551-572, 1938.
//
//    Ted Hill,
//    The First Digit Phenomenon,
//    American Scientist,
//    Volume 86, July/August 1998, pages 358 - 363.
//
//    R Raimi,
//    The Peculiar Distribution of First Digits,
//    Scientific American,
//    December 1969, pages 109-119.
//
//  Modified:
//
//    13 August 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the string of significant digits to be checked.
//    If X is 1, then we are asking for the Benford probability that
//    a value will have first digit 1.  If X is 123, we are asking for
//    the probability that the first three digits will be 123, and so on.
//
//    Output, double BENFORD_PDF, the Benford probability that an item taken
//    from a real world distribution will have the initial digits X.
//
{
  double pdf;

  if ( x <= 0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = log10 ( ( double ) ( x + 1 ) / ( double ) ( x ) );
  }

  return pdf;
}
//****************************************************************************80

double bernoulli_cdf ( int x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_CDF evaluates the Bernoulli CDF.
//
//
//  Modified:
//
//    22 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the number of successes on a single trial.
//    X = 0 or 1.
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Output, double BERNOULLI_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x < 0 )
  {
    cdf = 0.0;
  }
  else if ( x == 0 )
  {
    cdf = 1.0 - a;
  }
  else
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

int bernoulli_cdf_inv ( double cdf, double a )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_CDF_INV inverts the Bernoulli CDF.
//
//  Modified:
//
//    22 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, the parameter of the PDF.
//    0.0 <= A <= 1.0.
//
//    Output, int BERNOULLI_CDF_INV, the corresponding argument.
//
{
  int x;
//
  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "BERNOULLI_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf <= 1.0 - a )
  {
    x = 0;
  }
  else
  {
    x = 1;
  }

  return x;
}
//****************************************************************************80

bool bernoulli_check ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_CHECK checks the parameter of the Bernoulli CDF.
//
//  Modified:
//
//    22 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 <= A <= 1.0.
//
//    Output, bool BERNOULLI_CHECK, is TRUE if the data is acceptable.
{
  if ( a < 0.0 || 1.0 < a )
  {
    return false;
  }
  else
  {
    return true;
  }
}
//****************************************************************************80

double bernoulli_mean ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_MEAN returns the mean of the Bernoulli PDF.
//
//  Modified:
//
//    22 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of success.
//    0.0 <= A <= 1.0.
//
//    Output, double BERNOULLI_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double bernoulli_pdf ( int x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_PDF evaluates the Bernoulli PDF.
//
//  Formula:
//
//    PDF(X)(A) = A**X * ( 1.0 - A )**( X - 1 )
//
//    X = 0 or 1.
//
//  Discussion:
//
//    The Bernoulli PDF describes the simple case in which a single trial
//    is carried out, with two possible outcomes, called "success" and
//    "failure"; the probability of success is A.
//
//  Modified:
//
//    22 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the number of successes on a single trial.
//    X = 0 or 1.
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Output, double BERNOULLI_PDF, the value of the PDF.
//
{
  double pdf;
//
  if ( x < 0 )
  {
    pdf = 0.0;
  }
  else if ( x == 0 )
  {
    pdf = 1.0 - a;
  }
  else if ( x == 1 )
  {
    pdf = a;
  }
  else
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

int bernoulli_sample ( double a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_SAMPLE samples the Bernoulli PDF.
//
//  Modified:
//
//    22 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int BERNOULLI_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  int x;

  cdf = r8_uniform_01 ( seed );

  x = bernoulli_cdf_inv ( cdf, a );

  return x;
}
//****************************************************************************80

double bernoulli_variance ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_VARIANCE returns the variance of the Bernoulli PDF.
//
//  Modified:
//
//    22 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Output, double BERNOULLI_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = a * ( 1.0 - a );

  return variance;
}
//****************************************************************************80

double bessel_i0 ( double arg )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I0 evaluates the modified Bessel function I0.
//
//  Discussion:
//
//    The main computation evaluates slightly modified forms of
//    minimax approximations generated by Blair and Edwards, Chalk
//    River (Atomic Energy of Canada Limited) Report AECL-4928,
//    October, 1974.  This transportable program is patterned after
//    the machine dependent FUNPACK packet NATSI0, but cannot match
//    that version for efficiency or accuracy.  This version uses
//    rational functions that theoretically approximate I-SUB-0(X)
//    to at least 18 significant decimal digits.
//
//  Machine dependent constants:
//
//    beta   = Radix for the floating-point system
//    maxexp = Smallest power of beta that overflows
//    XSMALL = Positive argument such that 1.0 - X = 1.0 to
//             machine precision for all ABS(X) .LE. XSMALL.
//    XMAX =   Largest argument acceptable to BESI0;  Solution to
//             equation:
//               W(X) * (1+1/(8*X)+9/(128*X^2) = beta**maxexp
//             where  W(X) = EXP(X)/sqrt(2*PI*X)
//
//    Approximate values for some important machines are:
//
//                             beta       maxexp       XSMALL
//
//    CRAY-1        (S.P.)       2         8191       3.55D-15
//    Cyber 180/855
//      under NOS   (S.P.)       2         1070       3.55D-15
//    IEEE (IBM/XT,
//      SUN, etc.)  (S.P.)       2          128       2.98D-8
//    IEEE (IBM/XT,
//      SUN, etc.)  (D.P.)       2         1024       5.55D-17
//    IBM 3033      (D.P.)      16           63       6.95D-18
//    VAX           (S.P.)       2          127       2.98D-8
//    VAX D-Format  (D.P.)       2          127       6.95D-18
//    VAX G-Format  (D.P.)       2         1023       5.55D-17
//
//
//                                  XMAX
//
//    CRAY-1        (S.P.)       5682.810
//    Cyber 180/855
//      under NOS   (S.P.)       745.893
//    IEEE (IBM/XT,
//      SUN, etc.)  (S.P.)        91.900
//    IEEE (IBM/XT,
//      SUN, etc.)  (D.P.)       713.986
//    IBM 3033      (D.P.)       178.182
//    VAX           (S.P.)        91.203
//    VAX D-Format  (D.P.)        91.203
//    VAX G-Format  (D.P.)       713.293
//
//  Author:
//
//    W. J. Cody and L. Stoltz,
//    Mathematics and Computer Science Division,
//    Argonne National Laboratory,
//    Argonne, Illinois, 60439.
//
//    C++ translation by John Burkardt.
//
//  Parameters:
//
//    Input, double ARG, the argument.
//
//    Output, double BESSEL_I0, the value of the modified Bessel function
//    of the first kind.
//
{
  double a;
  double b;
  double exp40 = 2.353852668370199854E+17;
  int i;
  double p[15] = {
    -5.2487866627945699800E-18,
    -1.5982226675653184646E-14,
    -2.6843448573468483278E-11,
    -3.0517226450451067446E-08,
    -2.5172644670688975051E-05,
    -1.5453977791786851041E-02,
    -7.0935347449210549190,
    -2.4125195876041896775E+03,
    -5.9545626019847898221E+05,
    -1.0313066708737980747E+08,
    -1.1912746104985237192E+10,
    -8.4925101247114157499E+11,
    -3.2940087627407749166E+13,
    -5.5050369673018427753E+14,
    -2.2335582639474375249E+15 };
  double pp[8] = {
    -3.9843750000000000000E-01,
     2.9205384596336793945,
    -2.4708469169133954315,
     4.7914889422856814203E-01,
    -3.7384991926068969150E-03,
    -2.6801520353328635310E-03,
     9.9168777670983678974E-05,
    -2.1877128189032726730E-06 };
  double q[5] = {
    -3.7277560179962773046E+03,
     6.5158506418655165707E+06,
    -6.5626560740833869295E+09,
     3.7604188704092954661E+12,
    -9.7087946179594019126E+14 };
  double qq[7] = {
    -3.1446690275135491500E+01,
     8.5539563258012929600E+01,
    -6.0228002066743340583E+01,
     1.3982595353892851542E+01,
    -1.1151759188741312645,
     3.2547697594819615062E-02,
    -5.5194330231005480228E-04 };
  double rec15 = 6.6666666666666666666E-02;
  double sump;
  double sumq;
  double value;
  double x;
  double xmax = 91.9;
  double xsmall = 2.98E-08;
  double xx;

  x = r8_abs ( arg );

  if ( x < xsmall )
  {
    value = 1.0;
  }
  else if ( x < 15.0 )
  {
//
//  XSMALL <= ABS(ARG) < 15.0
//
    xx = x * x;
    sump = p[0];
    for ( i = 1; i < 15; i++ )
    {
      sump = sump * xx + p[i];
    }

    xx = xx - 225.0;
    sumq = ((((
           xx + q[0] )
         * xx + q[1] )
         * xx + q[2] )
         * xx + q[3] )
         * xx + q[4];

    value = sump / sumq;
  }
  else if ( 15.0 <= x )
  {
    if ( xmax < x )
    {
      value = r8_huge ( );
    }
    else
    {
//
//  15.0 <= ABS(ARG)
//
      xx = 1.0 / x - rec15;

      sump = ((((((
                  pp[0]
           * xx + pp[1] )
           * xx + pp[2] )
           * xx + pp[3] )
           * xx + pp[4] )
           * xx + pp[5] )
           * xx + pp[6] )
           * xx + pp[7];

      sumq = ((((((
             xx + qq[0] )
           * xx + qq[1] )
           * xx + qq[2] )
           * xx + qq[3] )
           * xx + qq[4] )
           * xx + qq[5] )
           * xx + qq[6];

      value = sump / sumq;
//
//  Calculation reformulated to avoid premature overflow.
//
      if ( x <= xmax - 15.0 )
      {
        a = exp ( x );
        b = 1.0;
      }
      else
      {
        a = exp ( x - 40.0 );
        b = exp40;
      }

      value = ( ( value * a - pp[0] * a ) / sqrt ( x ) ) * b;
    }
  }

  return value;
}
//****************************************************************************80

void bessel_i0_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I0_VALUES returns some values of the I0 Bessel function.
//
//  Discussion:
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    The modified Bessel function I0(Z) corresponds to N = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselI[0,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
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
//
# define N_MAX 20

  double fx_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1010025027795146E+01,
     0.1040401782229341E+01,
     0.1092045364317340E+01,
     0.1166514922869803E+01,
     0.1266065877752008E+01,
     0.1393725584134064E+01,
     0.1553395099731217E+01,
     0.1749980639738909E+01,
     0.1989559356618051E+01,
     0.2279585302336067E+01,
     0.3289839144050123E+01,
     0.4880792585865024E+01,
     0.7378203432225480E+01,
     0.1130192195213633E+02,
     0.1748117185560928E+02,
     0.2723987182360445E+02,
     0.6723440697647798E+02,
     0.4275641157218048E+03,
     0.2815716628466254E+04 };

  double x_vec[N_MAX] = {
     0.00,
     0.20,
     0.40,
     0.60,
     0.80,
     0.10E+01,
     0.12E+01,
     0.14E+01,
     0.16E+01,
     0.18E+01,
     0.20E+01,
     0.25E+01,
     0.30E+01,
     0.35E+01,
     0.40E+01,
     0.45E+01,
     0.50E+01,
     0.60E+01,
     0.80E+01,
     0.10E+02 };

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

double bessel_i1 ( double arg )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I1 evaluates the Bessel I function of order I.
//
//  Discussion:
//
//    The main computation evaluates slightly modified forms of
//    minimax approximations generated by Blair and Edwards.
//    This transportable program is patterned after the machine-dependent
//    FUNPACK packet NATSI1, but cannot match that version for efficiency
//    or accuracy.  This version uses rational functions that theoretically
//    approximate I-SUB-1(X) to at least 18 significant decimal digits.
//    The accuracy achieved depends on the arithmetic system, the compiler,
//    the intrinsic functions, and proper selection of the machine-dependent
//    constants.
//
//  Machine-dependent constants:
//
//    beta   = Radix for the floating-point system.
//    maxexp = Smallest power of beta that overflows.
//    XMAX =   Largest argument acceptable to BESI1;  Solution to
//             equation:
//               EXP(X) * (1-3/(8*X)) / SQRT(2*PI*X) = beta**maxexp
//
//
//    Approximate values for some important machines are:
//
//                            beta       maxexp    XMAX
//
//    CRAY-1        (S.P.)       2         8191    5682.810
//    Cyber 180/855
//      under NOS   (S.P.)       2         1070     745.894
//    IEEE (IBM/XT,
//      SUN, etc.)  (S.P.)       2          128      91.906
//    IEEE (IBM/XT,
//      SUN, etc.)  (D.P.)       2         1024     713.987
//    IBM 3033      (D.P.)      16           63     178.185
//    VAX           (S.P.)       2          127      91.209
//    VAX D-Format  (D.P.)       2          127      91.209
//    VAX G-Format  (D.P.)       2         1023     713.293
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2004
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz,
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Blair, Edwards,
//    Chalk River Report AECL-4928,
//    Atomic Energy of Canada, Limited,
//    October, 1974.
//
//  Parameters:
//
//    Input, double ARG, the argument.
//
//    Output, double BESSEL_I1, the value of the Bessel
//    I1 function.
//
{
  double a;
  double b;
  double exp40 = 2.353852668370199854E+17;
  double forty = 40.0;
  double half = 0.5;
  int j;
  double one = 1.0;
  double one5 = 15.0;
  double p[15] = {
    -1.9705291802535139930E-19,
    -6.5245515583151902910E-16,
    -1.1928788903603238754E-12,
    -1.4831904935994647675E-09,
    -1.3466829827635152875E-06,
    -9.1746443287817501309E-04,
    -4.7207090827310162436E-01,
    -1.8225946631657315931E+02,
    -5.1894091982308017540E+04,
    -1.0588550724769347106E+07,
    -1.4828267606612366099E+09,
    -1.3357437682275493024E+11,
    -6.9876779648010090070E+12,
    -1.7732037840791591320E+14,
    -1.4577180278143463643E+15 };
  double pbar = 3.98437500E-01;
  double pp[8] = {
    -6.0437159056137600000E-02,
     4.5748122901933459000E-01,
    -4.2843766903304806403E-01,
     9.7356000150886612134E-02,
    -3.2457723974465568321E-03,
    -3.6395264712121795296E-04,
     1.6258661867440836395E-05,
    -3.6347578404608223492E-07 };
  double q[5] = {
    -4.0076864679904189921E+03,
     7.4810580356655069138E+06,
    -8.0059518998619764991E+09,
     4.8544714258273622913E+12,
    -1.3218168307321442305E+15 };
  double qq[6] = {
    -3.8806586721556593450,
     3.2593714889036996297,
    -8.5017476463217924408E-01,
     7.4212010813186530069E-02,
    -2.2835624489492512649E-03,
     3.7510433111922824643E-05 };
  double rec15 = 6.6666666666666666666E-02;
  double sump;
  double sumq;
  double two25 = 225.0;
  double value;
  double x;
  double xmax = 713.987;
  double xx;
  double zero = 0.0;
//
  x = r8_abs ( arg );
//
//  ABS(ARG) < EPSILON ( ARG )
//
  if ( x < r8_epsilon ( ) )
  {
    value = half * x;
  }
//
//  EPSILON ( ARG ) <= ABS(ARG) < 15.0
//
  else if ( x < one5 )
  {
    xx = x * x;
    sump = p[0];
    for ( j = 1; j < 15; j++ )
    {
      sump = sump * xx + p[j];
    }

    xx = xx - two25;

    sumq = ((((
          xx + q[0]
      ) * xx + q[1]
      ) * xx + q[2]
      ) * xx + q[3]
      ) * xx + q[4];

    value = ( sump / sumq ) * x;
  }
  else if ( xmax < x )
  {
    value = r8_huge ( );
  }
//
//  15.0 <= ABS(ARG)
//
  else
  {
    xx = one / x - rec15;

    sump = ((((((
               pp[0]
        * xx + pp[1]
      ) * xx + pp[2]
      ) * xx + pp[3]
      ) * xx + pp[4]
      ) * xx + pp[5]
      ) * xx + pp[6]
      ) * xx + pp[7];

    sumq = (((((
          xx + qq[0]
      ) * xx + qq[1]
      ) * xx + qq[2]
      ) * xx + qq[3]
      ) * xx + qq[4]
      ) * xx + qq[5];

    value = sump / sumq;

    if ( xmax - one5 < x )
    {
      a = exp ( x - forty );
      b = exp40;
    }
    else
    {
      a = exp ( x );
      b = one;
    }
    value = ( ( value * a + pbar * a ) / sqrt ( x ) ) * b;
  }

  if ( arg < zero )
  {
    value = -value;
  }

  return value;
}
//****************************************************************************80

void bessel_i1_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I1_VALUES returns some values of the I1 Bessel function.
//
//  Discussion:
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselI[1,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
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
# define N_MAX 20

  double fx_vec[N_MAX] = {
     0.0000000000000000,
     0.1005008340281251,
     0.2040267557335706,
     0.3137040256049221,
     0.4328648026206398,
     0.5651591039924850,
     0.7146779415526431,
     0.8860919814143274,
     0.1084810635129880E+01,
     0.1317167230391899E+01,
     0.1590636854637329E+01,
     0.2516716245288698E+01,
     0.3953370217402609E+01,
     0.6205834922258365E+01,
     0.9759465153704450E+01,
     0.1538922275373592E+02,
     0.2433564214245053E+02,
     0.6134193677764024E+02,
     0.3998731367825601E+03,
     0.2670988303701255E+04 };

  double x_vec[N_MAX] = {
     0.00,
     0.20,
     0.40,
     0.60,
     0.80,
     0.10E+01,
     0.12E+01,
     0.14E+01,
     0.16E+01,
     0.18E+01,
     0.20E+01,
     0.25E+01,
     0.30E+01,
     0.35E+01,
     0.40E+01,
     0.45E+01,
     0.50E+01,
     0.60E+01,
     0.80E+01,
     0.10E+02 };

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

void bessel_ix_values ( int *n_data, double *nu, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_IX_VALUES returns some values of the Ix Bessel function.
//
//  Discussion:
//
//    This set of data considers the less common case in which the
//    index of the Bessel function In is actually not an integer.
//    We may suggest this case by occasionally replacing the symbol
//    "In" by "Ix".
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselI[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2007
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
//    Output, double *NU, the order of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 28

  double fx_vec[N_MAX] = {
    0.3592084175833614E+00,
    0.9376748882454876E+00,
    2.046236863089055E+00,
    3.053093538196718E+00,
    4.614822903407601E+00,
    26.47754749755907E+00,
    2778.784603874571E+00,
    4.327974627242893E+07,
    0.2935253263474798E+00,
    1.099473188633110E+00,
    21.18444226479414E+00,
    2500.906154942118E+00,
    2.866653715931464E+20,
    0.05709890920304825E+00,
    0.3970270801393905E+00,
    13.76688213868258E+00,
    2028.512757391936E+00,
    2.753157630035402E+20,
    0.4139416015642352E+00,
    1.340196758982897E+00,
    22.85715510364670E+00,
    2593.006763432002E+00,
    2.886630075077766E+20,
    0.03590910483251082E+00,
    0.2931108636266483E+00,
    11.99397010023068E+00,
    1894.575731562383E+00,
    2.716911375760483E+20 };

  double nu_vec[N_MAX] = {
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00 };

  double x_vec[N_MAX] = {
      0.2E+00,
      1.0E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      5.0E+00,
     10.0E+00,
     20.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *nu = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *nu = nu_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
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
//    BETA(X,Y) = ( GAMMA(X) * GAMMA(Y) ) / GAMMA(X+Y)
//
//    BETA(X,Y) = BETA(Y,X).
//    BETA(X,Y) = Integral ( 0 <= T <= 1 ) T**(X-1) (1-T)**(Y-1) dT.
//    BETA(X,Y) = GAMMA(X) * GAMMA(Y) / GAMMA(X+Y)
//
//    Both X and Y must be greater than 0.
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
  if ( x <= 0.0 || y <= 0.0 )
  {
    cout << "\n";
    cout << "BETA - Fatal error!\n";
    cout << "  Both X and Y must be greater than 0.\n";
    exit ( 1 );
  }

  return ( exp ( gamma_log ( x ) + gamma_log ( y ) - gamma_log ( x + y ) ) );
}
//****************************************************************************80

double beta_binomial_cdf ( int x, double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_BINOMIAL_CDF evaluates the Beta Binomial CDF.
//
//  Discussion:
//
//    A simple summing approach is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the argument of the CDF.
//
//    Input, double A, B, parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input, int C, a parameter of the PDF.
//    0 <= C.
//
//    Output, double BETA_BINOMIAL_CDF, the value of the CDF.
//
{
  double cdf;
  double pdf;
  int y;

  if ( x < 0 )
  {
    cdf = 0.0;
  }
  else if ( x < c )
  {
    cdf = 0.0;
    for ( y = 0; y <= x; y++ )
    {
      pdf = beta ( a + ( double ) y, b + ( double ) ( c - y) )
        / ( ( double ) ( c + 1 )
         * beta ( ( double ) ( y + 1 ), ( double ) ( c - y + 1 ) )
         * beta ( a, b ) );
      cdf = cdf + pdf;
    }
  }
  else if ( c <= x )
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

int beta_binomial_cdf_inv ( double cdf, double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_BINOMIAL_CDF_INV inverts the Beta Binomial CDF.
//
//  Discussion:
//
//    A simple discrete approach is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input, int C, a parameter of the PDF.
//    0 <= C.
//
//    Output, int BETA_BINOMIAL_CDF_INV, the smallest X whose cumulative
//    density function is greater than or equal to CDF.
//
{
  double cum;
  double pdf;
  int x;
  int y;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "BETA_BINOMIAL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  cum = 0.0;

  for ( y = 0; y <= c; y++ )
  {
    pdf = beta ( a + ( double ) ( y ),
      b + ( double ) ( c - y ) ) / ( ( double ) ( c + 1 )
      * beta ( ( double ) ( y + 1 ),
      ( double ) ( c - y + 1 ) ) * beta ( a, b ) );

    cum = cum + pdf;

    if ( cdf <= cum )
    {
      x = y;
      return x;
    }

  }

  x = c;

  return x;
}
//****************************************************************************80

bool beta_binomial_check ( double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_BINOMIAL_CHECK checks the parameters of the Beta Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input, int C, a parameter of the PDF.
//    0 <= C.
//
//    Output, bool BETA_BINOMIAL_CHECK, is TRUE if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << "\n";
    cout << "BETA_BINOMIAL_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "BETA_BINOMIAL_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  if ( c < 0 )
  {
    cout << "\n";
    cout << "BETA_BINOMIAL_CHECK - Fatal error!\n";
    cout << "  C < 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double beta_binomial_mean ( double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_BINOMIAL_MEAN returns the mean of the Beta Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input, int C, a parameter of the PDF.
//    0 <= N.
//
//    Output, double BETA_BINOMIAL_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = ( double ) ( c ) * a / ( a + b );

  return mean;
}
//****************************************************************************80

double beta_binomial_pdf ( int x, double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_BINOMIAL_PDF evaluates the Beta Binomial PDF.
//
//  Discussion:
//
//    PDF(X)(A,B,C) = Beta(A+X,B+C-X)
//      / ( (C+1) * Beta(X+1,C-X+1) * Beta(A,B) )  for 0 <= X <= C.
//
//    This PDF can be reformulated as:
//
//      The beta binomial probability density function for X successes
//      out of N trials is
//
//      PDF2(X)( N, MU, THETA ) =
//        C(N,X) * Product ( 0 <= R <= X - 1 ) ( MU + R * THETA )
//               * Product ( 0 <= R <= N - X - 1 ) ( 1 - MU + R * THETA )
//               / Product ( 0 <= R <= N - 1 )  ( 1 + R * THETA )
//
//      where
//
//        C(N,X) is the combinatorial coefficient;
//        MU is the expectation of the underlying Beta distribution;
//        THETA is a shape parameter.
//
//      A THETA value of 0 ( or A+B --> Infinity ) results in the binomial
//      distribution:
//
//        PDF2(X) ( N, MU, 0 ) = C(N,X) * MU**X * ( 1 - MU )**(N-X)
//
//    Given A, B, C for PDF, then the equivalent PDF2 has:
//
//      N     = C
//      MU    = A / ( A + B )
//      THETA = 1 / ( A + B )
//
//    Given N, MU, THETA for PDF2, the equivalent PDF has:
//
//      A = MU / THETA
//      B = ( 1 - MU ) / THETA
//      C = N
//
//    BETA_BINOMIAL_PDF(X)(1,1,C) = UNIFORM_DISCRETE_PDF(X)(0,C-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the argument of the PDF.
//
//    Input, double A, B, parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input, int C, a parameter of the PDF.
//    0 <= C.
//
//    Output, double BETA_BINOMIAL_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 0 )
  {

    pdf = 0.0;
  }
  else if ( x <= c )
  {
    pdf = beta ( a + ( double ) ( x ), b + ( double ) ( c - x ) )
      / ( ( double ) ( c + 1 )
      * beta ( ( double ) ( x + 1 ),
      ( double ) ( c - x + 1 ) ) * beta ( a, b ) );
  }
  else if ( c < x )
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

int beta_binomial_sample ( double a, double b, int c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_BINOMIAL_SAMPLE samples the Beta Binomial CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input, int C, a parameter of the PDF.
//    0 <= C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int X, a sample of the PDF.
//
{
  double cdf;
  int x;

  cdf = r8_uniform_01 ( seed );

  x = beta_binomial_cdf_inv ( cdf, a, b, c );

  return x;
}
//****************************************************************************80

double beta_binomial_variance ( double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_BINOMIAL_VARIANCE returns the variance of the Beta Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input, int C, a parameter of the PDF.
//    0 <= C.
//
//    Output, double BETA_BINOMIAL_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( ( double ) ( c ) * a * b )
    * ( a + b + ( double ) ( c ) )
    / ( ( a + b ) * ( a + b ) * ( a + b + 1.0 ) );

  return variance;
}
//****************************************************************************80

double beta_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_CDF evaluates the Beta CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double BETA_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else if ( x <= 1.0 )
  {
    cdf = beta_inc ( a, b, x );
  }
  else
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double beta_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_CDF_INV inverts the Beta CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2004
//
//  Author:
//
//    Abernathy and Smith.
//
//    C++ translation by John Burkardt.
//
//  Reference:
//
//    Abernathy, Smith,
//    Algorithm 724,
//    ACM Transactions on Mathematical Software,
//    Volume 19, Number 4, December 1993, pages 481-483.
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double BETA_CDF_INV, the argument of the CDF.
//
{
# define MAXK 20

  double bcoeff;
  double cdf_x;
  double d[MAXK * (MAXK-1)];
  double error = 0.0001;
  double errapp = 0.01;
  int i;
  int j;
  int k;
  int loopct;
  double pdf_x;
  double q;
  double s1;
  double s2;
  double sum2;
  double t;
  double tail;
  double x;
  double xold;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "BETA_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }
//
//  Estimate the solution.
//
  x = a / ( a + b );

  xold = 0.0;
  loopct = 2;

  while ( errapp <= r8_abs ( ( x - xold ) / x ) && loopct != 0 )
  {
    xold = x;
    loopct = loopct - 1;
//
//  CDF_X = PROB { BETA(A,B) <= X }.
//  Q = ( CDF - CDF_X ) / PDF_X.
//
    cdf_x = beta_cdf ( x, a, b );

    pdf_x = beta_pdf ( x, a, b );

    q = ( cdf - cdf_x ) / pdf_x;
//
//  D(N,K) = C(N,K) * Q**(N+K-1) / (N-1)!
//
    t = 1.0 - x;
    s1 = q * ( b - 1.0 ) / t;
    s2 = q * ( 1.0 - a ) / x;
    d[2-1+0*MAXK] = s1 + s2;
    tail = d[2-1+0*MAXK] * q / 2.0;
    x = x + q + tail;

    k = 3;

    while ( error < r8_abs ( tail / x ) && k <= MAXK )
    {
//
//  Find D(2,K-2).
//
      s1 = q * ( ( double ) ( k ) - 2.0 ) * s1 / t;
      s2 = q * ( 2.0 - ( double ) ( k ) ) * s2 / x;
      d[2-1+(k-2)*MAXK] = s1 + s2;
//
//  Find D(3,K-3), D(4,K-4), D(5,K-5), ... , D(K-1,1).
//
      for ( i = 3; i <= k-1; i++ )
      {
        sum2 = d[2-1+0*MAXK] * d[i-2+(k-i)*MAXK];
        bcoeff = 1.0;

        for ( j = 1; j <= k - i; j++ )
        {
          bcoeff = ( bcoeff * ( double ) ( k - i - j + 1 ) )
            / ( double ) ( j );
          sum2 = sum2 + bcoeff * d[2-1+j*MAXK] * d[i-2+(k-i-j)*MAXK];
        }
        d[i-1+(k-i)*MAXK] = sum2 + d[i-2+(k-i+1)*MAXK] / ( double ) ( i - 1 );
      }
//
//  Compute D(K,0) and use it to expand the series.
//
      d[k-1+0*MAXK] = d[2-1+0*MAXK] * d[k-2+0*MAXK] + d[k-2+1*MAXK]
        / ( double ) ( k - 1 );
      tail = d[k-1+0*MAXK] * q / ( double ) ( k );
      x = x + tail;
//
//  Check for divergence.
//
      if ( x <= 0.0 || 1.0 <= x )
      {
        cout << " \n";
        cout << "BETA_CDF_INV - Fatal error!\n";
        cout << "  The series has diverged.\n";
        x = -1.0;
        return x;
      }
      k = k + 1;
    }
  }
  return x;
# undef MAXK
}
//****************************************************************************80

bool beta_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_CHECK checks the parameters of the Beta PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, bool BETA_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << "\n";
    cout << "BETA_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "BETA_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

void beta_cdf_values ( int *n_data, double *a, double *b, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_CDF_VALUES returns some values of the Beta CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = BetaDistribution [ a, b ]
//      CDF [ dist, x ]
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
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
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
//    Output, double *A, B, the parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double a_vec[N_MAX] = {
      0.10E+01,
      0.10E+01,
      0.10E+01,
      0.10E+01,
      0.10E+01,
      0.10E+01,
      0.10E+01,
      0.10E+01,
      0.20E+01,
      0.30E+01,
      0.40E+01,
      0.50E+01 };

  double b_vec[N_MAX] = {
      0.50,
      0.50,
      0.50,
      0.50,
      0.20E+01,
      0.30E+01,
      0.40E+01,
      0.50E+01,
      0.20E+01,
      0.20E+01,
      0.20E+01,
      0.20E+01 };

  double fx_vec[N_MAX] = {
     0.5131670194948620E-01,
     0.1055728090000841,
     0.1633399734659245,
     0.2254033307585166,
     0.3600000000000000,
     0.4880000000000000,
     0.5904000000000000,
     0.6723200000000000,
     0.2160000000000000,
     0.8370000000000000E-01,
     0.3078000000000000E-01,
     0.1093500000000000E-01 };

  double x_vec[N_MAX] = {
     0.10,
     0.20,
     0.30,
     0.40,
     0.20,
     0.20,
     0.20,
     0.20,
     0.30,
     0.30,
     0.30,
     0.30 };

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
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double beta_inc ( double a, double b, double x )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC returns the value of the incomplete Beta function.
//
//  Discussion:
//
//    This calculation requires an iteration.  In some cases, the iteration
//    may not converge rapidly, or may become inaccurate.
//
//    BETA_INC(A,B,X)
//
//      =   Integral ( 0 <= T <= X ) T**(A-1) (1-T)**(B-1) dT
//        / Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT
//
//      =   Integral ( 0 <= T <= X ) T**(A-1) (1-T)**(B-1) dT
//        / BETA(A,B)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    Original FORTRAN77 version by Majumder, Bhattacharjee.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Majumder, Bhattacharjee,
//    Algorithm AS63,
//    Applied Statistics,
//    1973, volume 22, number 3.
//
//  Parameters:
//
//    Input, double A, B, the parameters of the function.
//    0.0 < A,
//    0.0 < B.
//
//    Input, double X, the argument of the function.
//    Normally, 0.0 <= X <= 1.0.
//
//    Output, double BETA_INC, the value of the function.
//
{
  double cx;
  int i;
  int it;
  int it_max = 1000;
  bool indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double tol = 1.0E-07;
  double value;
  double xx;

  if ( a <= 0.0 )
  {
    cout << "\n";
    cout << "BETA_INC - Fatal error!\n";
    cout << "  A <= 0.\n";
    exit ( 1 );
  }

  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "BETA_INC - Fatal error!\n";
    cout << "  B <= 0.\n";
    exit ( 1 );
  }

  if ( x <= 0.0 )
  {
    value = 0.0;
    return value;
  }
  else if ( 1.0 <= x )
  {
    value = 1.0;
    return value;
  }
//
//  Change tail if necessary and determine S.
//
  psq = a + b;

  if ( a < ( a + b ) * x )
  {
    xx = 1.0 - x;
    cx = x;
    pp = b;
    qq = a;
    indx = true;
  }
  else
  {
    xx = x;
    cx = 1.0 - x;
    pp = a;
    qq = b;
    indx = false;
  }

  term = 1.0;
  i = 1;
  value = 1.0;

  ns = ( int ) ( qq + cx * ( a + b ) );
//
//  Use Soper's reduction formulas.
//
  rx = xx / cx;

  temp = qq - ( double ) i;
  if ( ns == 0 )
  {
    rx = xx;
  }

  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    if ( it_max < it )
    {
      cout << "\n";
      cout << "BETA_INC - Fatal error!\n";
      cout << "  Maximum number of iterations exceeded!\n";
      cout << "  IT_MAX = " << it_max << "\n";
      exit ( 1 );
    }

    term = term * temp * rx / ( pp + ( double ) ( i ) );
    value = value + term;
    temp = r8_abs ( term );

    if ( temp <= tol && temp <= tol * value )
    {
      break;
    }

    i = i + 1;
    ns = ns - 1;

    if ( 0 <= ns )
    {
      temp = qq - ( double ) i;
      if ( ns == 0 )
      {
        rx = xx;
      }
    }
    else
    {
      temp = psq;
      psq = psq + 1.0;
    }
  }
//
//  Finish calculation.
//
  value = value * exp ( pp * log ( xx )
    + ( qq - 1.0 ) * log ( cx ) ) / ( beta ( a, b ) * pp );

  if ( indx )
  {
    value = 1.0 - value;
  }

  return value;
}
//****************************************************************************80

void beta_inc_values ( int *n_data, double *a, double *b, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC_VALUES returns some values of the incomplete Beta function.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
//                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
//
//    Thus,
//
//      BETA_INC(A,B,0.0) = 0.0;
//      BETA_INC(A,B,1.0) = 1.0
//
//    In Mathematica, the function can be evaluated by:
//
//      BETA[X,A,B] / BETA[A,B]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2004
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
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
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
//    Output, double *A, B, the parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 30

  double a_vec[N_MAX] = {
      0.5,
      0.5,
      0.5,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      5.5,
     10.0,
     10.0,
     10.0,
     10.0,
     20.0,
     20.0,
     20.0,
     20.0,
     20.0,
     30.0,
     30.0,
     40.0 };

  double b_vec[N_MAX] = {
      0.5,
      0.5,
      0.5,
      0.5,
      0.5,
      0.5,
      0.5,
      1.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      2.0,
      5.0,
      0.5,
      5.0,
      5.0,
     10.0,
      5.0,
     10.0,
     10.0,
     20.0,
     20.0,
     10.0,
     10.0,
     20.0 };

  double fx_vec[N_MAX] = {
     0.6376856085851985E-01,
     0.2048327646991335,
     0.1000000000000000E+01,
     0.0000000000000000,
     0.5012562893380045E-02,
     0.5131670194948620E-01,
     0.2928932188134525,
     0.5000000000000000,
     0.2800000000000000E-01,
     0.1040000000000000,
     0.2160000000000000,
     0.3520000000000000,
     0.5000000000000000,
     0.6480000000000000,
     0.7840000000000000,
     0.8960000000000000,
     0.9720000000000000,
     0.4361908850559777,
     0.1516409096347099,
     0.8978271484375000E-01,
     0.1000000000000000E+01,
     0.5000000000000000,
     0.4598773297575791,
     0.2146816102371739,
     0.9507364826957875,
     0.5000000000000000,
     0.8979413687105918,
     0.2241297491808366,
     0.7586405487192086,
     0.7001783247477069 };

  double x_vec[N_MAX] = {
     0.01,
     0.10,
     1.00,
     0.00,
     0.01,
     0.10,
     0.50,
     0.50,
     0.10,
     0.20,
     0.30,
     0.40,
     0.50,
     0.60,
     0.70,
     0.80,
     0.90,
     0.50,
     0.90,
     0.50,
     1.00,
     0.50,
     0.80,
     0.60,
     0.80,
     0.50,
     0.60,
     0.70,
     0.80,
     0.70 };

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
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double beta_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_MEAN returns the mean of the Beta PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double BETA_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a / ( a + b );

  return mean;
}
//****************************************************************************80

double beta_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_PDF evaluates the Beta PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = X**(A-1) * (1-X)**(B-1) / BETA(A,B).
//
//    A = B = 1 yields the Uniform distribution on [0,1].
//    A = B = 1/2 yields the Arcsin distribution.
//        B = 1 yields the power function distribution.
//    A = B -> Infinity tends to the Normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double BETA_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 0.0 || 1.0 < x )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = pow ( x, ( a - 1.0 ) ) * pow ( ( 1.0 - x ), ( b - 1.0 ) ) / beta ( a, b );
  }

  return pdf;
}
//****************************************************************************80

double beta_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_SAMPLE samples the Beta PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Kennedy, James Gentle,
//    Algorithm BN,
//    Statistical Computing,
//    Dekker, 1980.
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input/output, int SEED, a seed for the random number generator.
//
//    Output, double BETA_SAMPLE, a sample of the PDF.
//
{
  double mu;
  double stdev;
  double test;
  double u;
  double x;
  double y;

  mu = ( a - 1.0 ) / ( a + b - 2.0 );
  stdev = 0.5 / sqrt ( a + b - 2.0 );

  for ( ; ; )
  {
    y = normal_01_sample ( seed );

    x = mu + stdev * y;

    if ( x < 0.0 || 1.0 < x )
    {
      continue;
    }

    u = r8_uniform_01 ( seed );

    test =     ( a - 1.0 )     * log (         x   / ( a - 1.0 ) )
             + ( b - 1.0 )     * log ( ( 1.0 - x ) / ( b - 1.0 ) )
             + ( a + b - 2.0 ) * log ( a + b - 2.0 ) + 0.5 * y * y;

    if ( log ( u ) <= test )
    {
      break;
    }

  }

  return x;
}
//****************************************************************************80

double beta_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_VARIANCE returns the variance of the Beta PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double BETA_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( a * b ) / ( ( a + b ) * ( a + b ) * ( 1.0 + a + b ) );

  return variance;
}
//****************************************************************************80

double binomial_cdf ( double x, int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_CDF evaluates the Binomial CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    A sequence of trials with fixed probability of success on
//    any trial is known as a sequence of Bernoulli trials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the desired number of successes.
//    0 <= X <= A.
//
//    Input, int A, the number of trials.
//    1 <= A.
//
//    Input, double B, the probability of success on one trial.
//    0.0 <= B <= 1.0.
//
//    Output, double CDF, the value of the CDF.
//
{
  int cnk;
  double cdf;
  int j;
  double pr;

  if ( x < 0 )
  {
    cdf = 0.0;
  }
  else if ( a <= x )
  {
    cdf = 1.0;
  }
  else if ( b == 0.0 )
  {
    cdf = 1.0;
  }
  else if ( b == 1.0 )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 0.0;

    for ( j = 0; j <= x; j++ )
    {
      cnk = binomial_coef ( a, j );

      pr = ( double ) ( cnk ) * pow ( b, j ) * pow ( ( 1.0 - b ), ( a - j ) );

      cdf = cdf + pr;

    }

  }

  return cdf;
}
//****************************************************************************80

void binomial_cdf_values ( int *n_data, int *a, double *b, int *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = BinomialDistribution [ n, p ]
//      CDF [ dist, x ]
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
//    Output, int *A, a parameter of the function.
//
//    Output, double *B, a parameter of the function.
//
//    Output, int *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 17

  int a_vec[N_MAX] = {
     2,  2,  2,  2,
     2,  4,  4,  4,
     4, 10, 10, 10,
    10, 10, 10, 10,
    10 };

  double b_vec[N_MAX] = {
     0.05,
     0.05,
     0.05,
     0.50,
     0.50,
     0.25,
     0.25,
     0.25,
     0.25,
     0.05,
     0.10,
     0.15,
     0.20,
     0.25,
     0.30,
     0.40,
     0.50  };

  double fx_vec[N_MAX] = {
     0.9025000000000000,
     0.9975000000000000,
     0.1000000000000000E+01,
     0.2500000000000000,
     0.7500000000000000,
     0.3164062500000000,
     0.7382812500000000,
     0.9492187500000000,
     0.9960937500000000,
     0.9999363101685547,
     0.9983650626000000,
     0.9901259090013672,
     0.9672065024000000,
     0.9218730926513672,
     0.8497316674000000,
     0.6331032576000000,
     0.3769531250000000 };

  int x_vec[N_MAX] = {
     0, 1, 2, 0,
     1, 0, 1, 2,
     3, 4, 4, 4,
     4, 4, 4, 4,
     4 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0.0;
    *x = 0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int binomial_cdf_inv ( double cdf, int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_CDF_INV inverts the Binomial CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, int A, the number of trials.
//    1 <= A.
//
//    Input, double B, the probability of success on one trial.
//    0.0 <= B <= 1.0.
//
//    Output, int BINOMIAL_CDF_INV, the corresponding argument.
//
{
  double cdf2;
  double pdf;
  int x;
  int x2;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "BINOMIAL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  cdf2 = 0.0;

  for ( x2 = 0; x2 <= a; x2++ )
  {
    pdf = binomial_pdf ( x2, a, b );

    cdf2 = cdf2 + pdf;

    if ( cdf <= cdf2 )
    {
      x = x2;
      return x;
    }

  }

  return x;
}
//****************************************************************************80

bool binomial_check ( int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_CHECK checks the parameter of the Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the number of trials.
//    1 <= A.
//
//    Input, double B, the probability of success on one trial.
//    0.0 <= B <= 1.0.
//
//    Output, bool BINOMIAL_CHECK, is true if the parameters are legal.
//
{
  if ( a < 1 )
  {
    cout << "\n";
    cout << "BINOMIAL_CHECK - Fatal error!\n";
    cout << "  A < 1.\n";
    return false;
  }

  if ( b < 0.0 || 1.0 < b )
  {
    cout << "\n";
    cout << "BINOMIAL_CHECK - Fatal error!\n";
    cout << "  B < 0 or 1 < B.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

int binomial_coef ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_COEF computes the Binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
//
//    CNK = C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    M L Wolfson, H V Wright,
//    Combinatorial of M Things Taken N at a Time,
//    ACM Algorithm 160,
//    Communications of the ACM,
//    April, 1963.
//
//  Parameters:
//
//    Input, int N, K, are the values of N and K.
//
//    Output, int CNK, the number of combinations of N
//    things taken K at a time.
//
{
  int cnk;
  int i;
  int mn;
  int mx;

  mn = i4_min ( k, n-k );

  if ( mn < 0 )
  {
    cnk = 0;
  }
  else if ( mn == 0 )
  {
    cnk = 1;
  }
  else
  {
    mx = i4_max ( k, n-k );
    cnk = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      cnk = ( cnk * ( mx + i ) ) / i;
    }
  }

  return cnk;
}
//****************************************************************************80

double binomial_coef_log ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_COEF_LOG computes the logarithm of the Binomial coefficient.
//
//  Discussion:
//
//    CNK_LOG = LOG ( C(N,K) ) = LOG ( N! / ( K! * (N-K)! ) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, K, are the values of N and K.
//
//    Output, double CNK_LOG, the logarithm of C(N,K).
//
{
  double cnk_log;

  cnk_log = factorial_log ( n ) - factorial_log ( k ) - factorial_log ( n - k );

  return cnk_log;
}
//****************************************************************************80

double binomial_mean ( int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_MEAN returns the mean of the Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the number of trials.
//    1 <= A.
//
//    Input, double B, the probability of success on one trial.
//    0.0 <= B <= 1.0.
//
//    Output, double BINOMIAL_MEAN, the expected value of the number of
//    successes in A trials.
//
{
  double mean;

  mean = ( double ) ( a ) * b;

  return mean;
}
//****************************************************************************80

double binomial_pdf ( int x, int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_PDF evaluates the Binomial PDF.
//
//  Discussion:
//
//    PDF(A,B;X) is the probability of exactly X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    PDF(A,B;X) = C(N,X) * B**X * ( 1.0 - B )**( A - X )
//
//    Binomial_PDF(1,B;X) = Bernoulli_PDF(B;X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the desired number of successes.
//    0 <= X <= A.
//
//    Input, int A, the number of trials.
//    1 <= A.
//
//    Input, double B, the probability of success on one trial.
//    0.0 <= B <= 1.0.
//
//    Output, double BINOMIAL_PDF, the value of the PDF.
//
{
  int cnk;
  double pdf;

  if ( a < 1 )
  {
    pdf = 0.0;
  }
  else if ( x < 0 || a < x )
  {
    pdf = 0.0;
  }
  else if ( b == 0.0 )
  {
    if ( x == 0 )
    {
      pdf = 1.0;
    }
    else
    {
      pdf = 0.0;
    }
  }
  else if ( b == 1.0 )
  {
    if ( x == a )
    {
      pdf = 1.0;
    }
    else
    {
      pdf = 0.0;
    }
  }
  else
  {
    cnk = binomial_coef ( a, x );

    pdf = ( double ) ( cnk ) * pow ( b, x ) * pow ( ( 1.0 - b ), ( a - x ) );
  }

  return pdf;
}
//****************************************************************************80

int binomial_sample ( int a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_SAMPLE samples the Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Kennedy, James Gentle,
//    Algorithm BU,
//    Statistical Computing,
//    Dekker, 1980.
//
//  Parameters:
//
//    Input, int A, the number of trials.
//    1 <= A.
//
//    Input, double B, the probability of success on one trial.
//    0.0 <= B <= 1.0.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int BINOMIAL_SAMPLE, a sample of the PDF.
//
{
  int i;
  double u;
  int x;

  x = 0;

  for ( i = 1; i <= a; i++ )
  {
    u = r8_uniform_01 ( seed );

    if ( u <= b )
    {
      x = x + 1;
    }

  }

  return x;
}
//****************************************************************************80

double binomial_variance ( int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_VARIANCE returns the variance of the Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the number of trials.
//    1 <= A.
//
//    Input, double B, the probability of success on one trial.
//    0.0 <= B <= 1.0.
//
//    Output, double BINOMIAL_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( double ) ( a ) * b * ( 1.0 - b );

  return variance;
}
//****************************************************************************80

double birthday_cdf ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BIRTHDAY_CDF returns the Birthday Concurrence CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of people whose birthdays have been
//    disclosed.
//
//    Output, double BIRTHDAY_CDF, the probability that at least
//    two of the N people have matching birthays..
//
{
  double cdf;
  int i;

  if ( n < 1 )
  {
    cdf = 0.0;
    return cdf;
  }
  else if ( 365 < n )
  {
    cdf = 1.0;
    return cdf;
  }
//
//  Compute the probability that N people have distinct birthdays.
//
  cdf = 1.0;
  for ( i = 1; i <= n; i++ )
  {
    cdf = cdf * ( double ) ( 365 + 1 - i ) / 365.0;
  }
//
//  Compute the probability that it is NOT the case that N people
//  have distinct birthdays.  This is the cumulative probability
//  that person 2 matches person 1, or person 3 matches 1 or 2,
//  etc.
//
  cdf = 1.0 - cdf;

  return cdf;
}
//****************************************************************************80

int birthday_cdf_inv ( double cdf )

//****************************************************************************80
//
//  Purpose:
//
//    BIRTHDAY_CDF_INV inverts the Birthday Concurrence CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the probability that at least
//    two of the N people have matching birthays..
//
//    Output, int BIRTHDAY_CDF_INV, the corresponding number of people whose
//    birthdays need to be disclosed.
//
{
  double cdf_not;
  int i;
  int n;

  if ( cdf <= 0.0 )
  {
    n = 0;
    return n;
  }
  else if ( 1.0 <= cdf )
  {
    n = 365;
    return n;
  }
//
//  Compute the probability that N people have distinct birthdays.
//
  cdf_not = 1.0;

  for ( i = 1; i <= 365; i++ )
  {
    cdf_not = cdf_not * ( double ) ( 365 + 1 - i ) / 365.0;
    if ( cdf <= 1.0 - cdf_not )
    {
      n = i;
      return n;
    }
  }
  n = 365;
  return n;
}
//****************************************************************************80

double birthday_pdf ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BIRTHDAY_PDF returns the Birthday Concurrence PDF.
//
//  Discussion:
//
//    The probability is the probability that the N-th person is the
//    first one to match a birthday with someone earlier.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of people whose birthdays have been
//    disclosed.
//
//    Output, double BIRTHDAY_PDF, the probability that the N-th person
//    is the first to match a birthday with someone earlier.
//
{
  int i;
  double pdf;

  if ( n < 1 || 365 < n )
  {
    pdf = 0.0;
    return pdf;
  }
  pdf = 1.0;
//
//  Compute the probability that N-1 people have distinct birthdays.
//
  for ( i = 1; i <= n-1; i++ )
  {
    pdf = pdf * ( double ) ( 365 + 1 - i ) / 365.0;
  }
//
//  Compute the probability that person N has one of those N-1 birthdays.
//
  pdf = pdf * ( double ) ( n - 1 ) / 365.0;

  return pdf;
}
//****************************************************************************80

double bradford_cdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    BRADFORD_CDF evaluates the Bradford CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    A < B,
//    0.0 < C.
//
//    Output, double BRADFORD_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else if ( x <= b )
  {
    cdf = log ( 1.0 + c * ( x - a ) / ( b - a ) ) / log ( c + 1.0 );
  }
  else if ( b < x )
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double bradford_cdf_inv ( double cdf, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    BRADFORD_CDF_INV inverts the Bradford CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, C, the parameters of the PDF.
//    A < B,
//    0.0 < C.
//
//    Output, double BRADFORD_CDF_INV, the corresponding argument of the CDF.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "BRADFORD_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf <= 0.0 )
  {
    x = a;
  }
  else if ( cdf < 1.0 )
  {
    x = a + ( b - a ) * ( pow ( ( c + 1.0 ), cdf ) - 1.0 ) / c;
  }
  else if ( 1.0 <= cdf )
  {
    x = b;
  }

  return x;
}
//****************************************************************************80

bool bradford_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    BRADFORD_CHECK checks the parameters of the Bradford PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    A < B,
//    0.0 < C.
//
//    Output, bool BRADFORD_CHECK, is true if the parameters are legal.
//
{
  if ( b <= a )
  {
    cout << " \n";
    cout << "BRADFORD_CHECK - Fatal error!\n";
    cout << "  B <= A.\n";
    return false;
  }

  if ( c <= 0.0 )
  {
    cout << " \n";
    cout << "BRADFORD_CHECK - Fatal error!\n";
    cout << "  C <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double bradford_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    BRADFORD_MEAN returns the mean of the Bradford PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    A < B,
//    0.0 < C.
//
//    Output, double BRADFORD_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = ( c * ( b - a ) + log ( c + 1.0 ) * ( a * ( c + 1.0 ) - b ) )
    / ( c * log ( c + 1.0 ) );

  return mean;
}
//****************************************************************************80

double bradford_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    BRADFORD_PDF evaluates the Bradford PDF.
//
//  Discussion:
//
//    PDF(A,B,C;X) = C / ( ( C * ( X - A ) + B - A ) * log ( C + 1 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, C, the parameters of the PDF.
//    A < B,
//    0.0 < C.
//
//    Output, double BRADFORD_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else if ( x <= b )
  {
    pdf = c / ( ( c * ( x - a ) + b - a ) * log ( c + 1.0 ) );
  }
  else if ( b < x )
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

double bradford_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    BRADFORD_SAMPLE samples the Bradford PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    A < B,
//    0.0 < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double BRADFORD_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = a + ( b - a ) * ( pow ( ( c + 1.0 ), cdf ) - 1.0 ) / c;

  return x;
}
//****************************************************************************80

double bradford_variance ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    BRADFORD_VARIANCE returns the variance of the Bradford PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    A < B,
//    0.0 < C.
//
//    Output, double BRADFORD_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( b - a ) * ( b - a ) *
    ( c * ( log ( c + 1.0 ) - 2.0 ) + 2.0 * log ( c + 1.0 ) )
    / ( 2.0 * c * pow ( ( log ( c + 1.0 ) ), 2 ) );

  return variance;
}
//****************************************************************************80

double buffon_laplace_pdf ( double a, double b, double l )

//****************************************************************************80
//
//  Purpose:
//
//    BUFFON_LAPLACE_PDF evaluates the Buffon-Laplace PDF.
//
//  Discussion:
//
//    In the Buffon-Laplace needle experiment, we suppose that the plane has
//    been tiled into a grid of rectangles of width A and height B, and that a
//    needle of length L is dropped "at random" onto this grid.
//
//    We may assume that one end, the "eye" of the needle falls at the point
//    (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].
//
//    ANGLE, the angle that the needle makes is taken to be uniformly random.
//    The point of the needle, (X2,Y2), therefore lies at
//
//      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
//
//    The needle will have crossed at least one grid line if any of the
//    following are true:
//
//      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.
//
//    If L is larger than sqrt ( A*A + B*B ), then the needle will
//    cross every time, and the computation is uninteresting.  However, if
//    L is smaller than this limit, then the probability of a crossing on
//    a single trial is
//
//      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( PI * A * B )
//
//    and therefore, a record of the number of hits for a given number of
//    trials can be used as a very roundabout way of estimating PI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Sudarshan Raghunathan,
//    Making a Supercomputer Do What You Want: High Level Tools for
//    Parallel Programming,
//    Computing in Science and Engineering,
//    Volume 8, Number 5, September/October 2006, pages 70-80.
//
//  Parameters:
//
//    Input, double A, B, the horizontal and vertical dimensions
//    of each cell of the grid.  0 <= A, 0 <= B.
//
//    Input, double L, the length of the needle.
//    0 <= L <= min ( A, B ).
//
//    Output, double BUFFON_LAPLACE_PDF, the Buffon-Laplace PDF.
//
{
  double pdf;
  double pi = 3.141592653589793238462643;

  if ( a < 0.0 )
  {
    cout << "\n";
    cout << "BUFFON_LAPLACE_PDF - Fatal error!\n";
    cout << "  Input A < 0.\n";
    exit ( 1 );
  }
  else if ( a == 0.0 )
  {
    pdf = 1.0;
    return pdf;
  }

  if ( b < 0.0 )
  {
    cout << "\n";
    cout << "BUFFON_LAPLACE_PDF - Fatal error!\n";
    cout << "  Input B < 0.\n";
    exit ( 1 );
  }
  else if ( b == 0.0 )
  {
    pdf = 1.0;
    return pdf;
  }

  if ( l < 0.0 )
  {
    cout << "\n";
    cout << "BUFFON_LAPLACE_PDF - Fatal error!\n";
    cout << "  Input L < 0.\n";
    exit ( 1 );
  }
  else if ( l == 0.0 )
  {
    pdf = 0.0;
    return pdf;
  }
  else if ( r8_min ( a, b ) < l )
  {
    cout << "\n";
    cout << "BUFFON_LAPLACE_PDF - Fatal error!\n";
    cout << "  min ( A, B ) < L.\n";
    exit ( 1 );
  }

  pdf = ( 2.0 * l * ( a + b ) - l * l ) / ( pi * a * b );

  return pdf;
}
//****************************************************************************80

int buffon_laplace_simulate ( double a, double b, double l, int trial_num )

//****************************************************************************80
//
//  Purpose:
//
//    BUFFON_LAPLACE_SIMULATE simulates a Buffon-Laplace needle experiment.
//
//  Discussion:
//
//    In the Buffon-Laplace needle experiment, we suppose that the plane has
//    been tiled into a grid of rectangles of width A and height B, and that a
//    needle of length L is dropped "at random" onto this grid.
//
//    We may assume that one end, the "eye" of the needle falls at the point
//    (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].
//
//    ANGLE, the angle that the needle makes is taken to be uniformly random.
//    The point of the needle, (X2,Y2), therefore lies at
//
//      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
//
//    The needle will have crossed at least one grid line if any of the
//    following are true:
//
//      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.
//
//    This routine simulates the tossing of the needle, and returns the number
//    of times that the needle crossed at least one grid line.
//
//    If L is larger than sqrt ( A*A + B*B ), then the needle will
//    cross every time, and the computation is uninteresting.  However, if
//    L is smaller than this limit, then the probability of a crossing on
//    a single trial is
//
//      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( PI * A * B )
//
//    and therefore, a record of the number of hits for a given number of
//    trials can be used as a very roundabout way of estimating PI.
//    (Particularly roundabout, since we actually will use a good value of
//    PI in order to pick the random angles//)
//
//    Since this routine invokes the C++ random number generator,
//    the user should initialize the random number generator, particularly
//    if it is desired to control whether the sequence is to be varied
//    or repeated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Sudarshan Raghunathan,
//    Making a Supercomputer Do What You Want: High Level Tools for
//    Parallel Programming,
//    Computing in Science and Engineering,
//    Volume 8, Number 5, September/October 2006, pages 70-80.
//
//  Parameters:
//
//    Input, double A, B, the horizontal and vertical dimensions
//    of each cell of the grid.  0 <= A, 0 <= B.
//
//    Input, double L, the length of the needle.
//    0 <= L <= min ( A, B ).
//
//    Input, int TRIAL_NUM, the number of times the needle is
//    to be dropped onto the grid.
//
//    Output, int BUFFON_LAPLACE_SIMULATE, the number of times the needle
//    crossed at least one line of the grid of cells.
//
{
  double angle;
  int hits;
  double pi = 3.141592653589793238462643;
  int trial;
  double x1;
  double x2;
  double y1;
  double y2;

  hits = 0;

  for ( trial = 1; trial <= trial_num; trial++ )
  {
//
//  Randomly choose the location of the eye of the needle in [0,0]x[A,B],
//  and the angle the needle makes.
//
    x1 = a * ( double ) rand ( ) / ( double ) RAND_MAX;
    y1 = b * ( double ) rand ( ) / ( double ) RAND_MAX;
    angle = 2.0 * pi * ( double ) rand ( ) / ( double ) RAND_MAX;
//
//  Compute the location of the point of the needle.
//
    x2 = x1 + l * cos ( angle );
    y2 = y1 + l * sin ( angle );
//
//  Count the end locations that lie outside the cell.
//
    if ( x2 <= 0.0 || a <= x2 || y2 <= 0.0 || b <= y2 )
    {
      hits = hits + 1;
    }
  }
  return hits;
}
//****************************************************************************80

double buffon_pdf ( double a, double l )

//****************************************************************************80
//
//  Purpose:
//
//    BUFFON_PDF evaluates the Buffon PDF.
//
//  Discussion:
//
//    In the Buffon needle experiment, we suppose that the plane has been
//    ruled by vertical lines with a spacing of A units, and that a
//    needle of length L is dropped "at random" onto this grid.
//
//    Because of the various symmetries, we may assume that this eye of
//    this needle lands in the first infinite strip, and we may further
//    assume that its Y coordinate is 0.  Thus, we have
//    the eye as (X1,Y1) with 0 <= X1 <= A and Y1 = 0.
//
//    ANGLE, the angle that the needle makes is taken to be uniformly random.
//    The point of the needle, (X2,Y2), therefore lies at
//
//      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
//
//    The needle will have crossed at least one grid line if any of the
//    following are true:
//
//      X2 <= 0, A <= X2.
//
//    The probability of a crossing on a single trial is
//
//      P(A,L) = ( 2 * L ) / ( PI * A )
//
//    and therefore, a record of the number of hits for a given number of
//    trials can be used as a very roundabout way of estimating PI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the horizontal spacing between the
//    vertical grid lines.  0 <= A.
//
//    Input, double L, the length of the needle.
//
//    Output, double BUFFON_PDF, the Buffon PDF.
//
{
  double pdf;
  double pi = 3.141592653589793238462643;

  if ( a < 0.0 )
  {
    cout << "\n";
    cout << "BUFFON_PDF - Fatal error!\n";
    cout << "  Input A < 0.\n";
    exit ( 1 );
  }
  else if ( a == 0.0 )
  {
    pdf = 1.0;
    return pdf;
  }

  if ( l < 0.0 )
  {
    cout << "\n";
    cout << "BUFFON_PDF - Fatal error!\n";
    cout << "  Input L < 0.\n";
    exit ( 1 );
  }
  else if ( l == 0.0 )
  {
    pdf = 0.0;
    return pdf;
  }

  pdf = ( 2.0 * l ) / ( pi * a );

  return pdf;
}
//****************************************************************************80

int buffon_simulate ( double a, double l, int trial_num )

//****************************************************************************80
//
//  Purpose:
//
//    BUFFON_SIMULATE simulates a Buffon needle experiment.
//
//  Discussion:
//
//    In the Buffon needle experiment, we suppose that the plane has been
//    ruled by vertical lines with a spacing of A units, and that a
//    needle of length L is dropped "at random" onto this grid.
//
//    Because of the various symmetries, we may assume that this eye of
//    this needle lands in the first infinite strip, and we may further
//    assume that its Y coordinate is 0.  Thus, we have
//    the eye as (X1,Y1) with 0 <= X1 <= A and Y1 = 0.
//
//    ANGLE, the angle that the needle makes is taken to be uniformly random.
//    The point of the needle, (X2,Y2), therefore lies at
//
//      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
//
//    The needle will have crossed at least one grid line if any of the
//    following are true:
//
//      X2 <= 0, A <= X2.
//
//    The probability of a crossing on a single trial is
//
//      P(A,L) = ( 2 * L ) / ( PI * A )
//
//    and therefore, a record of the number of hits for a given number of
//    trials can be used as a very roundabout way of estimating PI.
//
//    Since this routine invokes the C++ random number generator,
//    the user should initialize the random number generator, particularly
//    if it is desired to control whether the sequence is to be varied
//    or repeated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the horizontal spacing between the
//    vertical grid lines.  0 <= A.
//
//    Input, double L, the length of the needle.
//
//    Input, int TRIAL_NUM, the number of times the needle is
//    to be dropped onto the grid.
//
//    Output, int BUFFON_SIMULATE, the number of times the needle
//    crossed at least one line of the grid of cells.
//
{
  double angle;
  int hits;
  double pi = 3.141592653589793238462643;
  int trial;
  double x1;
  double x2;

  hits = 0;

  for ( trial = 1; trial <= trial_num; trial++ )
  {
//
//  Randomly choose the location (X1,Y1) of the eye of the needle
//  in [0,0]x[A,0], and the angle the needle makes.
//
    x1 = a * ( double ) rand ( ) / ( double ) RAND_MAX;
    angle = 2.0 * pi * ( double ) rand ( ) / ( double ) RAND_MAX;
//
//  Compute the location of the point of the needle.
//  We only need to know the value of X2, not Y2!
//
    x2 = x1 + l * cos ( angle );
//
//  Count the end locations that lie outside the cell.
//
    if ( x2 <= 0.0 || a <= x2 )
    {
       hits = hits + 1;
    }
  }
  return hits;
}
//****************************************************************************80

double burr_cdf ( double x, double a, double b, double c, double d )

//****************************************************************************80
//
//  Purpose:
//
//    BURR_CDF evaluates the Burr CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, C, D, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double BURR_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 1.0 / pow ( ( 1.0 + pow ( ( b / ( x - a ) ), c )  ), d );
  }

  return cdf;
}
//****************************************************************************80

double burr_cdf_inv ( double cdf, double a, double b, double c, double d )

//****************************************************************************80
//
//  Purpose:
//
//    BURR_CDF_INV inverts the Burr CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, C, D, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double BURR_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "BURR_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a + b / pow ( ( pow ( ( 1.0 / cdf ), (1.0 / d ) ) - 1.0 ), ( 1.0 / c ) );

  return x;
}
//****************************************************************************80

bool burr_check ( double a, double b, double c, double d )

//****************************************************************************80
//
//  Purpose:
//
//    BURR_CHECK checks the parameters of the Burr CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, D, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, bool BURR_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "BURR_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  if ( c <= 0 )
  {
    cout << " \n";
    cout << "BURR_CHECK - Fatal error!\n";
    cout << "  C <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double burr_mean ( double a, double b, double c, double d )

//****************************************************************************80
//
//  Purpose:
//
//   BURR_MEAN returns the mean of the Burr PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, D, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double BURR_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a + b * r8_gamma ( 1.0 - 1.0 / c ) * r8_gamma ( d + 1.0 / c )
  / r8_gamma ( d );

  return mean;
}
//****************************************************************************80

double burr_pdf ( double x, double a, double b, double c, double d )

//****************************************************************************80
//
//  Purpose:
//
//    BURR_PDF evaluates the Burr PDF.
//
//  Discussion:
//
//    PDF(A,B,C,D;X) = ( C * D / B ) * ( ( X - A ) / B )**( - C - 1 )
//      * ( 1 + ( ( X - A ) / B )**( - C ) )**( - D - 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    M E Johnson,
//    Multivariate Statistical Simulation,
//    Wiley, New York, 1987.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, C, D, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double BURR_PDF, the value of the PDF.
//
{
  double pdf;
  double y;

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else
  {
    y = ( x - a ) / b;

    pdf = ( c * d / b ) * pow ( y, ( - c - 1.0 ) )
      * pow ( ( 1.0 + pow ( y, ( - c ) ) ), ( - d - 1.0 ) );

  }

  return pdf;
}
//****************************************************************************80

double burr_sample ( double a, double b, double c, double d, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    BURR_SAMPLE samples the Burr PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, D, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double BURR_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = burr_cdf_inv ( cdf, a, b, c, d );

  return x;
}
//****************************************************************************80

double burr_variance ( double a, double b, double c, double d )

//****************************************************************************80
//
//  Purpose:
//
//    BURR_VARIANCE returns the variance of the Burr PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, D, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double BURR_VARIANCE, the variance of the PDF.
//
{
  double k;
  double variance;

  if ( c <= 2.0 )
  {
    cout << " \n";
    cout << "BURR_VARIANCE - Warning!\n";
    cout << "  Variance undefined for C <= 2.\n";
    variance = r8_huge ( );
  }
  else
  {
    k = r8_gamma ( d ) * r8_gamma ( 1.0 - 2.0 / c ) * r8_gamma ( d + 2.0 / c )
      - pow ( ( r8_gamma ( 1.0 - 1.0 / c ) * r8_gamma ( d + 1.0 / c ) ), 2 );

    variance = k * b * b / pow ( r8_gamma ( d ), 2 );
  }

  return variance;
}
//****************************************************************************80

double cardioid_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CARDIOID_CDF evaluates the Cardioid CDF.
//
//  Discussion:
//
//    The angle X is assumed to lie between A - PI and A + PI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= B <= 0.5.
//
//    Output, double CDF, the value of the PDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;

  if ( x <= a - pi )
  {
    cdf = 0.0;
  }
  else if ( x < a + pi )
  {
    cdf = ( pi + x - a + 2.0 * b * sin ( x - a ) ) / ( 2.0 * pi );
  }
  else
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double cardioid_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CARDIOID_CDF_INV inverts the Cardioid CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0 <= CDF <= 1.
//
//    Input, double A, B, the parameters.
//    0.0 <= B <= 0.5.
//
//    Output, double X, the argument with the given CDF.
//    A - PI <= X <= A + PI.
//
{
  double fp;
  double fx;
  int it;
  const double pi = 3.14159265358979323;
  double tol = 0.000001;
  double x;

  if ( cdf <= 0.0 )
  {
    x = a - pi;
  }
  else if ( cdf < 1.0 )
  {
    x = a;

    it = 0;

    for ( ; ; )
    {
      fx = cdf - ( pi + x - a + 2.0 * b * sin ( x - a ) ) / ( 2.0 * pi );

      if ( r8_abs ( fx ) < tol )
      {
        break;
      }

      if ( 10 < it )
      {
        cout << "\n";
        cout << "CARDIOID_CDF_INV - Fatal error!\n";
        cout << "  Too many iterations.\n";
        exit ( 1 );
      }

      fp = - ( 1.0 + 2.0 * b * cos ( x - a ) ) / ( 2.0 * pi );

      x = x - fx / fp;
      x = r8_max ( x, a - pi );
      x = r8_min ( x, a + pi );

      it = it + 1;
    }
  }
  else
  {
    x = a + pi;
  }

  return x;
}
//****************************************************************************80

bool cardioid_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CARDIOID_CHECK checks the parameters of the Cardioid CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    -0.5 <= B <= 0.5.
//
//    Output, bool CARDIOID_CHECK, is true if the parameters are legal.
//
{
  bool value;

  if ( b < -0.5 || 0.5 < b )
  {
    cout << "\n";
    cout << "CARDIOID_CHECK - Fatal error!\n";
    cout << "  B < -0.5 or 0.5 < B.\n";
    value = false;
    return value;
  }

  value = true;

  return value;
}
//****************************************************************************80

double cardioid_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CARDIOID_MEAN returns the mean of the Cardioid PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= B <= 0.5.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double cardioid_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CARDIOID_PDF evaluates the Cardioid PDF.
//
//  Discussion:
//
//    The cardioid PDF can be thought of as being applied to points on
//    a circle.  Compare this distribution with the "Cosine PDF".
//
//    PDF(A,B;X) = ( 1 / ( 2 * PI ) ) * ( 1 + 2 * B * COS ( X - A ) )
//    for 0 <= B <= 1/2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    N I Fisher,
//    Statistical Analysis of Circular Data,
//    Cambridge, 1993.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= B <= 0.5.
//
//    Output, double CARDIOID_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  pdf = ( 1.0 + 2.0 * b * cos ( x - a ) ) / ( 2.0 * pi );

  return pdf;
}
//****************************************************************************80

double cardioid_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CARDIOID_SAMPLE samples the Cardioid PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= B <= 0.5.
//
//    Output, double X, a sample of the PDF.
//    A - PI <= X <= A + PI.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = cardioid_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double cardioid_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CARDIOID_VARIANCE returns the variance of the Cardioid PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= B <= 0.5.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = 0.0;

  return variance;
}
//****************************************************************************80

double cauchy_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_CDF evaluates the Cauchy CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;
  double y;

  y = ( x - a ) / b;

  cdf = 0.5 + atan ( y ) / pi;

  return cdf;
}
//****************************************************************************80

double cauchy_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_CDF_INV inverts the Cauchy CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double X, the corresponding argument.
//
{
  const double pi = 3.14159265358979323;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "CAUCHY_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a + b * tan ( pi * ( cdf - 0.5 ) );

  return x;
}
//****************************************************************************80

void cauchy_cdf_values ( int *n_data, double *mu, double *sigma, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_CDF_VALUES returns some values of the Cauchy CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = CauchyDistribution [ mu, sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
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
//    Output, double *MU, the mean of the distribution.
//
//    Output, double *SIGMA, the variance of the distribution.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double fx_vec[N_MAX] = {
     0.5000000000000000,
     0.8524163823495667,
     0.9220208696226307,
     0.9474315432887466,
     0.6475836176504333,
     0.6024163823495667,
     0.5779791303773693,
     0.5628329581890012,
     0.6475836176504333,
     0.5000000000000000,
     0.3524163823495667,
     0.2500000000000000 };

  double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  double sigma_vec[N_MAX] = {
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *mu = 0.0;
    *sigma = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *mu = mu_vec[*n_data-1];
    *sigma = sigma_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool cauchy_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_CHECK checks the parameters of the Cauchy CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, bool CAUCHY_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "CAUCHY_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double cauchy_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_MEAN returns the mean of the Cauchy PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double cauchy_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_PDF evaluates the Cauchy PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = 1 / ( PI * B * ( 1 + ( ( X - A ) / B )^2 ) )
//
//    The Cauchy PDF is also known as the Breit-Wigner PDF.  It
//    has some unusual properties.  In particular, the integrals for the
//    expected value and higher order moments are "singular", in the
//    sense that the limiting values do not exist.  A result can be
//    obtained if the upper and lower limits of integration are set
//    equal to +T and -T, and the limit as T=>INFINITY is taken, but
//    this is a very weak and unreliable sort of limit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;
  double y;

  y = ( x - a ) / b;

  pdf = 1.0 / ( pi * b * ( 1.0 + y * y ) );

  return pdf;
}
//****************************************************************************80

double cauchy_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_SAMPLE samples the Cauchy PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = cauchy_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double cauchy_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_VARIANCE returns the variance of the Cauchy PDF.
//
//  Discussion:
//
//    The variance of the Cauchy PDF is not well defined.  This routine
//    is made available for completeness only, and simply returns
//    a "very large" number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double VARIANCE, the mean of the PDF.
//
{
  double variance;

  variance = r8_huge ( );

  return variance;
}
//****************************************************************************80

double chi_cdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_CDF evaluates the Chi CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;
  double p2;
  double x2;
  double y;

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else
  {
    y = ( x - a ) / b;
    x2 = 0.5 * y * y;
    p2 = 0.5 * c;

    cdf = gamma_inc ( p2, x2 );
  }

  return cdf;
}
//****************************************************************************80

double chi_cdf_inv ( double cdf, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_CDF_INV inverts the Chi CDF.
//
//  Discussion:
//
//    A simple bisection method is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double X, the corresponding argument of the CDF.
//
{
  double cdf1;
  double cdf2;
  double cdf3;
  int it;
  int it_max = 100;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "CHI_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = a;
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = r8_huge ( );
    return x;
  }

  x1 = a;
  cdf1 = 0.0;

  x2 = a + 1.0;

  for ( ; ; )
  {
    cdf2 = chi_cdf ( x2, a, b, c );

    if ( cdf < cdf2 )
    {
      break;
    }

    x2 = a + 2.0 * ( x2 - a );
  }
//
//  Now use bisection.
//
  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = chi_cdf ( x3, a, b, c );

    if ( r8_abs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      break;
    }

    if ( it_max < it )
    {
      cout << " \n";
      cout << "CHI_CDF_INV - Fatal error!\n";
      cout << "  Iteration limit exceeded.\n";
      return x;
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
      cdf2 = cdf3;
    }

  }

  return x;
}
//****************************************************************************80

bool chi_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_CHECK checks the parameters of the Chi CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, bool CHI_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "CHI_CHECK - Fatal error!\n";
    cout << "  B <= 0.0.\n";
    return false;
  }

  if ( c <= 0.0 )
  {
    cout << " \n";
    cout << "CHI_CHECK - Fatal error!\n";
    cout << "  C <= 0.0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double chi_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_MEAN returns the mean of the Chi PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double MEAN, the mean value.
//
{
  double mean;

  mean = a + sqrt ( 2.0 ) * b * r8_gamma ( 0.5 * ( c + 1.0 ) )
  / r8_gamma ( 0.5 * c );

  return mean;
}
//****************************************************************************80

double chi_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_PDF evaluates the Chi PDF.
//
//  Discussion:
//
//    PDF(A,B,C;X) = EXP ( - 0.5 * ( ( X - A ) / B )^2 )
//      * ( ( X - A ) / B )**( C - 1 ) /
//      ( 2**( 0.5 * C - 1 ) * B * GAMMA ( 0.5 * C ) )
//
//    CHI(A,B,1) is the Half Normal PDF;
//    CHI(0,B,2) is the Rayleigh PDF;
//    CHI(0,B,3) is the Maxwell PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, C, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  double y;

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else
  {
    y = ( x - a ) / b;

    pdf = exp ( - 0.5 * y * y ) * pow ( y, ( c - 1.0 ) ) /
      ( pow ( 2.0, ( 0.5 * c - 1.0 ) ) * b * r8_gamma ( 0.5 * c ) );
  }

  return pdf;
}
//****************************************************************************80

double chi_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SAMPLE samples the Chi PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double x;

  x = chi_square_sample ( c, seed );

  x = a + b * sqrt ( x );

  return x;
}
//****************************************************************************80

double chi_variance ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_VARIANCE returns the variance of the Chi PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0 < B,
//    0 < C.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = b * b * ( c - 2.0 *
    pow ( ( r8_gamma ( 0.5 * ( c + 1.0 ) ) / r8_gamma ( 0.5 * c ) ), 2 ) );

  return variance;
}
//****************************************************************************80

double chi_square_cdf ( double x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_CDF evaluates the Chi squared CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value of the random deviate.
//
//    Input, double A, the parameter of the distribution, usually
//    the number of degrees of freedom.
//
//    Output, double CDF, the value of the CDF.
//
{
  double a2;
  double b2;
  double c2;
  double cdf;
  double x2;

  x2 = 0.5 * x;
  a2 = 0.0;
  b2 = 1.0;
  c2 = 0.5 * a;

  cdf = gamma_cdf ( x2, a2, b2, c2 );

  return cdf;
}
//****************************************************************************80

double chi_square_cdf_inv ( double cdf, double a )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_CDF_INV inverts the Chi squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    Original FORTRAN77 version by Best and Roberts.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Best, Roberts,
//    The Percentage Points of the Chi-Squared Distribution,
//    Algorithm AS 91,
//    Applied Statistics,
//    Volume 24, Number ?, pages 385-390, 1975.
//
//  Parameters:
//
//    Input, double CDF, a value of the chi-squared cumulative
//    probability density function.
//    0.000002 <= CDF <= 0.999998.
//
//    Input, double A, the parameter of the chi-squared
//    probability density function.  0 < A.
//
//    Output, double X, the value of the chi-squared random deviate
//    with the property that the probability that a chi-squared random
//    deviate with parameter A is less than or equal to X is CDF.
//
{
  double a2;
  double aa = 0.6931471806;
  double b;
  double c;
  double c1 = 0.01;
  double c2 = 0.222222;
  double c3 = 0.32;
  double c4 = 0.4;
  double c5 = 1.24;
  double c6 = 2.2;
  double c7 = 4.67;
  double c8 = 6.66;
  double c9 = 6.73;
  double c10 = 13.32;
  double c11 = 60.0;
  double c12 = 70.0;
  double c13 = 84.0;
  double c14 = 105.0;
  double c15 = 120.0;
  double c16 = 127.0;
  double c17 = 140.0;
  double c18 = 175.0;
  double c19 = 210.0;
  double c20 = 252.0;
  double c21 = 264.0;
  double c22 = 294.0;
  double c23 = 346.0;
  double c24 = 420.0;
  double c25 = 462.0;
  double c26 = 606.0;
  double c27 = 672.0;
  double c28 = 707.0;
  double c29 = 735.0;
  double c30 = 889.0;
  double c31 = 932.0;
  double c32 = 966.0;
  double c33 = 1141.0;
  double c34 = 1182.0;
  double c35 = 1278.0;
  double c36 = 1740.0;
  double c37 = 2520.0;
  double c38 = 5040.0;
  double cdf_max = 0.999998;
  double cdf_min = 0.000002;
  double ch;
  double e = 0.0000005;
  double g;
  int i;
  int it_max = 20;
  double p1;
  double p2;
  double q;
  double s1;
  double s2;
  double s3;
  double s4;
  double s5;
  double s6;
  double t;
  double x;
  double x2;
  double xx;
//
  if ( cdf < cdf_min )
  {
    x = -1.0;
    cout << " \n";
    cout << "CHI_SQUARE_CDF_INV - Fatal error!\n";
    cout << "  CDF < CDF_MIN.\n";
    exit ( 1 );
  }

  if ( cdf_max < cdf )
  {
    x = -1.0;
    cout << " \n";
    cout << "CHI_SQUARE_CDF_INV - Fatal error!\n";
    cout << "  CDF_MAX < CDF.\n";
    exit ( 1 );
  }

  xx = 0.5 * a;
  c = xx - 1.0;
//
//  Compute Log ( Gamma ( A/2 ) ).
//
  g = gamma_log ( a / 2.0 );
//
//  Starting approximation for small chi-squared.
//
  if ( a < -c5 * log ( cdf ) )
  {
    ch = pow ( ( cdf * xx * exp ( g + xx * aa ) ), ( 1.0 / xx ) );

    if ( ch < e )
    {
      x = ch;
      return x;
    }
  }
//
//  Starting approximation for A less than or equal to 0.32.
//
  else if ( a <= c3 )
  {
    ch = c4;
    a2 = log ( 1.0 - cdf );

    for ( ; ; )
    {
      q = ch;
      p1 = 1.0 + ch * ( c7 + ch );
      p2 = ch * ( c9 + ch * ( c8 + ch ) );

      t = - 0.5 + ( c7 + 2.0 * ch ) / p1 - ( c9 + ch * ( c10 + 3.0 * ch ) ) / p2;

      ch = ch - ( 1.0 - exp ( a2 + g + 0.5 * ch + c * aa ) * p2 / p1 ) / t;

      if ( r8_abs ( q / ch - 1.0 ) <= c1 )
      {
        break;
      }
    }
  }
//
//  Call to algorithm AS 111.
//  Note that P has been tested above.
//  AS 241 could be used as an alternative.
//
  else
  {
    x2 = normal_01_cdf_inv ( cdf );
//
//  Starting approximation using Wilson and Hilferty estimate.
//
    p1 = c2 / a;
    ch = a * pow ( ( x2 * sqrt ( p1 ) + 1.0 - p1 ), 3 );
//
//  Starting approximation for P tending to 1.
//
    if ( c6 * a + 6.0 < ch )
    {
      ch = - 2.0 * ( log ( 1.0 - cdf ) - c * log ( 0.5 * ch ) + g );
    }
  }
//
//  Call to algorithm AS 239 and calculation of seven term Taylor series.
//
  for ( i = 1; i <= it_max; i++ )
  {
    q = ch;
    p1 = 0.5 * ch;
    p2 = cdf - gamma_inc ( xx, p1 );
    t = p2 * exp ( xx * aa + g + p1 - c * log ( ch ) );
    b = t / ch;
    a2 = 0.5 * t - b * c;

    s1 = ( c19 + a2
       * ( c17 + a2
       * ( c14 + a2
       * ( c13 + a2
       * ( c12 + a2
       *   c11 ) ) ) ) ) / c24;

    s2 = ( c24 + a2
       * ( c29 + a2
       * ( c32 + a2
       * ( c33 + a2
       *   c35 ) ) ) ) / c37;

    s3 = ( c19 + a2
       * ( c25 + a2
       * ( c28 + a2
       *   c31 ) ) ) / c37;

    s4 = ( c20 + a2
       * ( c27 + a2
       *   c34 ) + c
       * ( c22 + a2
       * ( c30 + a2
       *   c36 ) ) ) / c38;

    s5 = ( c13 + c21 * a2 + c * ( c18 + c26 * a2 ) ) / c37;

    s6 = ( c15 + c * ( c23 + c16 * c ) ) / c38;

    ch = ch + t * ( 1.0 + 0.5 * t * s1 - b * c
      * ( s1 - b
      * ( s2 - b
      * ( s3 - b
      * ( s4 - b
      * ( s5 - b
      *   s6 ) ) ) ) ) );

    if ( e < r8_abs ( q / ch - 1.0 ) )
    {
      x = ch;
      return x;
    }

  }

  x = ch;
  cout << " \n";
  cout << "CHI_SQUARE_CDF_INV - Warning!\n";
  cout << "  Convergence not reached.\n";

  return x;
}
//****************************************************************************80

void chi_square_cdf_values ( int *n_data, int *a, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = ChiSquareDistribution [ df ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 August 2004
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
//    Output, int *A, the parameter of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
//
# define N_MAX 21

  int a_vec[N_MAX] = {
     1,  2,  1,  2,
     1,  2,  3,  4,
     1,  2,  3,  4,
     5,  3,  3,  3,
     3,  3, 10, 10,
    10 };

  double fx_vec[N_MAX] = {
     0.7965567455405796E-01,
     0.4987520807317687E-02,
     0.1124629160182849,
     0.9950166250831946E-02,
     0.4729107431344619,
     0.1812692469220181,
     0.5975750516063926E-01,
     0.1752309630642177E-01,
     0.6826894921370859,
     0.3934693402873666,
     0.1987480430987992,
     0.9020401043104986E-01,
     0.3743422675270363E-01,
     0.4275932955291202,
     0.6083748237289110,
     0.7385358700508894,
     0.8282028557032669,
     0.8883897749052874,
     0.1721156299558408E-03,
     0.3659846827343712E-02,
     0.1857593622214067E-01 };

  double x_vec[N_MAX] = {
     0.01,
     0.01,
     0.02,
     0.02,
     0.40,
     0.40,
     0.40,
     0.40,
     1.00,
     1.00,
     1.00,
     1.00,
     1.00,
     2.00,
     3.00,
     4.00,
     5.00,
     6.00,
     1.00,
     2.00,
     3.00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool chi_square_check ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_CHECK checks the parameter of the central Chi squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the distribution.
//    1 <= A.
//
//    Output, bool CHI_SQUARE_CHECK, is true if the parameters are legal.
//
{
  if ( a < 1.0 )
  {
    cout << " \n";
    cout << "CHI_SQUARE_CHECK - Fatal error!\n";
    cout << "  A < 1.0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double chi_square_mean ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_MEAN returns the mean of the central Chi squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the distribution.
//    1 <= A.
//
//    Output, double MEAN, the mean value.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double chi_square_pdf ( double x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_PDF evaluates the central Chi squared PDF.
//
//  Discussion:
//
//    PDF(A;X) =
//      EXP ( - X / 2 ) * X**((A-2)/2) / ( 2**(A/2) * GAMMA ( A/2 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X
//
//    Input, double A, the parameter of the PDF.
//    1 <= A.
//
//    Output, double PDF, the value of the PDF.
//
{
  double b;
  double pdf;

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    b = a / 2.0;
    pdf = exp ( -0.5 * x ) * pow ( x, ( b - 1.0 ) ) / ( pow ( 2.0, b )
    * r8_gamma ( b ) );
  }

  return pdf;
}
//****************************************************************************80

double chi_square_sample ( double a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_SAMPLE samples the central Chi squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    1 <= A.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double a2;
  double b2;
  double c2;
  int i;
  int it_max = 100;
  int n;
  double x;
  double x2;

  n = ( int ) a;

  if ( ( double ) ( n ) == a && n <= it_max )
  {
    x = 0.0;
    for ( i = 1; i <= n; i++ )
    {
      x2 = normal_01_sample ( seed );
      x = x + x2 * x2;
    }
  }
  else
  {
    a2 = 0.0;
    b2 = 1.0;
    c2 = a / 2.0;

    x = gamma_sample ( a2, b2, c2, seed );

    x = 2.0 * x;
  }

  return x;
}
//****************************************************************************80

double chi_square_variance ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_VARIANCE returns the variance of the central Chi squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the distribution.
//    1 <= A.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = 2.0 * a;

  return variance;
}
//****************************************************************************80

void chi_square_noncentral_cdf_values ( int *n_data, int *df, double *lambda,
  double *x, double *cdf )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NoncentralChiSquareDistribution [ df, lambda ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2004
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
//    Output, int *DF, the number of degrees of freedom.
//
//    Output, double *LAMBDA, the noncentrality parameter.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *CDF, the noncentral chi CDF.
//
{
# define N_MAX 28

  double cdf_vec[N_MAX] = {
     0.8399444269398261,
     0.6959060300435139,
     0.5350879697078847,
     0.7647841496310313,
     0.6206436532195436,
     0.4691667375373180,
     0.3070884345937569,
     0.2203818092990903,
     0.1500251895581519,
     0.3071163194335791E-02,
     0.1763982670131894E-02,
     0.9816792594625022E-03,
     0.1651753140866208E-01,
     0.2023419573950451E-03,
     0.4984476352854074E-06,
     0.1513252400654827E-01,
     0.2090414910614367E-02,
     0.2465021206048452E-03,
     0.2636835050342939E-01,
     0.1857983220079215E-01,
     0.1305736595486640E-01,
     0.5838039534819351E-01,
     0.4249784402463712E-01,
     0.3082137716021596E-01,
     0.1057878223400849,
     0.7940842984598509E-01,
     0.5932010895599639E-01,
     0.2110395656918684 };

  int df_vec[N_MAX] = {
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
     60,  80, 100,
      1,   2,   3,
     10,  10,  10,
     10,  10,  10,
     10,  10,  10,
      8 };

  double lambda_vec[N_MAX] = {
     0.5,
     0.5,
     0.5,
     1.0,
     1.0,
     1.0,
     5.0,
     5.0,
     5.0,
    20.0,
    20.0,
    20.0,
    30.0,
    30.0,
    30.0,
     5.0,
     5.0,
     5.0,
     2.0,
     3.0,
     4.0,
     2.0,
     3.0,
     4.0,
     2.0,
     3.0,
     4.0,
     0.5 };

  double x_vec[N_MAX] = {
      3.000,
      3.000,
      3.000,
      3.000,
      3.000,
      3.000,
      3.000,
      3.000,
      3.000,
      3.000,
      3.000,
      3.000,
     60.000,
     60.000,
     60.000,
      0.050,
      0.050,
      0.050,
      4.000,
      4.000,
      4.000,
      5.000,
      5.000,
      5.000,
      6.000,
      6.000,
      6.000,
      5.000 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *lambda = 0.0;
    *df = 0;
    *cdf = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *df = df_vec[*n_data-1];
    *cdf = cdf_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double chi_square_noncentral_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_NONCENTRAL_CHECK checks the parameters of the noncentral Chi Squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the parameter of the PDF.
//    1.0 <= A.
//
//    Input, double B, the noncentrality parameter of the PDF.
//    0.0 <= B.
//
//    Output, bool CHI_SQUARE_NONCENTRAL_CHECK, is true if the parameters
//    are legal.
//
{
  if ( a < 1.0 )
  {
    cout << " \n";
    cout << "CHI_SQUARE_NONCENTRAL_CHECK - Fatal error!\n";
    cout << "  A < 1.\n";
    return false;
  }

  if ( b < 0.0 )
  {
    cout << " \n";
    cout << "CHI_SQUARE_NONCENTRAL_CHECK - Fatal error!\n";
    cout << "  B < 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double chi_square_noncentral_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_NONCENTRAL_MEAN returns the mean of the noncentral Chi squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the parameter of the PDF.
//    1.0 <= A.
//
//    Input, double B, the noncentrality parameter of the PDF.
//    0.0 <= B.
//
//    Output, double MEAN, the mean value.
//
{
  double mean;

  mean = a + b;

  return mean;
}
//****************************************************************************80

double chi_square_noncentral_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_NONCENTRAL_SAMPLE samples the noncentral Chi squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 November 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the parameter of the PDF.
//    1.0 <= A.
//
//    Input, double B, the noncentrality parameter of the PDF.
//    0.0 <= B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double a1;
  double a2;
  double b2;
  double x;
  double x1;
  double x2;

  a1 = a - 1.0;

  x1 = chi_square_sample ( a1, seed );

  a2 = sqrt ( b );
  b2 = 1.0;
  x2 = normal_sample ( a2, b2, seed );

  x = x1 + x2 * x2;

  return x;
}
//****************************************************************************80

double chi_square_noncentral_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_NONCENTRAL_VARIANCE returns the variance of the noncentral Chi squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    1 <= A.
//
//    Input, double B, the noncentrality parameter of the PDF.
//    0.0 <= B.
//
//    Output, double VARIANCE, the variance value.
//
{
  double variance;

  variance = 2.0 * ( a + 2.0 * b );

  return variance;
}
//****************************************************************************80

double *circle_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_SAMPLE samples points from a circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the circle.
//    The circle is centered at (A,B) and has radius C.
//    0.0 < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double CIRCLE_SAMPLE[2], a sampled point of the circle.
//
{
  double angle;
  const double pi = 3.14159265358979323;
  double radius_frac;
  double *x;

  x = new double[2];

  radius_frac = sqrt ( r8_uniform_01 ( seed ) );

  angle = 2.0 * pi * r8_uniform_01 ( seed );

  x[0] = a + c * radius_frac * cos ( angle );
  x[1] = b + c * radius_frac * sin ( angle );

  return x;
}
//****************************************************************************80

double *circular_normal_01_mean ( )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_01_MEAN returns the mean of the Circular Normal 01 PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CIRCULAR_01_MEAN[2], the mean of the PDF.
//
{
  double *mean;

  mean = new double[2];

  mean[0] = 0.0;
  mean[1] = 0.0;

  return mean;
}
//****************************************************************************80

double circular_normal_01_pdf ( double x[2] )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_01_PDF evaluates the Circular Normal 01 PDF.
//
//  Discussion:
//
//    PDF(X) = EXP ( - 0.5 * ( X(1)^2 + X(2)^2 ) ) / ( 2 * PI )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X[2], the argument of the PDF.
//
//    Output, double CIRCULAR_NORMAL_01_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  pdf = exp ( - 0.5 * ( x[0] * x[0] + x[1] * x[1] ) ) / ( 2.0 * pi );

  return pdf;
}
//****************************************************************************80

double *circular_normal_01_sample ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_01_SAMPLE samples the Circular Normal 01 PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double CIRCULAR_NORMAL_01_SAMPLE[2], a sample of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double v1;
  double v2;
  double *x;

  x = new double[2];

  v1 = r8_uniform_01 ( seed );
  v2 = r8_uniform_01 ( seed );

  x[0] = sqrt ( - 2.0 * log ( v1 ) ) * cos ( 2.0 * pi * v2 );
  x[1] = sqrt ( - 2.0 * log ( v1 ) ) * sin ( 2.0 * pi * v2 );

  return x;
}
//****************************************************************************80

double *circular_normal_01_variance ( )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_01_VARIANCE returns the variance of the Circular Normal 01 PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CIRCULAR_NORMAL_01_VARIANCE[2], the variance of the PDF.
//
{
  double *variance;

  variance = new double[2];

  variance[0] = 1.0;
  variance[1] = 1.0;

  return variance;
}
//****************************************************************************80

double *circular_normal_mean ( double a[2], double b )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_MEAN returns the mean of the Circular Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[2], a parameter of the PDF, the mean value.
//
//    Input, double B, a parameter of the PDF, the standard deviation.
//
//    Output, double CIRCULAR_MEAN[2], the mean of the PDF.
//
{
  double *mean;

  mean = new double[2];

  mean[0] = a[0];
  mean[1] = a[1];

  return mean;
}
//****************************************************************************80

double circular_normal_pdf ( double x[2], double a[2], double b )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_PDF evaluates the Circular Normal PDF.
//
//  Discussion:
//
//    PDF(X) = EXP ( - 0.5D+00 * ( ( (X(1)-A(1))^2 + (X(2)-A(2))^2 ) / B^2 ) 
//      / ( 2 * PI * B^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X[2], the argument of the PDF.
//
//    Input, double A[2], a parameter of the PDF, the mean value.
//
//    Input, double B, a parameter of the PDF, the standard deviation.
//
//    Output, double CIRCULAR_NORMAL_PDF, the value of the PDF.
//
{
  double d;
  double pdf;
  const double pi = 3.14159265358979323;

  d = ( pow ( x[0] - a[0], 2 ) 
      + pow ( x[1] - a[1], 2 ) ) / pow ( b, 2 );

  pdf = exp ( - 0.5 * d ) / ( 2.0 * b * b * pi );

  return pdf;
}
//****************************************************************************80

double *circular_normal_sample ( double a[2], double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_SAMPLE samples the Circular Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[2], a parameter of the PDF, the mean value.
//
//    Input, double B, a parameter of the PDF, the standard deviation.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double CIRCULAR_NORMAL_SAMPLE[2], a sample of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double r;
  double v1;
  double v2;
  double *x;

  x = new double[2];

  v1 = r8_uniform_01 ( seed );
  v2 = r8_uniform_01 ( seed );

  r = sqrt ( - 2.0 * log ( v1 ) );

  x[0] = a[0] + b * r * cos ( 2.0 * pi * v2 );
  x[1] = a[1] + b * r * sin ( 2.0 * pi * v2 );

  return x;
}
//****************************************************************************80

double *circular_normal_variance ( double a[2], double b )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCULAR_NORMAL_VARIANCE returns the variance of the Circular Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[2], a parameter of the PDF, the mean value.
//
//    Input, double B, a parameter of the PDF, the standard deviation.
//
//    Output, double CIRCULAR_NORMAL_VARIANCE[2], the variance of the PDF.
//
{
  double *variance;

  variance = new double[2];

  variance[0] = b;
  variance[1] = b;

  return variance;
}
//****************************************************************************80

int combinatorial ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    COMBINATORIAL computes the binomial coefficient C(N,K).
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
//    24 January 2007
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
//    April, 1963.
//
//  Parameters:
//
//    Input, int N, K, are the values of N and K.
//
//    Output, int COMBINATORIAL, the number of combinations of N
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

double cosine_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_CDF evaluates the Cosine CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameter of the PDF.
//    0.0 < B.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;
  double y;

  if ( x <= a - pi * b )
  {
    cdf = 0.0;
  }
  else if ( x <= a + pi * b )
  {
    y = ( x - a ) / b;

    cdf = 0.5 + ( y + sin ( y ) ) / ( 2.0 * pi );
  }
  else if ( a + pi * b < x )
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double cosine_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_CDF_INV inverts the Cosine CDF.
//
//  Discussion:
//
//    A simple bisection method is used on the interval
//    [ A - PI * B, A + PI * B ].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double X, the corresponding argument of the CDF.
//
{
  double cdf1;
  double cdf2;
  double cdf3;
  int it;
  int it_max = 100;
  const double pi = 3.14159265358979323;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "COSINE_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = a - pi * b;
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = a + pi * b;
    return x;
  }

  x1 = a - pi * b;
  cdf1 = 0.0;

  x2 = a + pi * b;
  cdf2 = 1.0;
//
//  Now use bisection.
//
  it = 0;

  for ( it = 1; it <= it_max; it++ )
  {
    x3 = 0.5 * ( x1 + x2 );
    cdf3 = cosine_cdf ( x3, a, b );

    if ( r8_abs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      return x;
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
      cdf2 = cdf3;
    }

  }

  cout << " \n";
  cout << "COSINE_CDF_INV - Fatal error!\n";
  cout << "  Iteration limit exceeded.\n";

  exit ( 1 );
  return x;
}
//****************************************************************************80

bool cosine_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_CHECK checks the parameters of the Cosine CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameter of the PDF.
//    0.0 < B.
//
//    Output, bool COSINE_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "COSINE_CHECK - Fatal error!\n";
    cout << "  B <= 0.0\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double cosine_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_MEAN returns the mean of the Cosine PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double cosine_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_PDF evaluates the Cosine PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = ( 1 / ( 2 * PI * B ) ) * COS ( ( X - A ) / B )
//    for A - PI * B <= X <= A + PI * B
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;
  double y;

  if ( x < a - pi * b )
  {
    pdf = 0.0;
  }
  else if ( x <= a + pi * b )
  {
    y = ( x - a ) / b;

    pdf = 1.0 / ( 2.0 * pi * b ) * cos ( y );
  }
  else if ( a + pi * b < x )
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

double cosine_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_SAMPLE samples the Cosine PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = cosine_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double cosine_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_VARIANCE returns the variance of the Cosine PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double variance;

  variance = ( pi * pi / 3.0 - 2.0 ) * b * b;

  return variance;
}
//****************************************************************************80

double coupon_complete_pdf ( int type_num, int box_num )

//****************************************************************************80
//
//  Purpose:
//
//    COUPON_COMPLETE_PDF evaluates the Complete Coupon Collection PDF.
//
//  Discussion:
//
//    PDF(TYPE_NUM;BOX_NUM) is the probability that, given an inexhaustible
//    supply of boxes, inside each of which there is one of TYPE_NUM distinct
//    coupons, which are uniformly distributed among the boxes, that it will
//    require opening exactly BOX_NUM boxes to achieve at least one of each
//    kind of coupon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Herbert Wilf,
//    Some New Aspects of the Coupon Collector's Problem,
//    SIAM Review,
//    Volume 48, Number 3, September 2006, pages 549-565.
//
//  Parameters:
//
//    Input, int BOX_NUM, the number of boxes that had to be opened
//    in order to just get at least one of each coupon.
//    0 <= BOX_NUM.  If BOX_NUM < TYPE_NUM, then PDF is surely 0.
//
//    Input, int TYPE_NUM, the number of distinct coupons.
//    1 <= TYPE_NUM.
//
//    Output, double COUPON_COMPLETE_PDF, the value of the PDF.
//
{
  double factor;
  int i;
  double pdf;
//
//  Nonsense cases.
//
  if ( box_num < 0 )
  {
    pdf = 0.0;
  }
  else if ( type_num < 1 )
  {
    pdf = 0.0;
  }
//
//  Degenerate but meaningful case.
//
  else if ( type_num == 1 )
  {
    if ( box_num == 1 )
    {
      pdf = 1.0;
    }
    else
    {
      pdf = 0.0;
    }
  }
//
//  Easy cases.
//
  else if ( box_num < type_num )
  {
    pdf = 0.0;
  }
//
//  General case.
//
  else
  {
    factor = 1.0;
    for ( i = 1; i <= type_num; i++ )
    {
      factor = factor * ( double ) ( i ) / ( double ) ( type_num );
    }
    for ( i = type_num+1; i <= box_num; i++ )
    {
      factor = factor / ( double ) ( type_num );
    }
    pdf = factor * ( double ) ( stirling2_value ( box_num-1, type_num-1 ) );
  }
  return pdf;
}
//****************************************************************************80

double coupon_mean ( int j, int type_num )

//****************************************************************************80
//
//  Purpose:
//
//    COUPON_MEAN returns the mean of the Coupon PDF.
//
//  Discussion:
//
//    In this version of the coupon collector's problem, we assume
//    that each box contains 1 coupon, that there are TYPE_NUM distinct types
//    of coupon, uniformly distributed among an inexhaustible supply
//    of boxes, and that the collector's goal is to get J distinct
//    types of coupons by opening one box after another.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int J, the number of distinct coupons to be collected.
//    J must be between 1 and N.
//
//    Input, int TYPE_NUM, the number of distinct coupons.
//
//    Output, double COUPON_MEAN, the mean number of boxes that
//    must be opened in order to just get J distinct kinds.
//
{
  int i;
  double mean;

  if ( type_num < j )
  {
    cout << " \n";
    cout << "COUPON_MEAN - Fatal error!\n";
    cout << "  Number of distinct coupons desired must be no more\n";
    cout << "  than the total number of boxes opened.\n";
    exit ( 1 );
  }

  mean = 0.0;

  for ( i = 1; i <= j; i++ )
  {
    mean = mean + 1.0 / ( double ) ( type_num - i + 1 );
  }
  mean = mean * ( double ) ( type_num );

  return mean;
}
//****************************************************************************80

void coupon_simulate ( int type_num, int *seed, int coupon[], int *box_num )

//****************************************************************************80
//
//  Purpose:
//
//    COUPON_SIMULATE simulates the coupon collector's problem.
//
//  Discussion:
//
//    The coupon collector needs to collect one of each of TYPE_NUM
//    coupons.  The collector may draw one coupon (or, in some settings,
//    open one box) on each trial, and takes as many trials as necessary
//    to complete the task.  On each trial, the probability of picking
//    any particular type of coupon is always 1 / TYPE_NUM.
//
//    Interesting questions include;
//
//    * what is the expected number of drawings necessary to complete
//      the collection?
//
//    * How does the expected number of drawings necessary to complete
//      the collection vary as TYPE_NUM increases?
//
//    * What is the distribution of the numbers of each type of coupon
//      in a typical collection when it is just completed?
//
//    As TYPE_NUM increases, the number of coupons necessary to be
//    collected in order to get a complete set in any simulation
//    strongly tends to the value TYPE_NUM * LOG ( TYPE_NUM ).
//
//    If TYPE_NUM is 1, the simulation ends with a single drawing.
//
//    If TYPE_NUM is 2, then we may call the coupon taken on the first drawing
//    a "Head", say, and the process then is similar to the question of the
//    length, plus one, of a run of Heads or Tails in coin flipping.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TYPE_NUM, the number of types of coupons.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int COUPON[TYPE_NUM], the number of coupons of each type
//    that were collected during the simulation.
//
//    Output, int *BOX_NUM, the total number of boxes opened.
//
{
  int i;
  int box_max = 2000;
  int straight;

  for ( i = 0; i < type_num; i++ )
  {
    coupon[i] = 0;
  }

  straight = 0;
  *box_num = 0;
//
//  Draw another coupon.
//
  while ( *box_num < box_max )
  {
    i = i4_uniform ( 1, type_num, seed );
//
//  Increment the number of I coupons.
//
    coupon[i-1] = coupon[i-1] + 1;
    *box_num = *box_num + 1;
//
//  If I is the next one we needed, increase STRAIGHT by 1.
//
    if ( i == straight + 1 )
    {
      for ( ; ; )
      {
        straight = straight + 1;
//
//  If STRAIGHT = TYPE_NUM, we have all of them.
//
        if ( type_num <= straight )
        {
          return;
        }
//
//  If the next coupon has not been collected, our straight is over.
//
        if ( coupon[straight] <= 0 )
        {
          break;
        }
      }
    }
  }

  cout << " \n";
  cout << "COUPON_SIMULATE - Fatal error!\n";
  cout << "  Maximum number of coupons drawn without success.\n";

  exit ( 1 );
}
//****************************************************************************80

double coupon_variance ( int j, int type_num )

//****************************************************************************80
//
//  Purpose:
//
//    COUPON_VARIANCE returns the variance of the Coupon PDF.
//
//  Discussion:
//
//    In this version of the coupon collector's problem, we assume
//    that each box contains 1 coupon, that there are TYPE_NUM distinct types
//    of coupon, uniformly distributed among an inexhaustible supply
//    of boxes, and that the collector's goal is to get J distinct
//    types of coupons by opening one box after another.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int J, the number of distinct coupons to be collected.
//    J must be between 1 and N.
//
//    Input, int TYPE_NUM, the number of distinct coupons.
//
//    Output, double COUPON_VARIANCE, the variance of the number of
//    boxes that must be opened in order to just get J distinct kinds.
//
{
  int i;
  double variance;

  if ( type_num < j )
  {
    cout << " \n";
    cout << "COUPON_VARIANCE - Fatal error!\n";
    cout << "  Number of distinct coupons desired must be no more\n";
    cout << "  than the total number of distinct coupons.\n";
    exit ( 1 );
  }

  variance = 0.0;
  for ( i = 1; i <= j; i++ )
  {
    variance = variance + ( double ) ( i - 1 ) /
      pow ( double ( type_num - i + 1 ), 2 );
  }
  variance = variance * ( double ) ( type_num );

  return variance;
}
//****************************************************************************80

double deranged_cdf ( int x, int a )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_CDF evaluates the Deranged CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the maximum number of items in their correct places.
//    0 <= X <= A.
//
//    Input, int A, the number of items.
//    1 <= A.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;
  int cnk;
  int dnmk;
  int sum2;
  int x2;

  if ( x < 0 || a < x )
  {
    cdf = 0.0;
  }
  else
  {
    sum2 = 0;
    for ( x2 = 0; x2 <= x; x2++ )
    {
      cnk = binomial_coef ( a, x2 );
      dnmk = deranged_enum ( a - x2 );
      sum2 = sum2 + cnk * dnmk;
    }
    cdf = ( double ) ( sum2 ) / i4_factorial ( a );
  }

  return cdf;
}
//****************************************************************************80

int deranged_cdf_inv ( double cdf, int a )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_CDF_INV inverts the Deranged CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, int A, the number of items.
//    1 <= A.
//
//    Output, int X, the corresponding argument.
//
{
  double cdf2;
  double pdf;
  int x;
  int x2;
//
  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "DERANGED_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  cdf2 = 0.0;

  for ( x2 = 0; x2 <= a; x2++ )
  {
    pdf = deranged_pdf ( x2, a );

    cdf2 = cdf2 + pdf;

    if ( cdf <= cdf2 )
    {
      x = x2;
      return x;
    }

  }

  x = a;

  return x;
}
//****************************************************************************80

bool deranged_check ( int a )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_CHECK checks the parameter of the Deranged PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the total number of items.
//    1 <= A.
//
//    Output, bool DERANGED_CHECK, is true if the parameters are legal.
//
{
  if ( a < 1 )
  {
    cout << " \n";
    cout << "DERANGED_CHECK - Fatal error!\n";
    cout << "  A < 1.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

int deranged_enum ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_ENUM returns the number of derangements of N objects.
//
//  Discussion:
//
//    A derangement of N objects is a permutation with no fixed
//    points.  If we symbolize the permutation operation by "P",
//    then for a derangment, P(I) is never equal to I.
//
//  Recursion:
//
//      D(0) = 1
//      D(1) = 0
//      D(2) = 1
//      D(N) = (N-1) * ( D(N-1) + D(N-2) )
//
//    or
//
//      D(0) = 1
//      D(1) = 0
//      D(N) = N * D(N-1) + (-1)**N
//
//    D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
//
//    Based on the inclusion/exclusion law.
//
//    D(N) is the number of ways of placing N non-attacking rooks on
//    an N by N chessboard with one diagonal deleted.
//
//    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
//
//    The number of permutations with exactly K items in the right
//    place is COMB(N,K) * D(N-K).
//
//  First values:
//
//     N         D(N)
//     0           1
//     1           0
//     2           1
//     3           2
//     4           9
//     5          44
//     6         265
//     7        1854
//     8       14833
//     9      133496
//    10     1334961
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Output, int DERANGED_ENUM, the number of derangements of N objects.
//
{
  int dn;
  int dnm1;
  int dnm2;
  int i;

  if ( n < 0 )
  {
    dn = 0;
  }
  else if ( n == 0 )
  {
    dn = 1;
  }
  else if ( n == 1 )
  {
    dn = 0;
  }
  else if ( n == 2 )
  {
    dn = 1;
  }
  else
  {
    dnm1 = 0;
    dn = 1;

    for ( i = 3; i <= n; i++ )
    {
      dnm2 = dnm1;
      dnm1 = dn;
      dn = ( i - 1 ) * ( dnm1 + dnm2 );
    }

  }

  return dn;
}
//****************************************************************************80

double deranged_mean ( int a )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_MEAN returns the mean of the Deranged CDF.
//
//  Discussion:
//
//    The mean is computed by straightforward summation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the number of items.
//    1 <= A.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;
  double pdf;
  int x;

  mean = 0.0;
  for ( x = 1; x <= a; x++ )
  {
    pdf = deranged_pdf ( x, a );
    mean = mean + pdf * x;
  }

  return mean;
}
//****************************************************************************80

double deranged_pdf ( int x, int a )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_PDF evaluates the Deranged PDF.
//
//  Discussion:
//
//    PDF(A;X) is the probability that exactly X items will occur in
//    their proper place after a random permutation of A items.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the number of items in their correct places.
//    0 <= X <= A.
//
//    Input, int A, the total number of items.
//    1 <= A.
//
//    Output, double PDF, the value of the PDF.
//
{
  int cnk;
  int dnmk;
  double pdf;

  if ( x < 0 || a < x )
  {
    pdf = 0.0;
  }
  else
  {
    cnk = binomial_coef ( a, x );
    dnmk = deranged_enum ( a-x );
    pdf = ( double ) ( cnk * dnmk ) / i4_factorial ( a );
  }

  return pdf;
}
//****************************************************************************80

int deranged_sample ( int a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_SAMPLE samples the Deranged PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the number of items.
//    1 <= A.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int X, a sample of the PDF.
//
{
  double cdf;
  int x;

  cdf = r8_uniform_01 ( seed );

  x = deranged_cdf_inv ( cdf, a );

  return x;
}
//****************************************************************************80

double deranged_variance ( int a )

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_VARIANCE returns the variance of the Deranged CDF.
//
//  Discussion:
//
//    The variance is computed by straightforward summation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the number of items.
//    1 <= A.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double mean;
  double pdf;
  int x;
  double variance;

  mean = deranged_mean ( a );

  variance = 0.0;
  for ( x = 1; x <= a; x++ )
  {
    pdf = deranged_pdf ( x, a );
    variance = variance + pdf * pow ( ( x - mean ), 2 );
  }

  return variance;
}
//****************************************************************************80

double digamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    DIGAMMA calculates the digamma or Psi function.
//
//  Discussion:
//
//    DiGamma ( X ) = d ( log ( Gamma ( X ) ) ) / dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    Original FORTRAN77 by J Bernardo.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    J Bernardo,
//    Algorithm AS 103:
//    Psi ( Digamma ) Function,
//    Applied Statistics,
//    Volume 25, Number 3, pages 315-317, 1976.
//
//  Parameters:
//
//    Input, double X, the argument of the digamma function.
//    0 < X.
//
//    Output, double DIGAMMA, the value of the digamma function at X.
//
{
  double c = 8.5;
  double d1 = -0.5772156649;
  double r;
  double s = 0.00001;
  double s3 = 0.08333333333;
  double s4 = 0.0083333333333;
  double s5 = 0.003968253968;
  double value;
  double y;
//
//  The argument must be positive.
//
  if ( x <= 0.0 )
  {
    value = 0.0;
    cout << " \n";
    cout << "DIGAMMA - Fatal error!\n";
    cout << "  X <= 0.\n";
    exit ( 1 );
  }
//
//  Use approximation if argument <= S.
//
  else if ( x <= s )
  {
    value = d1 - 1.0 / x;
  }
//
//  Reduce the argument to DIGAMMA(X + N) where C <= (X + N).
//
  else
  {
    value = 0.0;
    y = x;

    while ( y < c )
    {
      value = value - 1.0 / y;
      y = y + 1.0;
    }
//
//  Use Stirling's (actually de Moivre's) expansion if C < argument.
//
    r = 1.0 / ( y * y );
    value = value + log ( y ) - 0.5 / y - r * ( s3 - r * ( s4 - r * s5 ) );
  }

  return value;
}
//****************************************************************************80

double dipole_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    DIPOLE_CDF evaluates the Dipole CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
//    is interesting, and -1.0 <= B <= 1.0.
//
//    Output, double DIPOLE_CDF, the value of the CDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;

  cdf = 0.5 + ( 1.0 / pi ) * atan ( x ) + b * b * ( x * cos ( 2.0 * a )
    - sin ( 2.0 * a ) ) / ( pi * ( 1.0 + x * x ) );

  return cdf;
}
//****************************************************************************80

double dipole_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    DIPOLE_CDF_INV inverts the Dipole CDF.
//
//  Discussion:
//
//    A simple bisection method is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    -1.0 <= B <= 1.0.
//
//    Output, double DIPOLE_CDF_INV, the corresponding argument of the CDF.
//
{
  double cdf1;
  double cdf2;
  double cdf3;
  int it;
  int it_max = 100;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "DIPOLE_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = - r8_huge ( );
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = r8_huge ( );
    return x;
  }
//
//  Seek X1 < X < X2.
//
  x1 = - 1.0;

  for ( ; ; )
  {
    cdf1 = dipole_cdf ( x1, a, b );

    if ( cdf1 <= cdf )
    {
      break;
    }
    x1 = 2.0 * x1;
  }

  x2 = 1.0;

  for ( ; ; )
  {
    cdf2 = dipole_cdf ( x2, a, b );

    if ( cdf <= cdf2 )
    {
      break;
    }
    x2 = 2.0 * x2;
  }
//
//  Now use bisection.
//
  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = dipole_cdf ( x3, a, b );

    if ( r8_abs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      break;
    }

    if ( it_max < it )
    {
      cout << " \n";
      cout << "DIPOLE_CDF_INV - Fatal error!\n";
      cout << "  Iteration limit exceeded.\n";
      exit ( 1 );
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
      cdf2 = cdf3;
    }

  }

  return x;
}
//****************************************************************************80

bool dipole_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    DIPOLE_CHECK checks the parameters of the Dipole CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
//    is interesting, and -1.0 <= B <= 1.0.
//
//    Output, bool DIPOLE_CHECK, is true if the parameters are legal.
//
{
  if ( b < -1.0 || 1.0 < b )
  {
    cout << " \n";
    cout << "DIPOLE_CHECK - Fatal error!\n";
    cout << "  -1.0 <= B <= 1.0 is required.\n";
    cout << "  The input B = " << b << "\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double dipole_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    DIPOLE_PDF evaluates the Dipole PDF.
//
//  Discussion:
//
//    PDF(A,B;X) =
//        1 / ( PI * ( 1 + X^2 ) )
//      + B^2 * ( ( 1 - X^2 ) * cos ( 2 * A ) + 2.0 * X * sin ( 2 * A ) )
//      / ( PI * ( 1 + X )^2 )
//
//    Densities of this kind commonly occur in the analysis of resonant
//    scattering of elementary particles.
//
//    DIPOLE_PDF(A,0;X) = CAUCHY_PDF(A;X)
//    A = 0, B = 1 yields the single channel dipole distribution.
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
//  Reference:
//
//    Robert Knop,
//    Algorithm 441,
//    ACM Transactions on Mathematical Software.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
//      is interesting,
//    and -1.0 <= B <= 1.0.
//
//    Output, double DIPOLE_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  pdf = 1.0 / ( pi * ( 1.0 + x * x ) )
    + b * b * ( ( 1.0 - x * x ) * cos ( 2.0 * a )
    + 2.0 * x * sin ( 2.0 * x ) )
    / ( pi * ( 1.0 + x * x ) * ( 1.0 + x * x ) );

  return pdf;
}
//****************************************************************************80

double dipole_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    DIPOLE_SAMPLE samples the Dipole PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Knop,
//    Algorithm 441,
//    ACM Transactions on Mathematical Software.
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
//      is interesting,
//    and -1.0 <= B <= 1.0.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double a2;
  double b2;
  double c2;
  double x;
  double *xc;
//
//  Find (X1,X2) at random in a circle.
//
  a2 = b * sin ( a );
  b2 = b * cos ( a );
  c2 = 1.0;

  xc = circle_sample ( a2, b2, c2, seed );
//
//  The dipole variate is the ratio X1 / X2.
//
  x = xc[0] / xc[1];

  delete [] xc;

  return x;
}
//****************************************************************************80

bool dirichlet_check ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_CHECK checks the parameters of the Dirichlet PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, double A[N], the probabilities for each component.
//    Each A[I] should be positive.
//
//    Output, bool DIRICHLET_CHECK, is true if the parameters are legal.
//
{
  int i;
  bool positive;

  positive = false;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] <= 0.0 )
    {
      cout << " \n";
      cout << "DIRICHLET_CHECK - Fatal error!\n";
      cout << "  A[" << i << "] <= 0.\n";
      cout << "  A[" << i << "] = " << a[i] << ".\n";
      return false;
    }
    else if ( 0.0 < a[i] )
    {
      positive = true;
    }
  }

  if ( !positive )
  {
    cout << " \n";
    cout << "DIRICHLET_CHECK - Fatal error!\n";
    cout << "  All parameters are zero!\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double *dirichlet_mean ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MEAN returns the means of the Dirichlet PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, double A[N], the probabilities for each component.
//    Each A[I] should be positive.
//
//    Output, double DIRICHLET_MEAN[N], the means of the PDF.
//
{
  int i;
  double *mean;

  mean = new double[n];

  for ( i = 0; i < n; i++ )
  {
    mean[i] = a[i];
  }

  r8vec_unit_sum ( n, mean );

  return mean;
}
//****************************************************************************80

double *dirichlet_moment2 ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MOMENT2 returns the second moments of the Dirichlet PDF.
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
//    Input, int N, the number of components.
//
//    Input, double A[N], the probabilities for each component.
//    Each A(I) should be positive.
//
//    Output, double DIRICHLET_MOMENT[N*N], the second moments of the PDF.
//
{
  double a_sum;
  double *m2;
  int i;
  int j;

  m2 = new double[n*n];

  a_sum = r8vec_sum ( n, a );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        m2[i+j*n] = a[i] * ( a[i] + 1.0 ) / ( a_sum * ( a_sum + 1.0 ) );
      }
      else
      {
        m2[i+j*n] = a[i] * a[j] / ( a_sum * ( a_sum + 1.0 ) );
      }
    }
  }

  return m2;
}
//****************************************************************************80

double dirichlet_pdf ( double x[], int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_PDF evaluates the Dirichlet PDF.
//
//  Discussion:
//
//    PDF(N,A;X) = Product ( 1 <= I <= N ) X(I)**( A(I) - 1 )
//      * Gamma ( A_SUM ) / A_PROD
//
//    where
//
//      0 < A(I) for all I;
//      0 <= X(I) for all I;
//      Sum ( 1 <= I <= N ) X(I) = 1;
//      A_SUM = Sum ( 1 <= I <= N ) A(I).
//      A_PROD = Product ( 1 <= I <= N ) Gamma ( A(I) )
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
//    Input, double X(N), the argument of the PDF.  Each X(I) should
//    be greater than 0.0, and the X(I)'s must add up to 1.0.
//
//    Input, int N, the number of components.
//
//    Input, double A(N), the probabilities for each component.
//    Each A(I) should be positive.
//
//    Output, double PDF, the value of the PDF.
//
{
  double a_prod;
  double a_sum;
  int i;
  double pdf;
  double tol = 0.0001;
  double x_sum;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] <= 0.0 )
    {
      cout << " \n";
      cout << "DIRICHLET_PDF - Fatal error!\n";
      cout << "  X(I) <= 0.\n";
      exit ( 1 );
    }
  }

  x_sum = r8vec_sum ( n, x );

  if ( tol < r8_abs ( x_sum - 1.0 ) )
  {
    cout << " \n";
    cout << "DIRICHLET_PDF - Fatal error!\n";
    cout << "  SUM X(I) =/= 1.\n";
    exit ( 1 );
  }

  a_sum = r8vec_sum ( n, a );

  a_prod = 1.0;
  for ( i = 0; i < n; i++ )
  {
    a_prod = a_prod * r8_gamma ( a[i] );
  }

  pdf = r8_gamma ( a_sum ) / a_prod;

  for ( i = 0; i < n; i++ )
  {
    pdf = pdf * pow ( x[i], a[i] - 1.0 );
  }

  return pdf;
}
//****************************************************************************80

double *dirichlet_sample ( int n, double a[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_SAMPLE samples the Dirichlet PDF.
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
//  Reference:
//
//    Jerry Banks, editor,
//    Handbook of Simulation,
//    Engineering and Management Press Books, 1998, page 169.
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, double A(N), the probabilities for each component.
//    Each A(I) should be positive.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double DIRICHLET_SAMPLE[N], a sample of the PDF.  The entries
//    of X should sum to 1.
//
{
  double a2;
  double b2;
  double c2;
  int i;
  double *x;

  x = new double[n];

  a2 = 0.0;
  b2 = 1.0;

  for ( i = 0; i < n; i++ )
  {
    c2 = a[i];
    x[i] = gamma_sample ( a2, b2, c2, seed );
  }
//
//  Rescale the vector to have unit sum.
//
  r8vec_unit_sum ( n, x );

  return x;
}
//****************************************************************************80

double *dirichlet_variance ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_VARIANCE returns the variances of the Dirichlet PDF.
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
//    Input, int N, the number of components.
//
//    Input, double A(N), the probabilities for each component.
//    Each A(I) should be positive.
//
//    Output, double DIRICHLET_VARIANCE(N), the variances of the PDF.
//
{
  double a_sum;
  int i;
  double *variance;

  variance = new double[n];

  a_sum = r8vec_sum ( n, a );

  for ( i = 0; i < n; i++ )
  {
    variance[i] = a[i] * ( a_sum - a[i] ) / ( a_sum * a_sum * ( a_sum + 1.0 ) );
  }

  return variance;
}
//****************************************************************************80

bool dirichlet_mix_check ( int comp_num, int elem_num, double a[],
  double comp_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_CHECK checks the parameters of a Dirichlet mixture PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int COMP_NUM, the number of components in the Dirichlet
//    mixture density, that is, the number of distinct Dirichlet PDF's
//    that are mixed together.
//
//    Input, int ELEM_NUM, the number of elements of an observation.
//
//    Input, double A[ELEM_NUM*COMP_NUM], the probabilities
//    for element ELEM_NUM in component COMP_NUM.
//    Each A[I,J] should be positive.
//
//    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of the densities.
//    These do not need to be normalized.  The weight of a given component is
//    the relative probability that that component will be used to generate
//    the sample.
//
//    Output, bool DIRICHLET_MIX_CHECK, is true if the parameters are legal.
//
{
  int comp_i;
  int elem_i;
  bool positive;

  for ( comp_i = 0; comp_i < comp_num; comp_i++ )
  {
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      if ( a[elem_i+comp_i*elem_num] <= 0.0 )
      {
        cout << " \n";
        cout << "DIRICHLET_MIX_CHECK - Fatal error!\n";
        cout << "  A(ELEM,COMP) <= 0.\n";
        cout << "  COMP = " << comp_i << "\n";
        cout << "  ELEM = " << elem_i << "\n";
        cout << "  A(COMP,ELEM) = " << a[elem_i+comp_i*elem_num] << "\n";
        return false;
      }
    }
  }

  positive = false;

  for ( comp_i = 0; comp_i < comp_num; comp_i++ )
  {
    if ( comp_weight[comp_i] < 0.0 )
    {
      cout << " \n";
      cout << "DIRICHLET_MIX_CHECK - Fatal error!\n";
      cout << "  COMP_WEIGHT(COMP) < 0.\n";
      cout << "  COMP = " << comp_i << "\n";
      cout << "  COMP_WEIGHT(COMP) = " << comp_weight[comp_i] << "\n";
      return false;
    }
    else if ( 0.0 < comp_weight[comp_i] )
    {
      positive = true;
    }
  }

  if ( !positive )
  {
    cout << " \n";
    cout << "DIRICHLET_MIX_CHECK - Fatal error!\n";
    cout << "  All component weights are zero.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double *dirichlet_mix_mean ( int comp_num, int elem_num, double a[],
  double comp_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_MEAN returns the means of a Dirichlet mixture PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int COMP_NUM, the number of components in the Dirichlet
//    mixture density, that is, the number of distinct Dirichlet PDF's
//    that are mixed together.
//
//    Input, int ELEM_NUM, the number of elements of an observation.
//
//    Input, double A[ELEM_NUM*COMP_NUM], the probabilities for
//    element ELEM_NUM in component COMP_NUM.
//    Each A[I,J] should be positive.
//
//    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of the densities.
//    These do not need to be normalized.  The weight of a given component is
//    the relative probability that that component will be used to generate
//    the sample.
//
//    Output, double DIRICHLET_MIX_MEAN[ELEM_NUM], the means for each element.
//
{
  double *a_vec;
  int comp_i;
  double *comp_mean;
  double comp_weight_sum;
  int elem_i;
  double *mean;

  a_vec = new double[elem_num];
  mean = new double[elem_num];

  comp_weight_sum = 0.0;
  for ( comp_i = 0; comp_i < comp_num; comp_i++ )
  {
    comp_weight_sum = comp_weight_sum + comp_weight[comp_i];
  }

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    mean[elem_i] = 0.0;
  }

  for ( comp_i = 0; comp_i < comp_num; comp_i++ )
  {
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      a_vec[elem_i] = a[elem_i+comp_i*elem_num];
    }
    comp_mean = dirichlet_mean ( elem_num, a_vec );
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      mean[elem_i] = mean[elem_i] + comp_weight[comp_i] * comp_mean[elem_i];
    }
    delete [] comp_mean;
  }

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    mean[elem_i] = mean[elem_i] / comp_weight_sum;
  }

  delete [] a_vec;

  return mean;
}
//****************************************************************************80

double dirichlet_mix_pdf ( double x[], int comp_num, int elem_num, double a[],
  double comp_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_PDF evaluates a Dirichlet mixture PDF.
//
//  Discussion:
//
//    The PDF is a weighted sum of Dirichlet PDF's.  Each PDF is a
//    "component", with an associated weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X[ELEM_NUM], the argument of the PDF.
//
//    Input, int COMP_NUM, the number of components in the Dirichlet
//    mixture density, that is, the number of distinct Dirichlet PDF's
//    that are mixed together.
//
//    Input, int ELEM_NUM, the number of elements of an observation.
//
//    Input, double A[ELEM_NUM*COMP_NUM], the probabilities for
//    element ELEM_NUM in component COMP_NUM.
//    Each A[I,J] should be positive.
//
//    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of the densities.
//    These do not need to be normalized.  The weight of a given component is
//    the relative probability that that component will be used to generate
//    the sample.
//
//    Output, double DIRICHLET_MIX_PDF, the value of the PDF.
//
{
  double *a_vec;
  double comp_pdf;
  double comp_weight_sum;
  int i;
  int j;
  double pdf;

  a_vec = new double[elem_num];

  comp_weight_sum = 0.0;
  for ( j = 0; j < comp_num; j++ )
  {
    comp_weight_sum = comp_weight_sum + comp_weight[j];
  }

  pdf = 0.0;
  for ( j = 0; j < comp_num; j++ )
  {
    for ( i = 0; i < elem_num; i++ )
    {
      a_vec[i] = a[i+j*elem_num];
    }
    comp_pdf = dirichlet_pdf ( x, elem_num, a_vec );

    pdf = pdf + comp_weight[j] * comp_pdf / comp_weight_sum;
  }

  delete [] a_vec;

  return pdf;
}
//****************************************************************************80

double *dirichlet_mix_sample ( int comp_num, int elem_num, double a[],
  double comp_weight[], int *seed, int *comp )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MIX_SAMPLE samples a Dirichlet mixture PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int COMP_NUM, the number of components in the Dirichlet
//    mixture density, that is, the number of distinct Dirichlet PDF's
//    that are mixed together.
//
//    Input, int ELEM_NUM, the number of elements of an observation.
//
//    Input, double A[ELEM_NUM*COMP_NUM], the probabilities for
//    element ELEM_NUM in component COMP_NUM.
//    Each A[I,J] should be positive.
//
//    Input, double COMP_WEIGHT[COMP_NUM], the mixture weights of the densities.
//    These do not need to be normalized.  The weight of a given component is
//    the relative probability that that component will be used to generate
//    the sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int *COMP, the index of the component of the Dirichlet
//    mixture that was chosen to generate the sample.
//
//    Output, double DIRICHLET_MIX_SAMPLE[ELEM_NUM], a sample of the PDF.
//
{
  double *a_vec;
  int elem_i;
  double *x;

  a_vec = new double[elem_num];
//
//  Choose a particular density component COMP.
//
  *comp = discrete_sample ( comp_num, comp_weight, seed );
//
//  Sample the density number COMP.
//
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    a_vec[elem_i] = a[elem_i+(*comp-1)*elem_num];
  }

  x = dirichlet_sample ( elem_num, a_vec, seed );

  delete [] a_vec;

  return x;
}
//****************************************************************************80

double dirichlet_multinomial_pdf ( int x[], int a, int b, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_MULTINOMIAL_PDF evaluates a Dirichlet Multinomial PDF.
//
//  Formula:
//
//    PDF(A,B,C;X) = Comb(A,B,X) * ( Gamma(C_Sum) / Gamma(C_Sum+A) )
//      Product ( 1 <= I <= B ) Gamma(C(I)+X(I)) / Gamma(C(I))
//
//    where:
//
//      Comb(A,B,X) is the multinomial coefficient C( A; X(1), X(2), ..., X(B) ),
//      C_Sum = Sum ( 1 <= I <= B ) C(I)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kenneth Lange,
//    Mathematical and Statistical Methods for Genetic Analysis,
//    Springer, 1997, page 45.
//
//  Parameters:
//
//    Input, int X[B]; X[I] counts the number of occurrences of
//    outcome I, out of the total of A trials.
//
//    Input, int A, the total number of trials.
//
//    Input, int B, the number of different possible outcomes on
//    one trial.
//
//    Input, double C[B]; C[I] is the Dirichlet parameter associated
//    with outcome I.
//
//    Output, double DIRICHLET_MULTINOMIAL_PDF, the value of the Dirichlet
//     multinomial PDF.
//
{
  double c_sum;
  int i;
  double pdf;
  double pdf_log;

  c_sum = r8vec_sum ( b, c );

  pdf_log = - gamma_log ( c_sum + ( double ) ( a ) ) + gamma_log ( c_sum )
    + gamma_log ( ( double ) ( a + 1 ) );

  for ( i = 0; i < b; i++ )
  {
    pdf_log = pdf_log + gamma_log ( c[i] + ( double ) ( x[i] ) )
      - gamma_log ( c[i] ) - gamma_log ( ( double ) ( x[i] + 1 ) );
  }

  pdf = exp ( pdf_log );

  return pdf;
}
//****************************************************************************80

double discrete_cdf ( int x, int a, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_CDF evaluates the Discrete CDF.
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
//    Input, int X, the item whose probability is desired.
//
//    Input, int A, the number of probabilities assigned.
//
//    Input, double B[A], the relative probabilities of outcomes
//    1 through A.  Each entry must be nonnegative.
//
//    Output, double DISCRETE_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x < 1 )
  {
    cdf = 0.0;
  }
  else if ( x < a )
  {
    cdf = r8vec_sum ( x, b ) / r8vec_sum ( a, b );
  }
  else if ( a <= x )
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

int discrete_cdf_inv ( double cdf, int a, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_CDF_INV inverts the Discrete CDF.
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
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, int A, the number of probabilities assigned.
//
//    Input, double B[A], the relative probabilities of outcomes
//    1 through A.  Each entry must be nonnegative.
//
//    Output, int DISCRETE_CDF_INV, the corresponding argument for which
//    CDF(X-1) < CDF <= CDF(X)
//
{
  double b_sum;
  double cum;
  int j;
  int x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "DISCRETE_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  b_sum = r8vec_sum ( a, b );

  cum = 0.0;

  for ( j = 1; j <= a; j++ )
  {
    cum = cum + b[j-1] / b_sum;

    if ( cdf <= cum )
    {
      x = j;
      return x;
    }
  }

  x = a;

  return x;
}
//****************************************************************************80

bool discrete_check ( int a, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_CHECK checks the parameters of the Discrete CDF.
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
//    Input, int A, the number of probabilities assigned.
//
//    Input, double B[A], the relative probabilities of
//    outcomes 1 through A.  Each entry must be nonnegative.
//
//    Output, bool DISCRETE_CHECK, is true if the parameters are legal.
//
{
  double b_sum;
  int j;

  for ( j = 0; j < a; j++ )
  {
    if ( b[j] < 0.0 )
    {
      cout << " \n";
      cout << "DISCRETE_CHECK - Fatal error!\n";
      cout << "  Negative probabilities not allowed.\n";
      return false;
    }
  }

  b_sum = r8vec_sum ( a, b );

  if ( b_sum == 0.0 )
  {
    cout << " \n";
    cout << "DISCRETE_CHECK - Fatal error!\n";
    cout << "  Total probablity is zero.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double discrete_mean ( int a, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_MEAN evaluates the mean of the Discrete PDF.
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
//    Input, int A, the number of probabilities assigned.
//
//    Input, double B[A], the relative probabilities of
//    outcomes 1 through A.  Each entry must be nonnegative.
//
//    Output, double DISCRETE_MEAN, the mean of the PDF.
//
{
  double b_sum;
  int j;
  double mean;

  b_sum = r8vec_sum ( a, b );

  mean = 0.0;
  for ( j = 0; j < a; j++ )
  {
    mean = mean + ( double ) ( j + 1 ) * b[j];
  }

  mean = mean / b_sum;

  return mean;
}
//****************************************************************************80

double discrete_pdf ( int x, int a, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_PDF evaluates the Discrete PDF.
//
//  Discusion:
//
//    PDF(A,B;X) = B(X) if 1 <= X <= A
//                = 0    otherwise
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
//    Input, int X, the item whose probability is desired.
//
//    Input, int A, the number of probabilities assigned.
//
//    Input, double B[A], the relative probabilities of
//    outcomes 1 through A.  Each entry must be nonnegative.
//
//    Output, double DISCRETE_PDF, the value of the PDF.
//
{
  double b_sum;
  double pdf;

  b_sum = r8vec_sum ( a, b );

  if ( 1 <= x && x <= a )
  {
    pdf = b[x-1] / b_sum;
  }
  else
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

int discrete_sample ( int a, double b[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_SAMPLE samples the Discrete PDF.
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
//    Input, int A, the number of probabilities assigned.
//
//    Input, double B[A], the relative probabilities of
//    outcomes 1 through A.  Each entry must be nonnegative.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int DISCRETE_SAMPLE, a sample of the PDF.
//
{
  double b_sum;
  double cdf;
  int x;

  b_sum = r8vec_sum ( a, b );

  cdf = r8_uniform_01 ( seed );

  x = discrete_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double discrete_variance ( int a, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_VARIANCE evaluates the variance of the Discrete PDF.
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
//    Input, int A, the number of probabilities assigned.
//
//    Input, double B[A], the relative probabilities of
//    outcomes 1 through A.  Each entry must be nonnegative.
//
//    Output, double DISCRETE_VARIANCE, the variance of the PDF.
//
{
  double b_sum;
  int j;
  double mean;
  double variance;

  b_sum = r8vec_sum ( a, b );

  mean = 0.0;
  for ( j = 1; j <= a; j++ )
  {
    mean = mean + ( double ) ( j ) * b[j-1];
  }

  mean = mean / b_sum;

  variance = 0.0;
  for ( j = 1; j <= a; j++ )
  {
    variance = variance + b[j-1] * pow ( ( double ) j - mean, 2 );
  }

  variance = variance / b_sum;

  return variance;
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
//    17 April 2004
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
  double value = 2.71828182845904523536028747135266249775724709369995;

  return value;
}
//****************************************************************************80

double empirical_discrete_cdf ( double x, int a, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_CDF evaluates the Empirical Discrete CDF.
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
//    Input, double X, the argument of the CDF.
//
//    Input, int A, the number of values.
//    0 < A.
//
//    Input, double B[A], the weights of each value.
//    0 <= B(1:A) and at least one value is nonzero.
//
//    Input, double C[A], the values.
//    The values must be distinct and in ascending order.
//
//    Output, double EMPIRICAL_DISCRETE_CDF, the value of the CDF.
//
{
  double bsum;
  double cdf;
  int i;

  cdf = 0.0;

  bsum = r8vec_sum ( a, b );

  for ( i = 1; i <= a; i++ )
  {
    if ( x < c[i-1] )
    {
      return cdf;
    }
    cdf = cdf + b[i-1] / bsum;
  }

  return cdf;
}
//****************************************************************************80

double empirical_discrete_cdf_inv ( double cdf, int a, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_CDF_INV inverts the Empirical Discrete CDF.
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
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, int A, the number of values.
//    0 < A.
//
//    Input, double B(A), the weights of each value.
//    0 <= B(1:A) and at least one value is nonzero.
//
//    Input, double C(A), the values.
//    The values must be distinct and in ascending order.
//
//    Output, double EMPIRICAL_DISCRETE_CDF_INV, the smallest argument
//    whose CDF is greater than or equal to CDF.
//
{
  double bsum;
  double cdf2;
  int i;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "EMPIRICAL_DISCRETE_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  bsum = r8vec_sum ( a, b );

  x = c[0];
  cdf2 = b[0] / bsum;

  for ( i = 1; i < a; i++ )
  {
    if ( cdf <= cdf2 )
    {
      return x;
    }

    x = c[i];
    cdf2 = cdf2 + b[i] / bsum;
  }

  return x;
}
//****************************************************************************80

bool empirical_discrete_check ( int a, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_CHECK checks the parameters of the Empirical Discrete CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the number of values.
//    0 < A.
//
//    Input, double B[A], the weights of each value.
//    0 <= B(1:A) and at least one value is nonzero.
//
//    Input, double C[A], the values.
//    The values must be distinct and in ascending order.
//
//    Output, bool EMPIRICAL_DISCRETE_CHECK, is true if the parameters
//    are legal.
//
{
  int i;
  int j;
  bool positive;

  if ( a <= 0 )
  {
    cout << " \n";
    cout << "EMPIRICAL_DISCRETE_CHECK - Fatal error!\n";
    cout << "  A must be positive.\n";
    cout << "  Input A = " << a << "\n";
    cout << "  A is the number of weights.\n";
    return false;
  }

  for ( i = 0; i < a; i++ )
  {
    if ( b[i] < 0.0 )
    {
      cout << " \n";
      cout << "EMPIRICAL_DISCRETE_CHECK - Fatal error!\n";
      cout << "  B[" << i << "] < 0.\n";
      cout << "  But all B values must be nonnegative.\n";
      return false;
    }
  }

  positive = false;

  for ( i = 0; i < a; i++ )
  {
    if ( 0.0 < b[i] )
    {
      positive = true;
    }
  }

  if ( !positive )
  {
    cout << " \n";
    cout << "EMPIRICAL_DISCRETE_CHECK - Fatal error!\n";
    cout << "  All B(*) = 0.\n";
    cout << "  But at least one B values must be nonzero.\n";
    return false;
  }

  for ( i = 0; i < a; i++ )
  {
    for ( j = i+1; j < a; j++ )
    {
      if ( c[i] == c[j] )
      {
        cout << " \n";
        cout << "EMPIRICAL_DISCRETE_CHECK - Fatal error!\n";
        cout << "  All values C must be unique.\n";
        cout << "  But at least two values are identical.\n";
        return false;
      }
    }
  }

  for ( i = 0; i < a-1; i++ )
  {
    if ( c[i+1] < c[i] )
    {
      cout << " \n";
      cout << "EMPIRICAL_DISCRETE_CHECK - Fatal error!\n";
      cout << "  The values in C must be in ascending order.\n";
      return false;
    }
  }

  return true;
}
//****************************************************************************80

double empirical_discrete_mean ( int a, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_MEAN returns the mean of the Empirical Discrete PDF.
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
//    Input, int A, the number of values.
//    0 < A.
//
//    Input, double B(A), the weights of each value.
//    0 <= B(1:A) and at least one value is nonzero.
//
//    Input, double C(A), the values.
//    The values must be distinct and in ascending order.
//
//    Output, double EMPIRICAL_DISCRETE_MEAN, the mean of the PDF.
//
{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < a; i++ )
  {
    mean = mean + b[i] * c[i];
  }
  mean = mean / r8vec_sum ( a, b );

  return mean;
}
//****************************************************************************80

double empirical_discrete_pdf ( double x, int a, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_PDF evaluates the Empirical Discrete PDF.
//
//  Discussion:
//
//    A set of A values C(1:A) are assigned nonnegative weights B(1:A),
//    with at least one B nonzero.  The probability of C(I) is the
//    value of B(I) divided by the sum of the weights.
//
//    The C's must be distinct, and given in ascending order.
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
//    Input, double X, the argument of the PDF.
//
//    Input, int A, the number of values.
//    0 < A.
//
//    Input, double B(A), the weights of each value.
//    0 <= B(1:A) and at least one value is nonzero.
//
//    Input, double C(A), the values.
//    The values must be distinct and in ascending order.
//
//    Output, double EMPIRICAL_DISCRETE_PDF, the value of the PDF.
//
{
  int i;
  double pdf;

  for ( i = 0; i <= a; i++ )
  {
    if ( x == c[i] )
    {
      pdf = b[i] / r8vec_sum ( a, b );
      return pdf;
    }
  }

  pdf = 0.0;

  return pdf;
}
//****************************************************************************80

double empirical_discrete_sample ( int a, double b[], double c[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_SAMPLE samples the Empirical Discrete PDF.
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
//    Input, int A, the number of values.
//    0 < A.
//
//    Input, double B(A), the weights of each value.
//    0 <= B(1:A) and at least one value is nonzero.
//
//    Input, double C(A), the values.
//    The values must be distinct and in ascending order.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double EMPIRICAL_DISCRETE_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = empirical_discrete_cdf_inv ( cdf, a, b, c );

  return x;
}
//****************************************************************************80

double empirical_discrete_variance ( int a, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_VARIANCE returns the variance of the Empirical Discrete PDF.
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
//    Input, int A, the number of values.
//    0 < A.
//
//    Input, double B(A), the weights of each value.
//    0 <= B(1:A) and at least one value is nonzero.
//
//    Input, double C(A), the values.
//    The values must be distinct and in ascending order.
//
//    Output, double EMPIRICAL_DISCRETE_VARIANCE, the variance of the PDF.
//
{
  double bsum;
  int i;
  double mean;
  double variance;

  bsum = r8vec_sum ( a, b );

  mean = empirical_discrete_mean ( a, b, c );

  variance = 0.0;

  for ( i = 0; i < a; i++ )
  {
    variance = variance + ( b[i] / bsum ) * pow ( c[i] - mean, 2 );
  }

  return variance;
}
//****************************************************************************80

double english_sentence_length_cdf ( int x )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_CDF evaluates the English Sentence Length CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Input, int X, the sentence length whose CDF is desired.
//
//    Output, double ENGLISH_SENTENCE_LENGTH_CDF, the value of the CDF.
//
{
# define SENTENCE_LENGTH_MAX 79

  double cdf;
  int i;
  double pdf_vec[SENTENCE_LENGTH_MAX] = {
    0.00806,
    0.01370,
    0.01862,
    0.02547,
    0.03043,
    0.03189,
    0.03516,
    0.03545,
    0.03286,
    0.03533,
    0.03562,
    0.03788,
    0.03669,
    0.03751,
    0.03518,
    0.03541,
    0.03434,
    0.03305,
    0.03329,
    0.03103,
    0.02867,
    0.02724,
    0.02647,
    0.02526,
    0.02086,
    0.02178,
    0.02128,
    0.01801,
    0.01690,
    0.01556,
    0.01512,
    0.01326,
    0.01277,
    0.01062,
    0.01051,
    0.00901,
    0.00838,
    0.00764,
    0.00683,
    0.00589,
    0.00624,
    0.00488,
    0.00477,
    0.00406,
    0.00390,
    0.00350,
    0.00318,
    0.00241,
    0.00224,
    0.00220,
    0.00262,
    0.00207,
    0.00174,
    0.00174,
    0.00128,
    0.00121,
    0.00103,
    0.00117,
    0.00124,
    0.00082,
    0.00088,
    0.00061,
    0.00061,
    0.00075,
    0.00063,
    0.00056,
    0.00052,
    0.00057,
    0.00031,
    0.00029,
    0.00021,
    0.00017,
    0.00021,
    0.00034,
    0.00031,
    0.00011,
    0.00011,
    0.00008,
    0.00006 };
  double pdf_sum = 0.99768;

  if ( x < 1 )
  {
    cdf = 0.0;
  }
  else if ( x < SENTENCE_LENGTH_MAX )
  {
    cdf = 0.0;
    for ( i = 0; i < x; i++ )
    {
      cdf = cdf + pdf_vec[i];
    }
    cdf = cdf / pdf_sum;
  }
  else if ( SENTENCE_LENGTH_MAX <= x )
  {
    cdf = 1.0;
  }

  return cdf;
# undef SENTENCE_LENGTH_MAX
}
//****************************************************************************80

int english_sentence_length_cdf_inv ( double cdf )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_CDF_INV inverts the English Sentence Length CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Output, int ENGLISH_SENTENCE_LENGTH_CDF_INV, the corresponding sentence
//    length for which CDF(X-1) < CDF <= CDF(X)
//
{
# define SENTENCE_LENGTH_MAX 79

  double cum;
  int j;
  double pdf_vec[SENTENCE_LENGTH_MAX] = {
    0.00806,
    0.01370,
    0.01862,
    0.02547,
    0.03043,
    0.03189,
    0.03516,
    0.03545,
    0.03286,
    0.03533,
    0.03562,
    0.03788,
    0.03669,
    0.03751,
    0.03518,
    0.03541,
    0.03434,
    0.03305,
    0.03329,
    0.03103,
    0.02867,
    0.02724,
    0.02647,
    0.02526,
    0.02086,
    0.02178,
    0.02128,
    0.01801,
    0.01690,
    0.01556,
    0.01512,
    0.01326,
    0.01277,
    0.01062,
    0.01051,
    0.00901,
    0.00838,
    0.00764,
    0.00683,
    0.00589,
    0.00624,
    0.00488,
    0.00477,
    0.00406,
    0.00390,
    0.00350,
    0.00318,
    0.00241,
    0.00224,
    0.00220,
    0.00262,
    0.00207,
    0.00174,
    0.00174,
    0.00128,
    0.00121,
    0.00103,
    0.00117,
    0.00124,
    0.00082,
    0.00088,
    0.00061,
    0.00061,
    0.00075,
    0.00063,
    0.00056,
    0.00052,
    0.00057,
    0.00031,
    0.00029,
    0.00021,
    0.00017,
    0.00021,
    0.00034,
    0.00031,
    0.00011,
    0.00011,
    0.00008,
    0.00006 };
  double pdf_sum = 0.99768;
  int x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "ENGLISH_SENTENCE_LENGTH_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  cum = 0.0;

  for ( j = 0; j < SENTENCE_LENGTH_MAX; j++ )
  {
    cum = cum + pdf_vec[j];

    if ( cdf <= cum / pdf_sum )
    {
      x = j + 1;
      return x;
    }
  }

  x = SENTENCE_LENGTH_MAX;

  return x;
# undef SENTENCE_LENGTH_MAX
}
//****************************************************************************80

double english_sentence_length_mean ( )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_MEAN evaluates the mean of the English Sentence Length PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Output, double ENGLISH_SENTENCE_LENGTH_MEAN, the mean of the PDF.
//
{
# define SENTENCE_LENGTH_MAX 79

  int j;
  double mean;
  double pdf_vec[SENTENCE_LENGTH_MAX] = {
    0.00806,
    0.01370,
    0.01862,
    0.02547,
    0.03043,
    0.03189,
    0.03516,
    0.03545,
    0.03286,
    0.03533,
    0.03562,
    0.03788,
    0.03669,
    0.03751,
    0.03518,
    0.03541,
    0.03434,
    0.03305,
    0.03329,
    0.03103,
    0.02867,
    0.02724,
    0.02647,
    0.02526,
    0.02086,
    0.02178,
    0.02128,
    0.01801,
    0.01690,
    0.01556,
    0.01512,
    0.01326,
    0.01277,
    0.01062,
    0.01051,
    0.00901,
    0.00838,
    0.00764,
    0.00683,
    0.00589,
    0.00624,
    0.00488,
    0.00477,
    0.00406,
    0.00390,
    0.00350,
    0.00318,
    0.00241,
    0.00224,
    0.00220,
    0.00262,
    0.00207,
    0.00174,
    0.00174,
    0.00128,
    0.00121,
    0.00103,
    0.00117,
    0.00124,
    0.00082,
    0.00088,
    0.00061,
    0.00061,
    0.00075,
    0.00063,
    0.00056,
    0.00052,
    0.00057,
    0.00031,
    0.00029,
    0.00021,
    0.00017,
    0.00021,
    0.00034,
    0.00031,
    0.00011,
    0.00011,
    0.00008,
    0.00006 };
  double pdf_sum = 0.99768;

  mean = 0.0;
  for ( j = 1; j <= SENTENCE_LENGTH_MAX; j++ )
  {
    mean = mean + ( double ) ( j ) * pdf_vec[j-1];
  }

  mean = mean / pdf_sum;

  return mean;
# undef SENTENCE_LENGTH_MAX
}
//****************************************************************************80

double english_sentence_length_pdf ( int x )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_PDF evaluates the English Sentence Length PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = B(X) if 1 <= X <= A
//                = 0    otherwise
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Input, int X, the sentence length whose probability is desired.
//
//    Output, double ENGLISH_SENTENCE_LENGTH_PDF, the value of the PDF.
//
{
# define SENTENCE_LENGTH_MAX 79

  double pdf;
  double pdf_vec[SENTENCE_LENGTH_MAX] = {
    0.00806,
    0.01370,
    0.01862,
    0.02547,
    0.03043,
    0.03189,
    0.03516,
    0.03545,
    0.03286,
    0.03533,
    0.03562,
    0.03788,
    0.03669,
    0.03751,
    0.03518,
    0.03541,
    0.03434,
    0.03305,
    0.03329,
    0.03103,
    0.02867,
    0.02724,
    0.02647,
    0.02526,
    0.02086,
    0.02178,
    0.02128,
    0.01801,
    0.01690,
    0.01556,
    0.01512,
    0.01326,
    0.01277,
    0.01062,
    0.01051,
    0.00901,
    0.00838,
    0.00764,
    0.00683,
    0.00589,
    0.00624,
    0.00488,
    0.00477,
    0.00406,
    0.00390,
    0.00350,
    0.00318,
    0.00241,
    0.00224,
    0.00220,
    0.00262,
    0.00207,
    0.00174,
    0.00174,
    0.00128,
    0.00121,
    0.00103,
    0.00117,
    0.00124,
    0.00082,
    0.00088,
    0.00061,
    0.00061,
    0.00075,
    0.00063,
    0.00056,
    0.00052,
    0.00057,
    0.00031,
    0.00029,
    0.00021,
    0.00017,
    0.00021,
    0.00034,
    0.00031,
    0.00011,
    0.00011,
    0.00008,
    0.00006 };
  double pdf_sum = 0.99768;

  if ( 1 <= x && x <= SENTENCE_LENGTH_MAX )
  {
    pdf = pdf_vec[x-1] / pdf_sum;
  }
  else
  {
    pdf = 0.0;
  }

  return pdf;
# undef SENTENCE_LENGTH_MAX
}
//****************************************************************************80

int english_sentence_length_sample ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_SAMPLE samples the English Sentence Length PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int ENGLISH_SENTENCE_LENGTH_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  int x;

  cdf = r8_uniform_01 ( seed );

  x = english_sentence_length_cdf_inv ( cdf );

  return x;
# undef SENTENCE_LENGTH_MAX
}
//****************************************************************************80

double english_sentence_length_variance ( )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_VARIANCE: variance of the English Sentence Length PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Output, double ENGLISH_SENTENCE_LENGTH_VARIANCE, the variance of the PDF.
//
{
# define SENTENCE_LENGTH_MAX 79

  int j;
  double mean;
  double pdf_vec[SENTENCE_LENGTH_MAX] = {
    0.00806,
    0.01370,
    0.01862,
    0.02547,
    0.03043,
    0.03189,
    0.03516,
    0.03545,
    0.03286,
    0.03533,
    0.03562,
    0.03788,
    0.03669,
    0.03751,
    0.03518,
    0.03541,
    0.03434,
    0.03305,
    0.03329,
    0.03103,
    0.02867,
    0.02724,
    0.02647,
    0.02526,
    0.02086,
    0.02178,
    0.02128,
    0.01801,
    0.01690,
    0.01556,
    0.01512,
    0.01326,
    0.01277,
    0.01062,
    0.01051,
    0.00901,
    0.00838,
    0.00764,
    0.00683,
    0.00589,
    0.00624,
    0.00488,
    0.00477,
    0.00406,
    0.00390,
    0.00350,
    0.00318,
    0.00241,
    0.00224,
    0.00220,
    0.00262,
    0.00207,
    0.00174,
    0.00174,
    0.00128,
    0.00121,
    0.00103,
    0.00117,
    0.00124,
    0.00082,
    0.00088,
    0.00061,
    0.00061,
    0.00075,
    0.00063,
    0.00056,
    0.00052,
    0.00057,
    0.00031,
    0.00029,
    0.00021,
    0.00017,
    0.00021,
    0.00034,
    0.00031,
    0.00011,
    0.00011,
    0.00008,
    0.00006 };
  double pdf_sum = 0.99768;
  double variance;

  mean = 0.0;
  for ( j = 1; j <= SENTENCE_LENGTH_MAX; j++ )
  {
    mean = mean + ( double ) ( j ) * pdf_vec[j-1];
  }

  mean = mean / pdf_sum;

  variance = 0.0;
  for ( j = 1; j <= SENTENCE_LENGTH_MAX; j++ )
  {
    variance = variance + pdf_vec[j-1] * pow ( ( double ) j - mean, 2 );
  }

  variance = variance / pdf_sum;

  return variance;
# undef SENTENCE_LENGTH_MAX
}
//****************************************************************************80

double english_word_length_cdf ( int x )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_CDF evaluates the English Word Length CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Input, int X, the word length whose CDF is desired.
//
//    Output, double ENGLISH_WORD_LENGTH_CDF, the value of the CDF.
//
{
# define WORD_LENGTH_MAX 27

  double cdf;
  int i;
  double pdf_vec[WORD_LENGTH_MAX] = {
    0.03160,
    0.16975,
    0.21192,
    0.15678,
    0.10852,
    0.08524,
    0.07724,
    0.05623,
    0.04032,
    0.02766,
    0.01582,
    0.00917,
    0.00483,
    0.00262,
    0.00099,
    0.00050,
    0.00027,
    0.00022,
    0.00011,
    0.00006,
    0.00005,
    0.00002,
    0.00001,
    0.00001,
    0.00001,
    0.00001,
    0.00001 };
  double pdf_sum = 0.99997;

  if ( x < 1 )
  {
    cdf = 0.0;
  }
  else if ( x < WORD_LENGTH_MAX )
  {
    cdf = 0.0;
    for ( i = 0; i < x; i++ )
    {
      cdf = cdf + pdf_vec[i];
    }
    cdf = cdf / pdf_sum;
  }
  else if ( WORD_LENGTH_MAX <= x )
  {
    cdf = 1.0;
  }

  return cdf;
# undef WORD_LENGTH_MAX
}
//****************************************************************************80

int english_word_length_cdf_inv ( double cdf )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_CDF_INV inverts the English Word Length CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Output, int ENGLISH_WORD_LENGTH_CDF_INV, the corresponding word
//    length for which CDF(X-1) < CDF <= CDF(X)
//
{
# define WORD_LENGTH_MAX 27

  double cum;
  int j;
  double pdf_vec[WORD_LENGTH_MAX] = {
    0.03160,
    0.16975,
    0.21192,
    0.15678,
    0.10852,
    0.08524,
    0.07724,
    0.05623,
    0.04032,
    0.02766,
    0.01582,
    0.00917,
    0.00483,
    0.00262,
    0.00099,
    0.00050,
    0.00027,
    0.00022,
    0.00011,
    0.00006,
    0.00005,
    0.00002,
    0.00001,
    0.00001,
    0.00001,
    0.00001,
    0.00001 };
  double pdf_sum = 0.99997;
  int x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "ENGLISH_WORD_LENGTH_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  cum = 0.0;

  for ( j = 0; j < WORD_LENGTH_MAX; j++ )
  {
    cum = cum + pdf_vec[j];

    if ( cdf <= cum / pdf_sum )
    {
      x = j + 1;
      return x;
    }
  }

  x = WORD_LENGTH_MAX;

  return x;
# undef WORD_LENGTH_MAX
}
//****************************************************************************80

double english_word_length_mean ( )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_MEAN evaluates the mean of the English Word Length PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Output, double ENGLISH_WORD_LENGTH_MEAN, the mean of the PDF.
//
{
# define WORD_LENGTH_MAX 27

  int j;
  double mean;
  double pdf_vec[WORD_LENGTH_MAX] = {
    0.03160,
    0.16975,
    0.21192,
    0.15678,
    0.10852,
    0.08524,
    0.07724,
    0.05623,
    0.04032,
    0.02766,
    0.01582,
    0.00917,
    0.00483,
    0.00262,
    0.00099,
    0.00050,
    0.00027,
    0.00022,
    0.00011,
    0.00006,
    0.00005,
    0.00002,
    0.00001,
    0.00001,
    0.00001,
    0.00001,
    0.00001 };
  double pdf_sum = 0.99997;

  mean = 0.0;
  for ( j = 1; j <= WORD_LENGTH_MAX; j++ )
  {
    mean = mean + ( double ) ( j ) * pdf_vec[j-1];
  }

  mean = mean / pdf_sum;

  return mean;
# undef WORD_LENGTH_MAX
}
//****************************************************************************80

double english_word_length_pdf ( int x )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_PDF evaluates the English Word Length PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = B(X) if 1 <= X <= A
//                = 0    otherwise
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Input, int X, the word length whose probability is desired.
//
//    Output, double ENGLISH_WORD_LENGTH_PDF, the value of the PDF.
//
{
# define WORD_LENGTH_MAX 27

  double pdf;
  double pdf_vec[WORD_LENGTH_MAX] = {
    0.03160,
    0.16975,
    0.21192,
    0.15678,
    0.10852,
    0.08524,
    0.07724,
    0.05623,
    0.04032,
    0.02766,
    0.01582,
    0.00917,
    0.00483,
    0.00262,
    0.00099,
    0.00050,
    0.00027,
    0.00022,
    0.00011,
    0.00006,
    0.00005,
    0.00002,
    0.00001,
    0.00001,
    0.00001,
    0.00001,
    0.00001 };
  double pdf_sum = 0.99997;

  if ( 1 <= x && x <= WORD_LENGTH_MAX )
  {
    pdf = pdf_vec[x-1] / pdf_sum;
  }
  else
  {
    pdf = 0.0;
  }

  return pdf;
# undef WORD_LENGTH_MAX
}
//****************************************************************************80

int english_word_length_sample ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_SAMPLE samples the English Word Length PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int ENGLISH_WORD_LENGTH_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  int x;

  cdf = r8_uniform_01 ( seed );

  x = english_word_length_cdf_inv ( cdf );

  return x;
# undef WORD_LENGTH_MAX
}
//****************************************************************************80

double english_word_length_variance ( )

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_VARIANCE: variance of the English Word Length PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Henry Kucera, Winthrop Francis,
//    Computational Analysis of Present-Day American English,
//    Brown University Press, 1967.
//
//  Parameters:
//
//    Output, double ENGLISH_WORD_LENGTH_VARIANCE, the variance of the PDF.
//
{
# define WORD_LENGTH_MAX 27

  int j;
  double mean;
  double pdf_vec[WORD_LENGTH_MAX] = {
    0.03160,
    0.16975,
    0.21192,
    0.15678,
    0.10852,
    0.08524,
    0.07724,
    0.05623,
    0.04032,
    0.02766,
    0.01582,
    0.00917,
    0.00483,
    0.00262,
    0.00099,
    0.00050,
    0.00027,
    0.00022,
    0.00011,
    0.00006,
    0.00005,
    0.00002,
    0.00001,
    0.00001,
    0.00001,
    0.00001,
    0.00001 };
  double pdf_sum = 0.99997;
  double variance;

  mean = 0.0;
  for ( j = 1; j <= WORD_LENGTH_MAX; j++ )
  {
    mean = mean + ( double ) ( j ) * pdf_vec[j-1];
  }

  mean = mean / pdf_sum;

  variance = 0.0;
  for ( j = 1; j <= WORD_LENGTH_MAX; j++ )
  {
    variance = variance + pdf_vec[j-1] * pow ( ( double ) j - mean, 2 );
  }

  variance = variance / pdf_sum;

  return variance;
# undef WORD_LENGTH_MAX
}
//****************************************************************************80

double erlang_cdf ( double x, double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_CDF evaluates the Erlang CDF.
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
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, int C, the parameters of the PDF.
//    0.0 < B.
//    0 < C.
//
//    Output, double ERLANG_CDF, the value of the CDF.
//
{
  double cdf;
  double p2;
  double x2;

  if ( x < a )
  {
    cdf = 0.0;
  }
  else
  {
    x2 = ( x - a ) / b;
    p2 = ( double ) ( c );

    cdf = gamma_inc ( p2, x2 );
  }

  return cdf;
}
//****************************************************************************80

double erlang_cdf_inv ( double cdf, double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_CDF_INV inverts the Erlang CDF.
//
//  Discussion:
//
//    A simple bisection method is used.
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
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, int C, the parameters of the PDF.
//    0.0 < B.
//    0 < C.
//
//    Output, double ERLANG_CDF_INV, the corresponding argument of the CDF.
//
{
  double cdf1;
  double cdf2;
  double cdf3;
  int it;
  int it_max = 100;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "ERLANG_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = a;
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = r8_huge ( );
    return x;
  }

  x1 = a;
  cdf1 = 0.0;

  x2 = a + 1.0;

  for ( ; ; )
  {
    cdf2 = erlang_cdf ( x2, a, b, c );

    if ( cdf < cdf2 )
    {
      break;
    }
    x2 = a + 2.0 * ( x2 - a );
  }
//
//  Now use bisection.
//
  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = erlang_cdf ( x3, a, b, c );

    if ( r8_abs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      break;
    }

    if ( it_max < it )
    {
      cout << " \n";
      cout << "ERLANG_CDF_INV - Warning!\n";
      cout << "  Iteration limit exceeded.\n";
      return x;
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
      cdf2 = cdf3;
    }
  }
  return x;
}
//****************************************************************************80

bool erlang_check ( double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_CHECK checks the parameters of the Erlang PDF.
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
//    Input, double A, B, int C, the parameters of the PDF.
//    0.0 < B.
//    0 < C.
//
//    Output, bool ERLANG_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "ERLANG_CHECK - Fatal error!\n";
    cout << "  B <= 0.0\n";
    return false;
  }

  if ( c <= 0 )
  {
    cout << " \n";
    cout << "ERLANG_CHECK - Fatal error!\n";
    cout << "  C <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double erlang_mean ( double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_MEAN returns the mean of the Erlang PDF.
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
//    Input, double A, B, int C, the parameters of the PDF.
//    0.0 < B.
//    0 < C.
//
//    Output, double ERLANG_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a + b * ( double ) ( c );

  return mean;
}
//****************************************************************************80

double erlang_pdf ( double x, double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_PDF evaluates the Erlang PDF.
//
//  Discussion:
//
//    PDF(A,B,C;X) = ( ( X - A ) / B )**( C - 1 )
//      / ( B * Gamma ( C ) * EXP ( ( X - A ) / B ) )
//
//    for 0 < B, 0 < C integer, A <= X.
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
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, int C, the parameters of the PDF.
//    0.0 < B.
//    0 < C.
//
//    Output, double ERLANG_PDF, the value of the PDF.
//
{
  double pdf;
  double y;

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else
  {
    y = ( x - a ) / b;

    pdf = pow ( y, c - 1 )
      / ( b * i4_factorial ( c - 1 ) * exp ( y ) );
  }

  return pdf;
}
//****************************************************************************80

double erlang_sample ( double a, double b, int c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_SAMPLE samples the Erlang PDF.
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
//    Input, double A, B, int C, the parameters of the PDF.
//    0.0 < B.
//    0 < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double a2;
  double b2;
  int i;
  double x;
  double x2;

  a2 = 0.0;
  b2 = b;
  x = a;

  for ( i = 1; i <= c; i++ )
  {
    x2 = exponential_sample ( a2, b2, seed );
    x = x + x2;
  }

  return x;
}
//****************************************************************************80

double erlang_variance ( double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_VARIANCE returns the variance of the Erlang PDF.
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
//    Input, double A, B, int C, the parameters of the PDF.
//    0.0 < B.
//    0 < C.
//
//    Output, double ERLANG_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance =  b * b * ( double ) ( c );

  return variance;
}
//****************************************************************************80

double error_f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    ERROR_F evaluates the error function ERF.
//
//  Discussion:
//
//    Since some compilers already supply a routine named ERF which evaluates
//    the error function, this routine has been given a distinct, if
//    somewhat unnatural, name.
//
//    The function is defined by:
//
//      ERF(X) = ( 2 / sqrt ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( - T^2 ) dT.
//
//    Properties of the function include:
//
//      Limit ( X -> -Infinity ) ERF(X) =          -1.0;
//                               ERF(0) =           0.0;
//                               ERF(0.476936...) = 0.5;
//      Limit ( X -> +Infinity ) ERF(X) =          +1.0.
//
//      0.5 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2006
//
//  Author:
//
//    Original FORTRAN77 versino by William Cody.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    "Rational Chebyshev Approximations for the Error Function",
//    Mathematics of Computation,
//    1969, pages 631-638.
//
//  Parameters:
//
//    Input, double X, the argument of the error function.
//
//    Output, double ERROR_F, the value of the error function.
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
  double erfxd;
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
  double xabs;
  double xbig = 26.543;
  double xden;
  double xnum;
  double xsmall = 1.11E-16;
  double xsq;

  xabs = r8_abs ( ( x ) );
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

    erfxd = x * ( xnum + a[3] ) / ( xden + b[3] );
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

    erfxd = ( xnum + c[7] ) / ( xden + d[7] );
    xsq = ( ( double ) ( ( int ) ( xabs * 16.0 ) ) ) / 16.0;
    del = ( xabs - xsq ) * ( xabs + xsq );
    erfxd = exp ( - xsq * xsq ) * exp ( -del ) * erfxd;

    erfxd = ( 0.5 - erfxd ) + 0.5;

    if ( x < 0.0 )
    {
      erfxd = -erfxd;
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
        erfxd = 1.0;
      }
      else
      {
        erfxd = - 1.0;
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

      erfxd = xsq * ( xnum + p[4] ) / ( xden + q[4] );
      erfxd = ( sqrpi - erfxd ) / xabs;
      xsq = ( ( double ) ( ( int ) ( xabs * 16.0 ) ) ) / 16.0;
      del = ( xabs - xsq ) * ( xabs + xsq );
      erfxd = exp ( - xsq * xsq ) * exp ( - del ) * erfxd;

      erfxd = ( 0.5 - erfxd ) + 0.5;

      if ( x < 0.0 )
      {
        erfxd = -erfxd;
      }
    }
  }

  return erfxd;
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
//      Gamma = limit ( M -> Infinity )
//        ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
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
  double value = 0.577215664901532860606512090082402431042;

  return value;
}
//****************************************************************************80

double exponential_01_cdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_CDF evaluates the Exponential 01 CDF.
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
//    Input, double X, the argument of the PDF.
//
//    Output, double EXPONENTIAL_01_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 1.0 - exp ( - x );
  }

  return cdf;
}
//****************************************************************************80

double exponential_01_cdf_inv ( double cdf )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_CDF_INV inverts the Exponential 01 CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Output, double EXPONENTIAL_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "EXPONENTIAL_01_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = - log ( 1.0 - cdf );

  return x;
}
//****************************************************************************80

double exponential_01_mean ( )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_MEAN returns the mean of the Exponential 01 PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double EXPONENTIAL_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = 1.0;

  return mean;
}
//****************************************************************************80

double exponential_01_pdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_PDF evaluates the Exponential 01 PDF.
//
//  Discussion:
//
//    PDF(X) = EXP ( - X )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = exp ( - x );
  }

  return pdf;
}
//****************************************************************************80

double exponential_01_sample ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_SAMPLE samples the Exponential PDF with parameter 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = - log ( 1.0 - cdf );

  return x;
}
//****************************************************************************80

double exponential_01_variance ( )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_VARIANCE returns the variance of the Exponential 01 PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = 1.0;

  return variance;
}
//****************************************************************************80

double exponential_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_CDF evaluates the Exponential CDF.
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
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameter of the PDF.
//    0.0 < B.
//
//    Output, double EXPONENTIAL_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 1.0 - exp ( ( a - x ) / b );
  }

  return cdf;
}
//****************************************************************************80

double exponential_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_CDF_INV inverts the Exponential CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXPONENTIAL_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "EXPONENTIAL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a - b * log ( 1.0 - cdf );

  return x;
}
//****************************************************************************80

void exponential_cdf_values ( int *n_data, double *lambda, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_CDF_VALUES returns some values of the Exponential CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = ExponentialDistribution [ lambda ]
//      CDF [ dist, x ]
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
//    Output, double *LAMBDA, the parameter of the distribution.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 9

  double fx_vec[N_MAX] = {
     0.3934693402873666,
     0.6321205588285577,
     0.7768698398515702,
     0.8646647167633873,
     0.8646647167633873,
     0.9816843611112658,
     0.9975212478233336,
     0.9996645373720975,
     0.9999546000702375 };

  double lambda_vec[N_MAX] = {
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *lambda = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *lambda = lambda_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool exponential_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_CHECK checks the parameters of the Exponential CDF.
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
//    Input, double A, B, the parameter of the PDF.
//    0.0 < B.
//
//    Output, bool EXPONENTIAL_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "EXPONENTIAL_CHECK - Fatal error!\n";
    cout << "  B <= 0.0\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double exponential_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_MEAN returns the mean of the Exponential PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXPONENTIAL_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a + b;

  return mean;
}
//****************************************************************************80

double exponential_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_PDF evaluates the Exponential PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = ( 1 / B ) * EXP ( ( A - X ) / B )
//
//    The time interval between two Poisson events is a random
//    variable with the Exponential PDF.  The parameter B is the
//    average interval between events.
//
//    In another context, the Exponential PDF is related to
//    the Boltzmann distribution, which describes the relative
//    probability of finding a system, which is in thermal equilibrium
//    at absolute temperature T, in a given state having energy E.
//    The relative probability is
//
//      Boltzmann_Relative_Probability(E,T) = exp ( - E / ( k * T ) ),
//
//    where k is the Boltzmann constant,
//
//      k = 1.38 * 10**(-23) joules / degree Kelvin
//
//    and normalization requires a determination of the possible
//    energy states of the system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXPONENTIAL_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < a )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = ( 1.0 / b ) * exp ( ( a - x ) / b );
  }

  return pdf;
}
//****************************************************************************80

double exponential_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_SAMPLE samples the Exponential PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double EXPONENTIAL_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = exponential_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double exponential_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_VARIANCE returns the variance of the Exponential PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXPONENTIAL_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = b * b;

  return variance;
}
//****************************************************************************80

double extreme_values_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_CDF evaluates the Extreme Values CDF.
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
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXTREME_VALUES_CDF, the value of the CDF.
//
{
  double cdf;
  double y;

  y = ( x - a ) / b;

  cdf = exp ( - exp ( - y ) );

  return cdf;
}
//****************************************************************************80

double extreme_values_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_CDF_INV inverts the Extreme Values CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXTREME_VALUES_CDF_INV, the corresponding argument of the CDF.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "EXTREME_VALUES_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a - b * log ( - log ( cdf ) );

  return x;
}
//****************************************************************************80

void extreme_values_cdf_values ( int *n_data, double *alpha, double *beta,
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_CDF_VALUES returns some values of the Extreme Values CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = ExtremeValuesDistribution [ alpha, beta ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
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
//    Output, double *ALPHA, the first parameter of the distribution.
//
//    Output, double *BETA, the second parameter of the distribution.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double alpha_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  double beta_vec[N_MAX] = {
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  double fx_vec[N_MAX] = {
     0.3678794411714423,
     0.8734230184931166,
     0.9818510730616665,
     0.9975243173927525,
     0.5452392118926051,
     0.4884435800065159,
     0.4589560693076638,
     0.4409910259429826,
     0.5452392118926051,
     0.3678794411714423,
     0.1922956455479649,
     0.6598803584531254E-01 };

  double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *alpha = 0.0;
    *beta = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *alpha = alpha_vec[*n_data-1];
    *beta = beta_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool extreme_values_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_CHECK checks the parameters of the Extreme Values CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, bool EXTREME_VALUES_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "EXTREME_VALUES_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double extreme_values_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_MEAN returns the mean of the Extreme Values PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXTREME_VALUES_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a + b * euler_constant ( );

  return mean;
}
//****************************************************************************80

double extreme_values_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_PDF evaluates the Extreme Values PDF.
//
//  Discussion:
//
//    PDF(A,B;X) =
//      ( 1 / B ) * exp ( ( A - X ) / B - exp ( ( A - X ) / B  ) ).
//
//    The Extreme Values PDF is also known as the Fisher-Tippet PDF
//    and the Log-Weibull PDF.
//
//    The special case A = 0 and B = 1 is the Gumbel PDF.
//
//    The Extreme Values PDF is the limiting distribution for the
//    smallest or largest value in a large sample drawn from
//    any of a great variety of distributions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein, editor,
//    CRC Concise Encylopedia of Mathematics,
//    CRC Press, 1998.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXTREME_VALUES_PDF, the value of the PDF.
//
{
  double pdf;

  pdf = ( 1.0 / b ) * exp ( ( a - x ) / b - exp ( ( a - x ) / b ) );

  return pdf;
}
//****************************************************************************80

double extreme_values_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_SAMPLE samples the Extreme Values PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double EXTREME_VALUES_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = extreme_values_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double extreme_values_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_VARIANCE returns the variance of the Extreme Values PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double EXTREME_VALUES_VARIANCE, the variance of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double variance;

  variance = pi * pi * b * b / 6.0;

  return variance;
}
//****************************************************************************80

double f_cdf ( double x, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F_CDF evaluates the F central CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
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
//    Input, double X, the argument of the CDF.
//
//    Input, int M, N, the parameters of the PDF.
//    1 <= M,
//    1 <= N.
//
//    Output, double F_CDF, the value of the CDF.
//
{
  double arg1;
  double arg2;
  double arg3;
  double cdf;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    arg1 = 0.5 * ( double ) ( n );
    arg2 = 0.5 * ( double ) ( m );

    arg3 = ( double ) ( n ) / ( ( double ) ( n ) + ( double ) ( m ) * x );

    cdf = 1.0 - beta_inc ( arg1, arg2, arg3 );
  }

  return cdf;
}
//****************************************************************************80

void f_cdf_values ( int *n_data, int *a, int *b, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    F_CDF_VALUES returns some values of the F CDF test function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = FRatioDistribution [ dfn, dfd ]
//      CDF [ dist, x ]
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
//    Output, int *A, int B, the parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 20

  int a_vec[N_MAX] = {
    1,
    1,
    5,
    1,
    2,
    4,
    1,
    6,
    8,
    1,
    3,
    6,
    1,
    1,
    1,
    1,
    2,
    3,
    4,
    5 };

  int b_vec[N_MAX] = {
     1,
     5,
     1,
     5,
    10,
    20,
     5,
     6,
    16,
     5,
    10,
    12,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5 };

  double fx_vec[N_MAX] = {
     0.5000000000000000,
     0.4999714850534485,
     0.4996034370170990,
     0.7496993658293228,
     0.7504656462757382,
     0.7514156325324275,
     0.8999867031372156,
     0.8997127554259699,
     0.9002845660853669,
     0.9500248817817622,
     0.9500574946122442,
     0.9501926400000000,
     0.9750133887312993,
     0.9900022327445249,
     0.9949977837872073,
     0.9989999621122122,
     0.5687988496283079,
     0.5351452100063650,
     0.5143428032407864,
     0.5000000000000000 };

  double x_vec[N_MAX] = {
      1.00,
      0.528,
      1.89,
      1.69,
      1.60,
      1.47,
      4.06,
      3.05,
      2.09,
      6.61,
      3.71,
      3.00,
     10.01,
     16.26,
     22.78,
     47.18,
      1.00,
      1.00,
      1.00,
      1.00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0;
    *b = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool f_check ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F_CHECK checks the parameters of the F PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the parameters of the PDF.
//    1 <= M,
//    1 <= N.
//
//    Output, bool F_CHECK, is true if the parameters are legal.
//
{
  if ( m <= 0 )
  {
    cout << "\n";
    cout << "F_CHECK - Fatal error!\n";
    cout << "  M <= 0.\n";
    return false;
  }

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "F_CHECK - Fatal error!\n";
    cout << "  N <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double f_mean ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F_MEAN returns the mean of the F central PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the parameters of the PDF.
//    1 <= M,
//    1 <= N.
//    Note, however, that the mean is not defined unless 3 <= N.
//
//    Output, double F_MEAN, the mean of the PDF.
//
{
  double mean;

  if ( n < 3 )
  {
    cout << " \n";
    cout << "F_MEAN - Fatal error!\n";
    cout << "  The mean is not defined for N < 3.\n";
    exit ( 1 );
  }

  mean = ( double ) ( n ) / ( double ) ( n - 2 );

  return mean;
}
//****************************************************************************80

double f_pdf ( double x, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F_PDF evaluates the F central PDF.
//
//  Discussion:
//
//    PDF(M,N;X) = M**(M/2) * X**((M-2)/2)
//      / ( Beta(M/2,N/2) * N**(M/2) * ( 1 + (M/N) * X )**((M+N)/2)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X
//
//    Input, int M, N, the parameters of the PDF.
//    1 <= M,
//    1 <= N.
//
//    Output, double F_PDF, the value of the PDF.
//
{
  double a;
  double b;
  double bot1;
  double bot2;
  double pdf;
  double top;

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    a = ( double ) ( m );
    b = ( double ) ( n );

    top = sqrt ( pow ( a, m ) * pow ( b, n ) * pow ( x, m - 2 ) );
    bot1 = beta ( a / 2.0, b / 2.0 ) ;
    bot2 = sqrt ( pow ( b + a * x, m + n ) );

    pdf = top / ( bot1 * bot2 );

  }

  return pdf;
}
//****************************************************************************80

double f_sample ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    F_SAMPLE samples the F central PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the parameters of the PDF.
//    1 <= M,
//    1 <= N.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double F_SAMPLE, a sample of the PDF.
//
{
  double a;
  double x;
  double xm;
  double xn;

  a = ( double ) ( m );
  xm = chi_square_sample ( a, seed );

  a = ( double ) ( n );
  xn = chi_square_sample ( a, seed );

  x = ( double ) ( n ) * xm / ( ( double ) ( m ) * xn );

  return x;
}
//****************************************************************************80

double f_variance ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F_VARIANCE returns the variance of the F central PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the parameters of the PDF.
//    1 <= M,
//    1 <= N.
//    Note, however, that the variance is not defined unless 5 <= N.
//
//    Output, double F_VARIANCE, the variance of the PDF.
//
{
  double variance;

  if ( n < 5 )
  {
    cout << " \n";
    cout << "F_VARIANCE - Fatal error!\n";
    cout << "  The variance is not defined for N < 5.\n";
    exit ( 1 );
  }

  variance = ( double ) ( 2 * n * n * ( m + n - 2 ) ) /
    ( double ) ( m * ( n - 2 ) * ( n - 2 ) * ( n - 4 ) );

  return variance;
}
//****************************************************************************80

void f_noncentral_cdf_values ( int *n_data, int *n1, int *n2, double *lambda,
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NoncentralFRatioDistribution [ n1, n2, lambda ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2004
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
//    Output, int *N1, int *N2, the numerator and denominator
//    degrees of freedom.
//
//    Output, double *LAMBDA, the noncentrality parameter.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 22

  double fx_vec[N_MAX] = {
     0.5000000000000000,
     0.6367825323508774,
     0.5840916116305482,
     0.3234431872392788,
     0.4501187879813550,
     0.6078881441188312,
     0.7059275551414605,
     0.7721782003263727,
     0.8191049017635072,
     0.3170348430749965,
     0.4327218008454471,
     0.4502696915707327,
     0.4261881186594096,
     0.6753687206341544,
     0.4229108778879005,
     0.6927667261228938,
     0.3632174676491226,
     0.4210054012695865,
     0.4266672258818927,
     0.4464016600524644,
     0.8445888579504827,
     0.4339300273343604 };

  double lambda_vec[N_MAX] = {
     0.00,
     0.00,
     0.25,
     1.00,
     1.00,
     1.00,
     1.00,
     1.00,
     1.00,
     2.00,
     1.00,
     1.00,
     1.00,
     2.00,
     1.00,
     1.00,
     0.00,
     1.00,
     1.00,
     1.00,
     1.00,
     1.00 };

  int n1_vec[N_MAX] = {
     1,  1,  1,  1,
     1,  1,  1,  1,
     1,  1,  2,  2,
     3,  3,  4,  4,
     5,  5,  6,  6,
     8, 16 };

  int n2_vec[N_MAX] = {
     1,  5,  5,  5,
     5,  5,  5,  5,
     5,  5,  5, 10,
     5,  5,  5,  5,
     1,  5,  6, 12,
    16,  8 };

  double x_vec[N_MAX] = {
     1.00,
     1.00,
     1.00,
     0.50,
     1.00,
     2.00,
     3.00,
     4.00,
     5.00,
     1.00,
     1.00,
     1.00,
     1.00,
     1.00,
     1.00,
     2.00,
     1.00,
     1.00,
     1.00,
     1.00,
     2.00,
     2.00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n1 = 0;
    *n2 = 0;
    *lambda = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *n1 = n1_vec[*n_data-1];
    *n2 = n2_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool f_noncentral_check ( double a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F_NONCENTRAL_CHECK checks the parameters of the F noncentral PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, a parameter of the PDF.
//
//    Input, int M, N, the parameters of the PDF.
//    1 <= M,
//    1 <= N.
//
//    Output, bool F_NONCENTRAL_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << " \n";
    cout << "F_NONCENTRAL_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    exit ( 1 );
  }

  if ( m <= 0 )
  {
    cout << "\n";
    cout << "F_NONCENTRAL_CHECK - Fatal error!\n";
    cout << "  M <= 0.\n";
    return false;
  }

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "F_NONCENTRAL_CHECK - Fatal error!\n";
    cout << "  N <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double f_noncentral_mean ( double a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F_NONCENTRAL_MEAN returns the mean of the F noncentral PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, a parameter of the PDF.
//
//    Input, int M, N, parameters of the PDF.
//    1 <= M,
//    1 <= N.
//    Note, however, that the mean is not defined unless 3 <= N.
//
//    Output, double F_NONCENTAL_MEAN, the mean of the PDF.
//
{
  double mean;

  if ( n < 3 )
  {
    cout << " \n";
    cout << "F_NONCENTRAL_MEAN - Fatal error!\n";
    cout << "  The mean is not defined for N < 3.\n";
    exit ( 1 );
  }

  mean = ( ( double ) ( m ) + a ) * ( double ) ( n )
    / ( ( double ) ( m ) * ( double ) ( n - 2 ) );

  return mean;
}
//****************************************************************************80

double f_noncentral_variance ( double a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F_NONCENTRAL_VARIANCE returns the variance of the F noncentral PDF.
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
//    Input, double A, a parameter of the PDF.
//
//    Input, int M, N, parameters of the PDF.
//    1 <= M,
//    1 <= N.
//    Note, however, that the variance is not defined unless 5 <= N.
//
//    Output, double F_NONCENTRAL_VARIANCE, the variance of the PDF.
//
{
  double mr;
  double nr;
  double variance;

  if ( n < 5 )
  {
    cout << " \n";
    cout << "F_NONCENTRAL_VARIANCE - Fatal error!\n";
    cout << "  The variance is not defined for N < 5.\n";
    exit ( 1 );
  }

  mr = ( double ) ( m );
  nr = ( double ) ( n );

  variance = ( ( mr + a ) * ( mr + a ) + 2.0 * ( mr + a ) * nr * nr )
    / ( ( nr - 2.0 ) * ( nr - 4.0 ) * mr * mr )
    - ( mr + a ) * ( mr + a ) * nr * nr
    / ( ( nr - 2.0 ) * ( nr - 2.0 ) * mr * mr );

  return variance;
}
//****************************************************************************80

double factorial_log ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    FACTORIAL_LOG returns the logarithm of N!.
//
//  Discussion:
//
//    N! = Product ( 1 <= I <= N ) I
//
//    N! = Gamma(N+1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the function.
//    0 <= N.
//
//    Output, double FACTORIAL_LOG, the logarithm of N!.
//
{
  int i;
  double value;

  if ( n < 0 )
  {
    cout << "\n";
    cout << "FACTORIAL_LOG - Fatal error!\n";
    cout << "  N < 0.\n";
    exit ( 1 );
  }

  value = 0.0;

  for ( i = 2; i <= n; i++ )
  {
    value = value + log ( ( double ) ( i ) );
  }

  return value;
}
//****************************************************************************80

double factorial_stirling ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    FACTORIAL_STIRLING computes Stirling's approximation to N!.
//
//  Discussion:
//
//    N! = Product ( 1 <= I <= N ) I
//
//    Stirling ( N ) = sqrt ( 2 * PI * N ) * ( N / E )**N * E**(1/(12*N) )
//
//    This routine returns the raw approximation for all nonnegative
//    values of N.  If N is less than 0, the value is returned as 0,
//    and if N is 0, the value of 1 is returned.  In all other cases,
//    Stirling's formula is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the function.
//
//    Output, double FACTORIAL_STIRLING, an approximation to N!.
//
{
  const double e_natural = 2.71828182845904523536028747135266249775724709369995;
  const double pi = 3.14159265358979323;
  double value;

  if ( n < 0 )
  {
    value = 0.0;
  }
  else if ( n == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = sqrt ( 2.0 * pi * ( double ) ( n ) )
      * pow ( ( double ) ( n ) / e_natural, n )
      * exp ( 1.0 / ( double ) ( 12 * n ) );
  }

  return value;
}
//****************************************************************************80

double fisher_pdf ( double x[3], double kappa, double mu[3] )

//****************************************************************************80
//
//  Purpose:
//
//    FISHER_PDF evaluates the Fisher PDF.
//
//  Discussion:
//
//    The formula for the PDF is:
//
//      PDF(KAPPA,MU;X) = C(KAPPA) * exp ( KAPPA * MU' * X )
//
//    where:
//
//      0 <= KAPPA is the concentration parameter,
//      MU is a point on the unit sphere, the mean direction,
//      X is any point on the unit sphere,
//      and C(KAPPA) is a normalization factor:
//
//      C(KAPPA) = sqrt ( KAPPA ) / ( ( 2 * pi )**(3/2) * I(0.5,KAPPA) )
//
//    where
//
//      I(nu,X) is the Bessel function of order NU and argument X.
//
//    For a fixed value of MU, the value of KAPPA determines the
//    tendency of sample points to tend to be near MU.  In particular,
//    KAPPA = 0 corresponds to a uniform distribution of points on the
//    sphere, but as KAPPA increases, the sample points will tend to
//    cluster more closely to MU.
//
//    The Fisher distribution for points on the unit sphere is
//    analogous to the normal distribution of points on a line,
//    and, more precisely, to the von Mises distribution on a circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kanti Mardia, Peter Jupp,
//    Directional Statistics,
//    Wiley, 2000,
//    LC: QA276.M335
//
//  Parameters:
//
//    Input, double X[3], the argument of the PDF.
//    X should have unit Euclidean norm, but this routine will
//    automatically work with a normalized version of X.
//
//    Input, double KAPPA, the concentration parameter.
//    0 <= KAPPA is required.
//
//    Input, double MU[3], the mean direction.
//    MU should have unit Euclidean norm, but this routine will
//    automatically work with a normalized version of MU.
//
//    Output, double FISHER_PDF, the value of the PDF.
//
{
# define NB 1

  double alpha;
  double arg;
  double b[NB];
  double cf;
  int ize;
  double mu_norm;
  int ncalc;
  double pdf;
  double pi = 3.141592653589793;
  double x_norm;

  if ( kappa < 0.0 )
  {
    cout << "\n";
    cout << "FISHER_PDF - Fatal error!\n";
    cout << "  KAPPA must be nonnegative.\n";
    cout << "  Input KAPPA = " << kappa << "\n";
    exit ( 1 );
  }

  if ( kappa == 0.0 )
  {
    pdf = 1.0 / ( 4.0 * pi );
    return pdf;
  }

  alpha = 0.5;
  ize = 1;

  ncalc = ribesl ( kappa, alpha, NB, ize, b );

  cf = sqrt ( kappa ) / ( sqrt ( pow ( 2.0 * pi, 3 ) ) * b[0] );

  mu_norm = r8vec_length ( 3, mu );

  if ( mu_norm == 0.0 )
  {
    pdf = cf;
    return pdf;
  }

  x_norm = r8vec_length ( 3, x );

  if ( x_norm == 0.0 )
  {
    pdf = cf;
    return pdf;
  }

  arg = kappa * r8vec_dot ( 3, x, mu ) / ( x_norm * mu_norm );

  pdf = cf * exp ( arg );

  return pdf;
# undef NB
}
//****************************************************************************80

double *fisher_sample ( double kappa, double mu[], int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    FISHER_SAMPLE samples the Fisher distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Nicholas Fisher, Toby Lewis, Brian Embleton,
//    Statistical Analysis of Spherical Data,
//    Cambridge, 2003,
//    ISBN13: 978-0521456999,
//    LC: QA276.F489.
//
//  Parameters:
//
//    Input, double KAPPA, the concentration parameter.
//
//    Input, double MU[3], the mean direction.
//    MU should have unit Euclidean norm, but this routine will
//    automatically work with a normalized version of MU.
//
//    Input, int N, the number of samples to choose.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double FISHER_SAMPLE[3*N], a sample of the Fisher distribution.
//
//  Local Parameters:
//
//    Local, double ALPHA, BETA, the colatitude (theta) and
//    longitude (phi) of the mean direction.
//
{
  double a[3*3];
  double alpha;
  double beta;
  int i;
  int j;
  int k;
  double lambda;
  double mu_norm;
  double *phi;
  double pi = 3.141592653589793;
  double rst[3];
  double *theta;
  double *xyz;

  mu_norm = r8vec_length ( 3, mu );

  if ( mu_norm == 0.0 )
  {
    cout << "\n";
    cout << "FISHER_SAMPLE - Fatal error!\n";
    cout << "  MU = 0.\n";
    exit ( 1 );
  }

  alpha = - acos ( mu[2] / mu_norm );
  beta = atan2 ( mu[1], mu[0] );

  lambda = exp ( - 2.0 * kappa );

  theta = r8vec_uniform_01 ( n, seed );

  for ( k = 0; k < n; k++ )
  {
    if ( kappa == 0.0 )
    {
      theta[k] = 2.0 * asin ( sqrt ( 1.0 - theta[k] ) );
    }
    else
    {
      theta[k] = 2.0 * asin ( sqrt (
        - log ( theta[k] * ( 1.0 - lambda ) + lambda )
        / ( 2.0 * kappa ) ) );
    }
  }

  phi = r8vec_uniform_01 ( n, seed );

  for ( k = 0; k < n; k++ )
  {
    phi[k] = 2.0 * pi * phi[k];
  }
//
//  Compute the rotation matrix.
//
  a[0+0*3] =   cos ( alpha ) * cos ( beta );
  a[1+0*3] =                 - sin ( beta );
  a[2+0*3] =   sin ( alpha ) * cos ( beta );

  a[0+1*3] =   cos ( alpha ) * sin ( beta );
  a[1+1*3] =                 + cos ( beta );
  a[2+1*3] =   sin ( alpha ) * sin ( beta );

  a[0+2*3] = - sin ( alpha );
  a[1+2*3] =   0.0;
  a[2+2*3] =   cos ( alpha );
//
//  Compute the unrotated points.
//
  xyz = new double[3*n];

  for ( k = 0; k < n; k++ )
  {
    rst[0] = sin ( theta[k] ) * cos ( phi[k] );
    rst[1] = sin ( theta[k] ) * sin ( phi[k] );
    rst[2] = cos ( theta[k] );
//
//  Rotate the points.
//
    for ( i = 0; i < 3; i++ )
    {
      xyz[i+k*3] = 0.0;
      for ( j = 0; j < 3; j++ )
      {
        xyz[i+k*3] = xyz[i+k*3] + a[i+j*3] * rst[j];
      }
    }
  }

  delete [] theta;
  delete [] phi;

  return xyz;
}
//****************************************************************************80

double fisk_cdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    FISK_CDF evaluates the Fisk CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double FISK_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 1.0 / ( 1.0 + pow ( ( b / ( x - a ) ), c ) );
  }

  return cdf;
}
//****************************************************************************80

double fisk_cdf_inv ( double cdf, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    FISK_CDF_INV inverts the Fisk CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double FISK_CDF_INV, the corresponding argument of the CDF.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "FISK_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf <= 0.0 )
  {
    x = a;
  }
  else if ( cdf < 1.0 )
  {
    x = a + b * pow ( cdf / ( 1.0 - cdf ), 1.0 / c );
  }
  else if ( 1.0 <= cdf )
  {
    x = r8_huge ( );
  }

  return x;
}
//****************************************************************************80

bool fisk_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    FISK_CHECK checks the parameters of the Fisk PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, bool FISK_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "FISK_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  if ( c <= 0.0 )
  {
    cout << " \n";
    cout << "FISK_CHECK - Fatal error!\n";
    cout << "  C <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double fisk_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    FISK_MEAN returns the mean of the Fisk PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double FISK_MEAN, the mean of the PDF.
//
{
  double mean;
  const double pi = 3.14159265358979323;

  if ( c <= 1.0 )
  {
    cout << " \n";
    cout << "FISK_MEAN - Fatal error!\n";
    cout << "  No mean defined for C <= 1.0\n";
    exit ( 1 );
  }

  mean = a + pi * ( b / c ) * r8_csc ( pi / c );

  return mean;
}
//****************************************************************************80

double fisk_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    FISK_PDF evaluates the Fisk PDF.
//
//  Discussion:
//
//    PDF(A,B,C;X) =
//      ( C / B ) * ( ( X - A ) / B )^( C - 1 ) /
//      ( 1 + ( ( X - A ) / B )**C )^2
//
//    The Fisk PDF is also known as the Log Logistic PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double FISK_PDF, the value of the PDF.
//
{
  double pdf;
  double y;

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else
  {
    y = ( x - a ) / b;

    pdf = ( c / b ) * pow ( y, ( c - 1.0 ) ) / pow ( ( 1.0 + pow ( y, c ) ), 2 );
  }

  return pdf;
}
//****************************************************************************80

double fisk_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    FISK_SAMPLE samples the Fisk PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double FISK_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = fisk_cdf_inv ( cdf, a, b, c );

  return x;
}
//****************************************************************************80

double fisk_variance ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    FISK_VARIANCE returns the variance of the Fisk PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double FISK_VARIANCE, the variance of the PDF.
//
{
  double g;
  const double pi = 3.14159265358979323;
  double variance;

  if ( c <= 2.0 )
  {
    cout << " \n";
    cout << "FISK_VARIANCE - Fatal error!\n";
    cout << "  No variance defined for C <= 2.0\n";
    exit ( 1 );
  }

  g = pi / c;

  variance = b * b * ( 2.0 * g * r8_csc ( 2.0 * g ) 
    - pow ( ( g * r8_csc ( g ) ), 2 ) );

  return variance;
}
//****************************************************************************80

double folded_normal_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_CDF evaluates the Folded Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//    0.0 <= X.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A,
//    0.0 < B.
//
//    Output, double FOLDED_NORMAL_CDF, the value of the CDF.
//
{
  double cdf;
  double cdf1;
  double cdf2;
  double x1;
  double x2;

  if ( x < 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    x1 = ( x - a ) / b;
    cdf1 = normal_01_cdf ( x1 );
    x2 = ( - x - a ) / b;
    cdf2 = normal_01_cdf ( x2 );
    cdf = cdf1 - cdf2;
  }

  return cdf;
}
//****************************************************************************80

double folded_normal_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_CDF_INV inverts the Folded Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A,
//    0.0 < B.
//
//    Output, double FOLDED_NORMAL_CDF_INV, the argument of the CDF.
//    0.0 <= X.
//
{
  double cdf1;
  double cdf2;
  double cdf3;
  int it;
  int it_max = 100;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;
  double xa;
  double xb;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "FOLDED_NORMAL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = 0.0;
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = r8_huge ( );
    return x;
  }
//
//  Find X1, for which the value of CDF will be too small.
//
  if ( 0.0 <= a )
  {
    x1 = normal_cdf_inv ( cdf, a, b );
  }
  else
  {
    x1 = normal_cdf_inv ( cdf, -a, b );
  }
  x1 = r8_max ( x1, 0.0 );
  cdf1 = folded_normal_cdf ( x1, a, b );
//
//  Find X2, for which the value of CDF will be too big.
//
  cdf2 = ( 1.0 - cdf ) / 2.0;

  xa = normal_cdf_inv ( cdf2, a, b );
  xb = normal_cdf_inv ( cdf2, -a, b );
  x2 = r8_max ( r8_abs ( xa ), r8_abs ( xb ) );
  cdf2 = folded_normal_cdf ( x2, a, b );
//
//  Now use bisection.
//
  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = folded_normal_cdf ( x3, a, b );

    if ( r8_abs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      break;
    }

    if ( it_max < it )
    {
      cout << " \n";
      cout << "FOLDED_NORMAL_CDF_INV - Fatal error!\n";
      cout << "  Iteration limit exceeded.\n";
      exit ( 1 );
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
      cdf2 = cdf3;
    }
  }

  return x;
}
//****************************************************************************80

bool folded_normal_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_CHECK checks the parameters of the Folded Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A,
//    0.0 < B.
//
//    Output, bool FOLDED_NORMAL_CHECK, is true if the parameters are legal.
//
{
  if ( a < 0.0 )
  {
    cout << " \n";
    cout << "FOLDED_NORMAL_CHECK - Fatal error!\n";
    cout << "  A < 0.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "FOLDED_NORMAL_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double folded_normal_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_MEAN returns the mean of the Folded Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A,
//    0.0 < B.
//
//    Output, double FOLDED_NORMAL_MEAN, the mean of the PDF.
//
{
  double a2;
  double cdf;
  double mean;
  const double pi = 3.14159265358979323;

  a2 = a / b;

  cdf = normal_01_cdf ( a2 );

  mean = b * sqrt ( 2.0 / pi ) * exp ( - 0.5 * a2 * a2 )
    - a * ( 1.0 - 2.0 * cdf );

  return mean;
}
//****************************************************************************80

double folded_normal_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_PDF evaluates the Folded Normal PDF.
//
//  Discussion:
//
//    PDF(A;X) = sqrt ( 2 / PI ) * ( 1 / B ) * COSH ( A * X / B^2 )
//      * EXP ( - 0.5 * ( X^2 + A^2 ) / B^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A,
//    0.0 < B.
//
//    Output, double FOLDED_NORMAL_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = sqrt ( 2.0 / pi ) * ( 1.0 / b ) * cosh ( a * x / b / b )
      * exp ( - 0.5 * ( x * x + a * a ) / b / b );
  }

  return pdf;
}
//****************************************************************************80

double folded_normal_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_SAMPLE samples the Folded Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A,
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double FOLDED_NORMAL_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = folded_normal_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double folded_normal_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_VARIANCE returns the variance of the Folded Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A,
//    0.0 < B.
//
//    Output, double FOLDED_NORMAL_VARIANCE, the variance of the PDF.
//
{
  double mean;
  double variance;

  mean = folded_normal_mean ( a, b );

  variance = a * a + b * b - mean * mean;

  return variance;
}
//****************************************************************************80

double frechet_cdf ( double x, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_CDF evaluates the Frechet CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the parameter.
//    It is required that 0.0 < ALPHA.
//
//    Input, double X, the argument of the CDF.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;

  if ( alpha <= 0.0 )
  {
    cout << "\n";
    cout << "FRECHET_CDF - Fatal error!\n";
    cout << "  ALPHA <= 0.0.\n";
    exit ( 1 );
  }

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = exp ( - 1.0 / pow ( x, alpha ) );
  }

  return cdf;
}
//****************************************************************************80

double frechet_cdf_inv ( double cdf, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_CDF_INV inverts the Frechet CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double ALPHA, the parameter.
//    It is required that 0.0 < ALPHA.
//
//    Output, double X, the corresponding argument of the CDF.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "FRECHET_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( alpha <= 0.0 )
  {
    cout << "\n";
    cout << "FRECHET_CDF_INV - Fatal error!\n";
    cout << "  ALPHA <= 0.0.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = 0.0;
  }
  else
  {
    x =  pow ( - 1.0 / log ( cdf ), 1.0 / alpha );
  }

  return x;
}
//****************************************************************************80

double frechet_mean ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_MEAN returns the mean of the Frechet PDF.
//
//  Discussion:
//
//    The distribution does not have a mean value unless 1 < ALPHA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the parameter.
//    It is required that 1.0 < ALPHA.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  if ( alpha <= 1.0 )
  {
    cout << "\n";
    cout << "FRECHET_MEAN - Fatal error!\n";
    cout << "  Mean does not exist if ALPHA <= 1.\n";
    exit ( 1 );
  }

  mean = r8_gamma ( ( alpha - 1.0 ) / alpha );

  return mean;
}
//****************************************************************************80

double frechet_pdf ( double x, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_PDF evaluates the Frechet PDF.
//
//  Discussion:
//
//    PDF(X) = ALPHA * exp ( -1 / X^ALPHA ) / X^(ALPHA+1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double ALPHA, the parameter.
//    It is required that 0.0 < ALPHA.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;

  if ( alpha <= 0.0 )
  {
    cout << "\n";
    cout << "FRECHET_PDF - Fatal error!\n";
    cout << "  ALPHA <= 0.0.\n";
    exit ( 1 );
  }

  pdf = alpha * exp ( - 1.0 / pow ( x, alpha ) ) / pow ( x, alpha + 1.0 );

  return pdf;
}
//****************************************************************************80

double frechet_sample ( double alpha, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_SAMPLE samples the Frechet PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the parameter.
//    It is required that 0.0 < ALPHA.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double cdf;
  double x;

  if ( alpha <= 0.0 )
  {
    cout << "\n";
    cout << "FRECHET_SAMPLE - Fatal error!\n";
    cout << "  ALPHA <= 0.0.\n";
    exit ( 1 );
  }

  cdf = r8_uniform_01 ( seed );

  x = frechet_cdf_inv ( cdf, alpha );

  return x;
}
//****************************************************************************80

double frechet_variance ( double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_VARIANCE returns the variance of the Frechet PDF.
//
//  Discussion:
//
//    The PDF does not have a variance unless 2 < ALPHA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ALPHA, the parameter.
//    It is required that 2.0 < ALPHA.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double mean;
  double variance;

  if ( alpha <= 2.0 )
  {
    cout << "\n";
    cout << "FRECHET_VARIANCE - Fatal error!\n";
    cout << "  Variance does not exist if ALPHA <= 2.\n";
    exit ( 1 );
  }

  mean = r8_gamma ( ( alpha - 1.0 ) / alpha );

  variance = r8_gamma ( ( alpha - 2.0 ) / alpha ) - mean * mean;

  return variance;
}
//****************************************************************************80

double gamma_cdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_CDF evaluates the Gamma CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double GAMMA_CDF, the value of the CDF.
//
{
  double cdf;
  double p2;
  double x2;

  x2 = ( x - a ) / b;
  p2 = c;

  cdf = gamma_inc ( p2, x2 );

  return cdf;
}
//****************************************************************************80

void gamma_cdf_values ( int *n_data, double *mu, double *sigma, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_CDF_VALUES returns some values of the Gamma CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = GammaDistribution [ mu, sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
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
//    Output, double *MU, the mean of the distribution.
//
//    Output, double *SIGMA, the variance of the distribution.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double fx_vec[N_MAX] = {
     0.8646647167633873,
     0.9816843611112658,
     0.9975212478233336,
     0.9996645373720975,
     0.6321205588285577,
     0.4865828809674080,
     0.3934693402873666,
     0.3296799539643607,
     0.4421745996289254,
     0.1911531694619419,
     0.6564245437845009E-01,
     0.1857593622214067E-01 };

  double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  double sigma_vec[N_MAX] = {
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *mu = 0.0;
    *sigma = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *mu = mu_vec[*n_data-1];
    *sigma = sigma_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool gamma_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_CHECK checks the parameters of the Gamma PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, bool GAMMA_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "GAMMA_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    cout << "  B = " << b << "\n";
    return false;
  }

  if ( c <= 0.0 )
  {
    cout << " \n";
    cout << "GAMMA_CHECK - Fatal error!\n";
    cout << "  C <= 0.\n";
    cout << "  C = " << c << "\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double gamma_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_MEAN returns the mean of the Gamma PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double GAMMA_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a + b * c;

  return mean;
}
//****************************************************************************80

double gamma_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_PDF evaluates the Gamma PDF.
//
//  Discussion:
//
//    PDF(A,B,C;X) = EXP ( - ( X - A ) / B ) * ( ( X - A ) / B )^(C-1)
//      / ( B * GAMMA ( C ) )
//
//    GAMMA_PDF(A,B,C;X), where C is an integer, is the Erlang PDF.
//    GAMMA_PDF(A,B,1;X) is the Exponential PDF.
//    GAMMA_PDF(0,2,C/2;X) is the Chi Squared PDF with C degrees of freedom.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, C, the parameters of the PDF.
//    A controls the location of the peak;  A is often chosen to be 0.0.
//    B is the "scale" parameter; 0.0 < B, and is often 1.0.
//    C is the "shape" parameter; 0.0 < C, and is often 1.0.
//
//    Output, double GAMMA_PDF, the value of the PDF.
//
{
  double pdf;
  double y;

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else
  {
    y = ( x - a ) / b;

    pdf = pow ( y, c - 1.0 ) / ( b * r8_gamma ( c ) * exp ( y ) );
  }

  return pdf;
}
//****************************************************************************80

double gamma_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_SAMPLE samples the Gamma PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    Original FORTRAN77 version by Ahrens and U Dieter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    J H Ahrens and U Dieter,
//    Generating Gamma Variates by a Modified Rejection Technique,
//    Communications of the ACM,
//    Volume 25, Number 1, January 1982, pages 47 - 54.
//
//    J H Ahrens and U Dieter,
//    Computer Methods for Sampling from Gamma, Beta, Poisson and
//      Binomial Distributions.
//    Computing, Volume 12, 1974, pages 223 - 246.
//
//    J H Ahrens, K D Kohrt, and U Dieter,
//    Algorithm 599,
//    ACM Transactions on Mathematical Software,
//    Volume 9, Number 2, June 1983, pages 255-257.
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double a1 =   0.3333333;
  double a2 = - 0.2500030;
  double a3 =   0.2000062;
  double a4 = - 0.1662921;
  double a5 =   0.1423657;
  double a6 = - 0.1367177;
  double a7 =   0.1233795;
  double bcoef;
  double co;
  double d;
  double e;
  double e1 = 1.0;
  double e2 = 0.4999897;
  double e3 = 0.1668290;
  double e4 = 0.0407753;
  double e5 = 0.0102930;
  double euler = 2.71828182845904;
  double p;
  double q;
  double q0;
  double q1 =  0.04166669;
  double q2 =  0.02083148;
  double q3 =  0.00801191;
  double q4 =  0.00144121;
  double q5 = -0.00007388;
  double q6 =  0.00024511;
  double q7 =  0.00024240;
  double r;
  double s;
  double si;
  double s2;
  double t;
  double u;
  double v;
  double w;
  double x;
//
//  Allow C = 0.
//
  if ( c == 0.0 )
  {
    x = a;
    return x;
  }
//
//  C < 1.
//
  if ( c < 1.0 )
  {
    for ( ; ; )
    {
      u = r8_uniform_01 ( seed );
      t = 1.0 + c / euler;
      p = u * t;

      s = exponential_01_sample ( seed );

      if ( p < 1.0 )
      {
        x = exp ( log ( p ) / c );
        if ( x <= s )
        {
          break;
        }
      }
      else
      {
        x = - log ( ( t - p ) / c );
        if ( ( 1.0 - c ) * log ( x ) <= s )
        {
          break;
        }
      }
    }

    x = a + b * x;
    return x;
  }
//
//  1 <= C.
//
  else
  {
    s2 = c - 0.5;
    s = sqrt ( c - 0.5 );
    d = sqrt ( 32.0 ) - 12.0 * sqrt ( c - 0.5 );

    t = normal_01_sample ( seed );
    x = pow ( ( sqrt ( c - 0.5 ) + 0.5 * t ), 2 );

    if ( 0.0 <= t )
    {
      x = a + b * x;
      return x;
    }

    u = r8_uniform_01 ( seed );

    if ( d * u <= t * t * t )
    {
      x = a + b * x;
      return x;
    }

    r = 1.0 / c;

    q0 = ( ( ( ( ( (
           q7   * r
         + q6 ) * r
         + q5 ) * r
         + q4 ) * r
         + q3 ) * r
         + q2 ) * r
         + q1 ) * r;

    if ( c <= 3.686 )
    {
      bcoef = 0.463 + s - 0.178 * s2;
      si = 1.235;
      co = 0.195 / s - 0.079 + 0.016 * s;
    }
    else if ( c <= 13.022 )
    {
      bcoef = 1.654 + 0.0076 * s2;
      si = 1.68 / s + 0.275;
      co = 0.062 / s + 0.024;
    }
    else
    {
      bcoef = 1.77;
      si = 0.75;
      co = 0.1515 / s;
    }

    if ( 0.0 < sqrt ( c - 0.5 ) + 0.5 * t )
    {
      v = 0.5 * t / s;

      if ( 0.25 < r8_abs ( v ) )
      {
        q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
      }
      else
      {
        q = q0 + 0.5 * t * t * ( ( ( ( ( (
               a7   * v
             + a6 ) * v
             + a5 ) * v
             + a4 ) * v
             + a3 ) * v
             + a2 ) * v
             + a1 ) * v;
      }

      if ( log ( 1.0 - u ) <= q )
      {
        x = a + b * x;
        return x;
      }
    }

    for ( ; ; )
    {
      e = exponential_01_sample ( seed );

      u = r8_uniform_01 ( seed );

      u = 2.0 * u - 1.0;
      t = bcoef + r8_abs ( si * e ) * r8_sign ( u );

      if ( -0.7187449 <= t )
      {
        v = 0.5 * t / s;

        if ( 0.25 < r8_abs ( v ) )
        {
          q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log ( 1.0 + v );
        }
        else
        {
          q = q0 + 0.5 * t * t * ( ( ( ( ( (
               a7   * v
             + a6 ) * v
             + a5 ) * v
             + a4 ) * v
             + a3 ) * v
             + a2 ) * v
             + a1 ) * v;
        }

        if ( 0.0 < q )
        {
          if ( 0.5 < q )
          {
            w = exp ( q ) - 1.0;
          }
          else
          {
            w = ( ( ( (
                    e5   * q
                  + e4 ) * q
                  + e3 ) * q
                  + e2 ) * q
                  + e1 ) * q;
          }

          if ( co * r8_abs ( u ) <= w * exp ( e - 0.5 * t * t ) )
          {
            x = a + b * pow ( s + 0.5 * t, 2 );
            return x;
          }
        }
      }
    }
  }
}
//****************************************************************************80

double gamma_variance ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_VARIANCE returns the variance of the Gamma PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = b * b * c;

  return variance;
}
//****************************************************************************80

double gamma_inc ( double p, double x )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC computes the incomplete Gamma function.
//
//  Discussion:
//
//    GAMMA_INC(P,       0) = 0,
//    GAMMA_INC(P,Infinity) = 1.
//
//    GAMMA_INC(P,X) = Integral ( 0 <= T <= X ) T**(P-1) EXP(-T) DT / GAMMA(P).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    Original FORTRAN77 version by B L Shea.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    B L Shea,
//    Chi-squared and Incomplete Gamma Integral,
//    Algorithm AS239,
//    Applied Statistics,
//    Volume 37, Number 3, 1988, pages 466-473.
//
//  Parameters:
//
//    Input, double P, the exponent parameter.
//    0.0 < P.
//
//    Input, double X, the integral limit parameter.
//    If X is less than or equal to 0, GAMMA_INC is returned as 0.
//
//    Output, double GAMMA_INC, the value of the function.
//
{
  double a;
  double arg;
  double b;
  double c;
  double exp_arg_min = -88.0;
  double overflow = 1.0E+37;
  double plimit = 1000.0;
  double pn1;
  double pn2;
  double pn3;
  double pn4;
  double pn5;
  double pn6;
  double rn;
  double value;
  double tol = 1.0E-07;
  double xbig = 1.0E+08;

  value = 0.0;

  if ( p <= 0.0 )
  {
    cout << " \n";
    cout << "GAMMA_INC - Fatal error!\n";
    cout << "  Parameter P <= 0.\n";
    exit ( 1 );
  }

  if ( x <= 0.0 )
  {
    value = 0.0;
    return value;
  }
//
//  Use a normal approximation if PLIMIT < P.
//
  if ( plimit < p )
  {
    pn1 = 3.0 * sqrt ( p ) * ( pow ( x / p, 1.0 / 3.0 )
      + 1.0 / ( 9.0 * p ) - 1.0 );
    value = normal_01_cdf ( pn1 );
    return value;
  }
//
//  Is X extremely large compared to P?
//
  if ( xbig < x )
  {
    value = 1.0;
    return value;
  }
//
//  Use Pearson's series expansion.
//  (P is not large enough to force overflow in the log of Gamma.
//
  if ( x <= 1.0 || x < p )
  {
    arg = p * log ( x ) - x - gamma_log ( p + 1.0 );
    c = 1.0;
    value = 1.0;
    a = p;

    for ( ; ; )
    {
      a = a + 1.0;
      c = c * x / a;
      value = value + c;

      if ( c <= tol )
      {
        break;
      }
    }

    arg = arg + log ( value );

    if ( exp_arg_min <= arg )
    {
      value = exp ( arg );
    }
    else
    {
      value = 0.0;
    }
  }
//
//  Use a continued fraction expansion.
//
  else
  {
    arg = p * log ( x ) - x - gamma_log ( p );
    a = 1.0 - p;
    b = a + x + 1.0;
    c = 0.0;
    pn1 = 1.0;
    pn2 = x;
    pn3 = x + 1.0;
    pn4 = x * b;
    value = pn3 / pn4;

    for ( ; ; )
    {
      a = a + 1.0;
      b = b + 2.0;
      c = c + 1.0;
      pn5 = b * pn3 - a * c * pn1;
      pn6 = b * pn4 - a * c * pn2;

      if ( 0.0 < r8_abs ( pn6 ) )
      {
        rn = pn5 / pn6;

        if ( r8_abs ( value - rn ) <= r8_min ( tol, tol * rn ) )
        {
          arg = arg + log ( value );

          if ( exp_arg_min <= arg )
          {
            value = 1.0 - exp ( arg );
          }
          else
          {
            value = 1.0;
          }

          return value;
        }
        value = rn;
      }
      pn1 = pn3;
      pn2 = pn4;
      pn3 = pn5;
      pn4 = pn6;
//
//  Rescale terms in continued fraction if terms are large.
//
      if ( overflow <= r8_abs ( pn5 ) )
      {
        pn1 = pn1 / overflow;
        pn2 = pn2 / overflow;
        pn3 = pn3 / overflow;
        pn4 = pn4 / overflow;
      }
    }
  }

  return value;
}
//****************************************************************************80

void gamma_inc_values ( int *n_data, double *a, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
//
//  Discussion:
//
//    The (normalized) incomplete Gamma function P(A,X) is defined as:
//
//      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
//
//    With this definition, for all A and X,
//
//      0 <= PN(A,X) <= 1
//
//    and
//
//      PN(A,INFINITY) = 1.0
//
//    In Mathematica, the function can be evaluated by:
//
//      1 - GammaRegularized[A,X]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 November 2004
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
//    Output, double *A, the parameter of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 20

  double a_vec[N_MAX] = {
     0.10,
     0.10,
     0.10,
     0.50,
     0.50,
     0.50,
     0.10E+01,
     0.10E+01,
     0.10E+01,
     0.11E+01,
     0.11E+01,
     0.11E+01,
     0.20E+01,
     0.20E+01,
     0.20E+01,
     0.60E+01,
     0.60E+01,
     0.11E+02,
     0.26E+02,
     0.41E+02  };

  double fx_vec[N_MAX] = {
     0.7382350532339351,
     0.9083579897300343,
     0.9886559833621947,
     0.3014646416966613,
     0.7793286380801532,
     0.9918490284064973,
     0.9516258196404043E-01,
     0.6321205588285577,
     0.9932620530009145,
     0.7205974576054322E-01,
     0.5891809618706485,
     0.9915368159845525,
     0.1018582711118352E-01,
     0.4421745996289254,
     0.9927049442755639,
     0.4202103819530612E-01,
     0.9796589705830716,
     0.9226039842296429,
     0.4470785799755852,
     0.7444549220718699 };

  double x_vec[N_MAX] = {
     0.30E-01,
     0.30,
     0.15E+01,
     0.75E-01,
     0.75,
     0.35E+01,
     0.10,
     0.10E+01,
     0.50E+01,
     0.10,
     0.10E+01,
     0.50E+01,
     0.15,
     0.15E+01,
     0.70E+01,
     0.25E+01,
     0.12E+02,
     0.16E+02,
     0.25E+02,
     0.45E+02 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double gamma_log ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
//
//  Discussion:
//
//    The program uses rational functions that theoretically approximate
//    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
//    approximation for 12 < X is from Hart et al, while approximations
//    for X < 12.0 are similar to those in Cody and Hillstrom, but are
//    unpublished.
//
//    The accuracy achieved depends on the arithmetic system, the compiler,
//    intrinsic functions, and proper selection of the machine-dependent
//    constants.
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
//    William Cody, Kenneth Hillstrom,
//    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
//    Mathematics of Computation,
//    Volume 21, 1967, pages 198-203.
//
//    Kenneth Hillstrom,
//    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
//    May 1969.
//
//    Hart, Ward Cheney, Charles Lawson, Maehly,
//    Charles Mesztenyi, John Rice, Thacher, Witzgall,
//    Computer Approximations,
//    Wiley and sons, New York, 1968.
//
//  Parameters:
//
//    Input, double X, the argument of the Gamma function.  X must be positive.
//
//    Output, double GAMMA_LOG, the logarithm of the Gamma function of X.
//    If X <= 0.0, or if overflow would occur, the program returns the
//    value XINF, the largest representable floating point number.
//
//  Machine-dependent constants:
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

double gamma_log_int ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LOG_INT computes the logarithm of Gamma of an integer N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the logarithm of the Gamma function.
//    0 < N.
//
//    Output, double GAMMA_LOG_INT, the logarithm of
//    the Gamma function of N.
//
{
  double value;

  if ( n <= 0 )
  {
    cout << " \n";
    cout << "GAMMA_LOG_INT - Fatal error!\n";
    cout << "  Illegal input value of N = " << n << "\n";
    cout << "  But N must be strictly positive.\n";
    exit ( 1 );
  }

  value = gamma_log ( ( double ) ( n ) );

  return value;
}
//****************************************************************************80

double genlogistic_cdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_CDF evaluates the Generalized Logistic CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;
  double y;

  y = ( x - a ) / b;

  cdf = 1.0 / pow ( ( 1.0 + exp ( - y ) ), c );

  return cdf;
}
//****************************************************************************80

double genlogistic_cdf_inv ( double cdf, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_CDF_INV inverts the Generalized Logistic CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double X, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "GENLOGISTIC_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = - r8_huge ( );
  }
  else if ( cdf < 1.0 )
  {
    x = a - b * log ( pow ( cdf, - 1.0 / c ) - 1.0 );
  }
  else if ( 1.0 == cdf )
  {
    x = r8_huge ( );
  }

  return x;
}
//****************************************************************************80

bool genlogistic_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_CHECK checks the parameters of the Generalized Logistic CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, bool GENLOGISTIC_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "GENLOGISTIC_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  if ( c <= 0.0 )
  {
    cout << " \n";
    cout << "GENLOGISTIC_CHECK - Fatal error!\n";
    cout << "  C <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double genlogistic_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_MEAN returns the mean of the Generalized Logistic PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a + b * ( euler_constant ( ) + digamma ( c ) );

  return mean;
}
//****************************************************************************80

double genlogistic_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_PDF evaluates the Generalized Logistic PDF.
//
//  Discussion:
//
//    PDF(A,B,C;X) = ( C / B ) * EXP ( ( A - X ) / B ) /
//      ( ( 1 + EXP ( ( A - X ) / B ) )**(C+1) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  double y;

  y = ( x - a ) / b;

  pdf = ( c / b ) * exp ( - y ) / pow ( 1.0 + exp ( - y ), c + 1.0 );

  return pdf;
}
//****************************************************************************80

double genlogistic_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_SAMPLE samples the Generalized Logistic PDF.
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
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double GENLOGISTIC_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = genlogistic_cdf_inv ( cdf, a, b, c );

  return x;
}
//****************************************************************************80

double genlogistic_variance ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_VARIANCE returns the variance of the Generalized Logistic PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double variance;

  variance = b * b * ( pi * pi / 6.0 + trigamma ( c ) );

  return variance;
}
//****************************************************************************80

double geometric_cdf ( int x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_CDF evaluates the Geometric CDF.
//
//  Discussion:
//
//    CDF(X,P) is the probability that there will be at least one
//    successful trial in the first X Bernoulli trials, given that
//    the probability of success in a single trial is P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the maximum number of trials.
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Output, double GEOMETRIC_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= 0 )
  {
    cdf = 0.0;
  }
  else if ( a == 0.0 )
  {
    cdf = 0.0;
  }
  else if ( a == 1.0 )
  {
    cdf = 1.0;
  }
  else
  {
    cdf = 1.0 - pow ( ( 1.0 - a ), x );
  }

  return cdf;
}
//****************************************************************************80

int geometric_cdf_inv ( double cdf, double a )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_CDF_INV inverts the Geometric CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Output, int GEOMETRIC_CDF_INV, the corresponding value of X.
//
{
  int x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "GEOMETRIC_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( a == 1.0 )
  {
    x = 1;
  }
  else if ( a == 0.0 )
  {
    x = i4_huge ( );
  }
  else
  {
    x = 1 + int ( log ( 1.0 - cdf ) / log ( 1.0 - a ) );
  }

  return x;
}
//****************************************************************************80

void geometric_cdf_values ( int *n_data, int *x, double *p, double *cdf )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_CDF_VALUES returns values of the geometric CDF.
//
//  Discussion:
//
//    The geometric or Pascal probability density function gives the
//    probability that the first success will happen on the X-th Bernoulli
//    trial, given that the probability of a success on a single trial is P.
//
//    The value of CDF ( X, P ) is the probability that the first success
//    will happen on or before the X-th trial.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = GeometricDistribution [ p ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2004
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
//    Daniel Zwillinger and Stephen Kokoska,
//    CRC Standard Probability and Statistics Tables and Formulae,
//    Chapman and Hall / CRC Press, 2000.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *X, the number of trials.
//
//    Output, double *P, the probability of success
//    on one trial.
//
//    Output, double *CDF, the cumulative density function.
//
{
# define N_MAX 14

  double cdf_vec[N_MAX] = {
     0.1900000000000000,
     0.2710000000000000,
     0.3439000000000000,
     0.6861894039100000,
     0.3600000000000000,
     0.4880000000000000,
     0.5904000000000000,
     0.9141006540800000,
     0.7599000000000000,
     0.8704000000000000,
     0.9375000000000000,
     0.9843750000000000,
     0.9995117187500000,
     0.9999000000000000 };

  double p_vec[N_MAX] = {
     0.1,
     0.1,
     0.1,
     0.1,
     0.2,
     0.2,
     0.2,
     0.2,
     0.3,
     0.4,
     0.5,
     0.5,
     0.5,
     0.9 };

  int x_vec[N_MAX] = {
    1,  2,  3, 10, 1,
    2,  3, 10,  3, 3,
    3,  5, 10,  3 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0;
    *p = 0.0;
    *cdf = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *p = p_vec[*n_data-1];
    *cdf = cdf_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool geometric_check ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_CHECK checks the parameter of the Geometric CDF.
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
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Output, bool GEOMETRIC_CHECK, is true if the parameters are legal.
//
{
  if ( a < 0.0 || 1.0 < a )
  {
    cout << " \n";
    cout << "GEOMETRIC_CHECK - Fatal error!\n";
    cout << "  A < 0 or 1 < A.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double geometric_mean ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_MEAN returns the mean of the Geometric PDF.
//
//  Discussion:
//
//    MEAN is the expected value of the number of trials required
//    to obtain a single success.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = 1.0 / a;

  return mean;
}
//****************************************************************************80

double geometric_pdf ( int x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_PDF evaluates the Geometric PDF.
//
//  Discussion:
//
//    PDF(A;X) = A * ( 1 - A )**(X-1)
//
//    PDF(A;X) is the probability that exactly X Bernoulli trials, each
//    with probability of success A, will be required to achieve
//    a single success.
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
//    Input, int X, the number of trials.
//    0 < X
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
//
//  Special cases.
//
  if ( x < 1 )
  {
    pdf = 0.0;
  }
  else if ( a == 0.0 )
  {
    pdf = 0.0;
  }
  else if ( a == 1.0 )
  {
    if ( x == 1 )
    {
      pdf = 1.0;
    }
    else
    {
      pdf = 0.0;
    }
  }
  else
  {
    pdf = a * pow ( ( 1.0 - a ), ( x - 1 ) );

  }

  return pdf;
}
//****************************************************************************80

int geometric_sample ( double a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_SAMPLE samples the Geometric PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int GEOMETRIC_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  int x;

  cdf = r8_uniform_01 ( seed );

  x = geometric_cdf_inv ( cdf, a );

  return x;
}
//****************************************************************************80

double geometric_variance ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_VARIANCE returns the variance of the Geometric PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of success on one trial.
//    0.0 <= A <= 1.0.
//
//    Output, double GEOMETRIC_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( 1.0 - a ) / ( a * a );

  return variance;
}
//****************************************************************************80

int get_seed ( )

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
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GET_SEED, a random seed value.
//
{
# define I_MAX 2147483647
  time_t clock;
  int ihour;
  int imin;
  int isec;
  int seed;
  struct tm *lt;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  clock = time ( &tloc );
  lt = localtime ( &clock );
//
//  Hours is 1, 2, ..., 12.
//
  ihour = lt->tm_hour;

  if ( 12 < ihour )
  {
    ihour = ihour - 12;
  }
//
//  Move Hours to 0, 1, ..., 11
//
  ihour = ihour - 1;

  imin = lt->tm_min;

  isec = lt->tm_sec;

  seed = isec + 60 * ( imin + 60 * ihour );
//
//  We want values in [1,43200], not [0,43199].
//
  seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,IMAX].
//
  seed = ( int )
    ( ( ( double ) seed )
    * ( ( double ) I_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;
#undef I_MAX
}
//****************************************************************************80

double gompertz_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    GOMPERTZ_CDF evaluates the Gompertz CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johnson, Kotz, and Balakrishnan,
//    Continuous Univariate Distributions, Volume 2, second edition,
//    Wiley, 1994, pages 25-26.
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    1 < A, 0 < B.
//
//    Output, double GOMPERTZ_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 1.0 - exp ( - b * ( pow ( a, x ) - 1.0 ) / log ( a ) );
  }

  return cdf;
}
//****************************************************************************80

double gompertz_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    GOMPERTZ_CDF_INV inverts the Gompertz CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johnson, Kotz, and Balakrishnan,
//    Continuous Univariate Distributions, Volume 2, second edition,
//    Wiley, 1994, pages 25-26.
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    1 < A, 0 < B.
//
//    Output, double GOMPERTZ_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "GOMPERTZ_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf < 1.0 )
  {
    x = log ( 1.0 - log ( 1.0 - cdf ) * log ( a ) / b  ) / log ( a );
  }
  else
  {
    x = r8_huge ( );
  }

  return x;
}
//****************************************************************************80

bool gompertz_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    GOMPERTZ_CHECK checks the parameters of the Gompertz PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johnson, Kotz, and Balakrishnan,
//    Continuous Univariate Distributions, Volume 2, second edition,
//    Wiley, 1994, pages 25-26.
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    1 < A, 0 < B.
//
//    Output, bool GOMPERTZ_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 1.0 )
  {
    cout << " \n";
    cout << "GOMPERTZ_CHECK - Fatal error!\n";
    cout << "  A <= 1.0!\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "GOMPERTZ_CHECK - Fatal error!\n";
    cout << "  B <= 0.0!\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double gompertz_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    GOMPERTZ_PDF evaluates the Gompertz PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = B * A**X / exp ( B * ( A**X - 1 ) / log ( A ) )
//
//    for
//
//      0.0 <= X
//      1.0 <  A
//      0.0 <  B
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johnson, Kotz, and Balakrishnan,
//    Continuous Univariate Distributions, Volume 2, second edition,
//    Wiley, 1994, pages 25-26.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    1 < A, 0 < B.
//
//    Output, double GOMPERTZ_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else if ( 1.0 < a )
  {
    pdf = exp ( log ( b ) + x * log ( a )
      - ( b / log ( a ) ) * ( pow ( a, x ) - 1.0 ) );
  }

  return pdf;
}
//****************************************************************************80

double gompertz_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    GOMPERTZ_SAMPLE samples the Gompertz PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    1 < A, 0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double GOMPERTZ_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = gompertz_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double gumbel_cdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_CDF evaluates the Gumbel CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Output, double GUMBEL_CDF, the value of the CDF.
//
{
  double cdf;

  cdf = exp ( - exp ( - x ) );

  return cdf;
}
//****************************************************************************80

double gumbel_cdf_inv ( double cdf )

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_CDF_INV inverts the Gumbel CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Output, double GUMBEL_CDF_INV, the corresponding argument of the CDF.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "GUMBEL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x =  - log ( - log ( cdf ) );

  return x;
}
//****************************************************************************80

double gumbel_mean ( )

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_MEAN returns the mean of the Gumbel PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double GUMBEL_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = euler_constant ( );

  return mean;
}
//****************************************************************************80

double gumbel_pdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_PDF evaluates the Gumbel PDF.
//
//  Discussion:
//
//    PDF(X) = EXP ( - X - EXP ( - X  ) ).
//
//    GUMBEL_PDF(X) = EXTREME_PDF(0,1;X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein, editor,
//    CRC Concise Encylopedia of Mathematics,
//    CRC Press, 1998.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Output, double GUMBEL_PDF, the value of the PDF.
//
{
  double pdf;

  pdf = exp ( - x - exp ( - x ) );

  return pdf;
}
//****************************************************************************80

double gumbel_sample ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_SAMPLE samples the Gumbel PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double GUMBEL_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = gumbel_cdf_inv ( cdf );

  return x;
}
//****************************************************************************80

double gumbel_variance ( )

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_VARIANCE returns the variance of the Gumbel PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double GUMBEL_VARIANCE, the variance of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double variance;

  variance = pi * pi / 6.0;

  return variance;
}
//****************************************************************************80

double half_normal_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_CDF evaluates the Half Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double HALF_NORMAL_CDF, the value of the CDF.
//
{
  double cdf;
  double cdf2;

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else
  {
    cdf2 = normal_cdf ( x, a, b );
    cdf = 2.0 * cdf2 - 1.0;
  }

  return cdf;
}
//****************************************************************************80

double half_normal_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_CDF_INV inverts the Half Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double HALF_NORMAL_CDF_INV, the corresponding argument.
//
{
  double cdf2;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "HALF_NORMAL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  cdf2 = 0.5 * ( cdf + 1.0 );

  x = normal_cdf_inv ( cdf2, a, b );

  return x;
}
//****************************************************************************80

bool half_normal_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_CHECK checks the parameters of the Half Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, bool HALF_NORMAL_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "HALF_NORMAL_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double half_normal_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_MEAN returns the mean of the Half Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double HALF_NORMAL_MEAN, the mean of the PDF.
//
{
  double mean;
  const double pi = 3.14159265358979323;

  mean = a + b * sqrt ( 2.0 / pi );

  return mean;
}
//****************************************************************************80

double half_normal_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_PDF evaluates the Half Normal PDF.
//
//  Discussion:
//
//    PDF(A,B;X) =
//      sqrt ( 2 / PI ) * ( 1 / B ) * exp ( - 0.5 * ( ( X - A ) / B )^2 )
//
//    for A <= X
//
//    The Half Normal PDF is a special case of both the Chi PDF and the
//    Folded Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double HALF_NORMAL_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;
  double y;

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else
  {
    y = ( x - a ) / b;

    pdf = sqrt ( 2.0 / pi ) * ( 1.0 / b ) * exp ( - 0.5 * y * y );

  }

  return pdf;
}
//****************************************************************************80

double half_normal_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_SAMPLE samples the Half Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double HALF_NORMAL_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = half_normal_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double half_normal_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_VARIANCE returns the variance of the Half Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double HALF_NORMAL_VARIANCE, the variance of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double variance;

  variance = b * b * ( 1.0 - 2.0 / pi );

  return variance;
}
//****************************************************************************80

double hypergeometric_cdf ( int x, int n, int m, int l )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CDF evaluates the Hypergeometric CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the argument of the CDF.
//
//    Input, int N, the number of balls selected.
//    0 <= N <= L.
//
//    Input, int M, the number of white balls in the population.
//    0 <= M <= L.
//
//    Input, int L, the number of balls to select from.
//    0 <= L.
//
//    Output, double HYPERGEOMETRIC_CDF, the value of the CDF.
//
{
  double cdf;
  double c1_log;
  double c2_log;
  double pdf;
  int x2;

  c1_log = binomial_coef_log ( l - m, n );
  c2_log = binomial_coef_log ( l, n );

  pdf = exp ( c1_log - c2_log );
  cdf = pdf;

  for ( x2 = 0; x2 <= x - 1; x2++ )
  {
    pdf = pdf * ( double ) ( ( m - x2 ) * ( n - x2 ) )
      / ( double ) ( ( x2 + 1 ) * ( l - m - n + x2 + 1 ) );

    cdf = cdf + pdf;
  }

  return cdf;
}
//****************************************************************************80

void hypergeometric_cdf_values ( int *n_data, int *sam, int *suc, int *pop,
  int *n, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = HypergeometricDistribution [ sam, suc, pop ]
//      CDF [ dist, n ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
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
//    Output, int *SAM, int *SUC, int *POP, the sample size,
//    success size, and population parameters of the function.
//
//    Output, int *N, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 16

  double fx_vec[N_MAX] = {
     0.6001858177500578E-01,
     0.2615284665839845,
     0.6695237889132748,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.5332595856827856,
     0.1819495964117640,
     0.4448047017527730E-01,
     0.9999991751316731,
     0.9926860896560750,
     0.8410799901444538,
     0.3459800113391901,
     0.0000000000000000,
     0.2088888139634505E-02,
     0.3876752992448843,
     0.9135215248834896 };

  int n_vec[N_MAX] = {
     7,  8,  9, 10,
     6,  6,  6,  6,
     6,  6,  6,  6,
     0,  0,  0,  0 };

  int pop_vec[N_MAX] = {
    100, 100, 100, 100,
    100, 100, 100, 100,
    100, 100, 100, 100,
    90,  200, 1000, 10000 };

  int sam_vec[N_MAX] = {
    10, 10, 10, 10,
     6,  7,  8,  9,
    10, 10, 10, 10,
    10, 10, 10, 10 };

  int suc_vec[N_MAX] = {
    90, 90, 90, 90,
    90, 90, 90, 90,
    10, 30, 50, 70,
    90, 90, 90, 90 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *sam = 0;
    *suc = 0;
    *pop = 0;
    *n = 0;
    *fx = 0.0;
  }
  else
  {
    *sam = sam_vec[*n_data-1];
    *suc = suc_vec[*n_data-1];
    *pop = pop_vec[*n_data-1];
    *n = n_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool hypergeometric_check ( int n, int m, int l )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CHECK checks the parameters of the Hypergeometric CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of balls selected.
//    0 <= N <= L.
//
//    Input, int M, the number of white balls in the population.
//    0 <= M <= L.
//
//    Input, int L, the number of balls to select from.
//    0 <= L.
//
//    Output, bool HYPERGEOMETRIC_CHECK, is true if the parameters are legal.
//
{
  if ( n < 0 || l < n )
  {
    cout << " \n";
    cout << "HYPERGEOMETRIC_CHECK - Fatal error!\n";
    cout << "  Input N is out of range.\n";
    return false;
  }

  if ( m < 0 || l < m )
  {
    cout << " \n";
    cout << "HYPERGEOMETRIC_CHECK - Fatal error!\n";
    cout << "  Input M is out of range.\n";
    return false;
  }

  if ( l < 0 )
  {
    cout << " \n";
    cout << "HYPERGEOMETRIC_CHECK - Fatal error!\n";
    cout << "  Input L is out of range.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double hypergeometric_mean ( int n, int m, int l )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_MEAN returns the mean of the Hypergeometric PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of balls selected.
//    0 <= N <= L.
//
//    Input, int M, the number of white balls in the population.
//    0 <= M <= L.
//
//    Input, int L, the number of balls to select from.
//    0 <= L.
//
//    Output, double HYPERGEOMETRIC_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = ( double ) ( n * m ) / ( double ) ( l );

  return mean;
}
//****************************************************************************80

double hypergeometric_pdf ( int x, int n, int m, int l )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_PDF evaluates the Hypergeometric PDF.
//
//  Discussion:
//
//    PDF(N,M,L;X) = C(M,X) * C(L-M,N-X) / C(L,N).
//
//    PDF(N,M,L;X) is the probability of drawing X white balls in a
//    single random sample of size N from a population containing
//    M white balls and a total of L balls.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the desired number of white balls.
//    0 <= X <= N, usually, although any value of X can be given.
//
//    Input, int N, the number of balls selected.
//    0 <= N <= L.
//
//    Input, int M, the number of white balls in the population.
//    0 <= M <= L.
//
//    Input, int L, the number of balls to select from.
//    0 <= L.
//
//    Output, double HYPERGEOMETRIC_PDF, the probability of exactly K white balls.
//
{
  double c1;
  double c2;
  double c3;
  double pdf;
  double pdf_log;
//
//  Special cases.
//
  if ( x < 0 )
  {
    pdf = 1.0;
  }
  else if ( n < x )
  {
    pdf = 0.0;
  }
  else if ( m < x )
  {
    pdf = 0.0;
  }
  else if ( l < x )
  {
    pdf = 0.0;
  }
  else if ( n == 0 )
  {
    if ( x == 0 )
    {
      pdf = 1.0;
    }
    else
    {
      pdf = 0.0;
    }
  }
  else
  {
    c1 = binomial_coef_log ( m, x );
    c2 = binomial_coef_log ( l-m, n-x );
    c3 = binomial_coef_log ( l, n );

    pdf_log = c1 + c2 - c3;

    pdf = exp ( pdf_log );

  }

  return pdf;
}
//****************************************************************************80

int hypergeometric_sample ( int n, int m, int l, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_SAMPLE samples the Hypergeometric PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jerry Banks, editor,
//    Handbook of Simulation,
//    Engineering and Management Press Books, 1998, page 165.
//
//  Parameters:
//
//    Input, int N, the number of balls selected.
//    0 <= N <= L.
//
//    Input, int M, the number of white balls in the population.
//    0 <= M <= L.
//
//    Input, int L, the number of balls to select from.
//    0 <= L.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int HYPERGEOMETRIC_SAMPLE, a sample of the PDF.
//
{
  double a;
  double b;
  double c1_log;
  double c2_log;
  double u;
  int x;

  c1_log = binomial_coef_log ( l - m, n );
  c2_log = binomial_coef_log ( l, n );

  a = exp ( c1_log - c2_log );
  b = a;

  u = r8_uniform_01 ( seed );

  x = 0;

  while ( a < u )
  {
    b = b * ( double ) ( ( m - x ) * ( n - x ) )
      / ( double ) ( ( x + 1 ) * ( l - m - n + x + 1 ) );

    a = a + b;

    x = x + 1;

  }

  return x;
}
//****************************************************************************80

double hypergeometric_variance ( int n, int m, int l )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_VARIANCE returns the variance of the Hypergeometric PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of balls selected.
//    0 <= N <= L.
//
//    Input, int M, the number of white balls in the population.
//    0 <= M <= L.
//
//    Input, int L, the number of balls to select from.
//    0 <= L.
//
//    Output, double HYPERGEOMETRIC_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( double ) ( n * m * ( l - m ) * ( l - n ) )
    / ( double ) ( l * l * ( l - 1 ) );

  return variance;
}
//****************************************************************************80

double i4_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL returns N!.
//
//  Discussion:
//
//    N! = Product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 1999
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
//    Output, double I4_FACTORIAL, the factorial of N.
//
{
  double fact;
  int i;
//
//  Check.
//
  if ( n < 0 )
  {
    cout << "\n";
    cout << "I4_FACTORIAL - Fatal error!\n";
    cout << "  N < 0.\n";
    return 0;
  }

  fact = 1.0;

  for ( i = 2; i <= n; i++ )
  {
    fact = fact * ( double ) i;
  }

  return fact;
}
//****************************************************************************80

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4
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
//    Output, int I4_HUGE, a "huge" integer.
//
{
  return 2147483647;
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
//
{
  if ( i2 < i1 )
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

int *i4row_max ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4ROW_MAX returns the maximums of an I4ROW.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    MAX =
//      3
//      7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int A[M*N], the array to be examined.
//
//    Output, int I4ROW_AMAX[M], the maximums of the rows.
//
{
  int i;
  int j;
  int *amax;

  amax = new int[m];

  for ( i = 0; i < m; i++ )
  {
    amax[i] = a[i+0*m];

    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[i+j*m] )
      {
        amax[i] = a[i+j*m];
      }
    }
  }

  return amax;
}
//****************************************************************************80

double *i4row_mean ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4ROW_MEAN returns the means of an I4ROW.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    MEAN =
//      2
//      5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int A[M*N], the array to be examined.
//
//    Output, double I4ROW_MEAN[M], the means, or averages, of the rows.
//
{
  int i;
  int j;
  double *mean;

  mean = new double[m];

  for ( i = 0; i < m; i++ )
  {
    mean[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      mean[i] = mean[i] + ( double ) a[i+j*m];
    }
    mean[i] = mean[i] / ( double ) ( n );
  }

  return mean;
}
//****************************************************************************80

int *i4row_min ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4ROW_MIN returns the minimums of an I4ROW.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    MIN =
//      1
//      2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int A[M*N], the array to be examined.
//
//    Output, int I4ROW_AMIN[M], the minimums of the rows.
//
{
  int i;
  int j;
  int *amin;

  amin = new int[m];

  for ( i = 0; i < m; i++ )
  {
    amin[i] = a[i+0*m];
    for ( j = 1; j < n; j++ )
    {
      if ( a[i+j*m] < amin[i] )
      {
        amin[i] = a[i+j*m];
      }
    }
  }

  return amin;
}
//****************************************************************************80

double *i4row_variance ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4ROW_VARIANCE returns the variances of an I4ROW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, int A[M*N], the array whose variances are desired.
//
//    Output, double I4ROW_VARIANCE[M], the variances of the rows.
//
{
  int i;
  int j;
  double mean;
  double *variance;

  variance = new double[m];

  for ( i = 0; i < m; i++ )
  {
    mean = 0.0;
    for ( j = 0; j < n; j++ )
    {
      mean = mean + ( double ) a[i+j*m];
    }
    mean = mean / ( double ) ( n );

    variance[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      variance[i] = variance[i] + pow ( ( ( double ) a[i+j*m] - mean ), 2 );
    }

    if ( 1 < n )
    {
      variance[i] = variance[i] / ( double ) ( n - 1 );
    }
    else
    {
      variance[i] = 0.0;
    }
  }

  return variance;
}
//****************************************************************************80

int i4vec_max ( int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MAX returns the maximum of an I4VEC
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int X[N], the vector whose maximum is desired.
//
//    Output, int I4VEC_MAX, the maximum of the vector entries.
//
{
  int i;
  int value;

  value = x[0];
  for ( i = 1; i < n; i++ )
  {
    if ( value < x[i] )
    {
      value = x[i];
    }
  }
  return value;
}
//****************************************************************************80

double i4vec_mean ( int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MEAN returns the mean of an I4VEC.
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
//    Input, int X[N], the vector whose mean is desired.
//
//    Output, double I4VEC_MEAN, the mean, or average, of the vector entries.
//
{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ )
  {
    mean = mean + ( double ) x[i];
  }

  mean = mean / ( double ) n;

  return mean;
}
//****************************************************************************80

int i4vec_min ( int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MIN returns the minimum of an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int X[N], the vector whose minimum is desired.
//
//    Output, int I4VEC_MIN, the minimum of the vector entries.
//
{
  int i;
  int value;

  value = x[0];
  for ( i = 1; i < n; i++ )
  {
    if ( x[i] < value )
    {
      value = x[i];
    }
  }
  return value;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2003
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
//    Input, char *TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( s_len_trim ( title ) != 0 )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(6) << i << "  " << setw(8) << a[i] << "\n";
  }
  return;
}
//****************************************************************************80

int i4vec_run_count ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_RUN_COUNT counts runs of equal values in an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of integer values.
//
//    A run is a sequence of equal values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be examined.
//
//    Output, int I4VEC_RUN_COUNT, the number of runs.
//
{
  int i;
  int run_count;
  int test;

  run_count = 0;

  if ( n < 1 )
  {
    return run_count;
  }

  test = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 || a[i] != test )
    {
      run_count = run_count + 1;
      test = a[i];
    }
  }

  return run_count;
}
//****************************************************************************80

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
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

double i4vec_variance ( int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_VARIANCE returns the variance of an I4VEC.
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
//    Input, int X[N], the vector whose variance is desired.
//
//    Output, double I4VEC_VARIANCE, the variance of the vector entries.
//
{
  int i;
  double mean;
  double variance;

  mean = i4vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ )
  {
    variance = variance + ( ( double ) x[i] - mean ) * ( ( double ) x[i] - mean );
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

double inverse_gaussian_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_CDF evaluates the Inverse Gaussian CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//    0.0 < X.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double INVERSE_GAUSSIAN_CDF, the value of the CDF.
//
{
  double cdf;
  double cdf1;
  double cdf2;
  double x1;
  double x2;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    x1 = sqrt ( b / x ) * ( x - a ) / a;
    cdf1 = normal_01_cdf ( x1 );

    x2 = - sqrt ( b / x ) * ( x + a ) / a;
    cdf2 = normal_01_cdf ( x2 );

    cdf = cdf1 + exp ( 2.0 * b / a ) * cdf2;
  }
  return cdf;
}
//****************************************************************************80

bool inverse_gaussian_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_CHECK checks the parameters of the Inverse Gaussian CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, bool INVERSE_GAUSSIAN_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << " \n";
    cout << "INVERSE_GAUSSIAN_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "INVERSE_GAUSSIAN_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double inverse_gaussian_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_MEAN returns the mean of the Inverse Gaussian PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double INVERSE_GAUSSIAN_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double inverse_gaussian_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_PDF evaluates the Inverse Gaussian PDF.
//
//  Discussion:
//
//    The Inverse Gaussian PDF is also known as the Wald PDF
//    and the Inverse Normal PDF.
//
//    PDF(A,B;X)
//      = sqrt ( B / ( 2 * PI * X^3 ) )
//        * exp ( - B * ( X - A )^2 / ( 2.0 * A^2 * X ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 < X
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double INVERSE_GAUSSIAN_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  if ( x <= 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = sqrt ( b / ( 2.0 * pi * pow ( x, 3 ) ) ) *
      exp ( - b * ( x - a ) * ( x - a ) / ( 2.0 * a * a * x ) );
  }

  return pdf;
}
//****************************************************************************80

double inverse_gaussian_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_SAMPLE samples the Inverse Gaussian PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double INVERSE_GAUSSIAN_SAMPLE, a sample of the PDF.
//
{
  double phi;
  double t;
  double u;
  double x;
  double y;
  double z;

  phi = b / a;
  z = normal_01_sample ( seed );
  y = z * z;

  t = 1.0 + 0.5 * ( y - sqrt ( 4.0 * phi * y + y * y ) ) / phi;
  u = r8_uniform_01 ( seed );

  if ( u * ( 1.0 + t ) <= 1.0 )
  {
    x = a * t;
  }
  else
  {
    x = a / t;
  }

  return x;
}
//****************************************************************************80

double inverse_gaussian_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_VARIANCE returns the variance of the Inverse Gaussian PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double INVERSE_GAUSSIAN_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = a * a * a / b;

  return variance;
}
//****************************************************************************80

void laplace_cdf_values ( int *n_data, double *mu, double *beta, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_CDF_VALUES returns some values of the Laplace CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = LaplaceDistribution [ mu, beta ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
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
//    Output, double *MU, the mean of the distribution.
//
//    Output, double *BETA, the shape parameter.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double beta_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  double fx_vec[N_MAX] = {
     0.5000000000000000,
     0.8160602794142788,
     0.9323323583816937,
     0.9751064658160680,
     0.6967346701436833,
     0.6417343447131054,
     0.6105996084642976,
     0.5906346234610091,
     0.5000000000000000,
     0.3032653298563167,
     0.1839397205857212,
     0.1115650800742149 };

  double mu_vec[N_MAX] = {
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01 };

  double x_vec[N_MAX] = {
     0.0000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *mu = 0.0;
    *beta = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *mu = mu_vec[*n_data-1];
    *beta = beta_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double laplace_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_CDF evaluates the Laplace CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LAPLACE_CDF, the value of the PDF.
//
{
  double cdf;
  double y;

  y = ( x - a ) / b;

  if ( x <= a )
  {
    cdf = 0.5 * exp ( y );
  }
  else
  {
    cdf = 1.0 - 0.5 * exp ( - y );
  }
  return cdf;
}
//****************************************************************************80

double laplace_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_CDF_INV inverts the Laplace CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LAPLACE_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "LAPLACE_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf <= 0.5 )
  {
    x = a + b * log ( 2.0 * cdf );
  }
  else
  {
    x = a - b * log ( 2.0 * ( 1.0 - cdf ) );
  }

  return x;
}
//****************************************************************************80

bool laplace_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_CHECK checks the parameters of the Laplace PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, bool LAPLACE_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "LAPLACE_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }
  return true;
}
//****************************************************************************80

double laplace_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_MEAN returns the mean of the Laplace PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LAPLACE_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double laplace_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_PDF evaluates the Laplace PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = exp ( - abs ( X - A ) / B ) / ( 2 * B )
//
//  Discussion:
//
//    The Laplace PDF is also known as the Double Exponential PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LAPLACE_PDF, the value of the PDF.
//
{
  double pdf;

  pdf = exp ( - r8_abs ( x - a ) / b ) / ( 2.0 * b );

  return pdf;
}
//****************************************************************************80

double laplace_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_SAMPLE samples the Laplace PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double LAPLACE_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = laplace_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double laplace_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_VARIANCE returns the variance of the Laplace PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LAPLACE_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = 2.0 * b * b;

  return variance;
}
//****************************************************************************80

double lerch ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    LERCH estimates the Lerch transcendent function.
//
//  Discussion:
//
//    The Lerch transcendent function is defined as:
//
//      LERCH ( A, B, C ) = Sum ( 0 <= K < Infinity ) A**K / ( C + K )**B
//
//    excluding any term with ( C + K ) = 0.
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
//  Reference:
//
//    Eric Weisstein, editor,
//    CRC Concise Encylopedia of Mathematics,
//    CRC Press, 1998.
//
//  Thanks:
//
//    Oscar van Vlijmen
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the function.
//
//    Output, double LERCH, an approximation to the Lerch
//    transcendent function.
//
{
  double a_k;
  int k;
  double sum2;
  double sum2_old;

  sum2 = 0.0;
  k = 0;
  a_k = 1.0;

  for ( ; ; )
  {
    sum2_old = sum2;

    if ( c + ( double ) ( k ) == 0.0 )
    {
      k = k + 1;
      a_k = a_k * a;
      continue;
    }

    sum2 = sum2 + a_k / pow ( c + ( double ) ( k ), b );

    if ( sum2 <= sum2_old )
    {
      break;
    }

    k = k + 1;
    a_k = a_k * a;
  }
  return sum2;
}
//****************************************************************************80

double levy_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LEVY_CDF evaluates the Levy CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//    Normally, A <= X.
//
//    Input, double A, B, the parameters of the PDF.
//    0 < B.
//
//    Output, double LEVY_CDF, the value of the PDF.
//
{
  double cdf;

  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "LEVY_CDF - Fatal error!\n";
    cout << "  B <= 0.0.\n";
    exit ( 1 );
  }

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 1.0 - error_f ( sqrt ( b / ( 2.0 * ( x - a ) ) ) );
  }

  return cdf;
}
//****************************************************************************80

double levy_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LEVY_CDF_INV inverts the Levy CDF.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0 < B.
//
//    Output, double X, the corresponding argument.
//
{
  double cdf1;
  double x;
  double x1;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "LEVY_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "LEVY_CDF_INV - Fatal error!\n";
    cout << "  B <= 0.0.\n";
    exit ( 1 );
  }

  cdf1 = 1.0 - 0.5 * cdf;
  x1 = normal_01_cdf_inv ( cdf1 );
  x = a + b / ( x1 * x1 );

  return x;
}
//****************************************************************************80

double levy_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LEVY_PDF evaluates the Levy PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = sqrt ( B / ( 2 * PI ) )
//               * exp ( - B / ( 2 * ( X - A ) )
//               / ( X - A )^(3/2)
//
//    for A <= X.
//
//    Note that the Levy PDF does not have a finite mean or variance.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    Normally, A <= X.
//
//    Input, double A, B, the parameters of the PDF.
//    0 < B.
//
//    Output, double LEVY_PDF, the value of the PDF.
//
{
  double pdf;
  double pi = 3.141592653589793;

  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "LEVY_PDF - Fatal error!\n";
    cout << "  B <= 0.0.\n";
    exit ( 1 );
  }

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = sqrt ( b / ( 2.0 * pi ) )
        * exp ( - b / ( 2.0 * ( x - a ) ) )
        / sqrt ( pow ( x - a, 3 ) );
  }

  return pdf;
}
//****************************************************************************80

double levy_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    LEVY_SAMPLE samples the Levy PDF.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double LEVY_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = levy_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double log_normal_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CDF evaluates the Lognormal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 < X.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;
  double logx;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    logx = log ( x );

    cdf = normal_cdf ( logx, a, b );
  }

  return cdf;
}
//****************************************************************************80

double log_normal_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CDF_INV inverts the Lognormal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input, double X, the corresponding argument.
//
{
  double logx;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "LOG_NORMAL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  logx = normal_cdf_inv ( cdf, a, b );

  x = exp ( logx );

  return x;
}
//****************************************************************************80

void log_normal_cdf_values ( int *n_data, double *mu, double *sigma,
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CDF_VALUES returns some values of the Log Normal CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = LogNormalDistribution [ mu, sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
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
//    Output, double *MU, the mean of the distribution.
//
//    Output, double *SIGMA, the shape parameter of the distribution.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double fx_vec[N_MAX] = {
     0.2275013194817921E-01,
     0.2697049307349095,
     0.5781741008028732,
     0.7801170895122241,
     0.4390310097476894,
     0.4592655190218048,
     0.4694258497695908,
     0.4755320473858733,
     0.3261051056816658,
     0.1708799040927608,
     0.7343256357952060E-01,
     0.2554673736161761E-01 };

  double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  double sigma_vec[N_MAX] = {
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *mu = 0.0;
    *sigma = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *mu = mu_vec[*n_data-1];
    *sigma = sigma_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool log_normal_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CHECK checks the parameters of the Lognormal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, bool LOG_NORMAL_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "LOG_NORMAL_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }
  return true;
}
//****************************************************************************80

double log_normal_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_MEAN returns the mean of the Lognormal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = exp ( a + 0.5 * b * b );

  return mean;
}
//****************************************************************************80

double log_normal_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_PDF evaluates the Lognormal PDF.
//
//  Discussion:
//
//    PDF(A,B;X)
//      = exp ( - 0.5 * ( ( log ( X ) - A ) / B )^2 )
//        / ( B * X * sqrt ( 2 * PI ) )
//
//    The Lognormal PDF is also known as the Cobb-Douglas PDF,
//    and as the Antilog_normal PDF.
//
//    The Lognormal PDF describes a variable X whose logarithm
//    is normally distributed.
//
//    The special case A = 0, B = 1 is known as Gilbrat's PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 < X
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;
  double y;

  if ( x <= 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    y = ( log ( x ) - a ) / b;
    pdf = exp ( -0.5 * y * y ) / ( b * x * sqrt ( 2.0 * pi ) );
  }

  return pdf;
}
//****************************************************************************80

double log_normal_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_SAMPLE samples the Lognormal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = log_normal_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double log_normal_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_VARIANCE returns the variance of the Lognormal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = exp ( 2.0 * a + b * b ) * ( exp ( b * b ) - 1.0 );

  return variance;
}
//****************************************************************************80

double log_series_cdf ( double x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_CDF evaluates the Logarithmic Series CDF.
//
//  Discussion:
//
//    Simple summation is used, with a recursion to generate successive
//    values of the PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Thanks:
//
//    Oscar van Vlijmen
//
//  Parameters:
//
//    Input, int X, the argument of the PDF.
//    0 < X
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A < 1.0.
//
//    Output, double LOG_SERIES_CDF, the value of the CDF.
//
{
  double cdf;
  double pdf;
  int x2;

  cdf = 0.0;

  for ( x2 = 1; x2 <= x; x2++ )
  {
    if ( x2 == 1 )
    {
      pdf = - a / log ( 1.0 - a );
    }
    else
    {
      pdf = ( double ) ( x2 - 1 ) * a * pdf / ( double ) ( x2 );
    }

    cdf = cdf + pdf;
  }
  return cdf;
}
//****************************************************************************80

int log_series_cdf_inv ( double cdf, double a )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_CDF_INV inverts the Logarithmic Series CDF.
//
//  Discussion:
//
//    Simple summation is used.  The only protection against an
//    infinite loop caused by roundoff is that X cannot be larger
//    than 1000.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A < 1.0.
//
//    Output, double LOG_SERIES_CDF_INV, the argument X of the CDF for which
//    CDF(X-1) <= CDF <= CDF(X).
//
{
  double cdf2;
  double pdf;
  int x;
  int xmax = 1000;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "LOG_SERIES_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  cdf2 = 0.0;
  x = 1;

  while ( cdf2 < cdf && x < xmax )
  {
    if ( x == 1 )
    {
      pdf = - a / log ( 1.0 - a );
    }
    else
    {
      pdf = ( double ) ( x - 1 ) * a * pdf / ( double ) ( x );
    }

    cdf2 = cdf2 + pdf;

    x = x + 1;
  }
  return x;
}
//****************************************************************************80

void log_series_cdf_values ( int *n_data, double *t, int *n, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_CDF_VALUES returns some values of the log series CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = LogSeriesDistribution [ t ]
//      CDF [ dist, n ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2004
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
//    Output, double *T, the parameter of the function.
//
//    Output, int *N, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 29

  double fx_vec[N_MAX] = {
     0.9491221581029903,
     0.9433541128559735,
     0.9361094611773272,
     0.9267370278044118,
     0.9141358246245129,
     0.8962840235449100,
     0.8690148741955517,
     0.8221011541254772,
     0.7213475204444817,
     0.6068261510845583,
     0.5410106403333613,
     0.4970679476476894,
     0.4650921887927060,
     0.4404842934597863,
     0.4207860535926143,
     0.4045507673897055,
     0.3908650337129266,
     0.2149757685421097,
     0.0000000000000000,
     0.2149757685421097,
     0.3213887739704539,
     0.3916213575531612,
     0.4437690508633213,
     0.4850700239649681,
     0.5191433267738267,
     0.5480569580144867,
     0.5731033910767085,
     0.5951442521714636,
     0.6147826594068904 };

  int n_vec[N_MAX] = {
     1, 1, 1, 1, 1,
     1, 1, 1, 1, 1,
     1, 1, 1, 1, 1,
     1, 1, 1, 0, 1,
     2, 3, 4, 5, 6,
     7, 8, 9, 10 };

  double t_vec[N_MAX] = {
     0.1000000000000000,
     0.1111111111111111,
     0.1250000000000000,
     0.1428571428571429,
     0.1666666666666667,
     0.2000000000000000,
     0.2500000000000000,
     0.3333333333333333,
     0.5000000000000000,
     0.6666666666666667,
     0.7500000000000000,
     0.8000000000000000,
     0.8333333333333333,
     0.8571485714857149,
     0.8750000000000000,
     0.8888888888888889,
     0.9000000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000,
     0.9900000000000000  };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *t = 0.0;
    *n = 0;
    *fx = 0.0;
  }
  else
  {
    *t = t_vec[*n_data-1];
    *n = n_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool log_series_check ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_CHECK checks the parameter of the Logarithmic Series PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A < 1.0.
//
//    Output, bool LOG_SERIES_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 || 1.0 <= a )
  {
    cout << " \n";
    cout << "LOG_SERIES_CHECK - Fatal error!\n";
    cout << "  A <= 0.0 or 1.0 <= A\n";
    return false;
  }
  return true;
}
//****************************************************************************80

double log_series_mean ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_MEAN returns the mean of the Logarithmic Series PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A < 1.0.
//
//    Output, double LOG_SERIES_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = - a / ( ( 1.0 - a ) * log ( 1.0 - a ) );

  return mean;
}
//****************************************************************************80

double log_series_pdf ( int x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_PDF evaluates the Logarithmic Series PDF.
//
//  Discussion:
//
//    PDF(A;X) = - A**X / ( X * log ( 1 - A ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the argument of the PDF.
//    0 < X
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A < 1.0.
//
//    Output, double LOG_SERIES_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x <= 0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = - pow ( a, x ) / ( ( double ) ( x ) * log ( 1.0 - a ) );
  }
  return pdf;
}
//****************************************************************************80

int log_series_sample ( double a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_SAMPLE samples the Logarithmic Series PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer-Verlag, New York, 1986, page 547.
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A < 1.0.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int LOG_SERIES_SAMPLE, a sample of the PDF.
//
{
  double u;
  double v;
  int x;

  u = r8_uniform_01 ( seed );
  v = r8_uniform_01 ( seed );

  x = ( int ) ( 1.0 + log ( v ) / ( log ( 1.0 - pow ( 1.0 - a, u ) ) ) );

  return x;
}
//****************************************************************************80

double log_series_variance ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_VARIANCE returns the variance of the Logarithmic Series PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A < 1.0.
//
//    Output, double LOG_SERIES_VARIANCE, the variance of the PDF.
//
{
  double alpha;
  double variance;

  alpha = - 1.0 / log ( 1.0 - a );

  variance = a * alpha * ( 1.0 - alpha * a ) / pow ( ( 1.0 - a ), 2 );

  return variance;
}
//****************************************************************************80

double log_uniform_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_CDF evaluates the Log Uniform CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else if ( x < b )
  {
    cdf = ( log ( x ) - log ( a ) ) / ( log ( b ) - log ( a ) );
  }
  else
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double log_uniform_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_CDF_INV inverts the Log Uniform CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double X, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "LOG_UNIFORM_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a * exp ( ( log ( b ) - log ( a ) ) * cdf );

  return x;
}
//****************************************************************************80

bool log_uniform_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_CHECK checks the parameters of the Log Uniform CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    1.0 < A < B.
//
//    Output, bool LOG_UNIFORM_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 1.0 )
  {
    cout << "\n";
    cout << "LOG_UNIFORM_CHECK - Fatal error!\n";
    cout << "  A <= 1.\n";
    return false;
  }

  if ( b <= a )
  {
    cout << "\n";
    cout << "LOG_UNIFORM_CHECK - Fatal error!\n";
    cout << "  B <= A.\n";
    return false;
  }
  return true;
}
//****************************************************************************80

double log_uniform_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_MEAN returns the mean of the Log Uniform PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    1.0 < A < B.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = ( b - a ) / ( log ( b ) - log ( a ) );

  return mean;
}
//****************************************************************************80

double log_uniform_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_PDF evaluates the Log Uniform PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = ( log ( B ) - log ( A ) ) / X  for A <= X <= B
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    1.0 < A < B.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < a )
  {
    pdf = 0.0;
  }
  else if ( x <= b )
  {
    pdf = 1.0 / ( x * ( log ( b ) - log ( a ) ) );
  }
  else
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

double log_uniform_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_SAMPLE samples the Log Uniform PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    1.0 < A < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = log_uniform_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double log_uniform_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_VARIANCE returns the variance of the Log Uniform PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    1.0 < A < B.
//
//    Output, double LOG_UNIFORM_VARIANCE, the variance of the PDF.
//
{
  double mean;
  double variance;

  mean = log_uniform_mean ( a, b );

  variance =
    ( ( 0.5 * b * b - 2.0 * mean * b + mean * mean * log ( b ) )
    - ( 0.5 * a * a - 2.0 * mean * a + mean * mean * log ( a ) ) )
    / ( log ( b ) - log ( a ) );

  return variance;
}
//****************************************************************************80

double logistic_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_CDF evaluates the Logistic CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LOGISTIC_CDF, the value of the CDF.
//
{
  double cdf;

  cdf = 1.0 / ( 1.0 + exp ( ( a - x ) / b ) );

  return cdf;
}
//****************************************************************************80

double logistic_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_CDF_INV inverts the Logistic CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LOGISTIC_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "LOGISTIC_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a - b * log ( ( 1.0 - cdf ) / cdf );

  return x;
}
//****************************************************************************80

void logistic_cdf_values ( int *n_data, double *mu, double *beta, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_CDF_VALUES returns some values of the Logistic CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = LogisticDistribution [ mu, beta ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2004
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
//    Output, double *MU, the mean of the distribution.
//
//    Output, double *BETA, the shape parameter of the distribution.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double beta_vec[N_MAX] = {
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  double fx_vec[N_MAX] = {
     0.5000000000000000,
     0.8807970779778824,
     0.9820137900379084,
     0.9975273768433652,
     0.6224593312018546,
     0.5825702064623147,
     0.5621765008857981,
     0.5498339973124779,
     0.6224593312018546,
     0.5000000000000000,
     0.3775406687981454,
     0.2689414213699951 };

  double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *mu = 0.0;
    *beta = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *mu = mu_vec[*n_data-1];
    *beta = beta_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool logistic_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_CHECK checks the parameters of the Logistic CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, bool LOGISTIC_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "LOGISTIC_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }
  return true;
}
//****************************************************************************80

double logistic_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_MEAN returns the mean of the Logistic PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LOGISTIC_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double logistic_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_PDF evaluates the Logistic PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = exp ( ( A - X ) / B ) /
//      ( B * ( 1 + exp ( ( A - X ) / B ) )^2 )
//
//    The Logistic PDF is also known as the Sech-Squared PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LOGISTIC_PDF, the value of the PDF.
//
{
  double pdf;
  double temp;

  temp = exp ( ( a - x ) / b );

  pdf = temp / ( b * ( 1.0 + temp ) * ( 1.0 + temp ) );

  return pdf;
}
//****************************************************************************80

double logistic_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_SAMPLE samples the Logistic PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double LOGISTIC_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = logistic_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double logistic_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_VARIANCE returns the variance of the Logistic PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double LOGISTIC_VARIANCE, the variance of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double variance;

  variance = pi * pi * b * b / 3.0;

  return variance;
}
//****************************************************************************80

double lorentz_cdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_CDF evaluates the Lorentz CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Output, double LORENTZ_CDF, the value of the CDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;

  cdf = 0.5 + atan ( x ) / pi;

  return cdf;
}
//****************************************************************************80

double lorentz_cdf_inv ( double cdf )

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_CDF_INV inverts the Lorentz CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Output, double LORENTZ_CDF_INV, the corresponding argument.
//
{
  const double pi = 3.14159265358979323;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "LORENTZ_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = tan ( pi * ( cdf - 0.5 ) );

  return x;
}
//****************************************************************************80

double lorentz_mean ( )

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_MEAN returns the mean of the Lorentz PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double LORENTZ_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = 0.0;

  return mean;
}
//****************************************************************************80

double lorentz_pdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_PDF evaluates the Lorentz PDF.
//
//  Discussion:
//
//    PDF(X) = 1 / ( PI * ( 1 + X^2 ) )
//
//    The chief interest of the Lorentz PDF is that it is easily
//    inverted, and can be used to dominate other PDF's in an
//    acceptance/rejection method.
//
//    LORENTZ_PDF(X) = CAUCHY_PDF(0,1;X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Output, double LORENTZ_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  pdf = 1.0 / ( pi * ( 1.0 + x * x ) );

  return pdf;
}
//****************************************************************************80

double lorentz_sample ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_SAMPLE samples the Lorentz PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double LORENTZ_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = lorentz_cdf_inv ( cdf );

  return x;
}
//****************************************************************************80

double lorentz_variance ( )

//****************************************************************************80
//
//  Purpose:
//
//    LORENTZ_VARIANCE returns the variance of the Lorentz PDF.
//
//  Discussion:
//
//    The variance of the Lorentz PDF is not well defined.  This routine
//    is made available for completeness only, and simply returns
//    a "very large" number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double LORENTZ_VARIANCE, the mean of the PDF.
//
{
  double variance;

  variance = r8_huge ( );

  return variance;
}
//****************************************************************************80

double maxwell_cdf ( double x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_CDF evaluates the Maxwell CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double MAXWELL_CDF, the value of the CDF.
//
{
  double cdf;
  double p2;
  double x2;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    x2 = x / a;
    p2 = 1.5;

    cdf = gamma_inc ( p2, x2 );
  }

  return cdf;
}
//****************************************************************************80

double maxwell_cdf_inv ( double cdf, double a )

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_CDF_INV inverts the Maxwell CDF.
//
//  Discussion:
//
//    A simple bisection method is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double MAXWELL_CDF_INV, the corresponding argument of the CDF.
//
{
  double cdf1;
  double cdf2;
  double cdf3;
  int it;
  int it_max = 100;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "MAXWELL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = 0.0;
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = r8_huge ( );
    return x;
  }

  x1 = 0.0;
  cdf1 = 0.0;

  x2 = 1.0;

  for ( ; ; )
  {
    cdf2 = maxwell_cdf ( x2, a );

    if ( cdf < cdf2 )
    {
      break;
    }

    x2 = 2.0 * x2;

    if ( 1000000.0 < x2 )
    {
      cout << " \n";
      cout << "MAXWELL_CDF_INV - Fatal error!\n";
      cout << "  Initial bracketing effort fails.\n";
      exit ( 1 );
    }
  }
//
//  Now use bisection.
//
  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = maxwell_cdf ( x3, a );

    if ( r8_abs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      break;
    }

    if ( it_max < it )
    {
      cout << " \n";
      cout << "MAXWELL_CDF_INV - Fatal error!\n";
      cout << "  Iteration limit exceeded.\n";
      exit ( 1 );
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
      cdf2 = cdf3;
    }
  }

  return x;
}
//****************************************************************************80

bool maxwell_check ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_CHECK checks the parameters of the Maxwell CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, bool MAXWELL_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << " \n";
    cout << "MAXWELL_CHECK - Fatal error!\n";
    cout << "  A <= 0.0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double maxwell_mean ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_MEAN returns the mean of the Maxwell PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double MAXWELL_MEAN, the mean value.
//
{
  double mean;

  mean = sqrt ( 2.0 ) * a * r8_gamma ( 2.0 ) / r8_gamma ( 1.5 );

  return mean;
}
//****************************************************************************80

double maxwell_pdf ( double x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_PDF evaluates the Maxwell PDF.
//
//  Discussion:
//
//    PDF(A;X) = EXP ( - 0.5 * ( X / A )^2 ) * ( X / A )^2 /
//      ( sqrt ( 2 ) * A * GAMMA ( 1.5 ) )
//
//    MAXWELL_PDF(A;X) = CHI_PDF(0,A,3;X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0 < X
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double MAXWELL_PDF, the value of the PDF.
//
{
  double pdf;
  double y;

  if ( x <= 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    y = x / a;

    pdf = exp ( - 0.5 * y * y ) * y * y / ( sqrt ( 2.0 ) * a * r8_gamma ( 1.5 ) );
  }

  return pdf;
}
//****************************************************************************80

double maxwell_sample ( double a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_SAMPLE samples the Maxwell PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double MAXWELL_SAMPLE, a sample of the PDF.
//
{
  double a2;
  double x;

  a2 = 3.0;
  x = chi_square_sample ( a2, seed );

  x = a * sqrt ( x );

  return x;
}
//****************************************************************************80

double maxwell_variance ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_VARIANCE returns the variance of the Maxwell PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double MAXWELL_VARIANCE, the variance of the PDF.
//
{
  double temp;
  double variance;

  temp = r8_gamma ( 2.0 ) / r8_gamma ( 1.5 );

  variance = a * a * ( 3.0 - 2.0 * temp * temp );

  return variance;
}
//****************************************************************************80

bool multicoef_check ( int nfactor, int factor[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTICOEF_CHECK checks the parameters of the multinomial coefficient.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NFACTOR, the number of factors.
//    1 <= NFACTOR.
//
//    Input, int FACTOR(NFACTOR), contains the factors.
//    0.0 <= FACTOR(I).
//
//    Output, bool MULTICOEF_CHECK, is true if the parameters are legal.
//
{
  int i;

  if ( nfactor < 1 )
  {
    cout << " \n";
    cout << "MULTICOEF_CHECK - Fatal error!\n";
    cout << "  NFACTOR < 1.\n";
    return false;
  }

  for ( i = 0; i < nfactor; i++ )
  {
    if ( factor[i] < 0 )
    {
      cout << " \n";
      cout << "MULTICOEF_CHECK - Fatal error\n";
      cout << "  Factor[" << i << "] = " << factor[i] << "\n";
      cout << "  But this value must be nonnegative.\n";
      return false;
    }

  }

  return true;
}
//****************************************************************************80

int multinomial_coef1 ( int nfactor, int factor[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_COEF1 computes a Multinomial coefficient.
//
//  Discussion:
//
//    The multinomial coefficient is a generalization of the binomial
//    coefficient.  It may be interpreted as the number of combinations of
//    N objects, where FACTOR(1) objects are indistinguishable of type 1,
//    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
//    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
//
//    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
//
//    The log of the gamma function is used, to avoid overflow.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NFACTOR, the number of factors.
//    1 <= NFACTOR.
//
//    Input, int FACTOR(NFACTOR), contains the factors.
//    0.0 <= FACTOR(I).
//
//    Output, int MULTINOMIAL_COEF1, the value of the multinomial coefficient.
//
{
  double facn;
  int i;
  int n;
  int ncomb;

  if (  !multicoef_check ( nfactor, factor ) )
  {
    cout << " \n";
    cout << "MULTINOMIAL_COEF1 - Fatal error!\n";
    cout << "  MULTICOEF_CHECK failed.\n";
    ncomb = -i4_huge ( );
    return ncomb;
  }
//
//  The factors sum to N.
//
  n = i4vec_sum ( nfactor, factor );

  facn = factorial_log ( n );

  for ( i = 0; i < nfactor; i++ )
  {
    facn = facn - factorial_log ( factor[i] );
  }

  ncomb = r8_nint ( exp ( facn ) );

  return ncomb;
}
//****************************************************************************80

int multinomial_coef2 ( int nfactor, int factor[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_COEF2 computes a Multinomial coefficient.
//
//  Discussion:
//
//    The multinomial coefficient is a generalization of the binomial
//    coefficient.  It may be interpreted as the number of combinations of
//    N objects, where FACTOR(1) objects are indistinguishable of type 1,
//    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
//    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
//
//    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
//
//    A direct method is used, which should be exact.  However, there
//    is a possibility of intermediate overflow of the result.
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
//    Input, int NFACTOR, the number of factors.
//    1 <= NFACTOR.
//
//    Input, int FACTOR[NFACTOR], contains the factors.
//    0 <= FACTOR(I).
//
//    Output, int MULTINOMIAL_COEF2, the value of the multinomial coefficient.
//
{
  int i;
  int j;
  int k;
  int ncomb;

  if ( !multicoef_check ( nfactor, factor ) )
  {
    cout << " \n";
    cout << "MULTINOMIAL_COEF2 - Fatal error!\n";
    cout << "  MULTICOEF_CHECK failed.\n";
    ncomb = - i4_huge ( );
    return ncomb;
  }

  ncomb = 1;
  k = 0;

  for ( i = 0; i < nfactor; i++ )
  {
    for ( j = 1; j <= factor[i]; j++ )
    {
      k = k + 1;
      ncomb = ( ncomb * k ) / j;
    }
  }

  return ncomb;
}
//****************************************************************************80

bool multinomial_check ( int a, int b, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_CHECK checks the parameters of the Multinomial PDF.
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
//    Input, int A, the number of trials.
//
//    Input, int B, the number of outcomes possible on one trial.
//    1 <= B.
//
//    Input, double C(B).  C(I) is the probability of outcome I on
//    any trial.
//    0.0 <= C(I) <= 1.0,
//    Sum ( 1 <= I <= B ) C(I) = 1.0.
//
//    Output, bool MULTINOMIAL_CHECK, is true if the parameters are legal.
//
{
  double c_sum;
  int i;

  if ( b < 1 )
  {
    cout << " \n";
    cout << "MULTINOMIAL_CHECK - Fatal error!\n";
    cout << "  B < 1.\n";
    return false;
  }

  for ( i = 0; i < b; i++ )
  {
    if ( c[i] < 0.0 || 1.0 < c[i] )
    {
      cout << " \n";
      cout << "MULTINOMIAL_CHECK - Fatal error!\n";
      cout << "  Input C(I) is out of range.\n";
      return false;
    }
  }

  c_sum = r8vec_sum ( b, c );

  if ( 0.0001 < r8_abs ( 1.0 - c_sum ) )
  {
    cout << " \n";
    cout << "MULTINOMIAL_CHECK - Fatal error!\n";
    cout << "  The probabilities do not sum to 1.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double *multinomial_covariance ( int a, int b, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_COVARIANCE returns the covariances of the Multinomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the number of trials.
//
//    Input, int B, the number of outcomes possible on one trial.
//    1 <= B.
//
//    Input, double C(B).  C(I) is the probability of outcome I on
//    any trial.
//    0.0 <= C(I) <= 1.0,
//    SUM ( 1 <= I <= B) C(I) = 1.0.
//
//    Output, double MULTINOMIAL_COVARIANCE[B*B], the covariance matrix.
//
{
  double *covariance;
  int i;
  int j;

  covariance = new double[b*b];

  for ( i = 0; i < b; i++)
  {
    for ( j = 0; j < b; j++ )
    {
      if ( i == j )
      {
        covariance[i+j*b] = ( double ) ( a ) * c[i] * ( 1.0 - c[i] );
      }
      else
      {
        covariance[i+j*b] = - ( double ) ( a ) * c[i] * c[j];
      }
    }
  }

  return covariance;
}
//****************************************************************************80

double *multinomial_mean ( int a, int b, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_MEAN returns the means of the Multinomial PDF.
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
//    Input, int A, the number of trials.
//
//    Input, int B, the number of outcomes possible on one trial.
//    1 <= B.
//
//    Input, double C(B).  C(I) is the probability of outcome I on
//    any trial.
//    0.0 <= C(I) <= 1.0,
//    SUM ( 1 <= I <= B) C(I) = 1.0.
//
//    Output, double MEAN(B), MEAN(I) is the expected value of the
//    number of outcome I in N trials.
//
{
  int i;
  double *mean;

  mean = new double[b];

  for ( i = 0; i < b; i++ )
  {
    mean[i] = ( double ) ( a ) * c[i];
  }

  return mean;
}
//****************************************************************************80

double multinomial_pdf ( int x[], int a, int b, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_PDF computes a Multinomial PDF.
//
//  Discussion:
//
//    PDF(A,B,C;X) = Comb(A,B,X) * Product ( 1 <= I <= B ) C(I)**X(I)
//
//    where Comb(A,B,X) is the multinomial coefficient
//      C( A; X(1), X(2), ..., X(B) )
//
//    PDF(A,B,C;X) is the probability that in A trials there
//    will be exactly X(I) occurrences of event I, whose probability
//    on one trial is C(I), for I from 1 to B.
//
//    As soon as A or B gets large, the number of possible X's explodes,
//    and the probability of any particular X can become extremely small.
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
//    Input, int X[B]; X(I) counts the number of occurrences of
//    outcome I, out of the total of A trials.
//
//    Input, int A, the total number of trials.
//
//    Input, int B, the number of different possible outcomes on
//    one trial.
//
//    Input, double C[B]; C(I) is the probability of outcome I on
//    any one trial.
//
//    Output, double MULTINOMIAL_PDF, the value of the multinomial PDF.
//
{
  int i;
  double pdf;
  double pdf_log;
//
//  To try to avoid overflow, do the calculation in terms of logarithms.
//  Note that Gamma(A+1) = A factorial.
//
  pdf_log = gamma_log ( ( double ) ( a + 1 ) );

  for ( i = 0; i < b; i++ )
  {
    pdf_log = pdf_log + x[i] * log ( c[i] )
      - gamma_log ( ( double ) ( x[i] + 1 ) );
  }

  pdf = exp ( pdf_log );

  return pdf;
}
//****************************************************************************80

int *multinomial_sample ( int a, int b, double c[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_SAMPLE samples the Multinomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer-Verlag, New York, 1986, page 559.
//
//  Parameters:
//
//    Input, int A, the total number of trials.
//    0 <= A.
//
//    Input, int B, the number of outcomes possible on one trial.
//    1 <= B.
//
//    Input, double C[B].  C(I) is the probability of outcome I on
//    any trial.
//    0.0 <= C(I) <= 1.0,
//    SUM ( 1 <= I <= B) C(I) = 1.0.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int X[B]; X(I) is the number of
//    occurrences of event I during the N trials.
//
{
  int i;
  int ifactor;
  int ntot;
  double prob;
  double sum2;
  int *x;

  x = new int[b];
  ntot = a;

  sum2 = 1.0;

  for ( i = 0; i < b; i++ )
  {
    x[i] = 0;
  }

  for ( ifactor = 0; ifactor < b-1; ifactor++ )
  {
    prob = c[ifactor] / sum2;
//
//  Generate a binomial random deviate for NTOT trials with
//  single trial success probability PROB.
//
    x[ifactor] = binomial_sample ( ntot, prob, seed );

    ntot = ntot - x[ifactor];
    if ( ntot <= 0 )
    {
      return x;
    }

    sum2 = sum2 - c[ifactor];
  }
//
//  The last factor gets what's left.
//
  x[b-1] = ntot;

  return x;
}
//****************************************************************************80

double *multinomial_variance ( int a, int b, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_VARIANCE returns the variances of the Multinomial PDF.
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
//    Input, int A, the number of trials.
//
//    Input, int B, the number of outcomes possible on one trial.
//    1 <= B.
//
//    Input, double C[B].  C(I) is the probability of outcome I on
//    any trial.
//    0.0 <= C(I) <= 1.0,
//    SUM ( 1 <= I <= B ) C(I) = 1.0.
//
//    Output, double VARIANCE(B), VARIANCE(I) is the variance of the
//    total number of events of type I.
//
{
  int i;
  double *variance;

  variance = new double[b];

  for ( i = 0; i < b; i++ )
  {
    variance[i] = ( double ) ( a ) * c[i] * ( 1.0 - c[i] );
  }

  return variance;
}
//****************************************************************************80

double *multivariate_normal_sample ( int n, double mean[], double covar_factor[],
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIVARIATE_NORMAL_SAMPLE samples the Multivariate Normal PDF.
//
//  Discussion:
//
//    PDF ( Mean(1:N), S(1:N,1:N); X(1:N) ) = 1 / ( 2 * PI * det ( S ) )^(N/2)
//      * exp ( - ( X - Mean )' * inverse ( S ) * ( X - Mean ) / 2 )
//
//    Here,
//
//      X is the argument vector of length N,
//      Mean is the mean vector of length N,
//      S is an N by N positive definite symmetric covariance matrix.
//
//    The properties of S guarantee that it has a lower triangular
//    matrix L, the Cholesky factor, such that S = L * L'.  It is the
//    matrix L, rather than S, that is required by this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jerry Banks, editor,
//    Handbook of Simulation,
//    Engineering and Management Press Books, 1998, page 167.
//
//  Parameters:
//
//    Input, int N, the spatial dimension.
//
//    Input, double MEAN[N], the mean vector.
//
//    Input, double COVAR_FACTOR[N*N], the lower triangular Cholesky
//    factor L of the covariance matrix S.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double MULTIVARIATE_NORMAL_SAMPLE[N], a sample point
//    of the distribution.
//
{
  int i;
  int j;
  double *x;
  double z;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    z = normal_01_sample ( seed );

    x[i] = mean[i];

    for ( j = 0; j <= i; j++ )
    {
      x[i] = x[i] + covar_factor[i+j*n] * z;
    }
  }

  return x;
}
//****************************************************************************80

double nakagami_cdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    NAKAGAMI_CDF evaluates the Nakagami CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B
//    0.0 < C.
//
//    Output, double NAKAGAMI_CDF, the value of the CDF.
//
{
  double cdf;
  double p2;
  double x2;
  double y;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else if ( 0.0 < x )
  {
    y = ( x - a ) / b;
    x2 = c * y * y;
    p2 = c;

    cdf = gamma_inc ( p2, x2 );
  }

  return cdf;
}
//****************************************************************************80

bool nakagami_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    NAKAGAMI_CHECK checks the parameters of the Nakagami PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, bool NAKAGAMI_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "NAKAGAMI_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  if ( c <= 0.0 )
  {
    cout << " \n";
    cout << "NAKAGAMI_CHECK - Fatal error!\n";
    cout << "  C <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double nakagami_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    NAKAGAMI_MEAN returns the mean of the Nakagami PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B
//    0.0 < C
//
//    Output, double NAKAGAMI_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a + b * r8_gamma ( c + 0.5 ) / ( sqrt ( c ) * r8_gamma ( c ) );

  return mean;
}
//****************************************************************************80

double nakagami_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    NAKAGAMI_PDF evaluates the Nakagami PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B
//    0.0 < C.
//
//    Output, double NAKAGAMI_PDF, the value of the PDF.
//
{
  double pdf;
  double y;

  if ( x <= 0.0 )
  {
    pdf = 0.0;
  }
  else if ( 0.0 < x )
  {
    y = ( x - a ) / b;

    pdf = 2.0 * pow ( c, c ) / ( b * r8_gamma ( c ) )
    * pow ( y, ( 2.0 * c - 1.0 ) )
    * exp ( - c * y * y );

  }

  return pdf;
}
//****************************************************************************80

double nakagami_variance ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    NAKAGAMI_VARIANCE returns the variance of the Nakagami PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B
//    0.0 < C
//
//    Output, double NAKAGAMI_VARIANCE, the variance of the PDF.
//
{
  double t1;
  double t2;
  double variance;

  t1 = r8_gamma ( c + 0.5 );
  t2 = r8_gamma ( c );

  variance = b * b * ( 1.0 - t1 * t1 / ( c * t2 * t2 ) );

  return variance;
}
//****************************************************************************80

double negative_binomial_cdf ( int x, int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CDF evaluates the Negative Binomial CDF.
//
//  Discussion:
//
//    A simple summing approach is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the argument of the CDF.
//
//    Input, int A, double B, parameters of the PDF.
//    0 <= A,
//    0 < B <= 1.
//
//    Output, double NEGATIVE_BINOMIAL_CDF, the value of the CDF.
//
{
  double cdf;
  int cnk;
  double pdf;
  int y;

  cdf = 0.0;

  for ( y = a; y <= x; y++ )
  {
    cnk = binomial_coef ( y-1, a-1 );

    pdf = ( double ) ( cnk ) * pow ( b, a ) * pow ( 1.0 - b, y - a );

    cdf = cdf + pdf;
  }

  return cdf;
}
//****************************************************************************80

int negative_binomial_cdf_inv ( double cdf, int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CDF_INV inverts the Negative Binomial CDF.
//
//  Discussion:
//
//    A simple discrete approach is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, int A, double B, parameters of the PDF.
//    0 <= A,
//    0 < B <= 1.
//
//    Output, int NEGATIVE_BINOMIAL_CDF_INV, the smallest X whose cumulative
//    density function is greater than or equal to CDF.
//
{
  double cum;
  double pdf;
  int x;
  int x_max = 1000;
//
  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "NEGATIVE_BINOMIAL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  cum = 0.0;

  x = a;

  for ( ; ; )
  {
    pdf = negative_binomial_pdf ( x, a, b );

    cum = cum + pdf;

    if ( cdf <= cum || x_max <= x )
    {
      break;
    }
    x = x + 1;
  }

  return x;
}
//****************************************************************************80

void negative_binomial_cdf_values ( int *n_data, int *f, int *s, double *p,
  double *cdf )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
//
//  Discussion:
//
//    Assume that a coin has a probability P of coming up heads on
//    any one trial.  Suppose that we plan to flip the coin until we
//    achieve a total of S heads.  If we let F represent the number of
//    tails that occur in this process, then the value of F satisfies
//    a negative binomial PDF:
//
//      PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
//
//    The negative binomial CDF is the probability that there are F or
//    fewer failures upon the attainment of the S-th success.  Thus,
//
//      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = NegativeBinomialDistribution [ s, p ]
//      CDF [ dist, f ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    F C Powell,
//    Statistical Tables for Sociology, Biology and Physical Sciences,
//    Cambridge University Press, 1982.
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
//    Output, int *F, the maximum number of failures.
//
//    Output, int *S, the number of successes.
//
//    Output, double *P, the probability of a success on one trial.
//
//    Output, double *CDF, the probability of at most F failures
//    before the S-th success.
//
{
# define N_MAX 27

  double cdf_vec[N_MAX] = {
     0.6367187500000000,
     0.3632812500000000,
     0.1445312500000000,
     0.5000000000000000,
     0.2265625000000000,
     0.6250000000000000E-01,
     0.3437500000000000,
     0.1093750000000000,
     0.1562500000000000E-01,
     0.1792000000000000,
     0.4096000000000000E-01,
     0.4096000000000000E-02,
     0.7047000000000000E-01,
     0.1093500000000000E-01,
     0.7290000000000000E-03,
     0.9861587127990000,
     0.9149749500510000,
     0.7471846521450000,
     0.8499053647030009,
     0.5497160941090026,
     0.2662040052146710,
     0.6513215599000000,
     0.2639010709000000,
     0.7019082640000000E-01,
     0.1000000000000000E+01,
     0.1990000000000000E-01,
     0.1000000000000000E-03 };

  int f_vec[N_MAX] = {
     4,  3,  2,
     3,  2,  1,
     2,  1,  0,
     2,  1,  0,
     2,  1,  0,
    11, 10,  9,
    17, 16, 15,
     9,  8,  7,
     2,  1,  0 };

  double p_vec[N_MAX] = {
     0.50,
     0.50,
     0.50,
     0.50,
     0.50,
     0.50,
     0.50,
     0.50,
     0.50,
     0.40,
     0.40,
     0.40,
     0.30,
     0.30,
     0.30,
     0.30,
     0.30,
     0.30,
     0.10,
     0.10,
     0.10,
     0.10,
     0.10,
     0.10,
     0.10E-01,
     0.10E-01,
     0.10E-01 };

  int s_vec[N_MAX] = {
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    0, 1, 2 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *f = 0;
    *s = 0;
    *p = 0.0;
    *cdf = 0.0;
  }
  else
  {
    *f = f_vec[*n_data-1];
    *s = s_vec[*n_data-1];
    *p = p_vec[*n_data-1];
    *cdf = cdf_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool negative_binomial_check ( int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CHECK checks the parameters of the Negative Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, double B, parameters of the PDF.
//    0 <= A,
//    0 < B <= 1.
//
//    Output, bool NEGATIVE_BINOMIAL_CHECK, is true if the
//    parameters are legal.
//
{
  if ( a < 0 )
  {
    cout << " \n";
    cout << "NEGATIVE_BINOMIAL_CHECK - Fatal error!\n";
    cout << "  A < 0.\n";
    return false;
  }

  if ( b <= 0.0 || 1.0 < b )
  {
    cout << " \n";
    cout << "NEGATIVE_BINOMIAL_CHECK - Fatal error!\n";
    cout << "  B <= 0 or 1 < B.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double negative_binomial_mean ( int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_MEAN returns the mean of the Negative Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, double B, parameters of the PDF.
//    0 <= A,
//    0 < B <= 1.
//
//    Output, double NEGATIVE_BINOMIAL_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = ( double ) ( a ) / b;

  return mean;
}
//****************************************************************************80

double negative_binomial_pdf ( int x, int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_PDF evaluates the Negative Binomial PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = C(X-1,A-1) * B**A * ( 1 - B )**(X-A)
//
//    PDF(A,B;X) is the probability that the A-th success will
//    occur on the X-th trial, given that the probability
//    of a success on a single trial is B.
//
//    The Negative Binomial PDF is also known as the Pascal PDF or
//    the "Polya" PDF.
//
//    NEGATIVE_BINOMIAL_PDF(1,B;X) = GEOMETRIC_PDF(B;X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the number of trials.
//    A <= X.
//
//    Input, int A, the number of successes required.
//    0 <= A <= X, normally.
//
//    Input, double B, the probability of a success on a single trial.
//    0.0 < B <= 1.0.
//
//    Output, double NEGATIVE_BINOMIAL_PDF, the value of the PDF.
//
{
  int cnk;
  double pdf;

  if ( x < a )
  {
    pdf = 0.0;
  }
  else
  {
    cnk = binomial_coef ( x-1, a-1 );

    pdf = ( double ) ( cnk ) * pow ( b, a ) * pow ( 1.0 - b, x - a );
  }

  return pdf;
}
//****************************************************************************80

int negative_binomial_sample ( int a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_SAMPLE samples the Negative Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, double B, parameters of the PDF.
//    0 <= A,
//    0 < B <= 1.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int NEGATIVE_BINOMIAL_SAMPLE, a sample of the PDF.
//
{
  int num_success;
  double r;
  int x;

  if ( b == 1.0 )
  {
    x = a;
    return x;
  }
  else if ( b == 0.0 )
  {
    x = i4_huge ( );
    return x;
  }

  x = 0;
  num_success = 0;

  while ( num_success < a )
  {
    x = x + 1;
    r = r8_uniform_01 ( seed );

    if ( r <= b )
    {
      num_success = num_success + 1;
    }

  }

  return x;
}
//****************************************************************************80

double negative_binomial_variance ( int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_VARIANCE returns the variance of the Negative Binomial PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, double B, parameters of the PDF.
//    0 <= A,
//    0 < B <= 1.
//
//    Output, double NEGATIVE_BINOMIAL_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( double ) ( a ) * ( 1.0 - b ) / ( b * b );

  return variance;
}
//****************************************************************************80

double normal_01_cdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF evaluates the Normal 01 CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    A G Adams,
//    Areas Under the Normal Curve,
//    Algorithm 39,
//    Computer j.,
//    Volume 12, pages 197-198, 1969.
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Output, double CDF, the value of the CDF.
//
{
  double a1 = 0.398942280444;
  double a2 = 0.399903438504;
  double a3 = 5.75885480458;
  double a4 = 29.8213557808;
  double a5 = 2.62433121679;
  double a6 = 48.6959930692;
  double a7 = 5.92885724438;
  double b0 = 0.398942280385;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364;
  double b5 = 0.151679116635;
  double b6 = 5.29330324926;
  double b7 = 4.8385912808;
  double b8 = 15.1508972451;
  double b9 = 0.742380924027;
  double b10 = 30.789933034;
  double b11 = 3.99019417011;
  double cdf;
  double q;
  double y;
//
//  |X| <= 1.28.
//
  if ( r8_abs ( x ) <= 1.28 )
  {
    y = 0.5 * x * x;

    q = 0.5 - r8_abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5
      + a6 / ( y + a7 ) ) ) );
//
//  1.28 < |X| <= 12.7
//
  }
  else if ( r8_abs ( x ) <= 12.7 )
  {
    y = 0.5 * x * x;

    q = exp ( - y ) * b0 / ( r8_abs ( x ) - b1
      + b2 / ( r8_abs ( x ) + b3
      + b4 / ( r8_abs ( x ) - b5
      + b6 / ( r8_abs ( x ) + b7
      - b8 / ( r8_abs ( x ) + b9
      + b10 / ( r8_abs ( x ) + b11 ) ) ) ) ) );
//
//  12.7 < |X|
//
  }
  else
  {
    q = 0.0;
  }
//
//  Take account of negative X.
//
  if ( x < 0.0 )
  {
    cdf = q;
  }
  else
  {
    cdf = 1.0 - q;
  }

  return cdf;
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

void normal_01_cdf_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NormalDistribution [ 0, 1 ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
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
# define N_MAX 17

  double fx_vec[N_MAX] = {
     0.5000000000000000,
     0.5398278372770290,
     0.5792597094391030,
     0.6179114221889526,
     0.6554217416103242,
     0.6914624612740131,
     0.7257468822499270,
     0.7580363477769270,
     0.7881446014166033,
     0.8159398746532405,
     0.8413447460685429,
     0.9331927987311419,
     0.9772498680518208,
     0.9937903346742239,
     0.9986501019683699,
     0.9997673709209645,
     0.9999683287581669 };

  double x_vec[N_MAX] = {
     0.0000000000000000,
     0.1000000000000000,
     0.2000000000000000,
     0.3000000000000000,
     0.4000000000000000,
     0.5000000000000000,
     0.6000000000000000,
     0.7000000000000000,
     0.8000000000000000,
     0.9000000000000000,
     0.1000000000000000E+01,
     0.1500000000000000E+01,
     0.2000000000000000E+01,
     0.2500000000000000E+01,
     0.3000000000000000E+01,
     0.3500000000000000E+01,
     0.4000000000000000E+01 };

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

double normal_01_mean ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_MEAN returns the mean of the Normal 01 PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = 0.0;

  return mean;
}
//****************************************************************************80

double normal_01_pdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_PDF evaluates the Normal 01 PDF.
//
//  Discussion:
//
//    The Normal 01 PDF is also called the "Standard Normal" PDF, or
//    the Normal PDF with 0 mean and variance 1.
//
//    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  pdf = exp ( -0.5 * x * x ) / sqrt ( 2.0 * pi );

  return pdf;
}
//**********************************************************************

double normal_01_sample ( int *seed )

//**********************************************************************
//
//  Purpose:
//
//    NORMAL_01_SAMPLE samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double NORMAL_01_SAMPLE, a normally distributed random value.
//
{
  const double pi = 3.14159265358979323;
  double r1;
  double r2;
  static int used = -1;
  double x;
  static double y = 0.0;

  if ( used == -1 )
  {
    used = 0;
  }
//
//  If we've used an even number of values so far, generate two more, return one,
//  and save one.
//
  if ( ( used % 2 )== 0 )
  {
    for ( ; ; )
    {
      r1 = r8_uniform_01 ( seed );
      if ( r1 != 0.0 )
      {
        break;
      }
    }

    r2 = r8_uniform_01 ( seed );

    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );
    y = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * pi * r2 );
  }
  else
  {

    x = y;

  }

  used = used + 1;

  return x;
}
//****************************************************************************80

double normal_01_variance ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_VARIANCE returns the variance of the Normal 01 PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = 1.0;

  return variance;
}
//****************************************************************************80

double *normal_01_vector ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_VECTOR samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
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
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X(N), a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, real R(N+1), is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int SAVED, is 0 or 1 depending on whether there is a
//    single saved value left over from the previous call.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.  This starts off as 1:N, but is adjusted
//    if we have a saved value that can be immediately stored in X(1),
//    and so on.
//
//    Local, real Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
  int i;
  int m;
  static int made = 0;
  const double pi = 3.14159265358979323;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;

  x = new double[n];
//
//  I'd like to allow the user to reset the internal data.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
//
//  Maybe we don't need any more values.
//
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
//
//  If we need just one new value, do that here to avoid null arrays.
//
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01 ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );
    y =         sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * pi * r[1] );

    saved = 1;

    made = made + 2;

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01 ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01 ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
    y           = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
}
//****************************************************************************80

double normal_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CDF evaluates the Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;
  double y;

  y = ( x - a ) / b;

  cdf = normal_01_cdf ( y );

  return cdf;
}
//****************************************************************************80

double normal_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CDF_INV inverts the Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Algorithm AS 111,
//    Applied Statistics,
//    Volume 26, pages 118-121, 1977.
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double X, the corresponding argument.
//
{
  double x;
  double x2;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "NORMAL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x2 = normal_01_cdf_inv ( cdf );

  x = a + b * x2;

  return x;
}
//****************************************************************************80

void normal_cdf_values ( int *n_data, double *mu, double *sigma, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CDF_VALUES returns some values of the Normal CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NormalDistribution [ mu, sigma ]
//      CDF [ dist, x ]
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
//    Output, double *MU, the mean of the distribution.
//
//    Output, double *SIGMA, the variance of the distribution.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double fx_vec[N_MAX] = {
     0.5000000000000000,
     0.9772498680518208,
     0.9999683287581669,
     0.9999999990134124,
     0.6914624612740131,
     0.6305586598182364,
     0.5987063256829237,
     0.5792597094391030,
     0.6914624612740131,
     0.5000000000000000,
     0.3085375387259869,
     0.1586552539314571 };

  double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  double sigma_vec[N_MAX] = {
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *mu = 0.0;
    *sigma = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *mu = mu_vec[*n_data-1];
    *sigma = sigma_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool normal_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CHECK checks the parameters of the Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, bool NORMAL_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "NORMAL_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double normal_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MEAN returns the mean of the Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double normal_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_PDF evaluates the Normal PDF.
//
//  Discussion:
//
//    PDF(A,B;X)
//      = exp ( - 0.5 * ( ( X - A ) / B )^2 )
//      / ( B * SQRT ( 2 * PI ) )
//
//    The normal PDF is also known as the Gaussian PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;
  double y;

  y = ( x - a ) / b;

  pdf = exp ( -0.5 * y * y )  / ( b * sqrt ( 2.0 * pi ) );

  return pdf;
}
//****************************************************************************80

double normal_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_SAMPLE samples the Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double NORMAL_SAMPLE, a sample of the PDF.
//
{
  double x;

  x = normal_01_sample ( seed );

  x = a + b * x;

  return x;
}
//****************************************************************************80

double normal_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_VARIANCE returns the variance of the Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = b * b;

  return variance;
}
//****************************************************************************80

double *normal_vector ( int n, double mean, double dev, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_VECTOR samples the normal probability distribution.
//
//  Discussion:
//
//    The normal probability distribution function (PDF) has
//    a user-specified mean and standard deviation.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
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
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input, double MEAN, the desired mean value.
//
//    Input, double DEV, the desired standard deviation.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double NORMAL_VECTOR[N], a sample of the normal PDF.
//
{
  int i;
  double *x;

  x = normal_01_vector ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    x[i] = mean + dev * x[i];
  }

  return x;
}
//****************************************************************************80

void owen_values ( int *n_data, double *h, double *a, double *t )

//****************************************************************************80
//
//  Purpose:
//
//    OWEN_VALUES returns some values of Owen's T function.
//
//  Discussion:
//
//    Owen's T function is useful for computation of the bivariate normal
//    distribution and the distribution of a skewed normal distribution.
//
//    Although it was originally formulated in terms of the bivariate
//    normal function, the function can be defined more directly as
//
//      T(H,A) = 1 / ( 2 * pi ) *
//        Integral ( 0 <= X <= A ) e^(H^2*(1+X^2)/2) / (1+X^2) dX
//
//    In Mathematica, the function can be evaluated by:
//
//      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2004
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
//    Output, double *H, a parameter.
//
//    Output, double *A, the upper limit of the integral.
//
//    Output, double *T, the value of the function.
//
{
# define N_MAX 22

  double a_vec[N_MAX] = {
    0.5000000000000000,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.1000000000000000E+02,
    0.1000000000000000E+03 };

  double h_vec[N_MAX] = {
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.5000000000000000,
    0.5000000000000000,
    0.5000000000000000,
    0.5000000000000000,
    0.2500000000000000,
    0.2500000000000000,
    0.2500000000000000,
    0.2500000000000000,
    0.1250000000000000,
    0.1250000000000000,
    0.1250000000000000,
    0.1250000000000000,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02 };

  double t_vec[N_MAX] = {
    0.4306469112078537E-01,
    0.6674188216570097E-01,
    0.7846818699308410E-01,
    0.7929950474887259E-01,
    0.6448860284750376E-01,
    0.1066710629614485,
    0.1415806036539784,
    0.1510840430760184,
    0.7134663382271778E-01,
    0.1201285306350883,
    0.1666128410939293,
    0.1847501847929859,
    0.7317273327500385E-01,
    0.1237630544953746,
    0.1737438887583106,
    0.1951190307092811,
    0.7378938035365546E-01,
    0.1249951430754052,
    0.1761984774738108,
    0.1987772386442824,
    0.2340886964802671,
    0.2479460829231492 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *h = 0.0;
    *a = 0.0;
    *t = 0.0;
  }
  else
  {
    *h = h_vec[*n_data-1];
    *a = a_vec[*n_data-1];
    *t = t_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double pareto_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_CDF evaluates the Pareto CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double PARETO_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x < a )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 1.0 - pow ( ( a / x ), b );
  }

  return cdf;
}
//****************************************************************************80

double pareto_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_CDF_INV inverts the Pareto CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double PARETO_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "PARETO_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a / pow ( 1.0 - cdf, 1.0 / b );

  return x;
}
//****************************************************************************80

bool pareto_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_CHECK checks the parameters of the Pareto CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, bool PARETO_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << " \n";
    cout << "PARETO_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "PARETO_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double pareto_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_MEAN returns the mean of the Pareto PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double PARETO_MEAN, the mean of the PDF.
//
{
  double mean;

  if ( b <= 1.0 )
  {
    cout << " \n";
    cout << "PARETO_MEAN - Fatal error!\n";
    cout << "  For B <= 1, the mean does not exist.\n";
    mean = 0.0;
    return mean;
  }

  mean = b * a / ( b - 1.0 );

  return mean;
}
//****************************************************************************80

double pareto_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_PDF evaluates the Pareto PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = B * A**B / X**(B+1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A.
//    0.0 < B.
//
//    Output, double PARETO_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < a )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = b * pow ( a, b ) / pow ( x, ( b + 1.0 ) );
  }

  return pdf;
}
//****************************************************************************80

double pareto_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_SAMPLE samples the Pareto PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double PARETO_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = pareto_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double pareto_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_VARIANCE returns the variance of the Pareto PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double PARETO_VARIANCE, the variance of the PDF.
//
{
  double variance;

  if ( b <= 2.0 )
  {
    cout << " \n";
    cout << "PARETO_VARIANCE - Warning!\n";
    cout << "  For B <= 2, the variance does not exist.\n";
    variance = 0.0;
    return variance;
  }

  variance = a * a * b / ( pow ( ( b - 1.0 ), 2 ) * ( b - 2.0 ) );

  return variance;
}
//****************************************************************************80

bool pearson_05_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_CHECK checks the parameters of the Pearson 5 PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, bool PEARSON_05_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << " \n";
    cout << "PEARSON_05_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "PEARSON_05_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double pearson_05_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_MEAN evaluates the mean of the Pearson 5 PDF.
//
//  Discussion:
//
//    The mean is undefined for B <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, double PEARSON_05_MEAN, the mean of the PDF.
//
{
  double mean;

  if ( b <= 1.0 )
  {
    cout << " \n";
    cout << "PEARSON_05_MEAN - Warning!\n";
    cout << "  MEAN undefined for B <= 1.\n";
    mean = c;
    return mean;
  }

  mean = c + a / ( b - 1.0 );

  return mean;
}
//****************************************************************************80

double pearson_05_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_PDF evaluates the Pearson 5 PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = A**B * ( X - C )**(-B-1)
//      * exp ( - A / ( X - C ) ) / Gamma ( B )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    C < X
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;

  if ( x <= c )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = pow ( a, b ) * pow ( x - c, - b - 1.0 )
      * exp ( - a / ( x - c ) ) / r8_gamma ( b );
  }

  return pdf;
}
//****************************************************************************80

double pearson_05_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_SAMPLE samples the Pearson 5 PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double PEARSON_05_SAMPLE, a sample of the PDF.
//
{
  double a2;
  double b2;
  double c2;
  double x;
  double x2;

  a2 = 0.0;
  b2 = b;
  c2 = 1.0 / a;

  x2 = gamma_sample ( a2, b2, c2, seed );

  x = c + 1.0 / x2;

  return x;
}
//****************************************************************************80

bool planck_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_CHECK checks the parameters of the Planck PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, bool PLANCK_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << "\n";
    cout << "PLANCK_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "PLANCK_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double planck_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_MEAN returns the mean of the Planck PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double PLANCK_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = ( b + 1.0 ) * zeta ( b + 2.0 ) / zeta ( b + 1.0 );

  return mean;
}
//****************************************************************************80

double planck_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_PDF evaluates the Planck PDF.
//
//  Discussion:
//
//    The Planck PDF has the form
//
//      PDF(A,B;X) = A**(B+1) * X**B / ( exp ( A * X ) - 1 ) / K
//
//    where K is the normalization constant, and has the value
//
//      K = Gamma ( B + 1 ) * Zeta ( B + 1 ).
//
//    The original Planck distribution governed the frequencies in
//    blackbody radiation at a given temperature T, and has the form
//
//      PDF(A;X) = K * X**3 / ( exp ( A * X ) - 1 )
//
//    where
//
//      K = 15 / PI**4.
//
//    Thus, in terms of the Planck PDF, the original Planck distribution
//    has A = 1, B = 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double PLANCK_PDF, the value of the PDF.
//
{
  double k;
  double pdf;

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    k = r8_gamma ( b + 1.0 ) * zeta ( b + 1.0 );
    pdf = pow ( a, b + 1.0 ) * pow ( x, b ) / ( exp ( a * x ) - 1.0 ) / k;
  }

  return pdf;
}
//****************************************************************************80

double planck_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_SAMPLE samples the Planck PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer Verlag, 1986, pages 552.
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double PLANCK_SAMPLE, a sample of the PDF.
//
{
  double a2;
  double b2;
  double c2;
  double g;
  double x;
  int z;
//
  a2 = 0.0;
  b2 = 1.0;
  c2 = b + 1.0;

  g = gamma_sample ( a2, b2, c2, seed );

  z = zipf_sample ( c2, seed );

  x = g / ( a * ( double ) ( z ) );

  return x;
}
//****************************************************************************80

double planck_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_VARIANCE returns the variance of the Planck PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double PLANCK_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = 0.0;

  return variance;
}
//****************************************************************************80

double point_distance_1d_pdf ( double x, int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_DISTANCE_1D_PDF evaluates the point distance PDF in 1D.
//
//  Discussion:
//
//    It is assumed that a set of points has been generated in 1D
//    according to a Poisson process.  The number of points in a region
//    of size LENGTH is a Poisson variate with mean value B * LENGTH.
//
//    For a point chosen at random, we may now find the nearest
//    Poisson point, the second nearest and so on.  We are interested
//    in the PDF that governs the expected behavior of the distances
//    of rank A = 1, 2, 3, ... with Poisson density B.
//
//    Note that this PDF is a form of the Gamma PDF.???
//
//    PDF(A,B;X) = B**A * X**( A - 1 ) * EXP ( - B * X ) / ( A - 1 )!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X.
//
//    Input, int A, indicates the degree of nearness of the point.
//    A = 1 means the nearest point, A = 2 the second nearest, and so on.
//    0 < A.
//
//    Input, double B, the point density.  0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;

  if ( a < 1 )
  {
    cout << " \n";
    cout << "POINT_DISTANCE_1D_PDF - Fatal error!\n";
    cout << "  Input parameter A < 1.\n";
    exit ( 1 );
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "POINT_DISTANCE_1D_PDF - Fatal error!\n";
    cout << "  Input parameter B <= 0.0.\n";
    exit ( 1 );
  }

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = pow ( b, a ) * pow ( x, a - 1 ) * exp ( - b * x )
      / i4_factorial ( a - 1 );
  }

  return pdf;
}
//****************************************************************************80

double point_distance_2d_pdf ( double x, int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_DISTANCE_2D_PDF evaluates the point distance PDF in 2D.
//
//  Discussion:
//
//    It is assumed that a set of points has been generated in 2D
//    according to a Poisson process.  The number of points in a region
//    of size AREA is a Poisson variate with mean value B * AREA.
//
//    For a point chosen at random, we may now find the nearest
//    Poisson point, the second nearest and so on.  We are interested
//    in the PDF that governs the expected behavior of the distances
//    of rank A = 1, 2, 3, ... with Poisson density B.
//
//    PDF(A,B;X) = 2 * ( B * PI )**A * X**( 2 * A - 1 )
//      * EXP ( - B * PI * X * X ) / ( A - 1 )!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
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
//    CRC Press, 1996, pages 579.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X.
//
//    Input, int A, indicates the degree of nearness of the point.
//    A = 1 means the nearest point, A = 2 the second nearest, and so on.
//    0 < A.
//
//    Input, double B, the point density.  0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  if ( a < 1 )
  {
    cout << " \n";
    cout << "POINT_DISTANCE_2D_PDF - Fatal error!\n";
    cout << "  Input parameter A < 1.\n";
    exit ( 1 );
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "POINT_DISTANCE_2D_PDF - Fatal error!\n";
    cout << "  Input parameter B <= 0.0.\n";
    exit ( 1 );
  }

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = 2.0 * pow ( b * pi, a ) * pow ( x, 2 * a - 1 )
      * exp ( - b * pi * x * x ) / i4_factorial ( a - 1 );
  }

  return pdf;
}
//****************************************************************************80

double point_distance_3d_pdf ( double x, int a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_DISTANCE_3D_PDF evaluates the point distance PDF in the 3D.
//
//  Discussion:
//
//    It is assumed that a set of points has been generated in 3D
//    according to a Poisson process.  The number of points in a region
//    of size VOLUME is a Poisson variate with mean value B * VOLUME.
//
//    For a point chosen at random, we may now find the nearest
//    Poisson point, the second nearest and so on.  We are interested
//    in the PDF that governs the expected behavior of the distances
//    of rank A = 1, 2, 3, ... with Poisson density B.
//
//    PDF(A,B;X) = 3 * ( (4/3) * B * PI )**A * X**( 3 * A - 1 )
//      * EXP ( - (4/3) * B * PI * X * X * X ) / ( A - 1 )!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
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
//    CRC Press, 1996, pages 580.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X.
//
//    Input, int A, indicates the degree of nearness of the point.
//    A = 1 means the nearest point, A = 2 the second nearest, and so on.
//    0 < A.
//
//    Input, double B, the Poisson point density.  0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  if ( a < 1 )
  {
    cout << " \n";
    cout << "POINT_DISTANCE_3D_PDF - Fatal error!\n";
    cout << "  Input parameter A < 1.\n";
    exit ( 1 );
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "POINT_DISTANCE_3D_PDF - Fatal error!\n";
    cout << "  Input parameter B <= 0.0.\n";
    exit ( 1 );
  }

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = 3.0 * pow ( ( ( 4.0 / 3.0 ) * b * pi ), a )
      * pow ( x, 3 * a - 1 ) * exp ( - ( 4.0 / 3.0 ) * b * pi * x * x * x )
      / i4_factorial ( a - 1 );
  }

  return pdf;
}
//****************************************************************************80

double poisson_cdf ( int k, double a )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_CDF evaluates the Poisson CDF.
//
//  Discussion:
//
//    CDF(K,A) is the probability that the number of events observed
//    in a unit time period will be no greater than K, given that the
//    expected number of events in a unit time period is A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int K, the argument of the CDF.
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double POISSON_CDF, the value of the CDF.
//
{
  double cdf;
  int i;
  double last;
  double next;
//
//  Check.
//
  if ( a <= 0.0 )
  {
    cout << "\n";
    cout << "POISSON_CDF - Fatal error!\n";
    cout << "  A <= 0.\n";
    return 0;
  }
//
//  Special cases.
//
  if ( k < 0 )
  {
    cdf = 0.0;
    return cdf;
  }
//
//  General case.
//
  next = exp ( - a );
  cdf = next;

  for ( i = 1; i <= k; i++ )
  {
    last = next;
    next = last * a / ( double ) i;
    cdf = cdf + next;
  }

  return cdf;
}
//****************************************************************************80

int poisson_cdf_inv ( double cdf, double a )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_CDF_INV inverts the Poisson CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, a value of the CDF.
//    0 <= CDF < 1.
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A.
//
//    Output, int POISSON_CDF_INV, the corresponding argument.
//
{
  int i;
  double last;
  double newval;
  double sum2;
  double sumold;
  int x;
  int xmax = 100;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "POISSON_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }
//
//  Now simply start at X = 0, and find the first value for which
//  CDF(X-1) <= CDF <= CDF(X).
//
  sum2 = 0.0;

  for ( i = 0; i <= xmax; i++ )
  {
    sumold = sum2;

    if ( i == 0 )
    {
      newval = exp ( - a );
      sum2 = newval;
    }
    else
    {
      last = newval;
      newval = last * a / ( double ) ( i );
      sum2 = sum2 + newval;
    }

    if ( sumold <= cdf && cdf <= sum2 )
    {
      x = i;
      return x;
    }
  }

  cout << " \n";
  cout << "POISSON_CDF_INV - Warning!\n";
  cout << "  Exceeded XMAX = " << xmax << "\n";

  x = xmax;

  return x;
}
//****************************************************************************80

void poisson_cdf_values ( int *n_data, double *a, int *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_CDF_VALUES returns some values of the Poisson CDF.
//
//  Discussion:
//
//    CDF(X)(A) is the probability of at most X successes in unit time,
//    given that the expected mean number of successes is A.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = PoissonDistribution [ a ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
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
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 653-658.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, the parameter of the function.
//
//    Output, int *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 21

  double a_vec[N_MAX] = {
     0.02,
     0.10,
     0.10,
     0.50,
     0.50,
     0.50,
     1.00,
     1.00,
     1.00,
     1.00,
     2.00,
     2.00,
     2.00,
     2.00,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00 };

  double fx_vec[N_MAX] = {
     0.9801986733067553,
     0.9048374180359596,
     0.9953211598395555,
     0.6065306597126334,
     0.9097959895689501,
     0.9856123220330293,
     0.3678794411714423,
     0.7357588823428846,
     0.9196986029286058,
     0.9810118431238462,
     0.1353352832366127,
     0.4060058497098381,
     0.6766764161830635,
     0.8571234604985470,
     0.6737946999085467E-02,
     0.4042768199451280E-01,
     0.1246520194830811,
     0.2650259152973617,
     0.4404932850652124,
     0.6159606548330631,
     0.7621834629729387 };

  int x_vec[N_MAX] = {
     0, 0, 1, 0,
     1, 2, 0, 1,
     2, 3, 0, 1,
     2, 3, 0, 1,
     2, 3, 4, 5,
     6 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *a = 0.0;
    *x = 0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool poisson_check ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_CHECK checks the parameter of the Poisson PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A.
//
//    Output, bool POISSON_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << " \n";
    cout << "POISSON_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double poisson_kernel ( double r, int n, double c[], double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_KERNEL evaluates the Poisson kernel.
//
//  Discussion:
//
//    P(X,Y) = ( R^2 - |X-C|^2 ) / ( R * A * |X-Y|^N )
//
//    where the N-dimensional ball has radius R and center C,
//    and A is the area of the unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the ball.
//
//    Input, int N, the spatial dimension.
//
//    Input, double C(N), the center of the ball.
//
//    Input, double X(N), a point inside the ball.
//
//    Input, double Y(N), a point on the surface of the ball.
//
//    Output, double POISSON_KERNEL, the Poisson kernel function P(X,Y).
//
{
  double area;
  double b;
  double p;
  double t;
  double xc_diff_norm;
  double xy_diff_norm;

  xc_diff_norm = r8vec_diff_norm ( n, x, c );
  xy_diff_norm = r8vec_diff_norm ( n, x, y );
  area = sphere_unit_area_nd ( n );

  t = ( r + xc_diff_norm ) * ( r - xc_diff_norm );
  b = r * area * pow ( xy_diff_norm, n );
  p = t / b;

  return p;
}
//****************************************************************************80

double poisson_mean ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_MEAN returns the mean of the Poisson PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double POISSON_MEAN, the mean of the PDF.
//
{
  if ( a <= 0.0 )
  {
    cout << "\n";
    cout << "POISSON_MEAN - Fatal error!\n" ;
    cout << "  A <= 0.\n";
    return 0;
  }

  return a;

}
//****************************************************************************80

double poisson_pdf ( int k, double a )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_PDF evaluates the Poisson PDF.
//
//  Discussion:
//
//    PDF(K,A) is the probability that the number of events observed
//    in a unit time period will be K, given the expected number
//    of events in a unit time.
//
//    The parameter A is the expected number of events per unit time.
//
//    The Poisson PDF is a discrete version of the exponential PDF.
//
//    The time interval between two Poisson events is a random
//    variable with the exponential PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int K, the argument of the PDF.
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double POISSON_PDF, the value of the PDF.
//
{
  double pdf;
//
//  Check.
//
  if ( a <= 0.0 )
  {
    cout << "\n";
    cout << "POISSON_PDF - Fatal error!\n";
    cout << "  A <= 0.\n";
    return (0.0);
  }

  pdf = exp ( -a ) * pow ( a, ( double ) k ) / i4_factorial ( k );

  return pdf;
}
//****************************************************************************80

int poisson_sample ( double a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_SAMPLE samples the Poisson PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Input/output, int *SEED, the random number generator seed.
//
//    Output, int POISSON_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  int i;
  int KMAX = 100;
  double last;
  double next;
  double sum;
  double sumold;
//
//  Check.
//
  if ( a <= 0.0 )
  {
    cout << "\n";
    cout << "POISSON_SAMPLE - Fatal error!\n";
    cout << "  A <= 0.\n";
    return 0;
  }
//
//  Pick a random value of CDF.
//
  cdf = uniform_01_sample ( seed );
//
//  Now simply start at K = 0, and find the first value for which
//  CDF(K-1) <= CDF <= CDF(K).
//
  sum = 0.0;

  for ( i = 0; i <= KMAX; i++ )
  {
    sumold = sum;

    if ( i == 0 )
    {
      next = exp ( -a );
      sum = next;
    }
    else
    {
      last = next;
      next = last * a / ( double ) i;
      sum = sum + next;
    }

    if ( sumold <= cdf && cdf <= sum )
    {
      return i;
    }

  }

  cout << "\n";
  cout << "POISSON_SAMPLE - Warning!\n";
  cout << "  Exceeded KMAX = " << KMAX << "\n";

  return KMAX;
}
//****************************************************************************80

double poisson_variance ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_VARIANCE returns the variance of the Poisson PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double POISSON_VARIANCE, the variance of the PDF.
//
{
//
//  Check.
//
  if ( a <= 0.0 )
  {
    cout << "\n";
    cout << "POISSON_VARIANCE - Fatal error!\n";
    cout << "  A <= 0.\n";
    return (0.0);
  }

  return a;
}
//****************************************************************************80

double power_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_CDF evaluates the Power CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A, 0.0 < B,
//
//    Output, double POWER_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else if ( x <= b )
  {
    cdf = pow ( ( x / b ), a );
  }
  else
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double power_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_CDF_INV inverts the Power CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, double POWER_CDF_INV, the argument of the CDF.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "POWER_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = 0.0;
  }
  else if ( cdf < 1.0 )
  {
    x = b * exp ( log ( cdf ) / a );
  }
  else
  {
    x = b;
  }

  return x;
}
//****************************************************************************80

bool power_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_CHECK checks the parameter of the Power PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, bool POWER_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << " \n";
    cout << "POWER_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "POWER_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double power_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_MEAN returns the mean of the Power PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, double POWER_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a * b / ( a + 1.0 );

  return mean;
}
//****************************************************************************80

double power_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_PDF evaluates the Power PDF.
//
//  Discussion:
//
//    PDF(A;X) = (A/B) * (X/B)**(A-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Daniel Zwillinger and Stephen Kokoska,
//    CRC Standard Probability and Statistics Tables and Formulae,
//    Chapman and Hall/CRC, 2000, pages 152-153.
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X <= B.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, double POWER_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 0.0 || b < x )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = ( a / b ) * pow ( x / b, a - 1.0 );
  }

  return pdf;
}
//****************************************************************************80

double power_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_SAMPLE samples the Power PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double POWER_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = power_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double power_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_VARIANCE returns the variance of the Power PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, double POWER_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = b * b * a / ( ( a + 1.0 ) * ( a + 1.0 ) * ( a + 2.0 ) );

  return variance;
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
     -0.5772156649015329,
     -0.4237549404110768,
     -0.2890398965921883,
     -0.1691908888667997,
     -0.6138454458511615E-01,
      0.3648997397857652E-01,
      0.1260474527734763,
      0.2085478748734940,
      0.2849914332938615,
      0.3561841611640597,
      0.4227843350984671 };

  double x_vec[N_MAX] = {
     1.0,
     1.1,
     1.2,
     1.3,
     1.4,
     1.5,
     1.6,
     1.7,
     1.8,
     1.9,
     2.0 };

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

double quasigeometric_cdf ( int x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_CDF evaluates the Quasigeometric CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the maximum number of trials.
//
//    Input, double A, the probability of 0 successes.
//    0.0 <= A <= 1.0.
//
//    Input, double B, the depreciation constant.
//    0.0 <= B < 1.0.
//
//    Output, double QUASIGEOMETRIC_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x < 0 )
  {
    cdf = 0.0;
  }
  else if ( x == 0 )
  {
    cdf = a;
  }
  else if ( b == 0.0 )
  {
    cdf = 1.0;
  }
  else
  {
    cdf = a + ( 1.0 - a ) * ( 1.0 - pow ( b, x ) );
  }

  return cdf;
}
//****************************************************************************80

int quasigeometric_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_CDF_INV inverts the Quasigeometric CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0
//
//    Input, double A, the probability of 0 successes.
//    0.0 <= A <= 1.0.
//
//    Input, double B, the depreciation constant.
//    0.0 <= B < 1.0.
//
//    Output, int QUASIGEOMETRIC_CDF_INV, the corresponding value of X.
//
{
  int x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cerr << "\n";
    cerr << "QUASIGEOMETRIC_CDF_INV - Fatal error!\n";
    cerr << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf < a )
  {
    x = 0;
  }
  else if ( b == 0.0 )
  {
    x = 1;
  }
  else
  {
    x = 1 + ( int ) ( ( log ( 1.0 - cdf ) - log ( 1.0 - a ) ) / log ( b ) );
  }

  return x;
}
//****************************************************************************80

bool quasigeometric_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_CHECK checks the parameters of the Quasigeometric CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of 0 successes.
//    0.0 <= A <= 1.0.
//
//    Input, double B, the depreciation constant.
//    0.0 <= B < 1.0.
//
//    Output, bool QUASIGEOMETRIC_CHECK, is true if the parameters are legal.
//
{
  bool check;

  if ( a < 0.0 || 1.0 < a )
  {
    cerr << "\n";
    cerr << "QUASIGEOMETRIC_CHECK - Fatal error!\n";
    cerr << "  A < 0 or 1 < A.\n";
    check = false;
    exit ( 1 );
  }

  if ( b < 0.0 || 1.0 <= b )
  {
    cerr << "\n";
    cerr << "QUASIGEOMETRIC_CHECK - Fatal error!\n";
    cerr << "  B < 0 or 1 <= B.\n";
    check = false;
    exit ( 1 );
  }

  check = true;

  return check;
}
//****************************************************************************80

double quasigeometric_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_MEAN returns the mean of the Quasigeometric PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of 0 successes.
//    0.0 <= A <= 1.0.
//
//    Input, double B, the depreciation constant.
//    0.0 <= B < 1.0.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = ( 1.0 - a  ) / ( 1.0 - b );

  return mean;
}
//****************************************************************************80

double quasigeometric_pdf ( int x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_PDF evaluates the Quasigeometric PDF.
//
//  Discussion:
//
//    PDF(A,B;X) =    A                     if 0  = X;
//               = (1-A) * (1-B) * B^(X-1)  if 1 <= X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Darren Glass, Philip Lowry,
//    Quasiquasigeometric Distributions and Extra Inning Baseball Games,
//    Mathematics Magazine,
//    Volume 81, Number 2, April 2008, pages 127-137.
//
//    Paul Nahin,
//    Digital Dice: Computational Solutions to Practical Probability Problems,
//    Princeton University Press, 2008,
//    ISBN13: 978-0-691-12698-2,
//    LC: QA273.25.N34.
//
//  Parameters:
//
//    Input, int X, the independent variable.
//    0 <= X
//
//    Input, double A, the probability of 0 successes.
//    0.0 <= A <= 1.0.
//
//    Input, double B, the depreciation constant.
//    0.0 <= B < 1.0.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 0 )
  {
    pdf = 0.0;
  }
  else if ( x == 0 )
  {
    pdf = a;
  }
  else if ( b == 0.0 )
  {
    if ( x == 1 )
    {
      pdf = 1.0;
    }
    else
    {
      pdf = 0.0;
    }
  }
  else
  {
    pdf = ( 1.0 - a ) * ( 1.0 - b ) * pow ( b, x - 1 );
  }
  return pdf;
}
//****************************************************************************80

int quasigeometric_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_SAMPLE samples the Quasigeometric PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of 0 successes.
//    0.0 <= A <= 1.0.
//
//    Input, double B, the depreciation constant.
//    0.0 <= B < 1.0.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, int QUASIGEOMETRIC_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  int x;

  cdf = r8_uniform_01 ( seed );

  x = quasigeometric_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double quasigeometric_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_VARIANCE returns the variance of the Quasigeometric PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the probability of 0 successes.
//    0.0 <= A <= 1.0.
//
//    Input, double B, the depreciation constant.
//    0.0 <= B < 1.0.
//
//    Output, double QUASIGEOMETRIC_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( 1.0 - a ) * ( a + b ) / ( 1.0 - b ) / ( 1.0 - b );

  return variance;
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

float r4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_01 returns a real pseudorandom number.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r4_uniform_01 = seed / ( 2**31 - 1 )
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
    value = -x;
  }
  return value;
}
//****************************************************************************80

int r8_ceiling ( double r )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CEILING rounds an R8 "up" to the nearest integer.
//
//  Example:
//
//    R     R8_CEILING
//
//   -1.1  -1
//   -1.0  -1
//   -0.9   0
//    0.0   0
//    5.0   5
//    5.1   6
//    5.9   6
//    6.0   6
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
//    Input, double R, the real value to be rounded up.
//
//    Output, int R8_CEILING, the rounded value.
//
{
  int value;

  value = int ( r );

  if ( ( double ) value < r )
  {
    value = value + 1;
  }

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
//    05 March 2012
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

  value = sin ( theta );

  if ( value == 0.0 )
  {
    cout << " \n";
    cout << "R8_CSC - Fatal error!\n";
    cout << "  Cosecant undefined for THETA = " << theta << "\n";
    exit ( 1 );
  }

  value = 1.0 / value;

  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 round off unit.
//
//  Discussion:
//
//    R8_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
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
//    Output, double R8_EPSILON, the floating point round-off unit.
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

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
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
//
//  Coefficients for minimax approximation over (12, INF).
//
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
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
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
  double twelve = 12.0;
  double two = 2.0;
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
  double zero = 0.0;;

  parity = false;
  fact = one;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + one;
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
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
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
        y = y + one;
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
      sum = sum + ( y - half ) * log ( y );
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

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

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
//    HUGE_VAL is the largest representable legal real number, and is usually
//    defined in math.h, or sometimes in stdlib.h.
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
  return HUGE_VAL;
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

double r8_modp ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MODP returns the nonnegative remainder of R8 division.
//
//  Discussion:
//
//    If
//      REM = R8_MODP ( X, Y )
//      RMULT = ( X - REM ) / Y
//    then
//      X = Y * RMULT + REM
//    where REM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360.0) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
//
//  Example:
//
//        I         J     MOD R8_MODP  R8_MODP Factorization
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
    cout << "\n";
    cout << "R8_MODP - Fatal error!\n";
    cout << "  R8_MODP ( X, Y ) called with Y = " << y << "\n";
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

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         R8_NINT
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
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int s;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }

  return ( s * ( int ) ( r8_abs ( x ) + 0.5 ) );
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
//    06 May 2003
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
  const double pi = 3.14159265358979323;

  return pi;
}
//****************************************************************************80

double r8_random ( double rlo, double rhi, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_RANDOM returns a scaled pseudorandom R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RLO, RHI, the minimum and maximum values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8_RANDOM, the randomly chosen value.
//
{
  double t;

  t = r8_uniform_01 ( seed );

  return ( 1.0 - t ) * rlo + t * rhi;
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
  if ( x < 0.0 )
  {
    return ( -1.0 );
  }
  else
  {
    return ( 1.0 );
  }
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
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
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
//    value should not be 0.  On output, SEED has been updated.
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
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT, with an optional title.
//
//  Discussion:
//
//    The doubly dimensioned array A is treated as a one dimensional vector,
//    stored by COLUMNS.  Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2003
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
//    Input, char *TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Modified:
//
//    09 April 2004
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
//    Input, char *TITLE, a title.
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
    cout << "  ---\n";
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
      cout << setw(5) << i << "  ";
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

double *r8row_max ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_MAX returns the maximums of an R8ROW.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    MAX =
//      3
//      7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, double A[M*N], the array to be examined.
//
//    Output, double R8ROW_MAX[M], the maximums of the rows.
//
{
  int i;
  int j;
  double *amax;

  amax = new double[m];

  for ( i = 0; i < m; i++ )
  {
    amax[i] = a[i+0*m];

    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[i+j*m] )
      {
        amax[i] = a[i+j*m];
      }
    }
  }

  return amax;
}
//****************************************************************************80

double *r8row_mean ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_MEAN returns the means of an R8ROW.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    MEAN =
//      2
//      5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, double A[M*N], the array to be examined.
//
//    Output, double R8ROW_MEAN[M], the means, or averages, of the rows.
//
{
  int i;
  int j;
  double *mean;

  mean = new double[m];

  for ( i = 0; i < m; i++ )
  {
    mean[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      mean[i] = mean[i] + a[i+j*m];
    }
    mean[i] = mean[i] / ( double ) ( n );
  }

  return mean;
}
//****************************************************************************80

double *r8row_min ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_MIN returns the minimums of an R8ROW.
//
//  Example:
//
//    A =
//      1  2  3
//      2  6  7
//
//    MIN =
//      1
//      2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, double A[M*N], the array to be examined.
//
//    Output, double R8ROW_MIN[M], the minimums of the rows.
//
{
  int i;
  int j;
  double *amin;

  amin = new double[m];

  for ( i = 0; i < m; i++ )
  {
    amin[i] = a[i+0*m];
    for ( j = 1; j < n; j++ )
    {
      if ( a[i+j*m] < amin[i] )
      {
        amin[i] = a[i+j*m];
      }
    }
  }

  return amin;
}
//****************************************************************************80

double *r8row_variance ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_VARIANCE returns the variances of an R8ROW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, double A[M*N], the array whose variances are desired.
//
//    Output, double R8ROW_VARIANCE[M], the variances of the rows.
//
{
  int i;
  int j;
  double mean;
  double *variance;

  variance = new double[m];

  for ( i = 0; i < m; i++ )
  {
    mean = 0.0;
    for ( j = 0; j < n; j++ )
    {
      mean = mean + a[i+j*m];
    }
    mean = mean / ( double ) ( n );

    variance[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      variance[i] = variance[i] + pow ( ( a[i+j*m] - mean ), 2 );
    }

    if ( 1 < n )
    {
      variance[i] = variance[i] / ( double ) ( n - 1 );
    }
    else
    {
      variance[i] = 0.0;
    }

  }

  return variance;
}
//****************************************************************************80

double r8vec_circular_variance ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CIRCULAR_VARIANCE returns the circular variance of an R8VEC
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
//    Input, double X(N), the vector whose variance is desired.
//
//    Output, double R8VEC_CIRCULAR VARIANCE, the circular variance
//    of the vector entries.
//
{
  int i;
  double mean;
  double sum_c;
  double sum_s;
  double value;

  mean = r8vec_mean ( n, x );

  sum_c = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum_c = sum_c + cos ( x[i] - mean );
  }

  sum_s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum_s = sum_s + sin ( x[i] - mean );
  }

  value = sqrt ( sum_c * sum_c + sum_s * sum_s ) / ( double ) n;

  value = 1.0 - value;

  return value;
}
//****************************************************************************80

double r8vec_diff_norm ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], B[N], the vectors.
//
//    Output, double R8VEC_DIFF_NORM, the L2 norm of A - B.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

double r8vec_dot ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
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
//    Output, double R8VEC_DOT, the dot product of the vectors.
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

double r8vec_length ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LENGTH returns the Euclidean length of an R8VEC
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the vector.
//
//    Output, double R8VEC_LENGTH, the Euclidean length of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    value = value + pow ( x[i], 2 );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

double r8vec_max ( int n, double *dvec )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double *RVEC, a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double rmax;
  double *r8vec_pointer;

  if ( n <= 0 )
  {
    return 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      rmax = *dvec;
      r8vec_pointer = dvec;
    }
    else
    {
      r8vec_pointer++;
      if ( rmax < *r8vec_pointer )
      {
        rmax = *r8vec_pointer;
      }
    }
  }

  return rmax;

}
//****************************************************************************80

double r8vec_mean ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MEAN returns the mean of an R8VEC.
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
//****************************************************************************80

double r8vec_min ( int n, double *dvec )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double *RVEC, a pointer to the first entry of the array.
//
//    Output, double R8VEC_MIN, the value of the minimum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double rmin;
  double *r8vec_pointer;

  if ( n <= 0 )
  {
    return 0.0;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      rmin = *dvec;
      r8vec_pointer = dvec;
    }
    else
    {
      r8vec_pointer++;
      if ( *r8vec_pointer < rmin )
      {
        rmin = *r8vec_pointer;
      }
    }
  }

  return rmin;

}
//****************************************************************************80

void r8vec_print ( int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC
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
//    Input, char *TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(6)  << i + 1 << "  "
         << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

double *r8vec_random ( int n, double alo, double ahi, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_RANDOM returns a scaled pseudorandom R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double ALO, AHI, the range allowed for the entries.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_RANDOM[N], the vector of randomly chosen integers.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = r8_random ( alo, ahi, seed );
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
  float sum;

  sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

double *r8vec_uniform_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
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
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
//
{
  int i;
  int k;
  double *r;

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void r8vec_unit_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIT_SUM normalizes an R8VEC to have unit sum.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, double A[N], the vector to be normalized.
//    On output, the entries of A should have unit sum.  However, if
//    the input vector has zero sum, the routine halts.
//
{
  double a_sum;
  int i;

  a_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    a_sum = a_sum + a[i];
  }

  if ( a_sum == 0.0 )
  {
    cout << "\n";
    cout << "R8VEC_UNIT_SUM - Fatal error!\n";
    cout << "  The vector entries sum to 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / a_sum;
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

double rayleigh_cdf ( double x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_CDF evaluates the Rayleigh CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//    0.0 <= X.
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A.
//
//    Output, double RAYLEIGH_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x < 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 1.0 - exp ( - x * x / ( 2.0 * a * a ) );
  }

  return cdf;
}
//****************************************************************************80

double rayleigh_cdf_inv ( double cdf, double a )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_CDF_INV inverts the Rayleigh CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A.
//
//    Output, double RAYLEIGH_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "RAYLEIGH_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = sqrt ( - 2.0 * a * a * log ( 1.0 - cdf ) );

  return x;
}
//****************************************************************************80

void rayleigh_cdf_values ( int *n_data, double *sigma, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_CDF_VALUES returns some values of the Rayleigh CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = RayleighDistribution [ sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2004
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
//    Output, double *SIGMA, the shape parameter of the distribution.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 9

  double fx_vec[N_MAX] = {
     0.8646647167633873,
     0.9996645373720975,
     0.9999999847700203,
     0.999999999999987,
     0.8646647167633873,
     0.3934693402873666,
     0.1992625970831920,
     0.1175030974154046,
     0.7688365361336422E-01 };

  double sigma_vec[N_MAX] = {
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *sigma = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *sigma = sigma_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool rayleigh_check ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_CHECK checks the parameter of the Rayleigh PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A.
//
//    Output, bool RAYLEIGH_CHECK, is true if the parameter is legal.
//
{
  if ( a <= 0.0 )
  {
    cout << " \n";
    cout << "RAYLEIGH_CHECK - Fatal error!\n";
    cout << "  A <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double rayleigh_mean ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_MEAN returns the mean of the Rayleigh PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A.
//
//    Output, double RAYLEIGH_MEAN, the mean of the PDF.
//
{
  double mean;
  const double pi = 3.14159265358979323;

  mean = a * sqrt ( 0.5 * pi );

  return mean;
}
//****************************************************************************80

double rayleigh_pdf ( double x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_PDF evaluates the Rayleigh PDF.
//
//  Discussion:
//
//    PDF(A;X) = ( X / A^2 ) * EXP ( - X^2 / ( 2 * A^2 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X
//
//    Input, double A, the parameter of the PDF.
//    0 < A.
//
//    Output, double RAYLEIGH_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = ( x / ( a * a ) ) * exp ( - x * x / ( 2.0 * a * a ) );
  }

  return pdf;
}
//****************************************************************************80

double rayleigh_sample ( double a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_SAMPLE samples the Rayleigh PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    0.0 < A.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double RAYLEIGH_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = rayleigh_cdf_inv ( cdf, a );

  return x;
}
//****************************************************************************80

double rayleigh_variance ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_VARIANCE returns the variance of the Rayleigh PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameters of the PDF.
//    0.0 < A.
//
//    Output, double RAYLEIGH_VARIANCE, the variance of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double variance;

  variance = 2.0 * a * a * ( 1.0 - 0.25 * pi );

  return variance;
}
//****************************************************************************80

double reciprocal_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_CDF evaluates the Reciprocal CDF.
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
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A <= B.
//
//    Output, double RECIPROCAL_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else if ( 0.0 < x )
  {
    cdf = log ( a / x ) / log ( a / b );
  }

  return cdf;
}
//****************************************************************************80

double reciprocal_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_CDF_INV inverts the Reciprocal CDF.
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
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A <= B.
//
//    Output, double RECIPROCAL_CDF_INV, the corresponding argument of the CDF.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "RECIPROCAL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = 0.0;
  }
  else if ( 0.0 < cdf )
  {
    x = pow ( b, cdf ) / pow ( a, ( cdf - 1.0 ) );
  }

  return x;
}
//****************************************************************************80

bool reciprocal_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_CHECK checks the parameters of the Reciprocal CDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A <= B.
//
//    Output, bool RECIPROCAL_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 0.0 )
  {
    cout << " \n";
    cout << "RECIPROCAL_CHECK - Fatal error!\n";
    cout << "  A <= 0.0\n";
    return false;
  }

  if ( b < a )
  {
    cout << " \n";
    cout << "RECIPROCAL_CHECK - Fatal error!\n";
    cout << "  B < A\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double reciprocal_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_MEAN returns the mean of the Reciprocal PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A <= B.
//
//    Output, double RECIPROCAL_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = ( a - b ) / log ( a / b );

  return mean;
}
//****************************************************************************80

double reciprocal_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_PDF evaluates the Reciprocal PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = 1.0 / ( X * LOG ( B / A ) )
//    for 0.0 <= X
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
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A <= B.
//
//    Output, double RECIPROCAL_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x <= 0.0 )
  {
    pdf = 0.0;
  }
  else if ( 0.0 < x )
  {
    pdf = 1.0 / ( x * log ( b / a ) );
  }

  return pdf;
}
//****************************************************************************80

double reciprocal_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_SAMPLE samples the Reciprocal PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A <= B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double RECIPROCAL_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = pow ( b, cdf ) / pow ( a, ( cdf - 1.0 ) );

  return x;
}
//****************************************************************************80

double reciprocal_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_VARIANCE returns the variance of the Reciprocal PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A <= B.
//
//    Output, double RECIPROCAL_VARIANCE, the variance of the PDF.
//
{
  double d;
  double variance;

  d = log ( a / b );

  variance = ( a - b )* ( a * ( d - 2.0 ) + b * ( d + 2.0 ) )
    / ( 2.0 * d * d );

  return variance;
}
//****************************************************************************80

int ribesl ( double x, double alpha, int nb, int ize, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    RIBESL calculates I Bessel function with non-integer orders.
//
//  Discussion:
//
//    This routine calculates Bessel functions I SUB(N+ALPHA) (X)
//    for non-negative argument X, and non-negative order N+ALPHA,
//    with or without exponential scaling.
//
//    This program is based on a program written by David
//    Sookne that computes values of the Bessel functions J or
//    I of real argument and integer order.  Modifications include
//    the restriction of the computation to the I Bessel function
//    of non-negative real argument, the extension of the computation
//    to arbitrary positive order, the inclusion of optional
//    exponential scaling, and the elimination of most underflow.
//
//    In case of an error, NCALC will not equal NB, and not all I's are
//    calculated to the desired accuracy.
//
//    If NCALC < 0:  An argument is out of range. For example,
//    NB <= 0, IZE is not 1 or 2, or IZE = 1 and EXPARG <= ABS(X)
//    In this case, the B-vector is not calculated, and NCALC is
//    set to MIN(NB,0)-1 so that NCALC /= NB.
//
//    If 0 < NCALC < NB, then not all requested function values could
//    be calculated accurately.  This usually occurs because NB is
//    much larger than ABS(X).  In this case, B(N) is calculated
//    to the desired accuracy for N <= NCALC, but precision
//    is lost for NCALC < N <= NB.  If B(N) does not vanish
//    for NCALC < N (because it is too small to be represented),
//    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
//    significant figures of B(N) can be trusted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    Original FORTRAN77 version by William Cody.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Frank Olver, David Sookne,
//    A Note on Backward Recurrence Algorithms,
//    Mathematics of Computation,
//    Volume 26, 1972, pages 941-947.
//
//    David Sookne,
//    Bessel Functions of Real Argument and Integer Order,
//    NBS Journal of Research B,
//    Volume 77B, 1973, pages 125-132.
//
//    William Cody,
//    Algorithm 597:
//    Sequence of Modified Bessel Functions of the First Kind,
//    ACM Transactions of Mathematical Software,
//    Volume 9, Number 2, June 1983, pages 242-245.
//
//  Parameters:
//
//    Input, double X, the argument for which the functions
//    are to be calculated.
//
//    Input, double ALPHA,the fractional part of the order
//    for which the functions are to be calculated.
//    0 <= ALPHA < 1.0.
//
//    Input, int NB, the number of functions to be calculated.
//    The first function calculated is of order ALPHA, and the
//    last is of order (NB - 1 + ALPHA).  1 <= NB.
//
//    Input, int IZE, scaling option.
//    1, unscaled I's are to calculated,
//    2, exponentially scaled I's are to be calculated.
//
//    Output, double B[NB], the values of the functions
//    I(ALPHA,X) through I(NB-1+ALPHA,X), with scaling if requested.
//
//    Output, int RIBESL, the value of NCALC, the error indicator.
//    If NCALC = NB, then all the requested values were calculated
//    to the desired accuracy.
//
//  Local Parameeters:
//
//    BETA, the radix for the floating-point system.
//
//    MINEXP, smallest representable power of BETA.
//
//    MAXEXP, smallest power of BETA that overflows
//
//    IT, number of bits in the mantissa of a working precision variable.
//
//    NSIG, decimal significance desired.  Should be set to
//    INT(LOG10(2)*IT+1).  Setting NSIG lower will result
//    in decreased accuracy while setting NSIG higher will
//    increase CPU time without increasing accuracy.  The
//    truncation error is limited to a relative error of
//    T=.5*10^(-NSIG).
//
//    ENTEN, 10.0^K, where K is the largest integer such that
//    ENTEN is machine-representable in working precision
//
//    ENSIG, 10.0^NSIG
//
//    RTNSIG, 10.0^(-K) for the smallest integer K such that
//    NSIG/4 <= K.
//
//    ENMTEN, smallest ABS(X) such that X/4 does not underflow
//
//    XLARGE, upper limit on the magnitude of X when IZE=2.  Bear
//    in mind that if ABS(X)=N, then at least N iterations
//    of the backward recursion will be executed.  The value
//    of 10.0^4 is used on every machine.
//
//    EXPARG, largest working precision argument that the library
//    EXP routine can handle and upper limit on the
//    magnitude of X when IZE=1; approximately log(BETA^MAXEXP).
//
//    Approximate values for some important machines are:
//
//                        beta       minexp      maxexp       it
//
//  CRAY-1        (S.P.)    2        -8193        8191        48
//  Cyber 180/855
//    under NOS   (S.P.)    2         -975        1070        48
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)    2         -126         128        24
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)    2        -1022        1024        53
//  IBM 3033      (D.P.)   16          -65          63        14
//  VAX           (S.P.)    2         -128         127        24
//  VAX D-Format  (D.P.)    2         -128         127        56
//  VAX G-Format  (D.P.)    2        -1024        1023        53
//
//
//                        NSIG       ENTEN       ENSIG      RTNSIG
//
// CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4
// Cyber 180/855
//   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4
// IEEE (IBM/XT,
//   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
// IEEE (IBM/XT,
//   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4
// IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2
// VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
// VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5
// VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4
//
//
//                         ENMTEN      XLARGE   EXPARG
//
// CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677
// Cyber 180/855
//   under NOS   (S.P.)   1.25E-293    1.0E+4     741
// IEEE (IBM/XT,
//   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88
// IEEE (IBM/XT,
//   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709
// IBM 3033      (D.P.)   2.16D-78     1.0D+4     174
// VAX           (S.P.)   1.17E-38     1.0E+4      88
// VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88
// VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709
//
{
  double constant = 1.585;
  double em;
  double empal;
  double emp2al;
  double en;
  double enmten = 8.9E-308;
  double ensig = 1.0E+16;
  double enten = 1.0E+308;
  double exparg = 709.0;
  bool flag;
  double half = 0.5;
  double halfx;
  int i;
  int k;
  int l;
  int magx;
  int n;
  int nbmx;
  int ncalc;
  int nend;
  int nsig = 16;
  int nstart;
  double one = 1.0;
  double p;
  double plast;
  double pold;
  double psave;
  double psavel;
  double rtnsig = 1.0E-04;
  double tempa;
  double tempb;
  double tempc;
  double test;
  double total;
  double tover;
  double two = 2.0;
  double xlarge = 1.0E+04;
  double zero = 0.0;
//
//  Check for X, NB, OR IZE out of range.
//
  if ( nb <= 0 )
  {
    ncalc = i4_min ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( x < 0.0 )
  {
    ncalc = i4_min ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( alpha < 0.0 )
  {
    ncalc = i4_min ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( 1.0 <= alpha )
  {
    ncalc = i4_min ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( ize == 1 && exparg < x )
  {
    ncalc = i4_min ( nb, 0 ) - 1;
    return ncalc;
  }

  if ( ize == 2 && xlarge < x )
  {
    ncalc = i4_min ( nb, 0 ) - 1;
    return ncalc;
  }
//
//  Use 2-term ascending series for small X.
//
  ncalc = nb;
  magx = ( int ) ( x );
//
//  Initialize the forward sweep, the P-sequence of Olver.
//
  if ( rtnsig <= x )
  {
    nbmx = nb - magx;
    n = magx + 1;
    en = ( double ) ( n + n ) + ( alpha + alpha );
    plast = one;
    p = en / x;
//
//  Calculate general significance test.
//
    test = ensig + ensig;

    if ( 5 * nsig < 2 * magx )
    {
      test = sqrt ( test * p );
    }
    else
    {
      test = test / pow ( constant, magx );
    }
//
//  Calculate P-sequence until N = NB-1.  Check for possible overflow.
//
    flag = false;

    if ( 3 <= nbmx )
    {
      tover = enten / ensig;
      nstart = magx + 2;
      nend = nb - 1;

      for ( k = nstart; k <= nend; k++ )
      {
        n = k;
        en = en + two;
        pold = plast;
        plast = p;
        p = en * plast / x + pold;
//
//  To avoid overflow, divide P-sequence by TOVER.  Calculate
//  P-sequence until 1 < ABS(P).
//
        if ( tover < p )
        {
          tover = enten;
          p = p / tover;
          plast = plast / tover;
          psave = p;
          psavel = plast;
          nstart = n + 1;

          for ( ; ; )
          {
            n = n + 1;
            en = en + two;
            pold = plast;
            plast = p;
            p = en * plast / x + pold;

            if ( 1.0 < p )
            {
              break;
            }
          }

          tempb = en / x;
//
//  Calculate backward test, and find NCALC, the highest N
//  such that the test is passed.
//
          test = pold * plast / ensig;
          test = test * ( half - half / ( tempb * tempb ) );
          p = plast * tover;
          n = n - 1;
          en = en - two;
          nend = i4_min ( nb, n );

          ncalc = nend + 1;

          for ( l = nstart; l <= nend; l++ )
          {
            pold = psavel;
            psavel = psave;
            psave = en * psavel / x + pold;

            if ( test < psave * psavel )
            {
              ncalc = l;
              break;
            }
          }
          ncalc = ncalc - 1;
          flag = true;
          break;
        }
      }

      if ( !flag )
      {
        n = nend;
        en = ( double ) ( n + n ) + ( alpha + alpha );
//
//  Calculate special significance test for 2 < NBMX.
//
        test = r8_max ( test, sqrt ( plast * ensig ) * sqrt ( p + p ) );
      }
    }
//
//  Calculate P-sequence until significance test passed.
//
    if ( !flag )
    {
      for ( ; ; )
      {
        n = n + 1;
        en = en + two;
        pold = plast;
        plast = p;
        p = en * plast / x + pold;

        if ( test <= p )
        {
          break;
        }
      }
    }
//
//  Initialize the backward recursion and the normalization sum.
//
    n = n + 1;
    en = en + two;
    tempb = zero;
    tempa = one / p;
    em = ( double ) ( n ) - one;
    empal = em + alpha;
    emp2al = ( em - one ) + ( alpha + alpha );
    total = tempa * empal * emp2al / em;
    nend = n - nb;
//
//  N < NB, so store B(N) and set higher orders to zero.
//
    if ( nend < 0 )
    {
      b[n-1] = tempa;
      nend = -nend;

      for ( l = 1; l <= nend; l++ )
      {
        b[n+l-1] = zero;
      }

      nend = n - 2;
//
//  Calculate via difference equation and store B(N), until N = 2.
//
      if ( 0 < nend )
      {
        for ( l = 1; l <= nend; l++ )
        {
          n = n - 1;
          en = en - two;
          b[n-1] = ( en * b[n] ) / x + b[n+1];
          em = em - one;
          emp2al = emp2al - one;
          if ( n == 2 )
          {
            emp2al = one;
          }
          empal = empal - one;
          total = ( total + b[n-1] * empal ) * emp2al / em;
        }
      }
//
//  Calculate B(1).
//
      b[0] = two * empal * b[1] / x + b[2];

      total = ( total + total ) + b[0];
    }
//
//  Recur backward via difference equation, calculating (but
//  not storing) B(N), until N = NB.
//
    else
    {
      if ( 0 < nend )
      {
        for ( l = 1; l <= nend; l++ )
        {
          n = n - 1;
          en = en - two;
          tempc = tempb;
          tempb = tempa;
          tempa = ( en * tempb ) / x + tempc;
          em = em - one;
          emp2al = emp2al - one;

          if ( n == 1 )
          {
            break;
          }

          if ( n == 2 )
          {
            emp2al = one;
          }
          empal = empal - one;
          total = ( total + tempa * empal ) * emp2al / em;
        }
      }
//
//  Store B(NB).
//
      b[n-1] = tempa;

      if ( nb <= 1 )
      {
        total = ( total + total ) + tempa;
      }
//
//  Calculate and Store B(NB-1).
//
      else
      {
        n = n - 1;
        en = en - two;
        b[n-1] = ( en * tempa ) / x + tempb;

        if ( 1 < n  )
        {
          em = em - one;
          emp2al = emp2al - one;

          if ( n == 2 )
          {
            emp2al = one;
          }
          empal = empal - one;
          total = ( total + b[n-1] * empal ) * emp2al / em;

          nend = n - 2;
//
//  Calculate via difference equation and store B(N), until N = 2.
//
          if ( 0 < nend )
          {
            for ( l = 1; l <= nend; l++ )
            {
              n = n - 1;
              en = en - two;
              b[n-1] = ( en * b[n] ) / x + b[n+1];
              em = em - one;
              emp2al = emp2al - one;
              if ( n == 2 )
              {
                emp2al = one;
              }
              empal = empal - one;
              total = ( total + b[n-1] * empal ) * emp2al / em;
            }
          }
//
//  Calculate B(1).
//
          b[0] = two * empal * b[1] / x + b[2];
        }
        total = ( total + total ) + b[0];
      }
    }
//
//  Normalize.  Divide all B(N) by TOTAL.
//
    if ( alpha != zero )
    {
       total = total * r8_gamma ( one + alpha ) * pow ( x * half, -alpha );
    }

    if ( ize == 1 )
    {
      total = total * exp ( -x );
    }

    tempa = enmten;

    if ( 1.0 < total )
    {
      tempa = tempa * total;
    }

    for ( n = 1; n <= nb; n++ )
    {
      if ( b[n-1] < tempa )
      {
        b[n-1] = zero;
      }
      b[n-1] = b[n-1] / total;
    }

    return ncalc;
  }
//
//  Two-term ascending series for small X.
//
  else
  {
    tempa = one;
    empal = one + alpha;
    halfx = zero;

    if ( enmten < x )
    {
      halfx = half * x;
    }

    if ( alpha != zero )
    {
      tempa = pow ( halfx, alpha ) / r8_gamma ( empal );
    }

    if ( ize == 2 )
    {
      tempa = tempa * exp ( - x );
    }

    tempb = zero;

    if ( one < x + one )
    {
      tempb = halfx * halfx;
    }

    b[0] = tempa + tempa * tempb / empal;

    if ( x != zero && b[0] == zero )
    {
      ncalc = 0;
    }

    if ( 1 < nb )
    {
      if ( x == zero )
      {
        for ( i = 1; i < nb; i++ )
        {
          b[i] = zero;
        }
      }
//
//  Calculate higher-order functions.
//
      else
      {
        tempc = halfx;
        tover = ( enmten + enmten ) / x;

        if ( tempb != zero )
        {
          tover = enmten / tempb;
        }

        for ( n = 2; n <= nb; n++ )
        {
          tempa = tempa / empal;
          empal = empal + one;
          tempa = tempa * tempc;

          if ( tempa <= tover * empal )
          {
            tempa = zero;
          }

          b[n-1] = tempa + tempa * tempb / empal;

          if ( b[n-1] == zero && n < ncalc )
          {
            ncalc = n - 1;
          }
        }
      }
    }
  }

  return ncalc;
}
//****************************************************************************80

double runs_mean ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    RUNS_MEAN returns the mean of the Runs PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the parameters of the PDF.
//
//    Output, double RUNS_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = ( double ) ( m + 2 * m * n + n )
       / ( double ) ( m             + n );

  return mean;
}
//****************************************************************************80

double runs_pdf ( int m, int n, int r )

//****************************************************************************80
//
//  Purpose:
//
//    RUNS_PDF evaluates the Runs PDF.
//
//  Discussion:
//
//    Suppose we have M symbols of one type and N of another, and we consider
//    the various possible permutations of these symbols.
//
//    Let "R" be the number of runs in a given permutation.  By a "run", we
//    mean a maximal sequence of identical symbols.  Thus, for instance,
//    the permutation
//
//      ABBBAAAAAAAA
//
//    has three runs.
//
//    The probability that a permutation of M+N symbols, with M of one kind
//    and N of another, will have exactly R runs is:
//
//      PDF(M,N)(R) = 2 * C(M-1,R/2-1) * C(N-1,R/2-1)
//                    / C(M+N,N) for R even;
//
//                  = ( C(M-1,(R-1)/2) * C(N-1,(R-3)/2 )
//                    + C(M-1,(R-3)/2) * C(N-1,(R-1)/2 )
//                    ) / C(M+N,N) for R odd.
//
//    Note that the maximum number of runs for a given M and N is:
//
//      M + N,                if M = N
//      2 * min ( M, N ) + 1  otherwise
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kalimutha Krishnamoorthy,
//    Handbook of Statistical Distributions with Applications,
//    Chapman and Hall, 2006,
//    ISBN: 1-58488-635-8,
//    LC: QA273.6.K75.
//
//  Parameters:
//
//    Input, int M, N, the parameters of the PDF.
//
//    Input, int R, the number of runs.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;

  if ( m < 0 )
  {
    cout << "\n";
    cout << "RUN_PDF - Fatal error!\n";
    cout << "  M must be at least 0.\n";
    cout << "  The input value of M = " << m << "\n";
    exit ( 1 );
  }

  if ( n < 0 )
  {
    cout << "\n";
    cout << "RUN_PDF - Fatal error!\n";
    cout << "  N must be at least 0.\n";
    cout << "  The input value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( n + m <= 0 )
  {
    cout << "\n";
    cout << "RUN_PDF - Fatal error!\n";
    cout << "  M+N must be at least 1.\n";
    cout << "  The input value of M+N = " << m + n << "\n";
    exit ( 1 );
  }
//
//  If all the symbols are of one type, there is always 1 run.
//
  if ( m == 0 || n == 0 )
  {
    if ( r == 1 )
    {
      pdf = 1.0;
    }
    else
    {
      pdf = 0.0;
    }
    return pdf;
  }
//
//  Take care of extreme values of R.
//
  if ( r < 2 || m + n < r )
  {
    pdf = 0.0;
    return pdf;
  }
//
//  The normal cases.
//
  if ( ( r % 2 ) == 0 )
  {
    pdf = ( double ) ( 2 * combinatorial ( m - 1, ( r / 2 ) - 1 )
                         * combinatorial ( n - 1, ( r / 2 ) - 1 ) )
        / ( double ) (     combinatorial ( m + n, n ) );
  }
  else
  {
    pdf = ( double ) (   combinatorial ( m - 1, ( r - 1 ) / 2 )
                       * combinatorial ( n - 1, ( r - 3 ) / 2 )
                       + combinatorial ( m - 1, ( r - 3 ) / 2 )
                       * combinatorial ( n - 1, ( r - 1 ) / 2 ) )
        / ( double ) (   combinatorial ( m + n, n ) );
  }

  return pdf;
}
//****************************************************************************80

int runs_sample ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RUNS_SAMPLE samples the Runs PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the parameters of the PDF.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int RUNS_SAMPLE, the number of runs.
//
{
  int *a;
  int r;

  a = runs_simulate ( m, n, seed );

  r = i4vec_run_count ( m+n, a );

  delete [] a;

  return r;
}
//****************************************************************************80

int *runs_simulate ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    RUNS_SIMULATE simulates a case governed by the Runs PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the parameters of the PDF.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int RUNS_SIMULATE[M+N], a sequence of M 0's and N 1's chosen
//    uniformly at random.
//
{
  int *a;
  int i;
  int j;
  int k;

  a = new int[m+n];

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0;
  }
  for ( i = m; i < m + n; i++ )
  {
    a[i] = 1;
  }

  for ( i = 1; i <= m+n-1; i++ )
  {
    j = i4_uniform ( i, m+n, seed );

    k      = a[i-1];
    a[i-1] = a[j-1];
    a[j-1] = k;
  }

  return a;
}
//****************************************************************************80

double runs_variance ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    RUNS_VARIANCE returns the variance of the Runs PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the parameters of the PDF.
//
//    Output, double RUNS_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( double ) ( 2 * m * n * ( 2 * m * n - m - n ) )
           / ( double ) ( ( m + n ) * ( m + n ) * ( m + n - 1 ) );

  return variance;
}
//****************************************************************************80

int s_len_trim ( char *s )

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

double sech ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    SECH returns the hyperbolic secant.
//
//  Discussion:
//
//    SECH ( X ) = 1.0 / COSH ( X ) = 2.0 / ( EXP ( X ) + EXP ( - X ) )
//
//    SECH is not a built-in function in FORTRAN, and occasionally it
//    is handier, or more concise, to be able to refer to it directly
//    rather than through its definition in terms of the sine function.
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
//    Input, double X, the argument.
//
//    Output, double SECH, the hyperbolic secant of X.
//
{
  double value;

  value = 1.0 / cosh ( x );

  return value;
}
//****************************************************************************80

double sech_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SECH_CDF evaluates the Hyperbolic Secant CDF.
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
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameter of the PDF.
//    0.0 < B.
//
//    Output, double SECH_CDF, the value of the CDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;
  double y;

  y = ( x - a ) / b;

  cdf = 2.0 * atan ( exp ( y ) ) / pi;

  return cdf;
}
//****************************************************************************80

double sech_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SECH_CDF_INV inverts the Hyperbolic Secant CDF.
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
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double SECH_CDF_INV, the corresponding argument of the CDF.
//
{
  const double pi = 3.14159265358979323;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "SECH_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = - r8_huge (  );
  }
  else if ( cdf < 1.0 )
  {
    x = a + b * log ( tan ( 0.5 * pi * cdf ) );
  }
  else if ( 1.0 == cdf )
  {
    x = r8_huge ( );
  }

  return x;
}
//****************************************************************************80

bool sech_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SECH_CHECK checks the parameters of the Hyperbolic Secant CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameter of the PDF.
//    0.0 < B.
//
//    Output, bool SECH_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "SECH_CHECK - Fatal error!\n";
    cout << "  B <= 0.0\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double sech_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SECH_MEAN returns the mean of the Hyperbolic Secant PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double SECH_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double sech_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SECH_PDF evaluates the Hypebolic Secant PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = sech ( ( X - A ) / B ) / ( PI * B )
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
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double SECH_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;
  double y;

  y = ( x - a ) / b;

  pdf = sech ( y ) / ( pi * b );

  return pdf;
}
//****************************************************************************80

double sech_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SECH_SAMPLE samples the Hyperbolic Secant PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double SECH_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = a + b * log ( tan ( 0.5 *pi * cdf ) );

  return x;
}
//****************************************************************************80

double sech_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SECH_VARIANCE returns the variance of the Hyperbolic Secant PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double SECH_VARIANCE, the variance of the PDF.
//
{
  const double pi = 3.14159265358979323;
  double variance;

  variance = 0.25 * pi * pi * b * b;

  return variance;
}
//****************************************************************************80

double semicircular_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_CDF evaluates the Semicircular CDF.
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
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameter of the PDF.
//    0.0 < B.
//
//    Output, double SEMICIRCULAR_CDF, the value of the CDF.
//
{
  double cdf;
  const double pi = 3.14159265358979323;
  double y;

  if ( x <= a - b )
  {
    cdf = 0.0;
  }
  else if ( x <= a + b )
  {
    y = ( x - a ) / b;

    cdf = 0.5 + ( y * sqrt ( 1.0 - y * y ) + asin ( y ) ) / pi;
  }
  else if ( a + b < x )
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double semicircular_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_CDF_INV inverts the Semicircular CDF.
//
//  Discussion:
//
//    A simple bisection method is used on the interval [ A - B, A + B ].
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
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double SEMICIRCULAR_CDF_INV, the corresponding argument
//    of the CDF.
//
{
  double cdf1;
  double cdf2;
  double cdf3;
  int it;
  int it_max = 100;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "SEMICIRCULAR_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = a - b;
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = a + b;
    return x;
  }

  x1 = a - b;
  cdf1 = 0.0;

  x2 = a + b;
  cdf2 = 1.0;
//
//  Now use bisection.
//
  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = semicircular_cdf ( x3, a, b );

    if ( r8_abs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      break;
    }

    if ( it_max < it )
    {
      cout << " \n";
      cout << "SEMICIRCULAR_CDF_INV - Fatal error!\n";
      cout << "  Iteration limit exceeded.\n";
      exit ( 1 );
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
      cdf2 = cdf3;
    }

  }

  return x;
}
//****************************************************************************80

bool semicircular_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_CHECK checks the parameters of the Semicircular CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameter of the PDF.
//    0.0 < B.
//
//    Output, bool SEMICIRCULAR_CHECK, is true if the parameters are legal.
//
{
  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "SEMICIRCULAR_CHECK - Fatal error!\n";
    cout << "  B <= 0.0\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double semicircular_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_MEAN returns the mean of the Semicircular PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double SEMICIRCULAR_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double semicircular_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_PDF evaluates the Semicircular PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = ( 2 / ( B * PI ) ) * SQRT ( 1 - ( ( X - A ) / B )^2 )
//    for A - B <= X <= A + B
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
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double SEMICIRCULAR_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;
  double y;

  if ( x < a - b )
  {
    pdf = 0.0;
  }
  else if ( x <= a + b )
  {
    y = ( x - a ) / b;

    pdf = 2.0 / ( b * pi ) * sqrt ( 1.0 - y * y );
  }
  else if ( a + b < x )
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

double semicircular_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_SAMPLE samples the Semicircular PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double SEMICIRCULAR_SAMPLE, a sample of the PDF.
//
{
  double angle;
  const double pi = 3.14159265358979323;
  double radius;
  double x;

  radius = r8_uniform_01 ( seed );
  radius = b * sqrt ( radius );
  angle = pi * r8_uniform_01 ( seed );
  x = a + radius * cos ( angle );

  return x;
}
//****************************************************************************80

double semicircular_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_VARIANCE returns the variance of the Semicircular PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    0.0 < B.
//
//    Output, double SEMICIRCULAR_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = b * b / 4.0;

  return variance;
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
//    Input, int N, the power of the sine function.
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
//
  if ( n < 0 )
  {
    cout << "\n";
    cout << "SIN_POWER_INT - Fatal error!\n";
    cout << "  Power N < 0.\n";
    value = 0.0;
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

double sphere_unit_area_nd ( int dim_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_AREA_ND computes the surface area of a unit sphere in ND.
//
//  Discussion:
//
//    The unit sphere in ND satisfies the equation:
//
//      Sum ( 1 <= I <= DIM_NUM ) X(I) * X(I) = 1
//
//    DIM_NUM   Area
//
//     2    2        * PI
//     3    4        * PI
//     4  ( 2 /   1) * PI^2
//     5  ( 8 /   3) * PI^2
//     6  ( 1 /   1) * PI^3
//     7  (16 /  15) * PI^3
//     8  ( 1 /   3) * PI^4
//     9  (32 / 105) * PI^4
//    10  ( 1 /  12) * PI^5
//
//    For the unit sphere, Area(DIM_NUM) = DIM_NUM * Volume(DIM_NUM)
//
//    Sphere_Unit_Area ( DIM_NUM ) = 2 * PI^(DIM_NUM/2) / Gamma ( DIM_NUM / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Output, double SPHERE_UNIT_AREA_ND, the area of the sphere.
//
{
  double area;
  int i;
  int m;
  double pi = 3.141592653589793;

  if ( ( dim_num % 2 ) == 0 )
  {
    m = dim_num / 2;
    area = 2.0 * pow ( pi, m );
    for ( i = 1; i <= m-1; i++ )
    {
      area = area / ( ( double ) i );
    }
  }
  else
  {
    m = ( dim_num - 1 ) / 2;
    area = pow ( 2.0, dim_num ) * pow ( pi, m );
    for ( i = m+1; i <= 2*m; i++ )
    {
      area = area / ( ( double ) i );
    }
  }

  return area;
}
//****************************************************************************80

int stirling2_value ( int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    STIRLING2_VALUE computes a Stirling number of the second kind.
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
//    X**N = sum ( 0 <= K <= N ) S2(N,K) X_K
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
//    Output, int STIRLING2_VALUE, the value of S2(N,M).
//
{
  int i;
  int j;
  int *s2;
  int value;

  s2 = new int[n*m];

  if ( n <= 0 )
  {
    value = 0;
    return value;
  }

  if ( m <= 0 )
  {
    value = 0;
    return value;
  }

  s2[0+0*n] = 1;
  for ( j = 2; j <= m; j++ )
  {
    s2[0+(j-1)*n] = 0;
  }

  for ( i = 2; i <= n; i++ )
  {
    s2[i-1+0*n] = 1;
    for ( j = 2; j <= m; j++ )
    {
      s2[i-1+(j-1)*n] = j * s2[i-2+(j-1)*n] + s2[i-2+(j-2)*n];
    }
  }

  value = s2[n-1+(m-1)*n];

  delete [] s2;

  return value;
}
//****************************************************************************80

double student_cdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_CDF evaluates the central Student T CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, shape parameters of the PDF,
//    used to transform the argument X to a shifted and scaled
//    value Y = ( X - A ) / B.  It is required that B be nonzero.
//    For the standard distribution, A = 0 and B = 1.
//
//    Input, double C, is usually called the number of
//    degrees of freedom of the distribution.  C is typically an
//    integer, but that is not essential.  It is required that
//    C be strictly positive.
//
//    Output, double STUDENT_CDF, the value of the CDF.
//
{
  double a2;
  double b2;
  double c2;
  double cdf;
  double y;

  y = ( x - a ) / b;

  a2 = 0.5 * c;
  b2 = 0.5;
  c2 = c / ( c + y * y );

  if ( y <= 0.0 )
  {
    cdf = 0.5 * beta_inc ( a2, b2, c2 );
  }
  else
  {
    cdf = 1.0 - 0.5 * beta_inc ( a2, b2, c2 );
  }

  return cdf;
}
//****************************************************************************80

void student_cdf_values ( int *n_data, double *c, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_CDF_VALUES returns some values of the Student CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = StudentTDistribution [ c ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
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
//    Output, double *C, is usually called the number of
//    degrees of freedom of the distribution.  C is typically an
//    integer, but that is not essential.  It is required that
//    C be strictly positive.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 13

  double c_vec[N_MAX] = {
    1.0, 2.0, 3.0, 4.0,
    5.0, 2.0, 5.0, 2.0,
    5.0, 2.0, 3.0, 4.0,
    5.0 };

  double fx_vec[N_MAX] = {
     0.6000231200328521,
     0.6001080279134390,
     0.6001150934648930,
     0.6000995134721354,
     0.5999341989834830,
     0.7498859393137811,
     0.7500879487671045,
     0.9500004222186464,
     0.9499969138365968,
     0.9900012348724744,
     0.9900017619355059,
     0.9900004567580596,
     0.9900007637471291 };

  double x_vec[N_MAX] = {
     0.325,
     0.289,
     0.277,
     0.271,
     0.267,
     0.816,
     0.727,
     2.920,
     2.015,
     6.965,
     4.541,
     3.747,
     3.365 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *c = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *c = c_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool student_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_CHECK checks the parameter of the central Student T CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, shape parameters of the PDF,
//    used to transform the argument X to a shifted and scaled
//    value Y = ( X - A ) / B.  It is required that B be nonzero.
//    For the standard distribution, A = 0 and B = 1.
//
//    Input, double C, is usually called the number of
//    degrees of freedom of the distribution.  C is typically an
//    integer, but that is not essential.  It is required that
//    C be strictly positive.
//
{
  if ( b == 0.0 )
  {
    cout << " \n";
    cout << "STUDENT_CHECK - Fatal error!\n";
    cout << "  B must be nonzero.\n";
    return false;
  }

  if ( c <= 0.0 )
  {
    cout << " \n";
    cout << "STUDENT_CHECK - Fatal error!\n";
    cout << "  C must be greater than 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double student_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_MEAN returns the mean of the central Student T PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, shape parameters of the PDF,
//    used to transform the argument X to a shifted and scaled
//    value Y = ( X - A ) / B.  It is required that B be nonzero.
//    For the standard distribution, A = 0 and B = 1.
//
//    Input, double C, is usually called the number of
//    degrees of freedom of the distribution.  C is typically an
//    integer, but that is not essential.  It is required that
//    C be strictly positive.
//
//    Output, double STUDENT_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double student_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_PDF evaluates the central Student T PDF.
//
//  Discussion:
//
//    PDF(A,B,C;X) = Gamma ( (C+1)/2 ) /
//      ( Gamma ( C / 2 ) * Sqrt ( PI * C )
//      * ( 1 + ((X-A)/B)^2/C )^(C + 1/2 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, shape parameters of the PDF,
//    used to transform the argument X to a shifted and scaled
//    value Y = ( X - A ) / B.  It is required that B be nonzero.
//    For the standard distribution, A = 0 and B = 1.
//
//    Input, double C, is usually called the number of
//    degrees of freedom of the distribution.  C is typically an
//    integer, but that is not essential.  It is required that
//    C be strictly positive.
//
//    Output, double STUDENT_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;
  double y;

  y = ( x - a ) / b;

  pdf = r8_gamma ( 0.5 * ( c + 1.0 ) ) / ( sqrt ( pi * c )
    * r8_gamma ( 0.5 * c )
    * sqrt ( pow ( ( 1.0 + y * y / c ), ( 2 * c + 1.0 ) ) ) );

  return pdf;
}
//****************************************************************************80

double student_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_SAMPLE samples the central Student T PDF.
//
//  Discussion:
//
//    For the sampling algorithm, it is necessary that 2 < C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, shape parameters of the PDF,
//    used to transform the argument X to a shifted and scaled
//    value Y = ( X - A ) / B.  It is required that B be nonzero.
//    For the standard distribution, A = 0 and B = 1.
//
//    Input, double C, is usually called the number of
//    degrees of freedom of the distribution.  C is typically an
//    integer, but that is not essential.  It is required that
//    C be strictly positive.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double STUDENT_SAMPLE, a sample of the PDF.
//
{
  double a2;
  double b2;
  double x;
  double x2;
  double x3;

  if ( c <= 2.0 )
  {
    cout << " \n";
    cout << "STUDENT_SAMPLE - Fatal error!\n";
    cout << "  Sampling fails for C <= 2.\n";
    exit ( 1 );
  }

  a2 = 0.0;
  b2 = c  / ( c - 2.0 );

  x2 = normal_sample ( a2, b2, seed );

  x3 = chi_square_sample ( c, seed );
  x3 = x3 * c / ( c - 2.0 );

  x = a + b * x2 * sqrt ( c ) / x3;

  return x;
}
//****************************************************************************80

double student_variance ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_VARIANCE returns the variance of the central Student T PDF.
//
//  Discussion:
//
//    The variance is not defined unless 2 < C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, shape parameters of the PDF,
//    used to transform the argument X to a shifted and scaled
//    value Y = ( X - A ) / B.  It is required that B be nonzero.
//    For the standard distribution, A = 0 and B = 1.
//
//    Input, double C, is usually called the number of
//    degrees of freedom of the distribution.  C is typically an
//    integer, but that is not essential.  It is required that
//    C be strictly positive.
//
//    Output, double STUDENT_VARIANCE, the variance of the PDF.
//
{
  double variance;

  if ( c <= 2.0 )
  {
    cout << " \n";
    cout << "STUDENT_VARIANCE - Fatal error!\n";
    cout << "  Variance not defined for C <= 2.\n";
    exit ( 1 );
  }

  variance = b * b * c / ( c - 2.0 );

  return variance;
}
//****************************************************************************80

double student_noncentral_cdf ( double x, int idf, double d )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_NONCENTRAL_CDF evaluates the noncentral Student T CDF.
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
//    Original FORTRAN77 version by B E Cooper.
//    C++ version by John Burkardt
//
//  Reference:
//
//    BE Cooper,
//    The Integral of the Noncentral T-Distribution,
//    Algorithm AS 5,
//    Applied Statistics,
//    Volume 17, 1968, page 193.
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, int IDF, the number of degrees of freedom.
//
//    Input, double D, the noncentrality parameter.
//
//    Output, double STUDENT_NONCENTRAL_CDF, the value of the CDF.
//
{
  double a;
  int a_max = 100;
  double ak;
  double b;
  double cdf;
  double cdf2;
  double drb;
  double emin = 12.5;
  double f;
  double fk;
  double fmkm1;
  double fmkm2;
  int k;
  const double pi = 3.14159265358979323;
  double sum2;
  double temp;

  f = ( double ) idf;

  if ( idf == 1 )
  {
    a = x / sqrt ( f );
    b = f / ( f + x * x );
    drb = d * sqrt ( b );

    cdf2 = normal_01_cdf ( drb );
    cdf = 1.0 - cdf2 + 2.0 * tfn ( drb, a );
  }
  else if ( idf <= a_max )
  {
    a = x / sqrt ( f );
    b = f / ( f + x * x );
    drb = d * sqrt ( b );
    sum2 = 0.0;

    fmkm2 = 0.0;
    if ( r8_abs ( drb ) < emin )
    {
      cdf2 = normal_01_cdf ( a * drb );
      fmkm2 = a * sqrt ( b ) * exp ( - 0.5 * drb * drb ) * cdf2
        / sqrt ( 2.0 * pi );
    }

    fmkm1 = b * d * a * fmkm2;
    if ( r8_abs ( d ) < emin )
    {
      fmkm1 = fmkm1 + 0.5 * b * a * exp ( - 0.5 * d * d ) / pi;
    }

    if ( idf % 2 == 0 )
    {
      sum2 = fmkm2;
    }
    else
    {
      sum2 = fmkm1;
    }

    ak = 1.0;

    for ( k = 2; k <= idf - 2; k = k + 2 )
    {
      fk = ( double ) ( k );

      fmkm2 = b * ( d * a * ak * fmkm1 + fmkm2 ) * ( fk - 1.0 ) / fk;

      ak = 1.0 / ( ak * ( fk - 1.0 ) );
      fmkm1 = b * ( d * a * ak * fmkm2 + fmkm1 ) * fk / ( fk + 1.0 );

      if ( idf % 2 == 0 )
      {
        sum2 = sum2 + fmkm2;
      }
      else
      {
        sum2 = sum2 + fmkm1;
      }

      ak = 1.0 / ( ak * fk );

    }

    if ( idf % 2 == 0 )
    {
      cdf2 = normal_01_cdf ( d );
      cdf = 1.0 - cdf2 + sum2 * sqrt ( 2.0 * pi );
    }
    else
    {
      cdf2 = normal_01_cdf ( drb );
      cdf = 1.0 - cdf2 + 2.0 * ( sum2 + tfn ( drb, a ) );
    }
  }
//
//  Normal approximation.
//
  else
  {
    a = sqrt ( 0.5 * f ) * exp ( gamma_log ( 0.5 * ( f - 1.0 ) )
      - gamma_log ( 0.5 * f ) ) * d;

    temp = ( x - a ) / sqrt ( f * ( 1.0 + d * d ) / ( f - 2.0 ) - a * a );

    cdf2 = normal_01_cdf ( temp );
    cdf = cdf2;
  }

  return cdf;
}
//****************************************************************************80

void student_noncentral_cdf_values ( int *n_data, int *df, double *lambda,
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NoncentralStudentTDistribution [ df, lambda ]
//      CDF [ dist, x ]
//
//    Mathematica seems to have some difficulty computing this function
//    to the desired number of digits.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2004
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
//    Output, int *DF, double *LAMBDA, the parameters of the
//    function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 30

  int df_vec[N_MAX] = {
     1,  2,  3,
     1,  2,  3,
     1,  2,  3,
     1,  2,  3,
     1,  2,  3,
    15, 20, 25,
     1,  2,  3,
    10, 10, 10,
    10, 10, 10,
    10, 10, 10 };

  double fx_vec[N_MAX] = {
     0.8975836176504333,
     0.9522670169,
     0.9711655571887813,
     0.8231218864,
     0.9049021510,
     0.9363471834,
     0.7301025986,
     0.8335594263,
     0.8774010255,
     0.5248571617,
     0.6293856597,
     0.6800271741,
     0.20590131975,
     0.2112148916,
     0.2074730718,
     0.9981130072,
     0.9994873850,
     0.9998391562,
     0.168610566972,
     0.16967950985,
     0.1701041003,
     0.9247683363,
     0.7483139269,
     0.4659802096,
     0.9761872541,
     0.8979689357,
     0.7181904627,
     0.9923658945,
     0.9610341649,
     0.8688007350 };

  double lambda_vec[N_MAX] = {
     0.0,
     0.0,
     0.0,
     0.5,
     0.5,
     0.5,
     1.0,
     1.0,
     1.0,
     2.0,
     2.0,
     2.0,
     4.0,
     4.0,
     4.0,
     7.0,
     7.0,
     7.0,
     1.0,
     1.0,
     1.0,
     2.0,
     3.0,
     4.0,
     2.0,
     3.0,
     4.0,
     2.0,
     3.0,
     4.0 };

  double x_vec[N_MAX] = {
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
      3.00,
     15.00,
     15.00,
     15.00,
      0.05,
      0.05,
      0.05,
      4.00,
      4.00,
      4.00,
      5.00,
      5.00,
      5.00,
      6.00,
      6.00,
      6.00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *df = 0;
    *lambda = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *df = df_vec[*n_data-1];
    *lambda = lambda_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double tfn ( double h, double a )

//****************************************************************************80
//
//  Purpose:
//
//    TFN calculates the T function of Owen.
//
//  Discussion:
//
//    Owen's T function is useful for computation of the bivariate normal
//    distribution and the distribution of a skewed normal distribution.
//
//    Although it was originally formulated in terms of the bivariate
//    normal function, the function can be defined more directly as
//
//      T(H,A) = 1 / ( 2 * PI ) *
//        Integral ( 0 <= X <= A ) e^( -H^2 * (1+X^2) / 2 ) / (1+X^2) dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2004
//
//  Author:
//
//    FORTRAN77 original version by J C Young, C E Minder.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    D B Owen,
//    Tables for computing the bivariate normal distribution,
//    Annals of Mathematical Statistics,
//    Volume 27, pages 1075-1090, 1956.
//
//    J C Young, C E Minder,
//    Algorithm AS 76,
//    An Algorithm Useful in Calculating Non-Central T and
//      Bivariate Normal Distributions,
//    Applied Statistics,
//    Volume 23, Number 3, 1974, pages 455-457.
//
//  Parameters:
//
//    Input, double H, A, the arguments of the T function.
//
//    Output, double TFN, the value of the T function.
//
{
# define NGAUSS 10

  double as;
  double h1;
  double h2;
  double hs;
  int i;
  double rt;
  double two_pi_inverse = 0.1591549430918953;
  double tv1 = 1.0E-35;
  double tv2 = 15.0;
  double tv3 = 15.0;
  double tv4 = 1.0E-05;
  double value;
  double weight[NGAUSS] = {
    0.666713443086881375935688098933E-01,
    0.149451349150580593145776339658,
    0.219086362515982043995534934228,
    0.269266719309996355091226921569,
    0.295524224714752870173892994651,
    0.295524224714752870173892994651,
    0.269266719309996355091226921569,
    0.219086362515982043995534934228,
    0.149451349150580593145776339658,
    0.666713443086881375935688098933E-01 };
  double x;
  double xtab[NGAUSS] = {
   -0.973906528517171720077964012084,
   -0.865063366688984510732096688423,
   -0.679409568299024406234327365115,
   -0.433395394129247190799265943166,
   -0.148874338981631210884826001130,
    0.148874338981631210884826001130,
    0.433395394129247190799265943166,
    0.679409568299024406234327365115,
    0.865063366688984510732096688423,
    0.973906528517171720077964012084 };
//
//  Test for H near zero.
//
  if ( r8_abs ( h ) < tv1 )
  {
    value = atan ( a ) * two_pi_inverse;
  }
//
//  Test for large values of abs(H).
//
  else if ( tv2 < r8_abs ( h ) )
  {
    value = 0.0;
  }
//
//  Test for A near zero.
//
  else if ( r8_abs ( a ) < tv1 )
  {
    value = 0.0;
  }
//
//  Test whether abs(A) is so large that it must be truncated.
//  If so, the truncated value is H2.
//
  else
  {
    hs = - 0.5 * h * h;
    h2 = a;
    as = a * a;
//
//  Computation of truncation point by Newton iteration.
//
    if ( tv3 <= log ( 1.0 + as ) - hs * as )
    {
      h1 = 0.5 * a;
      as = 0.25 * as;

      for ( ; ; )
      {
        rt = as + 1.0;
        h2 = h1 + ( hs * as + tv3 - log ( rt ) )
          / ( 2.0 * h1 * ( 1.0 / rt - hs ) );
        as = h2 * h2;

        if ( r8_abs ( h2 - h1 ) < tv4 )
        {
          break;
        }
        h1 = h2;
      }
    }
//
//  Gaussian quadrature on the interval [0,H2].
//
    rt = 0.0;
    for ( i = 0; i < NGAUSS; i++ )
    {
      x = 0.5 * h2 * ( xtab[i] + 1.0 );
      rt = rt + weight[i] * exp ( hs * ( 1.0 + x * x ) ) / ( 1.0 + x * x );
    }
    value = rt * 0.5 * h2 * two_pi_inverse;
  }

  return value;
# undef NGAUSS
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

double triangle_cdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CDF evaluates the Triangle CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    A <= B <= C and A < C.
//
//    Output, double TRIANGLE_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else if ( x <= b )
  {
    if ( a == b )
    {
      cdf = 0.0;
    }
    else
    {
      cdf = ( x - a ) * ( x - a ) / ( b - a ) / ( c - a );
    }
  }
  else if ( x <= c )
  {
    cdf = ( b - a ) / ( c - a )
        + ( 2.0 * c - b - x ) * ( x - b ) / ( c - b ) / ( c - a );
  }
  else
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double triangle_cdf_inv ( double cdf, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CDF_INV inverts the Triangle CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, C, the parameters of the PDF.
//    A <= B <= C and A < C.
//
//    Output, double TRIANGLE_X, the corresponding argument.
//
{
  double cdf_mid;
  double d;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "TRIANGLE_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  d = 2.0 / ( c - a );
  cdf_mid = 0.5 * d * ( b - a );

  if ( cdf <= cdf_mid )
  {
    x = a + sqrt ( cdf * ( b - a ) * ( c - a ) );
  }
  else
  {
    x = c - sqrt ( ( c - b ) * ( ( c - b ) - ( cdf - cdf_mid ) * ( c - a ) ) );
  }

  return x;
}
//****************************************************************************80

bool triangle_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CHECK checks the parameters of the Triangle CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    A <= B <= C and A < C.
//
//    Output, bool TRIANGLE_CHECK, is true if the parameters are legal.
//
{
  if ( b < a )
  {
    cout << "\n";
    cout << "TRIANGLE_CHECK - Fatal error!\n";
    cout << "  B < A.\n";
    return false;
  }

  if ( c < b )
  {
    cout << "\n";
    cout << "TRIANGLE_CHECK - Fatal error!\n";
    cout << "  C < B.\n";
    return false;
  }

  if ( a == c )
  {
    cout << "\n";
    cout << "TRIANGLE_CHECK - Fatal error!\n";
    cout << "  A == C.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double triangle_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_MEAN returns the mean of the Triangle PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    A <= B <= C and A < C.
//
//    Output, double TRIANGLE_MEAN, the mean of the discrete uniform PDF.
//
{
  double mean;

  mean = a + ( c + b - 2.0 * a ) / 3.0;

  return mean;
}
//****************************************************************************80

double triangle_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_PDF evaluates the Triangle PDF.
//
//  Discussion:
//
//    Given points A <= B <= C, the probability is 0 to the left of A,
//    rises linearly to a maximum of 2/(C-A) at B, drops linearly to zero
//    at C, and is zero for all values greater than C.
//
//    PDF(A,B,C;X)
//      = 2 * ( X - A ) / ( B - A ) / ( C - A ) for A <= X <= B
//      = 2 * ( C - X ) / ( C - B ) / ( C - A ) for B <= X <= C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, C, the parameters of the PDF.
//    A <= B <= C and A < C.
//
//    Output, double TRIANGLE_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else if ( x <= b )
  {
    if ( a == b )
    {
      pdf = 0.0;
    }
    else
    {
      pdf = 2.0 * ( x - a ) / ( b - a ) / ( c - a );
    }
  }
  else if ( x <= c )
  {
    if ( b == c )
    {
      pdf = 0.0;
    }
    else
    {
      pdf = 2.0 * ( c - x ) / ( c - b ) / ( c - a );
    }
  }
  else
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

double triangle_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SAMPLE samples the Triangle PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    A <= B <= C and A < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = triangle_cdf_inv ( cdf, a, b, c );

  return x;
}
//****************************************************************************80

double triangle_variance ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_VARIANCE returns the variance of the Triangle PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    A <= B <= C and A < C.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( ( c - a ) * ( c - a )
             - ( c - a ) * ( b - a )
             + ( b - a ) * ( b - a ) ) / 18.0;

  return variance;
}
//****************************************************************************80

double triangular_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULAR_CDF evaluates the Triangular CDF.
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
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double TRIANGULAR_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else if ( x <= 0.5 * ( a + b ) )
  {
    cdf = 2.0 * ( x * x - 2.0 * a * x + a * a ) / pow ( ( b - a ), 2 );
  }
  else if ( x <= b )
  {
    cdf = 0.5 + ( - 2.0 * x * x + 4.0 * b * x + 0.5 * a * a
      - a * b - 1.5 * b * b ) / pow ( ( b - a ), 2 );
  }
  else
  {
    cdf = 1.0;
  }

  return cdf;
}
//****************************************************************************80

double triangular_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULAR_CDF_INV inverts the Triangular CDF.
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
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double TRIANGULAR_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "TRIANGULAR_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf <= 0.5 )
  {
    x = a + 0.5 * ( b - a ) * sqrt ( 2.0 * cdf );
  }
  else
  {
    x = b - 0.5 * ( b - a ) * sqrt ( 2.0 * ( 1.0 - cdf ) );
  }

  return x;
}
//****************************************************************************80

bool triangular_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULAR_CHECK checks the parameters of the Triangular CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, bool TRIANGULAR_CHECK, is true if the parameters are legal.
//
{
  if ( b <= a )
  {
    cout << " \n";
    cout << "TRIANGULAR_CHECK - Fatal error!\n";
    cout << "  B <= A.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double triangular_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULAR_MEAN returns the mean of the Triangular PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double TRIANGULAR_MEAN, the mean of the discrete PDF.
//
{
  double mean;

  mean = 0.5 * ( a + b );

  return mean;
}
//****************************************************************************80

double triangular_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULAR_PDF evaluates the Triangular PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = 4 * ( X - A ) / ( B - A )^2 for A <= X <= (A+B)/2
//                = 4 * ( B - X ) / ( B - A )^2 for (A+B)/2 <= X <= B.
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
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double TRIANGULAR_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else if ( x <= 0.5 * ( a + b ) )
  {
    pdf = 4.0 * ( x - a ) / ( b - a ) / ( b - a );
  }
  else if ( x <= b )
  {
    pdf = 4.0 * ( b - x ) / ( b - a ) / ( b - a );
  }
  else
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

double triangular_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULAR_SAMPLE samples the Triangular PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double TRIANGULAR_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = triangular_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double triangular_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULAR_VARIANCE returns the variance of the Triangular PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double TRIANGULAR_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( b - a ) * ( b - a ) / 24.0;

  return variance;
}
//****************************************************************************80

double trigamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    TRIGAMMA calculates the TriGamma function.
//
//  Discussion:
//
//    TriGamma(x) = d^2 log ( Gamma ( x ) ) / dx^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2000
//
//  Author:
//
//    Original FORTRAN77 version by B Schneider
//    C++ version by John Burkardt
//
//  Reference:
//
//    B Schneider,
//    Trigamma Function,
//    Algorithm AS 121,
//    Applied Statistics,
//    Volume 27, Number 1, page 97-99, 1978.
//
//  Parameters:
//
//    Input, double X, the argument of the trigamma function.
//    0 < X.
//
//    Output, double TRIGAMMA, the value of the
//    trigamma function at X.
//
{
  double a = 0.0001;
  double b = 5.0;
  double b2 =   1.0 / 6.0;
  double b4 = - 1.0 / 30.0;
  double b6 =   1.0 / 42.0;
  double b8 = - 1.0 / 30.0;
  double value;
  double y;
  double z;
//
//  1): If X is not positive, fail.
//
  if ( x <= 0.0 )
  {
    value = 0.0;
    cout << " \n";
    cout << "TRIGAMMA - Fatal error!\n";
    cout << "  X <= 0.\n";
    exit ( 1 );
  }
//
//  2): If X is smaller than A, use a small value approximation.
//
  else if ( x <= a )
  {
    value = 1.0 / x / x;
  }
//
//  3): Otherwise, increase the argument to B <= ( X + I ).
//
  else
  {
    z = x;
    value = 0.0;

    while ( z < b )
    {
      value = value + 1.0 / z / z;
      z = z + 1.0;
    }
//
//  ...and then apply an asymptotic formula.
//
    y = 1.0 / z / z;

    value = value + 0.5 *
            y + ( 1.0
          + y * ( b2
          + y * ( b4
          + y * ( b6
          + y *   b8 )))) / z;
  }

  return value;
}
//****************************************************************************80

double uniform_01_cdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_CDF evaluates the Uniform 01 CDF.
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
//    Input, double X, the argument of the CDF.
//
//    Output, double UNIFORM_01_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x < 0.0 )
  {
    cdf = 0.0;
  }
  else if ( 1.0 < x )
  {
    cdf = 1.0;
  }
  else
  {
    cdf = x;
  }

  return cdf;
}
//****************************************************************************80

double uniform_01_cdf_inv ( double cdf )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_CDF_INV inverts the Uniform 01 CDF.
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
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Output, double UNIFORM_01_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "UNIFORM_01_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = cdf;

  return x;
}
//****************************************************************************80

double uniform_01_mean ( )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_MEAN returns the mean of the Uniform 01 PDF.
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
//    Output, double UNIFORM_01_MEAN, the mean of the discrete uniform PDF.
//
{
  double mean;

  mean = 0.5;

  return mean;
}
//****************************************************************************80

double uniform_01_pdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_PDF evaluates the Uniform 01 PDF.
//
//  Discussion:
//
//    PDF(X) = 1 for 0 <= X <= 1
//           = 0 otherwise
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X <= 1.0.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 0.0 || 1.0 < x )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = 1.0;
  }

  return pdf;
}
//**********************************************************************

double uniform_01_sample ( int *seed )

//**********************************************************************
//
//  Purpose:
//
//    UNIFORM_01_SAMPLE is a random number generator.
//
//  Discussion:
//
//    SEED = SEED * (7^5) mod (2^31 - 1)
//    UNIFORM_01_SAMPLE = SEED * / ( 2^31 - 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, the integer "seed" used to generate
//    the output random number, and updated in preparation for the
//    next one.  *SEED should not be zero.
//
//    Output, double UNIFORM_01_SAMPLE, a random value between 0 and 1.
//
//  Local Parameters:
//
//    Local, int IA = 7**5.
//
//    Local, int IB = 2**15.
//
//    Local, int IB16 = 2**16.
//
//    Local, int IP = 2**31 - 1.
//
{
  static int ia = 16807;
  static int ib15 = 32768;
  static int ib16 = 65536;
  static int ip = 2147483647;
  int iprhi;
  int ixhi;
  int k;
  int leftlo;
  int loxa;
  double temp;
//
//  Don't let SEED be 0.
//
  if ( *seed == 0 )
  {
    *seed = ip;
  }
//
//  Get the 15 high order bits of SEED.
//
  ixhi = *seed / ib16;
//
//  Get the 16 low bits of SEED and form the low product.
//
  loxa = ( *seed - ixhi * ib16 ) * ia;
//
//  Get the 15 high order bits of the low product.
//
  leftlo = loxa / ib16;
//
//  Form the 31 highest bits of the full product.
//
  iprhi = ixhi * ia + leftlo;
//
//  Get overflow past the 31st bit of full product.
//
  k = iprhi / ib15;
//
//  Assemble all the parts and presubtract IP.  The parentheses are essential.
//
  *seed = ( ( ( loxa - leftlo * ib16 ) - ip ) +
    ( iprhi - k * ib15 ) * ib16 ) + k;
//
//  Add IP back in if necessary.
//
  if ( *seed < 0 )
  {
    *seed = *seed + ip;
  }
//
//  Multiply by 1 / (2**31-1).
//
  temp = ( ( double ) *seed ) * 4.656612875E-10;

  return ( temp );
}
//****************************************************************************80

double uniform_01_variance ( )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_VARIANCE returns the variance of the Uniform 01 PDF.
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
//    Output, double UNIFORM_01_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = 1.0 / 12.0;

  return variance;
}
//****************************************************************************80

double *uniform_01_order_sample ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_01_ORDER_SAMPLE samples the Uniform 01 Order PDF.
//
//  Discussion:
//
//    In effect, this routine simply generates N samples of the
//    Uniform 01 PDF; but it generates them in order.  (Actually,
//    it generates them in descending order, but stores them in
//    the array in ascending order).  This saves the work of
//    sorting the results.  Moreover, if the order statistics
//    for another PDF are desired, and the inverse CDF is available,
//    then the desired values may be generated, presorted, by
//    calling this routine and using the results as input to the
//    inverse CDF routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jerry Banks, editor,
//    Handbook of Simulation,
//    Engineering and Management Press Books, 1998, page 168.
//
//  Parameters:
//
//    Input, int N, the number of elements in the sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_01_ORDER_SAMPLE[N], N samples of the Uniform 01 PDF, in
//    ascending order.
//
{
  int i;
  double u;
  double v;
  double *x;

  x = new double[n];

  v = 1.0;

  for ( i = n-1; 0 <= i; i-- )
  {
    u = r8_uniform_01 ( seed );
    v = v * pow ( u, 1.0 / ( double ) ( i + 1 ) );
    x[i] = v;
  }

  return x;
}
//****************************************************************************80

double uniform_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_CDF evaluates the Uniform CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double UNIFORM_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x < a )
  {
    cdf = 0.0;
  }
  else if ( b < x )
  {
    cdf = 1.0;
  }
  else
  {
    cdf = ( x - a ) / ( b - a );
  }

  return cdf;
}
//****************************************************************************80

double uniform_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_CDF_INV inverts the Uniform CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double UNIFORM_CDF_INV, the corresponding argument.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "UNIFORM_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a + ( b - a ) * cdf;

  return x;
}
//****************************************************************************80

bool uniform_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_CHECK checks the parameters of the Uniform CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, bool UNIFORM_CHECK, is true if the parameters are legal.
//
{
  if ( b <= a )
  {
    cout << " \n";
    cout << "UNIFORM_CHECK - Fatal error!\n";
    cout << "  B <= A.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double uniform_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_MEAN returns the mean of the Uniform PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double UNIFORM_MEAN, the mean of the discrete uniform PDF.
//
{
  double mean;

  mean = 0.5 * ( a + b );

  return mean;
}
//****************************************************************************80

double uniform_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_PDF evaluates the Uniform PDF.
//
//  Discussion:
//
//    The Uniform PDF is also known as the "Rectangular" or "de Moivre" PDF.
//
//    PDF(A,B;X) = 1 / ( B - A ) for A <= X <= B
//               = 0 otherwise
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double UNIFORM_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < a || b < x )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = 1.0 / ( b - a );
  }

  return pdf;
}
//****************************************************************************80

double uniform_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_SAMPLE samples the Uniform PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = uniform_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double uniform_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_VARIANCE returns the variance of the Uniform PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    A < B.
//
//    Output, double UNIFORM_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( b - a ) * ( b - a ) / 12.0;

  return variance;
}
//****************************************************************************80

double uniform_discrete_cdf ( int x, int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_CDF evaluates the Uniform Discrete CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the argument of the CDF.
//
//    Input, int A, B, the parameters of the PDF.
//    A <= B.
//
//    Output, double UNIFORM_DISCRETE_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x < a )
  {
    cdf = 0.0;
  }
  else if ( b < x )
  {
    cdf = 1.0;
  }
  else
  {
    cdf = ( double ) ( x + 1 - a ) / ( double ) ( b + 1 - a );
  }

  return cdf;
}
//****************************************************************************80

int uniform_discrete_cdf_inv ( double cdf, int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_CDF_INV inverts the Uniform Discrete CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, int A, B, the parameters of the PDF.
//    A <= B.
//
//    Output, int UNIFORM_DISCRETE_CDF_INV, the smallest argument whose
//    CDF is greater than or equal to CDF.
//
{
  double a2;
  double b2;
  int x;
  double x2;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "UNIFORM_DISCRETE_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  a2 = ( double ) ( a ) - 0.5;
  b2 = ( double ) ( b ) + 0.5;
  x2 = a + cdf * ( b2 - a2 );

  x = r8_nint ( x2 );

  x = i4_max ( x, a );
  x = i4_min ( x, b );

  return x;
}
//****************************************************************************80

bool uniform_discrete_check ( int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_CHECK checks the parameters of the Uniform discrete CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, B, the parameters of the PDF.
//    A <= B.
//
//    Output, bool UNIFORM_DISCRETE_CHECK, is true if the parameters
//    are legal.
//
{
  if ( b < a )
  {
    cout << " \n";
    cout << "UNIFORM_DISCRETE_CHECK - Fatal error!\n";
    cout << "  B < A.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double uniform_discrete_mean ( int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_MEAN returns the mean of the Uniform discrete PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, B, the parameters of the PDF.
//    A <= B.
//
//    Output, double UNIFORM_DISCRETE_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = 0.5 * ( double ) ( a + b );

  return mean;
}
//****************************************************************************80

double uniform_discrete_pdf ( int x, int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_PDF evaluates the Uniform discrete PDF.
//
//  Discussion:
//
//    The Uniform Discrete PDF is also known as the "Rectangular"
//    Discrete PDF.
//
//    PDF(A,B;X) = 1 / ( B + 1 - A ) for A <= X <= B.
//
//    The parameters define the interval of integers
//    for which the PDF is nonzero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the argument of the PDF.
//
//    Input, int A, B, the parameters of the PDF.
//    A <= B.
//
//    Output, double UNIFORM_DISCRETE_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < a || b < x )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = 1.0 / ( double ) ( b + 1 - a );
  }

  return pdf;
}
//****************************************************************************80

int uniform_discrete_sample ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_SAMPLE samples the Uniform discrete PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, B, the parameters of the PDF.
//    A <= B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int UNIFORM_DISCRETE_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  int x;

  cdf = r8_uniform_01 ( seed );

  x = uniform_discrete_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double uniform_discrete_variance ( int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_DISCRETE_VARIANCE returns the variance of the Uniform discrete PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, B, the parameters of the PDF.
//    A <= B.
//
//    Output, double UNIFORM_DISCRETE_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = ( double ) ( ( b + 1 - a ) * ( b + 1 - a ) - 1 ) / 12.0;

  return variance;
}
//****************************************************************************80

double *uniform_nsphere_sample ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_NSPHERE_SAMPLE samples the Uniform Unit Sphere PDF.
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
//  Reference:
//
//    Jerry Banks, editor,
//    Handbook of Simulation,
//    Engineering and Management Press Books, 1998, page 168.
//
//  Parameters:
//
//    Input, int N, the dimension of the sphere.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_NSPHERE_SAMPLE[N], a point on the unit N
//    sphere, chosen with a uniform probability.
//
{
  int i;
  double sum;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = normal_01_sample ( seed );
  }

  sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + x[i] * x[i];
  }
  sum = sqrt ( sum );

  for ( i = 0; i < n; i++ )
  {
    x[i] = x[i] / sum;
  }

  return x;
}
//****************************************************************************80

double von_mises_cdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_CDF evaluates the von Mises CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Geoffrey Hill
//    C++ version by John Burkardt
//
//  Reference:
//
//    Geoffrey Hill,
//    ACM TOMS Algorithm 518,
//    Incomplete Bessel Function I0: The von Mises Distribution,
//    ACM Transactions on Mathematical Software,
//    Volume 3, Number 3, September 1977, pages 279-284.
//
//    Kanti Mardia, Peter Jupp,
//    Directional Statistics,
//    Wiley, 2000, QA276.M335
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//    A - PI <= X <= A + PI.
//
//    Input, double A, B, the parameters of the PDF.
//    -PI <= A <= PI,
//    0.0 < B.
//
//    Output, double VON_MISES_CDF, the value of the CDF.
//
{
  double a1 = 12.0;
  double a2 = 0.8;
  double a3 = 8.0;
  double a4 = 1.0;
  double arg;
  double c;
  double c1 = 56.0;
  double cdf;
  double ck = 10.5;
  double cn;
  double erfx;
  int ip;
  int n;
  double p;
  const double pi = 3.14159265358979323;
  double r;
  double s;
  double sn;
  double u;
  double v;
  double y;
  double z;
//
//  We expect -PI <= X - A <= PI.
//
  if ( x - a <= -pi )
  {
    cdf = 0.0;
    return cdf;
  }

  if ( pi <= x - a )
  {
    cdf = 1.0;
    return cdf;
  }
//
//  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ).
//
  z = b;

  u = r8_modp ( x - a + pi, 2.0 * pi );

  if ( u < 0.0 )
  {
    u = u + 2.0 * pi;
  }

  y = u - pi;
//
//  For small B, sum IP terms by backwards recursion.
//
  if ( z <= ck )
  {
    v = 0.0;

    if ( 0.0 < z )
    {
      ip = ( int ) ( z * a2 - a3 / ( z + a4 ) + a1 );
      p = ( double ) ( ip );
      s = sin ( y );
      c = cos ( y );
      y = p * y;
      sn = sin ( y );
      cn = cos ( y );
      r = 0.0;
      z = 2.0 / z;

      for ( n = 2; n <= ip; n++ )
      {
        p = p - 1.0;
        y = sn;
        sn = sn * c - cn * s;
        cn = cn * c + y * s;
        r = 1.0 / ( p * z + r );
        v = ( sn / p + v ) * r;
      }
    }
    cdf = ( u * 0.5 + v ) / pi;
  }
//
//  For large B, compute the normal approximation and left tail.
//
  else
  {
    c = 24.0 * z;
    v = c - c1;
    r = sqrt ( ( 54.0 / ( 347.0 / v + 26.0 - c ) - 6.0 + c ) / 12.0 );
    z = sin ( 0.5 * y ) * r;
    s = 2.0 * z * z;
    v = v - s + 3.0;
    y = ( c - s - s - 16.0 ) / 3.0;
    y = ( ( s + 1.75 ) * s + 83.5 ) / v - y;
    arg = z * ( 1.0 - s / y / y );
    erfx = error_f ( arg );
    cdf = 0.5 * erfx + 0.5;
  }

  cdf = r8_max ( cdf, 0.0 );
  cdf = r8_min ( cdf, 1.0 );

  return cdf;
}
//****************************************************************************80

double von_mises_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_CDF_INV inverts the von Mises CDF.
//
//  Discussion:
//
//    A simple bisection method is used on the interval [ A - PI, A + PI ].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//
//    Input, double A, B, the parameters of the PDF.
//    -PI <= A <= PI,
//    0.0 < B.
//
//    Output, double VON_MISES_CDF_INV, the corresponding argument of the CDF.
//    A - PI <= X <= A + PI.
//
{
  double cdf1;
  double cdf2;
  double cdf3;
  int it;
  int it_max = 100;
  const double pi = 3.14159265358979323;
  double tol = 0.0001;
  double x;
  double x1;
  double x2;
  double x3;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "VON_MISES_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  if ( cdf == 0.0 )
  {
    x = a - pi;
    return x;
  }
  else if ( 1.0 == cdf )
  {
    x = a + pi;
    return x;
  }

  x1 = a - pi;
  cdf1 = 0.0;

  x2 = a + pi;
  cdf2 = 1.0;
//
//  Now use bisection.
//
  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    x3 = 0.5 * ( x1 + x2 );
    cdf3 = von_mises_cdf ( x3, a, b );

    if ( r8_abs ( cdf3 - cdf ) < tol )
    {
      x = x3;
      break;
    }

    if ( it_max < it )
    {
      cout << " \n";
      cout << "VON_MISES_CDF_INV - Fatal error!\n";
      cout << "  Iteration limit exceeded.\n";
      exit ( 1 );
    }

    if ( ( cdf3 <= cdf && cdf1 <= cdf ) || ( cdf <= cdf3 && cdf <= cdf1 ) )
    {
      x1 = x3;
      cdf1 = cdf3;
    }
    else
    {
      x2 = x3;
      cdf2 = cdf3;
    }
  }

  return x;
}
//****************************************************************************80

void von_mises_cdf_values ( int *n_data, double *a, double *b, double *x,
  double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_CDF_VALUES returns some values of the von Mises CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kanti Mardia, Peter Jupp,
//    Directional Statistics,
//    Wiley, 2000, QA276.M335
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *A, B, the parameters of the function.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 23

  double a_vec[N_MAX] = {
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.1E+01,
     0.1E+01,
     0.1E+01,
     0.1E+01,
     0.1E+01,
     0.1E+01,
    -0.2E+01,
    -0.1E+01,
     0.0E+01,
     0.1E+01,
     0.2E+01,
     0.3E+01,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0 };

  double b_vec[N_MAX] = {
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.3E+01,
      0.3E+01,
      0.3E+01,
      0.3E+01,
      0.3E+01,
      0.3E+01,
      0.0,
      0.1E+01,
      0.2E+01,
      0.3E+01,
      0.4E+01,
      0.5E+01 };

  double fx_vec[N_MAX] = {
     0.2535089956281180E-01,
     0.1097539041177346,
     0.5000000000000000,
     0.8043381312498558,
     0.9417460124555197,
     0.5000000000000000,
     0.6018204118446155,
     0.6959356933122230,
     0.7765935901304593,
     0.8410725934916615,
     0.8895777369550366,
     0.9960322705517925,
     0.9404336090170247,
     0.5000000000000000,
     0.5956639098297530E-01,
     0.3967729448207649E-02,
     0.2321953958111930E-03,
     0.6250000000000000,
     0.7438406999109122,
     0.8369224904294019,
     0.8941711407897124,
     0.9291058600568743,
     0.9514289900655436 };

  double x_vec[N_MAX] = {
     -0.2617993977991494E+01,
     -0.1570796326794897E+01,
      0.0000000000000000,
      0.1047197551196598E+01,
      0.2094395102393195E+01,
      0.1000000000000000E+01,
      0.1200000000000000E+01,
      0.1400000000000000E+01,
      0.1600000000000000E+01,
      0.1800000000000000E+01,
      0.2000000000000000E+01,
      0.0000000000000000,
      0.0000000000000000,
      0.0000000000000000,
      0.0000000000000000,
      0.0000000000000000,
      0.0000000000000000,
      0.7853981633974483,
      0.7853981633974483,
      0.7853981633974483,
      0.7853981633974483,
      0.7853981633974483,
      0.7853981633974483 };

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
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *a = a_vec[*n_data-1];
    *b = b_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool von_mises_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_CHECK checks the parameters of the von Mises PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    -PI <= A <= PI,
//    0.0 < B.
//
//    Output, bool VON_MISES_CHECK, is true if the parameters are legal.
//
{
  const double pi = 3.14159265358979323;

  if ( a < - pi || pi < a )
  {
    cout << " \n";
    cout << "VON_MISES_CHECK - Fatal error!\n";
    cout << "  A < -PI or PI < A.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "VON_MISES_CHECK - Fatal error!\n";
    cout << "  B <= 0.0\n";
    return false;
  }

  return true;
#
}
//****************************************************************************80

double von_mises_circular_variance ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_CIRCULAR_VARIANCE returns the circular variance of the von Mises PDF.
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
//    Input, double A, B, the parameters of the PDF.
//    -PI <= A <= PI,
//    0.0 < B.
//
//    Output, double VON_MISES_CIRCULAR_VARIANCE, the circular variance of the PDF.
//
{
  double value;

  value = 1.0 - bessel_i1 ( b ) / bessel_i0 ( b );

  return value;
}
//****************************************************************************80

double von_mises_mean ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_MEAN returns the mean of the von Mises PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    -PI <= A <= PI,
//    0.0 < B.
//
//    Output, double VON_MISES_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = a;

  return mean;
}
//****************************************************************************80

double von_mises_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_PDF evaluates the von Mises PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = EXP ( B * COS ( X - A ) ) / ( 2 * PI * I0(B) )
//
//    where:
//
//      I0(*) is the modified Bessel function of the first
//      kind of order 0.
//
//    The von Mises distribution for points on the unit circle is
//    analogous to the normal distribution of points on a line.
//    The variable X is interpreted as a deviation from the angle A,
//    with B controlling the amount of dispersion.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jerry Banks, editor,
//    Handbook of Simulation,
//    Engineering and Management Press Books, 1998, page 160.
//
//    D J Best, N I Fisher,
//    Efficient Simulation of the von Mises Distribution,
//    Applied Statistics,
//    Volume 28, Number 2, pages 152-157.
//
//    Kanti Mardia, Peter Jupp,
//    Directional Statistics,
//    Wiley, 2000, QA276.M335
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A - PI <= X <= A + PI.
//
//    Input, double A, B, the parameters of the PDF.
//    -PI <= A <= PI,
//    0.0 < B.
//
//    Output, double VON_MISES_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  if ( x < a - pi )
  {
    pdf = 0.0;
  }
  else if ( x <= a + pi )
  {
    pdf = exp ( b * cos ( x - a ) ) / ( 2.0 * pi * bessel_i0 ( b ) );
  }
  else
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80

double von_mises_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_SAMPLE samples the von Mises PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    D J Best, N I Fisher,
//    Efficient Simulation of the von Mises Distribution,
//    Applied Statistics,
//    Volume 28, Number 2, pages 152-157.
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    -PI <= A <= PI,
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double c;
  double f;
  const double pi = 3.14159265358979323;
  double r;
  double rho;
  double tau;
  double u1;
  double u2;
  double u3;
  double x;
  double z;

  tau = 1.0 + sqrt ( 1.0 + 4.0 * b * b );
  rho = ( tau - sqrt ( 2.0 * tau ) ) / ( 2.0 * b );
  r = ( 1.0 + rho * rho ) / ( 2.0 * rho );

  for ( ; ; )
  {
    u1 = r8_uniform_01 ( seed );
    z = cos ( pi * u1 );
    f = ( 1.0 + r * z ) / ( r + z );
    c = b * ( r - f );

    u2 = r8_uniform_01 ( seed );

    if ( u2 < c * ( 2.0 - c ) )
    {
      break;
    }

    if ( c <= log ( c / u2 ) + 1.0 )
    {
      break;
    }

  }

  u3 = r8_uniform_01 ( seed );

  x = a + r8_sign ( u3 - 0.5 ) * acos ( f );

  return x;
}
//****************************************************************************80

double weibull_cdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_CDF evaluates the Weibull CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//    A <= X.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;
  double y;

  if ( x < a )
  {
    cdf = 0.0;
  }
  else
  {
    y = ( x - a ) / b;
    cdf = 1.0 - 1.0 / exp ( pow ( y, c ) );
  }

  return cdf;
}
//****************************************************************************80

double weibull_cdf_inv ( double cdf, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_CDF_INV inverts the Weibull CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 < CDF < 1.0.
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double X, the corresponding argument of the CDF.
//
{
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << "\n";
    cout << "WEIBULL_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = a + b * pow ( ( - log ( 1.0 - cdf ) ), ( 1.0 / c ) );

  return x;
}
//****************************************************************************80

void weibull_cdf_values ( int *n_data, double *alpha, double *beta,
  double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_CDF_VALUES returns some values of the Weibull CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = WeibullDistribution [ alpha, beta ]
//      CDF [ dist, x ]
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
//    Output, double *ALPHA, the first parameter of the distribution.
//
//    Output, double *BETA, the second parameter of the distribution.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 12

  double alpha_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  double beta_vec[N_MAX] = {
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.5000000000000000,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  double fx_vec[N_MAX] = {
     0.8646647167633873,
     0.9816843611112658,
     0.9975212478233336,
     0.9996645373720975,
     0.6321205588285577,
     0.4865828809674080,
     0.3934693402873666,
     0.3296799539643607,
     0.8946007754381357,
     0.9657818816883340,
     0.9936702845725143,
     0.9994964109502630 };

  double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *alpha = 0.0;
    *beta = 0.0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *alpha = alpha_vec[*n_data-1];
    *beta = beta_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool weibull_check ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_CHECK checks the parameters of the Weibull CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, bool WEIBULL_CHECK, is true if the parameters are legal.
//
{

  if ( b <= 0.0 )
  {
    cout << "\n";
    cout << "WEIBULL_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  if ( c <= 0.0 )
  {
    cout << "\n";
    cout << "WEIBULL_CHECK - Fatal error!\n";
    cout << "  C <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double weibull_mean ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_MEAN returns the mean of the Weibull PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double MEAN, the mean of the PDF.
//
{
  double mean;

  mean = b * r8_gamma ( ( c + 1.0 ) / c ) + a;

  return mean;
}
//****************************************************************************80

double weibull_pdf ( double x, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_PDF evaluates the Weibull PDF.
//
//  Discussion:
//
//    PDF(A,B,C;X) = ( C / B ) * ( ( X - A ) / B )**( C - 1 )
//     * EXP ( - ( ( X - A ) / B )**C ).
//
//    The Weibull PDF is also known as the Frechet PDF.
//
//    WEIBULL_PDF(A,B,1;X) is the Exponential PDF.
//
//    WEIBULL_PDF(0,1,2;X) is the Rayleigh PDF.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A <= X
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  double y;

  if ( x < a )
  {
    pdf = 0.0;
  }
  else
  {
    y = ( x - a ) / b;

    pdf = ( c / b ) * pow ( y, ( c - 1.0 ) )  / exp ( pow ( y, c ) );
  }

  return pdf;
}
//****************************************************************************80

double weibull_sample ( double a, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_SAMPLE samples the Weibull PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = weibull_cdf_inv ( cdf, a, b, c );

  return x;
}
//****************************************************************************80

double weibull_variance ( double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_VARIANCE returns the variance of the Weibull PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Output, double VARIANCE, the variance of the PDF.
//
{
  double g1;
  double g2;
  double variance;

  g1 = r8_gamma ( ( c + 2.0 ) / c );
  g2 = r8_gamma ( ( c + 1.0 ) / c );

  variance = b * b * ( g1 - g2 * g2 );

  return variance;
}
//****************************************************************************80

double weibull_discrete_cdf ( int x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_CDF evaluates the Discrete Weibull CDF.
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
//    Input, int X, the argument of the CDF.
//    0 <= X.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A <= 1.0,
//    0.0 < B.
//
//    Output, double WEIBULL_DISCRETE_CDF, the value of the CDF.
//
{
  double cdf;

  if ( x < 0 )
  {
    cdf = 0.0;
  }
  else
  {
    cdf = 1.0 - pow ( 1.0 - a, pow ( x + 1, b ) );
  }

  return cdf;
}
//****************************************************************************80

int weibull_discrete_cdf_inv ( double cdf, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_CDF_INV inverts the Discrete Weibull CDF.
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
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A <= 1.0,
//    0.0 < B.
//
//    Output, int WEIBULL_DISCRETE_CDF_INV, the corresponding argument.
//
{
  int x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cout << " \n";
    cout << "WEIBULL_DISCRETE_CDF_INV - Fatal error!\n";
    cout << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x = r8_ceiling (
    pow ( log ( 1.0 - cdf ) / log ( 1.0 - a ), 1.0 / b ) - 1.0 );

  return x;
}
//****************************************************************************80

bool weibull_discrete_check ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_CHECK checks the parameters of the discrete Weibull CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A <= 1.0,
//    0.0 < B.
//
//    Output, bool WEIBULL_DISCRETE_CHECK, is true if the parameters
//    are legal.
//
{
  if ( a < 0.0 || 1.0 < a )
  {
    cout << " \n";
    cout << "WEIBULL_DISCRETE_CHECK - Fatal error!\n";
    cout << "  A < 0 or 1 < A.\n";
    return false;
  }

  if ( b <= 0.0 )
  {
    cout << " \n";
    cout << "WEIBULL_DISCRETE_CHECK - Fatal error!\n";
    cout << "  B <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double weibull_discrete_pdf ( int x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_PDF evaluates the discrete Weibull PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = ( 1 - A )**X**B - ( 1 - A )**(X+1)**B.
//
//    WEIBULL_DISCRETE_PDF(A,1;X) = GEOMETRIC_PDF(A;X)
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
//    Input, int X, the argument of the PDF.
//    0 <= X
//
//    Input, double A, B, the parameters that define the PDF.
//    0 <= A <= 1,
//    0 < B.
//
//    Output, double WEIBULL_DISCRETE_PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 0 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = pow ( 1.0 - a , pow ( x, b ) ) - pow ( 1.0 - a, pow ( x + 1, b ) );
  }

  return pdf;
}
//****************************************************************************80

int weibull_discrete_sample ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_SAMPLE samples the discrete Weibull PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 <= A <= 1.0,
//    0.0 < B.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int X, a sample of the PDF.
//
{
  double cdf;
  int x;

  cdf = r8_uniform_01 ( seed );

  x = weibull_discrete_cdf_inv ( cdf, a, b );

  return x;
}
//****************************************************************************80

double zeta ( double p )

//****************************************************************************80
//
//  Purpose:
//
//    ZETA estimates the Riemann Zeta function.
//
//  Discussion:
//
//    For 1 < P, the Riemann Zeta function is defined as:
//
//      ZETA ( P ) = Sum ( 1 <= N < Infinity ) 1 / N**P
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
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
//    P must be greater than 1.  For integral P up to 20, a
//    precomputed value of ZETA is returned; otherwise the infinite
//    sum is approximated.
//
//    Output, double ZETA, an approximation to the Riemann
//    Zeta function.
//
{
  int n;
  const double pi = 3.14159265358979323;
  double value;
  double zsum;
  double zsum_old;

  if ( p <= 1.0 )
  {
    cout << " \n";
    cout << "ZETA - Fatal error!\n";
    cout << "  Exponent P <= 1.0.\n";
    exit ( 1 );
  }
  else if ( p == 2.0 )
  {
    value = pow ( pi, 2 ) / 6.0;
  }
  else if ( p == 3.0 )
  {
    value = 1.2020569032;
  }
  else if ( p == 4.0 )
  {
    value = pow ( pi, 4 ) / 90.0;
  }
  else if ( p == 5.0 )
  {
    value = 1.0369277551;
  }
  else if ( p == 6.0 )
  {
    value = pow ( pi, 6 ) / 945.0;
  }
  else if ( p == 7.0 )
  {
    value = 1.0083492774;
  }
  else if ( p == 8.0 )
  {
    value = pow ( pi, 8 ) / 9450.0;
  }
  else if ( p == 9.0 )
  {
    value = 1.0020083928;
  }
  else if ( p == 10.0 )
  {
    value = pow ( pi, 10 ) / 93555.0;
  }
  else if ( p == 11.0 )
  {
    value = 1.0004941886;
  }
  else if ( p == 12.0 )
  {
    value = 1.0002460866;
  }
  else if ( p == 13.0 )
  {
    value = 1.0001227133;
  }
  else if ( p == 14.0 )
  {
    value = 1.0000612482;
  }
  else if ( p == 15.0 )
  {
    value = 1.0000305882;
  }
  else if ( p == 16.0 )
  {
    value = 1.0000152823;
  }
  else if ( p == 17.0 )
  {
    value = 1.0000076372;
  }
  else if ( p == 18.0 )
  {
    value = 1.0000038173;
  }
  else if ( p == 19.0 )
  {
    value = 1.0000019082;
  }
  else if ( p == 20.0 )
  {
    value = 1.0000009540;
  }
  else
  {
    zsum = 0.0;
    n = 0;

    for ( ; ; )
    {
      n = n + 1;
      zsum_old = zsum;
      zsum = zsum + 1.0 / pow ( ( double ) n, p );
      if ( zsum <= zsum_old )
      {
        break;
      }
    }
    value = zsum;
  }

  return value;
}
//****************************************************************************80

double zipf_cdf ( int x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_CDF evaluates the Zipf CDF.
//
//  Discussion:
//
//    Simple summation is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the argument of the PDF.
//    1 <= N
//
//    Input, double A, the parameter of the PDF.
//    1.0 < A.
//
//    Output, double CDF, the value of the CDF.
//
{
  double c;
  double cdf;
  double pdf;
  int y;

  if ( x < 1 )
  {
    cdf = 0.0;
  }
  else
  {
    c = zeta ( a );

    cdf = 0.0;
    for ( y = 1; y <= x; y++ )
    {
      pdf = 1.0 / pow ( y, a ) / c;
      cdf = cdf + pdf;
    }

  }

  return cdf;
}
//****************************************************************************80

bool zipf_check ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_CHECK checks the parameter of the Zipf PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    1.0 < A.
//
//    Output, bool ZIPF_CHECK, is true if the parameters are legal.
//
{
  if ( a <= 1.0 )
  {
    cout << " \n";
    cout << "ZIPF_CHECK - Fatal error!\n";
    cout << "  A <= 1.\n";
    return false;
  }

  return true;
}
//****************************************************************************80

double zipf_mean ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_MEAN returns the mean of the Zipf PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    1.0 < A.
//
//    Output, double ZIPF_MEAN, the mean of the PDF.
//    The mean is only defined for 2 < A.
//
{
  double mean;

  if ( a <= 2.0 )
  {
    cout << " \n";
    cout << "ZIPF_MEAN - Fatal error!\n";
    cout << "  No mean defined for A <= 2.\n";
    exit ( 1 );
  }

  mean = zeta ( a - 1.0 ) / zeta ( a );

  return mean;
}
//****************************************************************************80

double zipf_pdf ( int x, double a )

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_PDF evaluates the Zipf PDF.
//
//  Discussion:
//
//    PDF(A;X) = ( 1 / X**A ) / C
//
//    where the normalizing constant is chosen so that
//
//    C = Sum ( 1 <= I < Infinity ) 1 / I**A.
//
//    From observation, the frequency of different words in long
//    sequences of text seems to follow the Zipf PDF, with
//    parameter A slightly greater than 1.  The Zipf PDF is sometimes
//    known as the "discrete Pareto" PDF.
//
//    Lotka's law is a version of the Zipf PDF in which A is 2 or approximately
//    2.  Lotka's law describes the frequency of publications by authors in a
//    given field, and estimates that the number of authors with X papers is
//    about 1/X^A of the number of authors with 1 paper.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alfred Lotka,
//    The frequency distribution of scientific productivity,
//    Journal of the Washington Academy of Sciences,
//    Volume 16, Number 12, 1926, pages 317-324.
//
//  Parameters:
//
//    Input, int X, the argument of the PDF.
//    1 <= N
//
//    Input, double A, the parameter of the PDF.
//    1.0 < A.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;

  if ( x < 1 )
  {
    pdf = 0.0;
  }
  else
  {
    pdf = 1.0 / pow ( x, a ) / zeta ( a );
  }

  return pdf;
}
//****************************************************************************80

int zipf_sample ( double a, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_SAMPLE samples the Zipf PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer Verlag, 1986, pages 550-551.
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    1.0 < A.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int ZIPF_SAMPLE, a sample of the PDF.
//
{
  double b;
  double t;
  double u;
  double v;
  double w;
  int x;

  b = pow ( 2.0, ( a - 1.0 ) );

  for ( ; ; )
  {
    u = r8_uniform_01 ( seed );
    v = r8_uniform_01 ( seed );
    w = ( int ) ( 1.0 / pow ( u, 1.0 / ( a - 1.0 ) ) );

    t = pow ( ( w + 1.0 ) / w, a - 1.0 );

    if ( v * w * ( t - 1.0 ) * b <= t * ( b - 1.0 ) )
    {
      break;
    }

  }

  x = ( int ) w;

  return x;
}
//****************************************************************************80

double zipf_variance ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_VARIANCE returns the variance of the Zipf PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, the parameter of the PDF.
//    1.0 < A.
//
//    Output, double ZIPF_VARIANCE, the variance of the PDF.
//    The variance is only defined for 3 < A.
//
{
  double mean;
  double variance;

  if ( a <= 3.0 )
  {
    cout << " \n";
    cout << "ZIPF_VARIANCE - Fatal error!\n";
    cout << "  No variance defined for A <= 3.0.\n";
    exit ( 1 );
  }

  mean = zipf_mean ( a );

  variance = zeta ( a - 2.0 ) / zeta ( a ) - mean * mean;

  return variance;
}
