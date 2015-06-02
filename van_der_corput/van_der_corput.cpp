# include <cstdlib>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <cstring>
# include <ctime>
# include <cmath>

using namespace std;

# include "van_der_corput.hpp"

//
//  These variables, which define the van der corput base and seed,
//  are accessible to the user via some routines.
//
static int van_der_corput_BASE = 2;
static int van_der_corput_SEED = 1;

//****************************************************************************80

double *circle_unit_van_der_corput ( )

//****************************************************************************80
//
//  Purpose:
//
//    CIRCLE_UNIT_VAN_DER_CORPUT picks a van der Corput point on the unit circle.
//
//  Discussion:
//
//    This routine computes the "next" van der Corput number U, converts it
//    to an angle between 0 and 2 PI, and determines the corresponding
//    X and Y coordinates on the circle.
//
//    You can get or set the van der Corput seed, which determines the next
//    value, by calling VAN_DER_CORPUT_SEED_GET or VAN_DER_CORPUT_SEED_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CIRCLE_UNIT_VAN_DER_CORPUT[2], the next van der Corput
//    point on the circle.
//
{
  double angle;
  double pi = 3.141592653589793;
  double u;
  double *x;

  u = van_der_corput ( );

  angle = 2.0 * pi * u;

  x = new double[2];

  x[0] = cos ( angle );
  x[1] = sin ( angle );

  return x;
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
//    17 November 2004
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
  int i;
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
//  Remap ISEED from [1,43200] to [1,IMAX].
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
# undef I_MAX
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

double i4_to_van_der_corput ( int seed, int base )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_VAN_DER_CORPUT computes an element of a van der Corput sequence.
//
//  Discussion:
//
//    The van der Corput sequence is often used to generate a "subrandom"
//    sequence of points which have a better covering property
//    than pseudorandom points.
//
//    The van der Corput sequence generates a sequence of points in [0,1]
//    which (theoretically) never repeats.  Except for SEED = 0, the
//    elements of the van der Corput sequence are strictly between 0 and 1.
//
//    The van der Corput sequence writes an integer in a given base B,
//    and then its digits are "reflected" about the decimal point.
//    This maps the numbers from 1 to N into a set of numbers in [0,1],
//    which are especially nicely distributed if N is one less
//    than a power of the base.
//
//    Hammersley suggested generating a set of N nicely distributed
//    points in two dimensions by setting the first component of the
//    Ith point to I/N, and the second to the van der Corput
//    value of I in base 2.
//
//    Halton suggested that in many cases, you might not know the number
//    of points you were generating, so Hammersley's formulation was
//    not ideal.  Instead, he suggested that to generated a nicely
//    distributed sequence of points in M dimensions, you simply
//    choose the first M primes, P(1:M), and then for the J-th component of
//    the I-th point in the sequence, you compute the van der Corput
//    value of I in base P(J).
//
//    Thus, to generate a Halton sequence in a 2 dimensional space,
//    it is typical practice to generate a pair of van der Corput sequences,
//    the first with prime base 2, the second with prime base 3.
//    Similarly, by using the first K primes, a suitable sequence
//    in K-dimensional space can be generated.
//
//    The generation is quite simple.  Given an integer SEED, the expansion
//    of SEED in base BASE is generated.  Then, essentially, the result R
//    is generated by writing a decimal point followed by the digits of
//    the expansion of SEED, in reverse order.  This decimal value is actually
//    still in base BASE, so it must be properly interpreted to generate
//    a usable value.
//
//  Example:
//
//    BASE = 2
//
//    SEED     SEED      van der Corput
//    decimal  binary    binary   decimal
//    -------  ------    ------   -------
//        0  =     0  =>  .0     = 0.0
//        1  =     1  =>  .1     = 0.5
//        2  =    10  =>  .01    = 0.25
//        3  =    11  =>  .11    = 0.75
//        4  =   100  =>  .001   = 0.125
//        5  =   101  =>  .101   = 0.625
//        6  =   110  =>  .011   = 0.375
//        7  =   111  =>  .111   = 0.875
//        8  =  1000  =>  .0001  = 0.0625
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
//  Reference:
//
//    John Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, pages 84-90, 1960.
//
//    John Hammersley,
//    Monte Carlo methods for solving multivariable problems,
//    Proceedings of the New York Academy of Science,
//    Volume 86, pages 844-874, 1960.
//
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Input, int SEED, the index of the desired element.
//    SEED should be nonnegative.
//    SEED = 0 is allowed, and returns R = 0.
//
//    Input, int BASE, the van der Corput base, which is usually
//    a prime number.  BASE must be greater than 1.
//
//    Output, double VAN_DER_CORPUT, the SEED-th element of the van
//    der Corput sequence for base BASE.
//
{
  double base_inv;
  int digit;
  double r;

  if ( base <= 1 )
  {
    cout << "\n";
    cout << "I4_TO_VAN_DER_CORPUT - Fatal error!\n";
    cout << "  The input base BASE is <= 1!\n";
    cout << "  BASE = " << base << "\n";
    exit ( 1 );
  }

  if ( seed < 0 )
  {
    cout << "\n";
    cout << "I4_TO_VAN_DER_CORPUT - Fatal error!\n";
    cout << "  SEED < 0.";
    cout << "  SEED = " << seed << "\n";
    exit ( 1 );
  }

  r = 0.0;

  base_inv = 1.0 / ( ( double ) base );

  while ( seed != 0 )
  {
    digit = seed % base;
    r = r + ( ( double ) digit ) * base_inv;
    base_inv = base_inv / ( ( double ) base );
    seed = seed / base;
  }

  return r;
}
//****************************************************************************80

void i4_to_van_der_corput_sequence ( int seed, int base, int n, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_VAN_DER_CORPUT_SEQUENCE: next N elements of a van der Corput sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2003
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
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Input, int SEED, the index of the first desired element.
//    SEED should be nonnegative.
//    SEED = 0 is allowed, and returns R = 0.
//
//    Input, int BASE, the van der Corput base, which is typically a
//    prime number.
//
//    Input, int N, the number of elements desired.
//
//    Output, double R[N], the SEED-th through
//    (SEED+N-1)-th elements of the van der Corput sequence for
//    the given base.
//
{
  double base_inv;
  int digit;
  int i;
  int seed2;

  if ( base <= 1 )
  {
    cout << "\n";
    cout << "I4_TO_VAN_DER_CORPUT_SEQUENCE - Fatal error!\n";
    cout << "  The input base BASE is <= 1!\n";
    cout << "  BASE = " << base << "\n";
    exit ( 1 );
  }

  if ( seed < 0 )
  {
    cout << "\n";
    cout << "I4_TO_VAN_DER_CORPUT_SEQUENCE - Fatal error!\n";
    cout << "  The input base SEED is < 0!\n";
    cout << "  SEED = " << seed << "\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    r[i] = 0.0;
    seed2 = seed + i;
    base_inv = 1.0E+00 / ( ( double ) base );

    while ( seed2 != 0 )
    {
      digit = seed2 % base;
      r[i] = r[i] + ( ( double ) digit ) * base_inv;
      base_inv = base_inv / ( ( double ) base );
      seed2 = seed2 / base;
    }
  }

  return;
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
//    01 September 2012
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
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2003
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
//****************************************************************************80

double van_der_corput ( )

//****************************************************************************80
//
//  Purpose:
//
//    VAN_DER_CORPUT computes the next element in the van der Corput sequence.
//
//  Discussion:
//
//    The internal variables van_der_corput_SEED and van_der_corput_BASE
//    control the van der Corput sequence.
//
//    The value of van_der_corput_SEED is incremented by one on each call.
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Output, double VAN_DER_CORPUT, the next element of the van der Corput sequence.
//
{
  double r;

  r = i4_to_van_der_corput ( van_der_corput_SEED, van_der_corput_BASE );

  van_der_corput_SEED = van_der_corput_SEED + 1;

  return r;
}
//****************************************************************************80

int van_der_corput_base_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    VAN_DER_CORPUT_BASE_GET gets the base for a van der Corput sequence.
//
//  Modified:
//
//    28 August 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Output, int VAN_DER_CORPUT_BASE_GET, the base for the
//    van der Corput sequence.
//
{
  return van_der_corput_BASE;
}
//****************************************************************************80

void van_der_corput_base_set ( int base )

//****************************************************************************80
//
//  Purpose:
//
//    VAN_DER_CORPUT_BASE_SET sets the base for a van der Corput sequence.
//
//
//  Modified:
//
//    28 August 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Input, int BASE, the base for the van der Corput sequence.
//    BASE must be greater than 1.
//
{
  if ( base <= 1 )
  {
    cout << "\n";
    cout << "VAN_DER_CORPUT_BASE_SET - Fatal error!\n";
    cout << "  The input base BASE is <= 1!\n";
    cout << "  BASE = " << base << "\n";
    exit ( 1 );
  }

  van_der_corput_BASE = base;

  return;
}
//****************************************************************************80

int van_der_corput_seed_get ( )

//****************************************************************************80
//
//  Purpose:
//
//    VAN_DER_CORPUT_SEED_GET gets the "seed" for the van der Corput sequence.
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Output, int VAN_DER_CORPUT_SEED_GET, the current seed for
//    the van der Corput sequence.
//
{
  return van_der_corput_SEED;
}
//****************************************************************************80

void van_der_corput_seed_set ( int seed )

//****************************************************************************80
//
//  Purpose:
//
//    VAN_DER_CORPUT_SEED_SET sets the "seed" for the van der Corput sequence.
//
//  Discussion:
//
//    Calling VAN_DER_CORPUT repeatedly returns the elements of the
//    van der Corput sequence in order, starting with element number 1.
//    An internal counter, called SEED, keeps track of the next element
//    to return.  Each time the routine is called, the SEED-th element
//    is computed, and then SEED is incremented by 1.
//
//    To restart the van der Corput sequence, it is only necessary to reset
//    SEED to 1.  It might also be desirable to reset SEED to some other value.
//    This routine allows the user to specify any value of SEED.
//
//    The default value of SEED is 1, which restarts the van der Corput sequence.
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Input, int SEED, the seed for the van der Corput sequence.
//    SEED must be nonnegative.
//
{
  if ( seed < 0 )
  {
    cout << "\n";
    cout << "VAN_DER_CORPUT_SEED_SET - Fatal error!\n";
    cout << "  SEED < 0.";
    cout << "  SEED = " << seed << "\n";
    exit ( 1 );
  }

  van_der_corput_SEED = seed;

  return;
}
//****************************************************************************80

void van_der_corput_sequence ( int n, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    VAN_DER_CORPUT_SEQUENCE: next N elements in the van der Corput sequence.
//
//  Modified:
//
//    27 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Input, int N, the number of elements desired.
//
//    Output, double R[N], the next N elements of the van
//    der Corput sequence.
//
{
  int base;
  int seed;

  seed = van_der_corput_seed_get ( );

  base = van_der_corput_base_get ( );

  i4_to_van_der_corput_sequence ( seed, base, n, r );

  seed = seed + n;
  van_der_corput_seed_set ( seed );

  return;
}
//****************************************************************************80

int *vdc_numerator_sequence ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    VDC_NUMERATOR_SEQUENCE: van der Corput numerator sequence base 2.
//
//  Discussion:
//
//    The classical van der Corput sequence, base 2, can be considered
//    as a way of enumerating the dyadic fractions P/2^K in an order of
//    increasing denominator.
//
//    If we fix a value of K, then the first (2^K) - 1 items in the
//    sequence are fractions strictly between 0 and 1, which can be written
//    as P/2^K where 0 < P < 2^K.
//
//    This function determines the numerator sequence, that is, the values
//    P, which is interesting in its own right.  Note that the P sequence
//    is "nested" in the sense that if 2^(K-1) <= N1 < N2 < 2^(K), then the
//    sequence for N2 will begin with the sequence for N1.
//
//    The I-th value in the sequence can be determined by writing
//    the integer I in binary using K digits, and reversing the order.
//
//    N = 10
//
//    2^3 = 8 <= 10 < 16 = 2^4
//
//    1  0001  1000   8
//    2  0010  0100   4
//    3  0011  1100  12
//    4  0100  0010   2
//    5  0101  1010  10
//    6  0110  0110   6
//    7  0111  1110  14
//    8  1000  0001   1
//    9  1001  1001   9
//   10  1010  0101   5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 February 2011
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
//    Johannes van der Corput,
//    Verteilungsfunktionen I & II,
//    Nederl. Akad. Wetensch. Proc.,
//    Volume 38, 1935, pages 813-820, pages 1058-1066.
//
//  Parameters:
//
//    Input, int N, the number of elements to compute.
//
//    Output, int VDC_NUMERATOR_SEQUENCE[N], the elements of the van der Corput
//    numerator sequence base 2.
//
{
  int d;
  int i;
  int j;
  int n_log_2;
  int *p;
  int s;
//
//  Carry out the computation.
//
  p = new int[n];

  n_log_2 = i4_log_2 ( n ) + 1;

  for ( i = 0; i < n; i++ )
  {
    p[i] = 0;
    s = i + 1;
    for ( j = 0; j < n_log_2; j++ )
    {
      d = ( s % 2 );
      p[i] = 2 * p[i] + d;
      s = s / 2;
    }
  }

  return p;
}
