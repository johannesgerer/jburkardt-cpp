# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

# include "mpi.h"

int main ( int argc, char *argv[] );
int congruence ( int a, int b, int c, int *error );
int i4_gcd ( int i, int j );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_sign ( int i );
void lcrg_anbn ( int a, int b, int c, int n, int *an, int *bn );
int lcrg_evaluate ( int a, int b, int c, int x );
int power_mod ( int a, int n, int m );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for RANDOM_MPI.
//
//  Discussion:
//
//    This program demonstrates how P processors can generate the same
//    sequence of random numbers as 1 processor.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: Portable Parallel Programming with the
//    Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323,
//    LC: QA76.642.G76.
//
{
  int a;
  int an;
  int b;
  int bn;
  int c;
  int error;
  int id;
  int j;
  int k;
  int k_hi;
  int p;
  int u;
  int v;
//
//  Initialize MPI.
//
  MPI::Init ( argc, argv );
//
//  Get the number of processors.
//
  p = MPI::COMM_WORLD.Get_size ( );
//
//  Get the rank of this processor.
//
  id = MPI::COMM_WORLD.Get_rank ( );
//
//  Print a message.
//
  if ( id == 0 )
  {
    timestamp ( );
    cout << "\n";
    cout << "RANDOM_MPI - Master process:\n";
    cout << "  C++ version\n";
    cout << "  The number of processors is P = " << p << "\n";
    cout << "\n";
    cout << "  This program shows how a stream of random numbers\n";
    cout << "  can be computed 'in parallel' in an MPI program.\n";
    cout << "\n";
    cout << "  We assume we are using a linear congruential\n";
    cout << "  random number generator or LCRG, which takes\n";
    cout << "  an integer input and returns a new integer output:\n";
    cout << "\n";
    cout << "    U = ( A * V + B ) mod C\n";
    cout << "\n";
    cout << "  We assume that we want the MPI program to produce\n";
    cout << "  the same sequence of random values as a sequential\n";
    cout << "  program would - but we want each processor to compute\n";
    cout << "  one part of that sequence.\n";
    cout << "\n";
    cout << "  We do this by computing a new LCRG which can compute\n";
    cout << "  every P'th entry of the original one.\n";
    cout << "\n";
    cout << "  Our LCRG works with integers, but it is easy to\n";
    cout << "  turn each integer into a real number between [0,1].\n";
  }
//
//  A, B and C define the linear congruential random number generator.
//
  a = 16807;
  b = 0;
  c = 2147483647;

  if ( id == 0 )
  {
    cout << "\n";
    cout << "  LCRG parameters:\n";
    cout << "\n";
    cout << "  A  = " << a << "\n";
    cout << "  B  = " << b << "\n";
    cout << "  C  = " << c << "\n";
  }

  k_hi = p * 10;
//
//  Processor 0 generates 10 * P random values.
//
  if ( id == 0 )
  {
    cout << "\n";
    cout << "  Let processor 0 generate the entire random number sequence.\n";
    cout << "\n";
    cout << "     K    ID         Input        Output\n";
    cout << "\n";

    k = 0;
    v = 12345;
    cout << "  " << setw(4) << k
         << "  " << setw(4) << id
         << "  " << "            "
         << "  " << setw(12) << v << "\n";

    for ( k = 1; k <= k_hi; k++ )
    {
      u = v;
      v = lcrg_evaluate ( a, b, c, u );
      cout << "  " << setw(4) << k
           << "  " << setw(4) << id
           << "  " << setw(12) << u
           << "  " << setw(12) << v << "\n";
    }
  }
//
//  Processor P now participates by computing the P-th part of the sequence.
//
  lcrg_anbn ( a, b, c, p, &an, &bn );

  if ( id == 0 )
  {
    cout << "\n";
    cout << "  LCRG parameters for P processors:\n";
    cout << "\n";
    cout << "  AN = " << an << "\n";
    cout << "  BN = " << bn << "\n";
    cout << "  C  = " << c << "\n";
    cout << "\n";
    cout << "  Have ALL the processors participate in computing\n";
    cout << "  the same random number sequence.\n";
    cout << "\n";
    cout << "     K    ID         Input        Output\n";
    cout << "\n";
  }
//
//  On System X, just to try to keep the output of the various processors 
//  sorted, I had to resort to a SLEEP command to keep the processors from
//  competing.  Obviously, this is only to make the output readable
//  so we can verify it.  It's NOT what you want to do in a production
//  situation!
//
//  sleep ( id * 2 );
//
//  Use the basis LCRG to get the ID-th value in the sequence.
//
  v = 12345;
  for ( j = 1; j <= id; j++ )
  {
    u = v;
    v = lcrg_evaluate ( a, b, c, u );
  }
  k = id;

  cout << "  " << setw(4) << k
       << "  " << setw(4) << id
       << "  " << "            "
       << "  " << setw(12) << v << "\n";
//
//  Now use the "skipping" LCRG to compute the values with indices
//  ID, ID+P, ID+2P, ...,
//
  for ( k = id + p; k <= k_hi; k = k + p )
  {
    u = v;
    v = lcrg_evaluate ( an, bn, c, u );
    cout << "  " << setw(4) << k
         << "  " << setw(4) << id
         << "  " << setw(12) << u
         << "  " << setw(12) << v << "\n";
  }
//
//  Terminate MPI.
//
  MPI::Finalize ( );
//
//  Terminate.
//
  if ( id == 0 )
  {
    cout << "\n";
    cout << "RANDOM_MPI:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }
  return 0;
}
//****************************************************************************80

int congruence ( int a, int b, int c, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    CONGRUENCE solves a congruence of the form ( A * X = C ) mod B.
//
//  Discussion:
//
//    A, B and C are given integers.  The equation is solvable if and only
//    if the greatest common divisor of A and B also divides C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein, editor,
//    CRC Concise Encylopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input, int A, B, C, the coefficients of the Diophantine equation.
//
//    Output, bool *ERROR, error flag, is TRUE if an error occurred..
//
//    Output, int CONGRUENCE, the solution of the Diophantine equation.
//    X will be between 0 and B-1.
//
{
# define N_MAX 100

  int a_copy;
  int a_mag;
  int a_sign;
  int b_copy;
  int b_mag;
  int b_sign;
  int c_copy;
  int g;
  int k;
  int n;
  float norm_new;
  float norm_old;
  int q[N_MAX];
  bool swap;
  int temp;
  int x;
  int xnew;
  int y;
  int ynew;
  int z;
//
//  Defaults for output parameters.
//
  *error = false;
  x = 0;
  y = 0;
//
//  Special cases.
//
  if ( a == 0 && b == 0 && c == 0 )
  {
    x = 0;
    return x;
  }
  else if ( a == 0 && b == 0 && c != 0 )
  {
    *error = true;
    x = 0;
    return x;
  }
  else if ( a == 0 && b != 0 && c == 0 )
  {
    x = 0;
    return x;
  }
  else if ( a == 0 && b != 0 && c != 0 )
  {
    x = 0;
    if ( ( c % b ) != 0 )
    {
      *error = true;
    }
    return x;
  }
  else if ( a != 0 && b == 0 && c == 0 )
  {
    x = 0;
    return x;
  }
  else if ( a != 0 && b == 0 && c != 0 )
  {
    x = c / a;
    if ( ( c % a ) != 0 )
    {
      *error = true;
    }
    return x;
  }
  else if ( a != 0 && b != 0 && c == 0 )
  {
//  g = i4_gcd ( a, b );
//  x = b / g;
    x = 0;
    return x;
  }
//
//  Now handle the "general" case: A, B and C are nonzero.
//
//  Step 1: Compute the GCD of A and B, which must also divide C.
//
  g = i4_gcd ( a, b );

  if ( ( c % g ) != 0 )
  {
    *error = true;
    return x;
  }

  a_copy = a / g;
  b_copy = b / g;
  c_copy = c / g;
//
//  Step 2: Split A and B into sign and magnitude.
//
  a_mag = abs ( a_copy );
  a_sign = i4_sign ( a_copy );
  b_mag = abs ( b_copy );
  b_sign = i4_sign ( b_copy );
//
//  Another special case, A_MAG = 1 or B_MAG = 1.
//
  if ( a_mag == 1 )
  {
    x = a_sign * c_copy;
    return x;
  }
  else if ( b_mag == 1 )
  {
    x = 0;
    return x;
  }
//
//  Step 3: Produce the Euclidean remainder sequence.
//
  if ( b_mag <= a_mag )
  {
    swap = false;
    q[0] = a_mag;
    q[1] = b_mag;
  }
  else
  {
    swap = true;
    q[0] = b_mag;
    q[1] = a_mag;
  }

  n = 3;

  for ( ; ; )
  {
    q[n-1] = ( q[n-3] % q[n-2] );

    if ( q[n-1] == 1 )
    {
      break;
    }

    n = n + 1;

    if ( N_MAX < n )
    {
      *error = true;
      cout << "\n";
      cout << "CONGRUENCE - Fatal error!\n";
      cout << "  Exceeded number of iterations.\n";
      exit ( 1 );
    }
  }
//
//  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
//
  y = 0;
  for ( k = n; 2 <= k; k-- )
  {
    x = y;
    y = ( 1 - x * q[k-2] ) / q[k-1];
  }
//
//  Step 5: Undo the swapping.
//
  if ( swap )
  {
    z = x;
    x = y;
    y = z;
  }
//
//  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
//
  x = x * a_sign;
//
//  Step 7: Multiply by C, so that X * A + Y * B = C.
//
  x = x * c_copy;
//
//  Step 8: Now force 0 <= X < B.
//
  x = x % b;
//
//  Step 9: Force positivity.
//
  if ( x < 0 )
  {
    x = x + b;
  }

  return x;
# undef N_MAX
}
//****************************************************************************80

int i4_gcd ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_GCD finds the greatest common divisor of I and J.
//
//  Discussion:
//
//    Only the absolute values of I and J are considered, so that the 
//    result is always nonnegative.
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
//    Input, int I1, I2, two integers to be compared.
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

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
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

void lcrg_anbn ( int a, int b, int c, int n, int *an, int *bn )

//****************************************************************************80
//
//  Purpose:
//
//    LCRG_ANBN computes the "N-th power" of a linear congruential generator.
//
//  Discussion:
//
//    We are considering a linear congruential random number generator.
//    The LCRG takes as input an integer value called SEED, and returns
//    an updated value of SEED,
//
//      SEED(out) = ( a * SEED(in) + b ) mod c.
//
//    and an associated pseudorandom real value
//
//      U = SEED(out) / c.
//
//    In most cases, a user is content to call the LCRG repeatedly, with
//    the updating of SEED being taken care of automatically.
//
//    The purpose of this routine is to determine the values of AN and BN
//    that describe the LCRG that is equivalent to N applications of the
//    original LCRG.
//
//    One use for such a facility would be to do random number computations
//    in parallel.  If each of N processors is to compute many random values,
//    you can guarantee that they work with distinct random values
//    by starting with a single value of SEED, using the original LCRG to generate
//    the first N-1 "iterates" of SEED, so that you now have N "seed" values,
//    and from now on, applying the N-th power of the LCRG to the seeds.
//
//    If the K-th processor starts from the K-th seed, it will essentially
//    be computing every N-th entry of the original random number sequence,
//    offset by K.  Thus the individual processors will be using a random
//    number stream as good as the original one, and without repeating, and
//    without having to communicate.
//
//    To evaluate the N-th value of SEED directly, we start by ignoring
//    the modular arithmetic, and working out the sequence of calculations
//    as follows:
//
//      SEED(0)   =     SEED.
//      SEED(1)   = a * SEED      + b
//      SEED(2)   = a * SEED(1)   + b = a^2 * SEED           + a * b + b
//      SEED(3)   = a * SEED(2)   + b = a^3 * SEED + a^2 * b + a * b + b
//      ...
//      SEED(N-1) = a * SEED(N-2) + b
//
//      SEED(N) = a * SEED(N-1) + b = a^N * SEED
//                                    + ( a^(n-1) + a^(n-2) + ... + a + 1 ) * b
//
//    or, using the geometric series,
//
//      SEED(N) = a^N * SEED + ( a^N - 1) / ( a - 1 ) * b
//              = AN * SEED + BN
//
//    Thus, from any SEED, we can determine the result of N applications of the
//    original LCRG directly if we can solve
//
//      ( a - 1 ) * BN = ( a^N - 1 ) * b in modular arithmetic,
//
//    and evaluate:
//
//      AN = a^N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Barry Wilkinson, Michael Allen,
//    Parallel Programming:
//    Techniques and Applications Using Networked Workstations and Parallel Computers,
//    Prentice Hall,
//    ISBN: 0-13-140563-2,
//    LC: QA76.642.W54.
//
//  Parameters:
//
//    Input, int A, the multiplier for the LCRG.
//
//    Input, int  B, the added value for the LCRG.
//
//    Input, int  C, the base for the modular arithmetic.
//    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
//    required that 0 < C.
//
//    Input, int N, the "index", or number of times that the
//    LCRG is to be applied.  It is required that 0 <= N.
//
//    Output, int *AN, *BN, the multiplier and added value for
//    the LCRG that represent N applications of the original LCRG.
//
{
  int am1;
  int anm1tb;
  bool ierror;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "LCRG_ANBN - Fatal error!\n";
    cerr << "  Illegal input value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( c <= 0 )
  {
    cerr << "\n";
    cerr << "LCRG_ANBN - Fatal error!\n";
    cerr << "  Illegal input value of C = " << c << "\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    *an = 1;
    *bn = 0;
  }
  else if ( n == 1 )
  {
    *an = a;
    *bn = b;
  }
  else
  {
//
//  Compute A^N.
//
    *an = power_mod ( a, n, c );
//
//  Solve 
//    ( a - 1 ) * BN = ( a^N - 1 ) mod B
//  for BN.
//
    am1 = a - 1;
    anm1tb = ( *an - 1 ) * b;

    *bn = congruence ( am1, c, anm1tb, &ierror );

    if ( ierror )
    {
      cerr << "\n";
      cerr << "LCRG_ANBN - Fatal error!\n";
      cerr << "  An error occurred in the CONGRUENCE routine.\n";
      exit ( 1 );
    }
  }

  return;
}
//****************************************************************************80

int lcrg_evaluate ( int a, int b, int c, int x )

//****************************************************************************80
//
//  Purpose:
//
//    LCRG_EVALUATE evaluates an LCRG, y = ( A * x + B ) mod C.
//
//  Discussion:
//
//    This routine cannot be recommended for production use.  Because we want
//    to do modular arithmetic, but the base is not a power of 2, we need to
//    use "double precision" integers to keep accuracy.
//
//    If we knew the base C, we could try to avoid overflow while not changing
//    precision.
//
//    If the base C was a power of 2, we could rely on the usual properties of
//    integer arithmetic on computers, in which overflow bits, which are always
//    ignored, don ot actually matter.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int A, the multiplier for the LCRG.
//
//    Input, int B, the added value for the LCRG.
//
//    Input, int C, the base for the modular arithmetic.
//    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
//    required that 0 < C.
//
//    Input, int X, the value to be processed.
//
//    Output, int LCRG_EVALUATE, the processed value.
//
{
  long long int a8;
  long long int b8;
  long long int c8;
  long long int x8;
  int y;
  long long int y8;
//
//  To avoid roundoff issues, we need to go to "double precision" integers.
//  (Not available on all planets.)
//
  a8 = ( long long int ) a;
  b8 = ( long long int ) b;
  c8 = ( long long int ) c;
  x8 = ( long long int ) x;

  y8 = ( a8 * x8 + b8 ) % c8;

  y = ( int ) ( y8 );

  if ( y < 0 )
  {
    y = y + c;
  }

  return y;
}
//****************************************************************************80

int power_mod ( int a, int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    POWER_MOD computes ( A^N ) mod M.
//
//  Discussion:
//
//    Some programming tricks are used to speed up the computation, and to
//    allow computations in which A^N is much too large to store in a
//    real word.
//
//    First, for efficiency, the power A**N is computed by determining
//    the binary expansion of N, then computing A, A^2, A^4, and so on
//    by repeated squaring, and multiplying only those factors that
//    contribute to A^N.
//
//    Secondly, the intermediate products are immediately "mod'ed", which
//    keeps them small.
//
//    For instance, to compute mod ( A**13, 11 ), we essentially compute
//
//       13 = 1 + 4 + 8
//
//       A^13 = A * A^4 * A^8
//
//       A^13 mod ( 11 ) = A mod ( 11 ) * A^4 mod ( 11 ) * A^8 mod ( 11 ).
//
//    Fermat's little theorem says that if P is prime, and A is not divisible
//    by P, then ( A^(P-1) - 1 ) is divisible by P.
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
//    Input, int A, the base of the expression to be tested.
//    A should be nonnegative.
//
//    Input, int N, the power to which the base is raised.
//    N should be nonnegative.
//
//    Input, int M, the divisor against which the expression is tested.
//    M should be positive.
//
//    Output, int POWER_MOD, the remainder when A^N is divided by M.
//
{
  long long int a_square2;
  int d;
  long long int m2;
  int x;
  long long int x2;

  if ( a < 0 )
  {
    return -1;
  }

  if ( m <= 0 )
  {
    return -1;
  }

  if ( n < 0 )
  {
    return -1;
  }
//
//  A_SQUARE contains the successive squares of A.
//
  a_square2 = ( long long int ) a;
  m2 = ( long long int ) m;
  x2 = ( long long int ) 1;

  while ( 0 < n )
  {
    d = n % 2;

    if ( d == 1 )
    {
      x2 = ( x2 * a_square2 ) % m2;
    }

    a_square2 = ( a_square2 * a_square2 ) % m2;
    n = ( n - d ) / 2;
  }
//
//  Ensure that 0 <= X.
//
  while ( x2 < 0 )
  {
    x2 = x2 + m2;
  }

  x = ( int ) x2;

  return x;
}
//****************************************************************************80

void timestamp ( void )

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
//    02 October 2003
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
