# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <complex>

using namespace std;

# include "uniform.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test065 ( );
void test07 ( );
void test08 ( );
void test09 ( );

void test10 ( );
void test11 ( );
void test111 ( );
void test112 ( );
void test118 ( );
void test119 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );

void test20 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for UNIFORM_PRB.
//
//  Discussion:
//
//    UNIFORM_PRB calls sample problems for the UNIFORM routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "UNIFORM_PRB\n";
  cout << "  C++ version:\n";
  cout << "  Test the UNIFORM library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test065 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test111 ( );
  test112 ( );
  test118 ( );
  test119 ( );
  test12  ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test14 ( );
  test19 ( );

  test20 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "UNIFORM_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests C4_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  int seed_init = 123456789;
  complex <float> value;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  C4_UNIFORM_01 computes pseudorandom complex values\n";
  cout << "  uniformly distributed in the unit circle.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = c4_uniform_01 ( seed );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << real ( value )
         << "  " << setw(12) << imag ( value ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests C4VEC_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> *cvec;
  int i;
  int n;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  C4VEC_UNIFORM_01_NEW computes pseudorandom complex values\n";
  cout << "  uniformly distributed in the unit circle.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  n = 10;
  cvec = c4vec_uniform_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << real ( cvec[i] )
         << "  " << setw(12) << imag ( cvec[i] ) << "\n";
  }

  delete [] cvec;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests C8_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  int seed_init = 123456789;
  complex <float> value;

  cout << "\n";
  cout << "TEST93\n";
  cout << "  C8_UNIFORM_01 computes pseudorandom C8 values\n";
  cout << "  uniformly distributed in the unit circle.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = c8_uniform_01 ( seed );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << real ( value )
         << "  " << setw(12) << imag ( value ) << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests C8VEC_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int seed;
  int seed_init = 123456789;
  complex <double> *zvec;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  C8VEC_UNIFORM_01_NEW computes pseudorandom C8 values\n";
  cout << "  uniformly distributed in the unit circle.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  n = 10;
  zvec = c8vec_uniform_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << real ( zvec[i] )
         << "  " << setw(12) << imag ( zvec[i] ) << "\n";
  }

  delete [] zvec;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests CH_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  char chi;
  char clo;
  int i;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  CH_UNIFORM_AB computes pseudorandom characters\n";
  cout << "  in an interval [CLO,CHI].\n";

  clo = 'A';
  chi = 'J';
  seed = seed_init;

  cout << "\n";
  cout << "  The lower endpoint CLO = '" << clo << "'.\n";
  cout << "  The upper endpoint CHI = '" << chi << "'.\n";
  wcout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  '" << ch_uniform_ab ( clo, chi, seed ) << "'\n";
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests GET_SEED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int oops;
  int oops_max = 10000;
  int seed;
  int seed_old;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  GET_SEED picks an initial seed value for R8_UNIFORM_01.\n";
  cout << "  The value chosen should vary over time, because\n";
  cout << "  the seed is based on reading the clock.\n";
  cout << "\n";
  cout << "  This is just the \"calendar\" clock, which does\n";
  cout << "  not change very fast, so calling GET_SEED several\n";
  cout << "  times in a row may result in the same value.\n";

  oops = 0;
  seed = 12345678;
  seed_old = seed;

  cout << "\n";
  cout << "  Initial seed is " << seed << "\n";
  cout << "\n";
  cout << "  Next 3 values of R8_UNIFORM_01:\n";
  cout << "\n";

  for ( j = 1; j <= 3; j++ )
  {
    cout << "  " << setw(10) << r8_uniform_01 ( seed ) << "\n";
  }

  for ( i = 1; i <= 4; i++ )
  {
    while ( 1 )
    {
      seed = get_seed ( );

      if ( seed != seed_old )
      {
        seed_old = seed;
        break;
      }
      oops = oops + 1;
      if ( oops_max < oops ) 
      {
        cout << "\n";
        cout << "  Oops!\n";
        cout << "  Same seed returned for " << oops_max 
             << " calls to GET_SEED!\n";
        cout << "  Could be a bad algorithm, slow clock, or fast machine!\n";
        cout << "  To avoid infinite loops, we take what we have now.\n";
        oops = 0;
        oops_max = oops_max * 10;
        break;
      }
    }

    cout << "\n";
    cout << "  New seed from GET_SEED is " << seed << "\n";
    cout << "\n";
    cout << "  Next 3 values of R8_UNIFORM_01:\n";
    cout << "\n";

    for ( j = 1; j <= 3; j++ )
    {
      cout << "  " << setw(10) << r8_uniform_01 ( seed ) << "\n";
    }

  }
  return;
}
//****************************************************************************80

void test065 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST065 tests I4_SEED_ADVANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int seed;
  int seed_new;
  int step;

  cout << "\n";
  cout << "TEST065\n";
  cout << "  I4_SEED_ADVANCE advances the seed.\n";
  cout << "\n";
  cout << "  Step        SEED input       SEED output\n";
  cout << "\n";

  seed_new = 12345;

  for ( step = 1; step <= 10; step++)
  {
    seed = seed_new;
    seed_new = i4_seed_advance ( seed );

    cout << "  " << setw(4)  << step
         << "  " << setw(16) << seed
         << "  " << setw(16) << seed_new << "\n";
  }

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests I4_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
# define A 6
# define B 10

  int a = A;
  int b = B;
  int freq[B+1-A];
  int i;
  int j;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  I4_UNIFORM_AB computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = a; i <= b; i++ )
  {
    freq[i-a] = 0;
  }

  for ( i = 1; i <= 10000; i++ )
  {
    j = i4_uniform_ab ( a, b, seed );
    if ( j < a ) 
    {
      cout << "  Illegal value J = " << j << "\n";
    }
    else if ( j <= b )
    {
      freq[j-a] = freq[j-a] + 1;
    }
    else
    {
      cout << "  Illegal value J = " << j << "\n";
    }
  }

  cout << "\n";
  cout << "         I    Frequency\n";
  cout << "\n";
  for ( i = a; i <= b; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << freq[i-a] << "\n";
  }

  return;
# undef A
# undef B
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests I4_UNIFORM_AB.
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
{
  int a = -100;
  int b = 200;
  int i;
  int j;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  I4_UNIFORM_AB computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 20; i++ )
  {
    j = i4_uniform_ab ( a, b, seed );

    cout << "  " << setw(8) << i
         << "  " << setw(8) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests I4_UNIFORM_0I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 1000

  int i;
  double mean;
  int seed;
  double variance;
  int x[N];

  cout << "\n";
  cout << "TEST09\n";
  cout << "  I4_UNIFORM_0I samples a uniform random\n";
  cout << "  integer distribution in [0,2**31-1].\n";

  seed = 123456789;

  cout << "\n";
  cout << "  Starting with seed = " << seed << "\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = i4_uniform_0i ( seed );
  }

  cout << "\n";
  cout << "  First few values:\n";
  cout << "\n";
  for ( i = 0; i < 5; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << x[i] << "\n";
  }

  mean = i4vec_mean ( N, x );

  variance = i4vec_variance ( N, x );

  cout << "\n";
  cout << "  Number of values computed was N = " << N << "\n";
  cout << "  Average value was " << mean << "\n";
  cout << "  Minimum value was " << i4vec_min ( N, x ) << "\n";
  cout << "  Maximum value was " << i4vec_max ( N, x ) << "\n";
  cout << "  Variance was " << variance << "\n";

  return;
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests I4VEC_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 November 2006
//
//  Author:
//
//    John Burkardt
//
{
# define A 6
# define B 10
# define N 10000

  int a = A;
  int b = B;
  int freq[B+1-A];
  int i;
  int *i4vec;
  int j;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  I4VEC_UNIFORM_AB_NEW computes a vector of pseudorandom values\n";
  cout << "  in an interval [A,B].\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = a; i <= b; i++ )
  {
    freq[i-a] = 0;
  }

  i4vec = i4vec_uniform_ab_new ( N, a, b, seed );

  for ( i = 0; i < N; i++ )
  {
    j = i4vec[i];
    if ( j < a ) 
    {
      cout << "  Illegal value J = " << j << "\n";
    }
    else if ( j <= b )
    {
      freq[j-a] = freq[j-a] + 1;
    }
    else
    {
      cout << "  Illegal value J = " << j << "\n";
    }
  }

  cout << "\n";
  cout << "         I    Frequency\n";
  cout << "\n";
  for ( i = a; i <= b; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << freq[i-a] << "\n";
  }

  delete [] i4vec;

  return;
# undef A
# undef B
# undef N
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests I8_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2007
//
//  Author:
//
//    John Burkardt
//
{

  long long int a;
  long long int b;
  int i;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  I8_UNIFORM_AB computes pseudorandom values \n";
  cout << "  in an interval [A,B].\n";

//a = 100000;
  a = 1000000000LL;
//b = 800000;
  b = 8000000000LL;

  seed = seed_init;

  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is    " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(24) << i8_uniform_ab ( a, b, seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test111 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST111 tests L_UNIFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST111\n";
  cout << "  L_UNIFORM computes pseudorandom logical values.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
   cout << "  " << setw(8) << i
        << "  " << setw(1) <<  l_uniform ( seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test112 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST112 tests LVEC_UNIFORM_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  bool *lvec;
  int n = 10;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST112\n";
  cout << "  LVEC_UNIFORM_NEW computes a vector of\n";
  cout << "  pseudorandom logical values.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";

  lvec = lvec_uniform_new ( n, seed );

  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8) << i + 1
         << "  " << setw(1) << lvec[i] << "\n";
  }

  delete [] lvec;

  return;
}
//****************************************************************************80

void test118 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST118 tests LCRG_ANBN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int an;
  int b;
  int bn;
  int c;
  int n;

  cout << "\n";
  cout << "TEST118\n";
  cout << "  LCRG_ANBN determines a linear congruential random\n";
  cout << "  number generator equivalent to N steps of a given one.\n";
//
//  These parameters define the old (1969) IBM 360 random number generator:
//
  a = 16807;
  b = 0;
  c = 2147483647;

  cout << "\n";
  cout << "  LCRG parameters:\n";
  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "  C = " << c << "\n";
  cout << "\n";
  cout << "             N             A             B\n";
  cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    lcrg_anbn ( a, b, c, n, &an, &bn );
    cout << "  " << setw(12) << n
         << "  " << setw(12) << an
         << "  " << setw(12) << bn << "\n";
  }

  return;
}
//****************************************************************************80

void test119 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST119 tests LCRG_ANBN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int an;
  int b;
  int bn;
  int c;
  int j;
  int k;
  int n;
  int u;
  int v;
  int *x;
  int *y;

  cout << "\n";
  cout << "TEST119\n";
  cout << "  LCRG_ANBN determines a linear congruential random\n";
  cout << "  number generator equivalent to N steps of a given one.\n";
//
//  These parameters define the old (1969) IBM 360 random number generator:
//
  a = 16807;
  b = 0;
  c = 2147483647;

  cout << "\n";
  cout << " LCRG parameters:\n";
  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "  C = " << c << "\n";
  cout << "\n";
  cout << "                          N            In           Out\n";
  cout << "\n";

  k = 0;
  u = 12345;
  cout << "  " << "            "
       << "  " << setw(12) << k
       << "  " << "            "
       << "  " << setw(12) << u << "\n";
  for ( k = 1; k <= 11; k++ )
  {
    v = lcrg_evaluate ( a, b, c, u );
  cout << "  " << "            "
       << "  " << setw(12) << k
       << "  " << setw(12) << u 
       << "  " << setw(12) << v << "\n";
    u = v;
  }
//
//  Now try to replicate these results using N procesors.
//
  n = 4;
  x = new int[n];
  y = new int[n];

  lcrg_anbn ( a, b, c, n, &an, &bn );

  cout << "\n";
  cout << "  LCRG parameters:\n";
  cout << "\n";
  cout << "  AN = " << an << "\n";
  cout << "  BN = " << bn << "\n";
  cout << "  C  = " << c << "\n";
  cout << "\n";
  cout << "             J             N            In           Out\n";
  cout << "\n";

  x[0] = 12345;
  for ( j = 1; j < n; j++ )
  {
    x[j] = lcrg_evaluate ( a, b, c, x[j-1] );
  }

  for ( j = 0; j < n; j++ )
  {
    cout << "  " << setw(12) << j + 1
         << "  " << setw(12) << j
         << "  " << "            "
         << "  " << setw(12) << x[j] << "\n";
  }

  for ( k = n + 1; k <= 12; k = k + n )
  {
    for ( j = 0; j < n; j++ )
    {
      y[j] = lcrg_evaluate ( an, bn, c, x[j] );
      cout << "  " << setw(12) << j + 1
           << "  " << setw(12) << k + j - 1 
           << "  " << setw(12) << x[j]
           << "  " << setw(12) << y[j] << "\n";
      x[j] = y[j];
    }
  }

  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests LCRG_SEED
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  int c;
  int i;
  int seed;
  int seed_in;
  int seed_lcrg;
  int seed_out;
  int seed_start;
  float u;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  LCRG_SEED directly computes the updated value of a\n";
  cout << "  seed used by an linear congruential random number\n";
  cout << "  generator.\n";
  cout << "\n";
  cout << "       I          SEED          SEED          SEED    U\n";
  cout << "                 Input        Output          LCRG\n";
  cout << "\n";
//
//  These parameters define the old (1969) IBM 360 random number generator:
//
  a = 16807;
  b = 0;
  c = 2147483647;
//
//  This seed value was used in Pierre L'Ecuyer's article.
//
  seed_start = 12345;

  seed = seed_start;
//
//  Compute 1000 random numbers "the hard way", that is, sequentially.
//  Every now and then, call LCRG_SEED to compute SEED directly.
//
  for ( i = 1; i <= 1000; i++ )
  {
    seed_in = seed;
    u = r4_uniform_01 ( seed );
    seed_out = seed;

    if ( i <= 10 || i == 100 || i == 1000 )
    {
      seed_lcrg = lcrg_seed ( a, b, c, i, seed_start );

      cout << "  "
           << setw(6)  << i         << "  "
           << setw(12) << seed_in   << "  "
           << setw(12) << seed_out  << "  "
           << setw(12) << seed_lcrg << "  "
           << setw(14) << u         << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests R4_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  float a;
  float b;
  int i;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  R4_UNIFORM_AB computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";

  a = 5.0;
  b = 10.0;
  seed = seed_init;

  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r4_uniform_ab ( a, b, seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests R4_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  R4_UNIFORM_01 computes pseudorandom values \n";
  cout << "  in the interval [0,1].\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << r4_uniform_01 ( seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests R8_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  R8_UNIFORM_AB computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";

  a = 5.0;
  b = 10.0;
  seed = seed_init;

  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r8_uniform_ab ( a, b, seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests R8_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST16\n";
  cout << "  R8_UNIFORM_01 computes pseudorandom values \n";
  cout << "  in the interval [0,1].\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << r8_uniform_01 ( seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests R8_UNIFORM_01;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 1000

  int i;
  int seed;
  int seed_in;
  int seed_out;
  double u[N];
  double u_avg;
  double u_var;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  R8_UNIFORM_01 computes a sequence of uniformly distributed\n";
  cout << "  pseudorandom numbers.\n";

  seed = 12345;

  cout << "\n";
  cout << "  Initial SEED = " << seed << "\n";

  cout << "\n";
  cout << "  First 10 values:\n";
  cout << "\n";
  cout << "       I         Input        Output   R8_UNIFORM_01\n";
  cout << "                  SEED          SEED\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    seed_in = seed;
    u[i] = r8_uniform_01 ( seed );
    seed_out = seed;

    cout                         << "  "
         << setw(6)  << i + 1    << "  " 
         << setw(12) << seed_in  << "  " 
         << setw(12) << seed_out << "  " 
         << setw(10) << u[i]     << "\n";
  }

  cout << "\n";
  cout << "  Now call R8_UNIFORM_01 " << N << " times.\n";

  u_avg = 0.0;
  for ( i = 0; i < N; i++ )
  {
    u[i] = r8_uniform_01 ( seed );
    u_avg = u_avg + u[i];
  }

  u_avg = u_avg / ( ( double ) N );

  u_var = 0.0;
  for ( i = 0; i < N; i++ )
  {
    u_var = u_var + ( u[i] - u_avg ) * ( u[i] - u_avg );
  }
  u_var = u_var / ( ( double ) ( N - 1 ) );

  cout << "\n";
  cout << "  Average value = " << u_avg      << "\n";
  cout << "  Expecting       " << 0.5        << "\n";

  cout << "\n";
  cout << "  Variance =      " << u_var      << "\n";
  cout << "  Expecting       " << 1.0 / 12.0 << "\n";

  return;
# undef N
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests R8_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  int seed_in;
  int seed_out;
  int seed_save;
  double x;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  R8_UNIFORM_01 computes a sequence of pseudorandom numbers\n";
  cout << "  but all computations depend on the seed value.\n";
  cout << "  In this test, we show how a sequence of \"random\"\n";
  cout << "  values can be manipulated by accessing the seed.\n";

  seed = 1066;

  cout << "\n";
  cout << "  Set SEED to " << seed << "\n";
  cout << "\n";
  cout << "  Now call R8_UNIFORM_01 10 times, and watch SEED.\n";
  cout << "\n";
  cout << "       I         Input        Output   R8_UNIFORM_01\n";
  cout << "                  SEED          SEED\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;

    if ( i == 5 )
    {
      seed_save = seed;
    }
    x = r8_uniform_01 ( seed );
    seed_out = seed;
    cout                            << "  "
         << setw ( 6 )  << i        << "  "
         << setw ( 12 ) << seed_in  << "  "
         << setw ( 12 ) << seed_out << "  "
         << setw ( 10 ) << x        << "\n";

  }

  seed = seed_save;

  cout << "\n";
  cout << "  Reset SEED to its value at step 5, = " << seed << "\n";
  cout << "\n";
  cout << "  Now call R8_UNIFORM_01 10 times, and watch how SEED\n";
  cout << "  and R8_UNIFORM_01 restart themselves.\n";
  cout << "\n";
  cout << "       I         Input        Output   R8_UNIFORM_01\n";
  cout << "                  SEED          SEED\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = r8_uniform_01 ( seed );
    seed_out = seed;
    cout                            << "  "
         << setw ( 6 )  << i        << "  "
         << setw ( 12 ) << seed_in  << "  "
         << setw ( 12 ) << seed_out << "  "
         << setw ( 10 ) << x        << "\n";
  }

  seed = -12345678;

  cout << "\n";
  cout << "  What happens with an initial negative SEED?\n";
  cout << "\n";
  cout << "       I         Input        Output   R8_UNIFORM_01\n";
  cout << "                  SEED          SEED\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = r8_uniform_01 ( seed );
    seed_out = seed;
    cout                            << "  "
         << setw ( 6 )  << i        << "  "
         << setw ( 12 ) << seed_in  << "  "
         << setw ( 12 ) << seed_out << "  "
         << setw ( 10 ) << x        << "\n";
  }

  return;
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests R8_UNIFORM_01 and R8MAT_UNIFORM_01_NEW.
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
{
# define M 100
# define N 10

  double a[M*N];
  double *b;
  int i;
  int j;
  int k;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  R8_UNIFORM_01 computes pseudorandom values one at a time.\n";
  cout << "  R8MAT_UNIFORM_01_NEW computes a matrix of values.\n";
  cout << "\n";
  cout << "  For the same initial seed, the results should be identical,\n";
  cout << "  but R8MAT_UNIFORM_01_NEW might be faster.\n";
  cout << "\n";
  cout << "  Initial seed is " << seed_init << "\n";

  seed = seed_init;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < M; i++ )
    {
      a[i+j*M] = r8_uniform_01 ( seed );
    }
  }

  seed = seed_init;
  b = r8mat_uniform_01_new ( M, N, seed );

  cout << "\n";
  cout << "      I       J      A[I,J]        B[I,J]\n";
  cout << "                 (R8_UNIFORM_01)  (R8MAT_UNIFORM_01_NEW)\n";
  cout << "\n";

  for ( k = 0; k < 11; k++ )
  {
    i = ( k * ( M - 1 ) ) / 10;
    j = ( k * ( N - 1 ) ) / 10;

    cout << " "
         << setw(6) << i         << "  "
         << setw(6) << j         << "  "
         << setw(12) << a[i+j*M] << "  "
         << setw(12) << b[i+j*M] << "\n";
  }
  
  delete [] b;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests R8_UNIFORM_01 and R8VEC_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a[N];
  double *b;
  int i;
  int j;
  int k;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  R8_UNIFORM_01 computes pseudeorandom values one at a time.\n";
  cout << "  R8VEC_UNIFORM_01_NEW computes a vector of values.\n";
  cout << "\n";
  cout << "  For the same initial seed, the results should be identical,\n";
  cout << "  but R8VEC_UNIFORM_01 might be faster.\n";
  cout << "\n";
  cout << "  Initial seed is " << seed_init << "\n";

  seed = seed_init;
  for ( i = 0; i < N; i++ )
  {
    a[i] = r8_uniform_01 ( seed );
  }

  seed = seed_init;
  b = r8vec_uniform_01_new ( N, seed );

  cout << "\n";
  cout << "      I      A[I]          B[I]\n";
  cout << "         (R8_UNIFORM_01)  (R8VEC_UNIFORM_01)\n";
  cout << "\n";

  for ( i = 1; i < N; i++ )
  {
    cout << " "
         << setw(6) << i         << "  "
         << setw(12) << a[i] << "  "
         << setw(12) << b[i] << "\n";
  }
  
  delete [] b;

  return;
# undef N
}
