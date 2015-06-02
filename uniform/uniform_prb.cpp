# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <complex>

using namespace std;

# include "uniform.hpp"

int main ( );
void bvec_uniform_new_test ( );
void c4_uniform_01_test ( );
void c4mat_uniform_01_new_test ( );
void c4vec_uniform_01_new_test ( );
void c8_uniform_01_test ( );
void c8mat_uniform_01_new_test ( );
void c8vec_uniform_01_new_test ( );
void ch_uniform_ab_test ( );
void get_seed_test ( );
void i4_seed_advance_test ( );
void i4_uniform_0i_test ( );
void i4_uniform_ab_test ( );
void i4mat_uniform_ab_new_test ( );
void i4vec_uniform_ab_new_test ( );
void l4_uniform_test ( );
void l4mat_uniform_new_test ( );
void l4vec_uniform_new_test ( );
void lcrg_anbn_test ( );
void lcrg_seed_test ( );
void r4_uniform_01_test ( );
void r4_uniform_ab_test ( );
void r4mat_uniform_ab_new_test ( );
void r4vec_uniform_ab_new_test ( );
void r8_uniform_01_test ( );
void r8_uniform_ab_test ( );
void r8col_uniform_abvec_new_test ( );
void r8mat_uniform_ab_new_test ( );
void r8row_uniform_abvec_new_test ( );
void r8vec_uniform_01_new_test ( );
void r8vec_uniform_ab_new_test ( );
void r8vec_uniform_abvec_new_test ( );

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
//    UNIFORM_PRB tests the UNIFORM library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2014
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

  bvec_uniform_new_test ( );

  c4_uniform_01_test ( );
  c4mat_uniform_01_new_test ( );
  c4vec_uniform_01_new_test ( );

  c8_uniform_01_test ( );
  c8mat_uniform_01_new_test ( );
  c8vec_uniform_01_new_test ( );

  ch_uniform_ab_test ( );

  get_seed_test ( );

  i4_seed_advance_test ( );

  i4_uniform_0i_test ( );
  i4_uniform_ab_test ( );
  i4mat_uniform_ab_new_test ( );
  i4vec_uniform_ab_new_test ( );

  l4_uniform_test ( );
  l4mat_uniform_new_test ( );
  l4vec_uniform_new_test ( );

  lcrg_anbn_test ( );
  lcrg_seed_test ( );

  r4_uniform_01_test ( );
  r4_uniform_ab_test ( );
  r4mat_uniform_ab_new_test ( );
  r4vec_uniform_ab_new_test ( );

  r8_uniform_01_test ( );
  r8_uniform_ab_test ( );
  r8mat_uniform_ab_new_test ( );
  r8vec_uniform_01_new_test ( );
  r8vec_uniform_ab_new_test ( );

  r8col_uniform_abvec_new_test ( );
  r8row_uniform_abvec_new_test ( );
  r8vec_uniform_abvec_new_test ( );
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

void bvec_uniform_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_UNIFORM_NEW_TEST tests BVEC_UNIFORM_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  int i;
  int n = 10;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "BVEC_UNIFORM_NEW_TEST\n";
  cout << "  BVEC_UNIFORM_NEW computes a binary vector.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";
  
  for ( i = 0; i < 10; i++ )
  {
    b = bvec_uniform_new ( n, seed );
    bvec_print ( n, b, "" );
    delete [] b;
  }

  return;
}
//****************************************************************************80

void c4_uniform_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4_UNIFORM_01_TEST tests C4_UNIFORM_01.
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
  complex <float> value;

  seed = 123456789;

  cout << "\n";
  cout << "C4_UNIFORM_01_TEST\n";
  cout << "  C4_UNIFORM_01 computes pseudorandom complex values\n";
  cout << "  uniformly distributed in the unit circle.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

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

void c4mat_uniform_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_UNIFORM_01_NEW_TEST tests C4MAT_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> *c;
  int m;
  int n;
  int seed;

  m = 5;
  n = 2;
  seed = 123456789;

  cout << "\n";
  cout << "C4MAT_UNIFORM_01_NEW_TEST\n";
  cout << "  C4MAT_UNIFORM_01_NEW computes pseudorandom complex values\n";
  cout << "  uniformly distributed in the unit circle.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

  c = c4mat_uniform_01_new ( m, n, seed );

  c4mat_print ( m, n, c, "  Uniform C4MAT:" );

  delete [] c;

  return;
}
//****************************************************************************80

void c4vec_uniform_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNIFORM_01_NEW_TEST tests C4VEC_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> *c;
  int n;
  int seed;

  n = 10;
  seed = 123456789;

  cout << "\n";
  cout << "C4VEC_UNIFORM_01_NEW_TEST\n";
  cout << "  C4VEC_UNIFORM_01_NEW computes pseudorandom complex values\n";
  cout << "  uniformly distributed in the unit circle.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

  c = c4vec_uniform_01_new ( n, seed );

  c4vec_print ( n, c, "  Uniform C4VEC:" );

  delete [] c;

  return;
}
//****************************************************************************80

void c8_uniform_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8_UNIFORM_01_TEST tests C8_UNIFORM_01.
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
  complex <double> value;

  seed = 123456789;

  cout << "\n";
  cout << "C8_UNIFORM_01_TEST\n";
  cout << "  C8_UNIFORM_01 computes pseudorandom C8 values\n";
  cout << "  uniformly distributed in the unit circle.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

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

void c8mat_uniform_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_UNIFORM_01_NEW_TEST tests C8MAT_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  complex <double> *c;
  int m;
  int n;
  int seed;

  m = 5;
  n = 2;
  seed = 123456789;

  cout << "\n";
  cout << "C8MAT_UNIFORM_01_NEW_TEST\n";
  cout << "  C8MAT_UNIFORM_01_NEW computes pseudorandom complex values\n";
  cout << "  uniformly distributed in the unit circle.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

  c = c8mat_uniform_01_new ( m, n, seed );

  c8mat_print ( m, n, c, "  Uniform C8MAT:" );

  delete [] c;

  return;
}
//****************************************************************************80

void c8vec_uniform_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_UNIFORM_01_NEW_TEST tests C8VEC_UNIFORM_01_NEW.
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
  int n;
  int seed;
  complex <double> *c;

  n = 10;
  seed = 123456789;

  cout << "\n";
  cout << "C8VEC_UNIFORM_01_NEW_TEST\n";
  cout << "  C8VEC_UNIFORM_01_NEW computes pseudorandom C8 values\n";
  cout << "  uniformly distributed in the unit circle.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

  c = c8vec_uniform_01_new ( n, seed );

  c8vec_print ( n, c, "  Uniform C8VEC:" );

  delete [] c;

  return;
}
//****************************************************************************80

void ch_uniform_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CH_UNIFORM_AB_TEST tests CH_UNIFORM_AB.
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

  clo = 'A';
  chi = 'J';
  seed = 123456789;

  cout << "\n";
  cout << "CH_UNIFORM_AB_TEST\n";
  cout << "  CH_UNIFORM_AB computes pseudorandom characters\n";
  cout << "  in an interval [CLO,CHI].\n";
  cout << "\n";
  cout << "  The lower endpoint CLO = '" << clo << "'.\n";
  cout << "  The upper endpoint CHI = '" << chi << "'.\n";
  wcout << "  The initial seed is " << seed << "\n";

  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  '" << ch_uniform_ab ( clo, chi, seed ) << "'\n";
  }

  return;
}
//****************************************************************************80

void get_seed_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED_TEST tests GET_SEED.
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

  oops = 0;
  seed = 12345678;
  seed_old = seed;

  cout << "\n";
  cout << "GET_SEED_TEST\n";
  cout << "  GET_SEED picks an initial seed value for R8_UNIFORM_01.\n";
  cout << "  The value chosen should vary over time, because\n";
  cout << "  the seed is based on reading the clock.\n";
  cout << "\n";
  cout << "  This is just the \"calendar\" clock, which does\n";
  cout << "  not change very fast, so calling GET_SEED several\n";
  cout << "  times in a row may result in the same value.\n";
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

void i4_seed_advance_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SEED_ADVANCE_TEST tests I4_SEED_ADVANCE.
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

  seed_new = 12345;

  cout << "\n";
  cout << "I4_SEED_ADVANCE_TEST\n";
  cout << "  I4_SEED_ADVANCE advances the seed.\n";
  cout << "\n";
  cout << "  Step        SEED input       SEED output\n";
  cout << "\n";

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

void i4_uniform_0i_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_0I_TEST tests I4_UNIFORM_0I
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

  seed = 123456789;

  cout << "\n";
  cout << "I4_UNIFORM_0I_TEST\n";
  cout << "  I4_UNIFORM_0I samples a uniform random\n";
  cout << "  integer distribution in [0,2**31-1].\n";
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

void i4_uniform_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
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
  int seed = 123456789;

  cout << "\n";
  cout << "I4_UNIFORM_AB_TEST\n";
  cout << "  I4_UNIFORM_AB computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";
  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";

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

void i4mat_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_UNIFORM_AB_NEW_TEST tests I4MAT_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int a = -100;
  int b = 200;
  int m = 5;
  int n = 4;
  int seed = 123456789;
  int *v;

  cout << "\n";
  cout << "I4MAT_UNIFORM_AB_NEW_TEST\n";
  cout << "  I4MAT_UNIFORM_AB_NEW computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";
  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";

  v = i4mat_uniform_ab_new ( m, n, a, b, seed );

  i4mat_print ( m, n, v, "  Uniform I4MAT:" );

  delete [] v;

  return;
}
//****************************************************************************80

void i4vec_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIFORM_AB_NEW_TEST tests I4VEC_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int a = -100;
  int b = 200;
  int n = 20;
  int seed = 123456789;
  int *v;

  cout << "\n";
  cout << "I4VEC_UNIFORM_AB_NEW_TEST\n";
  cout << "  I4VEC_UNIFORM_AB_NEW computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";
  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";

  v = i4vec_uniform_ab_new ( n, a, b, seed );

  i4vec_print ( n, v, "  Uniform I4VEC:" );

  delete [] v;

  return;
}
//****************************************************************************80

void l4_uniform_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    L4_UNIFORM_TEST tests L4_UNIFORM.
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

  seed = 123456789;

  cout << "\n";
  cout << "L4_UNIFORM_TEST\n";
  cout << "  L4_UNIFORM computes pseudorandom logical values.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
   cout << "  " << setw(8) << i
        << "  " << setw(1) <<  l4_uniform ( seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void l4mat_uniform_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    L4MAT_UNIFORM_NEW_TEST tests L4MAT_UNIFORM_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  bool *l;
  int m = 5;
  int n = 4;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "L4MAT_UNIFORM_NEW_TEST\n";
  cout << "  L4MAT_UNIFORM_NEW computes a vector of\n";
  cout << "  pseudorandom logical values.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

  l = l4mat_uniform_new ( m, n, seed );

  l4mat_print ( m, n, l, "  Uniform L4MAT:" );

  delete [] l;

  return;
}
//****************************************************************************80

void l4vec_uniform_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    L4VEC_UNIFORM_NEW_TEST tests L4VEC_UNIFORM_NEW.
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
  bool *l;
  int n = 10;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "L4VEC_UNIFORM_NEW_TEST\n";
  cout << "  L4VEC_UNIFORM_NEW computes a vector of\n";
  cout << "  pseudorandom logical values.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

  l = l4vec_uniform_new ( n, seed );

  l4vec_print ( n, l, "  Uniform L4VEC:" );

  delete [] l;

  return;
}
//****************************************************************************80

void lcrg_anbn_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LCRG_ANBN_TEST tests LCRG_ANBN.
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
//
//  These parameters define the old (1969) IBM 360 random number generator:
//
  a = 16807;
  b = 0;
  c = 2147483647;

  cout << "\n";
  cout << "LCRG_ANBN_TEST\n";
  cout << "  LCRG_ANBN determines a linear congruential random\n";
  cout << "  number generator equivalent to N steps of a given one.\n";
  cout << "\n";
  cout << " LCRG parameters:\n";
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

void lcrg_seed_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LCRG_SEED_TEST tests LCRG_SEED
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

  cout << "\n";
  cout << "LCRG_SEED_TEST\n";
  cout << "  LCRG_SEED directly computes the updated value of a\n";
  cout << "  seed used by an linear congruential random number\n";
  cout << "  generator.\n";
  cout << "\n";
  cout << "       I          SEED          SEED          SEED    U\n";
  cout << "                 Input        Output          LCRG\n";
  cout << "\n";
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

void r4_uniform_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_01_TEST tests R4_UNIFORM_01.
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

  seed = 123456789;

  cout << "\n";
  cout << "R4_UNIFORM_01_TEST\n";
  cout << "  R4_UNIFORM_01 computes pseudorandom values \n";
  cout << "  in the interval [0,1].\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << r4_uniform_01 ( seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void r4_uniform_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_AB_TEST tests R4_UNIFORM_AB.
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

  a = 5.0;
  b = 10.0;
  seed = 123456789;

  cout << "\n";
  cout << "R4_UNIFORM_AB_TEST\n";
  cout << "  R4_UNIFORM_AB computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";
  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";

  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r4_uniform_ab ( a, b, seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void r4mat_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_UNIFORM_AB_NEW_TEST tests R4MAT_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  float a = -5.0;
  float b = 10.0;
  int m = 5;
  int n = 4;
  int seed = 123456789;
  float *v;

  cout << "\n";
  cout << "R4MAT_UNIFORM_AB_NEW_TEST\n";
  cout << "  R4MAT_UNIFORM_AB_NEW computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";
  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  v = r4mat_uniform_ab_new ( m, n, a, b, seed );

  r4mat_print ( m, n, v, "  Uniform R4MAT:" );

  delete [] v;

  return;
}
//****************************************************************************80

void r4vec_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_UNIFORM_AB_NEW_TEST tests R4VEC_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  float a = -5.0;
  float b = 10.0;
  int n = 20;
  int seed = 123456789;
  float *v;

  cout << "\n";
  cout << "R4VEC_UNIFORM_AB_NEW_TEST\n";
  cout << "  R4VEC_UNIFORM_AB_NEW computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";
  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  v = r4vec_uniform_ab_new ( n, a, b, seed );

  r4vec_print ( n, v, "  Uniform R4VEC:" );

  delete [] v;

  return;
}
//****************************************************************************80

void r8_uniform_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01_TEST tests R8_UNIFORM_01.
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

  seed = 123456789;

  cout << "\n";
  cout << "R8_UNIFORM_01_TEST\n";
  cout << "  R8_UNIFORM_01 computes pseudorandom values \n";
  cout << "  in the interval [0,1].\n";
  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";

  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << r8_uniform_01 ( seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void r8_uniform_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB_TEST tests R8_UNIFORM_AB.
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

  a = 5.0;
  b = 10.0;
  seed = 123456789;

  cout << "\n";
  cout << "R8_UNIFORM_AB_TEST\n";
  cout << "  R8_UNIFORM_AB computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";
  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";

  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r8_uniform_ab ( a, b, seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void r8mat_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_AB_NEW_TEST tests R8MAT_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a = -5.0;
  double b = 10.0;
  int m = 5;
  int n = 4;
  int seed = 123456789;
  double *v;

  cout << "\n";
  cout << "R8MAT_UNIFORM_AB_NEW_TEST\n";
  cout << "  R8MAT_UNIFORM_AB_NEW computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";
  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  v = r8mat_uniform_ab_new ( m, n, a, b, seed );

  r8mat_print ( m, n, v, "  Uniform R8MAT:" );

  delete [] v;

  return;
}
//****************************************************************************80

void r8vec_uniform_01_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW_TEST tests R8VEC_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n = 10;
  int seed;
  double *v;

  seed = 123456789;

  cout << "\n";
  cout << "R8VEC_UNIFORM_01_NEW_TEST\n";
  cout << "  R8VEC_UNIFORM_01_NEW computes a random R8VEC.\n";
  cout << "\n";
  cout << "  Initial seed is " << seed << "\n";

  v = r8vec_uniform_01_new ( n, seed );

  r8vec_print ( n, v, "  Uniform R8VEC:" );

  delete [] v;

  return;
# undef N
}
//****************************************************************************80

void r8vec_uniform_ab_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_AB_NEW_TEST tests R8VEC_UNIFORM_AB_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int n = 10;
  int seed;
  double *v;

  a = -1.0;
  b = 5.0;
  seed = 123456789;

  cout << "\n";
  cout << "R8VEC_UNIFORM_AB_NEW_TEST\n";
  cout << "  R8VEC_UNIFORM_AB_NEW computes a random R8VEC.\n";
  cout << "\n";
  cout << "  " << a << " <= X <= " << b << "\n";
  cout << "  Initial seed is " << seed << "\n";

  v = r8vec_uniform_ab_new ( n, a, b, seed );

  r8vec_print ( n, v, "  Uniform R8VEC:" );

  delete [] v;

  return;
# undef N
}
//****************************************************************************80

void r8col_uniform_abvec_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_UNIFORM_ABVEC_NEW_TEST tests R8COL_UNIFORM_ABVEC_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[5] = { 0.0, 0.20, 10.0, 52.0, -1.0 };
  double b[5] = { 1.0, 0.25, 20.0, 54.0, +1.0 };
  int i;
  int j;
  int m = 5;
  int n = 4;
  int seed;
  double *v;

  seed = 123456789;

  cout << "\n";
  cout << "R8COL_UNIFORM_ABVEC_NEW_TEST\n";
  cout << "  R8COL_UNIFORM_ABVEC_NEW computes a random R8COL.\n";
  cout << "\n";
  cout << "  Initial seed is " << seed << "\n";

  v = r8col_uniform_abvec_new ( m, n, a, b, seed );

  cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(8) << a[i] << ":  ";
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << v[i+j*m];
    }
    cout << "    :" << setw(8) <<  b[i] << "\n";
  }

  delete [] v;

  return;
}
//****************************************************************************80

void r8row_uniform_abvec_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8ROW_UNIFORM_ABVEC_NEW_TEST tests R8ROW_UNIFORM_ABVEC_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[5] = { 0.0, 0.20, 10.0, 52.0, -1.0 };
  double b[5] = { 1.0, 0.25, 20.0, 54.0, +1.0 };
  int i;
  int j;
  int m = 4;
  int n = 5;
  int seed;
  double *v;

  seed = 123456789;

  cout << "\n";
  cout << "R8ROW_UNIFORM_ABVEC_NEW_TEST\n";
  cout << "  R8ROW_UNIFORM_ABVEC_NEW computes a random R8ROW.\n";
  cout << "\n";
  cout << "  Initial seed is " << seed << "\n";
  cout << "\n";

  v = r8row_uniform_abvec_new ( m, n, a, b, seed );

  for ( j = 0; j < n; j++ )
  {
    cout << "  " << setw(8) << b[j];
  }
  cout << "\n";
  cout << "\n";

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(8) << v[i+j*m];
    }
    cout << "\n";
  }

  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout << "  " << setw(8) << a[j];
  }
  cout << "\n";;

  delete [] v;

  return;
}
//****************************************************************************80

void r8vec_uniform_abvec_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_ABVEC_NEW_TEST tests R8VEC_UNIFORM_ABVEC_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[5] = { 0.0, 0.20, 10.0, 52.0, -1.0 };
  double b[5] = { 1.0, 0.25, 20.0, 54.0, +1.0 };
  int i;
  int n = 5;
  int seed;
  double *v;

  seed = 123456789;

  cout << "\n";
  cout << "R8VEC_UNIFORM_ABVEC_NEW_TEST\n";
  cout << "  R8VEC_UNIFORM_ABVEC_NEW computes a random R8VEC.\n";
  cout << "\n";
  cout << "  Initial seed is " << seed << "\n";

  v = r8vec_uniform_abvec_new ( n, a, b, seed );

  cout << "\n";
  cout << "   I         A         X         B\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(8) << a[i]
         << "  " << setw(8) << v[i]
         << "  " << setw(8) << b[i] << "\n";
  }

  delete [] v;

  return;
# undef N
}
