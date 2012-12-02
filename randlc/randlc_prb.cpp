# include <cstdlib>
# include <iostream>
# include <iomanip>

# include "randlc.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    RANDLC_PRB calls sample problems for the RANDLC library.
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
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "RANDLC_PRB\n";
  cout << "  C++ version:\n";
  cout << "  Test the RANDLC library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "RANDLC_PRB\n";
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
//    TEST01 tests RANDLC.
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
//    John Burkardt
//
{
  int i;
  double seed;
  double seed_init = 123456789.0;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  RANDLC computes pseudorandom values \n";
  cout << "  in the interval [0,1].\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";
  cout << "         I          RANDLC\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(14) << randlc ( &seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests RANDLC;
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
//    John Burkardt
//
{
# define N 1000

  int i;
  double seed;
  double seed_in;
  double seed_out;
  double u[N];
  double u_avg;
  double u_var;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  RANDLC computes a sequence of uniformly distributed\n";
  cout << "  pseudorandom numbers.\n";

  seed = 123456789.0;

  cout << "\n";
  cout << "  Initial SEED = " << ( long long int ) seed << "\n";

  cout << "\n";
  cout << "  First 10 values:\n";
  cout << "\n";
  cout << "       I           Input          Output      RANDLC\n";
  cout << "                    SEED            SEED\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    seed_in = seed;
    u[i] = randlc ( &seed );
    seed_out = seed;
    cout << "  " << setw(6) << i + 1
         << "  " << setw(14) << ( long long int ) seed_in
         << "  " << setw(14) << ( long long int ) seed_out
         << "  " << setw(10) << u[i] << "\n";
  }

  cout << "\n";
  cout << "  Now call RANDLC " << N << " times.\n";

  u_avg = 0.0;
  for ( i = 0; i < N; i++ )
  {
    u[i] = randlc ( &seed );
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
  cout << "  Average value = " << u_avg << "\n";
  cout << "  Expecting       " << 0.5 << "\n";

  cout << "\n";
  cout << "  Variance =      " << u_var << "\n";
  cout << "  Expecting       " << 1.0 / 12.0 << "\n";

  return;
# undef N
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests RANDLC.
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
//    John Burkardt
//
{
  int i;
  double seed;
  double seed_in;
  double seed_out;
  double seed_save;
  double x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  RANDLC computes a sequence of pseudorandom numbers\n";
  cout << "  but all computations depend on the seed value.\n";
  cout << "  In this test, we show how a sequence of \"random\"\n";
  cout << "  values can be manipulated by accessing the seed.\n";

  seed = 1066.0;

  cout << "\n";
  cout << "  Set SEED to " << ( long long int ) seed << "\n";
  cout << "\n";
  cout << "  Now call RANDLC 10 times, and watch SEED.\n";
  cout << "\n";
  cout << "       I           Input          Output      RANDLC\n";
  cout << "                    SEED            SEED\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;

    if ( i == 5 )
    {
      seed_save = seed;
    }
    x = randlc ( &seed );
    seed_out = seed;
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << ( long long int ) seed_in
         << "  " << setw(14) << ( long long int ) seed_out
         << "  " << setw(10) << x << "\n";
  }

  seed = seed_save;

  cout << "\n";
  cout << "  Reset SEED to its value at step 5, = " << ( long long int ) seed << "\n";
  cout << "\n";
  cout << "  Now call RANDLC 10 times, and watch how SEED\n";
  cout << "  and RANDLC restart themselves.\n";
  cout << "\n";
  cout << "       I           Input          Output      RANDLC\n";
  cout << "                    SEED            SEED\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = randlc ( &seed );
    seed_out = seed;
    cout << "  " << setw(6) << i
         << "  " << setw(14) << ( long long int ) seed_in
         << "  " << setw(14) << ( long long int ) seed_out
         << "  " << setw(10) << x << "\n";
  }

  seed = 0.0;

  cout << "\n";
  cout << "  What happens with an initial zero SEED?\n";
  cout << "\n";
  cout << "       I           Input          Output      RANDLC\n";
  cout << "                    SEED            SEED\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = randlc ( &seed );
    seed_out = seed;
    cout << "  " << setw(6) << i
         << "  " << setw(14) << ( long long int ) seed_in
         << "  " << setw(14) << ( long long int ) seed_out
         << "  " << setw(10) << x << "\n";
  }

  seed = -123456789.0;

  cout << "\n";
  cout << "  What happens with an initial negative SEED?\n";
  cout << "\n";
  cout << "       I           Input          Output      RANDLC\n";
  cout << "                    SEED            SEED\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = randlc ( &seed );
    seed_out = seed;
    cout << "  " << setw(6) << i
         << "  " << setw(14) << ( long long int ) seed_in
         << "  " << setw(14) << ( long long int ) seed_out
         << "  " << setw(10) << x << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    RANDLC_TEST04 tests RANDLC_JUMP.
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
{
  int i;
  int k;
  int klog;
  double seed;
  double x1;
  double x2;

  cout << "\n";
  cout << "RANDLC_TEST04\n";
  cout << "  RANDLC_JUMP jumps directly to the K-th value\n";
  cout << "  returned by RANDLC.\n";
  cout << "\n";
  cout << "         K X(hard way)     X(jump)\n";
  cout << "\n";

  k = 1;

  for ( klog = 1; klog <= 10; klog++ )
  {
    seed = 123456789.0;
    for ( i = 1; i <= k; i++ )
    {
      x1 = randlc ( &seed );
    }

    seed = 123456789.0;
    x2 = randlc_jump ( seed, k );

    cout << "  " << setw(8) << k
         << "  " << setw(10) << x1
         << "  " << setw(10) << x2 << "\n";

    k = k * 2;
  }

  return;
}
