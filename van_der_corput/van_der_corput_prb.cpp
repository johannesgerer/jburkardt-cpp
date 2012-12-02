# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "van_der_corput.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test045 ( );
void test05 ( );
void test06 ( );
void test09 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for VAN_DER_CORPUT_PRB.
//
//  Discussion:
//
//    VAN_DER_CORPUT_PRB calls a set of problems for VAN_DER_CORPUT.
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
{
  timestamp ( );

  cout << "\n";
  cout << "VAN_DER_CORPUT_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the VAN_DER_CORPUT library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test045 ( );
  test05 ( );
  test06 ( );
  test09 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "VAN_DER_CORPUT_PRB\n";
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
//    TEST01 tests VAN_DER_CORPUT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
{
  int i;
  double r;
  int seed;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  VAN_DER_CORPUT computes the elements of a\n";
  cout << "  van der Corput sequence.\n";
  cout << "  Each call produces the next value.  By default,\n";
  cout << "  the base is 2, and the sequence starts at element 1.\n";
  cout << "\n";
  cout << "  In this test, we call VAN_DER_CORPUT several times.\n";
  cout << "\n";
  cout << "  Seed   van der Corput\n";
  cout << "\n";

  for ( i = 1; i <=10; i++ )
  {
    seed = van_der_corput_seed_get ( );

    r = van_der_corput ( );

    cout << setw(6)  << seed << "  "
         << setw(10) << r    << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests VAN_DER_CORPUT_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
{
# define N 10

  int base;
  int i;
  int n;
  double r[N];
  int seed;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  VAN_DER_CORPUT_SEQUENCE computes several elements of\n";
  cout << "  a van der Corput sequence on a single call.\n";
  cout << "\n";
  cout << "  In this test, we call VAN_DER_CORPUT_SEQUENCE once.\n";

  base = 2;
  van_der_corput_base_set ( base );

  seed = 0;
  van_der_corput_seed_set ( seed );

  van_der_corput_sequence ( N, r );

  cout << "\n";
  cout << "  Element   van der Corput\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(6)  << i    << "  "
         << setw(10) << r[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests VAN_DER_CORPUT_SEQUENCE, VAN_DER_CORPUT_SEED_GET, VAN_DER_CORPUT_SEED_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
{
# define NMAX 10

  int i;
  int n;
  double r[NMAX];
  int seed;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  VAN_DER_CORPUT_SEED_SET specifies the next element of\n";
  cout << "    the van der Corput sequence to compute.\n";
  cout << "  VAN_DER_CORPUT_SEED_GET reports the next element of the\n";
  cout << "    van der Corput sequence that will be computed.\n";
  cout << "\n";
  cout << "  By default, the sequence starts at element 1.\n";
  cout << "\n";
  cout << "  In this test, we demonstrate computing elements\n";
  cout << "  affects the seed, and how resetting the seed determines\n";
  cout << "  the next element computed.\n";
  cout << "\n";
  cout << "  We start at element 0 and compute 10 elements.\n";
  cout << "\n";

  seed = 0;
  van_der_corput_seed_set ( seed );

  n = NMAX;
  van_der_corput_sequence ( n, r );

  for ( i = 0; i < n; i++ )
  {
    cout << setw(6)  << seed+i << " "
         << setw(10) << r[i]   << "\n";
  }

  seed = van_der_corput_seed_get ( );

  cout << "\n";
  cout << "  The current seed is " << seed << "\n";

  cout << "\n";
  cout << "  We jump back to element 6 and compute 10 elements.\n";
  cout << "\n";

  seed = 6;
  van_der_corput_seed_set ( seed );

  n = NMAX;
  van_der_corput_sequence ( n, r );

  for ( i = 0; i < n; i++ )
  {
    cout << setw(6)  << seed+i << " "
         << setw(10) << r[i]   << "\n";
  }

  seed = van_der_corput_seed_get ( );

  cout << "\n";
  cout << "  The current seed is " << seed << "\n";

  cout << "\n";
  cout << "  We restart at element 0 and compute 6 elements.\n";
  cout << "\n";

  seed = 0;
  van_der_corput_seed_set ( seed );

  n = 6;
  van_der_corput_sequence ( n, r );

  for ( i = 0; i < n; i++ )
  {
    cout << setw(6)  << seed+i << " "
         << setw(10) << r[i]   << "\n";
  }

  seed = van_der_corput_seed_get ( );

  cout << "\n";
  cout << "  The current seed is " << seed << "\n";

  cout << "\n";
  cout << "  We jump to element 100 and compute 5 elements.\n";
  cout << "\n";

  seed = 100;
  van_der_corput_seed_set ( seed );

  n = 5;
  van_der_corput_sequence ( n, r );

  for ( i = 0; i < n; i++ )
  {
    cout << setw(6)  << seed+i << " "
         << setw(10) << r[i]   << "\n";
  }

  seed = van_der_corput_seed_get ( );

  cout << "\n";
  cout << "  The current seed is " << seed << "\n";

  return;
#undef N
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests VAN_DER_CORPUT, VAN_DER_CORPUT_BASE_GET, VAN_DER_CORPUT_BASE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
{
  int base;
  int i;
  double r;
  int seed;

  cout << "\n";
  cout << "TEST04";
  cout << "  VAN_DER_CORPUT_BASE_GET gets the current base.\n";
  cout << "  VAN_DER_CORPUT_BASE_SET sets the current base.\n";
  cout << "\n";
  cout << "  The van der Corput base is usually a prime, but this is\n";
  cout << "  not required.\n";
  cout << "\n";
  cout << "  In this test, we compute a van der Corput sequence\n";
  cout << "  with the default base, then change the base,\n";
  cout << "  reset the seed, and recompute the sequence.\n";

  seed = 0;
  van_der_corput_seed_set ( seed );

  base = van_der_corput_base_get ( );

  cout << "\n";
  cout << "  VAN_DER_CORPUT_BASE_GET: Current base is " << base << "\n";
  cout << "\n";
  cout << "  Seed   van der Corput\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed = van_der_corput_seed_get ( );
    r = van_der_corput ( );

    cout << setw(6)  << seed << "  "
         << setw(10) << r    << "\n";
  }

  base = 3;
  van_der_corput_base_set ( base );

  seed = 0;
  van_der_corput_seed_set ( seed );

  cout << "\n";
  cout << "  Reset base to " << base << "\n";
  cout << "  Reset seed to " << seed << "\n";

  cout << "\n";
  cout << "  Seed   van der Corput\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed = van_der_corput_seed_get ( );
    r = van_der_corput ( );

    cout << setw(6)  << seed << "  "
         << setw(10) << r    << "\n";
  }

  base = 4;
  van_der_corput_base_set ( base );

  seed = 0;
  van_der_corput_seed_set ( seed );

  cout << "\n";
  cout << "  Set BASE = " << base << "\n";
  cout << "  Set SEED = " << seed << "\n";

  cout << "\n";
  cout << "  Seed   van der Corput\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed = van_der_corput_seed_get ( );
    r = van_der_corput ( );

    cout << setw(6)  << seed << "  "
         << setw(10) << r    << "\n";
  }

  return;
}
//****************************************************************************80

void test045 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST045 tests VAN_DER_CORPUT, VAN_DER_CORPUT_SEED_GET, VAN_DER_CORPUT_SEED_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
{
  int base;
  int i;
  double r;
  int seed;

  cout << "\n";
  cout << "TEST045";
  cout << "  VAN_DER_CORPUT_SEED_GET gets the current seed.\n";
  cout << "  VAN_DER_CORPUT_SEED_SET sets the current seed.\n";
  cout << "\n";
  cout << "  The van der Corput base is usually a prime, but this is\n";
  cout << "  not required.\n";
  cout << "\n";
  cout << "  In this test, we compute a van der Corput sequence\n";
  cout << "  starting with the default seed, then check the seed,\n";
  cout << "  reset the seed, and recompute the sequence.\n";

  base = 2;
  van_der_corput_base_set ( base );

  seed = 0;
  van_der_corput_seed_set ( seed );

  cout << "\n";
  cout << "  All computations will use base " << base << ".\n";
  cout << "\n";

  seed = van_der_corput_seed_get ( );

  cout << "\n";
  cout << "  Set SEED = " << seed << "\n";
  cout << "\n";
  cout << "  Seed   van der Corput\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed = van_der_corput_seed_get ( );

    r = van_der_corput ( );

    cout << setw(6)  << seed << "  "
         << setw(10) << r    << "\n";
  }

  seed = van_der_corput_seed_get ( );

  cout << "\n";
  cout << "  Current seed is " << seed << "\n";
  cout << "\n";

  seed = 100;
  van_der_corput_seed_set ( seed );

  cout << "\n";
  cout << "  Set SEED = " << seed << "\n";
  cout << "\n";
  cout << "  Seed   van der Corput\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed = van_der_corput_seed_get ( );

    r = van_der_corput ( );

    cout << setw(6)  << seed << "  "
         << setw(10) << r    << "\n";
  }

  seed = van_der_corput_seed_get ( );

  cout << "\n";
  cout << "  Current seed is " << seed << "\n";
  cout << "\n";

  seed = 3;
  van_der_corput_seed_set ( seed );

  cout << "\n";
  cout << "  Reset seed to " << seed << "\n";
  cout << "\n";
  cout << "  Seed   van der Corput\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed = van_der_corput_seed_get ( );

    r = van_der_corput ( );

    cout << setw(6)  << seed << "  "
         << setw(10) << r    << "\n";
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests I4_TO_VAN_DER_CORPUT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
{
  int base;
  double r;
  int seed;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  I4_TO_VAN_DER_CORPUT returns the I-th element\n";
  cout << "  of a van der Corput sequence to a given base.\n";
  cout << "\n";
  cout << "\n";
  cout << "  Base    Seed   R\n";
  cout << "\n";

  for ( base = 2; base <= 5; base++ )
  {

    cout << "\n";
    cout << setw ( 6 ) << base << "\n";

    for ( seed = 0; seed <= 10; seed++ )
    {
      r = i4_to_van_der_corput ( seed, base );

      cout << "        "
           << setw ( 6 )  << seed << "  "
           << setw ( 10 ) << r    << "\n";
    }

  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests I4_TO_VAN_DER_CORPUT_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
{
# define N 10

  int base;
  int i;
  double r[N];
  int seed;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  I4_TO_VAN_DER_CORPUT_SEQUENCE returns N elements\n";
  cout << "  of a van der Corput sequence to a given base.\n";
  cout << "\n";
  cout << "\n";
  cout << "  Base    Seed   R\n";
  cout << "\n";

  for ( base = 2; base <= 5; base++ )
  {

    cout << "\n";
    cout << setw ( 6 ) << base << "\n";

    seed = 0;

    i4_to_van_der_corput_sequence ( seed, base, N, r );

    for ( i = 0; i < N; i++ )
    {

      cout << "        "
           << setw ( 6 ) << seed+i << "  "
           << setw ( 10 ) << r[i] << "\n";
    }

  }

  return;
# undef N
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests VDC_NUMERATOR_SEQUENCE.
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
{
  int i;
  int n;
  int *r;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  VDC_NUMERATOR_SEQUENCE returns N elements\n";
  cout << "  of a van der Corput numerator sequence in base 2.\n";
  cout << "\n";
  cout << "   N:  Sequence\n";
  cout << "\n";

  for ( n = 1; n <= 20; n++ )
  {
    r = vdc_numerator_sequence ( n );
    cout << "  " << setw(2) << n << ":";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(2) << r[i];
    }
    cout << "\n";

    delete [] r;
  }
  return;
}
