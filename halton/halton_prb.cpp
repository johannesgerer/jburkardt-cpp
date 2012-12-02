# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>

using namespace std;

#include "halton.hpp"

int main ( );
void test01 ( );
void test0125 ( );
void test0126 ( );
void test02 ( );
void test025 ( );
void test03 ( );
void test04 ( );
void test045 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    HALTON_PRB calls the HALTON tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "HALTON_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the HALTON library.\n";

  test01 ( );
  test0125 ( );
  test0126 ( );
  test02 ( );
  test025 ( );
  test03 ( );
  test04 ( );
  test045 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HALTON_PRB:\n";
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
//    TEST01 tests HALTON, HALTON_STEP_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_MAX 3
# define TEST_NUM 4

  int base[DIM_MAX];
  int i;
  int j;
  int n;
  int dim_num;
  double r[DIM_MAX];
  int seed[DIM_MAX];
  int step_vec[TEST_NUM] = { 0, 5, 1000, 1000000 };
  int step;
  int test;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  HALTON computes the next element of a Halton sequence.\n";
  cout << "  HALTON_STEP_SET sets the step.\n";
  cout << "\n";
  cout << "  In this test, we try several values of STEP.\n";
  cout << "  We repeat the test for several dimensions.\n";
  cout << "  We assume defaults for SEED, LEAP and BASE.\n";

  for ( dim_num = 1; dim_num <= DIM_MAX; dim_num++ )
  {
    for ( test = 0; test < TEST_NUM; test++ )
    {
      halton_dim_num_set ( dim_num );
      n = 11;
      step = step_vec[test];
      halton_step_set ( step );
      for ( i = 0; i < dim_num; i++ )
      {
        seed[i] = 0;
      }
      halton_seed_set ( seed );
      for ( i = 0; i < dim_num; i++ )
      {
        base[i] = prime ( i + 1 );
      }
      halton_base_set ( base );

      cout << "\n";
      cout << "  DIM_NUM = " << setw(12) << dim_num << "\n";
      cout << "  N =    " << setw(12) << n    << "\n";
      cout << "  STEP = " << setw(12) << step << "\n";
      i4vec_transpose_print ( dim_num, seed, "  SEED = " );
      i4vec_transpose_print ( dim_num, base, "  BASE = " );

      cout << "\n";
      cout << "    STEP   Halton\n";
      cout << "\n";
      for ( j = 0; j < n; j++ )
      {
        halton ( r );
        cout << setw(6) << step + j << "  ";
        for ( i = 0; i < dim_num; i++ )
        {
          cout << setw(12) << r[i] << "  ";
        }
        cout << "\n";
      }
    }
  }

  return;
# undef DIM_MAX
# undef TEST_NUM
}
//****************************************************************************80

void test0125 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0125 tests I4_TO_HALTON.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_MAX 3

  int base[DIM_MAX];
  int i;
  int j;
  int leap[DIM_MAX];
  int n;
  int dim_num;
  double r[DIM_MAX];
  int seed[DIM_MAX];
  int step;

  cout << "\n";
  cout << "TEST0125\n";
  cout << "  I4_TO_HALTON computes a Halton sequence.\n";
  cout << "  The user gives the seed and bases as explicit input.\n";
  cout << "\n";
  cout << "  In this test, we call I4_TO_HALTON repeatedly.\n";
  cout << "  We use distinct primes as bases.\n";

  for ( dim_num = 1; dim_num <= DIM_MAX; dim_num++ )
  {
    n = 11;
    step = 0;
    for ( i = 0; i < dim_num; i++ )
    {
      seed[i] = 0;
    }
    for ( i = 0; i < dim_num; i++ )
    {
      leap[i] = 1;
    }
    for ( i = 0; i < dim_num; i++ )
    {
      base[i] = prime ( i + 1 );
    }

    cout << "\n";
    cout << "  DIM_NUM = " << setw(12) << dim_num << "\n";
    cout << "  N =    " << setw(12) << n    << "\n";
    cout << "  STEP = " << setw(12) << step << "\n";
    i4vec_transpose_print ( dim_num, seed, "  SEED = " );
    i4vec_transpose_print ( dim_num, leap, "  LEAP = " );
    i4vec_transpose_print ( dim_num, base, "  BASE = " );

    cout << "\n";
    cout << "    STEP      Halton\n";
    cout << "\n";
    for ( j = 0; j < n; j++ )
    {
      step = j;
      i4_to_halton ( dim_num, step, seed, leap, base, r );
      cout                    << "  "
           << setw(6) << step << "  ";
      for ( i = 0; i < dim_num; i++ )
      {
        cout << setw(8) << r[i] << "  ";
      }
      cout << "\n";
    }

  }

  return;
# undef DIM_MAX
}
//****************************************************************************80

void test0126 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0126 tests I4_TO_HALTON.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int base[DIM_NUM];
  int i;
  int j;
  int leap[DIM_NUM];
  int n;
  int dim_num;
  double r[DIM_NUM];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST0126\n";
  cout << "  I4_TO_HALTON computes a Halton sequence.\n";
  cout << "  The user gives the seed and bases as explicit input.\n";
  cout << "\n";
  cout << "  In this test, we call I4_TO_HALTON repeatedly.\n";
  cout << "  We use the same value for all bases.\n";

  n = 11;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  for ( i = 0; i <  DIM_NUM; i++ )
  {
    leap[i] = 1;
  }
  for ( i = 0; i <  DIM_NUM; i++ )
  {
    base[i] = 2;
  }

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) <<  DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << n     << "\n";
  i4vec_transpose_print (  DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print (  DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print (  DIM_NUM, base, "  BASE = " );

  cout << "\n";
  cout << "    STEP      Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    step = j;
    i4_to_halton (  DIM_NUM, step, seed, leap, base, r );
    cout                 << "  "
         << setw(6) << j << "  ";
    for ( i = 0; i <  DIM_NUM; i++ )
    {
      cout << setw(8) << r[i] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests HALTON_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 3

  int i;
  int j;
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  HALTON_SEQUENCE computes N elements of \n";
  cout << "  a Halton sequence on a single call.\n";

  halton_dim_num_set ( DIM_NUM );
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << N    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );

  halton_sequence ( N, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 tests I4_TO_HALTON_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 3

  int base[DIM_NUM];
  int i;
  int j;
  int leap[DIM_NUM];
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;
  int test;

  cout << "\n";
  cout << "TEST025\n";
  cout << "  I4_TO_HALTON_SEQUENCE computes N elements of \n";
  cout << "  a Halton sequence on a single call.\n";
  cout << "  All arguments are specified explicitly.\n";

  for ( test = 1; test <= 4; test++ )
  {
    if ( test == 1 )
    {
      step = 0;
      for ( i = 0; i < DIM_NUM; i++ )
      {
        seed[i] = 0;
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        leap[i] = 1;
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        base[i] = prime ( i + 1 );
      }
    }
    else if ( test == 2 )
    {
      step = 0;
      for ( i = 0; i < DIM_NUM; i++ )
      {
        seed[i] = i + 1;
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        leap[i] = 1;
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        base[i] = prime ( i + 1 );
      }
    }
    else if ( test == 3 )
    {
      step = 0;
      for ( i = 0; i < DIM_NUM; i++ )
      {
        seed[i] = 1;
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        leap[i] = 3;
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        base[i] = prime ( i + 1 );
      }
    }
    else if ( test == 4 )
    {
      step = 0;
      for ( i = 0; i < DIM_NUM; i++ )
      {
        seed[i] = i + 1;
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        leap[i] = 1;
      }
      for ( i = 0; i < DIM_NUM; i++ )
      {
        base[i] = 2;
      }
    }

    cout << "\n";
    cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
    cout << "  N =    " << setw(12) << N    << "\n";
    cout << "  STEP = " << setw(12) << step << "\n";
    i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
    i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
    i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

    i4_to_halton_sequence ( DIM_NUM, N, step, seed, leap, base, r );

    cout << "\n";
    cout << "    STEP   Halton\n";
    cout << "\n";
    for ( j = 0; j < N; j++ )
    {
      cout                        << "  "
           << setw(6) << step + j << "  ";
      for ( i = 0; i < DIM_NUM; i++ )
      {
        cout << setw(12) << r[i+j*DIM_NUM] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests HALTON_SEQUENCE, HALTON_STEP_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 1
# define NMAX 11

  int base[DIM_NUM];
  int i;
  int j;
  int n;
  double r[DIM_NUM*NMAX];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  HALTON_STEP_SET specifies the next element of\n";
  cout << "    the Halton sequence to compute.\n";
  cout << "\n";
  cout << "  In this test, we demonstrate how resetting \n";
  cout << "  STEP determines the next element computed.\n";

  halton_dim_num_set ( DIM_NUM );
  n = 11;
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );
  base[0] = 2;
  halton_base_set ( base );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );

  halton_sequence ( n, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout                         << "  "
         << setw(6)  << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  n = 11;
  step = 6;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  halton_sequence ( n, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout                         << "  "
         << setw(6)  << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  n = 6;
  step = 0;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  halton_sequence ( n, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout                         << "  "
         << setw(6)  << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  n = 5;
  step = 100;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  halton_sequence ( n, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout                         << "  "
         << setw(6)  << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef NMAX
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests HALTON, HALTON_BASE_GET, HALTON_BASE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 1

  int *base;
  int i;
  int j;
  int n;
  double r[DIM_NUM];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  HALTON_BASE_GET gets the current Halton bases.\n";
  cout << "  HALTON_BASE_SET sets the current Halton bases.\n";
  cout << "\n";
  cout << "  In this test, we compute a Halton sequence\n";
  cout << "  with the default base, then change the base,\n";
  cout << "  reset the seed, and recompute the sequence.\n";

  halton_dim_num_set ( DIM_NUM );
  n = 10;
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );
  base = halton_base_get ( );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N    = " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    halton ( r );
    cout                         << "  "
         << setw(6)  << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << r[i] << "  ";
    }
    cout << "\n";
  }

  n = 10;
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = 3;
  }
  halton_base_set ( base );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N    = " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    halton ( r );
    cout                         << "  "
         << setw(6)  << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << r[i] << "  ";
    }
    cout << "\n";
  }

  n = 10;
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = 4;
  }
  halton_base_set ( base );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N    = " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    halton ( r );
    cout                         << "  "
         << setw(6)  << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << r[i] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test045 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST045 tests HALTON_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 101
# define DIM_NUM 2

  int base[DIM_NUM] = { 2, 3 };
  int i;
  int j;
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST045\n";
  cout << "  HALTON_SEQUENCE computes N elements of\n";
  cout << "  a Halton sequence on a single call.\n";

  halton_dim_num_set ( DIM_NUM );
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i + 1 );
  }
  halton_base_set ( base );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N    = " << setw(12) << N    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  halton_sequence ( N, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(7) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests HALTON.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 4

  int base[DIM_NUM];
  int i;
  int j;
  int n;
  double r[DIM_NUM];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  HALTON computes the elements of a vector Halton sequence.\n";
  cout << "\n";
  cout << "  Each call produces the next value.\n";
  cout << "\n";
  cout << "  In this test, we call HALTON several times.\n";

  halton_dim_num_set ( DIM_NUM );
  n = 11;
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i + 1 );
  }
  halton_base_set ( base );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    halton ( r );
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << r[i] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests HALTON_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 11
# define DIM_NUM 4

  int i;
  int j;
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  HALTON_SEQUENCE computes the next N elements\n";
  cout << "  of a vector Halton sequence.\n";
  cout << "\n";
  cout << "  Each call produces the next value.  By default,\n";
  cout << "  the bases are the first DIM_NUM primes.\n";
  cout << "\n";
  cout << "  In this test, we demonstrate how one call can compute\n";
  cout << "  many successive vector elements of the sequence.\n";

  halton_dim_num_set ( DIM_NUM );
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << N    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );

  halton_sequence ( N, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests HALTON_STEP_SET, HALTON_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 4
# define NMAX 10

  int *base;
  int i;
  int j;
  int *leap;
  int n;
  double r[DIM_NUM*NMAX];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  HALTON_STEP_SET specifies which element of the\n";
  cout << "    Halton subsequence to compute.\n";
  cout << "\n";
  cout << "  In this test, we show how STEP chooses the next element.\n";

  halton_dim_num_set ( DIM_NUM );
  n = 10;
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );
  leap = halton_leap_get ( );
  base = halton_base_get ( );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  halton_sequence ( n, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  n = 10;
  step = 6;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  halton_sequence ( n, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  n = 6;
  step = 0;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  halton_sequence ( n, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  n = 5;
  step = 100;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  halton_sequence ( n, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < n; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef DIM_NUM
# undef NMAX
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests HALTON_BASE_GET, HALTON_BASE_SET, HALTON_SEQUENCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 4

  int base[DIM_NUM];
  int i;
  int j;
  int n;
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  HALTON_BASE_GET gets the current bases.\n";
  cout << "  HALTON_BASE_SET sets the current bases.\n";
  cout << "  HALTON_SEQUENCE computes the next N elements\n";
  cout << "  of a vector Halton sequence.\n";
  cout << "\n";
  cout << "  In this test, we compute the first 10 elements of a\n";
  cout << "  sequence, then change bases, reset the seed\n";
  cout << "  and recompute.\n";

  halton_dim_num_set ( DIM_NUM );
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i + 1 );
  }
  halton_base_set ( base );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << N    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  halton_sequence ( N, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  n = 10;
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = prime ( 2 * ( i + 1 ) );
  }
  halton_base_set ( base );

  cout << "\n";
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  halton_sequence ( N, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests U1_TO_SPHERE_UNIT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 1
# define DIM_NUM2 2

  double average[DIM_NUM2];
  double average_dot;
  double dot_average;
  int i;
  int i2;
  int j;
  int j2;
  int n;
  int seed[DIM_NUM];
  int step;
  double u[DIM_NUM];
  double v[DIM_NUM2];
  double x[DIM_NUM2];

  cout << "\n";
  cout << "TEST09\n";
  cout << "  For the unit sphere in 2 dimensions (the circle):\n";
  cout << "  HALTON generates \"U1\" points,\n";
  cout << "  U1_TO_SPHERE_UNIT_2D maps U2 points to the circle;\n";

  halton_dim_num_set ( DIM_NUM );
  n = 5;
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u1_to_sphere_unit_2d ( u, x );
    cout << "  ";
    for ( i2 = 0; i2 < DIM_NUM2; i2++ )
    {
      cout << setw(10) << x[i2] << "  ";
    }
    cout << "\n";
  }

  n = 1000;
  step = 0;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  for ( i2 = 0; i2 < DIM_NUM2; i2++ )
  {
    average[i2] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u1_to_sphere_unit_2d ( u, x );
    for ( i2 = 0; i2 < DIM_NUM2; i2++ )
    {
      average[i2] = average[i2] + x[i2];
    }
  }

  for ( i2 = 0; i2 < DIM_NUM2; i2++ )
  {
    average[i2] = average[i2] / ( ( double ) n );
  }

  cout << "\n";
  cout << "  Average the points, which should get a value\n";
  cout << "  close to zero, and closer as N increases.\n";
  cout << "\n";

  cout << "  Average:        ";
  for ( i2 = 0; i2 < DIM_NUM2; i2++ )
  {
    cout << setw(10) << average[i2] << "  ";
  }
  cout << "\n";

  cout << "\n";
  cout << "  Now choose a random direction, sample the same\n";
  cout << "  number of points, and compute the dot product with\n";
  cout << "  the direction.\n";
  cout << "  Take the absolute value of each dot product\n";
  cout << "  and sum and average.\n";
  cout << "\n";
  cout << "  We expect a value near 2 / PI = 0.6366...\n";

  for ( j2 = 0; j2 < 5; j2++ )
  {
    step = get_seed ( ) + 111 * j2;
    halton_step_set ( step );

    halton ( u );
    u1_to_sphere_unit_2d ( u, v );

    step = 0;
    halton_step_set ( step );

    dot_average = 0.0;

    for ( j = 0; j < n; j++ )
    {
      halton ( u );
      u1_to_sphere_unit_2d ( u, x );
      dot_average = dot_average + fabs ( r8vec_dot_product ( DIM_NUM2, x, v ) );
    }

    dot_average = dot_average / ( ( double ) n );

    cout << "\n";
    cout << "  Random V:         ";
    for ( i2 = 0; i2 < DIM_NUM2; i2++ )
    {
      cout << setw(10) << v[i2] << "  ";
    }
    cout << "\n";
    cout << "  Average |(XdotV)| "
         << setw(10) << dot_average << "\n";
  }

  return;
# undef DIM_NUM
# undef DIM_NUM2
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests U2_TO_BALL_UNIT_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2

  double average[DIM_NUM];
  double average_r;
  double average_theta;
  int i;
  int j;
  int n = 1000;
  double r;
  int seed[DIM_NUM];
  int step;
  double theta;
  double u[DIM_NUM];
  double v[DIM_NUM];
  double x[DIM_NUM];

  cout << "\n";
  cout << "TEST10\n";
  cout << "  For the unit ball in 2 dimensions (the disk):\n";
  cout << "  U2_TO_BALL_UNIT_2D samples;\n";

  halton_dim_num_set ( DIM_NUM );
  n = 5;
  step = 0;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u2_to_ball_unit_2d ( u, x );
    cout << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << x[i] << "  ";
    }
    cout << "\n";
  }

  n = 1000;
  step = 0;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u2_to_ball_unit_2d ( u, x );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      average[i] = average[i] + x[i];
    }
  }

  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = average[i] / ( ( double ) n );
  }

  cout << "\n";
  cout << "  Average the points, which should get a value\n";
  cout << "  close to zero, and closer as N increases.\n";
  cout << "\n";
  cout << "  Average:        ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << setw(10) << average[i] << "  ";
  }
  cout << "\n";

  step = 0;
  halton_step_set ( step );

  average_r = 0.0;

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u2_to_ball_unit_2d ( u, x );

    r = r8vec_norm_l2 ( DIM_NUM, x );

    average_r = average_r + r;
  }

  average_r = average_r / ( ( double ) n );

  cout << "\n";
  cout << "  Average the distance of the points from\n";
  cout << "  the center, which should be  \n";
  cout << "  DIM_NUM / ( DIM_NUM + 1 ) = "
    << ( ( double ) DIM_NUM ) / ( ( double ) ( DIM_NUM + 1 ) ) << "\n";
  cout << "\n";
  cout << "  Average:        " << average_r << "\n";

  step = 0;
  halton_step_set ( step );

  average_theta = 0.0;

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u2_to_ball_unit_2d ( u, x );
    theta = atan4 ( x[1], x[0] );
    average_theta = average_theta + theta;
  }

  average_theta = average_theta / ( ( double ) n );

  cout << "\n";
  cout << "  Average the angle THETA, which should approach PI.\n";
  cout << "\n";
  cout << "  Average:        " << average_theta << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests U2_TO_SPHERE_UNIT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define DIM_NUM2 3

  double average[DIM_NUM2];
  double dot_average;
  int i;
  int i2;
  int j;
  int j2;
  int n;
  int seed[DIM_NUM];
  int step;
  double u[DIM_NUM];
  double v[DIM_NUM2];
  double x[DIM_NUM2];

  cout << "\n";
  cout << "TEST11\n";
  cout << "  For the unit sphere in 3 dimensions:\n";
  cout << "  U2_TO_SPHERE_UNIT_3D samples;\n";

  halton_dim_num_set ( DIM_NUM );
  n = 5;
  step = 123456789;
  halton_step_set ( step );
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  halton_seed_set ( seed );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u2_to_sphere_unit_3d ( u, x );
    cout << "  ";
    for ( i2 = 0; i2 < DIM_NUM2; i2++ )
    {
      cout << setw(10) << x[i2] << "  ";
    }
    cout << "\n";
  }

  n = 1000;
  step = 0;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  for ( i = 0; i < DIM_NUM2; i++ )
  {
    average[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u2_to_sphere_unit_3d ( u, x );
    for ( i2 = 0; i2 < DIM_NUM2; i2++ )
    {
      average[i2] = average[i2] + x[i2];
    }
  }

  for ( i2 = 0; i2 < DIM_NUM2; i2++ )
  {
    average[i2] = average[i2] / ( ( double ) n );
  }

  cout << "\n";
  cout << "  Average the points, which should get a value\n";
  cout << "  close to zero, and closer as N increases.\n";
  cout << "\n";
  cout << "  Average:        ";
  for ( i2 = 0; i2 < DIM_NUM2; i2++ )
  {
    cout << setw(10) << average[i2] << "  ";
  }
  cout << "\n";

  cout << "\n";
  cout << "  Now choose a random direction, sample the same\n";
  cout << "  number of points, and compute the dot product with\n";
  cout << "  the direction.\n";
  cout << "  Take the absolute value of each dot product\n";
  cout << "  and sum and average.\n";

  for ( j2 = 0; j2 < 5; j2++ )
  {
    step = get_seed ( ) + 111 * j2;
    halton_step_set ( step );

    halton ( u );
    u2_to_sphere_unit_3d ( u, v );

    step = 0;
    halton_step_set ( step );

    dot_average = 0.0;

    for ( j = 0; j < n; j++ )
    {
      halton ( u );
      u2_to_sphere_unit_3d ( u, x );
      dot_average = dot_average + fabs ( r8vec_dot_product ( DIM_NUM2, x, v ) );
    }

    dot_average = dot_average / ( ( double ) n );

    cout << "\n";
    cout << "  Random V:         ";
    for ( i2 = 0; i2 < DIM_NUM2; i2++ )
    {
      cout << setw(10) << v[i2] << "  ";
    }
    cout << "\n";
    cout << "  Average |(XdotV)| "
         << setw(10) << dot_average << "\n";

  }

  return;
# undef DIM_NUM
# undef DIM_NUM2
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests U3_TO_BALL_UNIT_3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  double average[DIM_NUM];
  double average_phi;
  double average_r;
  double average_theta;
  int i;
  int j;
  int n;
  double phi;
  double r;
  int seed[DIM_NUM];
  int step;
  double theta;
  double u[DIM_NUM];
  double v[DIM_NUM];
  double x[DIM_NUM];

  cout << "\n";
  cout << "TEST12\n";
  cout << "  For the unit ball in 3 dimensions:\n";
  cout << "  U3_TO_BALL_UNIT_3D samples;\n";

  halton_dim_num_set ( DIM_NUM );
  n = 5;
  step = 0;
  halton_step_set ( step );

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << n    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  cout << "\n";
  cout << "  A few sample values:\n";
  cout << "\n";

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u3_to_ball_unit_3d ( u, x );
    cout << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(10) << x[i] << "  ";
    }
    cout << "\n";
  }

  n = 1000;
  step = 0;
  halton_step_set ( step );

  cout << "\n";
  cout << "  N =    " << setw(12) << n << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";

  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u3_to_ball_unit_3d ( u, x );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      average[i] = average[i] + x[i];
    }
  }

  for ( i = 0; i < DIM_NUM; i++ )
  {
    average[i] = average[i] / ( ( double ) n );
  }

  cout << "\n";
  cout << "  Average the points, which should get a value\n";
  cout << "  close to zero, and closer as N increases.\n";

  cout << "\n";
  cout << "  Average:        ";
  for ( i = 0; i < DIM_NUM; i++ )
  {
    cout << setw(10) << average[i] << "  ";
  }
  cout << "\n";

  step = 0;
  halton_step_set ( step );

  average_r = 0.0;

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u3_to_ball_unit_3d ( u, x );
    r = r8vec_norm_l2 ( DIM_NUM, x );
    average_r = average_r + r;
  }

  average_r = average_r / ( ( double ) n );

  cout << "\n";
  cout << "  Average the distance of the points from\n";
  cout << "  the center, which should be  \n";
  cout << "  DIM_NUM / ( DIM_NUM + 1 ) = "
    << ( ( double ) DIM_NUM ) / ( ( double ) ( DIM_NUM + 1 ) ) << "\n";
  cout << "\n";
  cout << "  Average:        " << average_r << "\n";

  step = 0;
  halton_step_set ( step );

  average_theta = 0.0;

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u3_to_ball_unit_3d ( u, x );
    theta = atan4 ( x[1], x[0] );
    average_theta = average_theta + theta;
  }

  average_theta = average_theta / ( ( double ) n );

  cout << "\n";
  cout << "  Average the angle THETA,\n";
  cout << "  which should approach PI.\n";
  cout << "\n";
  cout << "  Average:        " << average_theta << "\n";

  step = 0;
  halton_step_set ( step );

  average_phi = 0.0;

  for ( j = 0; j < n; j++ )
  {
    halton ( u );
    u3_to_ball_unit_3d ( u, x );
    r = r8vec_norm_l2 ( DIM_NUM, x );
    if ( r == 0.0 )
    {
      phi = 0.0;
    }
    else
    {
      phi = acos ( x[2] / r );
    }
    average_phi = average_phi + phi;
  }

  average_phi = average_phi / ( ( double ) n );

  cout << "\n";
  cout << "  Average the angle PHI,\n";
  cout << "  which should approach PI/2.\n";
  cout << "\n";
  cout << "  Average:        " << average_phi << "\n";

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests HALHAM_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2008
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 3

  int base[DIM_NUM];
  char file_name[81] = "halton_03_00010.txt";
  int i;
  int j;
  int leap[DIM_NUM];
  double r[DIM_NUM*N];
  int seed[DIM_NUM];
  int step;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  HALHAM_WRITE writes a Halton or Hammersley dataset to a file\n";

  step = 0;
  for ( i = 0; i < DIM_NUM; i++ )
  {
    seed[i] = 0;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    leap[i] = 1;
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    base[i] = prime ( i + 1 );
  }

  cout << "\n";
  cout << "  DIM_NUM = " << setw(12) << DIM_NUM << "\n";
  cout << "  N =    " << setw(12) << N    << "\n";
  cout << "  STEP = " << setw(12) << step << "\n";
  i4vec_transpose_print ( DIM_NUM, seed, "  SEED = " );
  i4vec_transpose_print ( DIM_NUM, leap, "  LEAP = " );
  i4vec_transpose_print ( DIM_NUM, base, "  BASE = " );

  i4_to_halton_sequence ( DIM_NUM, N, step, seed, leap, base, r );

  cout << "\n";
  cout << "    STEP   Halton\n";
  cout << "\n";
  for ( j = 0; j < N; j++ )
  {
    cout                        << "  "
         << setw(6) << step + j << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << setw(12) << r[i+j*DIM_NUM] << "  ";
    }
    cout << "\n";
  }

  halham_write ( DIM_NUM, N, step, seed, leap, base, r, file_name );

  cout << "\n";
  cout << "  The data was written to \"" << file_name << "\".\n";

  return;
# undef N
# undef DIM_NUM
}
