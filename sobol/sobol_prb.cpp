# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "sobol.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test055 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SOBOL_PRB.
//
//  Discussion:
//
//    SOBOL_PRB tests the SOBOL library.
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
{
  timestamp ( );
  cout << "\n";
  cout << "SOBOL_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SOBOL library.\n";
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test055 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SOBOL_PRB\n";
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
//    TEST01 tests OR.
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
{
  int i;
  int j;
  int k;
  int seed;
  int test;

  seed = 123456789;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  The C++ function ^ computes the bitwise exclusive OR\n";
  cout << "  of two integers.\n";
  cout << "\n";
  cout << "       I       J     I^J\n";
  cout << "\n";

  for ( test = 1; test <= 10; test++ )
  {
    i = i4_uniform ( 0, 100, &seed );
    j = i4_uniform ( 0, 100, &seed );
    k = i ^ j;

    cout << "  "
         << setw(6) << i << "  "
         << setw(6) << j << "  "
         << setw(6) << k << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests I4_BIT_HI1.
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
{
  int i;
  int j;
  int k;
  int seed;
  int test;

  seed = 123456789;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  I4_BIT_HI1 returns the location of the high bit in an integer.\n";
  cout << "\n";
  cout << "       I  I4_BIT_HI1(I)\n";
  cout << "\n";

  for ( test = 1; test <= 10; test++ )
  {
    i = i4_uniform ( 0, 100, &seed );
    j = i4_bit_hi1 ( i );

    cout << "  "
         << setw(6) << i << "  "
         << setw(6) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests I4_BIT_LO0.
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
{
  int i;
  int j;
  int k;
  int seed;
  int test;

  seed = 123456789;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  I4_BIT_LO0 returns the location of the low 0 bit in an integer.\n";
  cout << "\n";
  cout << "       I  I4_BIT_LO0(I)\n";
  cout << "\n";

  for ( test = 1; test <= 10; test++ )
  {
    i = i4_uniform ( 0, 100, &seed );
    j = i4_bit_lo0 ( i );

    cout << "  "
         << setw(6) << i << "  "
         << setw(6) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests I4_SOBOL.
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
{
# define DIM_MAX 4

  int dim_num;
  int i;
  int j;
  float r[DIM_MAX];
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  I4_SOBOL computes the next element of a Sobol sequence.\n";
  cout << "\n";
  cout << "  In this test, we call I4_SOBOL repeatedly.\n";

  for ( dim_num = 2; dim_num <= DIM_MAX; dim_num++ )
  {

    seed = 0;

    cout << "\n";
    cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";
    cout << "\n";
    cout << "  Seed  Seed   I4_SOBOL\n";
    cout << "  In    Out\n";
    cout << "\n";

    for ( i = 0; i <= 110; i++ )
    {
      seed_in = seed;
      i4_sobol ( dim_num, &seed, r );
      seed_out = seed;

      if ( i <= 11 || 95 <= i )
      {
        cout << setw(6) << seed_in << "  "
             << setw(6) << seed_out << "  ";
        for ( j = 0; j < dim_num; j++ )
        {
          cout << setw(14) << r[j] << "  ";
        }
        cout << "\n";
      }
      else if ( i == 12 )
      {
        cout << "....................\n";
      }
    }

  }

  return;
# undef DIM_MAX
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests I4_SOBOL.
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
{
# define DIM_NUM 3

  int i;
  int j;
  float r[DIM_NUM];
  int seed;
  int seed_in;
  int seed_out;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  I4_SOBOL computes the next element of a Sobol sequence.\n";
  cout << "\n";
  cout << "  In this test, we demonstrate how the SEED can be\n";
  cout << "  manipulated to skip ahead in the sequence, or\n";
  cout << "  to come back to any part of the sequence.\n";
  cout << "\n";
  cout << "  Using dimension DIM_NUM =   " << DIM_NUM << "\n";

  seed = 0;

  cout << "\n";
  cout << "  Seed  Seed   I4_SOBOL\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 11; i++ )
  {
    seed_in = seed;
    i4_sobol ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(14) << r[j] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump ahead by increasing SEED:\n";

  seed = 100;

  cout << "\n";
  cout << "  Seed  Seed   I4_SOBOL\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    seed_in = seed;
    i4_sobol ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(14) << r[j] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump back by decreasing SEED:\n";

  seed = 3;

  cout << "\n";
  cout << "  Seed  Seed   I4_SOBOL\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 11; i++ )
  {
    seed_in = seed;
    i4_sobol ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(14) << r[j] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump ahead by increasing SEED:\n";

  seed = 98;

  cout << "\n";
  cout << "  Seed  Seed   I4_SOBOL\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    seed_in = seed;
    i4_sobol ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(14) << r[j] << "  ";
    }
    cout << "\n";
  }
  return;
# undef DIM_NUM
}
//****************************************************************************80

void test055 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST055 tests OR on long long ints.
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
{
  long long int i;
  long long int j;
  long long int k;
  int seed;
  int test;

  seed = 123456789;

  cout << "\n";
  cout << "TEST055\n";
  cout << "  The C++ function ^ computes the bitwise exclusive OR\n";
  cout << "  of two LONG LONG INT's.\n";
  cout << "\n";
  cout << "       I       J     I^J\n";
  cout << "\n";

  for ( test = 1; test <= 10; test++ )
  {
    i = ( long long int ) i4_uniform ( 0, 100, &seed );
    j = ( long long int ) i4_uniform ( 0, 100, &seed );
    k = i ^ j;

    cout << "  "
         << setw(6) << i << "  "
         << setw(6) << j << "  "
         << setw(6) << k << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests I8_BIT_HI1.
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
{
  long long int i;
  int j;
  int seed;
  int test;

  seed = 123456789;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  I8_BIT_HI1 returns the location of the high bit in an integer.\n";
  cout << "\n";
  cout << "       I  I8_BIT_HI1(I)\n";
  cout << "\n";

  for ( test = 1; test <= 10; test++ )
  {
    i = ( long long int ) i4_uniform ( 0, 100, &seed );
    j = i8_bit_hi1 ( i );

    cout << "  "
         << setw(6) << i << "  "
         << setw(6) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests I8_BIT_LO0.
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
{
  long long int i;
  int j;
  int k;
  int seed;
  int test;

  seed = 123456789;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  I8_BIT_LO0 returns the location of the low zero bit\n";
  cout << "  in an integer.\n";
  cout << "\n";
  cout << "       I  I8_BIT_LO0(I)\n";
  cout << "\n";

  for ( test = 1; test <= 10; test++ )
  {
    i = ( long long int ) i4_uniform ( 0, 100, &seed );
    j = i8_bit_lo0 ( i );

    cout << "  "
         << setw(6) << i << "  "
         << setw(6) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests I8_SOBOL.
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
{
# define DIM_MAX 4

  int dim_num;
  int i;
  int j;
  double r[DIM_MAX];
  long long int seed;
  long long int seed_in;
  long long int seed_out;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  I8_SOBOL computes the next element of a Sobol sequence.\n";
  cout << "\n";
  cout << "  In this test, we call I8_SOBOL repeatedly.\n";

  for ( dim_num = 2; dim_num <= DIM_MAX; dim_num++ )
  {

    seed = 0;

    cout << "\n";
    cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";
    cout << "\n";
    cout << "  Seed  Seed   I8_SOBOL\n";
    cout << "  In    Out\n";
    cout << "\n";

    for ( i = 0; i <= 110; i++ )
    {
      seed_in = seed;
      i8_sobol ( dim_num, &seed, r );
      seed_out = seed;

      if ( i <= 11 || 95 <= i )
      {
        cout << setw(6) << seed_in << "  "
             << setw(6) << seed_out << "  ";
        for ( j = 0; j < dim_num; j++ )
        {
          cout << setw(14) << r[j] << "  ";
        }
        cout << "\n";
      }
      else if ( i == 12 )
      {
        cout << "....................\n";
      }
    }

  }

  return;
# undef DIM_MAX
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests I8_SOBOL.
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
{
# define DIM_NUM 3

  int i;
  int j;
  double r[DIM_NUM];
  long long int seed;
  long long int seed_in;
  long long int seed_out;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  I8_SOBOL computes the next element of a Sobol sequence.\n";
  cout << "\n";
  cout << "  In this test, we demonstrate how the SEED can be\n";
  cout << "  manipulated to skip ahead in the sequence, or\n";
  cout << "  to come back to any part of the sequence.\n";
  cout << "\n";
  cout << "  Using dimension DIM_NUM =   " << DIM_NUM << "\n";

  seed = 0;

  cout << "\n";
  cout << "  Seed  Seed   I8_SOBOL\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 11; i++ )
  {
    seed_in = seed;
    i8_sobol ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(14) << r[j] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump ahead by increasing SEED:\n";

  seed = 100;

  cout << "\n";
  cout << "  Seed  Seed   I8_SOBOL\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    seed_in = seed;
    i8_sobol ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(14) << r[j] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump back by decreasing SEED:\n";

  seed = 3;

  cout << "\n";
  cout << "  Seed  Seed   I8_SOBOL\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 11; i++ )
  {
    seed_in = seed;
    i8_sobol ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(14) << r[j] << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  Jump ahead by increasing SEED:\n";

  seed = 98;

  cout << "\n";
  cout << "  Seed  Seed   I8_SOBOL\n";
  cout << "  In    Out\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    seed_in = seed;
    i8_sobol ( DIM_NUM, &seed, r );
    seed_out = seed;
    cout << setw(6) << seed_in << "  "
         << setw(6) << seed_out << "  ";
    for ( j = 0; j < DIM_NUM; j++ )
    {
      cout << setw(14) << r[j] << "  ";
    }
    cout << "\n";
  }
  return;
# undef DIM_NUM
}
