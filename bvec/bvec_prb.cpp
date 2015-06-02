# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "bvec.hpp"

int main ( );

void bvec_add_test ( );
void bvec_complement2_test ( );
void bvec_mul_test ( );
void bvec_next_test ( );
void bvec_next_grlex_test ( );
void bvec_print_test ( );
void bvec_sub_test ( );
void bvec_to_i4_test ( );
void bvec_uniform_new_test ( );
void i4_bclr_test ( );
void i4_bset_test ( );
void i4_btest_test ( );
void i4_to_bvec_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BVEC_PRB.
//
//  Discussion:
//
//    BVEC_PRB tests the BVEC library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BVEC_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BVEC library.\n";

  bvec_add_test ( );
  bvec_complement2_test ( );
  bvec_mul_test ( );
  bvec_next_test ( );
  bvec_next_grlex_test ( );
  bvec_print_test ( );
  bvec_sub_test ( );
  bvec_to_i4_test ( );
  bvec_uniform_new_test ( );
  i4_bclr_test ( );
  i4_bset_test ( );
  i4_btest_test ( );
  i4_to_bvec_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BVEC_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void bvec_add_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_ADD_TEST tests BVEC_ADD;
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
# define N 10

  int bvec1[N];
  int bvec2[N];
  int bvec3[N];
  int bvec4[N];
  int i;
  int j;
  int k;
  int l;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "BVEC_ADD_TEST\n";
  cout << "  BVEC_ADD adds binary vectors representing integers;\n";
  cout << "\n";
  cout << "        I        J        I + J   BVEC_ADD\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform_ab ( -100, 100, seed );
    j = i4_uniform_ab ( -100, 100, seed );

    k = i + j;

    i4_to_bvec ( i, N, bvec1 );
    i4_to_bvec ( j, N, bvec2 );
    bvec_add ( N, bvec1, bvec2, bvec3 );
    l = bvec_to_i4 ( N, bvec3 );

    cout << "  " << setw(8) << i
         << "  " << setw(8) << j
         << "  " << setw(8) << k
         << "  " << setw(8) << l << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void bvec_complement2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_COMPLEMENT2_TEST tests BVEC_COMPLEMENT2;
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
# define N 10

  int bvec1[N];
  int bvec2[N];
  int i;
  int j;
  int seed = 123456789;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "BVEC_COMPLEMENT2_TEST\n";
  cout << "  BVEC_COMPLEMENT2 returns the two's complement\n";
  cout << "  of a (signed) binary vector;\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( -100, 100, seed );

    i4_to_bvec ( i, N, bvec1 );

    bvec_complement2 ( N, bvec1, bvec2 );

    j = bvec_to_i4 ( N, bvec2 );

    cout << "\n";
    cout << "  I = " << "  " << i << "\n";
    cout << "  J = " << "  " << j << "\n";
    bvec_print ( N, bvec1, "" );
    bvec_print ( N, bvec2, "" );
  }

  return;
# undef N
}
//****************************************************************************80

void bvec_mul_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_MUL_TEST tests BVEC_MUL;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 January 2015
//
//  Author:
//
//    John Burkardt
//
{
# define N 15

  int bvec1[N];
  int bvec2[N];
  int bvec3[N];
  int i;
  int j;
  int k;
  int l;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "BVEC_MUL_TEST\n";
  cout << "  BVEC_MUL multiplies binary vectors\n";
  cout << "  representing integers;\n";
  cout << "\n";
  cout << "        I        J        I * J   BVEC_MUL\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform_ab ( -100, 100, seed );
    j = i4_uniform_ab ( -100, 100, seed );

    k = i * j;

    i4_to_bvec ( i, N, bvec1 );
    i4_to_bvec ( j, N, bvec2 );
    bvec_mul ( N, bvec1, bvec2, bvec3 );
    l = bvec_to_i4 ( N, bvec3 );

    cout << "  " << setw(8) << i
         << "  " << setw(8) << j
         << "  " << setw(8) << k
         << "  " << setw(8) << l << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void bvec_next_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_NEXT_TEST tests BVEC_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 January 2015
//
//  Author:
//
//    John Burkardt
//
{ 
  int *b;
  int i;
  int n = 4;

  cout << "\n";
  cout << "BVEC_NEXT_TEST\n";
  cout << "  BVEC_NEXT computes the 'next' BVEC.\n";
  cout << "\n";

  b = new int[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0;
  }

  for ( i = 0; i <= 16; i++ )
  {
    bvec_print ( n, b, "" );
    bvec_next ( n, b );
  }

  delete [] b;

  return;
}
//****************************************************************************80

void bvec_next_grlex_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_NEXT_GRLEX_TEST tests BVEC_NEXT_GRLEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2015
//
//  Author:
//
//    John Burkardt
//
{ 
  int *b;
  int i;
  int j;
  int n = 4;

  cout << "\n";
  cout << "BVEC_NEXT_GRLEX_TEST\n";
  cout << "  BVEC_NEXT_GRLEX computes binary vectors in GRLEX order.\n";
  cout << "\n";

  b = new int[n];

  for ( j = 0; j < n; j++ )
  {
    b[j] = 0;
  }

  for ( i = 0; i <= 16; i++ )
  {
    cout << "  " << setw(2) << i << ":  ";
    for ( j = 0; j < n; j++ )
    {
      cout << b[j];
    }
    cout << "\n";
    bvec_next_grlex ( n, b );
  }

  delete [] b;

  return;
}
//****************************************************************************80

void bvec_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_PRINT_TEST tests BVEC_PRINT.
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
  int n = 10;
  int bvec[10] = { 1, 0, 0, 1, 0, 1, 1, 1, 0, 0 };

  cout << "\n";
  cout << "BVEC_PRINT_TEST\n";
  cout << "  BVEC_PRINT prints a binary vector.\n";

  bvec_print ( n, bvec, "  BVEC:" );

  return;
}
//****************************************************************************80

void bvec_sub_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_SUB_TEST tests BVEC_SUB;
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
# define N 10

  int bvec1[N];
  int bvec2[N];
  int bvec3[N];
  int bvec4[N];
  int i;
  int j;
  int k;
  int l;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "BVEC_SUB_TEST\n";
  cout << "  BVEC_SUB subtracts binary vectors representing integers;\n";
  cout << "\n";
  cout << "        I        J        I - J   BVEC_SUB\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  { 
    i = i4_uniform_ab ( -100, 100, seed );
    j = i4_uniform_ab ( -100, 100, seed );

    k = i - j;

    i4_to_bvec ( i, N, bvec1 );
    i4_to_bvec ( j, N, bvec2 );
    bvec_sub ( N, bvec1, bvec2, bvec4 );
    l = bvec_to_i4 ( N, bvec4 );

    cout << "  " << setw(8) << i
         << "  " << setw(8) << j
         << "  " << setw(8) << k
         << "  " << setw(8) << l << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void bvec_to_i4_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BVEC_TO_I4_TEST tests BVEC_TO_I4;
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
# define N 10

  int bvec[N];
  int i;
  int i2;
  int j;

  cout << "\n";
  cout << "BVEC_TO_I4_TEST\n";
  cout << "  BVEC_TO_I4 converts a signed binary vector\n";
  cout << "  to an integer;\n";
  cout << "\n";
  cout << "  I --> BVEC  -->  I\n";
  cout << "\n";

  for ( i = -3; i <= 10; i++ )
  {
    i4_to_bvec ( i, N, bvec );
    i2 = bvec_to_i4 ( N, bvec );

    cout << setw(3) << i << "  ";
    for ( j = 0; j < N; j++ )
    {
      cout << setw(1) << bvec[j];
    }
    cout << "  ";
    cout << setw(3) << i2 << "\n";
  }

  return;
# undef N
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

void i4_bclr_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_BCLR_TEST tests I4_BCLR.
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
# define TEST_NUM 2

  int i4;
  int i4_test[TEST_NUM] = { 101, -31 };
  int ivec[32];
  int j1;
  int pos;
  int test;

  cout << "\n";
  cout << "I4_BCLR_TEST\n";
  cout << "  I4_BCLR sets a given bit to 0.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4 = i4_test[test];

    i4_to_bvec ( i4, 32, ivec );

    cout << "\n";
    cout << "  Working on I4 = " << i4 << "\n";
    cout << "\n";
    cout << "       Pos     Digit       I4_BCLR\n";
    cout << "\n";

    for ( pos = 0; pos <= 31; pos++ )
    {
      j1 = i4_bclr ( i4, pos );

      cout << "  " << setw(8) << pos
           << "  " << setw(8) << ivec[pos]
           << "  " << j1 << "\n";
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void i4_bset_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_BSET_TEST tests I4_BSET.
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
# define TEST_NUM 2

  int i4;
  int i4_test[TEST_NUM] = { 101, -31 };
  int ivec[32];
  int j1;
  int pos;
  int test;

  cout << "\n";
  cout << "I4_BSET_TEST\n";
  cout << "  I4_BSET sets a given bit to 0.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4 = i4_test[test];

    i4_to_bvec ( i4, 32, ivec );

    cout << "\n";
    cout << "  Working on I4 = " << i4 << "\n";
    cout << "\n";
    cout << "       Pos     Digit       I4_BSET\n";
    cout << "\n";

    for ( pos = 0; pos <= 31; pos++ )
    {
      j1 = i4_bset ( i4, pos );

      cout << "  " << setw(8) << pos
           << "  " << setw(8) << ivec[pos]
           << "  " << j1 << "\n";
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void i4_btest_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_BTEST_TEST tests I4_BTEST.
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
# define TEST_NUM 2

  int i4;
  int i4_test[TEST_NUM] = { 101, -31 };
  int ivec[32];
  int j1;
  int pos;
  int test;

  cout << "\n";
  cout << "I4_BTEST_TEST\n";
  cout << "  I4_BTEST reports whether a given bit is 0 or 1.\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    i4 = i4_test[test];

    i4_to_bvec ( i4, 32, ivec );

    cout << "\n";
    cout << "  Analyze the integer I4 = " << i4 << "\n";
    cout << "\n";
    cout << "       Pos     Digit  I4_BTEST\n";
    cout << "\n";

    for ( pos = 0; pos <= 31; pos++ )
    {
      j1 = i4_btest ( i4, pos );

      cout << "  " << setw(8) << pos
           << "  " << setw(8) << ivec[pos]
           << "  " << j1 << "\n";
    }
  }

  return;
# undef TEST_NUM
}
//****************************************************************************80

void i4_to_bvec_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_BVEC_TEST tests I4_TO_BVEC;
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
# define N 10

  int bvec[N];
  int i;
  int i2;
  int j;

  cout << "\n";
  cout << "I4_TO_BVEC_TEST\n";
  cout << "  I4_TO_BVEC converts an integer to a \n";
  cout << "  signed binary vector;\n";
  cout << "\n";
  cout << "  I --> BVEC  -->  I\n";
  cout << "\n";

  for ( i = -3; i <= 10; i++ )
  {
    i4_to_bvec ( i, N, bvec );
    i2 = bvec_to_i4 ( N, bvec );

    cout << setw(3) << i << "  ";
    for ( j = 0; j < N; j++ )
    {
      cout << setw(1) << bvec[j];
    }
    cout << "  ";
    cout << setw(3) << i2 << "\n";
  }

  return;
# undef N
}
