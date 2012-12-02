# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <complex>

using namespace std;

# include "normal.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for NORMAL_PRB.
//
//  Discussion:
//
//    NORMAL_PRB calls sample problems for the NORMAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "NORMAL_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the NORMAL library\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "NORMAL_PRB\n";
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
//    TEST01 tests C4_NORMAL_01.
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
  cout << "  C4_NORMAL_01 computes pseudorandom complex values\n";
  cout << "  normally distributed in the unit circle.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = c4_normal_01 ( seed );

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
//    TEST02 tests C8_NORMAL_01.
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
  complex <double> value;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  C8_NORMAL_01 computes pseudorandom complex values\n";
  cout << "  normally distributed in the unit circle.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = c8_normal_01 ( seed );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << real ( value )
         << "  " << setw(12) << imag ( value ) << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests I4_NORMAL.
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
  cout << "TEST03\n";
  cout << "  I4_NORMAL computes pseudonormal integer values\n";
  cout << "  with mean A and standard deviation B.\n";

  a = 70.0;
  b = 10.0;
  seed = seed_init;

  cout << "\n";
  cout << "  The mean A = " << a << "\n";
  cout << "  The standard deviation B = " << b << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(8) << i4_normal ( a, b, seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests I8_NORMAL.
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
  double a;
  double b;
  long long int i;
  long long int seed;
  long long int seed_init = 12345678987654321LL;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  I8_NORMAL computes pseudonormal integer values\n";
  cout << "  with mean A and standard deviation B.\n";

  a = 70.0;
  b = 10.0;
  seed = seed_init;

  cout << "\n";
  cout << "  The mean A = " << a << "\n";
  cout << "  The standard deviation B = " << b << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(8) << i8_normal ( a, b, seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests R4_NORMAL.
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
  cout << "TEST05\n";
  cout << "  R4_NORMAL computes pseudonormal values\n";
  cout << "  with mean A and standard deviation B.\n";

  a = 10.0;
  b = 2.0;
  seed = seed_init;

  cout << "\n";
  cout << "  The mean A = " << a << "\n";
  cout << "  The standard deviation B = " << b << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r4_normal ( a, b, seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests R4_NORMAL_01.
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
  cout << "TEST06\n";
  cout << "  R4_NORMAL_01 computes pseudonormal values\n";
  cout << "  with mean 0 and standard deviation 1.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r4_normal_01 ( seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests R8_NORMAL.
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
  cout << "TEST07\n";
  cout << "  R8_NORMAL computes pseudonormal values\n";
  cout << "  with mean A and standard deviation B.\n";

  a = 10.0;
  b = 2.0;
  seed = seed_init;

  cout << "\n";
  cout << "  The mean A = " << a << "\n";
  cout << "  The standard deviation B = " << b << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r8_normal ( a, b, seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests R8_NORMAL_01.
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
  cout << "TEST08\n";
  cout << "  R8_NORMAL_01 computes pseudonormal values\n";
  cout << "  with mean 0 and standard deviation 1.\n";

  seed = seed_init;

  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(6)  << i
         << "  " << setw(14) << r8_normal_01 ( seed ) << "\n";
  }

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests R8_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  int seed_init = 123456789;
  int seed_input;
  int seed_output;
  double value;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  R8_NORMAL_01 computes pseudonormal values\n";
  cout << "  with mean 0 and standard deviation 1.\n";
  cout << "\n";
  cout << "  Verify that we can change the seed\n";
  cout << "  and get the desired results.\n";
  cout << "\n";
  cout << "  The initial seed is " << seed_init << "\n";

  seed = seed_init;

  cout << "\n";
  cout << "         I    Seed(in)   Seed(out)   R8_NORMAL_01\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    seed_input = seed;
    value = r8_normal_01 ( seed );
    seed_output = seed;

    cout << "  " << setw(8)  << i
         << "  " << setw(12) << seed_input
         << "  " << setw(12) << seed_output
         << "  " << setw(14) << value << "\n";
  }

  seed = seed_init;

  cout << "\n";
  cout << "  Resetting seed to repeat, after an ODD number of steps.\n";
  cout << "\n";

  for ( i = 6; i <= 10; i++ )
  {
    seed_input = seed;
    value = r8_normal_01 ( seed );
    seed_output = seed;

    cout << "  " << setw(8)  << i
         << "  " << setw(12) << seed_input
         << "  " << setw(12) << seed_output
         << "  " << setw(14) << value << "\n";
  }

  seed = seed_init;

  cout << "\n";
  cout << "  Resetting seed to repeat, after an EVEN number of steps.\n";
  cout << "\n";

  for ( i = 11; i <= 15; i++ )
  {
    seed_input = seed;
    value = r8_normal_01 ( seed );
    seed_output = seed;

    cout << "  " << setw(8)  << i
         << "  " << setw(12) << seed_input
         << "  " << setw(12) << seed_output
         << "  " << setw(14) << value << "\n";
  }

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//   TEST10 tests R8_NORMAL_01;
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
  cout << "TEST10\n";
  cout << "  R8_NORMAL_01 computes a sequence of normally distributed\n";
  cout << "  pseudorandom numbers.\n";

  seed = 12345;

  cout << "\n";
  cout << "  Initial SEED = " << seed << "\n";

  cout << "\n";
  cout << "  First 10 values:\n";
  cout << "\n";
  cout << "       I         Input        Output   R8_NORMAL_01\n";
  cout << "                  SEED          SEED\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    seed_in = seed;
    u[i] = r8_normal_01 ( seed );
    seed_out = seed;

    cout                         << "  "
         << setw(6)  << i + 1    << "  " 
         << setw(12) << seed_in  << "  " 
         << setw(12) << seed_out << "  " 
         << setw(10) << u[i]     << "\n";
  }

  cout << "\n";
  cout << "  Now call R8_NORMAL_01 " << N << " times.\n";

  u_avg = 0.0;
  for ( i = 0; i < N; i++ )
  {
    u[i] = r8_normal_01 ( seed );
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
  cout << "  Expecting       " << 0.0        << "\n";

  cout << "\n";
  cout << "  Variance =      " << u_var      << "\n";
  cout << "  Expecting       " << 1.0        << "\n";

  return;
# undef N
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests R8_NORMAL_01 and R8MAT_NORMAL_01_NEW.
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
  cout << "TEST11\n";
  cout << "  R8_NORMAL_01 computes pseudorandom normal values one at a time.\n";
  cout << "  R8MAT_NORMAL_01_NEW computes a matrix of values.\n";
  cout << "\n";
  cout << "  For the same initial seed, the results should be identical,\n";
  cout << "  but R8MAT_NORMAL_01 might be faster.\n";
  cout << "\n";
  cout << "  Initial seed is " << seed_init << "\n";

  seed = seed_init;
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < M; i++ )
    {
      a[i+j*M] = r8_normal_01 ( seed );
    }
  }

  seed = seed_init;
  b = r8mat_normal_01_new ( M, N, seed );

  cout << "\n";
  cout << "      I       J      A[I,J]        B[I,J]\n";
  cout << "                 (R8_NORMAL_01)  (R8MAT_NORMAL_01_NEW)\n";
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

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//   TEST12 tests R8_NORMAL_01 and R8VEC_NORMAL_01_NEW.
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
# define N 10

  double a[N];
  double *b;
  int i;
  int j;
  int k;
  int seed;
  int seed_init = 123456789;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  R8_NORMAL_01 computes pseudeorandom normal values one at a time.\n";
  cout << "  R8VEC_NORMAL_01_NEW computes a vector of values.\n";
  cout << "\n";
  cout << "  For the same initial seed, the results should be identical,\n";
  cout << "  but R8VEC_NORMAL_01_NEW might be faster.\n";
  cout << "\n";
  cout << "  Initial seed is " << seed_init << "\n";

  seed = seed_init;
  for ( i = 0; i < N; i++ )
  {
    a[i] = r8_normal_01 ( seed );
  }

  seed = seed_init;
  b = r8vec_normal_01_new ( N, seed );

  cout << "\n";
  cout << "      I      A[I]          B[I]\n";
  cout << "         (R8_NORMAL_01)  (R8VEC_NORMAL_01_NEW)\n";
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


