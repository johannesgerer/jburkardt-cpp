# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "c4lib.hpp"

int main ( );

void test0061 ( );
void test0062 ( );
void test0063 ( );
void test0064 ( );
void test0065 ( );
void test0066 ( );
void test0067 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for C4LIB_PRB.
//
//  Discussion:
//
//    C4LIB_PRB tests the C4LIB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "C4LIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the C4LIB library.\n";

  test0061 ( );
  test0062 ( );
  test0063 ( );
  test0064 ( );
  test0065 ( );
  test0066 ( );
  test0067 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "C4LIB_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test0061 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0061 tests C4_ARG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  float argument;
  complex <float> c1;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST0061\n";
  cout << "  C4_ARG computes the argument of a C4.\n";
  cout << "\n";
  cout << 
    "            C1=random            ARG=C4_ARG(C1)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c4_uniform_01 ( &seed );
    argument = c4_arg ( c1 );

    cout << "  (" << setw(10) << real ( c1 )
         << ",  " << setw(10) << imag ( c1 ) << ")"
         << "  " << setw(10) << argument << "\n";
  }

  return;
}
//****************************************************************************80

void test0062 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0062 tests C4_CUBE_ROOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> c1;
  complex <float> c2;
  complex <float> c3;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST0062\n";
  cout << "  C4_CUBE_ROOT computes the principal cube root of a C4.\n";
  cout << "\n";
  cout << 
    "            C1=random            C2=C4_CUBE_ROOT(C1)         C3=C2*C2*C2\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c4_uniform_01 ( &seed );
    c2 = c4_cube_root ( c1 );
    c3 = c2 * c2 * c2;

    cout << "  (" << setw(10) << real ( c1 )
         << ",  " << setw(10) << imag ( c1 ) << ")"
         << "  (" << setw(10) << real ( c2 )
         << ",  " << setw(10) << imag ( c2 ) << ")"
         << "  (" << setw(10) << real ( c3 )
         << ",  " << setw(10) << imag ( c3 ) << ")\n";
  }

  return;
}
//****************************************************************************80

void test0063 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0063 tests C4_MAG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> c1;
  float magnitude;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST0063\n";
  cout << "  C4_MAG computes the magnitude of a C4.\n";
  cout << "\n";
  cout << 
    "            C1=random            MAG=C4_MAG(C1)\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c4_uniform_01 ( &seed );
    magnitude = c4_mag ( c1 );

    cout << "  (" << setw(10) << real ( c1 )
         << ",  " << setw(10) << imag ( c1 ) << ")"
         << "  " << setw(10) << magnitude << "\n";
  }

  return;
}
//****************************************************************************80

void test0064 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0064 tests C4_NORMAL_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int seed = 123456789;
  int test;
  int test_num = 20;
  complex <float> x;

  cout << "\n";
  cout << "TEST0064\n";
  cout << "  C4_NORMAL_01 generates unit pseudonormal\n";
  cout << "    complex values.\n";
  cout << "  Using initial random number seed = " << seed << "\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    x = c4_normal_01 ( &seed );
    cout << "  " << setw(10) << real ( x ) 
         << "  " << setw(10) << imag ( x ) << "\n";
  }

  return;
}
//****************************************************************************80

void test0065 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0065 tests C4_SQRT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> c1;
  complex <float> c2;
  complex <float> c3;
  int seed = 123456789;
  int test;
  int test_num = 10;

  cout << "\n";
  cout << "TEST0065\n";
  cout << "  C4_SQRT computes the principal square root of a C4.\n";
  cout << "\n";
  cout << 
    "            C1=random            C2=C4_SQRT(C1)              C3=C2*C2\n";
  cout << "\n";

  for ( test = 0; test < test_num; test++ )
  {
    c1 = c4_uniform_01 ( &seed );
    c2 = c4_sqrt ( c1 );
    c3 = c2 * c2;

    cout << "  (" << setw(10) << real ( c1 )
         << ",  " << setw(10) << imag ( c1 ) << ")"
         << "  (" << setw(10) << real ( c2 )
         << ",  " << setw(10) << imag ( c2 ) << ")"
         << "  (" << setw(10) << real ( c3 )
         << ",  " << setw(10) << imag ( c3 ) << ")\n";
  }

  return;
}
//****************************************************************************80

void test0066 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0066 tests C4MAT_UNIFORM_01_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> *a;
  int m = 5;
  int n = 4;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST0066\n";
  cout << "  C4MAT_UNIFORM_01_NEW computes a random complex matrix.\n";

  a = c4mat_uniform_01_new ( m, n, &seed );

  c4mat_print ( m, n, a, "  The matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void test0067 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0067 tests C4VEC_INDICATOR_NEW;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  complex <float> *a;
  int n = 10;

  cout << "\n";
  cout << "TEST0067\n";
  cout << "  C4VEC_INDICATOR_NEW sets A = (1-1i,2-2i,...,N-Ni)\n";

  a = c4vec_indicator_new ( n );
 
  c4vec_print ( n, a, "  The indicator vector:" );

  delete [] a;

  return;
}
