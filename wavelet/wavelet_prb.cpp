# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "wavelet.hpp"

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

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WAVELET_PRB.
//
//  Discussion:
//
//    WAVELET_PRB tests the WAVELET library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "WAVELET_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the WAVELET library.\n";

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
//
//  Terminate.
//
  cout << "\n";
  cout << "WAVELET_PRB\n";
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
//    TEST01 tests DAUB2_TRANSFORM and DAUB2_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  DAUB2_TRANSFORM computes the DAUB2 transform of a vector.\n";
  cout << "  DAUB2_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub2_transform ( n, u );

  w = daub2_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D2(U)    D2inv(D2(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub2_transform ( n, u );

  w = daub2_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D2(U)    D2inv(D2(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub2_transform ( n, u );

  w = daub2_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D2(U)    D2inv(D2(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub2_transform ( n, u );

  w = daub2_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D2(U)    D2inv(D2(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests DAUB4_TRANSFORM and DAUB4_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  DAUB4_TRANSFORM computes the DAUB4 transform of a vector.\n";
  cout << "  DAUB4_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub4_transform ( n, u );

  w = daub4_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D4(U)    D4inv(D4(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub4_transform ( n, u );

  w = daub4_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D4(U)    D4inv(D4(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub4_transform ( n, u );

  w = daub4_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D4(U)    D4inv(D4(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub4_transform ( n, u );

  w = daub4_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D4(U)    D4inv(D4(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests DAUB6_TRANSFORM and DAUB6_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  DAUB6_TRANSFORM computes the DAUB6 transform of a vector.\n";
  cout << "  DAUB6_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub6_transform ( n, u );

  w = daub6_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D6(U)    D6inv(D6(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub6_transform ( n, u );

  w = daub6_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D6(U)    D6inv(D6(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub6_transform ( n, u );

  w = daub6_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D6(U)    D6inv(D6(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub6_transform ( n, u );

  w = daub6_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D6(U)    D6inv(D6(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests DAUB8_TRANSFORM and DAUB8_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  DAUB8_TRANSFORM computes the DAUB8 transform of a vector.\n";
  cout << "  DAUB8_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub8_transform ( n, u );

  w = daub8_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D8(U)    D8inv(D8(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub8_transform ( n, u );

  w = daub8_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D8(U)    D8inv(D8(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub8_transform ( n, u );

  w = daub8_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D8(U)    D8inv(D8(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub8_transform ( n, u );

  w = daub8_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D8(U)    D8inv(D8(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests DAUB10_TRANSFORM and DAUB10_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  DAUB10_TRANSFORM computes the DAUB10 transform of a vector.\n";
  cout << "  DAUB10_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub10_transform ( n, u );

  w = daub10_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D10(U)    D10inv(D10(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub10_transform ( n, u );

  w = daub10_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D10(U)    D10inv(D10(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub10_transform ( n, u );

  w = daub10_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D10(U)    D10inv(D10(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub10_transform ( n, u );

  w = daub10_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D10(U)    D10inv(D10(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests DAUB12_TRANSFORM and DAUB12_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  DAUB12_TRANSFORM computes the DAUB12 transform of a vector.\n";
  cout << "  DAUB12_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub12_transform ( n, u );

  w = daub12_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D12(U)    D12inv(D12(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub12_transform ( n, u );

  w = daub12_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D12(U)    D12inv(D12(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub12_transform ( n, u );

  w = daub12_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D12(U)    D12inv(D12(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub12_transform ( n, u );

  w = daub12_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D12(U)    D12inv(D12(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests DAUB14_TRANSFORM and DAUB14_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  DAUB14_TRANSFORM computes the DAUB14 transform of a vector.\n";
  cout << "  DAUB14_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub14_transform ( n, u );

  w = daub14_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D14(U)    D14inv(D14(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub14_transform ( n, u );

  w = daub14_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D14(U)    D14inv(D14(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub14_transform ( n, u );

  w = daub14_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D14(U)    D14inv(D14(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub14_transform ( n, u );

  w = daub14_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D14(U)    D14inv(D14(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests DAUB16_TRANSFORM and DAUB16_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  DAUB16_TRANSFORM computes the DAUB16 transform of a vector.\n";
  cout << "  DAUB16_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub16_transform ( n, u );

  w = daub16_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D16(U)    D16inv(D16(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub16_transform ( n, u );

  w = daub16_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D16(U)    D16inv(D16(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub16_transform ( n, u );

  w = daub16_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D16(U)    D16inv(D16(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub16_transform ( n, u );

  w = daub16_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D16(U)    D16inv(D16(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests DAUB18_TRANSFORM and DAUB18_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  DAUB18_TRANSFORM computes the DAUB18 transform of a vector.\n";
  cout << "  DAUB18_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub18_transform ( n, u );

  w = daub18_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D18(U)    D18inv(D18(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub18_transform ( n, u );

  w = daub18_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D18(U)    D18inv(D18(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub18_transform ( n, u );

  w = daub18_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D18(U)    D18inv(D18(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub18_transform ( n, u );

  w = daub18_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D18(U)    D18inv(D18(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests DAUB20_TRANSFORM and DAUB20_TRANSFORM_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a_first;
  double a_last;
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  DAUB20_TRANSFORM computes the DAUB20 transform of a vector.\n";
  cout << "  DAUB20_TRANSFORM_INVERSE inverts it.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;

  u = r8vec_uniform_01_new ( n, seed );

  v = daub20_transform ( n, u );

  w = daub20_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D20(U)    D20inv(D20(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = new double[8];
  for ( i = 0; i < n; i++ )
  {
    u[i] = 1.0;
  }

  v = daub20_transform ( n, u );

  w = daub20_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D20(U)    D20inv(D20(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  a_first = 1.0;
  a_last = ( double ) ( n );
  u = r8vec_linspace_new ( n, a_first, a_last );

  v = daub20_transform ( n, u );

  w = daub20_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D20(U)    D20inv(D20(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = new double[n];
  for ( i = 0; i < n; i++ )
  {
    u[i] = pow ( ( double ) ( i - 5 ), 2 );
  }

  v = daub20_transform ( n, u );

  w = daub20_transform_inverse ( n, v );

  cout << "\n";
  cout << "   i      U          D20(U)    D20inv(D20(U))\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
