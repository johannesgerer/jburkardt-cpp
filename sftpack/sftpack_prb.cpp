# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "sftpack.hpp"

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

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SFTPACK_PRB.
//
//  Discussion:
//
//    SFTPACK_PRB tests the SFTPACK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//   22 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SFTPACK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SFTPACK library.\n";

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
//
//  Terminate.
//
  cout << "\n";
  cout << "SFTPACK_PRB\n";
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
//    TEST01 tests R8VEC_SCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double ahi = 5.0;
  double alo = 0.0;
  double *c;
  double *d;
  double *e;
  int i;
  int n = 256;
  int seed;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For slow cosine transforms,\n";
  cout << "  R8VEC_SCT does a forward or backward transform.\n";
  cout << "\n";
  cout << "  The number of data items is N = " << n << "\n";
//
//  Set the data values.
//
  seed = 123456789;

  c = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, c, 10, "  The original data:" );
//
//  Compute the coefficients.
//
  d = r8vec_sct ( n, c );

  r8vec_print_part ( n, d, 10, "  The cosine coefficients:" );
//
//  Now compute inverse transform of coefficients.  Should get back the
//  original data.
//
  e = r8vec_sct ( n, d );

  for ( i = 0; i < n; i++ )
  {
    e[i] = e[i] / ( double ) ( 2 * n );
  }

  r8vec_print_part ( n, e, 10, "  The retrieved data:" );

  delete [] c;
  delete [] d;
  delete [] e;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R8VEC_SFTB and R8VEC_SFTF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double ahi = 5.0;
  double alo = 0.0;
  double azero;
  double *b;
  int i;
  int n = 36;
  int seed;
  double *x;
  double *z;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For real slow Fourier transforms,\n";
  cout << "  R8VEC_SFTF computes the forward transform.\n";
  cout << "  R8VEC_SFTB computes the backward transform.\n";
  cout << "\n";
  cout << "  The number of data values, N = " << n << "\n";

  seed = 123456789;

  x = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, x, 10, "  The original data:" );
//
//  Compute the slow Fourier transform of the data.
//
  a = new double[n/2];
  b = new double[n/2];

  r8vec_sftf ( n, x, &azero, a, b );

  cout << "\n";
  cout << "  A (cosine) coefficients:\n";
  cout << "\n";

  cout << "  " << setw(4) << 0
       << "  " << setw(14) << azero << "\n";

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << a[i] << "\n";
  }

  cout << "\n";
  cout << "  B (sine) coefficients:\n";
  cout << "\n";

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << b[i] << "\n";
  }
//
//  Now try to retrieve the data from the coefficients.
//
  z = r8vec_sftb ( n, azero, a, b );

  r8vec_print_part ( n, z, 10, "  The retrieved data:" );

  delete [] a;
  delete [] b;
  delete [] x;
  delete [] z;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R8VEC_SHT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double ahi = 5.0;
  double alo = 0.0;
  double *c;
  double *d;
  double *e;
  int n = 17;
  int seed;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For real slow Hartley transforms,\n";
  cout << "  R8VEC_SHT does a forward or backward transform.\n";
  cout << "\n";
  cout << "  The number of data items is N = " << n << "\n";
//
//  Set the data values.
//
  seed = 123456789;

  c = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, c, 10, "  The original data:" );
//
//  Compute the coefficients.
//
  d = r8vec_sht ( n, c );

  r8vec_print_part ( n, d, 10, "  The Hartley coefficients:" );
//
//  Now compute inverse transform of coefficients.  Should get back the
//  original data.
//
  e = r8vec_sht ( n, d );

  r8vec_print_part ( n, e, 10, "  The retrieved data:" );

  delete [] c;
  delete [] d;
  delete [] e;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests R8VEC_SQCTB and R8VEC_SQCTF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double ahi = 5.0;
  double alo = 0.0;
  int n = 256;
  int seed;
  double *x;
  double *y;
  double *z;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For real slow quarter wave cosine transforms,\n";
  cout << "  R8VEC_SQCTF does a forward transform;\n";
  cout << "  R8VEC_SQCTB does a backward transform.\n";
  cout << "\n";
  cout << "  The number of data items is N = " << n << "\n";
//
//  Set the data values.
//
  seed = 123456789;

  x = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, x, 10, "  The original data:" );
//
//  Compute the coefficients.
//
  y = r8vec_sqctf ( n, x );

  r8vec_print_part ( n, y, 10, "  The cosine coefficients:" );
//
//  Now compute inverse transform of coefficients.  Should get back the
//  original data.
//
  z = r8vec_sqctb ( n, y );

  r8vec_print_part ( n, z, 10, "  The retrieved data:" );

  delete [] x;
  delete [] y;
  delete [] z;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests R8VEC_SQSTB and R8VEC_SQSTF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double ahi = 5.0;
  double alo = 0.0;
  int n = 256;
  int seed;
  double *x;
  double *y;
  double *z;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For real slow quarter wave sine transforms,\n";
  cout << "  R8VEC_SQSTF does a forward transform;\n";
  cout << "  R8VEC_SQSTB does a backward transform.\n";
  cout << "\n";
  cout << "  The number of data items is N = " << n << "\n";
//
//  Set the data values.
//
  seed = 123456789;

  x = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, x, 10, "  The original data:" );
//
//  Compute the coefficients.
//
  y = r8vec_sqstf ( n, x );

  r8vec_print_part ( n, y, 10, "  The sine coefficients:" );
//
//  Now compute inverse transform of coefficients.  Should get back the
//  original data.
//
  z = r8vec_sqstb ( n, y );

  r8vec_print_part ( n, z, 10, "  The retrieved data:" );

  delete [] x;
  delete [] y;
  delete [] z;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests R8VEC_SST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double ahi = 5.0;
  double alo = 0.0;
  double *c;
  double *d;
  double *e;
  int i;
  int n = 256;
  int seed;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  For slow sine transforms,\n";
  cout << "  R8VEC_SST does a forward or backward transform.\n";
  cout << "\n";
  cout << "  The number of data items is N = " << n << "\n";
//
//  Set the data values.
//
  seed = 123456789;

  c = r8vec_uniform_new ( n, alo, ahi, &seed );

  r8vec_print_part ( n, c, 10, "  The original data:" );
//
//  Compute the coefficients;
//
  d = r8vec_sst ( n, c );

  r8vec_print_part ( n, d, 10, "  The sine coefficients:" );
//
//  Now compute inverse transform of coefficients.  Should get back the
//  original data.
//
  e = r8vec_sst ( n, d );

  for ( i = 0; i < n; i++ )
  {
    e[i] = e[i] / ( double ) ( 2 * ( n + 1 ) );
  }

  r8vec_print_part ( n, e, 10, "  The retrieved data:" );

  delete [] c;
  delete [] d;
  delete [] e;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests C4VEC_SFTB and C4VEC_SFTF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 36;
  int seed;
  complex <float> *x;
  complex <float> *x2;
  complex <float> *y;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  For complex slow Fourier transforms,\n";
  cout << "  C4VEC_SFTF computes the forward transform.\n";
  cout << "  C4VEC_SFTB computes the backward transform.\n";
  cout << "\n";
  cout << "  The number of data values, N = " << n << "\n";

  seed = 123456789;

  x = c4vec_uniform_01_new ( n, &seed );

  c4vec_print_part ( n, x, 10, "  The original data:" );
//
//  Compute the slow Fourier transform of the data.
//
  y = c4vec_sftf ( n, x );

  c4vec_print_part ( n, y, 10, "  The Fourier coefficients:" );
//
//  Now try to retrieve the data from the coefficients.
//
  x2 = c4vec_sftb ( n, y );

  c4vec_print_part ( n, x2, 10, "  The retrieved data:" );

  delete [] x;
  delete [] x2;
  delete [] y;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests C8VEC_SFTB and C8VEC_SFTF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 36;
  int seed;
  complex <double> *x;
  complex <double> *x2;
  complex <double> *y;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  For complex slow Fourier transforms,\n";
  cout << "  C8VEC_SFTF computes the forward transform.\n";
  cout << "  C8VEC_SFTB computes the backward transform.\n";
  cout << "\n";
  cout << "  The number of data values, N = " << n << "\n";

  seed = 123456789;

  x = c8vec_uniform_01_new ( n, &seed );

  c8vec_print_part ( n, x, 10, "  The original data:" );
//
//  Compute the slow Fourier transform of the data.
//
  y = c8vec_sftf ( n, x );

  c8vec_print_part ( n, y, 10, "  The Fourier coefficients:" );
//
//  Now try to retrieve the data from the coefficients.
//
  x2 = c8vec_sftb ( n, y );

  c8vec_print_part ( n, x2, 10, "  The retrieved data:" );

  delete [] x;
  delete [] x2;
  delete [] y;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests R4VEC_SFTB and R4VEC_SFTF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  float *a;
  float ahi = 5.0;
  float alo = 0.0;
  float azero;
  float *b;
  int i;
  int n = 36;
  int seed;
  float *x;
  float *z;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  For real slow Fourier transforms,\n";
  cout << "  R4VEC_SFTF computes the forward transform.\n";
  cout << "  R4VEC_SFTB computes the backward transform.\n";
  cout << "\n";
  cout << "  The number of data values, N = " << n << "\n";

  seed = 123456789;

  x = r4vec_uniform_new ( n, alo, ahi, &seed );

  r4vec_print_part ( n, x, 10, "  The original data:" );
//
//  Compute the slow Fourier transform of the data.
//
  a = new float[n/2];
  b = new float[n/2];

  r4vec_sftf ( n, x, &azero, a, b );

  cout << "\n";
  cout << "  A (cosine) coefficients:\n";
  cout << "\n";

  cout << "  " << setw(4) << 0
       << "  " << setw(14) << azero << "\n";

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << a[i] << "\n";
  }

  cout << "\n";
  cout << "  B (sine) coefficients:\n";
  cout << "\n";

  for ( i = 0; i < ( n / 2 ); i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << b[i] << "\n";
  }
//
//  Now try to retrieve the data from the coefficients.
//
  z = r4vec_sftb ( n, azero, a, b );

  r4vec_print_part ( n, z, 10, "  The retrieved data:" );

  delete [] a;
  delete [] b;
  delete [] x;
  delete [] z;

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests C4MAT_SFTB and C4MAT_SFTF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n1 = 10;
  int n2 = 4;
  int seed;
  complex <float> *x;
  complex <float> *x2;
  complex <float> *y;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  For complex slow Fourier transforms,\n";
  cout << "  C4MAT_SFTF computes the forward transform.\n";
  cout << "  C4MAT_SFTB computes the backward transform.\n";
  cout << "\n";
  cout << "  The data has dimensions N1 = " << n1 << ", N2 = " << n2 << "\n";

  seed = 123456789;

  x = c4mat_uniform_01_new ( n1, n2, &seed );

  c4mat_print_some ( n1, n2, x, 1, 1, 10, 10, "  The original data:" );
//
//  Compute the slow Fourier transform of the data.
//
  y = c4mat_sftf ( n1, n2, x );

  c4mat_print_some ( n1, n2, y, 1, 1, 10, 10, "  The Fourier coefficients:" );
//
//  Now try to retrieve the data from the coefficients.
//
  x2 = c4mat_sftb ( n1, n2, y );

  c4mat_print_some ( n1, n2, x2, 1, 1, 10, 10, "  The retrieved data:" );

  delete [] x;
  delete [] x2;
  delete [] y;

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests C8MAT_SFTB and C8MAT_SFTF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n1 = 10;
  int n2 = 4;
  int seed;
  complex <double> *x;
  complex <double> *x2;
  complex <double> *y;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  For complex slow Fourier transforms,\n";
  cout << "  C8MAT_SFTF computes the forward transform.\n";
  cout << "  C8MAT_SFTB computes the backward transform.\n";
  cout << "\n";
  cout << "  The data has dimensions N1 = " << n1 << ", N2 = " << n2 << "\n";

  seed = 123456789;

  x = c8mat_uniform_01_new ( n1, n2, &seed );

  c8mat_print_some ( n1, n2, x, 1, 1, 10, 10, "  The original data:" );
//
//  Compute the slow Fourier transform of the data.
//
  y = c8mat_sftf ( n1, n2, x );

  c8mat_print_some ( n1, n2, y, 1, 1, 10, 10, "  The Fourier coefficients:" );
//
//  Now try to retrieve the data from the coefficients.
//
  x2 = c8mat_sftb ( n1, n2, y );

  c8mat_print_some ( n1, n2, x2, 1, 1, 10, 10, "  The retrieved data:" );

  delete [] x;
  delete [] x2;
  delete [] y;

  return;
}
