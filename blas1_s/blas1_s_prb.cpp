# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "blas1_s.hpp"

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
float r4_uniform_01 ( int *seed );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLAS1_S_PRB.
//
//  Discussion:
//
//    BLAS1_S_PRB tests the BLAS1 single precision real routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "BLAS1_S_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the BLAS1_S library.\n";

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
  cout << "BLAS1_S_PRB:\n";
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
//    TEST01 demonstrates ISAMAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int i1;
  int incx;
  float x[N];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  ISAMAX returns the index of maximum magnitude;\n";

  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = ( float ) ( ( 7 * i ) % 11 ) - ( float ) ( N / 2 );
  }

  cout << "\n";
  cout << "  The vector X:\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(8) << x[i-1] << "\n";
  }

  incx = 1;

  i1 = isamax ( N, x, incx );

  cout << "\n";
  cout << "  The index of maximum magnitude = " << i1 << "\n";

  return;
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests ISAMAX, SAXPY and SSCAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define LDA N

  float a[LDA*N];
  float b[N];
  int i;
  int info;
  int ipvt[N];
  int j;
  int k;
  int l;
  float t;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Use ISAMAX, SAXPY and SSCAL\n";
  cout << "  in a Gauss elimination routine.\n";
//
//  Set the matrix.
//
  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      if ( i == j )
      {
        a[i-1+(j-1)*LDA] = 2.0;
      }
      else if ( i == j + 1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else if ( i == j - 1 )
      {
        a[i-1+(j-1)*LDA] = -1.0;
      }
      else
      {
        a[i-1+(j-1)*LDA] = 0.0;
      }
    }
  }
//
//  Set the right hand side.
//
  for ( i = 1; i <= N-1; i++ )
  {
    b[i-1] = 0.0;
  }
  b[N-1] = ( float ) ( N + 1 );

  info = 0;

  for ( k = 1; k <= N-1; k++ )
  {
    l = isamax ( N-k+1, a+(k-1)+(k-1)*LDA, 1 ) + k - 1;
    ipvt[k-1] = l;

    if ( a[l-1+(k-1)*LDA] == 0.0 )
    {
      info = k;
    }
    else
    {
      if ( l != k )
      {
        t = a[l-1+(k-1)*LDA];
        a[l-1+(k-1)*LDA] = a[k-1+(k-1)*LDA];
        a[k-1+(k-1)*LDA] = t;
      }

      t = -1.0 / a[k-1+(k-1)*LDA];
      sscal ( N-k, t, a+k+(k-1)*LDA, 1 );

      for ( j = k+1; j <= N; j++ )
      {
        t = a[l-1+(j-1)*LDA];
        if ( l != k )
        {
          a[l-1+(j-1)*LDA] = a[k-1+(j-1)*LDA];
          a[k-1+(j-1)*LDA] = t;
        }
        saxpy ( N-k, t, a+k+(k-1)*LDA, 1, a+k+(j-1)*LDA, 1 );
      }
    }
  }

  ipvt[N-1] = N;
  if ( a[N-1+(N-1)*LDA] == 0.0 )
  {
    info = N;
  }

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  The matrix is singular.\n";
    return;
  }

  for ( k = 1; k <= N-1; k++ )
  {
    l = ipvt[k-1];
    t = b[l-1];
    if ( l != k )
    {
      b[l-1] = b[k-1];
      b[k-1] = t;
    }
    saxpy ( N-k, t, a+k+(k-1)*LDA, 1, b+k, 1 );
  }

  for ( k = N; 1 <= k; k-- )
  {
    b[k-1] = b[k-1] / a[k-1+(k-1)*LDA];
    t = -b[k-1];
    saxpy ( k-1, t, a+0+(k-1)*LDA, 1, b, 1 );
  }

  cout << "\n";
  cout << "  First five entries of solution:\n";
  cout << "\n";
  for ( i = 1; i <= 5; i++ )
  {
    cout << "  " << setw(14) << b[i-1];
  }
  cout << "\n";

  return;
# undef LDA
# undef N
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests SASUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define LDA 6
# define MA 5
# define NA 4
# define NX 10

  float a[LDA*NA];
  int i;
  int j;
  float x[NX];

  for ( i = 0; i < NX; i++ )
  {
    x[i] = pow ( -1.0, i + 1 ) * ( float ) ( 2 * ( i + 1 ) );
  }

  cout << "\n";
  cout << "TEST03\n";
  cout << "  SASUM adds the absolute values of elements of a vector.\n";
  cout << "\n";
  cout << "  X = \n";
  cout << "\n";
  for ( i = 0; i < NX; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "\n";
  }

  cout << "\n";
  cout << "  SASUM ( NX,   X, 1 ) =    " << sasum ( NX,   x, 1 )    << "\n";
  cout << "  SASUM ( NX/2, X, 2 ) =    " << sasum ( NX/2, x, 2 )    << "\n";
  cout << "  SASUM ( 2,    X, NX/2 ) = " << sasum ( 2,    x, NX/2 ) << "\n";

  for ( i = 0; i < MA; i++ )
  {
    for ( j = 0; j < NA; j++ )
    {
      a[i+j*LDA] = pow ( -1.0, i + 1 + j + 1)
        * ( float ) ( 10 * ( i + 1 ) + j + 1 );
    }
  }

  cout << "\n";
  cout << "  Demonstrate with a matrix A:\n";
  cout << "\n";
  for ( i = 0; i < MA; i++ )
  {
    for ( j = 0; j < NA; j++ )
    {
      cout << "  " << setw(14) << a[i+j*LDA];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  SASUM(MA,A(1,2),1) =   " << sasum ( MA, a+0+1*LDA, 1 )   << "\n";
  cout << "  SASUM(NA,A(2,1),LDA) = " << sasum ( NA, a+1+0*LDA, LDA ) << "\n";

  return;
# undef LDA
# undef MA
# undef NA
# undef NX
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests SAXPY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  float da;
  int i;
  float x[N];
  float y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  cout << "\n";
  cout << "TEST04\n";
  cout << "  SAXPY adds a multiple of vector X to vector Y.\n";
  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "\n";
  }
  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << y[i]  << "\n";
  }

  da = 1.0;
  saxpy ( N, da, x, 1, y, 1 );
  cout << "\n";
  cout << "  SAXPY ( N, " << da << ", X, 1, Y, 1 )\n";
  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << y[i]  << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  da = -2.0;
  saxpy ( N, da, x, 1, y, 1 );
  cout << "\n";
  cout << "  SAXPY ( N, " << da << ", X, 1, Y, 1 )\n";
  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << y[i]  << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  da = 3.0;
  saxpy ( 3, da, x, 2, y, 1 );
  cout << "\n";
  cout << "  SAXPY ( 3, " << da << ", X, 2, Y, 1 )\n";
  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << y[i]  << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  da = -4.0;
  saxpy ( 3, da, x, 1, y, 2 );
  cout << "\n";
  cout << "  SAXPY ( 3, " << da << ", X, 1, Y, 2 )\n";
  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << y[i]  << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 demonstrates SCOPY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
  float a[5*5];
  int i;
  int j;
  float x[10];
  float y[10];

  cout << "\n";
  cout << "TEST05\n";
  cout << "  SCOPY copies one vector into another.\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < 10; i++ )
  {
    y[i] = ( float ) ( 10 * ( i + 1 ) );
  }

  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      a[i+j*5] = ( float ) ( 10 * ( i + 1 ) + j + 1 );
    }
  }

  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "\n";
  }
  cout << "\n";
  cout << "  Y =\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << y[i]  << "\n";
  }
  cout << "\n";
  cout << "  A =\n";
  cout << "\n";
  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      cout << "  " << setw(14) << a[i+j*5];
    }
      cout << "\n";
  }

  scopy ( 5, x, 1, y, 1 );
  cout << "\n";
  cout << "  SCOPY ( 5, X, 1, Y, 1 )\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << y[i]  << "\n";
  }

  for ( i = 0; i < 10; i++ )
  {
    y[i] = ( float ) ( 10 * ( i + 1 ) );
  }

  scopy ( 3, x, 2, y, 3 );
  cout << "\n";
  cout << "  SCOPY ( 3, X, 2, Y, 3 )\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << y[i]  << "\n";
  }

  scopy ( 5, x, 1, a, 1 );
  cout << "\n";
  cout << "  SCOPY ( 5, X, 1, A, 1 )\n";
  cout << "\n";
  cout << "  A =\n";
  cout << "\n";
  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      cout << "  " << setw(14) << a[i+j*5];
    }
      cout << "\n";
  }

  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      a[i+j*5] = ( float ) ( 10 * ( i + 1 ) + j + 1 );
    }
  }

  scopy ( 5, x, 2, a, 5 );
  cout << "\n";
  cout << "  SCOPY ( 5, X, 2, A, 5 )\n";
  cout << "\n";
  cout << "  A =\n";
  cout << "\n";
  for ( i = 0; i < 5; i++ )
  {
    for ( j = 0; j < 5; j++ )
    {
      cout << "  " << setw(14) << a[i+j*5];
    }
      cout << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 demonstrates SDOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define LDA 10
# define LDB 7
# define LDC 6

  float a[LDA*LDA];
  float b[LDB*LDB];
  float c[LDC*LDC];
  int i;
  int j;
  float sum1;
  float x[N];
  float y[N];

  cout << "\n";
  cout << "TEST06\n";
  cout << "  SDOT computes the dot product of vectors.\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = - ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*LDA] = ( float ) ( i + 1 + j + 1 );
    }
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      b[i+j*LDB] = ( float ) ( ( i + 1 ) - ( j + 1 ) );
    }
  }

  sum1 = sdot ( N, x, 1, y, 1 );

  cout << "\n";
  cout << "  Dot product of X and Y is " << sum1 << "\n";
//
//  To multiply a ROW of a matrix A times a vector X, we need to
//  specify the increment between successive entries of the row of A:
//
  sum1 = sdot ( N, a+1+0*LDA, LDA, x, 1 );

  cout << "\n";
  cout << "  Product of row 2 of A and X is " << sum1 << "\n";
//
//  Product of a column of A and a vector is simpler:
//
  sum1 = sdot ( N, a+0+1*LDA, 1, x, 1 );

  cout << "\n";
  cout << "  Product of column 2 of A and X is " << sum1 << "\n";
//
//  Here's how matrix multiplication, c = a*b, could be done
//  with SDOT:
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      c[i+j*LDC] = sdot ( N, a+i, LDA, b+0+j*LDB, 1 );
    }
  }

  cout << "\n";
  cout << "  Matrix product computed with SDOT:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(14) << c[i+j*LDC];
    }
    cout << "\n";
  }

  return;
# undef N
# undef LDA
# undef LDB
# undef LDC
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 demonstrates SMACH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
  int job;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  SMACH returns some approximate machine numbers.\n";
  cout << "\n";
  job = 1;
  cout << "  SMACH(1) = EPS =  " << smach ( job ) << "\n";
  job = 2;
  cout << "  SMACH(2) = TINY = " << smach ( job ) << "\n";
  job = 3;
  cout << "  SMACH(3) = HUGE = " << smach ( job ) << "\n";

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 demonstrates SNRM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define LDA 10
//
//  These parameters illustrate the fact that matrices are typically
//  dimensioned with more space than the user requires.
//
  float a[LDA*LDA];
  int i;
  int incx;
  int j;
  float sum1;
  float x[N];

  cout << "\n";
  cout << "TEST08\n";
  cout << "  SNRM2 computes the Euclidean norm of a vector.\n";
  cout << "\n";
//
//  Compute the euclidean norm of a vector:
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "\n";
  }
  cout << "\n";
  cout << "  The 2-norm of X is " << snrm2 ( N, x, 1 ) << "\n";
//
//  Compute the euclidean norm of a row or column of a matrix:
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*LDA] = ( float ) ( i + 1 + j + 1 );
    }
  }

  cout << "\n";
  cout << "  The 2-norm of row 2 of A is "
       << snrm2 ( N, a+1+0*LDA, LDA ) << "\n";

  cout << "\n";
  cout << "  The 2-norm of column 2 of A is "
       << snrm2 ( N, a+0+1*LDA, 1 ) << "\n";

  return;
# undef N
# undef LDA
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests SROT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  float c;
  int i;
  float s;
  float x[N];
  float y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( ( i + 1 ) * ( i + 1 ) - 12 );
  }

  cout << "\n";
  cout << "TEST09\n";
  cout << "  SROT carries out a Givens rotation.\n";
  cout << "\n";
  cout << "  X and Y\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "  "
         << setw(14) << y[i]  << "\n";
  }

  c = 0.5;
  s = sqrt ( 1.0 - c * c );
  srot ( N, x, 1, y, 1, c, s );
  cout << "\n";
  cout << "  SROT ( N, X, 1, Y, 1, " << c << "," << s << " )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "  "
         << setw(14) << y[i]  << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( ( i + 1 ) * ( i + 1 ) - 12 );
  }

  c = x[0] / sqrt ( x[0] * x[0] + y[0] * y[0] );
  s = y[0] / sqrt ( x[0] * x[0] + y[0] * y[0] );
  srot ( N, x, 1, y, 1, c, s );
  cout << "\n";
  cout << "  SROT ( N, X, 1, Y, 1, " << c << "," << s << " )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "  "
         << setw(14) << y[i]  << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests SROTG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  float a;
  float b;
  float c;
  float r;
  float s;
  float sa;
  float sb;
  int seed;
  int test;
  int test_num = 5;
  float z;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  SROTG generates a real Givens rotation\n";
  cout << "    (  C  S ) * ( A ) = ( R )\n";
  cout << "    ( -S  C )   ( B )   ( 0 )\n";
  cout << "\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r4_uniform_01 ( &seed );
    b = r4_uniform_01 ( &seed );

    sa = a;
    sb = b;

    srotg ( &sa, &sb, &c, &s );

    r = sa;
    z = sb;

    cout << "\n";
    cout << "  A =  " << a << "  B =  " << b << "\n";
    cout << "  C =  " << c << "  S =  " << s << "\n";
    cout << "  R =  " << r << "  Z =  " << z << "\n";
    cout << "   C*A+S*B = " <<  c * a + s * b << "\n";
    cout << "  -S*A+C*B = " << -s * a + c * b << "\n";
  }

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests SSCAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  float da;
  int i;
  float x[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  cout << "\n";
  cout << "TEST11\n";
  cout << "  SSCAL multiplies a vector by a scalar.\n";
  cout << "\n";
  cout << "  X =\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "\n";
  }

  da = 5.0;
  sscal ( N, da, x, 1 );
  cout << "\n";
  cout << "  SSCAL ( N, " << da << ", X, 1 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  da = -2.0;
  sscal ( 3, da, x, 2 );
  cout << "\n";
  cout << "  SSCAL ( 3, " << da << ", X, 2 )\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests SSWAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int i;
  float x[N];
  float y[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  cout << "\n";
  cout << "TEST12\n";
  cout << "  SSWAP swaps two vectors.\n";
  cout << "\n";
  cout << "  X and Y:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "  "
         << setw(14) << y[i]  << "\n";
  }

  sswap ( N, x, 1, y, 1 );
  cout << "\n";
  cout << "  SSWAP ( N, X, 1, Y, 1 )\n";
  cout << "\n";
  cout << "  X and Y:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "  "
         << setw(14) << y[i]  << "\n";
  }

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( float ) ( i + 1 );
  }

  for ( i = 0; i < N; i++ )
  {
    y[i] = ( float ) ( 100 * ( i + 1 ) );
  }

  sswap ( 3, x, 2, y, 1 );
  cout << "\n";
  cout << "  SSWAP ( 3, X, 2, Y, 1 )\n";

  cout << "\n";
  cout << "  X and Y:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i + 1 << "  "
         << setw(14) << x[i]  << "  "
         << setw(14) << y[i]  << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

float r4_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_01 returns a real pseudorandom number.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r4_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    P A Lewis, A S Goodman, J M Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  float r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( float ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
