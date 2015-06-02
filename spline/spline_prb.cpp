# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "spline.hpp"

int main ( );
void test001 ( );
void test002 ( );
void test003 ( );
void test004 ( );
void test005 ( );
void test006 ( );

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
void test115 ( );
void test116 ( );
void test12 ( );
void test125 ( );
void test126 ( );
void test127 ( );
void test13 ( );
void test14 ( );
void test145 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );

void test20 ( );
void test205 ( );
void test21 ( );
void test215 ( );
void test22 ( );
void test225 ( );
void test23 ( );
void test235 ( );
void test24 ( );

double frunge ( double x );
double fprunge ( double x );
double fpprunge ( double x );
double fcube ( double x );
double fpcube ( double x );
double fppcube ( double x );
void parabola_formula ( double x, double *y, double *yp, double *ypp );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPLINE_PRB.
//
//  Discussion:
//
//    SPLINE_PRB tests the SPLINE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPLINE_PRB\n";
  cout << "  C++ version:\n";
  cout << "  Test the SPLINE library.\n";

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );

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
  test115 ( );
  test116 ( );
  test12 ( );
  test125 ( );
  test126 ( );
  test127 ( );
  test13 ( );
  test14 ( );
  test145 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test205 ( );
  test21 ( );
  test215 ( );
  test22 ( );
  test225 ( );
  test23 ( );
  test235 ( );
  test24 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPLINE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test001 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST001 tests PARABOLA_VAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2004
//
//  Author:
//
//    John Burkardt
//
# define NDIM 1
# define NDATA 5
{
  int i;
  int left;
  double xdata[NDATA];
  double xval;
  double ydata[NDIM*NDATA];
  double yval[NDIM];
  double zdata[NDATA];
  double zval[NDIM];

  cout << "\n";
  cout << "TEST001\n";
  cout << "  PARABOLA_VAL2 evaluates parabolas through\n";
  cout << "    3 points in a table\n";
  cout << "\n";
  cout << "  Our data tables will actually be parabolas:\n";
  cout << "    A: 2*x**2 + 3 * x + 1.\n";
  cout << "    B: 4*x**2 - 2 * x + 5.\n";
  cout << "\n";

  for ( i = 0; i < NDATA; i++ )
  {
    xval = 2.0 * ( double ) ( i + 1 );
    xdata[i] = xval;
    ydata[0+i*NDIM] = 2.0 * xval * xval + 3.0 * xval + 1.0;
    zdata[i] = 4.0 * xval * xval - 2.0 * xval + 5.0;
    cout << setw(6)  << i        << "  "
         << setw(10) << xdata[i] << "  "
         << setw(10) << ydata[i] << "  "
         << setw(10) << zdata[i] << "\n";
  }

  cout << "\n";
  cout << "  Interpolated data:\n";
  cout << "\n";
  cout << "  LEFT, X, Y1, Y2\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    xval = ( double ) ( 2 * i - 1 );
    left = i;
    if ( NDATA - 2 < left )
    {
      left = NDATA - 2;
    }
    if ( left < 1 )
    {
      left = 1;
    }
    parabola_val2 ( NDIM, NDATA, xdata, ydata, left, xval, yval );
    parabola_val2 ( NDIM, NDATA, xdata, zdata, left, xval, zval );

    cout                        << "  "
         << setw(6)  << left    << "  "
         << setw(10) << xval    << "  "
         << setw(10) << yval[0] << "  "
         << setw(10) << zval[0] << "\n";
  }

  return;
# undef NDATA
# undef NDIM
}
//****************************************************************************80

void test002 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST002 tests R8VEC_BRACKET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
# define N 10
# define NTEST 6
{
  int i;
  int itest;
  int left;
  int right;
  double x[N];
  double xtest[NTEST] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n" ;
  cout << "TEST002\n";
  cout << "  R8VEC_BRACKET finds a pair of entries in a\n";
  cout << "    sorted real array which bracket a value.\n";

  for ( i = 1; i <= N; i++ )
  {
    x[i-1] = ( double ) i;
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  Sorted array:" );

  for ( itest = 0; itest < NTEST; itest++ )
  {
    xval = xtest[itest];

    cout << "\n";
    cout << "  Search for XVAL = " << xval << "\n";

    r8vec_bracket ( N, x, xval, &left, &right );

    cout << "  X["
         << setw(2)  << left      << "-1] = "
         << setw(12) << x[left-1] << "\n";

    cout << "  X["
         << setw(2)  << right      << "-1] = "
         << setw(12) << x[right-1] << "\n";

  }

  return;

# undef N
# undef NTEST
}
//****************************************************************************80

void test003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST003 tests R8VEC_BRACKET3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
# define N 10
# define NTEST 6
{
  int i;
  int itest;
  int left;
  double x[N];
  double xtest[NTEST] = { -10.0, 1.0, 4.5, 5.0, 10.0, 12.0 };
  double xval;

  cout << "\n";
  cout << "TEST003\n";
  cout << "  R8VEC_BRACKET3 finds a pair of entries in a\n";
  cout << "    sorted real array which bracket a value.\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  x[5] = x[4];

  r8vec_print ( N, x, "  Sorted array:" );

  left = ( N + 1 ) / 2;

  for ( itest = 0; itest < NTEST; itest++ )
  {
    xval = xtest[itest];

    cout << "\n";
    cout << "  Search for XVAL = " << xval << "\n";

    cout << "  Starting guess for interval is = " << left << "\n";

    r8vec_bracket3 ( N, x, xval, &left );

    cout << "  Nearest interval:\n";
    cout << "   X[" << left   << "-1]= " << x[left-1] << "\n";
    cout << "   X[" << left+1 << "-1]= " << x[left  ] << "\n";

  }

  return;
# undef N
# undef NTEST
}
//****************************************************************************80

void test004 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST004 tests R8VEC_ORDER_TYPE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
# define N 4
# define NTEST 6
{
  int itest;
  int j;
  int order;
  double x[N];

  cout << "\n";
  cout << "TEST004\n";
  cout << "  R8VEC_ORDER_TYPE classifies a real vector as\n";
  cout << "  -1: no order\n";
  cout << "   0: all equal;\n";
  cout << "   1: ascending;\n";
  cout << "   2: strictly ascending;\n";
  cout << "   3: descending;\n";
  cout << "   4: strictly descending.\n";
  cout << "\n";

  for ( itest = 1; itest <= NTEST; itest++ )
  {
    if ( itest == 1 )
    {
      x[0] = 1.0;
      x[1] = 3.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 2 )
    {
      x[0] = 2.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 2.0;
    }
    else if ( itest == 3 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 2.0;
      x[3] = 4.0;
    }
    else if ( itest == 4 )
    {
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 3.0;
      x[3] = 4.0;
    }
    else if ( itest == 5 )
    {
      x[0] = 4.0;
      x[1] = 4.0;
      x[2] = 3.0;
      x[3] = 1.0;
    }
    else if ( itest == 6 )
    {
      x[0] = 9.0;
      x[1] = 7.0;
      x[2] = 3.0;
      x[3] = 0.0;
    }

    order = r8vec_order_type ( N, x );

    cout << "\n";
    cout << "  Vector of order type " << order << ":\n";
    cout << "\n";

    for ( j = 0; j < N; j++ )
    {
      cout                     << "  "
           << setw(6)  << j    << "  "
           << setw(10) << x[j] << "\n";
    }

  }

  return;
# undef N
# undef NTEST
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests D3_NP_FS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
# define N 10
{
  double *a;
  double *b;
  int i;
  int seed;
  double *x;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  D3_NP_FS factors and solves a tridiagonal\n";
  cout << "    linear system.\n";
//
//  Set the matrix.
//
  seed = 123456789;
  a = d3_uniform ( N, &seed );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute b = A * x.
//
  b = d3_mxv ( N, a, x );
//
//  Wipe out the solution.
//  Solve the system.
//
  delete [] x;

  x = d3_np_fs ( N, a, b );
//
//  Print the solution.
//
  r8vec_print ( N, x, "  Computed solution:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests DATA_TO_DIF and DIF_VAL.
//
//  Discussion:
//
//    This test demonstrates how divided difference approximation
//    improves with N.
//
//    Evaluate these polynomials at 2.5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXTAB 8

  double diftab[MAXTAB];
  double error;
  int j;
  int ntab;
  double true_value;
  double xtab[MAXTAB];
  double xval;
  double ytab[MAXTAB];
  double yval;

  cout << "\n";
  cout << "TEST006\n";
  cout << "  Approximate Y = EXP(X) using orders 1 to " << MAXTAB << ".\n";

  cout << "\n";
  cout << "  Original data:\n";
  cout << "\n";
  cout << "       X          Y\n";
  cout << "\n";
  for ( j = 0; j < MAXTAB; j++ )
  {
    xtab[j] = ( double ) j;
    ytab[j] = exp ( xtab[j] );

    cout                        << "  "
         << setw(10) << xtab[j] << "  "
         << setw(10) << ytab[j] << "\n";
  }

  xval = 2.5;
  true_value = exp ( xval );
  cout << "\n";
  cout << "  Evaluate at X = " << xval << " where EXP(X) = "
    << true_value << "\n";
  cout << "\n";
  cout << "  Order  Approximate Y     Error\n";
  cout << "\n";

  for ( ntab = 1; ntab <= MAXTAB; ntab++ )
  {

    for ( j = 0; j < ntab; j++ )
    {
      xtab[j] = ( double ) j;
      ytab[j] = exp ( xtab[j] );
    }

    data_to_dif ( ntab, xtab, ytab, diftab );

    yval = dif_val ( ntab, xtab, diftab, xval );

    error = yval - true_value;

    cout                      << "  "
         << setw(6)  << ntab  << "  "
         << setw(10) << yval  << "  "
         << setw(10) << error << "\n";

  }

  return;
# undef MAXTAB
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests BASIS_FUNCTION_B_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 5

  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double tdata[NDATA] = { 0.0, 1.0, 4.0, 6.0, 10.0 };
  double thi;
  double tlo;
  double tval;
  double yval;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  BASIS_FUNCTION_B_VAL evaluates the \n";
  cout << "    B spline basis function.\n";
  cout << "\n";
  cout << "           T            B(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {

      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_function_b_val ( tdata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  return;
# undef NDATA
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BASIS_FUNCTION_BETA_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 5

  double beta1;
  double beta2;
  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double tdata[NDATA] = { 0.0, 1.0, 4.0, 6.0, 10.0 };
  double thi;
  double tlo;
  double tval;
  double yval;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  BASIS_FUNCTION_BETA_VAL evaluates the \n";
  cout << "    Beta spline basis function.\n";

  beta1 = 1.0;
  beta2 = 0.0;

  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";
  cout << "\n";
  cout << "              T           B(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_function_beta_val ( beta1, beta2, tdata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  beta1 = 1.0;
  beta2 = 100.0;

  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";
  cout << "\n";
  cout << "              T           B(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_function_beta_val ( beta1, beta2, tdata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }

      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  beta1 = 100.0;
  beta2 = 0.0;

  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";
  cout << "\n";
  cout << "              T           B(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_function_beta_val ( beta1, beta2, tdata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  return;
# undef NDATA
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests BASIS_MATRIX_B_UNI and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { -1.0, 0.0, 1.0, 2.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 4.0, 7.0, 12.0, 19.0 };
  double yval;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  BASIS_MATRIX_B_UNI sets up the basis matrix\n";
  cout << "    for the uniform B spline.\n";

  mbasis = basis_matrix_b_uni ( );

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  left = 2;

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";

    }
  }

  delete [] mbasis;

  return;
# undef N
# undef NDATA
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests BASIS_MATRIX_BETA_UNI and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  double beta1;
  double beta2;
  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { -1.0, 0.0, 1.0, 2.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 4.0, 7.0, 12.0, 19.0 };
  double yval;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  BASIS_MATRIX_BETA_UNI sets up the basis matrix\n";
  cout << "    for the uniform beta spline.\n";
//
//  First test
//
  beta1 = 1.0;
  beta2 = 0.0;

  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";

  mbasis = basis_matrix_beta_uni ( beta1, beta2 );

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }
//
//  Second test
//
  beta1 = 1.0;
  beta2 = 100.0;

  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";

  delete [] mbasis;
  mbasis = basis_matrix_beta_uni ( beta1, beta2 );

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }
//
//  Third test
//
  beta1 = 100.0;
  beta2 = 0.0;

  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";

  delete [] mbasis;
  mbasis = basis_matrix_beta_uni ( beta1, beta2 );

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  delete [] mbasis;

  return;
# undef N
# undef NDATA
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests BASIS_MATRIX_BEZIER and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { 0.0, 0.0, 1.0, 1.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 7.0,  8.3333333,   10.0, 12.0 };
  double yval;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  BASIS_MATRIX_BEZIER sets up the basis\n";
  cout << "    matrix for the uniform Bezier spline.\n";

  mbasis = basis_matrix_bezier ( );

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";


  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  delete [] mbasis;

  return;
# undef N
# undef NDATA
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests BASIS_MATRIX_HERMITE and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { 0.0, 0.0, 1.0, 1.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 7.0, 12.0, 4.0, 6.0 };
  double yval;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  BASIS_MATRIX_HERMITE sets up the basis matrix\n";
  cout << "    for the Hermite spline.\n";

  mbasis = basis_matrix_hermite ( );

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";


  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  delete [] mbasis;

  return;
# undef N
# undef NDATA
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests BASIS_MATRIX_OVERHAUSER_UNI and BASIS_MATRIX_TMP.
//
//  Discussion:
//
//    YDATA(1:NDATA) = ( TDATA(1:NDATA) + 2 )**2 + 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA] = { -1.0, 0.0, 1.0, 2.0 };
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA] = { 4.0, 7.0, 12.0, 19.0 };
  double yval;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  BASIS_MATRIX_OVERHAUSER_UNI sets up the basis\n";
  cout << "    matrix for the uniform Overhauser spline.\n";

  mbasis = basis_matrix_overhauser_uni ( );

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout                         << "  "
           << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  delete [] mbasis;

  return;
# undef N
# undef NDATA
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
//
//  Discussion:
//
//    YDATA(1:NDATA) = ( TDATA(1:NDATA) - 2 )**2 + 3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  double alpha;
  double beta;
  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the\n";
  cout << "    basis matrix for the nonuniform Overhauser\n";
  cout << "    spline.\n";

  tdata[0] = 0.0;
  tdata[1] = 1.0;
  tdata[2] = 2.0;
  tdata[3] = 3.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  cout << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA  = " << beta  << "\n";

  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] - 2.0 ), 2 ) + 3.0;
  }

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout                         << "  "
           << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  tdata[0] = 0.0;
  tdata[1] = 1.0;
  tdata[2] = 2.0;
  tdata[3] = 5.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  cout << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA  = " << beta  << "\n";

  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] - 2.0 ), 2 ) + 3.0;
  }

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout                         << "  "
           << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  tdata[0] = 0.0;
  tdata[1] = 3.0;
  tdata[2] = 4.0;
  tdata[3] = 5.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  cout << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA  = " << beta  << "\n";

  delete [] mbasis;
  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] - 2.0 ), 2 ) + 3.0;
  }

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  delete [] mbasis;

  return;
# undef N
# undef NDATA
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4
# define NDATA 4

  double alpha;
  double beta;
  int i;
  int j;
  int jhi;
  int left;
  char mark;
  double *mbasis;
  int nsample = 4;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the\n";
  cout << "    basis matrix for the nonuniform Overhauser \n";
  cout << "    spline.\n";
  cout << "\n";
  cout << "  First test that the nonuniform code can match\n";
  cout << "  the uniform code.  Compare these results with\n";
  cout << "  the uniform output.\n";
  cout << "\n";

  tdata[0] = -1.0;
  tdata[1] =  0.0;
  tdata[2] =  1.0;
  tdata[3] =  2.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  cout << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA  = " << beta  << "\n";

  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] + 2.0 ), 2 ) + 3.0;
  }

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  cout << "\n";
  cout << "  Now test that the nonuniform code on a\n";
  cout << "  nonuniform grid.\n";
  cout << "\n";

  tdata[0] = -4.0;
  tdata[1] = -3.0;
  tdata[2] = -1.0;
  tdata[3] =  2.0;

  alpha = ( tdata[2] - tdata[1] ) / ( tdata[2] - tdata[0] );
  beta =  ( tdata[2] - tdata[1] ) / ( tdata[3] - tdata[1] );

  cout << "\n";
  cout << "  ALPHA = " << alpha << "\n";
  cout << "  BETA  = " << beta  << "\n";

  delete [] mbasis;
  mbasis = basis_matrix_overhauser_nonuni ( alpha, beta );

  for ( i = 0; i < NDATA; i++ )
  {
    ydata[i] = pow ( ( tdata[i] + 2.0 ), 2 ) + 3.0;
  }

  cout << "\n";
  cout << "    TDATA, YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
      cout << setw(12) << tdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
  }

  left = 2;

  cout << "\n";
  cout << "              T      Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = basis_matrix_tmp ( left, N, mbasis, NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  delete [] mbasis;

  return;
# undef N
# undef NDATA
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests BC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  int i;
  int nsample = 101;
  double t;
  double xcon[N+1] = { 0.0, 0.75, 1.0 };
  double xval;
  double ycon[N+1] = { 1.0, 0.0,  1.0 };
  double yval;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  BC_VAL evaluates a general Bezier function.\n";
//
//  One point on the curve should be about (0.75, 0.536).
//
  cout << "\n";
  cout << "        T             X(T)           Y(T)\n";
  cout << "\n";

  for ( i = 1; i <= nsample; i++ )
  {
    t = ( double ) ( i - 1 ) / ( double ) ( nsample - 1 );

    bc_val ( N, t, xcon, ycon, &xval, &yval );

      cout << setw(12) << t    << "  "
           << setw(12) << xval << "  "
           << setw(12) << yval << "\n";
  }

  cout << "\n";
  cout << "  The point ( 0.75, 0.536 ) should be on the curve.\n";

  return;
# undef N
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests BEZ_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  double a = 0.0;
  double b = 1.0;
  double bval;
  int i;
  int nsample = 21;
  double x;
  double y[N+1] = { 1.0, 0.0, 1.0 };

  cout << "\n";
  cout << "TEST11\n";
  cout << "  BEZ_VAL evaluates a Bezier function.\n";
//
//  One point on the curve should be (0.75, 20/32).
//
  cout << "\n";
  cout << "        I             X           B\n";
  cout << "\n";

  for ( i = 1; i <= nsample; i++ )
  {
    x = ( ( double ) ( nsample - i     ) * a
        + ( double ) (           i - 1 ) * b )
        / ( double ) ( nsample     - 1 );

    bval = bez_val ( N, x, a, b, y );

    cout << setw(6)  << i    << "  "
         << setw(12) << x    << "  "
         << setw(12) << bval << "\n";
  }

  cout << "\n";
  cout << "  When X = " << 0.75 << "\n";
  cout << "  BEZ_VAL(X) should be " << 0.625 << "\n";

  return;
# undef N
}
//****************************************************************************80

void test115 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST115 tests BP01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  double a = 0.0;
  double b = 1.0;
  double *bern;
  int i;
  int j;
  int n;
  int nsample = 11;
  double x;

  cout << "\n";
  cout << "TEST115\n";
  cout << "  BP01 evaluates the Bernstein basis polynomials\n";
  cout << "  for the interval [0,1].\n";

  for ( n = 0; n <= N_MAX; n++ )
  {
    cout << "\n";
    cout << "  Degree N = " << n << "\n";
    cout << "\n";
    cout << "   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)\n";
    cout << "\n";

    for ( i = 1; i <= nsample; i++ )
    {
      x = ( ( double ) ( nsample - i     ) * a
          + ( double ) (           i - 1 ) * b )
          / ( double ) ( nsample     - 1 );

      bern = bp01 ( n, x );

      cout << setw(12) << x << "  ";
      for ( j = 0; j <= n; j++ )
      {
        cout << setw(12) << bern[j] << "  ";
      }
      cout << "\n";

      delete [] bern;
    }
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test116 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST116 tests BPAB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 3

  double a = 1.0;
  double b = 3.0;
  double *bern;
  int i;
  int j;
  int n;
  int nsample = 11;
  double x;

  cout << "\n";
  cout << "TEST116\n";
  cout << "  BPAB evaluates the Bernstein basis polynomials\n";
  cout << "  for the interval [A,B].\n";

  for ( n = 0; n <= N_MAX; n++ )
  {
    cout << "\n";
    cout << "  Degree N = " << n << "\n";
    cout << "\n";
    cout << "   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)\n";
    cout << "\n";

    for ( i = 1; i <= nsample; i++ )
    {
      x = ( ( double ) ( nsample - i     ) * a
          + ( double ) (           i - 1 ) * b )
          / ( double ) ( nsample     - 1 );

      bern = bpab ( n, a, b, x );

      cout << setw(12) << x << "  ";
      for ( j = 0; j <= n; j++ )
      {
        cout << setw(12) << bern[j] << "  ";
      }
      cout << "\n";

      delete [] bern;
    }
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests BPAB_APPROX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXDATA 10

  double a;
  double b;
  int i;
  int ndata;
  int nsample;
  double xdata[MAXDATA+1];
  double xval;
  double ydata[MAXDATA+1];
  double yval;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  BPAB_APPROX evaluates the Bernstein polynomial\n";
  cout << "  approximant to a function F(X).\n";

  a = 1.0;
  b = 3.0;

  for ( ndata = 0; ndata <= 9; ndata = ndata + 3 )
  {
    for ( i = 0; i <= ndata; i++ )
    {
      if ( ndata == 0 )
      {
        xdata[i] = 0.5 * ( a + b );
      }
      else
      {
        xdata[i] = ( ( double ) ( ndata - i ) * a
                   + ( double ) (         i ) * b )
                   / ( double ) ( ndata );
      }

      ydata[i] = sin ( xdata[i] );

    }

    cout << "\n";
    cout << "    XDATA    YDATA\n";
    cout << "\n";
    for ( i = 0; i <= ndata; i++ )
    {
      cout << setw(12) << xdata[i] << "  "
           << setw(12) << ydata[i] << "\n";
    }

    cout << "\n";
    cout << "  Bernstein approximant of degree N = " << ndata << "\n";
    cout << "\n";
    cout << "    X      F(X)     BERN(X)    ERROR\n";
    cout << "\n";

    nsample = 2 * ndata + 1;

    for ( i = 1; i <= nsample; i++ )
    {
      if ( nsample == 1 )
      {
        xval = 0.5 * ( a + b );
      }
      else
      {
        xval = ( ( double ) ( nsample - i     ) * a
               + ( double ) (           i - 1 ) * b )
               / ( double ) ( nsample     - 1 );
      }

      yval = bpab_approx ( ndata, a, b, ydata, xval );

      cout << setw(12) <<              xval   << "  "
           << setw(12) <<        sin ( xval ) << "  "
           << setw(12) << yval                << "  "
           << setw(12) << yval - sin ( xval ) << "\n";
    }
  }

  return;
# undef MAXDATA
}
//****************************************************************************80

void test125 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST125 tests LEAST_SET_OLD and LEAST_VAL_OLD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXDEG 6
# define NTAB 21

  double b[MAXDEG];
  double c[MAXDEG+1];
  double d[MAXDEG-1];
  double eps;
  double error;
  int i;
  int ierror;
  int j;
  int jhi;
  int ndeg;
  double ptab[NTAB];
  double xtab[NTAB];
  double xval;
  double ytab[NTAB];
  double ytrue;
  double yval;

  cout << "\n";
  cout << "TEST125\n";
  cout << "  LEAST_SET_OLD sets a least squares polynomial,\n";
  cout << "  LEAST_VAL_OLD evaluates it.\n";

  for ( i = 0; i < NTAB; i++ )
  {
    xtab[i] = ( ( double ) ( NTAB - i - 1 ) * ( -1.0 )
              + ( double ) (        i     ) * ( +1.0 ) )
              / ( double ) ( NTAB     - 1 );
    ytab[i] = ( double ) ( ( int ) ( exp ( xtab[i] ) * 100.0 + 0.5 ) )
      / 100.0;
  }

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Number of data values = " << NTAB << "\n";
  cout << "\n";
  cout << "       X             Y\n";
  cout << "\n";

    cout << "\n";
    cout << "    XTAB    YTAB\n";
    cout << "\n";
    for ( i = 0; i < NTAB; i++ )
    {
      cout << setw(12) << xtab[i] << "  "
           << setw(12) << ytab[i] << "\n";
    }

  for ( ndeg = 1; ndeg <= MAXDEG; ndeg++ )
  {
    cout << "\n";
    cout << "  Use a polynomial of degree: " << ndeg << "\n";
    cout << "\n";

    least_set_old ( NTAB, xtab, ytab, ndeg, ptab, b, c, d, &eps, &ierror );

    cout << "\n";
    cout << "  Total approximation error = " << eps << "\n";
    cout << "\n";
    cout << "       X            F(X)          P(X)          Error\n";
    cout << "\n";

    for ( i = 1; i <= NTAB; i++ )
    {
      if ( i < NTAB )
      {
        jhi = 2;
      }
      else
      {
        jhi = 0;
      }

      for ( j = 0; j <= jhi; j++ )
      {
        if ( i < NTAB )
        {
          xval = ( ( double ) ( 3 - j ) * xtab[i-1]
                 + ( double ) (     j ) * xtab[i]   )
                 / ( double ) ( 3     );
        }
        else
        {
          xval = xtab[i-1];
        }

        yval = least_val_old ( xval, ndeg, b, c, d );

        ytrue = ( double ) ( ( int ) ( exp ( xval ) * 100.0 + 0.5 ) )
          / 100.0;

        error = yval - ytrue;

        cout << setw(12) << xval  << "  "
             << setw(12) << yval  << "  "
             << setw(12) << ytrue << "  "
             << setw(12) << error << "\n";
      }
    }
  }

  return;
# undef MAXDEG
}
//****************************************************************************80

void test126 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST126 tests LEAST_SET and LEAST_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define POINT_NUM 21
# define NTERMS 4

  double b[NTERMS];
  double c[NTERMS];
  double d[NTERMS];
  double f[POINT_NUM];
  double fp[POINT_NUM];
  int i;
  int nterms2;
  double px;
  double w[POINT_NUM];
  double x[POINT_NUM];

  for ( i = 0; i < POINT_NUM; i++ )
  {
    w[i] = 1.0;
    x[i] = -1.0 + ( double ) ( i ) / 10.0;
    f[i] = x[i] * x[i] - x[i] - 6.0;
    fp[i] = 2.0 * x[i] - 1.0;
  }

  least_set ( POINT_NUM, x, f, w, NTERMS, b, c, d );

  cout << "\n";
  cout << "TEST126\n";
  cout << "  LEAST_SET sets a least squares polynomial,\n";
  cout << "  LEAST_VAL evaluates it.\n";
  cout << "\n";
  cout << "  X, F(X), P(X), Error\n";
  cout << "\n";
  for ( nterms2 = 1; nterms2 <= NTERMS; nterms2++ )
  {
    cout << "\n";
    cout << "  Using polynomial order = " << nterms2 << "\n";
    cout << "\n";
    for ( i = 0; i < POINT_NUM; i++ )
    {
      px = least_val ( nterms2, b, c, d, x[i] );
      cout << "  " << setw(12) << x[i]
           << "  " << setw(12) << f[i]
           << "  " << setw(12) << px
           << "  " << setw(12) << px - f[i] << "\n";
    }
  }

  return;
# undef POINT_NUM
# undef NTERMS
}
//****************************************************************************80

void test127 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST127 tests LEAST_SET and LEAST_VAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define POINT_NUM 21
# define NTERMS 4

  double b[NTERMS];
  double c[NTERMS];
  double d[NTERMS];
  double f[POINT_NUM];
  double fp[POINT_NUM];
  int i;
  int nterms2;
  double px;
  double pxp;
  double w[POINT_NUM];
  double x[POINT_NUM];

  for ( i = 0; i < POINT_NUM; i++ )
  {
    w[i] = 1.0;
    x[i] = -1.0 + ( double ) ( i ) / 10.0;
    f[i] = x[i] * x[i] - x[i] - 6.0;
    fp[i] = 2.0 * x[i] - 1.0;
  }

  least_set ( POINT_NUM, x, f, w, NTERMS, b, c, d );

  cout << "\n";
  cout << "TEST127\n";
  cout << "  LEAST_SET sets a least squares polynomial,\n";
  cout << "  LEAST_VAL2 evaluates it.\n";
  cout << "\n";
  cout << "  X, F(X), P(X), FP(X), PP(X)\n";
  cout << "\n";

  for ( nterms2 = 1; nterms2 <= NTERMS; nterms2++ )
  {
    cout << "\n";
    cout << "  Using polynomial order = " << nterms2 << "\n";
    cout << "\n";
    for ( i = 0; i < POINT_NUM; i++ )
    {
      least_val2 ( nterms2, b, c, d, x[i], &px, &pxp );
      cout << "  " << setw(12) << x[i]
           << "  " << setw(12) << f[i]
           << "  " << setw(12) << px
           << "  " << setw(12) << fp[i]
           << "  " << setw(12) << pxp << "\n";
    }
  }

  return;
# undef POINT_NUM
# undef NTERMS
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests SPLINE_B_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 11

  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  SPLINE_B_VAL evaluates the\n";
  cout << "    B spline.\n";
  cout << "\n";
  cout << "  TDATA   YDATA\n";
  cout << "\n";

  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = ( double ) i;
    ydata[i] = sin ( 2.0 * pi * tdata[i] / ( double ) ( NDATA - 1) );
    cout << setw(12) << tdata[i] << "  "
         << setw(12) << ydata[i] << "\n";
  }

  cout << "\n" ;
  cout << "    T, Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[1] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_b_val ( NDATA, tdata, ydata, tval );

      if ( 0 < i && j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }

      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";

    }

  }

  return;
# undef NDATA
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests SPLINE_BETA_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 11

  double beta1;
  double beta2;
  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  SPLINE_BETA_VAL evaluates the BETA spline.\n";
  cout << "\n";
  cout << "       TDATA         YDATA\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = ( double ) i;
    ydata[i] = sin ( 2.0 * pi * tdata[i] / ( double ) ( NDATA - 1) );
    cout << setw(12) << tdata[i] << "  "
         << setw(12) << ydata[i] << "\n";
  }

  beta1 = 1.0;
  beta2 = 0.0;

  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";
  cout << "\n";
  cout << "    T, Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {

    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_beta_val ( beta1, beta2, NDATA, tdata, ydata, tval );

      if ( 0 < i && j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  beta1 = 1.0;
  beta2 = 100.0;

  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";
  cout << "\n";
  cout << "    T, Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_beta_val ( beta1, beta2, NDATA, tdata, ydata, tval );

      if ( 0 < i && j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }

  }

  beta1 = 100.0;
  beta2 = 0.0;

  cout << "\n";
  cout << "  BETA1 = " << beta1 << "\n";
  cout << "  BETA2 = " << beta2 << "\n";
  cout << "\n";
  cout << "    T, Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_beta_val ( beta1, beta2, NDATA, tdata, ydata, tval );

      if ( 0 < i && j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";
    }
  }

  return;
# undef NDATA
}
//****************************************************************************80

void test145 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST145 tests SPLINE_CONSTANT_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 12
# define N_TEST 20

  double fval;
  int i;
  int j;
  int seed = 123456789;
  double *tdata;
  double thi;
  double tlo;
  double *t_test;
  double tval;
  double ydata[NDATA];
  double yval;

  cout << "\n";
  cout << "TEST145\n";
  cout << "  SPLINE_CONSTANT_VAL evaluates a piecewise\n";
  cout << "  constant spline.\n";
  cout << "\n";
  cout << "  Runge's function, evenly spaced knots.\n";
//
//  Set the data.
//
  tlo = -1.0;
  thi = +1.0;

  tdata = r8vec_even_new ( NDATA-1, tlo, thi );

  for ( i = 0; i < NDATA; i++ )
  {
    if ( i == 0 )
    {
      tval = tdata[0];
    }
    else if ( i < NDATA - 1 )
    {
      tval = 0.5 * ( tdata[i-1] + tdata[i] );
    }
    else if ( i == NDATA - 1 )
    {
      tval = tdata[i-1];
    }

    ydata[i] = frunge ( tval );

  }

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Number of data values = " << NDATA << "\n";
  cout << "\n";
  cout << "       T             Y\n";
  cout << "\n";

  for ( i = 0; i < NDATA; i++ )
  {
    cout << "  *"
         << "              "
         << setw(12) << ydata[i] << "\n";
    if ( i < NDATA -1 )
    {
      cout << "  *"
            << setw(12) << tdata[i] << "\n";
    }
  }
//
//  Sample the spline.
//
  t_test = r8vec_uniform_new ( N_TEST, tlo-1.0, thi+1.0, &seed );

  r8vec_sort_bubble_a ( N_TEST, t_test );

  cout << "\n";
  cout << "     T     Y(interp)    Y(exact)\n";
  cout << "\n";

  j = 0;
  cout << "  *"
       << "              "
       << setw(12) << ydata[j] << "\n";
  j = j + 1;

  for ( i = 0; i < N_TEST; i++ )
  {
    tval = t_test[i];

    yval = spline_constant_val ( NDATA, tdata, ydata, tval );

    if ( j <= NDATA - 1 )
    {
      while ( tdata[j-1] <= tval )
      {
        fval = frunge ( tdata[j-1] );
        cout << "  *"
             << setw(12) << tdata[j-1] << "  "
             << "            "         << "  "
             << setw(12) << fval       << "\n";
        cout << "  *"
             << "            "         << "  "
             << setw(12) << ydata[j]   << "\n";
        j = j + 1;
        if ( NDATA <= j )
        {
          break;
        }
      }
    }

    fval = frunge ( tval );

    cout << "   "
         << setw(12) << tval << "  "
         << setw(12) << yval << "  "
         << setw(12) << fval << "\n";

  }

  delete [] tdata;
  delete [] t_test;

  return;
# undef NDATA
# undef N_TEST
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2013
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  cout << "\n";
  cout << "TEST15\n";
  cout << "  SPLINE_CUBIC_SET sets up a cubic spline;\n";
  cout << "  SPLINE_CUBIC_VAL evaluates it.\n";
  cout << "\n";
  cout << "  Runge's function, evenly spaced knots.\n";
  cout << "\n";
  cout << "     I     T         Y\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( ( double ) ( N - i - 1 ) * (-1.0)
            + ( double ) (     i     ) * (+1.0) )
            / ( double ) ( N     - 1 );
    y[i] =  frunge ( t[i] );
    cout << setw(6)  << i    << "  "
         << setw(10) << t[i] << "  "
         << setw(10) << y[i] << "\n";
  }
//
//  Try various boundary conditions.
//
  for ( k = 0; k <= 4; k++ )
  {
    if ( k == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      ibcend = 0;
      ybcend = 0.0;

      cout << "\n";
      cout << "  Boundary condition 0 at both ends:\n";
      cout << "  Spline is quadratic in boundary intervals.\n";
    }
    else if ( k == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fprunge ( t[0] );

      ibcend = 1;
      ybcend = fprunge ( t[N-1] );

      cout << "\n";
      cout << "  Boundary condition 1 at both ends:\n";
      cout << "  Y'(left) =  " << ybcbeg << "\n";
      cout << "  Y'(right) = " << ybcend << "\n";

    }
    else if ( k == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fpprunge ( t[0] );

      ibcend = 2;
      ybcend = fpprunge ( t[N-1] );

      cout << "\n";
      cout << "  Boundary condition 2 at both ends:\n";
      cout << "  YP''(left) =  " << ybcbeg << "\n";
      cout << "  YP''(right) = " << ybcend << "\n";
    }
    else if ( k == 3 )
    {
      ibcbeg = 2;
      ybcbeg = 0.0;

      ibcend = 2;
      ybcend = 0.0;

      cout << "\n";
      cout << "  Natural spline:\n";
      cout << "  YP''(left) =  " << ybcbeg << "\n";
      cout << "  YP''(right) = " << ybcend << "\n";
    }
    else if ( k == 4 )
    {
      ibcbeg = 3;
      ibcend = 3;

      cout << "\n";
      cout << "  \"Not-a-knot\" spline:\n";
    }

    ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

    cout << "\n";
    cout << "  SPLINE''(T), F''(T):\n";
    cout << "\n";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(10) << ypp[i] << "  "
           << setw(10) << fpprunge ( t[i] ) << "\n";
    }

    cout << "\n";
    cout << "  T, SPLINE(T), F(T)\n";
    cout << "\n";

    for ( i = 0; i <= N; i++ )
    {
      if ( i == 0 )
      {
        jhi = 1;
      }
      else if ( i < N )
      {
        jhi = 2;
      }
      else
      {
        jhi = 2;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        if ( i == 0 )
        {
          tval = t[0] - 1.0;
        }
        else if ( i < N )
        {
          tval = (
              ( double ) ( jhi - j + 1 ) * t[i-1]
            + ( double ) (       j - 1 ) * t[i] )
            / ( double ) ( jhi         );
        }
        else
        {
          if ( j == 1 )
          {
            tval = t[N-1];
          }
          else
          {
            tval = t[N-1] + 1.0;
          }
        }

        yval = spline_cubic_val ( N, t, y, ypp, tval, &ypval, &yppval );

        cout << setw(10) << tval << "  "
             << setw(10) << yval << "  "
             << setw(10) << frunge ( tval ) << "\n";
      }
    }
    delete [] ypp;
  }

  return;
# undef N
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k;
  int left;
  int left_in;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  cout << "\n";
  cout << "TEST16\n";
  cout << "  SPLINE_CUBIC_SET sets up a cubic spline;\n";
  cout << "  SPLINE_CUBIC_VAL2 evaluates it.\n";
  cout << "\n";
  cout << "  Runge's function, evenly spaced knots.\n";
  cout << "\n";
  cout << "     I      T       Y\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( ( double ) ( N - i - 1 ) * (-1.0)
            + ( double ) (     i     ) * (+1.0) )
            / ( double ) ( N     - 1 );
    y[i] =  frunge ( t[i] );
    cout << setw(6)  << i    << "  "
         << setw(12) << t[i] << "  "
         << setw(12) << y[i] << "\n";
  }
//
//  Try all three types of boundary condition.
//
  for ( k = 0; k < 3; k++ )
  {
    if ( k == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      ibcend = 0;
      ybcend = 0.0;

      cout << "\n";
      cout << "  Boundary condition 0 at both ends:\n";
      cout << "  Spline is quadratic in boundary intervals.\n";
    }
    else if ( k == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fprunge ( t[0] );

      ibcend = 1;
      ybcend = fprunge ( t[N-1] );

      cout << "\n";
      cout << "  Boundary condition 1 at both ends:\n";
      cout << "  Y'(left) =  " << ybcbeg << "\n";
      cout << "  Y'(right) = " << ybcend << "\n";
    }
    else if ( k == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fpprunge ( t[0] );

      ibcend = 2;
      ybcend = fpprunge ( t[N-1] );

      cout << "\n";
      cout << "  Boundary condition 2 at both ends:\n";
      cout << "  YP''(left) =  " << ybcbeg << "\n";
      cout << "  YP''(right) = " << ybcend << "\n";
    }

    ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

    cout << "\n";
    cout << "  SPLINE''(T), F''(T):\n";
    cout << "\n";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(12) << ypp[i] << "  "
           << setw(12) << fpprunge(t[i]) << "\n";
    }

    left = 0;

    cout << "\n";
    cout << "  T, SPLINE(T), F(T), LEFT_IN, LEFT_OUT\n";
    cout << "\n";

    for ( i = 0; i <= N; i++ )
    {
      if ( i == 0 )
      {
        jhi = 1;
      }
      else if ( i < N )
      {
        jhi = 2;
      }
      else
      {
        jhi = 2;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        if ( i == 0 )
        {
          tval = t[0] - 1.0;
        }
        else if ( i < N )
        {
          tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
                 + ( double ) (       j - 1 ) * t[i] )
                 / ( double ) ( jhi         );
        }
        else
        {
          if ( j == 1 )
          {
            tval = t[N-1];
          }
          else
          {
            tval = t[N-1] + 1.0;
          }
        }

        left_in = left;
        spline_cubic_val2 ( N, t, tval, &left, y, ypp, &yval, &ypval,
          &yppval );

        cout << setw(12) << tval << "  "
             << setw(12) << yval << "  "
             << setw(12) << frunge ( tval ) << "  "
             << setw(6)  << left_in << "  "
             << setw(6)  << left << "\n";
      }
    }
    delete [] ypp;
  }

  return;

# undef N
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  cout << "\n";
  cout << "TEST17\n";
  cout << "  SPLINE_CUBIC_SET sets up a cubic spline;\n";
  cout << "  SPLINE_CUBIC_VAL evaluates it.\n";
  cout << "\n";
  cout << "  Cubic function, unevenly spaced knots.\n";
  cout << "\n";
  cout << "  T, Y\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( double ) ( i ) / ( double ) ( N - 1 );
    t[i] = t[i] * t[i];
    y[i] =  fcube ( t[i] );
    cout << setw(6)  << i    << "  "
         << setw(10) << t[i] << "  "
         << setw(10) << y[i] << "\n";
  }
//
//  Try all three types of boundary condition.
//
  for ( k = 0; k < 3; k++ )
  {
    if ( k == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      ibcend = 0;
      ybcend = 0.0;

      cout << "\n";
      cout << "  Boundary condition 0 at both ends:\n";
      cout << "  Spline is quadratic in boundary intervals.\n";

    }
    else if ( k == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fpcube ( t[0] );

      ibcend = 1;
      ybcend = fpcube ( t[N-1] );

      cout << "\n";
      cout << "  Boundary condition 1 at both ends:\n";
      cout << "  Y'(left) =  " << ybcbeg << "\n";
      cout << "  Y'(right) = " << ybcend << "\n";

    }
    else if ( k == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fppcube ( t[0] );

      ibcend = 2;
      ybcend = fppcube ( t[N-1] );

      cout << "\n";
      cout << "  Boundary condition 2 at both ends:\n";
      cout << "  YP''(left) =  " << ybcbeg << "\n";
      cout << "  YP''(right) = " << ybcend << "\n";

    }

    ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

    cout << "\n";
    cout << "  SPLINE''(T), F''(T):\n";
    cout << "\n";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(12) << ypp[i]        << "  "
           << setw(12) << fppcube(t[i]) << "\n";
    }

    cout << "\n";
    cout << "  T, SPLINE(T), F(T)\n";
    cout << "\n";

    for ( i = 0; i <= N; i++ )
    {
      if ( i == 0 )
      {
        jhi = 1;
      }
      else if ( i < N )
      {
        jhi = 2;
      }
      else
      {
        jhi = 2;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        if ( i == 0 )
        {
          tval = t[0] - 1.0;
        }
        else if ( i < N )
        {
          tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
                 + ( double ) (       j - 1 ) * t[i] )
                 / ( double ) ( jhi         );
        }
        else
        {
          if ( j == 1 )
          {
            tval = t[N-1];
          }
          else
          {
            tval = t[N-1] + 1.0;
          }
        }

        yval = spline_cubic_val ( N, t, y, ypp, tval, &ypval, &yppval );

        cout << setw(12) << tval << "  "
             << setw(12) << yval << "  "
             << setw(12) << fcube ( tval ) << "\n";
      }
    }
    delete [] ypp;
  }

  return;
# undef N
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  cout << "\n";
  cout << "TEST18\n";
  cout << "  SPLINE_CUBIC_SET sets up a cubic spline;\n";
  cout << "  SPLINE_CUBIC_VAL evaluates it.\n";
  cout << "\n";
  cout << "  Cubic function, evenly spaced knots.\n";
  cout << "\n";
  cout << "        T            Y\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( double ) ( i ) / ( double ) ( N - 1 );
    y[i] =  fcube ( t[i] );
    cout << setw(12) << t[i] << "  "
         << setw(12) << y[i] << "\n";
  }
//
//  Try all three types of boundary condition.
//
  for ( k = 0; k < 3; k++ )
  {
    if ( k == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      ibcend = 0;
      ybcend = 0.0;

      cout << "\n";
      cout << "  Boundary condition 0 at both ends:\n";
      cout << "  Spline is quadratic in boundary intervals.\n";

    }
    else if ( k == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fpcube ( t[0] );

      ibcend = 1;
      ybcend = fpcube ( t[N-1] );

      cout << "\n";
      cout << "  Boundary condition 1 at both ends:\n";
      cout << "  Y'(left) =  " << ybcbeg << "\n";
      cout << "  Y'(right) = " << ybcend << "\n";

    }
    else if ( k == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fppcube ( t[0] );

      ibcend = 2;
      ybcend = fppcube ( t[N-1] );

      cout << "\n";
      cout << "  Boundary condition 2 at both ends:\n";
      cout << "  YP''(left) =  " << ybcbeg << "\n";
      cout << "  YP''(right) = " << ybcend << "\n";

    }

    ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

    cout << "\n";
    cout << "     SPLINE''(T)       F''(T):\n";
    cout << "\n";
    for ( i = 0; i < N; i++ )
    {
      cout << setw(12) << ypp[i]        << "  "
           << setw(12) << fppcube(t[i]) << "\n";
    }

    cout << "\n";
    cout << "        T   SPLINE(T)     F(T)\n";
    cout << "\n";

    for ( i = 0; i <= N; i++ )
    {
      if ( i == 0 )
      {
        jhi = 1;
      }
      else if ( i < N )
      {
        jhi = 2;
      }
      else
      {
        jhi = 2;
      }

      for ( j = 1; j <= jhi; j++ )
      {
        if ( i == 0 )
        {
          tval = t[0] - 1.0;
        }
        else if ( i < N )
        {
          tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
                 + ( double ) (       j - 1 ) * t[i] )
                 / ( double ) ( jhi         );
        }
        else
        {
          if ( j == 1 )
          {
            tval = t[N-1];
          }
          else
          {
            tval = t[N-1] + 1.0;
          }
        }

        yval = spline_cubic_val ( N, t, y, ypp, tval, &ypval, &yppval );

        cout << "\n";
        cout << setw(12) << tval
             << setw(12) << yval << "  "
             << setw(12) << fcube ( tval ) << "\n";
        cout << "            "
             << setw(12) << ypval << "  "
             << setw(12) << fpcube ( tval ) << "\n";
        cout << "            "
             << setw(12) << yppval << "  "
             << setw(12) << fppcube ( tval ) << "\n";
      }
    }
    delete [] ypp;
  }

  return;
# undef N
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  int i;
  int ibcbeg;
  int ibcend;
  int j;
  int jhi;
  int k1;
  int k2;
  double t[N];
  double tval;
  double y[N];
  double ybcbeg;
  double ybcend;
  double *ypp;
  double yppval;
  double ypval;
  double yval;
//
//  Set up the data.
//
  cout << "\n";
  cout << "TEST19\n";
  cout << "  SPLINE_CUBIC_SET sets up a cubic spline;\n";
  cout << "  SPLINE_CUBIC_VAL evaluates it.\n";
  cout << "\n";
  cout << "  Cubic function, evenly spaced knots.\n";
  cout << "  ONLY TWO KNOTS!\n";
  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Number of data values = " << N << "\n";
  cout << "\n";
  cout << "           T             Y\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( double ) ( i ) / ( double ) ( N - 1 );
    y[i] =  fcube ( t[i] );
    cout << setw(12) << t[i] << "  "
         << setw(12) << y[i] << "\n";
  }
//
//  Try all nine pairs of boundary condition.
//
  for ( k1 = 0; k1 < 3; k1++ )
  {
    if ( k1 == 0 )
    {
      ibcbeg = 0;
      ybcbeg = 0.0;

      cout << "\n";
      cout << "  Boundary condition 0 at left end.\n";
    }
    else if ( k1 == 1 )
    {
      ibcbeg = 1;
      ybcbeg = fpcube ( t[0] );

      cout << "\n";
      cout << "  Boundary condition 1 at left end.\n";
      cout << "  Y'(left) =  " << ybcbeg << "\n";
    }
    else if ( k1 == 2 )
    {
      ibcbeg = 2;
      ybcbeg = fppcube ( t[0] );

      cout << "\n";
      cout << "  Boundary condition 2 at left end.\n";
      cout << "  YP''(left) =  " << ybcbeg << "\n";
    }

    for ( k2 = 0; k2 < 3; k2++ )
    {
      if ( k2 == 0 )
      {
        ibcend = 0;
        ybcend = 0.0;

        cout << "  Boundary condition 0 at right end.\n";
      }
      else if ( k2 == 1 )
      {
        ibcend = 1;
        ybcend = fpcube ( t[N-1] );

        cout << "\n";
        cout << "  Boundary condition 1 at right end.\n";
        cout << "  Y'(right) = " << ybcend << "\n";
      }
      else if ( k2 == 2 )
      {
        ibcend = 2;
        ybcend = fppcube ( t[N-1] );

        cout << "\n";
        cout << "  Boundary condition 2 at right end.\n";
        cout << "  YP''(right) = " << ybcend << "\n";

      }

      ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );

      cout << "\n";
      cout << "  SPLINE''(T)        F''(T)\n";
      cout << "\n";
      for ( i = 0; i < N; i++ )
      {
        cout << setw(12) << ypp[i]        << "  "
             << setw(12) << fppcube(t[i]) << "\n";
      }

      cout << "\n";
      cout << "           T    SPLINE(T)         F(T)\n";
      cout << "\n";

      for ( i = 0; i <= N; i++ )
      {

        if ( i == 0 )
        {
          jhi = 1;
        }
        else if ( i < N )
        {
          jhi = 2;
        }
        else
        {
          jhi = 2;
        }

        for ( j = 1; j <= jhi; j++ )
        {

          if ( i == 0 )
          {
            tval = t[0] - 1.0;
          }
          else if ( i < N )
          {
            tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
                   + ( double ) (       j - 1 ) * t[i] )
                   / ( double ) ( jhi         );
          }
          else
          {
            if ( j == 1 )
            {
              tval = t[N-1];
            }
            else
            {
              tval = t[N-1] + 1.0;
            }
          }

          yval = spline_cubic_val ( N, t, y, ypp, tval, &ypval, &yppval );

          cout << "\n";
          cout << setw(12) << tval
               << setw(12) << yval << "  "
               << setw(12) << fcube ( tval ) << "\n";
          cout << "            "
               << setw(12) << ypval << "  "
               << setw(12) << fpcube ( tval ) << "\n";
          cout << "            "
               << setw(12) << yppval << "  "
               << setw(12) << fppcube ( tval ) << "\n";
        }
      }
    }
    delete [] ypp;
  }

  return;
# undef N
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests SPLINE_HERMITE_SET and SPLINE_HERMITE_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 4

  double *c;
  double fpval;
  double fval;
  int i;
  int j;
  int jhi;
  char mark;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double tval;
  double ydata[NDATA];
  double ypdata[NDATA];
  double ypval;
  double yval;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  SPLINE_HERMITE_SET sets up a Hermite spline;\n";
  cout << "  SPLINE_HERMITE_VAL evaluates it.\n";
//
//  Set the data.
//
  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = 0.5 * ( double ) i * pi / ( double ) ( NDATA - 1 );
    ydata[i] = sin ( tdata[i] );
    ypdata[i] = cos ( tdata[i] );
  }

  cout << "\n";
  cout << "  Data\n";
  cout << "\n";
  cout << "     TDATA(I)     YDATA[I]     Y'DATA[I]\n";
  cout << "\n";

  for ( i = 0; i < NDATA; i++ )
  {
    cout << setw(12) << tdata[i]  << "  "
         << setw(12) << ydata[i]  << "  "
         << setw(12) << ypdata[i] << "\n";
  }
//
//  Set up the spline.
//
  c = spline_hermite_set ( NDATA, tdata, ydata, ypdata );
//
//  Now evaluate the spline all over the place.
//
  cout << "\n";
  cout << "              T     Y(hermite)     Y(exact)   Y'(hermite)     Y'(exact)\n";
  cout << "\n";

  for ( i = 0; i < NDATA; i++ )
  {
    if ( i == NDATA-1 )
    {
      jhi = 0;
    }
    else
    {
      jhi = 2;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = 0.5 * ( double ) ( 3 * i + j ) * pi
        / ( double ) ( 3 * ( NDATA - 1 ) );

      fval = sin ( tval );
      fpval = cos ( tval );

      spline_hermite_val ( NDATA, tdata, c, tval, &yval, &ypval );

      if ( j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                      << "  "
                       << mark  << "  "
           << setw(12) << tval  << "  "
           << setw(12) << yval  << "  "
           << setw(12) << fval  << "  "
           << setw(12) << ypval << "  "
           << setw(12) << fpval << "\n";
    }

  }

  delete [] c;

  return;
# undef NDATA
}
//****************************************************************************80

void test205 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST205 tests SPLINE_LINEAR_INT and SPLINE_LINEAR_INTSET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a;
  double b;
  double data_x[N];
  double data_y[N];
  int i;
  double int_x[N+1] = { 0.0, 1.0, 4.0, 5.0, 10.0 };
  double int_v[N] = { 10.0, 2.0, 8.0, 27.5 };
  double value;

  cout << "\n";
  cout << "TEST205\n";
  cout << "  SPLINE_LINEAR_INTSET is given some interval endpoints,\n";
  cout << "  and a value associated with each interval.\n";
  cout << "\n";
  cout << "  It determines a linear spline, with breakpoints\n";
  cout << "  at the centers of each interval, whose integral\n";
  cout << "  over each interval is equal to the given value.\n";

  r8vec_print ( N+1, int_x, "  The interval end points:" );
  r8vec_print ( N, int_v, "  The desired interval integral values:" );

  spline_linear_intset ( N, int_x, int_v, data_x, data_y );

  r8vec_print ( N, data_x, "  The spline break points:" );
  r8vec_print ( N, data_y, "  The spline data values: " );

  cout << "\n";
  cout << "  As a check, call SPLINE_LINEAR_INT to compute\n";
  cout << "  the integral of the spline over each interval,\n";
  cout << "  and compare to the desired value.\n";
  cout << "\n";
  cout << "       A         B       Desired      Computed\n";
  cout << "\n";

  for ( i = 1; i <= N; i++ )
  {
    a = int_x[i-1];
    b = int_x[i];
    value = spline_linear_int ( N, data_x, data_y, a, b );
    cout << setw(8) << a << "  "
         << setw(8) << b << "  "
         << setw(12) << int_v[i-1] << "  "
         << setw(12) << value << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests SPLINE_LINEAR_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
#define N 11

  double fval;
  int i;
  int j;
  int jhi;
  double t[N];
  double tval;
  double y[N];
  double ypval;
  double yval;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  SPLINE_LINEAR_VAL evaluates a linear spline.\n";
  cout << "\n";
  cout << "  Runge's function, evenly spaced knots.\n";

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( ( double ) ( N - i - 1 ) * (-1.0)
            + ( double ) (     i     ) * (+1.0) )
            / ( double ) ( N     - 1 );
    y[i] =  frunge ( t[i] );
  }

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Number of data values = " << N << "\n";
  cout << "\n";
  cout << "        T             Y\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(12) << t[i] << "  "
         << setw(12) << y[i] << "\n";
  }

  cout << "\n";
  cout << "  Interpolation:\n";
  cout << "\n";
  cout << "       T             Y            Yexact\n";
  cout << "\n";
  for ( i = 0; i <= N; i++ )
  {

    if ( i == 0 )
    {
      jhi = 1;
    }
    else if ( i < N )
    {
      jhi = 2;
    }
    else
    {
      jhi = 2;
    }

    for ( j = 1; j <= jhi; j++ )
    {
      if ( i == 0 )
      {
        tval = t[0] - 1.0;
      }
      else if ( i < N )
      {
        tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
               + ( double ) (       j - 1 ) * t[i] )
               / ( double ) ( jhi         );
      }
      else
      {
        if ( j == 1 )
        {
          tval = t[N-1];
        }
        else
        {
          tval = t[N-1] + 1.0;
        }
      }

      spline_linear_val ( N, t, y, tval, &yval, &ypval );

      fval = frunge ( tval );

      cout << setw(12) << tval << "  "
           << setw(12) << yval << "  "
           << setw(12) << fval << "\n";

    }

  }

  return;
#undef N
}
//****************************************************************************80

void test215 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST215 tests SPLINE_LINEAR_INT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a;
  double b;
  int i;
  double int_val;
  double t[N] = { 2.0, 4.5,  7.5 };
  double y[N] = { 3.0, 3.75, 5.5 };

  cout << "\n";
  cout << "TEST215\n";
  cout << "  SPLINE_LINEAR_INT computes the integral of a linear spline.\n";
  cout << "\n";

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Number of data values = " << N << "\n";
  cout << "\n";
  cout << "	  T		Y\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << setw(12) << t[i] << "  "
         << setw(12) << y[i] << "\n";
  }

  cout << "\n";
  cout << "    A             B           Integral\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    if ( i == 1 )
    {
      a = 0.0;
      b = 4.0;
    }
    else if ( i == 2 )
    {
      a = 4.0;
      b = 5.0;
    }
    else if ( i == 3 )
    {
      a = 5.0;
      b = 10.0;
    }
    else if ( i == 4 )
    {
      a = 0.0;
      b = 10.0;
    }
    else
    {
      a = 10.0;
      b = 0.0;
    }

    int_val = spline_linear_int ( N, t, y, a, b );

    cout                        << "  "
         << setw(12) << a       << "  "
         << setw(12) << b       << "  "
         << setw(12) << int_val << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests SPLINE_OVERHAUSER_UNI_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 11

  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  SPLINE_OVERHAUSER_UNI_VAL evaluates the\n";
  cout << "    uniform Overhauser spline.\n";

  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = ( double ) ( i );
    ydata[i] = sin ( 2.0 * pi * tdata[i] / ( double ) ( NDATA - 1) );
  }

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Number of data values = " << NDATA << "\n";
  cout << "\n";
  cout << "       T             Y\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
    cout                         << "  "
         << setw(12) << tdata[i] << "  "
         << setw(12) << ydata[i] << "\n";
  }

  cout << "\n";
  cout << "    T, Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_overhauser_uni_val ( NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";

    }

  }

  return;
# undef NDATA
}
//****************************************************************************80

void test225 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST225 tests SPLINE_OVERHAUSER_NONUNI_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 11

  int i;
  int j;
  int jhi;
  char mark;
  int nsample = 4;
  double pi = 3.141592653589793;
  double tdata[NDATA];
  double thi;
  double tlo;
  double tval;
  double ydata[NDATA];
  double yval;

  cout << "\n";
  cout << "TEST225\n";
  cout << "  SPLINE_OVERHAUSER_NONUNI_VAL evaluates the\n";
  cout << "    nonuniform Overhauser spline.\n";
  cout << "\n";
  cout << "  In this draft of a test, we simply repeat the work\n";
  cout << "  for the uniform test.\n";

  for ( i = 0; i < NDATA; i++ )
  {
    tdata[i] = ( double ) ( i );
    ydata[i] = sin ( 2.0 * pi * tdata[i] / ( double ) ( NDATA - 1) );
  }

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Number of data values = " << NDATA << "\n";
  cout << "\n";
  cout << "       T             Y\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
    cout                         << "  "
         << setw(12) << tdata[i] << "  "
         << setw(12) << ydata[i] << "\n";
  }

  cout << "\n";
  cout << "    T, Spline(T)\n";
  cout << "\n";

  for ( i = 0; i <= NDATA; i++ )
  {
    if ( i == 0 )
    {
      tlo = tdata[0] - 0.5 * ( tdata[1] - tdata[0] );
      thi = tdata[0];
    }
    else if ( i < NDATA )
    {
      tlo = tdata[i-1];
      thi = tdata[i];
    }
    else if ( NDATA <= i )
    {
      tlo = tdata[NDATA-1];
      thi = tdata[NDATA-1] + 0.5 * ( tdata[NDATA-1] - tdata[NDATA-2] );
    }

    if ( i < NDATA )
    {
      jhi = nsample - 1;
    }
    else
    {
      jhi = nsample;
    }

    for ( j = 0; j <= jhi; j++ )
    {
      tval = ( ( double ) ( nsample - j ) * tlo
             + ( double ) (           j ) * thi )
             / ( double ) ( nsample     );

      yval = spline_overhauser_nonuni_val ( NDATA, tdata, ydata, tval );

      if ( 0 < i & j == 0 )
      {
        mark = '*';
      }
      else
      {
        mark = ' ';
      }
      cout                     << "  "
                       << mark << "  "
           << setw(12) << tval << "  "
           << setw(12) << yval << "\n";

    }

  }

  return;
# undef NDATA
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests SPLINE_OVERHAUSER_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NDATA 4
# define NDIM 1

  int i;
  double tdata[NDATA];
  double tval;
  double ydata[NDIM*NDATA];
  double yval[NDIM];
  double zdata[NDATA];
  double zval[NDIM];

  cout << "\n";
  cout << "TEST23\n";
  cout << "  SPLINE_OVERHAUSER_VAL evaluates the\n";
  cout << "    Overhauser spline.\n";
//
//  Set the data.
//
  tdata[0] = 1.0;
  ydata[0+0*NDIM] =   0.0;
  zdata[0+0*NDIM] =   0.0;

  tdata[1] = 2.0;
  ydata[0+1*NDIM] =   1.0;
  zdata[0+1*NDIM] =   1.0;

  tdata[2] = 3.0;
  ydata[0+2*NDIM] =   2.0;
  zdata[0+2*NDIM] = - 1.0;

  tdata[3] = 4.0;
  ydata[0+3*NDIM] =   3.0;
  zdata[0+3*NDIM] =   0.0;

  cout << "\n";
  cout << "  Data\n";
  cout << "  TDATA[I], YDATA[I], ZDATA[I]\n";
  cout << "\n";
  for ( i = 0; i < NDATA; i++ )
  {
    cout                                << "  "
         << setw(12) << tdata[i]        << "  "
         << setw(12) << ydata[0+i*NDIM] << "  "
         << setw(12) << zdata[0+i*NDIM] << "\n";
  }
//
//  Now evaluate the spline all over the place.
//
  cout << "\n";
  cout << "  T, Spline value\n";
  cout << "\n";

  for ( i = 0; i <= 6 * NDATA + 3; i++ )
  {
    tval = ( ( double ) i ) / 6.0;
    spline_overhauser_val ( NDIM, NDATA, tdata, ydata, tval, yval );
    spline_overhauser_val ( NDIM, NDATA, tdata, zdata, tval, zval );

    cout                        << "  "
         << setw(12) << tval    << "  "
         << setw(12) << yval[0] << "  "
         << setw(12) << zval[0] << "\n";
  }

  return;
# undef NDATA
# undef NDIM
}
//****************************************************************************80

void test235 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST235 tests SPLINE_PCHIP_SET and SPLINE_PCHIP_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 21
# define NE 101

  double d[N];
  double diff;
  double f[N];
  double fd[NE];
  double fe[NE];
  int i;
  double x[N];
  double xe[NE];

  cout << "\n";
  cout << "TEST235\n";
  cout << "  SPLINE_PCHIP_SET sets up a piecewise cubic\n";
  cout << "    Hermite interpolant.\n";
  cout << "  SPLINE_PCHIP_VAL evaluates the interpolant.\n";
  cout << "\n";
//
//  Compute Runge's function at N points in [-1,1].
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = -1.0 + ( double ) ( i ) / 10.0;
    f[i] = frunge ( x[i] );
  }
//
//  SPLINE_PCHIP_SET takes the data in X and F, and constructs a table in D
//  that defines the interpolant.
//
  spline_pchip_set ( N, x, f, d );
//
//  Evaluate the interpolant and derivative at NE points from -1 to 0.
//
  for ( i = 0; i < NE; i++ )
  {
    xe[i] = -1.0 + ( double ) ( i ) / ( double ) ( NE - 1 );
  }

  spline_pchip_val ( N, x, f, d, NE, xe, fe );
//
//  Print the table of X, F(exact) and F(interpolated)
//
  for ( i = 0; i < NE; i++ )
  {
    diff = fe[i] - frunge ( xe[i] );

    cout << "  " << setw(8)  << xe[i]
         << "  " << setw(10) << frunge ( xe[i] )
         << "  " << setw(10) << fe[i]
         << "  " << setw(14) << diff << "\n";
  }

  return;
# undef N
# undef NE
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests SPLINE_QUADRATIC_VAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 11

  double fval;
  int i;
  int j;
  int jhi;
  double t[N];
  double tval;
  double y[N];
  double ypval;
  double yval;

  cout << "\n";
  cout << "TEST24\n";
  cout << "  SPLINE_QUADRATIC_VAL evaluates a\n";
  cout << "    quadratic spline.\n";
  cout << "\n";
  cout << "  Runge''s function, evenly spaced knots.\n";

  for ( i = 0; i < N; i++ )
  {
    t[i] =  ( ( double ) ( N - i - 1 ) * (-1.0)
            + ( double ) (     i     ) * (+1.0) )
            / ( double ) ( N     - 1 );
    y[i] =  frunge ( t[i] );
  }

  cout << "\n";
  cout << "  The data to be interpolated:\n";
  cout << "\n";
  cout << "  Number of data values = " << N << "\n";
  cout << "\n";
  cout << "       T             Y\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout                     << "  "
         << setw(12) << t[i] << "  "
         << setw(12) << y[i] << "\n";
  }

  cout << "\n";
  cout << "  Interpolated values\n";
  cout << "\n";
  cout << "       T             Y           Y(exact)\n";
  cout << "\n";

  for ( i = 0; i <= N; i++ )
  {
    if ( i == 0 )
    {
      jhi = 1;
    }
    else if ( i < N )
    {
      jhi = 2;
    }
    else
    {
      jhi = 2;
    }

    for ( j = 1; j <= jhi; j++ )
    {
      if ( i == 0 )
      {
        tval = t[0] - 1.0;
      }
      else if ( i < N )
      {
        tval = ( ( double ) ( jhi - j + 1 ) * t[i-1]
               + ( double ) (       j - 1 ) * t[i] )
               / ( double ) ( jhi         );
      }
      else
      {
        if ( j == 1 )
        {
          tval = t[N-1];
        }
        else
        {
          tval = t[N-1] + 1.0;
        }
      }

      spline_quadratic_val ( N, t, y, tval, &yval, &ypval );

      fval = frunge ( tval );

      cout << setw(12) << tval << "  "
           << setw(12) << yval << "  "
           << setw(12) << fval << "\n";
    }
  }

  return;
# undef N
}
//****************************************************************************80

void parabola_formula ( double x, double *y, double *yp, double *ypp )

//****************************************************************************80
//
//  Purpose:
//
//    PARABOLA_FORMULA evaluates a parabola for us.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  *y = 2.0 * x * x + 3.0 * x + 1.0;
  *yp = 4.0 * x + 3.0;
  *ypp = 4.0;

  return;
}
//****************************************************************************80

double frunge ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FRUNGE sets the Runge data values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;

  fx = 1.0 / ( 1.0 + 25.0 * x * x );

  return fx;
}
//****************************************************************************80

double fprunge ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FPRUNGE sets the Runge derivative values at the endpoints.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bot;
  double fx;

  bot = 1.0 + 25.0 * x * x;
  fx = -50.0 * x / ( bot * bot );

  return fx;
}
//****************************************************************************80

double fpprunge ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FPPRUNGE sets the Runge second derivative values at the endpoints.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bot;
  double fx;

  bot = 1.0 + 25.0 * x * x;
  fx = ( -50.0 + 3750.0 * x * x ) / ( bot * bot * bot );

  return fx;
}
//****************************************************************************80

double fcube ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FCUBE evaluates a cubic function.
//
//  Discussion:
//
//    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which the function is evaluated.
//
//    Output, double FCUBE, the value of the function.
//
{
  double fx;

  fx = ( ( x + 2.0 ) * x + 3.0 ) * x + 4.0;

  return fx;
}
//****************************************************************************80

double fpcube ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FPCUBE evaluates the derivative of a cubic function.
//
//  Discussion:
//
//    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which the function is evaluated.
//
//    Output, double FPCUBE, the value of the derivative of the function.
//
{
  double fx;

  fx = ( 3.0 * x + 4.0 ) * x + 3.0;

  return fx;
}
//****************************************************************************80

double fppcube ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FPPCUBE evaluates the second derivative of a cubic function.
//
//  Discussion:
//
//    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which the function is evaluated.
//
//    Output, double FPPCUBE, the value of the second derivative of the function.
//
{
  double fx;

  fx = 6.0 * x + 4.0;

  return fx;
}
