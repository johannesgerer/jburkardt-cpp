# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>

# include "qls.h"
# include "qr_solve.h"
# include "r8lib.h"

int main ( );
void test01 ( );
void test02 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for QLS_PRB.

  Discussion:

    QLS_PRB tests the QUADRATURE_LEAST_SQUARES library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 March 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "QLS_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the QUADRATURE_LEAST_SQUARES library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "QLS_PRB\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    QLS_TEST01 shows that we can compute the Newton-Cotes rules.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 April 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int d;
  int i;
  int n;
  double *w1;
  double *w2;
  double *x;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  WEIGHTS_LS computes the weights for a\n" );
  printf ( "  least squares quadrature rule.\n" );
/*
  Demonstrate the 5 point Newton-Cotes closed rule.
*/
  printf ( "\n" );
  printf ( "  W1 = classical Newton Cotes weights, N = 5\n" );
  printf ( "  W2 = least squares weights, D = 4, N = 5\n" );

  n = 5;
  x = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );

  ncc_set ( n, x, w1 );
/*
  Using the same points, compute the least squares weights
  for polynomial approximation up to degree 4.
*/
  d = n - 1;
  a = -1.0;
  b = +1.0;

  w2 = weights_ls ( d, a, b, n, x );

  printf ( "\n" );
  printf ( "   I        X(i)          W1(i)           W2(i)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10.4f  %14.6g  %14.6g\n", i, x[i], w1[i], w2[i] );
  }

  free ( w1 );
  free ( w2 );
  free ( x );
/*
  Look at a 9 point rule.
  Note that Newton Cotes rules soon have many negative weights.
*/
  printf ( "\n" );
  printf ( "  W1 = classical Newton Cotes weights, N = 9\n" );
  printf ( "  W2 = least squares weights, D = 4, N = 9\n" );

  n = 9;
  x = ( double * ) malloc ( n * sizeof ( double ) );
  w1 = ( double * ) malloc ( n * sizeof ( double ) );
 
  ncc_set ( n, x, w1 );
/*
  Using the same points, compute the least squares weights
  for polynomial approximation up to degree 4.
*/
  d = 4;
  a = -1.0;
  b = +1.0;
  w2 = weights_ls ( d, a, b, n, x );

  printf ( "\n" );
  printf ( "   I        X(i)          W1(i)           W2(i)\n" );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %2d  %10.4f  %14.6g  %14.6g\n", i, x[i], w1[i], w2[i] );
  }

  free ( w1 );
  free ( w2 );
  free ( x );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 uses random points as abscissas.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 April 2014

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int d;
  double e;
  double exact;
  double *f;
  int i;
  int n = 50;
  double q;
  int seed;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  WEIGHTS_LS computes the weights for a\n" );
  printf ( "  least squares quadrature rule.\n" );
  printf ( "\n" );
  printf ( "  Pick 50 random values in [-1,+1].\n" );
  printf ( "  Compare Monte Carlo (equal weight) integral estimate\n" );
  printf ( "  to least squares estimates of degree D = 0, 1, 2, 3, 4.\n" );
  printf ( "  For low values of D, the least squares estimate improves.\n" );
  printf ( "\n" );
  printf ( "  As D increases, the estimate can deteriorate.\n" );
/*
  Define the integration interval.
*/
  a = -5.0;
  b = +5.0;
/*
  Get random values.
*/
  seed = 123456789;
  x = r8vec_uniform_ab_new ( n, a, b, &seed );
/*
  Evaluate the function.
*/
  f = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0 / ( 1.0 + x[i] * x[i] );
  }
  exact = atan ( b ) - atan ( a );

  printf ( "\n" );
  printf ( "  Rule         Estimate          Error\n" );
/*
  Get the MC estimate.
*/
  q = ( b - a ) * r8vec_sum ( n, f ) / ( double ) ( n );
  e = fabs ( q - exact );

  printf ( "\n" );
  printf ( "  MC     %14.6g  %14.6g\n", q, e );
  printf ( "\n" );
/*
  Using the same points, compute the least squares weights
  for polynomial approximation of degree D.
*/
  for ( d = 0; d <= 15; d++ )
  {
    w = weights_ls ( d, a, b, n, x );
    q = r8vec_dot_product ( n, w, f );
    e = fabs ( q - exact );
    printf ( "  LS%2d  %14.6g  %14.6g\n", d, q, e );
    free ( w );
  }

  q = exact;
  e = fabs ( q - exact );
  printf ( "\n" );
  printf ( "  EXACT %14.6g  %14.6g\n", q, e );

  free ( f );
  free ( x );

  return;
}
