# include <stdlib.h>
# include <stdio.h>

# include "toeplitz_cholesky.h"

int main ( );
void toeplitz_cholesky_prb01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN tests TOEPLITZ_CHOLESKY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "TOEPLITZ_CHOLESKY_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TOEPLITZ_CHOLESKY library.\n" );

  toeplitz_cholesky_prb01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TOEPLITZ_CHOLESKY_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void toeplitz_cholesky_prb01 ( )

/******************************************************************************/
/*
  Purpose:

    TOEPLITZ_CHOLESKY_PRB01 tests TOEPLITZ_CHOLESKY.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 November 2012

  Author:

    John Burkardt
*/
{
  double a[3*3] = {  
    1.0,    0.5, -0.375,
    0.5,    1.0,  0.5,
   -0.375,  0.5,  1.0 };
  double *b;
  double *g;
  double g_save[2*3] = {
    1.0,    0.0,
    0.5,    0.5,
   -0.375, -0.375 };
  double *l;
  int n = 3;
  double *r;

  printf ( "\n" );
  printf ( "TOEPLITZ_CHOLESKY_PRB01\n" );
  printf ( "  Test the factorization of a simple example.\n" );
/*
  TOEPLITZ_CHOLESKY_UPPER.
*/
  printf ( "\n" );
  printf ( "TOEPLITZ_CHOLESKY_UPPER:\n" );

  r8mat_print ( n, n, a, "  Toeplitz matrix A:" );

  r = toeplitz_cholesky_upper ( n, a );
  r8mat_print ( n, n, r, "  Computed upper Cholesky factor R:" );

  b = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, b, "  Product R'R:" );

  free ( b );
  free ( r );
/*
  TOEP_CHOLESKY_UPPER.
*/
  printf ( "\n" );
  printf ( "TOEP_CHOLESKY_UPPER:\n" );
  
  g = r8mat_copy_new ( 2, n, g_save );

  r8mat_print ( 2, n, g, "  Compressed Toeplitz matrix G:" );

  r = toep_cholesky_upper ( n, g );
  r8mat_print ( n, n, r, "  Computed upper Cholesky factor R:" );

  b = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, b, "  Product R'R:" );

  free ( b );
  free ( g );
  free ( r );
/*
  TOEPLITZ_CHOLESKY_LOWER.
*/
  printf ( "\n" );
  printf ( "TOEPLITZ_CHOLESKY_LOWER:\n" );

  r8mat_print ( n, n, a, "  Toeplitz matrix A:" );

  l = toeplitz_cholesky_lower ( n, a );
  r8mat_print ( n, n, l, "  Computed lower Cholesky factor L:" );

  b = r8mat_mmt_new ( n, n, n, l, l );
  r8mat_print ( n, n, b, "  Product LL':" );

  free ( b );
  free ( l );
/*
  TOEP_CHOLESKY_LOWER.
*/
  printf ( "\n" );
  printf ( "TOEP_CHOLESKY_LOWER:\n" );

  g = r8mat_copy_new ( 2, n, g_save );

  r8mat_print ( 2, n, g, "  Compressed Toeplitz matrix G:" );

  l = toep_cholesky_lower ( n, g );
  r8mat_print ( n, n, l, "  Computed lower Cholesky factor L:" );

  b = r8mat_mmt_new ( n, n, n, l, l );
  r8mat_print ( n, n, b, "  Product LL':" );

  free ( b );
  free ( g );
  free ( l );

  return;
}



