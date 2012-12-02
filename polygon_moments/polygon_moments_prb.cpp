# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "polygon_moments.h"

int main ( );
void test01 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    POLYGON_MOMENTS_PRB tests POLYGON_MOMENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 October 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "POLYGON_MOMENTS_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test POLYGON_MOMENTS library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "POLYGON_MOMENTS_PRB:\n" );
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

    TEST01 carries out a test on a rectangle.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 October 2012

  Author:

    John Burkardt
*/
{
  double alpha_exact[6] = {
    1.0, 
    5.0, 4.0, 
    30.66666666666667, 22.0, 18.66666666666666 };
  double alpha_pq;
  int k;
  double mu_exact[6] = {
    1.0, 
    0.0, 0.0, 
    5.666666666666667, 2.0, 2.666666666666667 };
  double mu_pq;
  int n = 4;
  double nu_exact[6] = {
    40.0, 
    200.0, 160.0, 
    1226.66666666666667, 880.0, 746.66666666666666 };
  double nu_pq;
  int p;
  int q;
  int s;
  double x[4] = {
    2.0, 10.0, 8.0, 0.0 };
  double y[4] = {
    0.0,  4.0, 8.0, 4.0 };

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Check normalized moments of a rectangle.\n" );
  printf ( "\n" );
  printf ( "   P   Q             Nu(P,Q)\n" );
  printf ( "            Computed         Exact\n" );
  printf ( "\n" );
  k = 0;
  for ( s = 0; s <= 2; s++ )
  {
    for ( p = s; 0 <= p; p-- )
    {
      q = s - p;
      nu_pq = moment ( n, x, y, p, q );
      printf ( "  %2d  %2d  %14.6g  %14.6g\n", p, q, nu_pq, nu_exact[k] );
      k = k + 1;
    }
  }

  printf ( "\n" );
  printf ( "   P   Q           Alpha(P,Q)\n" );
  printf ( "            Computed         Exact\n" );
  printf ( "\n" );
  k = 0;
  for ( s = 0; s <= 2; s++ )
  {
    for ( p = s; 0 <= p; p-- )
    {
      q = s - p;
      alpha_pq = moment_normalized ( n, x,y, p, q );
      printf ( "  %2d  %2d  %14.6g  %14.6g\n", p, q, alpha_pq, alpha_exact[k] );
      k = k + 1;
    }
  }

  printf ( "\n" );
  printf ( "   P   Q             Mu(P,Q)\n" );
  printf ( "            Computed         Exact\n" );
  printf ( "\n" );
  k = 0;
  for ( s = 0; s <= 2; s++ )
  {
    for ( p = s; 0 <= p; p-- )
    {
      q = s - p;
      mu_pq = moment_central ( n, x, y , p, q );
      printf ( "  %2d  %2d  %14.6g  %14.6g\n", p, q, mu_pq, mu_exact[k] );
      k = k + 1;
    }
  }

  return;
}
