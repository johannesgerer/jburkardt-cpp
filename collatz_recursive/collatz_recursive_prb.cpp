# include <stdlib.h>
# include <stdio.h>

# include "collatz_recursive.h"

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for COLLATZ_RECURSIVE_PRB.

  Discussion:

    COLLATZ_RECURSIVE_PRB tests the COLLATZ_RECURSIVE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "COLLATZ_RECURSIVE_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the COLLATZ_RECURSIVE library.\n" );

  test01 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "COLLATZ_RECURSIVE_PRB\n" );
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

    TEST01 tests COLLATZ_PATH;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2012

  Author:

    John Burkardt
*/
{
  int n;
  int n_test[5] = { 7, 8, 9, 10, 600 };
  int test;
  int test_num = 5;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  COLLATZ_PATH prints the members of a Collatz path.\n" );

  for ( test = 0; test < test_num; test++ )
  {
    n = n_test[test];
    printf ( "\n" );
    printf ( "  %d is the starting point.\n", n );
    collatz_path ( n );
  }

  return;
}
