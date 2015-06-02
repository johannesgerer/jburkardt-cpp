# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "asa159.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA159_PRB.
//
//  Discussion:
//
//    ASA159_PRB tests the ASA159 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA159_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA159 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA159_PRB\n";
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
//    TEST01 tests RCONT2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5

  int a[M*N];
  int c[N] = { 2, 2, 2, 2, 1 };
  int i;
  int ierror;
  bool key = false;
  int m = M;
  int n = N;
  int ntest = 10;
  int r[M] = { 3, 2, 2, 1, 1 };
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  RCONT2 constructs a random matrix with\n";
  cout << "  given row and column sums.\n";

  i4vec_print ( m, r, "  The rowsum vector:" );
  i4vec_print ( n, c, "  The columnsum vector:" );

  for ( i = 1; i <= ntest; i++ )
  {
    rcont2 ( m, n, r, c, &key, &seed, a, &ierror );

    if ( ierror != 0 )
    {
      cout << "\n";
      cout << "  RCONT2 returned error flag IERROR = " << ierror << "\n";
      return;
    }

    i4mat_print ( m, n, a, "  The rowcolsum matrix:" );
  }

  return;
}
