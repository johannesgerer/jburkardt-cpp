# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "asa299.hpp"

int main ( void );

void test01 ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA299_PRB.
//
//  Discussion:
//
//    ASA299_PRB calls the ASA299 test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "ASA299_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA299 library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA299_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests SIMPLEX_LATTICE_POINT_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  int i;
  int j;
  bool more;
  int t = 4;
  int x[N];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  SIMPLEX_LATTICE_POINT_NEXT generates lattice points\n";
  cout << "  in the simplex\n";
  cout << "    0 <= X\n";
  cout << "    sum ( X(1:N) ) <= T\n";
  cout << "  Here N = " << N << "\n";
  cout << "  and T =  " << t << "\n";
  cout << "\n";
  cout << "     Index        X(1)      X(2)      X(3)      X(4)\n";
  cout << "\n";

  more = false;

  i = 0;

  for ( ; ; )
  {
    simplex_lattice_point_next ( N, t, &more, x );

    i = i + 1;

    cout << "  " << setw(8) << i;
    cout << "  ";
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(8) << x[j];
    }
    cout << "\n";

    if ( !more )
    {
      break;
    }

  }
 
  return;
# undef N
}
