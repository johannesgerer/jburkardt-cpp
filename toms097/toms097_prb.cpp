# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "toms097.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TOMS097_PRB.
//
//  Discussion:
//
//    TOMS097_PRB tests the TOMS097 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TOMS097_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TOMS097 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TOMS097_PRB\n";
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
//    TEST01 tests I4MAT_SHORTEST_PATH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2014
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  int a[N*N] = {
     0, -1, -1, -1, -1, -1, 
     2,  0, -1, -1, -1,  5, 
     5,  7,  0, -1,  2, -1, 
    -1,  1,  4,  0, -1,  2, 
    -1, -1, -1,  3,  0,  4, 
    -1,  8, -1, -1,  3,  0  };
  int i;
  int j;
  int n = N;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  I4MAT_SHORTEST_PATH uses Floyd''s algorithm to find the\n";
  cout << "  shortest distance between all pairs of nodes\n";
  cout << "  in a directed graph, starting from the initial array\n";
  cout << "  of direct node-to-node distances.\n";

  cout << "\n";
  cout << "  In the initial direct distance array, if\n";
  cout << "    A(I,J) = HUGE,\n";
  cout << "  this indicates there is NO directed link from\n";
  cout << "  node I to node J.  In that case, the value of\n";
  cout << "  of A(I,J) is essentially 'infinity'.\n";

  cout << "\n";
  cout << "  Initial direct-link distance matrix:\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << setw(6) << a[i+j*n];
    }
    cout << "\n";
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == -1 )
      {
        a[i+j*n] = i4_huge ( );
      }
    }
  } 

  i4mat_shortest_path ( n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == i4_huge ( ) )
      {
        a[i+j*n] = -1;
      }
    }
  }

  cout << "\n";
  cout << "  In the final shortest distance array, if\n";
  cout << "    A(I,J) = -1,\n";
  cout << "  this indicates there is NO directed path from\n";
  cout << "  node I to node J.\n";

  cout << "\n";
  cout << "  Final distance matrix:\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << setw(6) << a[i+j*n];
    }
    cout << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R8MAT_SHORTEST_PATH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2014
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  double a[N*N] = {
     0.0, -1.0, -1.0, -1.0, -1.0, -1.0, 
     2.0,  0.0, -1.0, -1.0, -1.0,  5.0, 
     5.0,  7.0,  0.0, -1.0,  2.0, -1.0, 
    -1.0,  1.0,  4.0,  0.0, -1.0,  2.0, 
    -1.0, -1.0, -1.0,  3.0,  0.0,  4.0, 
    -1.0,  8.0, -1.0, -1.0,  3.0,  0.0  };
  int i;
  int j;
  int n = N;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R8MAT_SHORTEST_PATH uses Floyd''s algorithm to find the\n";
  cout << "  shortest distance between all pairs of nodes\n";
  cout << "  in a directed graph, starting from the initial array\n";
  cout << "  of direct node-to-node distances.\n";

  cout << "\n";
  cout << "  In the initial direct distance array, if\n";
  cout << "    A(I,J) = -1,\n";
  cout << "  this indicates there is NO directed link from\n";
  cout << "  node I to node J.  In that case, the value of\n";
  cout << "  of A(I,J) is essentially 'infinity'.\n";

  cout << "\n";
  cout << "  Initial direct-link distance matrix:\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << setw(10) << a[i+j*n];
    }
    cout << "\n";
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == -1.0 )
      {
        a[i+j*n] = r8_huge ( );
      }
    }
  } 

  r8mat_shortest_path ( n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == r8_huge ( ) )
      { 
        a[i+j*n] = -1.0;
      }
    }
  }

  cout << "\n";
  cout << "  In the final shortest distance array, if\n";
  cout << "    A(I,J) = -1,\n";
  cout << "  this indicates there is NO directed path from\n";
  cout << "  node I to node J.\n";

  cout << "\n";
  cout << "  Final distance matrix:\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      cout << setw(10) << a[i+j*n];
    }
    cout << "\n";
  }

  return;
}
