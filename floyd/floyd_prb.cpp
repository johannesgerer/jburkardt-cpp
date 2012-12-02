# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "floyd.hpp"

int main ( );

void test01 ( );
void test02 ( );
double test03 ( int n );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FLOYD_PRB.
//
//  Discussion:
//
//    FLOYD_PRB calls a set of problems for FLOYD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2009
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double wtime;

  timestamp ( );

  cout << "\n";
  cout << "FLOYD_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FLOYD library.\n";

  test01 ( );
  test02 ( );

  cout << "\n";
  cout << "FLOYD_TEST03\n";
  cout << "  Test I4MAT_FLOYD on the MOD(I,J) matrix.\n";
  cout << "  The work is roughly N^3.\n";
  cout << "\n";
  cout << "         N   Time (seconds)  Time/N^3\n";
  cout << "\n";

  n = 1;
  while ( n <= 2048 )
  {
    wtime = test03 ( n );
    cout << "  " << setw(8) << n
         << "  " << setw(14) << wtime
         << "  " << setw(14) << 1000000.0 * wtime / ( double ) ( n )
           / ( double ) ( n ) / ( double ) ( n )
         << "\n";
    n = n * 2;
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "FLOYD_PRB\n";
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
//    TEST01 tests I4MAT_FLOYD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 November 2008
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
  int huge;
  int i;
  int j;
  int n = N;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  I4MAT_FLOYO uses Floyd's algorithm to find the\n";
  cout << "  shortest distance between all pairs of nodes\n";
  cout << "  in a directed graph, starting from the initial array\n";
  cout << "  of direct node-to-node distances.\n";

  cout << "\n";
  cout << "  In the initial direct distance array, if\n";
  cout << "    A(I,J) = -1,\n";
  cout << "  this indicates there is NO directed link from\n";
  cout << "  node I to node J.  In that case, the value of\n";
  cout << "  of A(I,J) is essentially \"infinity\".\n";

  i4mat_print ( n, n, a, "  Initial direct distance array:" );

  huge = i4_huge ( ) / 2;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == - 1 )
      {
        a[i+j*n] = huge;
      }
    }
  }

  i4mat_floyd ( n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == huge )
      {
        a[i+j*n] = - 1;
      }
    }
  }

  cout << "\n";
  cout << "  In the final shortest distance array, if\n";
  cout << "    A(I,J) = -1,\n";
  cout << "  this indicates there is NO directed path from\n";
  cout << "  node I to node J.\n";

  i4mat_print ( n, n, a, "  Final shortest distance array:" );

  return;
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R8MAT_FLOYD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 November 2008
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
  double huge;
  int i;
  int j;
  int n = N;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R8MAT_FLOYO uses Floyd's algorithm to find the\n";
  cout << "  shortest distance between all pairs of nodes\n";
  cout << "  in a directed graph, starting from the initial array\n";
  cout << "  of direct node-to-node distances.\n";

  cout << "\n";
  cout << "  In the initial direct distance array, if\n";
  cout << "    A(I,J) = -1,\n";
  cout << "  this indicates there is NO directed link from\n";
  cout << "  node I to node J.  In that case, the value of\n";
  cout << "  of A(I,J) is essentially \"infinity\".\n";

  r8mat_print ( n, n, a, "  Initial direct distance array:" );

  huge = r8_huge ( );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == - 1.0 )
      {
        a[i+j*n] = huge;
      }
    }
  }

  r8mat_floyd ( n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*n] == huge )
      {
        a[i+j*n] = - 1.0;
      }
    }
  }

  cout << "\n";
  cout << "  In the final shortest distance array, if\n";
  cout << "    A(I,J) = -1,\n";
  cout << "  this indicates there is NO directed path from\n";
  cout << "  node I to node J.\n";

  r8mat_print ( n, n, a, "  Final shortest distance array:" );

  return;
# undef N
}
//****************************************************************************80

double test03 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests I4MAT_FLOYD.
//
//  Discussion:
//
//    The matrix size is input by the user.
//
//    The matrix A has the property that
//
//      A(I,J) = 1 if I is divisible by J.
//
//  Example:
//
//    N = 6
//
//    1 0 0 0 0 0
//    1 1 0 0 0 0
//    1 0 1 0 0 0
//    1 1 0 1 0 0
//    1 0 0 0 1 0
//    1 1 1 0 0 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the matrix.
//
//    Output, double TEST03, the CPU time required by I4MAT_FLOYD.
//
{
  int *a;
  int huge;
  int i;
  int j;
  double time1;
  double time2;
  double wtime;

  a = new int[n*n];

  huge = i4_huge ( ) / 2;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( ( i + 1 ) % ( j + 1 ) == 0 )
      {
        a[i+j*n] = 1;
      }
      else
      {
        a[i+j*n] = huge;
      }
    }
  }

  time1 = cpu_time ( );

  i4mat_floyd ( n, a );

  time2 = cpu_time ( );

  wtime = time2 - time1;

  delete [] a;

  return wtime;
}
