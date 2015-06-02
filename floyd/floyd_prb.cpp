# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "floyd.hpp"

int main ( );

void test01 ( );
void test02 ( );
void test03 ( );
double test03_sub ( int n );
void test04 ( );

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
//    FLOYD_PRB tests the FLOYD library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FLOYD_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FLOYD library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
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

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 applies Floyd's method to problems of increasing size.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double wtime;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Test I4MAT_FLOYD on the MOD(I,J) matrix.\n";
  cout << "  The work is roughly N^3.\n";
  cout << "\n";
  cout << "         N   Time (seconds)  Time/N^3\n";
  cout << "\n";

  n = 1;
  while ( n <= 2048 )
  {
    wtime = test03_sub ( n );
    cout << "  " << setw(8) << n
         << "  " << setw(14) << wtime
         << "  " << setw(14) << 1000000.0 * wtime / ( double ) ( n )
           / ( double ) ( n ) / ( double ) ( n )
         << "\n";
    n = n * 2;
  }

  return ;
}
//****************************************************************************80

double test03_sub ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03_SUB tests I4MAT_FLOYD.
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
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 uses Floyd's method for a triangulation.
//
//  Discussion:
//
//     8  41--42--43--44  45--46--47--48
//     |   | \ | \ | \ |   | \ | \ | \ |
//     7  33--34--35--36  37--38--39--40
//     |   | \ |                   | \ |
//     6  29--30                  31--32
//     |   | \ |                   | \ |
//     5  25--26                  27--28
//     |   | \ |                   | \ |
//     4  21--22                  23--24
//     |   | \ |                   | \ |
//     3  17--18                  19--20
//     |   | \ |                   | \ |
//     2   9--10--11--12--13--14--15--16
//     |   | \ | \ | \ | \ | \ | \ | \ |
//     1   1---2---3---4---5---6---7---8
//     |    
//     +---1---2---3---4---5---6---7---8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2014
//
//  Author:
//
//    John Burkardt
//
{
# define ELEMENT_NUM 46
# define NODE_NUM 48
  
  double d[NODE_NUM*NODE_NUM];
  int element;
  int element_node[3*ELEMENT_NUM] = {
     1,  2,  9, 
     2, 10,  9, 
     2,  3, 10, 
     3, 11, 10, 
     3,  4, 11, 
     4, 12, 11, 
     4,  5, 12, 
     5, 13, 12, 
     5,  6, 13, 
     6, 14, 13, 
     6,  7, 14, 
     7, 15, 14, 
     7,  8, 15, 
     8, 16, 15, 
     9, 10, 17, 
    10, 18, 17, 
    15, 16, 19, 
    16, 20, 19, 
    17, 18, 21, 
    18, 22, 21, 
    19, 20, 23, 
    20, 24, 23, 
    21, 22, 25, 
    22, 26, 25, 
    23, 24, 27, 
    24, 28, 27, 
    25, 26, 29, 
    26, 30, 29, 
    27, 28, 31, 
    28, 32, 31, 
    29, 30, 33, 
    30, 34, 33, 
    31, 32, 39, 
    32, 40, 39, 
    33, 34, 41, 
    34, 42, 41, 
    34, 35, 42, 
    35, 43, 42, 
    35, 36, 43, 
    36, 44, 43, 
    37, 38, 45, 
    38, 46, 45, 
    38, 39, 46, 
    39, 47, 46, 
    39, 40, 47, 
    40, 48, 47 };
  int element_num = ELEMENT_NUM;
  int i;
  int j;
  int n1;
  int n2;
  int node_num = NODE_NUM;
  double xy[2*NODE_NUM] = {
    1.0, 1.0, 
    2.0, 1.0, 
    3.0, 1.0, 
    4.0, 1.0, 
    5.0, 1.0, 
    6.0, 1.0, 
    7.0, 1.0, 
    8.0, 1.0, 
    1.0, 2.0, 
    2.0, 2.0, 
    3.0, 2.0, 
    4.0, 2.0, 
    5.0, 2.0, 
    6.0, 2.0, 
    7.0, 2.0, 
    8.0, 2.0, 
    1.0, 3.0,  
    2.0, 3.0, 
    7.0, 3.0, 
    8.0, 3.0, 
    1.0, 4.0, 
    2.0, 4.0, 
    7.0, 4.0, 
    8.0, 4.0, 
    1.0, 5.0, 
    2.0, 5.0, 
    7.0, 5.0, 
    8.0, 5.0, 
    1.0, 6.0, 
    2.0, 6.0, 
    7.0, 6.0, 
    8.0, 6.0, 
    1.0, 7.0, 
    2.0, 7.0, 
    3.0, 7.0, 
    4.0, 7.0, 
    5.0, 7.0, 
    6.0, 7.0, 
    7.0, 7.0, 
    8.0, 7.0, 
    1.0, 8.0,  
    2.0, 8.0, 
    3.0, 8.0, 
    4.0, 8.0, 
    5.0, 8.0, 
    6.0, 8.0, 
    7.0, 8.0, 
    8.0, 8.0 };

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Start with a triangulation, and use R8_FLOYD\n";
  cout << "  to determine the pairwise distance matrix.\n";
//
//  Must initialize distances to -1
//
  for ( j = 0; j < node_num; j++ )
  {
    for ( i = 0; i < node_num; i++ )
    {
      d[i+j*node_num] = -1.0;
    }
  }
//
//  Diagonals are 0.
//
  for ( i = 0; i < node_num; i++ )
  {
    d[i+i*node_num] = 0.0;
  }
//
//  Initialize D to the one-step distance.
//
  for ( element = 0; element < element_num; element++ )
  {
    n1 = element_node[2+element*3] - 1;
    for ( i = 0; i < 3; i++ )
    {
      n2 = element_node[i+element*3] - 1;
      d[n1+n2*node_num] = sqrt ( pow ( xy[0+n1*2] - xy[0+n2*2], 2 )
                               + pow ( xy[1+n1*2] - xy[1+n2*2], 2 ) );
      d[n2+n1*node_num] = d[n1+n2*node_num];
      n1 = n2;
    }
  }
//
//  Reset -1 values to R8_HUGE.
//
  for ( j = 0; j < node_num; j++ )
  {
    for ( i = 0; i < node_num; i++ )
    {
      if ( d[i+j*node_num] == -1.0 )
      {
        d[i+j*node_num] = r8_huge ( );
      }
    }
  }
//
//  Update D to the N-1 step distance.
//
  r8mat_floyd ( node_num, d );

  r8mat_print ( node_num, node_num, d, "  Distance matrix" );

  return;
# undef ELEMENT_NUM
# undef NODE_NUM
}
