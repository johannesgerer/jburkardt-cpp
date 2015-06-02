# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "bellman_ford.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BELLMAN_FORD_PRB.
//
//  Discussion:
//
//    BELLMAN_FORD_PRB tests the BELLMAN_FORD library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BELLMAN_FORD_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BELLMAN_FORD library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BELLMAN_FORD_PRB\n";
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
//    TEST01 runs a simple test.
//
//  Discussion:
//
//    The correct distances are { 0, -6, -2, 3, 0, 0 }.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  int e[2*10] = 
  {
    1, 0,
    4, 1,
    1, 2,
    2, 4,
    4, 0,
    2, 5,
    5, 0,
    3, 2,
    5, 3,
    3, 0,
  };
  int e_num = 10;
  double e_weight[10] = 
  {
    -3.0,
     6.0,
    -4.0,
    -1.0,
     4.0,
    -2.0,
     2.0,
     8.0,
    -3.0,
     3.0
  };
  int predecessor[6];
  int source = 0;
  int v_num = 6;
  double v_weight[6];
  
  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Bellman-Ford shortest path algorithm.\n";

  cout << "\n";
  cout << "  Number of vertices = " << v_num << "\n";
  cout << "  Number of edges = " << e_num << "\n";
  cout << "  The reference vertex is " << source << "\n";
  i4mat_transpose_print ( 2, e_num, e, "  The edge array:" );
  r8vec_print ( e_num, e_weight, "  The edge weights:" );

  bellman_ford ( v_num, e_num, source, e, e_weight, v_weight, predecessor );

  r8vec_print ( v_num, v_weight, "  The shortest distances:" );

  i4vec_print ( v_num, predecessor, "  The vertex predecessor parents for the shortest paths:" );

  return;
}