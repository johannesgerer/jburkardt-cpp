# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "solve.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Discussion:
//
//    MAIN is the main program for SOLVE_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SOLVE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SOLVE library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SOLVE_PRB\n";
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
//    TEST01 demonstrates how a 3X3 linear system can be set up and solved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double **a;
  double *b;
  int n;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Set up a linear system, and solve it by calling a function.\n";
  cout << "  The linear system is to be accessed using the A[i][j] notation.\n";
  cout << "  It is usually difficult to use this notation and still be able\n";
  cout << "  to pass the array A to a function.  An ordinary doubly indexed\n";
  cout << "  array would require the receiving function to know, IN ADVANCE,\n";
  cout << "  the exact fixed value of the second dimension of A, which defeats\n";
  cout << "  the goal of writing general usable library software.\n";
  cout << "\n";
  cout << "  Here, I think I have made it easy, at the cost of:\n";
  cout << "  * declaring the array as a 'double **a'\n";
  cout << "  * creating the array with r8rmat_new() or r8rmat_zero()\n";
  cout << "  * solving the system with r8rmat_fs_new()\n";
  cout << "  * deleting the array with r8rmat_delete()\n";
//
//  Define the array size.
//
  n = 3;
//
//  Create the array that will contain the matrix.
//
  a = r8rmat_new ( n, n );
//
//  Set the array values.
//
  a[0][0] = 1;
  a[0][1] = 2;
  a[0][2] = 3;

  a[1][0] = 4;
  a[1][1] = 5;
  a[1][2] = 6;

  a[2][0] = 7;
  a[2][1] = 8;
  a[2][2] = 0;
//
//  Create the right hand side.
//
  b = new double[n];
//
//  Set the right hand side values.
//
  b[0] = 14; 
  b[1] = 32;
  b[2] = 23;
//
//  Request the solution of A*x=b.
//
  x = r8rmat_fs_new ( n, a, b );

  r8vec_print ( n, x, "  Solution:" );
//
//  Free memory.
//
  r8rmat_delete ( n, n, a );
  delete [] b;
  delete [] x;

  return;
}
