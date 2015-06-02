# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "umfpack.h"

int main ( );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for UMFPACK_SIMPLE.
//
//  Discussion:
//
//    This program uses UMFPACK to solve the 5x5 linear system A*X=B:
//
//        2  3  0  0  0        1.0         8.0
//        3  0  4  0  6        2.0        45.0
//    A = 0 -1 -3  2  0    X = 3.0    B = -3.0
//        0  0  1  0  0        4.0         3.0
//        0  4  2  0  1        5.0        10.0
//
//    The matrix contains 12 nonzero values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2014
//
//  Reference:
//
//    Timothy Davis,
//    UMFPACK User Guide,
//    Version 5.6.2, 25 April 2013
//    http://suitesparse.com
//
{
# define N 5
# define NCC 12

  int Ai[NCC] = { 
    0, 1, 
    0, 2, 4, 
    1, 2, 3, 4, 
    2, 
    1, 4 };
  int Ap[N+1] = { 0, 2, 5, 9, 10, 12 };
  double Ax[NCC] = { 
    2.0,  3.0, 
    3.0, -1.0, 4.0, 
    4.0, -3.0, 1.0, 2.0, 
    2.0, 
    6.0, 1.0 };
  double b[N] = { 8.0, 45.0, -3.0, 3.0, 19.0 };
  int i;
  int n = 5;
  double *null = ( double * ) NULL;
  void *Numeric;
  int status;
  void *Symbolic;
  double x[N];

  timestamp ( );
  cout << "\n";
  cout << "UMFPACK_SIMPLE:\n";
  cout << "  C++ version\n";
  cout << "  Use UMFPACK to solve the sparse linear system A*x=b.\n";
//
//  Carry out the symbolic factorization.
//
  status = umfpack_di_symbolic ( n, n, Ap, Ai, Ax, &Symbolic, null, null );
//
//  Use the symbolic factorization to carry out the numeric factorization.
//
  status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null );
//
//  Free the memory associated with the symbolic factorization.
//
  umfpack_di_free_symbolic ( &Symbolic );
//
//  Solve the linear system.
//
  status = umfpack_di_solve ( UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null );
//
//  Free the memory associated with the numeric factorization.
//
  umfpack_di_free_numeric ( &Numeric );
//
//  Print the solution.
//
  cout << "\n";
  cout << "  Computed solution:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  x[" << i << "] = " << x[i] << "\n";
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "UMFPACK_SIMPLE:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
# undef N
# undef NCC
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
