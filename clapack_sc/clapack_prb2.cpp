# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "clapack.h"

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CLAPACK_PRB2.
//
//    We want to use the LAPACK routine DGESV to solve a linear
//    system of the form A*x=b.  Our 4x4 version of this problem is:
//
//      1.0 -1.0  2.0 -1.0   x(0)     -8.0
//      2.0 -2.0  3.0 -3.0 * x(1)  = -20.0
//      1.0  1.0  1.0  0.0   x(2)     -2.0
//      1.0 -1.0  4.0  3.0   x(3)      4.0
//
//    The solution is x = ( -7, 3, 2, 2 ).
//
//    We want to store the matrix A as a doubly-indexed array.
//
//    In order for this C program to call the LAPACK routine, we have to
//    use the CLAPACK interface.  This requires the following things:
//
//    1) we must "include" clapack.h.
//    
//    2) all integer variables that are passed to the CLAPACK routines must
//       be declared using the "long int" data type.
//
//    3) all double-indexed arrays must be converted to column-major vectors.
//       copy a[i][j] into A[i+j*N], where N is the size of the linear system.
//
//    4) all scalar arguments to the CLAPACK routines must be passed as
//       addresses, that is, their names must be preceded by the ampersand.
//
//    5) the CLAPACK routines must be called in lower case, and with an
//       underscore at the end.
//
//    6) to compile the program, you might use a command like this:
//       gcc -c myprog.c -I/usr/common/clapack
//
//    7) to load the program, you might use a command like this:
//       gcc myprog.o -L/usr/common/clapack -lclapack -lm
//
//    8) to compile and load in one step, use the command
//       gcc myprog.c -I/usr/common/clapack -L/usr/common/clapack -lclapack -lm
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[4][4] = {
    {  1.0, -1.0,  2.0, -1.0 },
    {  2.0, -2.0,  3.0, -3.0 },
    {  1.0,  1.0,  1.0,  0.0 },
    {  1.0, -1.0,  4.0,  3.0 } };
  double *A;
  double B[4] = {
    -8.0, -20.0, -2.0, 4.0 };
  int i;
  long int INFO;
  long int IPIV[4];
  int j;
  long int LDA = 4;
  long int LDB = 4;
  long int N = 4;
  long int NRHS = 1;

  cout << "\n";
  cout << "DGESV_TEST\n";
  cout << "  Demonstrate the use of DGESV to solve a linear system\n";
  cout << "  using double precision real arithmetic.\n";
//
//  Print the coefficient matrix "a" in its original form.
//
  cout << "\n";
  cout << "  Coefficient matrix:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << "  " << setw(14) << a[i][j];
    }
    cout << "\n";
  }
//
//  Print the right hand side B.
//
  cout << "\n";
  cout << "  Right hand side:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(14) << B[i] << "\n";
  }
//
//  Convert the matrix "a" to a FORTRAN vector "A".
//
  A = new double[N*N];

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      A[i+j*N] = a[i][j];
    }
  }
//
//  Call DGESV to compute the solution.
//  We used capitals for the names of all the arguments.
//  We must use addresses (&prefix) for all scalar arguments.
//  We must use a lower case version of the function name,
//  followed by an underscore.
//
  dgesv_ ( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );
//
//  If INFO is not zero, the algorithm failed.
//
  if ( ( int ) INFO != 0 )
  {
    cout << "\n";
    cout << "  Return value of error flag INFO = " << ( int ) INFO << "\n";
    return 1;
  }
//
//  Print the solution.
//
  cout << "\n";
  cout << "  Computed solution:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(14) << B[i] << "\n";
  }

  return 0;
}

