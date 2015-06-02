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
//    MAIN is the main program for CLAPACK_PRB3.
//
//    We want to use the LAPACK routines DGBTRF and DGBTRS to solve a banded 
//    linear system of the form A*x=b.  Our 4x4 version of this problem is:
//
//      2.0 -1.0  0.0  0.0   x(0)      0.0
//     -1.0  2.0 -1.0  0.0 * x(1)  =   0.0
//      0.0 -1.0  2.0 -1.0   x(2)      0.0
//      0.0  0.0 -1.0  2.0   x(3)      5.0
//
//    The solution is x = ( 1, 2, 3, 4 ).
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
//       The C copy of the matrix can be initialized as a collection of rows.
//       This information must be converted to a FORTRAN vector, and the
//       indexing is pretty horrible, so look at the example below.
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
//    09 January 2014
//
//  Author:
//
//    John Burkardt
//
{
//
//  The 4x4 matrix is stored as 4 rows, containing the 3 nonzero elements
//  of each band.  The first and last rows have one extra dummy element.
//
  double a[4][3] = {
    {  0.0,  2.0, -1.0 },
    { -1.0,  2.0, -1.0 },
    { -1.0,  2.0, -1.0 },
    { -1.0,  2.0,  0.0 } };
  double *AB;
  double B[4] = {
    0.0, 0.0, 0.0, 5.0 };
  int i;
  long int INFO;
  long int IPIV[4];
  int j;
  int k;
  long int KL = 1;
  long int KU = 1;
  long int LDAB = 4;
  long int LDB = 4;
  long int M = 4;
  long int N = 4;
  long int NRHS = 1;
  char TRANS = 'N';

  cout << "\n";
  cout << "CLAPACK_PRB3\n";
  cout << "  Demonstrate the use of DGBTRF to factor a banded matrix\n";
  cout << "  and DGBTRS to solve an associated linear system\n";
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
      if ( j < i - KL || i + KU < j  )
      {
        cout << "  " << setw(14) << 0.0;
      }
      else
      {
        cout << "  " << setw(14) << a[i][j+KL-i];
      }
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
  AB = new double[ ( 2 * KL + KU + 1 ) * N];

  for ( k = 0; k < ( 2 * KL + KU + 1 ) * N; k++ )
  {
    AB[k] = 0.0;
  }

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i - KL <= j && j <= i + KU )
      {
        AB[KL+KU+i-j+j*N] = a[i][j+KL-i];
      }
    }
  }
//
//  Call DGBTRF to factor the matrix.
//
  dgbtrf_ ( &M, &N, &KL, &KU, AB, &LDAB, IPIV, &INFO );
//
//  If INFO is not zero, the factor algorithm failed.
//
  if ( ( int ) INFO != 0 )
  {
    cout << "\n";
    cout << "  DGBTRF returned error flag INFO = " << ( int ) INFO << "\n";
    return 1;
  }
/*
  Call DGBTRS to compute the solution.
*/
  dgbtrs_ ( &TRANS, &N, &KL, &KU, &NRHS, AB, &LDAB, IPIV, B, &LDB, &INFO );
/*
  If INFO is not zero, the solve algorithm failed.
*/
  if ( ( int ) INFO != 0 )
  {
    cout << "\n";
    cout << "  DGBTRS returned error flag INFO = " << ( int ) INFO << "\n";
    return 1;
  }
/*
  Print the solution.
*/
  cout << "\n";
  cout << "  Computed solution:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << "  " << setw(14) << B[i] << "\n";
  }
//
//  Free memory.
//
  delete [] AB;

  return 0;
}

