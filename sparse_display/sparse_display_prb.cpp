# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "sparse_display.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
double *wathen_ge ( int nx, int ny, int n, int &seed );
int wathen_order ( int nx, int ny );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_DISPLAY_PRB.
//
//  Discussion:
//
//    SPARSE_DISPLAY_PRB tests the SPARSE_DISPLAY library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPARSE_DISPLAY_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPARSE_DISPLAY library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPARSE_DISPLAY_PRB\n";
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
//    TEST01 tests SPY_GE for a general storage matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  string header = "wathen_ge";
  int n;
  int nx;
  int ny;
  int seed;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  SPY_GE generates a sparsity plot for a matrix stored\n";
  cout << "  in general (GE) format.\n";

  nx = 5;
  ny = 5;
  n = wathen_order ( nx, ny );
  seed = 123456789;
  a = wathen_ge ( nx, ny, n, seed );

  spy_ge ( n, n, a, header );

  delete [] a;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests SPY_FILE in cases where the indices are in a file.
//
//  Discussion:
//
//    The files used in this example actually use negative column indices
//    because they were output by DEAL.II and intended to be passed directly
//    into GNUPLOT without any clever commands.
//
//    So my own "SPY_FILE" is currently set up to deal exactly with such
//    files, and hence, given sensible data will actually show a spy plot
//    that is transposed - still readable and all, but wrong way round.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST02\n";
  cout << "  SPY_FILE generates a sparsity plot for a matrix for\n";
  cout << "  which there exists a file containing all the pairs\n";
  cout << "  (I,J) corresponding to nonzero entries.\n";

  spy_file ( "before", "before_data.txt" );
  spy_file ( "after", "after_data.txt" );

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests SPY_GE for a general storage matrix.
//
//  Discussion:
//
//    It is not clear to me whether the plot being created is correctly
//    oriented.  I might be seeing the transpose of the matrix.
//    One way to check is to set up a moderate sized, highly asymmetric matrix.
//    In this case, I will create a certain upper triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  string header = "20x30";
  int i;
  int j;
  int m;
  int n;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  SPY_GE generates a sparsity plot for a matrix stored\n";
  cout << "  in general (GE) format.\n";
  cout << "  Just to orient ourselves, generate an upper triangular matrix.\n";

  m = 20;
  n = 30;
  a = new double[m*n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( j < i || ( ( i - j ) % ( i + 1 ) ) != 0 )
      {
        a[i+j*m] = 0.0;
      }
      else
      {
        a[i+j*m] = 1.0;
      }
    }
  }

  spy_ge ( m, n, a, header );

  delete [] a;

  return;
}
//****************************************************************************80

double *r8mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
//    in column-major order.
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }
      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

double *wathen_ge ( int nx, int ny, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    WATHEN_GE returns the Wathen matrix as a general storage (GE) matrix.
//
//  Discussion:
//
//    The Wathen matrix is a finite element matrix which is sparse.
//
//    The entries of the matrix depend in part on a physical quantity
//    related to density.  That density is here assigned random values between
//    0 and 100.
//
//    The matrix order N is determined by the input quantities NX and NY,
//    which would usually be the number of elements in the X and Y directions.
//    The value of N is
//
//      N = 3*NX*NY + 2*NX + 2*NY + 1,
//
//    The matrix is the consistent mass matrix for a regular NX by NY grid
//    of 8 node serendipity elements.
//
//    Here is an illustration for NX = 3, NY = 2:
//
//     23-24-25-26-27-28-29
//      |     |     |     |
//     19    20    21    22
//      |     |     |     |
//     12-13-14-15-16-17-18
//      |     |     |     |
//      8     9    10    11
//      |     |     |     |
//      1--2--3--4--5--6--7
//
//    For this example, the total number of nodes is, as expected,
//
//      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
//
//    The matrix is symmetric positive definite for any positive values of the
//    density RHO(X,Y).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Nicholas Higham,
//    Algorithm 694: A Collection of Test Matrices in MATLAB,
//    ACM Transactions on Mathematical Software,
//    Volume 17, Number 3, September 1991, pages 289-305.
//
//    Andrew Wathen,
//    Realistic eigenvalue bounds for the Galerkin mass matrix,
//    IMA Journal of Numerical Analysis,
//    Volume 7, Number 4, October 1987, pages 449-457.
//
//  Parameters:
//
//    Input, int NX, NY, values which determine the size 
//    of the matrix.
//
//    Input, int N, the number of rows and columns.
//
//    Input/output, int &SEED, the random number seed.
//
//    Output, double WATHEN_GE[N*N], the matrix.
//
{
  double *a;
  const double em[8*8] = {
     6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, 
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, 
     2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, 
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, 
     3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, 
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, 
     2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, 
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 };
  int i;
  int ii;
  int j;
  int jj;
  int kcol;
  int krow;
  int node[8];
  double *rho;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }

  rho = r8mat_uniform_01_new ( nx, ny, seed );

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      rho[i+j*nx] = 100.0 * rho[i+j*nx];
    }
  }

  for ( j = 0; j < nx; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      node[0] = 3 * ( j + 1 ) * nx + 2 * ( j + 1 ) + 2 * ( i + 1 );
      node[1] = node[0] - 1;
      node[2] = node[0] - 2;
      node[3] = ( 3 * ( j + 1 ) - 1 ) * nx + 2 * ( j + 1 ) + ( i + 1 ) - 2;
      node[4] = ( 3 * ( j + 1 ) - 3 ) * nx + 2 * ( j + 1 ) + 2 * ( i + 1 ) - 4;
      node[5] = node[4] + 1;
      node[6] = node[4] + 2;
      node[7] = node[3] + 1;

      for ( krow = 0; krow < 8; krow++ )
      {
        ii = node[krow];
        for ( kcol = 0; kcol < 8; kcol++ )
        {
          jj = node[kcol];
          a[ii+jj*n] = a[ii+jj*n] + rho[i+j*nx] * em[krow+kcol*8];
        }
      }
    }
  }

  delete [] rho;

  return a;
}
//****************************************************************************80

int wathen_order ( int nx, int ny )

//****************************************************************************80
//
//  Purpose:
//
//    WATHEN_ORDER returns the order of the WATHEN matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Nicholas Higham,
//    Algorithm 694: A Collection of Test Matrices in MATLAB,
//    ACM Transactions on Mathematical Software,
//    Volume 17, Number 3, September 1991, pages 289-305.
//
//    Andrew Wathen,
//    Realistic eigenvalue bounds for the Galerkin mass matrix,
//    IMA Journal of Numerical Analysis,
//    Volume 7, 1987, pages 449-457.
//
//  Parameters:
//
//    Input, int NX, NY, values which determine the size of A.
//
//    Output, int WATHEN_ORDER, the order of the matrix,
//    as determined by NX and NY.
//
{
  int n;

  n = 3 * nx * ny + 2 * nx + 2 * ny + 1;

  return n;
}
