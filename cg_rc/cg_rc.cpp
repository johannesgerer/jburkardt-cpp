# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "cg_rc.hpp"

//****************************************************************************80

int cg_rc ( int n, double b[], double x[], double r[], double z[], 
  double p[], double q[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    CG_RC is a reverse communication conjugate gradient routine.
//
//  Discussion:
//
//    This routine seeks a solution of the linear system A*x=b
//    where b is a given right hand side vector, A is an n by n
//    symmetric positive definite matrix, and x is an unknown vector
//    to be determined.
//
//    Under the assumptions that the matrix A is large and sparse,
//    the conjugate gradient method may provide a solution when
//    a direct approach would be impractical because of excessive
//    requirements of storage or even of time.
//
//    The conjugate gradient method presented here does not require the 
//    user to store the matrix A in a particular way.  Instead, it only 
//    supposes that the user has a way of calculating
//      y = alpha * A * x + b * y
//    and of solving the preconditioned linear system
//      M * x = b
//    where M is some preconditioning matrix, which might be merely
//    the identity matrix, or a diagonal matrix containing the
//    diagonal entries of A.
//
//    This routine was extracted from the "templates" package.
//    There, it was not intended for direct access by a user;
//    instead, a higher routine called "cg()" was called once by
//    the user.  The cg() routine then made repeated calls to 
//    cgrevcom() before returning the result to the user.
//
//    The reverse communication feature of cgrevcom() makes it, by itself,
//    a very powerful function.  It allows the user to handle issues of
//    storage and implementation that would otherwise have to be
//    mediated in a fixed way by the function argument list.  Therefore,
//    this version of cgrecom() has been extracted from the templates
//    library and documented as a stand-alone procedure.
//
//    The user sets the value of JOB to 1 before the first call,
//    indicating the beginning of the computation, and to the value of
//    2 thereafter, indicating a continuation call.  
//    The output value of JOB is set by cgrevcom(), which
//    will return with an output value of JOB that requests a particular
//    new action from the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
//    Charles Romine, Henk van der Vorst,
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Methods,
//    SIAM, 1994,
//    ISBN: 0898714710,
//    LC: QA297.8.T45.
//
//  Parameters:
//
//    Input, int N, the dimension of the matrix.
//
//    Input, double B[N], the right hand side vector.
//
//    Input/output, double X[N].  On first call, the user 
//    should store an initial guess for the solution in X.  On return with
//    JOB = 4, X contains the latest solution estimate.
//
//    Input/output, double R[N], Z[N], P[N], Q[N],
//    information used by the program during the calculation.  The user
//    does not need to initialize these vectors.  However, specific
//    return values of JOB may require the user to carry out some computation
//    using data in some of these vectors.
//
//    Input/output, int JOB, communicates the task to be done.
//    The user needs to set the input value of JOB to 1, before the first call,
//    and then to 2 for every subsequent call for the given problem.
//    The output value of JOB indicates the requested user action.  
//    * JOB = 1, compute Q = A * P;
//    * JOB = 2: solve M*Z=R, where M is the preconditioning matrix;
//    * JOB = 3: compute R = R - A * X;
//    * JOB = 4: check the residual R for convergence.  
//               If satisfactory, terminate the iteration.
//               If too many iterations were taken, terminate the iteration.
//
{
  double alpha;
  double beta;
  int i;
  static int iter;
  int job_next;
  double pdotq;
  static double rho;
  static double rho_old;
  static int rlbl;
//
//  Initialization.
//  Ask the user to compute the initial residual.
//
  if ( job == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      r[i] = b[i];
    }
    job_next = 3;
    rlbl = 2;
  }
//
//  Begin first conjugate gradient loop.
//  Ask the user for a preconditioner solve.
//
  else if ( rlbl == 2 )
  {
    iter = 1;

    job_next = 2;
    rlbl = 3;
  }
//
//  Compute the direction.
//  Ask the user to compute ALPHA.
//  Save A*P to Q.
//
  else if ( rlbl == 3 )
  {
    rho = 0.0;
    for ( i = 0; i < n; i++ )
    {
      rho = rho + r[i] * z[i];
    }

    if ( 1 < iter )
    {
      beta = rho / rho_old;
      for ( i = 0; i < n; i++ )
      {
        z[i] = z[i] + beta * p[i];
      }
    }

    for ( i = 0; i < n; i++ )
    {
      p[i] = z[i];
    }

    job_next = 1;
    rlbl = 4;
  }
//
//  Compute current solution vector.
//  Ask the user to check the stopping criterion.
//
  else if ( rlbl == 4 )
  {
    pdotq = 0.0;
    for ( i = 0; i < n; i++ )
    {
      pdotq = pdotq + p[i] * q[i];
    }
    alpha = rho / pdotq;
    for ( i = 0; i < n; i++ )
    {
      x[i] = x[i] + alpha * p[i];
    }
    for ( i = 0; i < n; i++ )
    {
      r[i] = r[i] - alpha * q[i];
    }
    job_next = 4;
    rlbl = 5;
  }
//
//  Begin the next step.
//  Ask for a preconditioner solve.
//
  else if ( rlbl == 5 )
  {
    rho_old = rho;
    iter = iter + 1;

    job_next = 2;
    rlbl = 3;
  }

  return job_next;
}
//****************************************************************************80

void r8mat_mv ( int m, int n, double a[], double x[], double ax[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV multiplies a matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as an argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double AX[M], the product A*X.
//
{
  int i;
  int j;

  for ( i = 0; i < m; i++ )
  {
    ax[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      ax[i] = ax[i] + a[i+j*m] * x[j];
    }
  }

  return;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
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
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
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
//****************************************************************************80

double *wathen ( int nx, int ny, int n )

//****************************************************************************80
//
//  Purpose:
//
//    WATHEN returns the WATHEN matrix.
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
//    and sufficient storage in A must have been set aside to hold
//    the matrix.
//
//    A is the consistent mass matrix for a regular NX by NY grid
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
//  Properties:
//
//    A is symmetric positive definite for any positive values of the
//    density RHO(NX,NY), which is here given the value 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2013
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
//    Input, int N, the order of the matrix.
//
//    Output, double WATHEN[N*N], the matrix.
//
{
  double *a;
  static double em[8*8] = {
     6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, 
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, 
     2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, 
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, 
     3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, 
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, 
     2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, 
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 };
  int i;
  int j;
  int kcol;
  int krow;
  int node[8];
  double rho;
  
  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }

  for ( j = 1; j <= ny; j++ )
  {
    for ( i = 1; i <= nx; i++ )
    {
//
//  For the element (I,J), determine the indices of the 8 nodes.
//
      node[0] = 3 * j * nx + 2 * j + 2 * i;
      node[1] = node[0] - 1;
      node[2] = node[0] - 2;
      node[3] = ( 3 * j - 1 ) * nx + 2 * j + i - 2;
      node[4] = ( 3 * j - 3 ) * nx + 2 * j + 2 * i - 4;
      node[5] = node[4] + 1;
      node[6] = node[4] + 2;
      node[7] = node[3] + 1;
//
//  The density RHO can also be set to a random positive value.
//
      for ( krow = 0; krow < 8; krow++ )
      {
        for ( kcol = 0; kcol < 8; kcol++ )
        {
          rho = 1.0;

          if ( node[krow] < 0 || n <= node[krow] ||
               node[kcol] < 0 || n <= node[kcol] )
          {
            cerr << "\n";
            cerr << "WATHEN - Fatal error!\n";
            cerr << "  I = " << i << "  J = " << j << "\n";
            cerr << "  KROW = " << krow << "\n";
            cerr << "  KCOL = " << kcol << "\n";
            cerr << "  NODE[KROW] = " << node[krow] << "\n";
            cerr << "  NODE[KCOL] = " << node[kcol] << "\n";
            exit ( 1 );
          }
          a[node[krow]+node[kcol]*n] = a[node[krow]+node[kcol]*n]
            + 20.0 * rho * em[krow+kcol*8] / 9.0;
        }
      }
    }
  }
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
//  Discussion:
//
//    N = 3*NX*NY + 2*NX + 2*NY + 1,
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2011
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
//    Output, int WATHEN_ORDER, the order of the matrix.
//
{
  int n;

  n = 3 * nx * ny + 2 * nx + 2 * ny + 1;

  return n;
}
