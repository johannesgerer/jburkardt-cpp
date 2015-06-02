# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "stochastic_heat2d.hpp"

void boundary ( int nx, int ny, double x[], double y[], 
  int n, double a[], double u[] );

//****************************************************************************80

double diffusivity_2d_bnt ( double dc0, double omega[], double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFUSIVITY_2D_BNT evaluates a 2D stochastic diffusivity function.
//
//  Discussion:
//
//    The 2D diffusion equation has the form
//
//      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
//
//    where DC(X,Y) is a function called the diffusivity.
//
//    In the stochastic version of the problem, the diffusivity function
//    includes the influence of stochastic parameters:
//
//      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
//
//    In this function, the domain is the rectangle [-1.5,0]x[-0.4,0.8].
//
//    The four stochastic parameters OMEGA(1:4) are assumed to be independent
//    identically distributed random variables with mean value zero and 
//    variance 1.  The distribution is typically taken to be Gaussian or
//    uniform.
//
//    A collocation approach to this problem would then use the roots of
//    Hermite or Legendre polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ivo Babuska, Fabio Nobile, Raul Tempone,
//    A stochastic collocation method for elliptic partial differential equations
//    with random input data,
//    SIAM Journal on Numerical Analysis,
//    Volume 45, Number 3, 2007, pages 1005-1034.
//
//  Parameters:
//
//    Input, double DC0, the constant term in the expansion of the 
//    diffusion coefficient.  Take DC0 = 10.
//
//    Input, double OMEGA[4], the stochastic parameters.
//
//    Input, double X, Y, the points where the diffusion 
//    coefficient is to be evaluated.
//
//    Output, double DIFFUSIVITY_2D_BNT, the value of the diffusion
//    coefficient at (X,Y).
//
{
  double arg;
  double dc;
  double pi = 3.141592653589793;

  arg = omega[0] * cos ( pi * x ) 
      + omega[1] * sin ( pi * x ) 
      + omega[2] * cos ( pi * y ) 
      + omega[3] * sin ( pi * y );

  arg = exp ( - 0.125 ) * arg;

  dc = dc0 + exp ( arg );

  return dc;
}
//****************************************************************************80

void interior ( double omega[], int nx, int ny, double x[], double y[], 
  double f ( double x, double y ), int n, double a[], double rhs[] )

//****************************************************************************80
//
//  Purpose:
//
//    INTERIOR sets up the matrix and right hand side at interior nodes.
//
//  Discussion:
//
//    Nodes are assigned a single index K, which increases as:
//
//    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
//           ....         ....  ...    .....
//           NX+1         NX+2  ...   2 * NX
//              1            2  ...       NX
//
//    Therefore, the neighbors of an interior node numbered C are
//
//             C+NY
//              |
//      C-1 --- C --- C+1
//              |
//             C-NY
//
//    If we number rows from bottom I = 1 to top I = NY
//    and columns from left J = 1 to right J = NX, then the relationship
//    between the single index K and the row and column indices I and J is:
//      K = ( I - 1 ) * NX + J
//    and
//      J = 1 + mod ( K - 1, NX )
//      I = 1 + ( K - J ) / NX
//      
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double OMEGA[4], the stochastic coefficients.
//
//    Input, int NX, NY, the number of grid points in X and Y.
//
//    Input, double X[NX], Y[NY], the coordinates of grid lines.
//
//    Input, double function F ( double X, double Y ), evaluates the heat 
//    source term.
//
//    Input, int N, the number of nodes.
//
//    Output, double A[N*N], the system matrix, with the entries for 
//    the interior nodes filled in.
//
//    Output, double RHS[N], the system right hand side, with the 
//    entries for the interior nodes filled in.
//
{
  double dc0;
  double dce;
  double dcn;
  double dcs;
  double dcw;
  double dx;
  double dy;
  int ic;
  int in;
  int is;
  int jc;
  int je;
  int jw;
  int kc;
  int ke;
  int kn;
  int ks;
  int kw;
  double xce;
  double xcw;
  double ycn;
  double ycs;

  dc0 = 1.0;
//
//  For now, assume X and Y are equally spaced.
//
  dx = x[1] - x[0];
  dy = y[1] - y[0];

  for ( ic = 1; ic < ny - 1; ic++ )
  {
    for ( jc = 1; jc < nx - 1; jc++ )
    {
      in = ic + 1;
      is = ic - 1;
      je = jc + 1;
      jw = jc - 1;

      kc = ic * nx + jc;
      ke = kc + 1;
      kw = kc - 1;
      kn = kc + nx;
      ks = kc - nx;

      xce = 0.5 * ( x[jc] + x[je] );
      dce = diffusivity_2d_bnt ( dc0, omega, xce,   y[ic] );
      xcw = 0.5 * ( x[jc] + x[jw] );
      dcw = diffusivity_2d_bnt ( dc0, omega, xcw,   y[ic] );
      ycn = 0.5 * ( y[ic] + y[in] );
      dcn = diffusivity_2d_bnt ( dc0, omega, x[jc], ycn );
      ycs = 0.5 * ( y[ic] + y[is] );
      dcs = diffusivity_2d_bnt ( dc0, omega, x[jc], ycs );

      a[kc+kc*n] = ( dce + dcw ) / dx / dx + ( dcn + dcs ) / dy / dy;
      a[kc+ke*n] = - dce         / dx / dx;
      a[kc+kw*n] =       - dcw   / dx / dx;
      a[kc+kn*n] =                           - dcn         / dy / dy;
      a[kc+ks*n] =                                 - dcs   / dy / dy;

      rhs[kc] = f ( x[jc], y[ic] );
    }
  }

  return;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
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
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8mat_fs ( int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_FS factors and solves a system with one right hand side.
//
//  Discussion:
//
//    This routine differs from R8MAT_FSS in two ways:
//    * only one right hand side is allowed;
//    * the input matrix A is not modified.
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[N*N], the coefficient matrix of the linear system.
//
//    Input/output, double X[N], on input, the right hand side of the
//    linear system.  On output, the solution of the linear system.
//
{
  double *a2;
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  a2 = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }

  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = fabs ( a2[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < fabs ( a2[i-1+(jcol-1)*n] ) )
      {
        piv = fabs ( a2[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cerr << "\n";
      cerr << "R8MAT_FS - Fatal error!\n";
      cerr << "  Zero pivot on step " << jcol << "\n";
      exit ( 1 );
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                  = a2[jcol-1+(j-1)*n];
        a2[jcol-1+(j-1)*n] = a2[ipiv-1+(j-1)*n];
        a2[ipiv-1+(j-1)*n] = t;
      }
      t         = x[jcol-1];
      x[jcol-1] = x[ipiv-1];
      x[ipiv-1] = t;
    }
//
//  Scale the pivot row.
//
    t = a2[jcol-1+(jcol-1)*n];
    a2[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a2[jcol-1+(j-1)*n] = a2[jcol-1+(j-1)*n] / t;
    }
    x[jcol-1] = x[jcol-1] / t;
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a2[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a2[i-1+(jcol-1)*n];
        a2[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a2[i-1+(j-1)*n] = a2[i-1+(j-1)*n] + t * a2[jcol-1+(j-1)*n];
        }
        x[i-1] = x[i-1] + t * x[jcol-1];
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      x[i-1] = x[i-1] - a2[i-1+(jcol-1)*n] * x[jcol-1];
    }
  }

  delete [] a2;

  return;
}
//****************************************************************************80

double r8mat_max ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MAX returns the maximum entry of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Output, double R8MAT_MAX, the maximum entry of A.
//
{
  int i;
  int j;
  double value;

  value = a[0+0*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < a[i+j*m] )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
//****************************************************************************80

double r8mat_mean ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MEAN returns the mean of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Output, double R8MAT_MEAN, the mean of A.
//
{
  int i;
  int j;
  double value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + a[i+j*m];
    }
  }
  value = value / ( double ) ( m * n );

  return value;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
//
//    In other words, the interval is divided into N-1 even subintervals,
//    and the endpoints of intervals are used as the points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

void r8vec_mesh_2d ( int nx, int ny, double xvec[], double yvec[], 
  double xmat[], double ymat[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MESH_2D creates a 2D mesh from X and Y vectors.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    NX = 2
//    XVEC = ( 1, 2, 3 )
//    NY = 3
//    YVEC = ( 4, 5 )
//
//    XMAT = (
//      1, 2, 3
//      1, 2, 3 )
//
//    YMAT = (
//      4, 4, 4
//      5, 5, 5 ) 
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 July 2013
//
//  Parameters:
//
//    Input, int NX, NY, the number of X and Y values.
//
//    Input, double XVEC[NX], YVEC[NY], the X and Y coordinate
//    values.
//
//    Output, double XMAT[NX*NY], YMAT[NX*NY], the coordinate
//    values of points on an NX by NY mesh.
//
{
  int i;
  int j;

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      xmat[i+j*nx] = xvec[i];
    }
  }

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      ymat[i+j*nx] = yvec[j];
    }
  }

 return;
}
//****************************************************************************80

double *r8vec_normal_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, double R[N+1], is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.
//
{
  int i;
  int m;
  const double pi = 3.141592653589793;
  double *r;
  double *x;
  int x_hi;
  int x_lo;

  x = new double[n];
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  If we need just one new value, do that here to avoid null arrays.
//
  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );

    delete [] r;
  }

  return x;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
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

double *stochastic_heat2d ( double omega[], int nx, int ny, double x[], 
  double y[], double f ( double x, double y ) )

//****************************************************************************80
//
//  Purpose:
//
//    STOCHASTIC_HEAT2D solves the steady 2D heat equation.
//
//  Discussion:
//
//    Nodes are assigned a singled index K, which increases as:
//
//    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
//           ....         ....  ...    .....
//           NX+1         NX+2  ...   2 * NX
//              1            2  ...       NX
//
//    Therefore, the neighbors of an interior node numbered C are
//
//             C+NY
//              |
//      C-1 --- C --- C+1
//              |
//             C-NY
//
//    Nodes on the lower boundary satisfy:
//      1 <= K <= NX
//    Nodes on the upper boundary satisfy:
//      (NY-1)*NX+1 <= K <= NY * NX
//    Nodes on the left boundary satisfy:
//      mod ( K, NX ) = 1
//    Nodes on the right boundary satisfy:
//      mod ( K, NX ) = 0
//
//    If we number rows from bottom I = 1 to top I = NY
//    and columns from left J = 1 to right J = NX, we have
//      K = ( I - 1 ) * NX + J
//    and
//      J = 1 + mod ( K - 1, NX )
//      I = 1 + ( K - J ) / NX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double OMEGA[4], the stochastic coefficients.
//
//    Input, int NX, NY, the number of grid points in X and Y.
//
//    Input, double X[NX], Y[NY], the coordinates of grid lines.
//
//    Input, double F ( double X, double Y ), evaluates the heat 
//    source term.
//
//    Output, double U[NX*NY], the approximation to the solution at 
//    the grid points.
//
{
  double *a;
  int n;
  double *u;
//
//  Set the total number of unknowns.
//
  n = nx * ny;
//
//  Set up the matrix and right hand side.
//
  a = new double[n*n];
  u = new double[n];
//
//  Define the matrix at interior points.
//
  interior ( omega, nx, ny, x, y, f, n, a, u );
//
//  Handle boundary conditions.
//
  boundary ( nx, ny, x, y, n, a, u );
//
//  Solve the linear system.
//
  r8mat_fs ( n, a, u );
//
//  Free memory.
//
  delete [] a;

  return u;
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
