# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "stochastic_diffusion.hpp"

//****************************************************************************80

double *diffusivity_1d_xk ( double dc0, int m, double omega[], int n,
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFUSIVITY_1D_XK evaluates a 1D stochastic diffusivity function.
//
//  Discussion:
//
//    The 1D diffusion equation has the form
//
//      - d/dx ( DC(X) Del U(X) ) = F(X)
//
//    where DC(X) is a function called the diffusivity.
//
//    In the stochastic version of the problem, the diffusivity function
//    includes the influence of stochastic parameters:
//
//      - d/dx ( DC(X;OMEGA) d/dx U(X) ) = F(X).
//
//    In this function, the domain is assumed to be the unit interval [0.1].
//
//
//    For DC0 = 1 and F(X) = 0, with boundary conditions U(0:OMEGA) = 0,
//    U(1;OMEGA) = 1, the exact solution is
//
//    If OMEGA ~= 0:
//
//      U(X;OMEGA) = log ( 1 + OMEGA * X ) / log ( 1 + OMEGA )
//
//    If OMEGA = 0:
//
//      U(X;OMEGA) = X
//
//    In the numerical experiments described in the paper, OMEGA was taken
//    to be a random variable with a Beta, or Uniform, or Gaussian or 
//    Poisson or Binomial distribution.
//
//    For the Gaussian and Poisson distributions, the positivity requirement 
//    could not be guaranteed, and the experiments were simply made with a 
//    "small" variance of 0.1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu, George Karniadakis,
//    Modeling uncertainty in steady state diffusion problems via
//    generalized polynomial chaos,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 191, 2002, pages 4927-4948.
//
//  Parameters:
//
//    Input, double DC0, the constant term in the expansion of the 
//    diffusion coefficient.
//
//    Input, int M, the number of stochastic parameters.
//
//    Input, double OMEGA[M], the stochastic parameters.  
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], the point where the diffusion coefficient 
//    is to be evaluated.
//
//    Output, double DIFFUSIVITY_1D_XK[N], the value of the diffusion coefficient 
//    at X.
//
{
  double *dc;
  int j;
  int k;
  double pi = 3.141592653589793;
  double w;

  k = 0;
  w = 1.0;

  dc = new double[n];

  for ( j = 0; j < n; j++ )
  {
    dc[j] = 0.0;
  }

  while ( k < m )
  {
    if ( k < m )
    {
      k = k + 1;
      for ( j = 0; j < n; j++ )
      {
        dc[j] = dc[j] + omega[k-1] * sin ( w * pi * x[j] );
      }
    }

    if ( k < m )
    {
      k = k + 1;
      for ( j = 0; j < n; j++ )
      {
        dc[j] = dc[j] + omega[k-1] * cos ( w * pi * x[j] );    
      }
    }

    w = w + 1.0;

  }

  for ( j = 0; j < n; j++ )
  {
    dc[j] = exp ( - 0.125 ) * dc[j];
  }

  for ( j = 0; j < n; j++ )
  {
    dc[j] = dc0 + exp ( dc[j] );
  }

  return dc;
}
//****************************************************************************80

double *diffusivity_2d_bnt ( double dc0, double omega[], int n, double x[], 
  double y[] )

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
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the points where the diffusion 
//    coefficient is to be evaluated.
//
//    Output, double DIFFUSIVITY_2D_BNT[N], the value of the diffusion
//    coefficient at (X,Y).
//
{
  double *arg;
  double *dc;
  int j;
  double pi = 3.141592653589793;

  arg = new double[n];

  for ( j = 0; j < n; j++ )
  {
    arg[j] = omega[0] * cos ( pi * x[j] ) 
           + omega[1] * sin ( pi * x[j] ) 
           + omega[2] * cos ( pi * y[j] ) 
           + omega[3] * sin ( pi * y[j] );
  }

  for ( j = 0; j < n; j++ )
  {
    arg[j] = exp ( - 0.125 ) * arg[j];
  }

  dc = new double[n];
  for ( j = 0; j < n; j++ )
  {
    dc[j] = dc0 + exp ( arg[j] );
  }

  delete [] arg;

  return dc;
}
//****************************************************************************80

double *diffusivity_2d_elman ( double a, double cl, double dc0, int m_1d, 
  double omega[], int n1, int n2, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFUSIVITY_2D_ELMAN evaluates a 2D stochastic diffusivity function.
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
//    In this function, the domain is assumed to be the square [-A,+A]x[-A,+A].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Howard Elman, Darran Furnaval,
//    Solving the stochastic steady-state diffusion problem using multigrid,
//    IMA Journal on Numerical Analysis,
//    Volume 27, Number 4, 2007, pages 675-688.
//
//    Roger Ghanem, Pol Spanos,
//    Stochastic Finite Elements: A Spectral Approach,
//    Revised Edition,
//    Dover, 2003,
//    ISBN: 0486428184,
//    LC: TA347.F5.G56.
//
//  Parameters:
//
//    Input, double A, the "radius" of the square region.  The region
//    is assumed to be [-A,+A]x[-A,+A].
//    0 < A.
//
//    Input, double CL, the correlation length.
//    0 < CL.
//
//    Input, double DC0, the constant term in the expansion of the 
//    diffusion coefficient.  Take DC0 = 10.
//
//    Input, int M_1D, the first and second dimensions of the
//    stochastic parameter array.
//
//    Input, double OMEGA[M_1D*M_1D], the stochastic parameters.
//
//    Input, int N1, N2, the dimensions of the X and Y arrays.
//
//    Input, double X[N1*N2], Y[N1*N2], the points where the diffusion 
//    coefficient is to be evaluated.
//
//    Output, double DIFFUSIVITY_2D_ELMAN[N1*N2], the value of the diffusion 
//    coefficient at X.
//
{
  double *c_1dx;
  double *c_1dy;
  double *dc;
  int i;
  int i1;
  int i2;
  int j;
  int k;
  double *lambda_1d;
  int m;
  double *theta_1d;

  m = m_1d * m_1d;
//
//  Compute THETA.
//
  theta_1d = theta_solve ( a, cl, m_1d );
//
//  Compute LAMBDA_1D.
//
  lambda_1d = new double[m_1d];

  for ( i = 0; i < m_1d; i++ )
  {
    lambda_1d[i] = 2.0 * cl / ( 1.0 + cl * cl * theta_1d[i] * theta_1d[i] );
  }
//
//  Compute C_1DX(1:M1D) and C_1DY(1:M1D) at (X,Y).
//
  c_1dx = new double[m_1d*n1*n2];
  c_1dy = new double[m_1d*n1*n2];

  for ( k = 0; k < n2; k++ )
  {
    for ( j = 0; j < n1; j++ )
    {
      for ( i = 0; i < m_1d; i++ )
      {
        c_1dx[i+j*m_1d+k*m_1d*n1] = 0.0;
        c_1dy[i+j*m_1d+k*m_1d*n1] = 0.0;
      }
    }
  }

  i = 0;

  for ( ; ; )
  {
    if ( m_1d <= i )
    {
      break;
    }

    for ( k = 0; k < n2; k++ )
    {
      for ( j = 0; j < n1; j++ )
      {
        c_1dx[i+j*m_1d+k*m_1d*n1] = cos ( theta_1d[i] * a * x[j+k*n1] ) 
          / sqrt ( a + sin ( 2.0 * theta_1d[i] * a ) 
          / ( 2.0 * theta_1d[i] ) );

        c_1dy[i+j*m_1d+k*m_1d*n1] = cos ( theta_1d[i] * a * y[j+k*n1] ) 
          / sqrt ( a + sin ( 2.0 * theta_1d[i] * a ) 
          / ( 2.0 * theta_1d[i] ) );
      }
    }

    i = i + 1;

    if ( m_1d <= i )
    {
      break;
    }

    for ( k = 0; k < n2; k++ )
    {
      for ( j = 0; j < n1; j++ )
      {
        c_1dx[i+j*m_1d+k*m_1d*n1] = sin ( theta_1d[i] * a * x[j+k*n1] ) 
          / sqrt ( a - sin ( 2.0 * theta_1d[i] * a ) 
          / ( 2.0 * theta_1d[i] ) );

        c_1dy[i+j*m_1d+k*m_1d*n1] = sin ( theta_1d[i] * a * y[j+k*n1] ) 
          / sqrt ( a - sin ( 2.0 * theta_1d[i] * a ) 
          / ( 2.0 * theta_1d[i] ) );
      }
    }

    i = i + 1;
  }
//
//  Evaluate the diffusion coefficient DC at (X,Y).
//
  dc = new double[n1*n2];

  for ( k = 0; k < n2; k++ )
  {
    for ( j = 0; j < n1; j++ )
    {
      dc[j+k*n1] = dc0;
      for ( i2 = 0; i2 < m_1d; i2++ )
      {
        for ( i1 = 0; i1 < m_1d; i1++ )
        {
          dc[j+k*n1] = dc[j+k*n1] + sqrt ( lambda_1d[i1] * lambda_1d[i2] ) 
            * c_1dx[i1+j*m_1d+k*m_1d*n1] * c_1dy[i2+j*m_1d+k*m_1d*n1] 
            * omega[i1+i2*m_1d];
        }
      }
    }
  }

  delete [] c_1dx;
  delete [] c_1dy;
  delete [] lambda_1d;
  delete [] theta_1d;

  return dc;
}
//****************************************************************************80

double *diffusivity_2d_ntw ( double cl, double dc0, int m, double omega[], 
  int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIFFUSIVITY_2D_NTW evaluates a 2D stochastic diffusivity function.
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
//    In this function, the domain is the rectangle [0,D]x[0,D] where D = 1.
//
//    Note that in this problem the diffusivity has a one-dimensional
//    spatial dependence on X, but not on Y
//
//    The random variables OMEGA are independent, have zero mean and unit
//    variance, and are uniformly distributed in [-sqrt(3),+sqrt(3)].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Xiang Ma, Nicholas Zabaras,
//    An adaptive hierarchical sparse grid collocation algorithm for the solution
//    of stochastic differential equations,
//    Journal of Computational Physics,
//    Volume 228, pages 3084-3113, 2009.
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, double CL, the desired physical correlation length for 
//    the coefficient.
//
//    Input, double DC0, the constant term in the expansion of the 
//    diffusion coefficient.  Take DC0 = 0.5.
//
//    Input, int M, the number of terms in the expansion.
//
//    Input, double OMEGA[M], the stochastic parameters.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the points where the diffusion 
//    coefficient is to be evaluated.
//
//    Output, double DIFFUSIVITY_2D_NTW[N], the value of the diffusion coefficient
//    at (X,Y).
//
{
  double d;
  double *dc;
  double *dc_arg;
  int i;
  double ihalf_r8;
  int j;
  double l;
  double lj;
  double lp;
  double *phi;
  double pi = 3.141592653589793;
  double zeta;
  double zeta_arg;

  d = 1.0;
  lp = r8_max ( d, 2.0 * cl );
  l = cl / lp;

  dc_arg = new double[n];

  for ( j = 0; j < n; j++ )
  {
    dc_arg[j] = 1.0 + omega[0] * sqrt ( sqrt ( pi ) * l / 2.0 );
  }

  dc = new double[n];
  phi = new double[n];

  for ( i = 2; i <= m; i++ )
  {
    ihalf_r8 = ( double ) ( i / 2 );
    zeta_arg = - pow ( ihalf_r8 * pi * l, 2 ) / 8.0;
    zeta = sqrt ( sqrt ( pi ) * l ) * exp ( zeta_arg );

    if ( ( i % 2 ) == 0 )
    {
      for ( j = 0; j < n; j++ )
      {
        phi[j] = sin ( ihalf_r8 * pi * x[j] / lp );
      }
    }
    else
    {
      for ( j = 0; j < n; j++ )
      {
        phi[j] = cos ( ihalf_r8 * pi * x[j] / lp );
      }
    }

    for ( j = 0; j < n; j++ )
    {
      dc_arg[j] = dc_arg[j] + zeta * phi[j] * omega[i-1];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    dc[j] = dc0 + exp ( dc_arg[j] );
  }

  delete [] dc_arg;
  delete [] phi;

  return dc;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
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

double r8vec_max ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
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
//    22 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
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

void r8vec_uniform_01 ( int n, int &seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
//    Output, double R[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

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

double *theta_solve ( double a, double cl, int m )

//****************************************************************************80
//
//  Purpose:
//
//    THETA_SOLVE solves a pair of transcendental equations.
//
//  Discussion:
//
//    The vector THETA returned by this function is needed in order to define
//    the terms in a Karhunen-Loeve expansion of a diffusion coefficient.
//
//    The two equations are:
//
//      1/CL - THETA * TAN ( A * THETA ) = 0
//      THETA - 1/CL * TAN ( A * THETA ) = 0
//
//    A and CL are taken to be positive.  Over each open interval 
//
//      ( n - 1/2 pi, n + 1/2 pi ) / A, for N = 0, 1, ...
//
//    the function TAN ( A * THETA ) monotonically rises from -oo to +00; 
//    therefore, it can be shown that there is one root of each equation 
//    in every interval of this form.  Moreover, because of the positivity
//    of A and CL, we can restrict our search to the interval 
//
//      [ n pi, n + 1/2 pi ) / A, for N = 0, 1, ...
//
//    This function computes K such roots, starting in the first interval,
//    finding those two roots, moving to the next interval, and so on, until
//    the requested number of roots have been found.  Odd index roots will
//    correspond to the first equation, and even index roots to the second.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Howard Elman, Darran Furnival,
//    Solving the Stochastic Steady-State Diffusion Problem Using Multigrid,
//    University of Maryland Department of Computer Science,
//    Technical Report TR-4786.
//
//  Parameters:
//
//    Input, double A, the "radius" of the domain, D = (-A,A)x(-A,A).
//    0 < A.
//
//    Input, double CL, the correlation length.
//    0 < CL.
//
//    Input, int M, the number of values to compute.
//
//    Output, double THETA_SOLVE[M], the values of Theta.
//
{
  double bmatol;
  double eps;
  double fa;
  double fb;
  double fc;
  double ftol;
  int k;
  double pi = 3.141592653589793;
  double *theta;
  double xa;
  double xa_init;
  double xb;
  double xb_init;
  double xc;

  theta = new double[m];
  for ( k = 0; k < m; k++ )
  {
    theta[k] = 0.0;
  }
//
//  [ XA_INIT, XB_INIT] = [ n * pi, n+1/2 pi ] / a, n = 0, 1, 2, ...
//
  xa_init = 0.0;
  xb_init = ( pi / 2.0 ) / a;
  eps = r8_epsilon ( );

  k = 0;
  for ( ; ; )
  {
//
//  Seek root of equation 1 in interval.
//
    if ( m <= k )
    {
      break;
    }

    k = k + 1;
    xa = xa_init;
    fa = 1.0 / cl - xa * tan ( a * xa );
    ftol = eps * ( fabs ( fa ) + 1.0 );
    xb = xb_init;
    fb = - fa;
    fc = fa;
    bmatol = 100.0 * eps * ( fabs ( xa ) + fabs ( xb ) );

    while ( bmatol < xb - xa )
    {
      xc = ( xa + xb ) / 2.0;
      fc = 1.0 / cl - xc * tan ( a * xc );

      if ( fabs ( fc ) <= ftol )
      {
        break;
      }
      else if ( 0.0 < fc )
      {
        xa = xc;
      }
      else
      {
        xb = xc;
      }
    }

    theta[k-1] = xc;
//
//  Seek root of equation 2 in interval.
//
    if ( m <= k )
    {
      break;
    }

    k = k + 1;
//
//  In the first interval, we need to skip the zero root of equation 2.
//
    if ( k == 2 )
    {
      k = k - 1;
    }
    else
    {
      xa = xa_init;
      fa = xa - tan ( a * xa ) / cl;
      ftol = eps * ( fabs ( fa ) + 1.0 );
      xb = xb_init;
      fb = - fa;

      while ( bmatol < xb - xa )
      {
        xc = ( xa + xb ) / 2.0;
        fc = xc - tan ( a * xc ) / cl;

        if ( fabs ( fc ) <= ftol )
        {
          break;
        }
        else if ( 0.0 < fc )
        {
          xa = xc;
        }
        else
        {
          xb = xc;
        }
      }
      theta[k-1] = xc;
    }
//
//  Advance the interval.
//
    xa_init = xa_init + pi / a;
    xb_init = xb_init + pi / a;
  }

  return theta;
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
