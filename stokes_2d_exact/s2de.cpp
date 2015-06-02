# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>
# include <cstring>
# include <ctime>
# include <fstream>

using namespace std;

# include "s2de.hpp"

//****************************************************************************80

void grid_2d ( int x_num, double x_lo, double x_hi, int y_num, double y_lo, 
  double y_hi, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_2D returns a regular 2D grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of X values to use.
//
//    Input, double X_LO, X_HI, the range of X values.
//
//    Input, int Y_NUM, the number of Y values to use.
//
//    Input, double Y_LO, Y_HI, the range of Y values.
//
//    Output, double X[X_NUM*Y_NUM], Y[X_NUM*Y_NUM], 
//    the coordinates of the grid.
//
{
  int i;
  int j;
  double xi;
  double yj;

  if ( x_num == 1 )
  {
    for ( j = 0; j < y_num; j++ )
    {
      for ( i = 0; i < x_num; i++ )
      {
        x[i+j*x_num] = ( x_lo + x_hi ) / 2.0;
      }
    }
  }
  else
  {
    for ( i = 0; i < x_num; i++ )
    {
      xi = ( ( double ) ( x_num - i - 1 ) * x_lo   
           + ( double ) (         i     ) * x_hi ) 
           / ( double ) ( x_num     - 1 );
      for ( j = 0; j < y_num; j++ )
      {
        x[i+j*x_num] = xi;
      }
    }
  }

  if ( y_num == 1 )
  {
    for ( j = 0; j < y_num; j++ )
    {
      for ( i = 0; i < x_num; i++ )
      {
        y[i+j*x_num] = ( y_lo + y_hi ) / 2.0;
      }
    }
  }
  else
  {
    for ( j = 0; j < y_num; j++ )
    {
      yj = ( ( double ) ( y_num - j - 1 ) * y_lo   
           + ( double ) (         j     ) * y_hi ) 
           / ( double ) ( y_num     - 1 );
      for ( i = 0; i < x_num; i++ )
      {
        y[i+j*x_num] = yj;
      }
    }
  }

  return;
}
//****************************************************************************80

double r8vec_amax ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_AMAX returns the maximum absolute value in an R8VEC.
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
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], the array.
//
//    Output, double AMAX, the value of the entry
//    of largest magnitude.
//
{
  double amax;
  int i;

  amax = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( amax < fabs ( a[i] ) )
    {
      amax = fabs ( a[i] );
    }
  }

  return amax;
}
//****************************************************************************80

double r8vec_amin ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_AMIN returns the minimum absolute value in an R8VEC.
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
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], the array.
//
//    Output, double R8VEC_AMIN, the value of the entry
//    of smallest magnitude.
//
{
  int i;
  const double r8_huge = 1.79769313486231571E+308;
  double value;

  value = r8_huge;
  for ( i = 0; i < n; i++ )
  {
    if ( fabs ( a[i] ) < value )
    {
      value = fabs ( a[i] );
    }
  }

  return value;
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

double r8vec_min ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
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
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], the array to be checked.
//
//    Output, double R8VEC_MIN, the value of the minimum element.
//
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

double r8vec_norm_l2 ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], the vector whose L2 norm is desired.
//
//    Output, double R8VEC_NORM_L2, the L2 norm of A.
//
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
//****************************************************************************80

double *r8vec_uniform_ab_new ( int n, double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.
//
//  Discussion:
//
//    Each dimension ranges from A to B.
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
//    Input, int N, the number of entries in the vector.
//
//    Input, double A, B, the lower and upper limits of the pseudorandom values.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_AB_NEW - Fatal error!\n";
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

    r[i] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void resid_stokes1 ( int n, double x[], double y[], double ur[], double vr[], 
  double pr[] )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_STOKES1 returns residuals of the exact Stokes solution #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Junping Wang, Yanqiu Wang, Xiu Ye,
//    A robust numerical method for Stokes equations based on divergence-free
//    H(div) finite element methods,
//    SIAM Journal on Scientific Computing,
//    Volume 31, Number 4, 2009, pages 2784-2802.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the points.
//
//    Output, double UR[N], VR[N], PR[N], the residuals in the U, 
//    V and P equations.
//
{
  double *f;
  double *g;
  double *h;
  int i;
  double p;
  double px;
  double py;
  double u;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double v;
  double vx;
  double vxx;
  double vy;
  double vyy;
//
//  Get the right hand sides.
//
  f = new double[n];
  g = new double[n];
  h = new double[n];

  rhs_stokes1 ( n, x, y, f, g, h );
//
//  Form the functions and derivatives.
//
  for ( i = 0; i < n; i++ )
  {
    u = - 2.0 
          * pow ( x[i], 2 ) * pow ( x[i] - 1.0, 2 ) 
          * y[i] * ( y[i] - 1.0 ) * ( 2.0 * y[i] - 1.0 );

    ux = - 2.0 
          * ( 4.0 * pow ( x[i], 3 ) - 6.0 * pow ( x[i], 2 ) 
          + 2.0 * x[i] ) 
          * y[i] * ( y[i] - 1.0 ) * ( 2.0 * y[i] - 1.0 );

    uxx = - 2.0 
          * ( 12.0 * pow ( x[i], 2 ) - 12.0 * x[i] + 2.0 ) 
          * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] );

    uy = - 2.0  
          * pow ( x[i], 2 ) * pow ( x[i] - 1.0, 2 )  
          * ( 6.0 * pow ( y[i], 2 ) - 3.0 * y[i] + 1.0 );

    uyy = - 2.0 
          * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) ) 
          * ( 12.0 * y[i] - 6.0 );

    v =   2.0 
          * x[i] * ( x[i] - 1.0 ) * ( 2.0 * x[i] - 1.0 ) 
          * pow ( y[i], 2 ) * pow ( y[i] - 1.0, 2 );

    vx =   2.0 
          * ( 6.0 * pow ( x[i], 2 ) - 6.0 * x[i] + 1.0 ) 
          * pow ( y[i], 2 ) * pow ( y[i] - 1.0, 2 );

    vxx =   2.0 
          * ( 12.0 * x[i] - 6.0 ) 
          * pow ( y[i], 2 ) * pow ( y[i] - 1.0, 2 );

    vy =   2.0 
          * x[i] * ( x[i] - 1.0 ) * ( 2.0 * x[i] - 1.0 ) 
          * ( 4.0 * pow ( y[i], 3 ) - 6.0 * pow ( y[i], 2 )  
          + 2.0 * y[i] );

    vyy =   2.0 
          * x[i] * ( x[i] - 1.0 ) * ( 2.0 * x[i] - 1.0 ) 
          * ( 12.0 * pow ( y[i], 2 ) - 12.0 * y[i] + 2.0 );

    p = 0.0;
    px = 0.0;
    py = 0.0;

    ur[i] = px - ( uxx + uyy ) - f[i];
    vr[i] = py - ( vxx + vyy ) - g[i];
    pr[i] = ux + vy - h[i];
  }
//
//  Deallocate memory.
//
  delete [] f;
  delete [] g;
  delete [] h;

  return;
}
//****************************************************************************80

void resid_stokes2 ( int n, double x[], double y[], double ur[], double vr[], 
  double pr[] )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_STOKES2 returns residuals of the exact Stokes solution #2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Junping Wang, Yanqiu Wang, Xiu Ye,
//    A robust numerical method for Stokes equations based on divergence-free
//    H(div) finite element methods,
//    SIAM Journal on Scientific Computing,
//    Volume 31, Number 4, 2009, pages 2784-2802.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the points.
//
//    Output, double UR[N], VR[N], PR[N], the residuals in the U, 
//    V and P equations.
//
{
  double *f;
  double *g;
  double *h;
  int i;
  double p;
  double px;
  double py;
  const double r8_pi = 3.141592653589793;
  double u;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double v;
  double vx;
  double vxx;
  double vy;
  double vyy;
//
//  Get the right hand sides.
//
  f = new double[n];
  g = new double[n];
  h = new double[n];

  rhs_stokes2 ( n, x, y, f, g, h );
//
//  Form the functions and derivatives.
//
  for ( i = 0; i < n; i++ )
  {
    u =   2.0 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    ux =   4.0 * r8_pi 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    uxx = - 8.0 * pow ( r8_pi, 2 )
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    uy = - 4.0 * r8_pi 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    uyy = - 8.0 * pow ( r8_pi, 2 )
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    v = - 2.0 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vx =   4.0 * r8_pi 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vxx =   8.0 * pow ( r8_pi, 2 )
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vy = - 4.0 * r8_pi 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    vyy =   8.0 * pow ( r8_pi, 2 )
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    p = pow ( x[i], 2 ) + pow ( y[i], 2 );

    px = 2.0 * x[i];
    py = 2.0 * y[i];

    ur[i] = px - ( uxx + uyy ) - f[i];
    vr[i] = py - ( vxx + vyy ) - g[i];
    pr[i] = ux + vy - h[i];
  }
//
//  Deallocate memory.
//
  delete [] f;
  delete [] g;
  delete [] h;

  return;
}
//****************************************************************************80

void resid_stokes3 ( int n, double x[], double y[], double ur[], double vr[], 
  double pr[] )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_STOKES3 returns residuals of the exact Stokes solution #3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Howard Elman, Alison Ramage, David Silvester,
//    Finite Elements and Fast Iterative Solvers with
//    Applications in Incompressible Fluid Dynamics,
//    Oxford, 2005,
//    ISBN: 978-0198528678,
//    LC: QA911.E39.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the points.
//
//    Output, double UR[N], VR[N], PR[N], the residuals in the U, 
//    V and P equations.
//
{
  double *f;
  double *g;
  double *h;
  int i;
  double p;
  double px;
  double py;
  double u;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double v;
  double vx;
  double vxx;
  double vy;
  double vyy;
//
//  Get the right hand sides.
//
  f = new double[n];
  g = new double[n];
  h = new double[n];

  rhs_stokes3 ( n, x, y, f, g, h );
//
//  Form the functions and derivatives.
//
  for ( i = 0; i < n; i++ )
  {
    u =   20.0 * x[i] * pow ( y[i], 3 );
    ux = 20.0 * pow ( y[i], 3 );
    uxx = 0.0;
    uy = 60.0 * x[i] * pow ( y[i], 2 );
    uyy = 120.0 * x[i] * y[i];

    v = 5.0 * ( pow ( x[i], 4 )  - pow ( y[i], 4 ) );
    vx = 20.0 * pow ( x[i], 3 );
    vxx = 60.0 * pow ( x[i], 2 );
    vy = - 20.0 * pow ( y[i], 3 );
    vyy = - 60.0 * pow ( y[i], 2 );

    p =   60.0 * pow ( x[i], 2 ) * y[i] - 20.0 * pow ( y[i], 3 ) + 10.0;
    px = 120.0 * x[i] * y[i];
    py =  60.0 * pow ( x[i], 2 ) - 60.0 * pow ( y[i], 2 );

    ur[i] = px - ( uxx + uyy ) - f[i];
    vr[i] = py - ( vxx + vyy ) - g[i];
    pr[i] = ux + vy - h[i];
  }
//
//  Deallocate memory.
//
  delete [] f;
  delete [] g;
  delete [] h;

  return;
}
//****************************************************************************80

void rhs_stokes1 ( int n, double x[], double y[], double f[], double g[], 
  double h[] )

//****************************************************************************80
//
//  Purpose:
//
//    RHS_STOKES1 returns the right hand sides of the exact Stokes solution #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Junping Wang, Yanqiu Wang, Xiu Ye,
//    A robust numerical method for Stokes equations based on divergence-free
//    H(div) finite element methods,
//    SIAM Journal on Scientific Computing,
//    Volume 31, Number 4, 2009, pages 2784-2802.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the points.
//
//    Output, double F[N], G[N], H[N], the right hand sides in the U,
//    V and P equations.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = + 2.0 
          * ( 12.0 * pow ( x[i], 2 ) - 12.0 * x[i] + 2.0 ) 
          * ( 2.0 * pow ( y[i], 3 ) - 3.0 * pow ( y[i], 2 ) + y[i] ) 
          + 2.0 
          * ( pow ( x[i], 4 ) - 2.0 * pow ( x[i], 3 ) + pow ( x[i], 2 ) ) 
          * ( 12.0 * y[i] - 6.0 );

    g[i] = - 2.0 
          * ( 12.0 * x[i] - 6.0 ) 
          * ( pow ( y[i], 4 ) - 2.0 * pow ( y[i], 3 ) + pow ( y[i], 2 ) ) 
          - 2.0 
          * ( 2.0 * pow ( x[i], 3 ) - 3.0 * pow ( x[i], 2 ) + x[i] ) 
          * ( 12.0 * pow ( y[i], 2 ) - 12.0 * y[i] + 2.0 );

    h[i] = 0.0;
  }

  return;
}
//****************************************************************************80

void rhs_stokes2 ( int n, double x[], double y[], double f[], double g[], 
  double h[] )

//****************************************************************************80
//
//  Purpose:
//
//    RHS_STOKES2 returns the right hand sides of the exact Stokes solution #2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Junping Wang, Yanqiu Wang, Xiu Ye,
//    A robust numerical method for Stokes equations based on divergence-free
//    H(div) finite element methods,
//    SIAM Journal on Scientific Computing,
//    Volume 31, Number 4, 2009, pages 2784-2802.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the points.
//
//    Output, double F[N], G[N], H[N], the right hand sides in the U,
//    V and P equations.
//
{
  int i;
  double p;
  double px;
  double py;
  const double r8_pi = 3.141592653589793;
  double u;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double v;
  double vx;
  double vxx;
  double vy;
  double vyy;

  for ( i = 0; i < n; i++ )
  {
    u =   2.0 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    ux =   4.0 * r8_pi 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    uxx = - 8.0 * pow ( r8_pi, 2 )
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    uy = - 4.0 * r8_pi 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    uyy = - 8.0 * pow ( r8_pi, 2 )
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    v = - 2.0 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vx =   4.0 * r8_pi 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vxx =   8.0 * pow ( r8_pi, 2 )
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    vy = - 4.0 * r8_pi 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );

    vyy =   8.0 * pow ( r8_pi, 2 )
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    p = pow ( x[i], 2 ) + pow ( y[i], 2 );

    px = 2.0 * x[i];
    py = 2.0 * y[i];

    f[i] = px - ( uxx + uyy );
    g[i] = py - ( vxx + vyy );
    h[i] = ux + vy;
  }

  return;
}
//****************************************************************************80

void rhs_stokes3 ( int n, double x[], double y[], double f[], double g[], 
  double h[] )

//****************************************************************************80
//
//  Purpose:
//
//    RHS_STOKES3 returns the right hand sides of the exact Stokes solution #3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Howard Elman, Alison Ramage, David Silvester,
//    Finite Elements and Fast Iterative Solvers with
//    Applications in Incompressible Fluid Dynamics,
//    Oxford, 2005,
//    ISBN: 978-0198528678,
//    LC: QA911.E39.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the points.
//
//    Output, double F[N], G[N], H[N], the right hand sides in the U,
//    V and P equations.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f[i] = 0.0;
    g[i] = 0.0;
    h[i] = 0.0;
  }

  return;
}
//****************************************************************************80

void stokes_gnuplot ( string header, int n, double x[], double y[], double u[], 
  double v[], double s )

//****************************************************************************80
//
//  Purpose:
//
//    STOKES_GNUPLOT writes the Stokes velocity vector field to files for GNUPLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string HEADER, a header to be used to name the files.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the evaluation points.
//
//    Input, double U[N], V[N], the velocity components.
//
//    Input, double S, a scale factor for the velocity vectors.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  string plot_filename;
//
//  Write the data file.
//
  data_filename = header + "_data.txt";

  data_unit.open ( data_filename.c_str ( ) );

  for ( i = 0; i < n; i++ )
  {
    data_unit << "  " << x[i]
              << "  " << y[i]
              << "  " << s * u[i]
              << "  " << s * v[i] << "\n";
  }

  data_unit.close ( );

  cout << "\n";
  cout << "  Data written to '" << data_filename << "'\n";
//
//  Write the command file.
//
  command_filename = header + "_commands.txt";
  plot_filename = header + ".png";

  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "#  " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output '" << plot_filename << "'\n";
  command_unit << "#\n";
  command_unit << "#  Add titles and labels.\n";
  command_unit << "#\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set title 'Stokes velocity field'\n";
  command_unit << "unset key\n";
  command_unit << "#\n";
  command_unit << "#  Add grid lines.\n";
  command_unit << "#\n";
  command_unit << "set grid\n";
  command_unit << "set size ratio -1\n";
  command_unit << "#\n";
  command_unit << "#  Timestamp the plot.\n";
  command_unit << "#\n";
  command_unit << "set timestamp\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2:3:4 with vectors \\\n";
  command_unit << "  head filled lt 2 linecolor rgb 'blue'\n";
  command_unit << "quit\n";

  command_unit.close ( );

  cout << "  Commands written to '" << command_filename << "'\n";

  return;
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

void uvp_stokes1 ( int n, double x[], double y[], double u[], double v[], 
  double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_STOKES1 evaluates the exact Stokes solution #1.
//
//  Discussion:
//
//    The solution is defined over the unit square [0,1]x[0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Junping Wang, Yanqiu Wang, Xiu Ye,
//    A robust numerical method for Stokes equations based on divergence-free
//    H(div) finite element methods,
//    SIAM Journal on Scientific Computing,
//    Volume 31, Number 4, 2009, pages 2784-2802.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the points.
//
//    Output, double U[N], V[N], P[N], the velocity components and
//    pressure at each of the points.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {

    u[i] = - 2.0 
          * pow ( x[i], 2 ) * pow ( x[i] - 1.0, 2 )
          * y[i] * ( y[i] - 1.0 ) * ( 2.0 * y[i] - 1.0 );

    v[i] =   2.0 
          * x[i] * ( x[i] - 1.0 ) * ( 2.0 * x[i] - 1.0 ) 
          * pow ( y[i], 2 ) * pow ( y[i] - 1.0, 2 );

    p[i] = 0.0;
  }

  return;
}
//****************************************************************************80

void uvp_stokes2 ( int n, double x[], double y[], double u[], double v[], 
  double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_STOKES2 evaluates the exact Stokes solution #2.
//
//  Discussion:
//
//    The solution is defined over the unit square [0,1]x[0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Junping Wang, Yanqiu Wang, Xiu Ye,
//    A robust numerical method for Stokes equations based on divergence-free
//    H(div) finite element methods,
//    SIAM Journal on Scientific Computing,
//    Volume 31, Number 4, 2009, pages 2784-2802.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the points.
//
//    Output, double U[N], V[N], P[N], the velocity components and
//    pressure at each of the points.
//
{
  int i;
  const double r8_pi = 3.141592653589793;

  for ( i = 0; i < n; i++ )
  {
    u[i] =   2.0 
          * sin ( 2.0 * r8_pi * x[i] ) 
          * cos ( 2.0 * r8_pi * y[i] );


    v[i] = - 2.0 
          * cos ( 2.0 * r8_pi * x[i] ) 
          * sin ( 2.0 * r8_pi * y[i] );

    p[i] = pow ( x[i], 2 ) + pow ( y[i], 2 );
  }

  return;
}
//****************************************************************************80

void uvp_stokes3 ( int n, double x[], double y[], double u[], double v[], 
  double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    UVP_STOKES3 evaluates the exact Stokes solution #3.
//
//  Discussion:
//
//    The solution is defined over the unit square [-1,+1]x[-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Howard Elman, Alison Ramage, David Silvester,
//    Finite Elements and Fast Iterative Solvers with
//    Applications in Incompressible Fluid Dynamics,
//    Oxford, 2005,
//    ISBN: 978-0198528678,
//    LC: QA911.E39.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the points.
//
//    Output, double U[N], V[N], P[N], the velocity components and
//    pressure at each of the points.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    u[i] =   20.0 * x[i] * pow ( y[i], 3 );
    v[i] =    5.0 * ( pow ( x[i], 4 )  - pow ( y[i], 4 ) );
    p[i] =   60.0 * pow ( x[i], 2 ) * y[i] - 2.0 * pow ( y[i], 3 ) + 10.0;
  }

  return;
}
