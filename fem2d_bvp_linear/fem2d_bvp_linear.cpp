# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "fem2d_bvp_linear.hpp"

//****************************************************************************80

double *fem2d_bvp_linear ( int nx, int ny, double a ( double x, double y ), 
  double c ( double x, double y ), double f ( double x, double y ), 
  double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_BVP_LINEAR solves a boundary value problem on a rectangle.
//
//  Discussion:
//
//    The procedure uses the finite element method, with piecewise linear basis
//    functions to solve a 2D boundary value problem over a rectangle
//
//    The following differential equation is imposed inside the region:
//
//      - d/dx a(x,y) du/dx - d/dy a(x,y) du/dy + c(x,y) * u(x,y) = f(x,y)
//
//    where a(x,y), c(x,y), and f(x,y) are given functions.
//
//    On the boundary, the solution is constrained to have the value 0.
//
//    The finite element method will use a regular grid of NX nodes in X, and 
//    NY nodes in Y.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of X and Y grid values.
//
//    Input, double A ( double X, double Y ), evaluates a(x,y);
//
//    Input, double C ( double X, double Y ), evaluates c(x,y);
//
//    Input, double F ( double X, double Y ), evaluates f(x,y);
//
//    Input, double X[NX], Y[NY], the grid coordinates.
//
//    Output, double FEM1D_BVP_LINEAR[NX*NY], the finite element coefficients, 
//    which are also the value of the computed solution at the mesh points.
//
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956, 
     0.000000000000000000000000000000, 
     0.774596669241483377035853079956 };
  double *amat;
  double aq;
  double *b;
  double cq;
  int e;
  int ex;
  int ey;
  double fq;
  int i;
  int ierror;
  int ii;
  int j;
  int jj;
  int k;
  int mn;
  int n;
  int ne;
  int nw;
  int quad_num = QUAD_NUM;
  int qx;
  int qy;
  int s;
  int se;
  int sw;
  int w;
  double weight[QUAD_NUM] = {
    0.555555555555555555555555555556, 
    0.888888888888888888888888888889, 
    0.555555555555555555555555555556 };
  double wq;
  double *u;
  double vne;
  double vnex;
  double vney;
  double vnw;
  double vnwx;
  double vnwy;
  double vse;
  double vsex;
  double vsey;
  double vsw;
  double vswx;
  double vswy;
  double xe;
  double xq;
  double xw;
  double yn;
  double yq;
  double ys;

  mn = nx * ny;

  amat = new double[mn*mn];
  b = new double[mn];

  for ( jj = 0; jj < mn; jj++ )
  {
    for ( ii = 0; ii < mn; ii++ )
    {
      amat[ii+jj*mn] = 0.0;
    }
  }
  for ( ii = 0; ii < mn; ii++ )
  {
    b[ii] = 0.0;
  }

  for ( ex = 0; ex < nx - 1; ex++ )
  {
    w = ex;
    e = ex + 1;

    xw = x[w];
    xe = x[e];

    for ( ey = 0; ey < ny - 1; ey++ )
    {
      s = ey;
      n = ey + 1;

      ys = y[s];
      yn = y[n];

      sw =   ey       * nx + ex;
      se =   ey       * nx + ex + 1;
      nw = ( ey + 1 ) * nx + ex;
      ne = ( ey + 1 ) * nx + ex + 1;

      for ( qx = 0; qx < quad_num; qx++ )
      {
        xq = ( ( 1.0 - abscissa[qx] ) * xw   
             + ( 1.0 + abscissa[qx] ) * xe ) 
               / 2.0;

        for ( qy = 0; qy < quad_num; qy++ )
        {
          yq = ( ( 1.0 - abscissa[qy] ) * ys   
               + ( 1.0 + abscissa[qy] ) * yn ) 
                 / 2.0;

          wq = weight[qx] * ( xe - xw ) / 2.0 
             * weight[qy] * ( yn - ys ) / 2.0;

          aq = a ( xq, yq );
          cq = c ( xq, yq );
          fq = f ( xq, yq );

          vsw  = ( xe - xq ) / ( xe - xw ) * ( yn - yq ) / ( yn - ys );
          vswx = (-1.0 ) / ( xe - xw ) * ( yn - yq ) / ( yn - ys );
          vswy = ( xe - xq ) / ( xe - xw ) * (-1.0 ) / ( yn - ys );

          vse  = ( xq - xw ) / ( xe - xw ) * ( yn - yq ) / ( yn - ys );
          vsex = ( 1.0 ) / ( xe - xw ) * ( yn - yq ) / ( yn - ys );
          vsey = ( xq - xw ) / ( xe - xw ) * (-1.0 ) / ( yn - ys );

          vnw  = ( xe - xq ) / ( xe - xw ) * ( yq - ys ) / ( yn - ys );
          vnwx = (-1.0 ) / ( xe - xw ) * ( yq - ys ) / ( yn - ys );
          vnwy = ( xe - xq ) / ( xe - xw ) * ( 1.0 ) / ( yn - ys );

          vne  = ( xq - xw ) / ( xe - xw ) * ( yq - ys ) / ( yn - ys );
          vnex = ( 1.0 ) / ( xe - xw ) * ( yq - ys ) / ( yn - ys );
          vney = ( xq - xw ) / ( xe - xw ) * ( 1.0 ) / ( yn - ys );

          amat[sw+sw*mn] = amat[sw+sw*mn] + wq * ( vswx * aq * vswx 
                                                 + vswy * aq * vswy 
                                                 + vsw  * cq * vsw );
          amat[sw+se*mn] = amat[sw+se*mn] + wq * ( vswx * aq * vsex 
                                                 + vswy * aq * vsey 
                                                 + vsw  * cq * vse );
          amat[sw+nw*mn] = amat[sw+nw*mn] + wq * ( vswx * aq * vnwx 
                                                 + vswy * aq * vnwy 
                                                 + vsw  * cq * vnw );
          amat[sw+ne*mn] = amat[sw+ne*mn] + wq * ( vswx * aq * vnex 
                                                 + vswy * aq * vney 
                                                 + vsw  * cq * vne );
          b[sw] =    b[sw]                + wq * ( vsw * fq );

          amat[se+sw*mn] = amat[se+sw*mn] + wq * ( vsex * aq * vswx 
                                                 + vsey * aq * vswy 
                                                 + vse  * cq * vsw );
          amat[se+se*mn] = amat[se+se*mn] + wq * ( vsex * aq * vsex 
                                                 + vsey * aq * vsey 
                                                 + vse  * cq * vse );
          amat[se+nw*mn] = amat[se+nw*mn] + wq * ( vsex * aq * vnwx 
                                                 + vsey * aq * vnwy 
                                                 + vse  * cq * vnw );
          amat[se+ne*mn] = amat[se+ne*mn] + wq * ( vsex * aq * vnex 
                                                 + vsey * aq * vney 
                                                 + vse  * cq * vne );
          b[se] =    b[se]                + wq * ( vse * fq );

          amat[nw+sw*mn] = amat[nw+sw*mn] + wq * ( vnwx * aq * vswx 
                                                 + vnwy * aq * vswy 
                                                 + vnw  * cq * vsw );
          amat[nw+se*mn] = amat[nw+se*mn] + wq * ( vnwx * aq * vsex 
                                                 + vnwy * aq * vsey 
                                                 + vnw  * cq * vse );
          amat[nw+nw*mn] = amat[nw+nw*mn] + wq * ( vnwx * aq * vnwx 
                                                 + vnwy * aq * vnwy 
                                                 + vnw  * cq * vnw );
          amat[nw+ne*mn] = amat[nw+ne*mn] + wq * ( vnwx * aq * vnex 
                                                 + vnwy * aq * vney 
                                                 + vnw  * cq * vne );
          b[nw] =    b[nw]                + wq * ( vnw * fq );

          amat[ne+sw*mn] = amat[ne+sw*mn] + wq * ( vnex * aq * vswx 
                                                 + vney * aq * vswy 
                                                 + vne  * cq * vsw );
          amat[ne+se*mn] = amat[ne+se*mn] + wq * ( vnex * aq * vsex 
                                                 + vney * aq * vsey 
                                                 + vne  * cq * vse );
          amat[ne+nw*mn] = amat[ne+nw*mn] + wq * ( vnex * aq * vnwx 
                                                 + vney * aq * vnwy 
                                                 + vne  * cq * vnw );
          amat[ne+ne*mn] = amat[ne+ne*mn] + wq * ( vnex * aq * vnex 
                                                 + vney * aq * vney 
                                                 + vne  * cq * vne );
          b[ne] =    b[ne]                + wq * ( vne * fq );
        }
      }
    }
  }
//
//  Where a node is on the boundary, 
//  replace the finite element equation by a boundary condition.
//
  k = 0;
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        for ( jj = 0; jj < mn; jj++ )
        {
          amat[k+jj*mn] = 0.0;
        }
        for ( ii = 0; ii < mn; ii++ )
        {
          amat[ii+k*mn] = 0.0;
        }
        amat[k+k*mn] = 1.0;
        b[k] = 0.0;
      }
      k = k + 1;
    }
  }
//
//  Solve the linear system.
//
  u = r8mat_solve2 ( mn, amat, b, ierror );

  delete [] amat;
  delete [] b;

  return u;
# undef QUAD_NUM
}
//****************************************************************************80

double fem2d_h1s_error_linear ( int nx, int ny, double x[], double y[],
  double u[], double exact_ux ( double x, double y ),
  double exact_uy ( double x, double y ) )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_H1S_ERROR_LINEAR: seminorm error of a finite element solution.
//
//  Discussion:
//
//    The finite element method has been used, over a rectangle,
//    involving a grid of NX*NY nodes, with piecewise linear elements used 
//    for the basis.
//
//    The finite element solution U(x,y) has been computed, and formulas for the
//    exact derivatives Vx(x,y) and Vy(x,y) are known.
//
//    This function estimates the H1 seminorm of the error:
//
//      H1S = sqrt ( integral ( x, y )   ( Ux(x,y) - Vx(x,y) )^2 
//                                     + ( Uy(x,y) - Vy(x,y) )^2 dx dy )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of X and Y grid values.
//
//    Input, double X[NX], Y[NY], the grid coordinates.
//
//    Input, double U[NX*NY], the finite element coefficients.
//
//    Input, function EXACT_UX(X,Y), EXACT_UY(X,Y) return the 
//    value of the derivatives of the exact solution with respect to
//    X and Y, respectively, at the point (X,Y).
//
//    Output, double FEM2D_H1S_ERROR_LINEAR, the estimated seminorm of 
//    the error.
//
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956, 
     0.000000000000000000000000000000, 
     0.774596669241483377035853079956 };
  int e;
  int ex;
  int ey;
  double exq;
  double eyq;
  double h1s;
  int mn;
  int n;
  int ne;
  int nw;
  int quad_num = QUAD_NUM;
  int qx;
  int qy;
  int s;
  int se;
  int sw;
  double uxq;
  double uyq;
  double vnex;
  double vney;
  double vnwx;
  double vnwy;
  double vsex;
  double vsey;
  double vswx;
  double vswy;
  int w;
  double weight[QUAD_NUM] = {
    0.555555555555555555555555555556, 
    0.888888888888888888888888888889, 
    0.555555555555555555555555555556 };
  double wq;
  double xe;
  double xq;
  double xw;
  double yn;
  double yq;
  double ys;

  mn = nx * ny;

  h1s = 0.0;

  for ( ex = 0; ex < nx - 1; ex++ )
  {
    w = ex;
    e = ex + 1;

    xw = x[w];
    xe = x[e];

    for ( ey = 0; ey < ny - 1; ey++ )
    {
      s = ey;
      n = ey + 1;

      ys = y[s];
      yn = y[n];

      sw =   ey       * nx + ex;
      se =   ey       * nx + ex + 1;
      nw = ( ey + 1 ) * nx + ex;
      ne = ( ey + 1 ) * nx + ex + 1;

      for ( qx = 0; qx < quad_num; qx++ )
      {
        xq = ( ( 1.0 - abscissa[qx] ) * xw   
             + ( 1.0 + abscissa[qx] ) * xe ) 
               / 2.0;

        for ( qy = 0; qy < quad_num; qy++ )
        {
          yq = ( ( 1.0 - abscissa[qy] ) * ys   
               + ( 1.0 + abscissa[qy] ) * yn ) 
                 / 2.0;

          wq = weight[qx] * ( xe - xw ) / 2.0 
             * weight[qy] * ( yn - ys ) / 2.0;

          vswx = (-1.0 ) / ( xe - xw ) * ( yn - yq ) / ( yn - ys );
          vswy = ( xe - xq ) / ( xe - xw ) * (-1.0 ) / ( yn - ys );

          vsex = ( 1.0 ) / ( xe - xw ) * ( yn - yq ) / ( yn - ys );
          vsey = ( xq - xw ) / ( xe - xw ) * (-1.0 ) / ( yn - ys );

          vnwx = (-1.0 ) / ( xe - xw ) * ( yq - ys ) / ( yn - ys );
          vnwy = ( xe - xq ) / ( xe - xw ) * ( 1.0 ) / ( yn - ys );

          vnex = ( 1.0 ) / ( xe - xw ) * ( yq - ys ) / ( yn - ys );
          vney = ( xq - xw ) / ( xe - xw ) * ( 1.0 ) / ( yn - ys );
//
//  Note that the south-west component of U is stored in U(W,S), not U(S,W)!
//
          uxq = u[w+s*nx] * vswx + u[e+s*nx] * vsex 
              + u[w+n*nx] * vnwx + u[e+n*nx] * vnex;
          uyq = u[w+s*nx] * vswy + u[e+s*nx] * vsey 
              + u[w+n*nx] * vnwy + u[e+n*nx] * vney;

          exq = exact_ux ( xq, yq );
          eyq = exact_uy ( xq, yq );

          h1s = h1s + wq * ( pow ( uxq - exq, 2 ) + pow ( uyq - eyq, 2 ) );
        }
      }
    }
  }

  h1s = sqrt ( h1s );

  return h1s;
# undef QUAD_NUM
}
//****************************************************************************80

double fem2d_l1_error ( int nx, int ny, double x[], double y[], double u[], 
  double exact ( double x, double y ) )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_L1_ERROR estimates the l1 error norm of a finite element solution.
//
//  Discussion:
//
//    The finite element method has been used, over a rectangle,
//    involving a grid of NX*NY nodes, with piecewise linear elements used 
//    for the basis.
//
//    The finite element coefficients have been computed, and a formula for the
//    exact solution is known.
//
//    This function estimates the little l1 norm of the error:
//      E1 = sum ( 1 <= I <= NX, 1 <= J <= NY ) 
//        abs ( U(i,j) - EXACT(X(i),Y(j)) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of X and Y grid values.
//
//    Input, double X[NX], Y[NY], the grid coordinates.
//
//    Input, double U[NX*NY], the finite element coefficients.
//
//    Input, function EQ = EXACT(X,Y), returns the value of the exact
//    solution at the point (X,Y).
//
//    Output, double FEM2D_L1_ERROR, the little l1 norm of the error.
//
{
  int i;
  int j;
  double e1;
  int mn;

  mn = nx * ny;

  e1 = 0.0;

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      e1 = e1 + fabs ( u[i+j*nx] - exact ( x[i], y[j] ) );
    }
  }

  e1 = e1 / ( double ) ( nx ) / ( double ) ( ny );

  return e1;
}
//****************************************************************************80

double fem2d_l2_error_linear ( int nx, int ny, double x[], double y[],
  double u[], double exact ( double x, double y ) )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_L2_ERROR_LINEAR: L2 error norm of a finite element solution.
//
//  Discussion:
//
//    The finite element method has been used, over a rectangle,
//    involving a grid of NX*NY nodes, with piecewise linear elements used 
//    for the basis.
//
//    The finite element coefficients have been computed, and a formula for the
//    exact solution is known.
//
//    This function estimates E2, the L2 norm of the error:
//
//      E2 = Integral ( X, Y ) ( U(X,Y) - EXACT(X,Y) )^2 dX dY
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of X and Y grid values.
//
//    Input, double X[NX], Y[NY], the grid coordinates.
//
//    Input, double U[NX*NY], the finite element coefficients.
//
//    Input, function EQ = EXACT(X,Y), returns the value of the exact
//    solution at the point (X,Y).
//
//    Output, double FEM2D_L2_ERROR_LINEAR, the estimated L2 norm of the error.
//
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956, 
     0.000000000000000000000000000000, 
     0.774596669241483377035853079956 };
  int e;
  double e2;
  double eq;
  int ex;
  int ey;
  int mn;
  int n;
  int ne;
  int nw;
  int quad_num = QUAD_NUM;
  int qx;
  int qy;
  int s;
  int se;
  int sw;
  double uq;
  double vne;
  double vnw;
  double vse;
  double vsw;
  int w;
  double weight[QUAD_NUM] = {
    0.555555555555555555555555555556, 
    0.888888888888888888888888888889, 
    0.555555555555555555555555555556 };
  double wq;
  double xe;
  double xq;
  double xw;
  double yn;
  double yq;
  double ys;

  mn = nx * ny;

  e2 = 0.0;
//
//  Integrate over each interval.
//
  for ( ex = 0; ex < nx - 1; ex++ )
  {
    w = ex;
    e = ex + 1;

    xw = x[w];
    xe = x[e];

    for ( ey = 0; ey < ny - 1; ey++ )
    {
      s = ey;
      n = ey + 1;

      ys = y[s];
      yn = y[n];

      sw =   ey       * nx + ex;
      se =   ey       * nx + ex + 1;
      nw = ( ey + 1 ) * nx + ex;
      ne = ( ey + 1 ) * nx + ex + 1;

      for ( qx = 0; qx < quad_num; qx++ )
      {
        xq = ( ( 1.0 - abscissa[qx] ) * xw  
             + ( 1.0 + abscissa[qx] ) * xe ) 
               / 2.0;

        for ( qy = 0; qy < quad_num; qy++ )
        {
          yq = ( ( 1.0 - abscissa[qy] ) * ys   
               + ( 1.0 + abscissa[qy] ) * yn ) 
                 / 2.0;

          wq = weight[qx] * ( xe - xw ) / 2.0
             * weight[qy] * ( yn - ys ) / 2.0;

          vsw  = ( xe - xq ) / ( xe - xw ) * ( yn - yq ) / ( yn - ys );
          vse  = ( xq - xw ) / ( xe - xw ) * ( yn - yq ) / ( yn - ys );
          vnw  = ( xe - xq ) / ( xe - xw ) * ( yq - ys ) / ( yn - ys );
          vne  = ( xq - xw ) / ( xe - xw ) * ( yq - ys ) / ( yn - ys );
//
//  Note that the south-west component of U is stored in U(W,S), not U(S,W)!
//
          uq = u[w+s*nx] * vsw + u[e+s*nx] * vse
             + u[w+n*nx] * vnw + u[e+n*nx] * vne;
          eq = exact ( xq, yq );

          e2 = e2 + wq * pow ( uq - eq, 2 );
        }
      }
    }
  }

  e2 = sqrt ( e2 );

  return e2;
# undef QUAD_NUM
}
//****************************************************************************80

int *i4vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO_NEW creates and zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
}
//****************************************************************************80

double *r8mat_solve2 ( int n, double a[], double b[], int &ierror )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE2 computes the solution of an N by N linear system.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The linear system may be represented as
//
//      A*X = B
//
//    If the linear system is singular, but consistent, then the routine will
//    still produce a solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of equations.
//
//    Input/output, double A[N*N].
//    On input, A is the coefficient matrix to be inverted.
//    On output, A has been overwritten.
//
//    Input/output, double B[N].
//    On input, B is the right hand side of the system.
//    On output, B has been overwritten.
//
//    Output, double R8MAT_SOLVE2[N], the solution of the linear system.
//
//    Output, int &IERROR.
//    0, no error detected.
//    1, consistent singularity.
//    2, inconsistent singularity.
//
{
  double amax;
  int i;
  int imax;
  int j;
  int k;
  int *piv;
  double *x;

  ierror = 0;

  piv = i4vec_zero_new ( n );
  x = r8vec_zero_new ( n );
//
//  Process the matrix.
//
  for ( k = 1; k <= n; k++ )
  {
//
//  In column K:
//    Seek the row IMAX with the properties that:
//      IMAX has not already been used as a pivot;
//      A(IMAX,K) is larger in magnitude than any other candidate.
//
    amax = 0.0;
    imax = 0;
    for ( i = 1; i <= n; i++ )
    {
      if ( piv[i-1] == 0 )
      {
        if ( amax < fabs ( a[i-1+(k-1)*n] ) )
        {
          imax = i;
          amax = fabs ( a[i-1+(k-1)*n] );
        }
      }
    }
//
//  If you found a pivot row IMAX, then,
//    eliminate the K-th entry in all rows that have not been used for pivoting.
//
    if ( imax != 0 )
    {
      piv[imax-1] = k;
      for ( j = k+1; j <= n; j++ )
      {
        a[imax-1+(j-1)*n] = a[imax-1+(j-1)*n] / a[imax-1+(k-1)*n];
      }
      b[imax-1] = b[imax-1] / a[imax-1+(k-1)*n];
      a[imax-1+(k-1)*n] = 1.0;

      for ( i = 1; i <= n; i++ )
      {
        if ( piv[i-1] == 0 )
        {
          for ( j = k+1; j <= n; j++ )
          {
            a[i-1+(j-1)*n] = a[i-1+(j-1)*n] - a[i-1+(k-1)*n] * a[imax-1+(j-1)*n];
          }
          b[i-1] = b[i-1] - a[i-1+(k-1)*n] * b[imax-1];
          a[i-1+(k-1)*n] = 0.0;
        }
      }
    }
  }
//
//  Now, every row with nonzero PIV begins with a 1, and
//  all other rows are all zero.  Begin solution.
//
  for ( j = n; 1 <= j; j-- )
  {
    imax = 0;
    for ( k = 1; k <= n; k++ )
    {
      if ( piv[k-1] == j )
      {
        imax = k;
      }
    }

    if ( imax == 0 )
    {
      x[j-1] = 0.0;

      if ( b[j-1] == 0.0 )
      {
        ierror = 1;
        cout << "\n";
        cout << "R8MAT_SOLVE2 - Warning:\n";
        cout << "  Consistent singularity, equation = " << j << "\n";
      }
      else
      {
        ierror = 2;
        cout << "\n";
        cout << "R8MAT_SOLVE2 - Warning:\n";
        cout << "  Inconsistent singularity, equation = " << j << "\n";
      }
    }
    else
    {
      x[j-1] = b[imax-1];

      for ( i = 1; i <= n; i++ )
      {
        if ( i != imax )
        {
          b[i-1] = b[i-1] - a[i-1+(j-1)*n] * x[j-1];
        }
      }
    }
  }

  delete [] piv;

  return x;
}
//****************************************************************************80

double *r8mat_zero_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_ZERO_NEW returns a new zeroed R8MAT.
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
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, double R8MAT_ZERO_NEW[M*N], the new zeroed matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

double *r8vec_even_new ( int n, double alo, double ahi )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EVEN_NEW returns an R8VEC of values evenly spaced between ALO and AHI.
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
//    18 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values.
//
//    Input, double ALO, AHI, the low and high values.
//
//    Output, double R8VEC_EVEN_NEW[N], N evenly spaced values.
//    Normally, A[0] = ALO and A[N-1] = AHI.
//    However, if N = 1, then A[0] = 0.5*(ALO+AHI).
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - i - 1 ) * alo
             + ( double ) (     i     ) * ahi )
             / ( double ) ( n     - 1 );
    }
  }

  return a;
}
//****************************************************************************80

double *r8vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO_NEW creates and zeroes an R8VEC.
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
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return a;
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
