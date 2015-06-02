# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "fem2d_bvp_serene.hpp"

//****************************************************************************80

double *basis_serene ( double xq, double yq, double xw, double ys, double xe, 
  double yn, double xx[], double yy[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_SERENE evaluates the serendipity basis functions.
//
//  Discussion:
//
//    This procedure assumes that a serendipity element has been defined,
//    whose sides are parallel to coordinate axes.
//
//    The local element numbering is
//
//      YN  3--2--1
//       |  |     |
//       |  4     8
//       |  |     |
//      YS  5--6--7
//       |
//       +--XW---XE--
//
//    We note that each basis function can be written as the product of
//    three linear terms, which never result in an x^2y^2 term.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XQ, YQ, the evaluation point.
//
//    Input, double XW, YS, the coordinates of the lower left corner.
//
//    Input, double XE, YN, the coordinates of the upper right corner.
//
//    Input, double XX[8], YY[8], the coordinates of the 8 nodes.
//
//    Output, double BASIS_SERENE[8], the value of the basis functions 
//    at (XQ,YQ).
//
{
  double *v;

  v = new double[8];

  v[0] = not1 ( xq, xw, xx[0] ) 
       * not1 ( yq, ys, yy[0] ) 
       * not2 ( xq, yq, xx[7], yy[7], xx[1], yy[1], xx[0], yy[0] );

  v[1] = not1 ( xq, xw, xx[1] ) 
       * not1 ( xq, xe, xx[1] ) 
       * not1 ( yq, ys, yy[1] );

  v[2] = not1 ( xq, xe, xx[2] ) 
       * not1 ( yq, ys, yy[2] ) 
       * not2 ( xq, yq, xx[1], yy[1], xx[3], yy[3], xx[2], yy[2] );

  v[3] = not1 ( xq, xe, xx[3] ) 
       * not1 ( yq, yn, yy[3] ) 
       * not1 ( yq, ys, yy[3] );

  v[4] = not1 ( xq, xe, xx[4] ) 
       * not1 ( yq, yn, yy[4] ) 
       * not2 ( xq, yq, xx[3], yy[3], xx[5], yy[5], xx[4], yy[4] );

  v[5] = not1 ( xq, xe, xx[5] ) 
       * not1 ( xq, xw, xx[5] ) 
       * not1 ( yq, yn, yy[5] );

  v[6] = not1 ( xq, xw, xx[6] ) 
       * not1 ( yq, yn, yy[6] ) 
       * not2 ( xq, yq, xx[5], yy[5], xx[7], yy[7], xx[6], yy[6] );

  v[7] = not1 ( yq, ys, yy[7] ) 
       * not1 ( yq, yn, yy[7] ) 
       * not1 ( xq, xw, xx[7] );

  return v;
}
//****************************************************************************80

double *basis_dx_serene ( double xq, double yq, double xw, double ys, 
  double xe, double yn, double xx[], double yy[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_DX_SERENE differentiates the serendipity basis functions for X.
//
//  Discussion:
//
//    This procedure assumes that a serendipity element has been defined,
//    whose sides are parallel to coordinate axes.
//
//    The local element numbering is
//
//      YN  3--2--1
//       |  |     |
//       |  4     8
//       |  |     |
//      YS  5--6--7
//       |
//       +--XW---XE--
//
//    We note that each basis function can be written as the product of
//    three linear terms, which never result in an x^2y^2 term.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XQ, YQ, the evaluation point.
//
//    Input, double XW, YS, the coordinates of the lower left corner.
//
//    Input, double XE, YN, the coordinates of the upper right corner.
//
//    Input, double XX[8], YY[8], the coordinates of the 8 nodes.
//
//    Output, double BASIS_DX_SERENE[8], the derivatives of the basis 
//    functions at (XQ,YQ) with respect to X.
//
{
  double *vx;

  vx = new double[8];

  vx[0] = 
      not1d ( xw, xx[0] ) 
    * not1 ( yq, ys, yy[0] ) 
    * not2 ( xq, yq, xx[7], yy[7], xx[1], yy[1], xx[0], yy[0] ) 
    + not1 ( xq, xw, xx[0] ) 
    * not1 ( yq, ys, yy[0] ) 
    * not2dx ( xx[7], yy[7], xx[1], yy[1], xx[0], yy[0] );

  vx[1] = 
      not1d ( xw, xx[1] ) 
    * not1 ( xq, xe, xx[1] ) 
    * not1 ( yq, ys, yy[1] ) 
    + not1 ( xq, xw, xx[1] ) 
    * not1d ( xe, xx[1] ) 
    * not1 ( yq, ys, yy[1] );

  vx[2] = 
      not1d ( xe, xx[2] ) 
    * not1 ( yq, ys, yy[2] ) 
    * not2 ( xq, yq, xx[1], yy[1], xx[3], yy[3], xx[2], yy[2] ) 
    + not1 ( xq, xe, xx[2] ) 
    * not1 ( yq, ys, yy[2] ) 
    * not2dx ( xx[1], yy[1], xx[3], yy[3], xx[2], yy[2] );

  vx[3] = 
      not1d ( xe, xx[3] ) 
    * not1 ( yq, yn, yy[3] ) 
    * not1 ( yq, ys, yy[3] );

  vx[4] = 
      not1d ( xe, xx[4] ) 
    * not1 ( yq, yn, yy[4] ) 
    * not2 ( xq, yq, xx[3], yy[3], xx[5], yy[5], xx[4], yy[4] ) 
    + not1 ( xq, xe, xx[4] ) 
    * not1 ( yq, yn, yy[4] ) 
    * not2dx ( xx[3], yy[3], xx[5], yy[5], xx[4], yy[4] );

  vx[5] = 
      not1d ( xe, xx[5] ) 
    * not1 ( xq, xw, xx[5] ) 
    * not1 ( yq, yn, yy[5] ) 
    + not1 ( xq, xe, xx[5] ) 
    * not1d ( xw, xx[5] ) 
    * not1 ( yq, yn, yy[5] );

  vx[6] = 
      not1d ( xw, xx[6] ) 
    * not1 ( yq, yn, yy[6] ) 
    * not2 ( xq, yq, xx[5], yy[5], xx[7], yy[7], xx[6], yy[6] ) 
    + not1 ( xq, xw, xx[6] ) 
    * not1 ( yq, yn, yy[6] ) 
    * not2dx ( xx[5], yy[5], xx[7], yy[7], xx[6], yy[6] );

  vx[7] = 
      not1 ( yq, ys, yy[7] ) 
    * not1 ( yq, yn, yy[7] ) 
    * not1d ( xw, xx[7] );

  return vx;
}
//****************************************************************************80

double *basis_dy_serene ( double xq, double yq, double xw, double ys, 
  double xe, double yn, double xx[], double yy[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_DY_SERENE differentiates the serendipity basis functions for Y.
//
//  Discussion:
//
//    This procedure assumes that a serendipity element has been defined,
//    whose sides are parallel to coordinate axes.
//
//    The local element numbering is
//
//      YN  3--2--1
//       |  |     |
//       |  4     8
//       |  |     |
//      YS  5--6--7
//       |
//       +--XW---XE--
//
//    We note that each basis function can be written as the product of
//    three linear terms, which never result in an x^2y^2 term.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XQ, YQ, the evaluation point.
//
//    Input, double XW, YS, the coordinates of the lower left corner.
//
//    Input, double XE, YN, the coordinates of the upper right corner.
//
//    Input, double XX[8], YY[8], the coordinates of the 8 nodes.
//
//    Output, double BASIS_DY_SERENE[8], the derivatives of the basis 
//    functions at (XQ,YQ) with respect to Y.
//
{
  double *vy;

  vy = new double[8];

  vy[0] = 
      not1 ( xq, xw, xx[0] ) 
    * not1d ( ys, yy[0] ) 
    * not2 ( xq, yq, xx[7], yy[7], xx[1], yy[1], xx[0], yy[0] ) 
    + not1 ( xq, xw, xx[0] ) 
    * not1 ( yq, ys, yy[0] ) 
    * not2dy ( xx[7], yy[7], xx[1], yy[1], xx[0], yy[0] );

  vy[1] = 
      not1 ( xq, xw, xx[1] ) 
    * not1 ( xq, xe, xx[1] ) 
    * not1d ( ys, yy[1] );

  vy[2] = not1 ( xq, xe, xx[2] ) 
    * not1d ( ys, yy[2] ) 
    * not2 ( xq, yq, xx[1], yy[1], xx[3], yy[3], xx[2], yy[2] ) 
    + not1 ( xq, xe, xx[2] ) 
    * not1 ( yq, ys, yy[2] ) 
    * not2dy ( xx[1], yy[1], xx[3], yy[3], xx[2], yy[2] );

  vy[3] = 
      not1 ( xq, xe, xx[3] ) 
    * not1d ( yn, yy[3] ) 
    * not1 ( yq, ys, yy[3] ) 
    + not1 ( xq, xe, xx[3] ) 
    * not1 ( yq, yn, yy[3] ) 
    * not1d ( ys, yy[3] );

  vy[4] = 
      not1 ( xq, xe, xx[4] ) 
    * not1d ( yn, yy[4] ) 
    * not2 ( xq, yq, xx[3], yy[3], xx[5], yy[5], xx[4], yy[4] ) 
    + not1 ( xq, xe, xx[4] ) 
    * not1 ( yq, yn, yy[4] ) 
    * not2dy ( xx[3], yy[3], xx[5], yy[5], xx[4], yy[4] );

  vy[5] = 
      not1 ( xq, xe, xx[5] ) 
    * not1 ( xq, xw, xx[5] ) 
    * not1d ( yn, yy[5] );

  vy[6] = 
      not1 ( xq, xw, xx[6] ) 
    * not1d ( yn, yy[6] ) 
    * not2 ( xq, yq, xx[5], yy[5], xx[7], yy[7], xx[6], yy[6] ) 
    + not1 ( xq, xw, xx[6] ) 
    * not1 ( yq, yn, yy[6] ) 
    * not2dy ( xx[5], yy[5], xx[7], yy[7], xx[6], yy[6] );

  vy[7] = 
      not1d ( ys, yy[7] ) 
    * not1 ( yq, yn, yy[7] ) 
    * not1 ( xq, xw, xx[7] ) 
    + not1 ( yq, ys, yy[7] ) 
    * not1d ( yn, yy[7] ) 
    * not1 ( xq, xw, xx[7] );

  return vy;
}
//****************************************************************************80

double *fem2d_bvp_serene ( int nx, int ny, double a ( double x, double y ), 
  double c ( double x, double y ), double f ( double x, double y ), 
  double x[], double y[], bool show11 )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_BVP_SERENE solves boundary value problem on a rectangle.
//
//  Discussion:
//
//    The program uses the finite element method, with piecewise 
//    serendipity basis functions to solve a 2D boundary value problem 
//    over a rectangle.
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
//    NY nodes in Y.  Both NX and NY must be odd.
//
//    The local element numbering is
//
//      3--2--1
//      |     |
//      4     8
//      |     |
//      5--6--7
//
//    The serendipity element mass matrix is a multiple of:
//
//       6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0
//      -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0
//       2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0
//      -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0
//       3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0
//      -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0
//       2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0
//      -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of X and Y grid values.
//    NX and NY must be odd and at least 3.
//
//    Input, function A(X,Y), evaluates a(x,y)
//
//    Input, function C(X,Y), evaluates c(x,y)
//
//    Input, function F(X,Y), evaluates f(x,y)
//
//    Input, double X[NX], Y[NY], the mesh points.
//
//    Input, bool SHOW11, is true to print out the element matrix
//    for the element in row 1, column 1.
//
//    Output, double FEM2D_BVP_SERENE[MN], the finite element coefficients, which 
//    are also the value of the computed solution at the mesh points.
//
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956, 
     0.000000000000000000000000000000, 
     0.774596669241483377035853079956 };
  double *ae;
  double *amat;
  double aq;
  double *b;
  double *be;
  int cc;
  double cq;
  int e;
  int ex;
  int ex_num;
  int ey;
  int ey_num;
  double fq;
  int i;
  int ierror;
  int ihi;
  int ii;
  int inc;
  int j;
  int jj;
  int k;
  int mm;
  int mn;
  int n;
  int node[8];
  int quad_num = QUAD_NUM;
  int qx;
  int qy;
  int s;
  double scale;
  double *u;
  double *v;
  double *vx;
  double *vy;
  int w;
  double weight[QUAD_NUM] = {
    0.555555555555555555555555555556, 
    0.888888888888888888888888888889, 
    0.555555555555555555555555555556 };
  double wq;
  double xc;
  double xe;
  double xq;
  double xw;
  double xx[8];
  double ym;
  double yn;
  double yq;
  double ys;
  double yy[8];
//
//  Make room for the matrix A and right hand side b.
//
  mn = fem2d_bvp_serene_node_num ( nx, ny );

  amat = r8mat_zero_new ( mn, mn );
  b = r8vec_zero_new ( mn );
//
//  Compute the matrix entries by integrating over each element.
//
  ex_num = ( nx - 1 ) / 2;
  ey_num = ( ny - 1 ) / 2;

  for ( ey = 0; ey < ey_num; ey++ )
  {
    s = 2 * ey;
    mm = 2 * ey + 1;
    n = 2 * ey + 2;

    ys = y[s];
    ym = y[mm];
    yn = y[n];

    yy[0] = y[n];
    yy[1] = y[n];
    yy[2] = y[n];
    yy[3] = y[mm];
    yy[4] = y[s];
    yy[5] = y[s];
    yy[6] = y[s];
    yy[7] = y[mm];

    for ( ex = 0; ex < ex_num; ex++ )
    {
      w = 2 * ex;
      cc = 2 * ex + 1;
      e = 2 * ex + 2;

      xe = x[e];
      xc = x[cc];
      xw = x[w];

      xx[0] = x[e];
      xx[1] = x[cc];
      xx[2] = x[w];
      xx[3] = x[w];
      xx[4] = x[w];
      xx[5] = x[cc];
      xx[6] = x[e];
      xx[7] = x[e];
//
//  Node indices
//
//  3  2  1  wn  cn  en
//  4     8  wm      em
//  5  6  7  ws  cs  es
//
      node[0] = ( 3 * ey + 3 ) * ey_num + 2 * ey + 2 * ex + 4;
      node[1] = ( 3 * ey + 3 ) * ey_num + 2 * ey + 2 * ex + 3;
      node[2] = ( 3 * ey + 3 ) * ey_num + 2 * ey + 2 * ex + 2;
      node[3] = ( 3 * ey + 2 ) * ey_num + 2 * ey +     ex + 1;
      node[4] = ( 3 * ey     ) * ey_num + 2 * ey + 2 * ex;
      node[5] = ( 3 * ey     ) * ey_num + 2 * ey + 2 * ex + 1;
      node[6] = ( 3 * ey     ) * ey_num + 2 * ey + 2 * ex + 2;
      node[7] = ( 3 * ey + 2 ) * ey_num + 2 * ey +     ex + 2;

      if ( show11 )
      {
        if ( ey == 0 && ex == 0 )
        {
          ae = r8mat_zero_new ( 8, 8 );
          be = r8vec_zero_new ( 8 );
        }
      }

      for ( qx = 0; qx < quad_num; qx++ )
      {
        xq = ( ( 1.0 - abscissa[qx] ) * x[e]   
             + ( 1.0 + abscissa[qx] ) * x[w] ) 
               / 2.0;

        for ( qy = 0; qy < quad_num; qy++ )
        {
          yq = ( ( 1.0 - abscissa[qy] ) * y[n]   
               + ( 1.0 + abscissa[qy] ) * y[s] ) 
                 / 2.0;

          wq = weight[qx] * ( x[e] - x[w] ) / 2.0 
             * weight[qy] * ( y[n] - y[s] ) / 2.0;

          v = basis_serene ( xq, yq, xw, ys, xe, yn, xx, yy );
          vx = basis_dx_serene ( xq, yq, xw, ys, xe, yn, xx, yy );
          vy = basis_dy_serene ( xq, yq, xw, ys, xe, yn, xx, yy );

          aq = a ( xq, yq );
          cq = c ( xq, yq );
          fq = f ( xq, yq );
//
//  Build the element matrix.
//
          if ( show11 )
          {
            if ( ey == 0 && ex == 0 )
            {
              for ( i = 0; i < 8; i++ )
              {
                for ( j = 0; j < 8; j++ )
                {
                  ae[i+j*8] = ae[i+j*8] + wq * ( vx[i] * aq * vx[j] 
                                               + vy[i] * aq * vy[j] 
                                               + v[i]  * cq * v[j] );
                }
                be[i] = be[i] + wq * ( v[i] * fq );
              }  
            }
          }

          for ( i = 0; i < 8; i++ )
          {
            ii = node[i];
            for ( j = 0; j < 8; j++ )
            {
              jj = node[j];
              amat[ii+jj*mn] = amat[ii+jj*mn] + wq * ( vx[i] * aq * vx[j] 
                                                     + vy[i] * aq * vy[j] 
                                                     + v[i]  * cq * v[j] );
            }
            b[ii] = b[ii] + wq * ( v[i] * fq );
          }

          delete [] v;
          delete [] vx;
          delete [] vy;
        }
      }
//
//  Print a sample element matrix.
//
      if ( show11 )
      {
        if ( ey == 0 && ex == 0 )
        {
          scale = 0.5 * ae[0+2*8];
          for ( j = 0; j < 8; j++ )
          {
            for ( i = 0; i < 8; i++ )
            {
              ae[i+j*8] = ae[i+j*8] / scale;
            }
          }
          r8mat_print ( 8, 8, ae, "  Wathen elementary mass matrix:" );
          delete [] ae;
          delete [] be;
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
    if ( ( j % 2 ) == 0 )
    {
      inc = 1;
    }
    else
    {
      inc = 2;
    }

    for ( i = 0; i < nx; i = i + inc )
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

int fem2d_bvp_serene_node_num ( int nx, int ny )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_BVP_SERENE_NODE_NUM counts the number of nodes.
//
//  Discussion:
//
//    The program uses the finite element method, with piecewise serendipity 
//    basis functions to solve a 2D boundary value problem over a rectangle.
//
//    The grid uses NX nodes in the X direction and NY nodes in the Y direction.
//
//    Both NX and NY must be odd.
//
//    Because of the peculiar shape of the serendipity elements, counting the
//    number of nodes and variables is a little tricky.  Here is a grid for
//    the case when NX = 7 and NY = 5, for which there are 29 nodes 
//    and variables.
//
//     23 24 25 26 27 28 29
//     19    20    21    22
//     12 13 14 15 16 17 18
//      8     9    10    11
//      1  2  3  4  5  6  7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of X and Y grid values.
//    NX and NY must be odd and at least 3.
//
//    Output, int FEM2D_BVP_SERENE_NODE_NUM, the number of nodes and variables.
//
{
  int value;

  value = 
        nx           * ( ny + 1 ) / 2 
    + ( nx + 1 ) / 2 * ( ny - 1 ) / 2;

  return value;
}
//****************************************************************************80

double fem2d_h1s_error_serene ( int nx, int ny, double x[], double y[], 
  double u[], double exact_ux ( double x, double y ), 
  double exact_uy ( double x, double y ) )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_H1S_ERROR_SERENE: seminorm error of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over a product region
//    involving a grid of NX*NY nodes, with serendipity elements used 
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
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of nodes.
//
//    Input, double X[NX], Y[NY], the grid coordinates.
//
//    Input, double U[*], the finite element coefficients.
//
//    Input, function EQX = EXACT_UX(X,Y), function EQY = EXACT_UY(X,Y)
//    returns the exact derivatives with respect to X and Y.
//
//    Output, double FEM2D_BVP_H1S, the estimated seminorm of the error.
//
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956, 
     0.000000000000000000000000000000, 
     0.774596669241483377035853079956 };
  int cc;
  int e;
  int ex;
  int ex_num;
  int ey;
  int ey_num;
  double exq;
  double eyq;
  double h1s;
  int k;
  int mm;
  int n;
  int node[8];
  int quad_num = QUAD_NUM;
  int qx;
  int qy;
  int s;
  double uxq;
  double uyq;
  double *vx;
  double *vy;
  int w;
  double weight[QUAD_NUM] = {
    0.555555555555555555555555555556, 
    0.888888888888888888888888888889, 
    0.555555555555555555555555555556 };
  double wq;
  double xc;
  double xe;
  double xq;
  double xw;
  double xx[8];
  double ym;
  double yn;
  double yq;
  double ys;
  double yy[8];

  h1s = 0.0;
//
//  Quadrature definitions.
//
  ex_num = ( nx - 1 ) / 2;
  ey_num = ( ny - 1 ) / 2;

  for ( ey = 0; ey < ey_num; ey++ )
  {
    s = 2 * ey;
    mm = 2 * ey + 1;
    n = 2 * ey + 2;

    ys = y[s];
    ym = y[mm];
    yn = y[n];

    yy[0] = y[n];
    yy[1] = y[n];
    yy[2] = y[n];
    yy[3] = y[mm];
    yy[4] = y[s];
    yy[5] = y[s];
    yy[6] = y[s];
    yy[7] = y[mm];

    for ( ex = 0; ex < ex_num; ex++ )
    {
      w = 2 * ex;
      cc = 2 * ex + 1;
      e = 2 * ex + 2;

      xe = x[e];
      xc = x[cc];
      xw = x[w];

      xx[0] = x[e];
      xx[1] = x[cc];
      xx[2] = x[w];
      xx[3] = x[w];
      xx[4] = x[w];
      xx[5] = x[cc];
      xx[6] = x[e];
      xx[7] = x[e];
//
//  Node indices
//
//  3  2  1  wn  cn  en
//  4     8  wm      em
//  5  6  7  ws  cs  es
//
      node[0] = ( 3 * ey + 3 ) * ey_num + 2 * ey + 2 * ex + 4;
      node[1] = ( 3 * ey + 3 ) * ey_num + 2 * ey + 2 * ex + 3;
      node[2] = ( 3 * ey + 3 ) * ey_num + 2 * ey + 2 * ex + 2;
      node[3] = ( 3 * ey + 2 ) * ey_num + 2 * ey +     ex + 1;
      node[4] = ( 3 * ey     ) * ey_num + 2 * ey + 2 * ex;
      node[5] = ( 3 * ey     ) * ey_num + 2 * ey + 2 * ex + 1;
      node[6] = ( 3 * ey     ) * ey_num + 2 * ey + 2 * ex + 2;
      node[7] = ( 3 * ey + 2 ) * ey_num + 2 * ey +     ex + 2;

      for ( qx = 0; qx < quad_num; qx++ )
      {
        xq = ( ( 1.0 - abscissa[qx] ) * x[e]   
             + ( 1.0 + abscissa[qx] ) * x[w] ) 
               / 2.0;

        for ( qy = 0; qy < quad_num; qy++ )
        {
          yq = ( ( 1.0 - abscissa[qy] ) * y[n]   
               + ( 1.0 + abscissa[qy] ) * y[s] ) 
                 / 2.0;

          wq = weight[qx] * ( x[e] - x[w] ) / 2.0 
             * weight[qy] * ( y[n] - y[s] ) / 2.0;

          vx = basis_dx_serene ( xq, yq, xw, ys, xe, yn, xx, yy );
          vy = basis_dy_serene ( xq, yq, xw, ys, xe, yn, xx, yy );

          uxq = 0.0;
          uyq = 0.0;
          for ( k = 0; k < 8; k++ )
          {
            uxq = uxq + u[node[k]] * vx[k];
            uyq = uyq + u[node[k]] * vy[k];
          }

          exq = exact_ux ( xq, yq );
          eyq = exact_uy ( xq, yq );

          h1s = h1s + wq * ( pow ( uxq - exq, 2 ) + pow ( uyq - eyq, 2 ) );
 
          delete [] vx;
          delete [] vy;
        }
      }
    }
  }

  h1s = sqrt ( h1s );

  return h1s;
# undef QUAD_NUM
}
//****************************************************************************80

double fem2d_l1_error_serene ( int nx, int ny, double x[], double y[], 
  double u[], double exact ( double x, double y ) )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_L1_ERROR_SERENE: l1 error norm of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over a product
//    region with NX*NY nodes and the serendipity element.
//
//    The coefficients U(*) have been computed, and a formula for the
//    exact solution is known.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
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
//    Input, double U[*], the finite element coefficients.
//
//    Input, function EQ = EXACT(X,Y), returns the value of the exact
//    solution at the point (X,Y).
//
//    Output, double FEM2D_L1_ERROR_SERENE, the little l1 norm of the error.
//
{
  int i;
  int inc;
  int j;
  int k;
  double e1;
  int mn;

  e1 = 0.0;
  k = 0;

  for ( j = 0; j < ny; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      inc = 1;
    }
    else
    {
      inc = 2;
    }
    for ( i = 0; i < nx; i = i + inc )
    {
      e1 = e1 + fabs ( u[k] - exact ( x[i], y[j] ) );
      k = k + 1;
    }
  }

  e1 = e1 / ( double ) ( k );

  return e1;
}
//****************************************************************************80

double fem2d_l2_error_serene ( int nx, int ny, double x[], double y[], 
  double u[], double exact ( double x, double y ) )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_L2_ERROR_SERENE: L2 error norm of a finite element solution.
//
//  Discussion:
//
//    The finite element method has been used, over a rectangle,
//    involving a grid of NXxNY nodes, with serendipity elements used 
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
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of nodes in the X 
//    and Y directions.
//
//    Input, double X[NX], Y[NY], the grid coordinates.
//
//    Input, double U[*], the finite element coefficients.
//
//    Input, function EQ = EXACT(X,Y), returns the value of the exact
//    solution at the point (X,Y).
//
//    Output, double E2, the estimated L2 norm of the error.
//
{
# define QUAD_NUM 3

  double abscissa[QUAD_NUM] = {
    -0.774596669241483377035853079956, 
     0.000000000000000000000000000000, 
     0.774596669241483377035853079956 };
  int cc;
  int e;
  double e2;
  double eq;
  int ex;
  int ex_num;
  int ey;
  int ey_num;
  int k;
  int mm;
  int n;
  int node[8];
  int quad_num = QUAD_NUM;
  int qx;
  int qy;
  int s;
  double uq;
  double *v;
  int w;
  double weight[QUAD_NUM] = {
    0.555555555555555555555555555556, 
    0.888888888888888888888888888889, 
    0.555555555555555555555555555556 };
  double wq;
  double xc;
  double xe;
  double xq;
  double xw;
  double xx[8];
  double ym;
  double yn;
  double yq;
  double ys;
  double yy[8];

  e2 = 0.0;
//
//  Compute the matrix entries by integrating over each element.
//
  ex_num = ( nx - 1 ) / 2;
  ey_num = ( ny - 1 ) / 2;

  for ( ey = 0; ey < ey_num; ey++ )
  {
    s = 2 * ey;
    mm = 2 * ey + 1;
    n = 2 * ey + 2;

    ys = y[s];
    ym = y[mm];
    yn = y[n];

    yy[0] = y[n];
    yy[1] = y[n];
    yy[2] = y[n];
    yy[3] = y[mm];
    yy[4] = y[s];
    yy[5] = y[s];
    yy[6] = y[s];
    yy[7] = y[mm];

    for ( ex = 0; ex < ex_num; ex++ )
    {
      w = 2 * ex;
      cc = 2 * ex + 1;
      e = 2 * ex + 2;

      xe = x[e];
      xc = x[cc];
      xw = x[w];

      xx[0] = x[e];
      xx[1] = x[cc];
      xx[2] = x[w];
      xx[3] = x[w];
      xx[4] = x[w];
      xx[5] = x[cc];
      xx[6] = x[e];
      xx[7] = x[e];
//
//  Node indices
//
//  3  2  1  wn  cn  en
//  4     8  wm      em
//  5  6  7  ws  cs  es
//
      node[0] = ( 3 * ey + 3 ) * ey_num + 2 * ey + 2 * ex + 4;
      node[1] = ( 3 * ey + 3 ) * ey_num + 2 * ey + 2 * ex + 3;
      node[2] = ( 3 * ey + 3 ) * ey_num + 2 * ey + 2 * ex + 2;
      node[3] = ( 3 * ey + 2 ) * ey_num + 2 * ey +     ex + 1;
      node[4] = ( 3 * ey     ) * ey_num + 2 * ey + 2 * ex;
      node[5] = ( 3 * ey     ) * ey_num + 2 * ey + 2 * ex + 1;
      node[6] = ( 3 * ey     ) * ey_num + 2 * ey + 2 * ex + 2;
      node[7] = ( 3 * ey + 2 ) * ey_num + 2 * ey +     ex + 2;

      for ( qx = 0; qx < quad_num; qx++ )
      {
        xq = ( ( 1.0 - abscissa[qx] ) * x[e]   
             + ( 1.0 + abscissa[qx] ) * x[w] ) 
               / 2.0;

        for ( qy = 0; qy < quad_num; qy++ )
        {
          yq = ( ( 1.0 - abscissa[qy] ) * y[n]   
               + ( 1.0 + abscissa[qy] ) * y[s] ) 
                 / 2.0;

          wq = weight[qx] * ( x[e] - x[w] ) / 2.0 
             * weight[qy] * ( y[n] - y[s] ) / 2.0;

          v = basis_serene ( xq, yq, xw, ys, xe, yn, xx, yy );

          uq = 0.0;
          for ( k = 0; k < 8; k++ )
          {
            uq = uq + u[node[k]] * v[k];
          }

          eq = exact ( xq, yq );
          e2 = e2 + wq * pow ( uq - eq, 2 );

          delete [] v;
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

double not1 ( double x1, double x2, double x3 )

//****************************************************************************80
//
//  Purpose:
//
//    NOT1 evaluates a factor for serendipity basis functions.
//
//  Discussion:
//
//    not1(x1,x2,x3) evaluates at the point x1, the basis factor that
//    is 0 at x2 and 1 at x3:
//
//      not1(x1,x2,x3) = ( x1 - x2 ) / ( x3 - x2 )  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X1, the evaluation point.
//
//    Input, double X2, X3, values that define the factor.
//
//    Output, double NOT1, the value of the basis function factor.
//
{
  double value;

  value = ( x1 - x2 ) / ( x3 - x2 );

  return value;
}
//****************************************************************************80

double not1d ( double x2, double x3 )

//****************************************************************************80
//
//  Purpose:
//
//    NOT1D differentiates a factor for serendipity basis functions.
//
//  Discussion:
//
//    not1(x1,x2,x3) evaluates at the point x1, the basis factor that
//    is 0 at x2 and 1 at x3:
//
//      not1(x1,x2,x3) = ( x1 - x2 ) / ( x3 - x2 )  
//
//    This function returns the derivative of the factor with respect to x1:
//
//      not1d(x1,x2,x3) = 1 / ( x3 - x2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X2, X3, values that define the factor.
//
//    Output, double NOT1D, the derivative of the basis function 
//    factor.
//
{
  double value;

  value = 1.0 / ( x3 - x2 );

  return value;
}
//****************************************************************************80

double not2 ( double x1, double y1, double x2, double y2, double x3, double y3, 
  double x4, double y4 )

//****************************************************************************80
//
//  Purpose:
//
//    NOT2 evaluates a factor for serendipity basis functions.
//
//  Discussion:
//
//    not2(x1,y1,x2,y2,x3,y3,x4,y4) evaluates at the point (x1,y1), the basis 
//    factor that is 0 at (x2,y2) and (x3,y3) and 1 at (x4,y4):
//
//          ( ( x1 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y1 - y2 ) )
//        / ( ( x4 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y4 - y2 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X1, Y2, the evaluation point.
//
//    Input, double X2, Y2, X3, Y3, values that define the factor.
//
//    Output, double NOT2, the value of the basis function factor.
//
{
  double value;

  value = ( ( x1 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y1 - y2 ) ) 
        / ( ( x4 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y4 - y2 ) );

  return value;
}
//****************************************************************************80

double not2dx ( double x2, double y2, double x3, double y3, double x4, double y4 )

//****************************************************************************80
//
//  Purpose:
//
//    NOT2DX evaluates a factor for serendipity basis functions.
//
//  Discussion:
//
//    not2(x1,y1,x2,y2,x3,y3,x4,y4) evaluates at the point (x1,y1), the basis 
//    factor that is 0 at (x2,y2) and (x3,y3) and 1 at (x4,y4):
//
//          ( ( x1 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y1 - y2 ) )
//        / ( ( x4 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y4 - y2 ) )
//
//    not2dx returns the derivative of this function with respect to X1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X2, Y2, X3, Y3, values that define the factor.
//
//    Output, double NOT2DX, the derivative of the basis function 
//    factor with respect to X1.
//
{
  double value;

  value = (   1.0       * ( y3 - y2 ) +                 0.0       ) 
        / ( ( x4 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y4 - y2 ) );

  return value;
}
//****************************************************************************80

double not2dy ( double x2, double y2, double x3, double y3, double x4, double y4 )

//****************************************************************************80
//
//  Purpose:
//
//    NOT2DY evaluates a factor for serendipity basis functions.
//
//  Discussion:
//
//    not2(x1,y1,x2,y2,x3,y3,x4,y4) evaluates at the point (x1,y1), the basis 
//    factor that is 0 at (x2,y2) and (x3,y3) and 1 at (x4,y4):
//
//          ( ( x1 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y1 - y2 ) )
//        / ( ( x4 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y4 - y2 ) )
//
//    not2dy returns the derivatives of this function with respect to Y1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X2, Y2, X3, Y3, values that define the factor.
//
//    Output, double NOT2DY, the derivative of the basis function 
//    factor with respect to Y1.
//
{
  double value;

  value = (   0.0                     - ( x3 - x2 ) *   1.0       ) 
        / ( ( x4 - x2 ) * ( y3 - y2 ) - ( x3 - x2 ) * ( y4 - y2 ) );

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
  const int i4_huge = 2147483647;
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

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
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
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
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
//    26 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }
    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( ihi < m )
    {
      i2hi = ihi;
    }
    else
    {
      i2hi = m;
    }

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
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
//
//    The value of N is
//      N = 3*NX*NY + 2*NX + 2*NY + 1,
//
//    and sufficient storage in A must have been set aside to hold
//    the matrix.
//
//    A is the consistent mass matrix for a regular NX by NY grid
//    of 8 node serendipity elements.  
//
//    The local element numbering is
//
//      3--2--1
//      |     |
//      4     8
//      |     |
//      5--6--7
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
//    The problem could be reprogrammed so that RHO is nonconstant,
//    but positive.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2014
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
  int seed;
  
  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }

  seed = 123456789;

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

      rho = 100.0 * r8_uniform_01 ( seed );

      for ( krow = 0; krow < 8; krow++ )
      {
        for ( kcol = 0; kcol < 8; kcol++ )
        {
          a[node[krow]+node[kcol]*n] = a[node[krow]+node[kcol]*n]
            + rho * em[krow+kcol*8];
        }
      }
    }
  }
  return a;
}
