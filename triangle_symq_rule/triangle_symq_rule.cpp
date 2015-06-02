# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "triangle_symq_rule.hpp"

//****************************************************************************80

void kjacopols ( double x, double a, double b, int n, double pols[] )

//****************************************************************************80
//
//  Purpose:
//
//    KJACOPOLS evaluates Jacobi polynomials.
//
//  Discussion:
//
//    This routine evaluates Jacobi polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double A, B, the parameter values.
//
//    Input, int N, the highest degree to be evaluated.
//
//    Output, double POLS[N+1], the polynomial values.
//
{
  double alpha;
  double beta;
  int k;
  double pk;
  double pkm1;
  double pkp1;

  pkp1 = 1.0;
  pols[0] = pkp1;

  if ( n == 0 )
  {
    return;
  }

  pk = pkp1;
  pkp1 = ( a / 2.0 - b / 2.0 ) 
    + ( 1.0 + a / 2.0 + b / 2.0 ) * x;
  pols[1] = pkp1;

  if ( n == 1 )
  {
    return;
  }

  for ( k = 2; k <= n; k++ )
  {
    pkm1 = pk;
    pk = pkp1;

    alpha = ( 2.0 * k + a + b - 1.0 ) 
      * ( a * a - b * b + ( 2.0 * k + a + b - 2.0 ) 
      * ( 2.0 * k + a + b ) * x );

    beta = 2.0 * ( k + a - 1.0 ) * ( k + b - 1.0 ) 
      * ( 2.0 * k + a + b );

    pkp1 = ( alpha * pk - beta * pkm1 ) 
      / ( 2.0 * k * ( k + a + b ) 
      * ( 2.0 * k + a + b - 2.0 ) );

    pols[k] = pkp1;
  }
  return;
}
//****************************************************************************80

void kjacopols2 ( double x, double a, double b, int n, double pols[], 
  double ders[] )

//****************************************************************************80
//
//  Purpose:
//
//    KJACOPOLS2 evaluates Jacobi polynomials and derivatives.
//
//  Discussion:
//
//    This routine evaluates Jacobi polynomials and  derivatives.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double A, B, the parameter values.
//
//    Input, int N, the highest degree to be evaluated.
//
//    Output, double POLS[N+1], the polynomial values.
//
//    Output, double DERS[N+1], the polynomial derivative values.
//
{
  double alpha1;
  double alpha2;
  double beta;
  double dk;
  double dkm1;
  double dkp1;
  double gamma;
  int k;
  double pk;
  double pkm1;
  double pkp1;

  pkp1 = 1.0;
  pols[0] = pkp1;

  dkp1 = 0.0;
  ders[0] = dkp1;

  if ( n == 0 )
  {
    return;
  }

  pk = pkp1;
  pkp1 = ( a / 2.0 - b / 2.0 ) 
    + ( 1.0 + a / 2.0 + b / 2.0 ) * x;
  pols[1] = pkp1;

  dk = dkp1;
  dkp1 = ( 1.0 + a / 2.0 + b / 2.0 );
  ders[1] = dkp1;

  if ( n == 1 )
  {
    return;
  }

  for ( k = 2; k <= n; k++ )
  {
    pkm1 = pk;
    pk = pkp1;
    dkm1 = dk;
    dk = dkp1;

    alpha1 = ( 2.0 * k + a + b - 1.0 ) * ( a * a - b * b );
    alpha2 = ( 2.0 * k + a + b - 1.0 ) 
      * ( ( 2.0 * k + a + b - 2.0 ) 
      * ( 2.0 * k + a + b ) );
    beta = 2.0 * ( k + a - 1.0 ) * ( k + b - 1.0 ) 
      * ( 2.0 * k + a + b );
    gamma = ( 2.0 * k * ( k + a + b ) 
      * ( 2.0 * k + a + b - 2.0 ) );
    pkp1 = ( ( alpha1 + alpha2 * x ) * pk - beta * pkm1 ) / gamma;
    dkp1 = ( ( alpha1 + alpha2 * x ) * dk 
      - beta * dkm1 + alpha2 * pk ) / gamma;

    pols[k] = pkp1;
    ders[k] = dkp1;
  }

  return;
}
//****************************************************************************80

double *klegeypols ( double x, double y, int n )

//****************************************************************************80
//
//  Purpose:
//
//    KLEGEYPOLS evaluates scaled Legendre polynomials.
//
//  Discussion:
//
//    This routine evaluate a sequence of scaled Legendre polynomials
//    P_n(x/y) y^n, with the parameter y in [0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double Y, the parameter.
//
//    Input, int N, the highest degree to be evaluated.
//
//    Output, double KLEGEYPOLS[N+1], the polynomial values.
//
{
  int k;
  double pk;
  double pkm1;
  double pkp1;
  double *pols;

  pols = new double[n+1];

  pkp1 = 1.0;
  pols[0] = pkp1;
  if ( n == 0 )
  {
    return pols;
  }

  pk = pkp1;
  pkp1 = x;
  pols[1] = pkp1;
  if ( n == 1 )
  {
    return pols;
  }

  for ( k = 1; k < n; k++ )
  {
    pkm1 = pk;
    pk = pkp1;
    pkp1 = ( ( 2.0 * k + 1.0 ) * x * pk - k * pkm1 * y * y ) / ( k + 1.0 );
    pols[k+1] = pkp1;
  }

  return pols;
}
//****************************************************************************80

void klegeypols3 ( double x, double y, int n, double pols[], double dersx[], 
  double dersy[] )

//****************************************************************************80
//
//  Purpose:
//
//    KLEGEYPOLS3 evaluates scaled Legendre polynomials and derivatives.
//
//  Discussion:
//
//    This routine evaluates a sequence of scaled Legendre polynomials
//    P_n(x/y) y^n, with the parameter y in [0,1], together with their
//    derivatives with respect to the parameters x and y.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double Y, the parameter value.
//
//    Input, int N, the highest degree to be evaluated.
//
//    Output, double POLS[N+1], the polynomial values.
//
//    Output, double DERSX[N+1], the derivatives with respect to X.
//
//    Output, double DERSY[N+1], the derivatives with respect to Y.
//
{
  double dk;
  double dkm1;
  double dkp1;
  int k;
  double pk;
  double pkm1;
  double pkp1;
  double yk;
  double ykm1;
  double ykp1;

  pkp1 = 1.0;
  pols[0] = pkp1;
  dkp1 = 0.0;
  dersx[0] = dkp1;
  ykp1 = 0.0;
  dersy[0] = ykp1;

  if ( n == 0 )
  {
    return;
  }

  pk = pkp1;
  pkp1 = x;
  pols[1] = pkp1;
  dk = dkp1;
  dkp1 = 1.0;
  dersx[1] = dkp1;
  yk = ykp1;
  ykp1 = 0.0;
  dersy[1] = ykp1;

  if ( n == 1 )
  {
    return;
  }

  for ( k = 1; k <= n - 1; k++ )
  {
    pkm1 = pk;
    pk = pkp1;
    dkm1 = dk;
    dk = dkp1;
    ykm1 = yk;
    yk = ykp1;
    pkp1 = ( ( 2.0 * k + 1.0 ) * x * pk - k * pkm1 * y * y ) / ( k + 1.0 );
    dkp1 = ( ( 2.0 * k + 1.0 ) * ( x * dk + pk ) 
      - k * dkm1 * y * y ) / ( k + 1.0 );
    ykp1 = ( ( 2.0 * k + 1.0 ) * ( x * yk ) 
      - k * ( pkm1 * 2.0 * y + ykm1 * y * y ) ) 
      / ( k + 1.0 );
    pols[k+1] = pkp1;
    dersx[k+1] = dkp1;
    dersy[k+1] = ykp1;
  }

  return;
}
//****************************************************************************80

double *ortho2eva ( int mmax, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    ORTHO2EVA evaluates orthogonal polynomials on the reference triangle.
//
//  Discussion:
//
//    This routine evaluates at the user-supplied point Z
//    a collection of polynomials (of X,Y) orthogonal on the
//    reference triangle.
//
//    The reference triangle has the vertices
//      (0,2/sqrt(3)), (-1,-1/sqrt(3)), (1,-1/sqrt(3))
//
//    The polynomials evaluated by this routine are all
//    orthogonal polynomials up to order mmax, arranged in the
//    increasing order.
//
//    This routine is based on the Koornwinder's representation
//    of the orthogonal polynomials on the right triangle
//      (-1,-1), (-1,1), (1,-1)
//    given by
//      K_mn(x,y) = P_m ((2*x+1+y)/(1-y)) ((1-y)/2)^m P_n^{2m+1,0} (y)
//    where P_m are the Legendre polynomials or order m
//    and P_n^{2m+1,0} are the Jacobi polynomials of order n with
//    the parameters alpha=2*m+1 and beta=0.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int MMAX, the maximum order to which the polynomials
//    are to be evaluated;
//
//    Input, double Z[2], the location where the polynomials are
//    to be evaluated; normally, expected to be inside (including boundary)
//    the reference triangle.
//
//    Output, double POLS[(mmax+1)*(mmax+2)/2], the orthogonal
//    polynomials evaluated at the point Z.
//
{
  double *pols;

  if ( mmax == 0 )
  {
    pols = new double[1];

    pols[0] = 1.0 / sqrt ( 3.0 ) * sqrt ( sqrt ( 3.0 ) );
  }
  else if ( mmax == 1 )
  {
    pols = new double[3];

    pols[0] = 1.0 / sqrt ( 3.0 ) * sqrt ( sqrt ( 3.0 ) );
    pols[1] = z[0] * sqrt ( 2.0 ) * sqrt ( sqrt ( 3.0 ) );
    pols[2] = z[1] * sqrt ( 2.0 ) * sqrt ( sqrt ( 3.0 ) );
  }
  else
  {
    pols = ortho2eva0 ( mmax, z );
  }
  return pols;
}
//****************************************************************************80

double *ortho2eva0 ( int mmax, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    ORTHO2EVA0 evaluates the orthonormal polynomials on the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int MMAX, the maximum order to which the polynomials are
//    to be evaluated.
//
//    Input, double Z[2], the coordinates of the evaluation point.
//
//    Output, double ORTHO2EVA0[(mmax+1)*(mmax+2)/2], the orthogonal
//    polynomials evaluated at the point Z.
//
{
  double a;
  double b;
  double *f1;
  double *f2;
  int kk;
  int m;
  int n;
  int npols;
  double par1;
  double par2;
  double *pols;
  double r11;
  double r12;
  double r21;
  double r22;
  double scale;
  double sqrt2;
  double sqrt3;
  double x;
  double y;
  double zero;

  npols = ( ( mmax + 1 ) * ( mmax + 2 ) ) / 2;
  pols = new double[npols];

  zero = 0.0;
  sqrt2 = sqrt ( 2.0 );
  sqrt3 = sqrt ( 3.0 );
  r11 = -1.0 / 3.0;
  r12 = -1.0 / sqrt3;
  r21 = - 1.0 / 3.0;
  r22 = 2.0 / sqrt3;

  a = z[0];
  b = z[1];
//
//  Map the reference triangle to the right
//  triangle with the vertices (-1,-1), (1,-1), (-1,1)
//
  x = r11 + r12 * b + a;
  y = r21 + r22 * b;
//
//  Evaluate the Koornwinder's polynomials via the three term recursion.
//
  par1 = ( 2.0 * x + 1.0 + y ) / 2.0;
  par2 = ( 1.0 - y ) / 2.0;

  f1 = klegeypols ( par1, par2, mmax );

  f2 = new double[(mmax+1)*(mmax+1)];

  for ( m = 0; m <= mmax; m++ )
  {
    par1 = 2 * m + 1;
    kjacopols ( y, par1, zero, mmax - m, f2+m*(mmax+1) );
  }

  kk = 0;
  for ( m = 0; m <= mmax; m++ )
  {
    for ( n = 0; n <= m; n++ )
    {
//
//  Evaluate the polynomial (m-n, n)
//
      pols[kk] = f1[m-n] * f2[n+(m-n)*(mmax+1)];
//
//  Normalize.
//
      scale = sqrt 
      ( ( double ) 
        ( 
          ( 1 + ( m - n ) + n ) * ( 1 + ( m - n ) + ( m - n ) )
        ) / sqrt3 
      );

      pols[kk] = pols[kk] * scale;
      kk = kk + 1;
    }
  }

  delete [] f1;
  delete [] f2;

  return pols;
}
//****************************************************************************80

void quaecopy ( int nf, double xs[], double ys[], double ws[], double z[], 
  double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAECOPY copies a quadrature rule into user arrays Z and W.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int NF, the number of values to copy.
//
//    Input, double XS[NF], YS[NF], the point coordinates to copy.
//
//    Input, double WS[NF], the weights to copy.
//
//    Output, double Z[2*NF], the copied point coordinates.
//
//    Output, double W[NF], the copied weights.
//
{
  int j;

  for ( j = 0; j < nf; j++ )
  {
    z[0+j*2] = xs[j];
    z[1+j*2] = ys[j];
  }
  r8vec_copy ( nf, ws, w );

  return;
}
//****************************************************************************80

void quaecopy2 ( double xs[], double ys[], double ws[], double xnew[], 
  double ynew[], double w[], int kk )

//****************************************************************************80
//
//  Purpose:
//
//    QUAECOPY2 copies a quadrature rule into a user arrays X, Y, and W.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double XS[KK], YS[KK], the point coordinates to copy.
//
//    Input, double WS[KK], the weights to copy.
//
//    Output, double XNEW[KK], YNEW[KK], the copied point coordinates.
//
//    Output, double W[KK], the copied weights.
//
//    Input, int KK, the number of values to copy.
//
{
  r8vec_copy ( kk, xs, xnew );
  r8vec_copy ( kk, ys, ynew );
  r8vec_copy ( kk, ws, w );

  return;
}
//****************************************************************************80

int quaeinside ( int iitype, double xsout, double ysout )

//****************************************************************************80
//
//  Purpose:
//
//    QUAEINSIDE checks whether a point is inside a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int IITYPE, indicates the check to perform.
//    * 0 check whether the point is inside the whole triangle
//    * 1 check whether the point is inside the bottom 1/3 of the triangle
//    * 2 check whether the point is inside the lower-left 1/6 of the triangle
//
//    Input, double XSOUT, YSOUT, the coordinates of the point.
//
//    Output, int QUAEINSIDE.
//    * 1, the point is inside.
//    * 2, the point is outside.
//
{
  int nbool;
  double s;

  s = sqrt ( 3.0 );
//
//  The 1/6 of triangle.
//
  if ( iitype == 2 )
  {
    nbool = 1;
    if ( xsout < -1.0 || 0.0 < xsout )
    {
      nbool = 0;
    }
//
//  hx
//
    if ( ysout < -1.0 / s - 1.0E-30 || xsout / s < ysout )
    {
      nbool = 0;
    }
  }
//
//  The 1/3 of triangle.
//
  else if ( iitype == 1 )
  {
    nbool = 0;

    if ( xsout <= 0.0 && -1.0 <= xsout &&  -1.0 / s <= ysout && ysout <= xsout / s )
    {
      nbool = 1;
    }
    if ( iitype == 1 )
    {
      if  ( 0.0 <= xsout && xsout <= 1.0 &&
            -1.0 / s <= ysout && ysout <= - xsout / s )
      {
        nbool = 1;
      }
    }
  }
//
//  The entire triangle.
//
  else if ( iitype == 0 )
  {
    nbool = 1;
    if ( ysout < -1.0 / s )
    {
      nbool = 0;
    }
    if ( s * xsout + 2.0 / s < ysout )
    {
      nbool = 0;
    }
    if ( - s * xsout + 2.0 / s < ysout )
    {
      nbool = 0;
    }
  }
  return nbool;
}
//****************************************************************************80

void quaenodes ( int nptsout, double xsout[], double ysout[], double wsout[], 
  int nptsoutout, double xs2[], double ys2[], double ws2[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAENODES expands nodes to the reference triangle.
//
//  Discussion:
//
//    This routine expands nodes to the reference triangle
//    assuming that the points are already in the lower-left 1/6 of the
//    triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int NPTSOUT, the number of points expanded.
//
//    Input, double XSOUT[NPTSOUT], YSOUT[NPTSOUT], WSOUT[NPTSOUT],
//    the coordinates and weights of the points to be expanded.
//
//    Output, int NPTSOUTPUT, the number of points in the expanded set.
//
//    Output, double XS2[NPTSOUTOUT], YS2[NPTSOUTOUT], WS2[NPTSOUTOUT],
//    the coordinates and weights of the expanded set of points.
//
{
  const double eps = 1.0E-12;
  int i;
  int ntot;
  double w0;
  double x0;
  double x1;
  double x2;
  double y0;
  double y1;
  double y2;

  ntot = 0;

  for ( i = 0; i < nptsout; i++ )
  {
    if ( pow ( xsout[i], 2 ) + pow ( ysout[i], 2 ) < eps )
    {
      xs2[ntot] = xsout[i];
      ys2[ntot] = ysout[i];
      ws2[ntot] = wsout[i];
      ntot = ntot + 1;
    }
    else if ( pow ( xsout[i], 2 ) < eps ||
      fabs ( ysout[i] - xsout[i] / sqrt ( 3.0 ) ) < sqrt ( eps ) )
    {
      x0 = xsout[i];
      y0 = ysout[i];
      w0 = wsout[i] / 3.0;

      xs2[ntot] = x0;
      ys2[ntot] = y0;
      ws2[ntot] = w0;
      ntot = ntot + 1;

      quaerotate ( x0, y0, x1, y1 );
      xs2[ntot] = x1;
      ys2[ntot] = y1;
      ws2[ntot] = w0;
      ntot = ntot + 1;

      quaerotate ( x1, y1, x2, y2 );
      xs2[ntot] = x2;
      ys2[ntot] = y2;
      ws2[ntot] = w0;
      ntot = ntot + 1;
    }
    else
    {
      x0 = xsout[i];
      y0 = ysout[i];
      w0 = wsout[i] / 6.0;

      xs2[ntot] = x0;
      ys2[ntot] = y0;
      ws2[ntot] = w0;
      ntot = ntot + 1;

      quaerotate ( x0, y0, x1, y1 );
      xs2[ntot] = x1;
      ys2[ntot] = y1;
      ws2[ntot] = w0;
      ntot = ntot + 1;

      quaerotate ( x1, y1, x2, y2 );
      xs2[ntot] = x2;
      ys2[ntot] = y2;
      ws2[ntot] = w0;
      ntot = ntot + 1;

      xs2[ntot] = -x0;
      ys2[ntot] = y0;
      ws2[ntot] = w0;
      ntot = ntot + 1;

      quaerotate ( -x0, y0, x1, y1 );
      xs2[ntot] = x1;
      ys2[ntot] = y1;
      ws2[ntot] = w0;
      ntot = ntot + 1;

      quaerotate ( x1, y1, x2, y2 );
      xs2[ntot] = x2;
      ys2[ntot] = y2;
      ws2[ntot] = w0;
      ntot = ntot + 1;
    }
  }

  return;
}
//****************************************************************************80

void quaenodes2 ( int nptsout, double xsout[], double ysout[], double wsout[], 
  int nptsoutout, double xs2[], double ys2[], double ws2[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAENODES2 expands nodes from 1/6 to 1/3 of the triangle.
//
//  Discussion:
//
//    This routine only expands to 1/3 of the triangle, assuming the points are
//    already in the lower-left 1/6 of the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int NPTSOUT, the number of points expanded.
//
//    Input, double XSOUT[NPTSOUT], YSOUT[NPTSOUT], WSOUT[NPTSOUT],
//    the coordinates and weights of the points to be expanded.
//
//    Output, int NPTSOUTPUT, the number of points in the expanded set.
//
//    Output, double XS2[NPTSOUTOUT], YS2[NPTSOUTOUT], WS2[NPTSOUTOUT],
//    the coordinates and weights of the expanded set of points.
//
{
  const double eps = 1.0E-12;
  int i;
  int ntot;

  ntot = 0;

  for ( i = 0; i < nptsout; i++ )
  {
    if ( pow ( xsout[i], 2 ) + pow ( ysout[i], 2 ) < eps )
    {
      xs2[ntot] = xsout[i];
      ys2[ntot] = ysout[i];
      ws2[ntot] = wsout[i];
      ntot = ntot + 1;
    }
    else if ( pow ( xsout[i], 2 ) < eps )
    {
      xs2[ntot] = xsout[i];
      ys2[ntot] = ysout[i];
      ws2[ntot] = wsout[i];
      ntot = ntot + 1;
    }
    else
    {
      xs2[ntot] = -xsout[i];
      ys2[ntot] = ysout[i];
      ws2[ntot] = wsout[i] / 2.0;
      ntot = ntot + 1;

      xs2[ntot] = xsout[i];
      ys2[ntot] = ysout[i];
      ws2[ntot] = wsout[i] / 2.0;
      ntot = ntot + 1;
    }
  }

  return;
}
//***************************************************************************80

void quaequad ( int itype, int mmax, double zs[], double whts[], int numnodes )

//****************************************************************************80
//
//  Purpose:
//
//    QUAEQUAD returns a symmetric quadrature formula for a reference triangle.
//
//  Discussion:
//
//    This routine constructs (or rather, retrieves)
//    D_3 symmetric quadrature formulae for smooth functions
//    on the triangle with vertices
//      (-1,-1/sqrt(3)), (1,-1/sqrt(3)), (0,2/sqrt(3)).
//
//    All quadratures are obtained to the extended precision.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int ITYPE, the configuration of the quadrature nodes.
//    * 0, all nodes on the reference triangle
//    * 1, quadrature nodes on the bottom 1/3 of the reference triangle
//    * 2, quadrature nodes on lower-left 1/6 of the reference triangle
//
//    Input, int MMAX, the degree of the quadrature (the maximum degree of
//    the polynomials of two variables that are integrated
//    exactly. 1 <= MMAX <= 50.
//
//    Output, double ZS[2*NUMNODES], the nodes.
//
//    Output, double WHTS[NUMNODES], the weights.
//
//    Input, int NUMNODES, the number of nodes in the quadrature rule.
//
{
  int nc;
  int nf;
  double *wc;
  double *wf;
  double *xc;
  double *xf;
  double *yc;
  double *yf;

  if ( mmax < 1 || 50 < mmax )
  {
    cerr << "\n";
    cerr << "QUAEQUAD - Fatal error!\n";
    cerr << "  1 <= MMAX <= 50 is required.\n";
    exit ( 1 );
  }

  if ( itype < 0 || 2 < itype )
  {
    cerr << "\n";
    cerr << "QUAEQUAD - Fatal error!\n";
    cerr << "  0 <= ITYPE <= 2 is required.\n";
    exit ( 1 );
  }
//
//  Get the size of the compressed rule.
//
  nc = rule_compressed_size ( mmax );
//
//  Retrieve the compressed rule.
//
  xc = new double[nc];
  yc = new double[nc];
  wc = new double[nc];

  quaequad0 ( mmax, nc, xc, yc, wc );
//
//  Expand the nodes to the entire triangle.
//
  nf = rule_full_size ( mmax );
  xf = new double[nf];
  yf = new double[nf];
  wf = new double[nf];

  if ( itype == 0 )
  {
    quaenodes ( nc, xc, yc, wc, nf, xf, yf, wf );
  }
//
//  Expand the nodes to the lower 1/3 of the triangle.
//
  else if ( itype == 1 )
  {
    quaenodes2 ( nc, xc, yc, wc, nf, xf, yf, wf );
  }
//
//  Simply copy the nodes; they are already in the lower-left
//  1/6 of the triangle.
//
  else if ( itype == 2 )
  {
    r8vec_copy ( nc, xc, xf );
    r8vec_copy ( nc, yc, yf );
    r8vec_copy ( nc, wc, wf );
  }

  quaecopy ( nf, xf, yf, wf, zs, whts );

  delete [] xc;
  delete [] xf;
  delete [] yc;
  delete [] yf;
  delete [] wc;
  delete [] wf;

  return;
}
//****************************************************************************80

void quaequad0 ( int mmax, int kk, double xnew[], double ynew[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAEQUAD0 returns the requested quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int MMAX, the degree of the quadrature (the 
//    maximum degree of the polynomials of two variables that are integrated
//    exactly.  1 <= MMAX <= 50.
//
//    Input, int KK, the size of the compressed rule.
//
//    Output, double XNEW[KK], YNEW[KK], the coordinates of the nodes.
//
//    Output, double W[KK], the weights.
//
{
  int i;
  int iitype;
  int nbool0;
  int nbool1;
  int nbool2;
  double x0;
  double x1;
  double y0;
  double y1;
//
//  Copy the arrays defining the compressed rule.
//
  if ( mmax == 1 )
  {
    rule01 ( xnew, ynew, w );
  }
  else if ( mmax == 2 )
  {
    rule02 ( xnew, ynew, w );
  }
  else if ( mmax == 3 )
  {
    rule03 ( xnew, ynew, w );
  }
  else if ( mmax == 4 )
  {
    rule04 ( xnew, ynew, w );
  }
  else if ( mmax == 5 )
  {
    rule05 ( xnew, ynew, w );
  }
  else if ( mmax == 6 )
  {
    rule06 ( xnew, ynew, w );
  }
  else if ( mmax == 7 )
  {
    rule07 ( xnew, ynew, w );
  }
  else if ( mmax == 8 )
  {
    rule08 ( xnew, ynew, w );
  }
  else if ( mmax == 9 )
  {
    rule09 ( xnew, ynew, w );
  }
  else if ( mmax == 10 )
  {
    rule10 ( xnew, ynew, w );
  }
  else if ( mmax == 11 )
  {
    rule11 ( xnew, ynew, w );
  }
  else if ( mmax == 12 )
  {
    rule12 ( xnew, ynew, w );
  }
  else if ( mmax == 13 )
  {
    rule13 ( xnew, ynew, w );
  }
  else if ( mmax == 14 )
  {
    rule14 ( xnew, ynew, w );
  }
  else if ( mmax == 15 )
  {
    rule15 ( xnew, ynew, w );
  }
  else if ( mmax == 16 )
  {
    rule16 ( xnew, ynew, w );
  }
  else if ( mmax == 17 )
  {
    rule17 ( xnew, ynew, w );
  }
  else if ( mmax == 18 )
  {
    rule18 ( xnew, ynew, w );
  }
  else if ( mmax == 19 )
  {
    rule19 ( xnew, ynew, w );
  }
  else if ( mmax == 20 )
  {
    rule20 ( xnew, ynew, w );
  }
  else if ( mmax == 21 )
  {
    rule21 ( xnew, ynew, w );
  }
  else if ( mmax == 22 )
  {
    rule22 ( xnew, ynew, w );
  }
  else if ( mmax == 23 )
  {
    rule23 ( xnew, ynew, w );
  }
  else if ( mmax == 24 )
  {
    rule24 ( xnew, ynew, w );
  }
  else if ( mmax == 25 )
  {
    rule25 ( xnew, ynew, w );
  }
  else if ( mmax == 26 )
  {
    rule26 ( xnew, ynew, w );
  }
  else if ( mmax == 27 )
  {
    rule27 ( xnew, ynew, w );
  }
  else if ( mmax == 28 )
  {
    rule28 ( xnew, ynew, w );
  }
  else if ( mmax == 29 )
  {
    rule29 ( xnew, ynew, w );
  }
  else if ( mmax == 30 )
  {
    rule30 ( xnew, ynew, w );
  }
  else if ( mmax == 31 )
  {
    rule31 ( xnew, ynew, w );
  }
  else if ( mmax == 32 )
  {
    rule32 ( xnew, ynew, w );
  }
  else if ( mmax == 33 )
  {
    rule33 ( xnew, ynew, w );
  }
  else if ( mmax == 34 )
  {
    rule34 ( xnew, ynew, w );
  }
  else if ( mmax == 35 )
  {
    rule35 ( xnew, ynew, w );
  }
  else if ( mmax == 36 )
  {
    rule36 ( xnew, ynew, w );
  }
  else if ( mmax == 37 )
  {
    rule37 ( xnew, ynew, w );
  }
  else if ( mmax == 38 )
  {
    rule38 ( xnew, ynew, w );
  }
  else if ( mmax == 39 )
  {
    rule39 ( xnew, ynew, w );
  }
  else if ( mmax == 40 )
  {
    rule40 ( xnew, ynew, w );
  }
  else if ( mmax == 41 )
  {
    rule41 ( xnew, ynew, w );
  }
  else if ( mmax == 42 )
  {
    rule42 ( xnew, ynew, w );
  }
  else if ( mmax == 43 )
  {
    rule43 ( xnew, ynew, w );
  }
  else if ( mmax == 44 )
  {
    rule44 ( xnew, ynew, w );
  }
  else if ( mmax == 45 )
  {
    rule45 ( xnew, ynew, w );
  }
  else if ( mmax == 46 )
  {
    rule46 ( xnew, ynew, w );
  }
  else if ( mmax == 47 )
  {
    rule47 ( xnew, ynew, w );
  }
  else if ( mmax == 48 )
  {
    rule48 ( xnew, ynew, w );
  }
  else if ( mmax == 49 )
  {
    rule49 ( xnew, ynew, w );
  }
  else if ( mmax == 50 )
  {
    rule50 ( xnew, ynew, w );
  }
  else
  {
    cerr << "\n";
    cerr << "QUAEQUAD0 - Fatal error!\n";
    cerr << "  Illegal input value of MMAX.\n";
    cerr << "  1 <= MMAX <= 50 required.\n";
    exit ( 1 );
  }

  for ( i = 0; i < kk; i++ )
  {
//
//  The lower-left 1/6
//
    iitype = 2;
    nbool2 = quaeinside ( iitype, xnew[i], ynew[i] );
//
//  The lower 1/3
//
    iitype = 1;
    nbool1 = quaeinside ( iitype, xnew[i], ynew[i] );
//
//  The whole triangle
//
    iitype = 0;
    nbool0 = quaeinside ( iitype, xnew[i], ynew[i] );

    if ( nbool2 == 1 )
    {
    }
    else if ( nbool1 == 1 )
    {
      xnew[i] = -xnew[i];
      ynew[i] = ynew[i];
    }
    else if ( nbool0 == 1 )
    {
      x0 = xnew[i];
      y0 = ynew[i];
      quaerotate ( x0, y0, x1, y1 );
      xnew[i] = x1;
      ynew[i] = y1;
    }
    else
    {
      cerr << "\n";
      cerr << "QUAEQUAD0 - Fatal error!\n";
      cerr << "  Point does not lie inside triangle.\n";
      exit ( 1 );
    }
  }

  return;
}
//****************************************************************************80

void quaerotate ( double xin, double yin, double &xout, double &yout )

//****************************************************************************80
//
//  Purpose:
//
//    QUAEROTATE applies a rotation.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double XIN, YIN, the coordinates of the point.
//
//    Output, double &XOUT, &YOUT, the coordinates of the point
//    after rotation.
//
{
  double a11;
  double a12;
  double a21;
  double a22;
  const double r8_pi = 3.141592653589793;
  double theta;
//
//  Initialize the matrix of rotation.
//
  theta = 2.0 * r8_pi / 3.0;
  a11 = cos ( theta );
  a22 = cos ( theta );
  a12 = - sin ( theta );
  a21 = -a12;
//
//  Apply the rotation matrix to the input vector.
//
  xout = a11 * xin + a12 * yin;
  yout = a21 * xin + a22 * yin;

  return;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
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
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
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
  const int i4_huge = 2147483647;
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

double *ref_to_koorn ( double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    REF_TO_KOORN maps points from the reference to Koornwinder's triangle.
//
//  Discussion:
//
//    The reference triangle has vertices:
//
//      ( -1, -1/sqrt(3) )
//      ( +1, -1/sqrt(3) )
//      (  0, +2/sqrt(3) )
//
//    Koornwinder's triangle has vertices:
//
//      ( -1, -1 )
//      ( +1, -1 )
//      ( -1, +1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R[2], the coordinates of a point in the
//    reference triangle.
//
//    Output, double REF_TO_KOORN[2], the coordinates of the point in
//    the Koornwinder triangle.
//
{
  double a10;
  double a11;
  double a12;
  double a20;
  double a21;
  double a22;
  double *u;

  u = new double[2];

  a10 = -1.0 / 3.0;
  a11 =  1.0;
  a12 = -1.0 / sqrt ( 3.0 );

  a20 = - 1.0 / 3.0;
  a21 =   0.0;
  a22 =   2.0 * sqrt ( 3.0 ) / 3.0;

  u[0] = a10 + r[0] + a12 * r[1];
  u[1] = a20 +        a22 * r[1];

  return u;
}
//****************************************************************************80

double *ref_to_triangle ( double tvert1[], double tvert2[], double tvert3[], 
  double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    REF_TO_TRIANGLE maps points from the reference triangle to a triangle.
//
//  Discussion:
//
//    The reference triangle has vertices:
//
//      ( -1, -1/sqrt(3) )
//      ( +1, -1/sqrt(3) )
//      (  0, +2/sqrt(3) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TVERT1[2], TVERT2[2], TVERT3[2], the coordinates
//    of the vertices of the triangle.  These vertices will be taken
//    to be the images of (0,0), (1,0) and (0,1) respectively.
//
//    Input, double R[2], the coordinates of a point in the
//    reference triangle.
//
//    Output, double T(2), the coordinates of the point in
//    the triangle.
//
{
  int i;
  double *s;
  double *t;
  double rvert1[2];
  double rvert2[2];
  double rvert3[2];

  rvert1[0] = - 1.0;
  rvert1[1] = - 1.0 / sqrt ( 3.0 );
  rvert2[0] = + 1.0;
  rvert2[1] = - 1.0 / sqrt ( 3.0 );
  rvert3[0] =   0.0;
  rvert3[1] =   2.0 / sqrt ( 3.0 );

  s = triangle_to_simplex ( rvert1, rvert2, rvert3, r );

  t = new double[2];

  for ( i = 0; i < 2; i++ )
  {
    t[i] = tvert1[i] * ( 1.0 - s[0] - s[1] ) 
         + tvert2[i] * s[0]
         + tvert3[i] * s[1];
  }

  delete [] s;

  return t;
}
//****************************************************************************80

int rule_compressed_size ( int mmax )

//****************************************************************************80
//
//  Purpose:
//
//    RULE_COMPRESSED_SIZE returns the compressed size of the requested quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int MMAX, the degree of the quadrature (the maximum degree of
//    the polynomials of two variables that are integrated
//    exactly.  1 <= MMAX <= 50.
//
//    Output, int RULE_COMPRESSED_SIZE, the number of nodes in the compressed rule.
//
{
  const int nnodes[50] = {
     1,  1,  2,  2,  3, 
     3,  4,  5,  6,  7, 
     8,  8, 10, 10, 12, 
    13, 13, 15, 17, 18, 
    19, 20, 23, 24, 25, 
    27, 29, 30, 32, 34, 
    37, 39, 41, 44, 46, 
    46, 49, 51, 54, 58, 
    58, 62, 65, 67, 70, 
    73, 75, 78, 82, 84 };
  int npts;

  if ( 1 <= mmax && mmax <= 50 )
  {
    npts = nnodes[mmax-1];
  }
  else
  {
    cerr << "\n";
    cerr << "RULE_COMPRESSED_SIZE - Fatal error!\n";
    cerr << "  Degree MMAX must be between 1 and 50.\n";
    exit ( 1 );
  }

  return npts;
}
//****************************************************************************80

int rule_full_size ( int mmax )

//****************************************************************************80
//
//  Purpose:
//
//    RULE_FULL_SIZE returns the full size of the requested quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int MMAX, the degree of the quadrature (the maximum degree of
//    the polynomials of two variables that are integrated
//    exactly.  1 <= MMAX <= 50.
//
//    Output, int RULE_FULL_SIZE, the number of nodes in the full rule.
//
{
  const int nnodes[50] = {
      1,   3,   6,   6,   7,  12,  15,  16,  19,  25, 
     28,  33,  37,  42,  49,  55,  60,  67,  73,  79, 
     87,  96, 103, 112, 120, 130, 141, 150, 159, 171, 
    181, 193, 204, 214, 228, 243, 252, 267, 282, 295, 
    309, 324, 339, 354, 370, 385, 399, 423, 435, 453 };
  int npts;

  if ( 1 <= mmax && mmax <= 50 )
  {
    npts = nnodes[mmax-1];
  }
  else
  {
    cerr << "\n";
    cerr << "RULE_FULL_SIZE - Fatal error!\n";
    cerr << "  Degree MMAX must be between 1 and 50.\n";
    exit ( 1 );
  }

  return npts;
}
//****************************************************************************80

void rule01 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE01 returns the rule of degree 1.
//
//  Discussion:
//
//    Order 1 (1 pts)
//    1/6 data for 2-th order quadrature with 1 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 1;
  static double xs[1] = { 
       0.00000000000000000000000000000000 };
  static double ys[1] = { 
       0.00000000000000000000000000000000 };
  static double ws[1] = { 
       0.21934566882541541013653648363283 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule02 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE02 returns the rule of degree 2.
//
//  Discussion:
//
//    Order 2 (3 pts)
//    1/6 data for 2-th order quadrature with 1 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 1;
  static double xs[1] = { 
       0.00000000000000000000000000000000 };
  static double ys[1] = { 
       0.57735026918962576450914878050196 };
  static double ws[1] = { 
       0.21934566882541541013653648363283 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule03 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE03 returns the rule of degree 3.
//
//  Discussion:
//
//    Order 3 (6 pts)
//    1/6 data for 4-th order quadrature with 2 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 2;

  static double xs[2] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[2] = { 
      -0.39011034927098671317419951524681, 
       0.83747122925185483484415621413487 };
  static double ws[2] = { 
       0.14699335257362380894667320570121, 
       0.72352316251791601189863277931619E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule04 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE04 returns the rule of degree 4.
//
//  Discussion:
//
//    Order 4 (6 pts)
//    1/6 data for 4-th order quadrature with 2 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 2;
  static double xs[2] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[2] = { 
      -0.39011034927098671317419951524681, 
       0.83747122925185483484415621413487 };
  static double ws[2] = { 
       0.14699335257362380894667320570121, 
       0.72352316251791601189863277931619E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule05 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE05 returns the rule of degree 5.
//
//  Discussion:
//
//    Order 5 (7 pts)
//    1/6 data for 5-th order quadrature with 3 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 3;
  static double xs[3] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[3] = { 
       0.80383378476840441740523498549521, 
      -0.47391934523147540911429282520838, 
       0.00000000000000000000000000000000 };
  static double ws[3] = { 
       0.82872641363789568276918148034109E-01, 
       0.87120251975907374578897626781336E-01, 
       0.49352775485718467280720708817387E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule06 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE06 returns the rule of degree 6.
//
//  Discussion:
//
//    Order 6 (12 pts)
//    1/6 data for 6-th order quadrature with 3 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 3;
  static double xs[3] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.69739024379082289479659372931424 };
  static double ys[3] = { 
       0.39457278141889907274877084543840, 
      -0.50854615859082585320530997387282, 
      -0.54379745836573696334891625304668 };
  static double ws[3] = { 
       0.11274353612785067755726196993061, 
       0.53124044525363740527637434536885E-01, 
       0.53478088172200992051637079165337E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule07 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE07 returns the rule of degree 7.
//
//  Discussion:
//
//    Order 7 (15 pts)
//    1/6 data for 7-th order quadrature with 4 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 4;
  static double xs[4] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.43435077013834314810201324540358, 
       0.00000000000000000000000000000000 };
  static double ys[4] = { 
      -0.48449728984184861759282345604996, 
       0.95448364011528521109785290893117, 
      -0.49599375367952699580224212194439, 
       0.31755324913853224106767424614158 };
  static double ws[4] = { 
       0.34994956344512530405766865129365E-01, 
       0.26925670356590145219766443119660E-01, 
       0.73377101909709748747894757835629E-01, 
       0.84047940214602985763108417548179E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule08 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE08 returns the rule of degree 8.
//
//  Discussion:
//
//    Order 8 (16 pts)
//    1/6 data for 8-th order quadrature with 5 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 5;
  static double xs[5] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.46537956332076616781921459289142 };
  static double ys[5] = { 
       0.56383112390345027361169431475509, 
      -0.43633565854637050936665616521528, 
       0.00000000000000000000000000000000, 
       0.97959980312548768236416519395659, 
      -0.56281008819734772609629802285497 };
  static double ws[5] = { 
       0.67920849523015500098940229347144E-01, 
       0.62573814354178010612492306599248E-01, 
       0.31655003488030481682328160339460E-01, 
       0.21358892610685618050556600653540E-01, 
       0.35837108849505799692219186693442E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule09 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE09 returns the rule of degree 9.
//
//  Discussion:
//
//    Order 9 (19 pts)
//    1/6 data for 9-th order quadrature with 6 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 6;
  static double xs[6] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.51923560962373232501497734583023, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[6] = { 
      -0.54160946728182000656873212814730, 
       0.00000000000000000000000000000000, 
       0.50274436666672432312125361006830, 
      -0.51354426784066472011973483470318, 
      -0.35942222147133163342323862565198, 
       0.99975295878520206908663626641606 };
  static double ws[6] = { 
       0.20619392336297146785884881219859E-01, 
       0.21306316202539810241408079330466E-01, 
       0.52411159696263022262999116690283E-01, 
       0.56964341363056457388652374913231E-01, 
       0.51213402104188971492967950562937E-01, 
       0.16831057123070001964624080916057E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule10 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE10 returns the rule of degree 10.
//
//  Discussion:
//
//    Order 10 (25 pts)
//    1/6 data for   10-th order quadrature with   7 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 7;
  static double xs[7] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.69780686931593427582366555189730, 
       0.00000000000000000000000000000000, 
      -0.30903100009613455447142535490429, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[7] = { 
      -0.56063064349133316993278274030104, 
       0.10883996591237330573849553975247E+01, 
      -0.51720719429149531106975531479846, 
       0.00000000000000000000000000000000, 
      -0.51225507594767383123122215407598, 
       0.51562570796758001502147445483555, 
      -0.32874839651011214404597359718045 };
  static double ws[7] = { 
       0.64438869372269120316756176676570E-02, 
       0.42018026730627469477523487256964E-02, 
       0.38116505989607355372802824839739E-01, 
       0.18340560543312397707554978722912E-01, 
       0.50983455788600478135692116331631E-01, 
       0.51743930451848519084272818465418E-01, 
       0.49515526441757000856785778879781E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule11 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE11 returns the rule of degree 11.
//
//  Discussion:
//
//    Order 11 (28 pts)
//    1/6 data for   11-th order quadrature with   8 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 8;
  static double xs[8] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.66702609775505741425256941082402, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.33107350040092300374508907880435 };
  static double ys[8] = { 
       0.10478437573860252650592594105291E+01, 
      -0.57312463741413046743394091946546, 
      -0.55246647968395743746723714345900, 
       0.00000000000000000000000000000000, 
       0.76253712102917906246868919741961, 
      -0.35791680916635253698118135523597, 
       0.41170804295590891080211946164601, 
      -0.49479368349849501640395735410577 };
  static double ws[8] = { 
       0.80604906968824787275184319645950E-02, 
       0.82027549569428763887170584691844E-02, 
       0.19158909765241473284488780746507E-01, 
       0.17864637545398712988860649788563E-01, 
       0.26406526528755833416145171420307E-01, 
       0.41518760800101172427427091335633E-01, 
       0.44644591603719541882313987831736E-01, 
       0.53488996928373321021065312076308E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule12 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE12 returns the rule of degree 12.
//
//  Discussion:
//
//    Order 12 (33 pts)
//    1/6 data for   12-th order quadrature with   8 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 8;
  static double xs[8] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.37279552304503872030654041088416, 
      -0.72405807527665067283016173568468, 
      -0.39365448416805093945053087420319, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[8] = { 
       0.21432682937950204373207766849052, 
       0.77622032111803989421823022245158, 
      -0.36989093457992084889029420576697, 
      -0.37591965438942697512621254516818, 
      -0.54031470967359184488712033281188, 
      -0.53745394007281752834331765730183, 
      -0.53648686378750904365565925207548, 
       0.10693230309921692958873012567969E+01 };
  static double ws[8] = { 
       0.41154432712824560996161389102816E-01, 
       0.18744876429730660409284009808805E-01, 
       0.32848111684339866575373457398258E-01, 
       0.56890409960601999633893263927781E-01, 
       0.19851236078200935116627045158968E-01, 
       0.28668810178252124180808712582749E-01, 
       0.15968477487762470187363091359494E-01, 
       0.52193142937027930370255142939621E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule13 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE13 returns the rule of degree 13.
//
//  Discussion:
//
//    Order 13 (37 pts)
//    1/6 data for 13-th order quadrature with 10 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 10;
  static double xs[10] = { 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.39685461846296817096661395897363, 
      0.00000000000000000000000000000000, 
     -0.36877346231328110106712038942974, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.72443410422579609569939260421083, 
      0.00000000000000000000000000000000 };
  static double ys[10] = { 
     -0.56396461592102123502624022053391, 
     -0.47207168193213434618577598142010, 
      0.35411102701218903723318259570554, 
     -0.54446208086261457427222068671521, 
      0.00000000000000000000000000000000, 
     -0.40806649765315498179968814806897, 
     -0.28109188226279360944191830212271, 
      0.76131746182322280252086625652347, 
     -0.53930344496737094791136532996013, 
      0.10684585018901699613124809207464E+01 };
  static double ws[10] = { 
      0.65418593445945714253573279005246E-02, 
      0.21571270093488444532084979604155E-01, 
      0.30310770119514528295026716361971E-01, 
      0.23854497740070562467907639645715E-01, 
      0.11323203959116968208057898102191E-01, 
      0.48973694128817658616454407740346E-01, 
      0.30892926213314122881388117756484E-01, 
      0.20335382082811117457514058128897E-01, 
      0.20258422938614600267531787728451E-01, 
      0.52836422050728359852135506640993E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule14 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE14 returns the rule of degree 14.
//
//  Discussion:
//
//    Order 14 (42 pts)
//    1/6 data for 14-th order quadrature with 10 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 10;
  static double xs[10] = { 
       0.00000000000000000000000000000000, 
      -0.38860728567183008438779563596070, 
       0.00000000000000000000000000000000, 
      -0.23336083105033817175364101627343, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.59834186695364090421374367255377, 
      -0.76078267367321428411638967504470 };
  static double ys[10] = { 
      -0.29206320844176911462814202555440, 
      -0.55198100751730853202008670921318, 
       0.94061946354883534183314930438733, 
      -0.41641460445461153990144627293196, 
       0.20734659086072252500578277182211, 
       0.54084256733761408785791946134147, 
       0.10875282781985526175124163597842E+01, 
      -0.53912013325044374969918965513440, 
      -0.47840728699646114957331822352154, 
      -0.57515345557308018596067002721924 };
  static double ws[10] = { 
       0.21575950013461064401435712745522E-01, 
       0.18999249951197107599659818693116E-01, 
       0.94979085230770220958819598950727E-02, 
       0.50762962987167304291061144573638E-01, 
       0.34069276742966482732002423831004E-01, 
       0.27744543677779980654066228474803E-01, 
       0.32397817681977165704923716215682E-02, 
       0.14400206375318340262255497721148E-01, 
       0.32461956812954507173175490962641E-01, 
       0.65938319732958843565058351143197E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule15 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE15 returns the rule of degree 15.
//
//  Discussion:
//
//    Order 15 (49 pts)
//    1/6 data for 15-th order quadrature with 12 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 12;
  static double xs[12] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.55076221170342555936584321474523, 
      -0.68357214208317707300355657115991, 
      -0.25612692724233233019551175197799, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.36565038512802113430475383248411, 
       0.00000000000000000000000000000000, 
      -0.92281600612376773678051247314954 };
  static double ys[12] = { 
       0.70444274213533005064746095394944, 
       0.00000000000000000000000000000000, 
      -0.43905276476834389171283475266302, 
      -0.43082877955573503227054242095029, 
      -0.54959053538220730475284121644503, 
      -0.40821474151886640425221018734748, 
      -0.54854909315008697510886410926313, 
       0.38728999882555268441925180574268, 
      -0.22031826248214061224887524925114, 
      -0.54538656727512647368330462485919, 
       0.95952641028823274876646693877390, 
      -0.57542156951901391210911429306642 };
  static double ws[12] = { 
       0.48678314316748716941104750027452E-02, 
       0.65212388041010424283650376293385E-02, 
       0.14209708983278779348806400026132E-01, 
       0.31888484893082254935417221561584E-01, 
       0.14777542712078757890433156238493E-01, 
       0.40897290108008783478410013863721E-01, 
       0.10418223735073013293430828867470E-01, 
       0.30458747186574088002369116260931E-01, 
       0.30490829969029447439964565952835E-01, 
       0.21631995447453411916722584709018E-01, 
       0.99261422781568007641581796089061E-02, 
       0.32576332769041589443489039116581E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule16 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE16 returns the rule of degree 16.
//
//  Discussion:
//
//    Order 16 (55 pts)
//    1/6 data for 16-th order quadrature with 13 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 13;
  static double xs[13] = { 
       -0.16279607394216944053888791563128, 
       -0.36133516018585422888297687826246, 
        0.00000000000000000000000000000000, 
       -0.80996884917848457646874211083899, 
       -0.30011609466308557056235128866156, 
        0.00000000000000000000000000000000, 
        0.00000000000000000000000000000000, 
       -0.60911897435540735910318766737015, 
        0.00000000000000000000000000000000, 
       -0.56820549744094988938273289446974, 
        0.00000000000000000000000000000000, 
        0.00000000000000000000000000000000, 
        0.00000000000000000000000000000000 };
  static double ys[13] = { 
       -0.56061007710964974635124810013984, 
       -0.52485883552713864174217026883955, 
        0.92373339140338472305543117791951, 
       -0.55862165096055508852660603004778, 
       -0.39262157635321675549452096003032, 
        0.31873771449384928263441293513356, 
       -0.27527401254502863816148429232788, 
       -0.48840198582588110096280723102788, 
        0.63486450609449659829895244998069, 
       -0.57094697658781472615514574003317, 
       -0.47186155418111734634796824747685, 
        0.00000000000000000000000000000000, 
        0.10956666024303233991933608013003E+01 };
  static double ws[13] = { 
        0.10768394677601291424698896830942E-01, 
        0.18403461944093881198776232104775E-01, 
        0.81763898630962949174826141110659E-02, 
        0.75698865940531301551762664596891E-02, 
        0.41648559391949793671705497882676E-01, 
        0.27100623100161070449309097892289E-01, 
        0.26969791338286996589221138520990E-01, 
        0.23232761214636807134379750252095E-01, 
        0.18940606006195797095777111782639E-01, 
        0.60732744287626264990846070830020E-02, 
        0.17828637150989055090637459591942E-01, 
        0.10139891906267619230864844708608E-01, 
        0.24933912093210466794229664121214E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule17 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE17 returns the rule of degree 17.
//
//  Discussion:
//
//    Order 17 (60 pts)
//    1/6 data for 17-th order quadrature with 13 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 13;
  static double xs[13] = { 
       0.00000000000000000000000000000000, 
      -0.84341388249881453210290659667990, 
       0.00000000000000000000000000000000, 
      -0.15581940864945499700831493674102, 
      -0.44328038905528761799845288018270, 
      -0.24405663711918945323690771490419, 
       0.00000000000000000000000000000000, 
      -0.32008743864026572053341280894294, 
      -0.58451263244119866196418566874803, 
       0.00000000000000000000000000000000, 
      -0.66559778269229534571098852849491, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[13] = { 
      -0.29018817691328482198192482970493, 
      -0.55730147641876997108891034072276, 
       0.52992169651771413561569476668012, 
      -0.55443580380164528313259557057725, 
      -0.55459827350347093500389453629695, 
      -0.30454277894986876852327724535219, 
       0.16498418183331297460767706117349, 
      -0.46069772486520133352908482893009, 
      -0.44217697018376683763624710581761, 
       0.92380408942408538210867783233523, 
      -0.54960689880096757714733164550461, 
       0.11035860158850820902562773489079E+01, 
      -0.45817780070044727282559523202757 };
  static double ws[13] = { 
       0.17971600336645011505239938032544E-01, 
       0.60333417978448673307205126288851E-02, 
       0.17314684664654665424326078868145E-01, 
       0.13685116601127277244126428127471E-01, 
       0.11439597619776484984723997813429E-01, 
       0.34443796770210477782263821006565E-01, 
       0.24818679791573586362196055551041E-01, 
       0.29595573057886602948646133219328E-01, 
       0.27055715740469881665013596685831E-01, 
       0.81984835916342231806329716122509E-02, 
       0.10500033568557456872117406813490E-01, 
       0.18253206778903200987505636654074E-02, 
       0.16463724607144554737778979608446E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule18 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE18 returns the rule of degree 18.
//
//  Discussion:
//
//    Order 18 (67 pts)
//    1/6 data for 18-th order quadrature with 15 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 15;
  static double xs[15] = { 
      -0.13948489081933204088461180178994, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.89294878385120099482738372845036, 
       0.00000000000000000000000000000000, 
      -0.34186434923523414986986916455812, 
      -0.47648266163227498302335864700291, 
       0.00000000000000000000000000000000, 
      -0.57691929083982134519411643537248, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.31378647708219846613556828441497, 
      -0.16617683797609618890338301653712, 
      -0.72011081467827938658866318769494, 
       0.00000000000000000000000000000000 };
  static double ys[15] = { 
      -0.42072604091782665082295117027375, 
      -0.49046440452670875509062398243114, 
       0.62940934145837335576901873813141, 
      -0.55570148308783861437134955410860, 
       0.00000000000000000000000000000000, 
      -0.48379919937747237041994289927636, 
      -0.55915504286023151615285736575981, 
      -0.26927767315911711440578206048817, 
      -0.46282281236309048295631522905634, 
       0.23458453920186008350135824349789, 
       0.11416791732161437494841337947029E+01, 
      -0.31915880712448193443612017530713, 
      -0.55709943481991491343285204410962, 
      -0.55252734012256597819639469738327, 
       0.90376550142496547719360163488612 };
  static double ws[15] = { 
       0.20173122273677479015897043873638E-01, 
       0.86249091344656329557281505991774E-02, 
       0.13370218870435472743618008520977E-01, 
       0.55505642264323721863858554667335E-02, 
       0.67445549565863583016393815247161E-02, 
       0.21538746762008261153428035301757E-01, 
       0.10173035336419548391103790096248E-01, 
       0.22025810771933006688763669533528E-01, 
       0.22256988236841747085623568927686E-01, 
       0.20475740472311754752083054066923E-01, 
       0.35007938360486604839471917528428E-03, 
       0.36314280849967402171385914715287E-01, 
       0.12616049305635046682659603780667E-01, 
       0.10057049329246148211330983846490E-01, 
       0.90745189158503137484947042037209E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule19 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE19 returns the rule of degree 19.
//
//  Discussion:
//
//    Order 19 (73 pts)
//    1/6 data for 19-th order quadrature with 17 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 17;
  static double xs[17] = { 
      -0.71015029250539571205725188857831, 
       0.00000000000000000000000000000000, 
      -0.87005513863591847539846627491196, 
      -0.69945621064432203805136450119586, 
      -0.24805042378404735475371987811852, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.48257092691741797006456094245205, 
       0.00000000000000000000000000000000, 
      -0.25105983153559812324481383767489, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.50128815307877630097134071409913, 
      -0.26869006747347756381758551053812 };
  static double ys[17] = { 
      -0.56868110833624004481988611063840, 
       0.97274416749947380836266322101212, 
      -0.56041590202911501559258672039248, 
      -0.50955355803618757246691504082806, 
      -0.35337391260199717622188423519488, 
       0.76863314856603408380627827197299, 
       0.11143817650139625458794243310849E+01, 
       0.27079298080151709690311725774220, 
      -0.44820650104172157467171108779575, 
      -0.24469161409484972069659993886505, 
      -0.50652963110043023623110954119541, 
       0.53749806844809532876819461539117, 
      -0.43599548606828934030698875208946, 
       0.00000000000000000000000000000000, 
      -0.55141263467657578807862062425378, 
      -0.55150176836813067274425858670454, 
      -0.57376647619684422448809207050633 };
  static double ws[17] = { 
       0.38504278531892868609886548778606E-02, 
       0.46782440974053020444480928944566E-02, 
       0.43790899840937610846456372822265E-02, 
       0.12760020705410626922975774421129E-01, 
       0.34673634319836242706583072061103E-01, 
       0.10025165180245639919338696393626E-01, 
       0.11615619347983174558257617949020E-02, 
       0.20894553379853770048532357879375E-01, 
       0.23831566393070936109285638209472E-01, 
       0.20752749075081240697538867154710E-01, 
       0.21191699422660145703833059227548E-01, 
       0.16221915816210101079240533113508E-01, 
       0.15124040243920516850764130389157E-01, 
       0.75606611406926425730513610639911E-02, 
       0.67921804525194115915872354653919E-02, 
       0.11128622936210414622021183322419E-01, 
       0.43195358902170538658764280819587E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule20 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE20 returns the rule of degree 20.
//
//  Discussion:
//
//    Order 20 (79 pts)
//    1/6 data for 20-th order quadrature with 18 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 18;
  static double xs[18] = { 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.86696389117550812661190351318543, 
     -0.46255868049954115157854976204512, 
      0.00000000000000000000000000000000, 
     -0.67416180418116902260753583861473, 
     -0.22447168033665606945391317277106, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.55640337063475932672618703937903, 
      0.00000000000000000000000000000000, 
     -0.76173172264834809386157300085862, 
     -0.15012093407474928012860293037678, 
     -0.27874288623783821051052255615910, 
     -0.42809996429665844299450237790398, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000 };
  static double ys[18] = { 
      0.50935573580030290124442500215952, 
      0.10254518566344431484070667698816E+01, 
     -0.49506265376045875810405496904400, 
     -0.56894127058564452950719739998646, 
     -0.39335935346809758365690922062188, 
     -0.38873359764810130895317005413822, 
     -0.56423729270253943031731484496373, 
     -0.33519558519158101515864503661923, 
      0.27281208605145075547218518170954, 
      0.00000000000000000000000000000000, 
     -0.49670535155060417023204843137510, 
     -0.20816484443009696237324462486882, 
     -0.51090241799312035370196362374051, 
     -0.56032152802942467258347810157577, 
     -0.48210916153383866244702005141785, 
     -0.55875287099133282930781054909746, 
      0.11166780705147990599713407600049E+01, 
      0.77578464434062421187747427697944 };
  static double ws[18] = { 
      0.12072956229196140184720818378917E-01, 
      0.28443984028101927395806965629882E-02, 
      0.93465277263442986124162626755397E-02, 
      0.29739840427656477081802036308844E-02, 
      0.20327046933776882823607497409423E-01, 
      0.12440057912161099228836551604472E-01, 
      0.57983520915299366578479709959397E-02, 
      0.30774405447413411002155380888934E-01, 
      0.18534535260005961010495899016617E-01, 
      0.61022450704916040183649592581025E-02, 
      0.15757087201875994349887096272192E-01, 
      0.18146095122192897098014042659411E-01, 
      0.10912126413380354963231665469691E-01, 
      0.97275807652705556020653743902143E-02, 
      0.22813420666829580712213934309259E-01, 
      0.94183526939491428417826897656981E-02, 
      0.10513336056091899919123590224643E-02, 
      0.10305163239812520591223081322085E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule21 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE21 returns the rule of degree 21.
//
//  Discussion:
//
//    Order 21 (87 pts)
//    1/6 data for 21-th order quadrature with 19 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 19;
  static double xs[19] = { 
       0.00000000000000000000000000000000, 
      -0.21632544850764896966019269484948, 
       0.00000000000000000000000000000000, 
      -0.51732142577254394519688365999304, 
       0.00000000000000000000000000000000, 
      -0.23848997800485075950680771153334, 
      -0.49726138663597772473398444337837, 
      -0.72412950658754887280604342699955, 
      -0.25333184867896101624130032052433, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.46634273354592047897176529761299, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.69288931850735639864308358232905, 
       0.00000000000000000000000000000000, 
      -0.26783854808923110943600296777461, 
       0.00000000000000000000000000000000, 
      -0.87825896621183219776432484250545 };
  static double ys[19] = { 
       0.11915504280142053470339944811982, 
      -0.22176792985285491406799950307586, 
      -0.56698524576800618707986473265874, 
      -0.56534402375875828330832096864561, 
      -0.24347179507408359903990857698094, 
      -0.36295805875690062325727639767329, 
      -0.50981524806986605698152573307478, 
      -0.56083300395461387753780504376756, 
      -0.48546117771023431955108992405636, 
       0.74251201445344514748786620220395, 
       0.49552108269676082383881281439860, 
      -0.40335187439348537560379590829971, 
      -0.51360341626925635313747432894532, 
      -0.40349668011940382629268981404988, 
      -0.49169851113649364618692601028995, 
       0.96892916731392889311013905422251, 
      -0.55958871884460214682886011537049, 
       0.11174875776997437235819388085953E+01, 
      -0.55950684866453675407246048332841 };
  static double ws[19] = { 
       0.14115632059803205472519666814894E-01, 
       0.23025262548389042223555118865658E-01, 
       0.29202561691086176915750601020522E-02, 
       0.55355656065861078854118874771594E-02, 
       0.15135314836891179901866111095885E-01, 
       0.24278255412699264374036359934954E-01, 
       0.13779168816408564476674209048439E-01, 
       0.58970817066310603066620906137064E-02, 
       0.19083475799818969202994305068534E-01, 
       0.89861747152080313144225770269461E-02, 
       0.12802269114319752046734087041489E-01, 
       0.20930889409056329487663117084041E-01, 
       0.80375338997864963118800061347869E-02, 
       0.12907050562519806085002873135871E-01, 
       0.12913281131480700366519063232496E-01, 
       0.47063366701120542614104893201110E-02, 
       0.90017947131145907457116462987061E-02, 
       0.99277995286896342401230585292059E-03, 
       0.42975457006126745578855094841825E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule22 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE22 returns the rule of degree 22.
//
//  Discussion:
//
//    Order 22 (96 pts)
//    1/6 data for 22-th order quadrature with 20 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 20;
  static double xs[20] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.85243937884353040180383114561341, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.77445005331884859316935899829027, 
      -0.13904123725801872358993428337677, 
      -0.23304071399817585574874678659693, 
      -0.26744719790035475459155300021176, 
      -0.22878505379994225637259252841697, 
       0.00000000000000000000000000000000, 
      -0.38014474895820048454509648055133, 
      -0.45773619140744569728808634071512, 
      -0.69910784887798561027380170772980, 
      -0.59915020623635852909092766037930, 
      -0.50285479867282866077430913574272, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.96233289290517364585870888247669 };
  static double ys[20] = { 
      -0.17961779550825402665053150536356, 
      -0.43105921890016050589097298460315, 
      -0.56370814820709365122621895231987, 
       0.13432079752143621239079354899726, 
       0.50168092900112124608440590878512, 
      -0.30708760414931325796548069627474, 
      -0.49983703894353862614316233613125, 
      -0.51105561801610634168636959252531, 
      -0.39938711325112050072985537594951, 
      -0.56453267538284060850226818603628, 
      -0.24638775273505893357364918660516, 
      -0.56405220111847180239219000037482, 
      -0.50115450425960692839217225970244, 
      -0.38909376761177045913072941229458, 
      -0.56151116446763500432595204778315, 
      -0.49377028890768631789742714088225, 
      -0.56147791277581875684420297893483, 
       0.10538658381143172910281642949327E+01, 
       0.75483396039626819990343725694118, 
      -0.57423523645313343793449840203234 };
  static double ws[20] = { 
       0.88789485269040415184442000093073E-02, 
       0.91213138484983310796492167899650E-02, 
       0.34157891281528518722433813845768E-02, 
       0.13868632623761798604503963871564E-01, 
       0.10542607716173479246541995276511E-01, 
       0.12406032586949329015430983875961E-01, 
       0.98936888063393894436023103307196E-02, 
       0.14736493275333048326268550862834E-01, 
       0.23319638474620094676086399806944E-01, 
       0.64543696619767795008356091295137E-02, 
       0.28567254691249129006535922249255E-01, 
       0.34805812400404550009501818114158E-02, 
       0.15348348448970428911081050216648E-01, 
       0.20675736766822369708334848331883E-01, 
       0.54047041342285612785845222377716E-02, 
       0.13902459659667285747407760303210E-01, 
       0.66515392643746496117585264091036E-02, 
       0.23486059112870622145372854661800E-02, 
       0.94860727130382766679283796797328E-02, 
       0.84285134702804870581139558973652E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule23 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE23 returns the rule of degree 23.
//
//  Discussion:
//
//    Order 23 (103 pts)
//    1/6 data for 23-th order quadrature with 23 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 23;
  static double xs[23] = { 
      -0.65712214849613188354910369404174, 
      -0.76660745759442637356707284950207, 
       0.00000000000000000000000000000000, 
      -0.77617921446782377203133884585557, 
      -0.37522786404062973902300110350405, 
      -0.58004024629858734710994966045544, 
      -0.89831420119030494322671376626199, 
      -0.21445385867353683733519417564443, 
       0.00000000000000000000000000000000, 
      -0.41730612660269037167833945853218, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.59578381850403663485092059957580, 
      -0.21256122928680759054743806564239, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.18972387666403954451075595735093, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.40941795217926264839064707861320, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[23] = { 
      -0.53600577707072872244385605785502, 
      -0.56836123421761655516206642728715, 
       0.10195753956759042232454257963319E+01, 
      -0.52064114314269073816779046417287, 
      -0.57311090188249385547279923552875, 
      -0.56223762431711984631227409503690, 
      -0.56494438615384257743957669387858, 
      -0.45865810284942809229483887487305, 
      -0.50920750140222330666237718651679, 
      -0.40115163422567045202527649710931, 
       0.85387432302293394609857524739803, 
      -0.21127615568211163281540431623842, 
       0.23237891808812345451062467122205, 
      -0.47628237617657111987911856906806, 
      -0.30917657348348663047591112852387, 
       0.67967040631012229365709649358046, 
      -0.57374563316209184490726633368261, 
      -0.55178698374998971587781023092642, 
      -0.38575926863496184163982523488357, 
       0.46621101303269434423430967614670, 
      -0.52020349027005747473328130414494, 
       0.11234666733002454051771461134785E+01, 
       0.00000000000000000000000000000000 };
  static double ws[23] = { 
       0.33272536459172081243423482336248E-02, 
       0.29282906445970122313947765632691E-02, 
       0.25767019981925537648569711133608E-02, 
       0.70120823907155191278558743427503E-02, 
       0.30021012691706992785886911937141E-02, 
       0.54153159980733361962405203321805E-02, 
       0.25697547045533841474189815315258E-02, 
       0.19716254021883959296530794843124E-01, 
       0.75002329339313453980338046569360E-02, 
       0.21216747175005973181144357747315E-01, 
       0.58959569777444837996575869478525E-02, 
       0.15578768482573967095704720947014E-01, 
       0.15666454825087261920671374524318E-01, 
       0.13779632479555575444967581228253E-01, 
       0.27432767705683395855549801326363E-01, 
       0.95806564943421797133558151503143E-02, 
       0.15842534442935738142299813918141E-02, 
       0.93412022737984016156870664924120E-02, 
       0.12471084885337250207089191625388E-01, 
       0.13118150579496436724006564496384E-01, 
       0.13391809372550109677210592618043E-01, 
       0.70104711646684561385279668637298E-03, 
       0.55391494064449379081462896405052E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule24 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE24 returns the rule of degree 24.
//
//  Discussion:
//
//    Order 24 (112 pts)
//    1/6 data for 24-th order quadrature with 24 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 24;
  static double xs[24] = { 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.34667319739036827052339346563820, 
     -0.17228858629081128185759902607603, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.77534695461921828552255298716799, 
     -0.11830295887885374209920762360797, 
     -0.63345645296256978775220168640372, 
     -0.45325106659208302305278781769339, 
     -0.23690281718952127487435067945971, 
     -0.34221988285107932085443071544544, 
     -0.22606667538238729829778396624439, 
     -0.44495329180516622777269400394083, 
      0.00000000000000000000000000000000, 
     -0.80452275594748579655555201318815, 
      0.00000000000000000000000000000000, 
     -0.54364142669910333813418299173702, 
      0.00000000000000000000000000000000, 
     -0.64253738555354065050641405694698, 
      0.00000000000000000000000000000000, 
     -0.91734214502720875776048093000479, 
      0.00000000000000000000000000000000 };
  static double ys[24] = { 
      0.00000000000000000000000000000000, 
     -0.29638036437519615283072347461721, 
      0.59226680488565760168026231308527, 
     -0.28226548001723516079013789395034, 
     -0.28331767767631024291018810948384, 
      0.10127221547590450105942364741960E+01, 
      0.11313827320245876948911869224423E+01, 
     -0.51098115521172256268950056659125, 
     -0.41686453554510700393115126619431, 
     -0.50600932341283035759161511854566, 
     -0.50881145127983229474141858425667, 
     -0.51050575168543342224214849423552, 
     -0.41361091237025662602034980078080, 
     -0.56455388712044611888415763422386, 
     -0.56428020831842462418864564259509, 
     -0.56437817053789745666612468626219, 
     -0.56481964637042001394869703493070, 
      0.23937641179955802312096963113378, 
     -0.41182463487840787556899445189961, 
     -0.51019021318165292778033437868107, 
     -0.56351478264510604651378681421503, 
      0.82100884119022006576877128851229, 
     -0.56336411329583398015443759191052, 
     -0.14556014125612880230802044190807 };
  static double ws[24] = { 
      0.27518427300594242408167672948018E-02, 
      0.86272156924574954117622248939788E-02, 
      0.68297766559224770708029514748518E-02, 
      0.18615927197937334900149347500352E-01, 
      0.20102296969937806928616412243917E-01, 
      0.25227162946694301850912768013948E-02, 
      0.40617630703098676858295809545372E-03, 
      0.70624104072789811560052517607139E-02, 
      0.19783032876181461967862187704075E-01, 
      0.94811737732738431203731642290050E-02, 
      0.11719477113696286147199954244934E-01, 
      0.13091319693878567354029144297298E-01, 
      0.18888756936770716280013255717068E-01, 
      0.55829824091059095164775149600824E-02, 
      0.53712600687036400647412622640952E-02, 
      0.28580170714449744209873764706737E-02, 
      0.34075951304887990048003455619043E-02, 
      0.13502925079066951153738813940684E-01, 
      0.15587004356588926954110113389875E-01, 
      0.68123246685728385281460756274493E-02, 
      0.48789848045465228958508578200054E-02, 
      0.65983957168252596123510923456370E-02, 
      0.23649160198439673642978641242145E-02, 
      0.12499140851132809089730270870368E-01 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule25 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE25 returns the rule of degree 25.
//
//  Discussion:
//
//    Order 25 (120 pts)
//    1/6 data for 25-th order quadrature with 25 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 25;
  static double xs[25] = { 
      0.0000000000000000000000000000000000, 
      0.0000000000000000000000000000000000, 
     -0.1173472784070389885314628951496541, 
      0.0000000000000000000000000000000000, 
     -0.6450240460256727525710559747446056, 
      0.0000000000000000000000000000000000, 
     -0.5664343384798589665228955002024030, 
      0.0000000000000000000000000000000000, 
     -0.3910189989265832790334726502623686, 
     -0.2012181263167401097079073382399585, 
     -0.2038066062734541936977294977052958, 
     -0.7946055588153782250820470421122720, 
     -0.1972942521563999602645401846428416, 
      0.0000000000000000000000000000000000, 
     -0.6355041410448193873418685582169463, 
     -0.2629485129580234408600628587611446, 
     -0.7847262192370077797696941981117842, 
     -0.3957476856294985173592517287113179, 
      0.0000000000000000000000000000000000, 
     -0.4622691203344359376572172601919067, 
      0.0000000000000000000000000000000000, 
     -0.9106735335147376498998363730530587, 
      0.0000000000000000000000000000000000, 
      0.0000000000000000000000000000000000, 
     -0.4402761583839805030337834115210799 };
    static double ys[25] = { 
     -0.1881308452404751474791770362222539, 
      0.4237594812020736794258174198266686, 
     -0.5741998997668708578748305284336055, 
      0.1172287234795072247760841771164101, 
     -0.5133334261186957616800512329664083, 
      0.1025756540328684422937411013290944E+01, 
     -0.4407640888171437954387566505919770, 
      0.6520273733414071456865438634979216, 
     -0.4581028502386777788723229469606161, 
     -0.3764327086872663044277454612227560, 
     -0.4936620181125184332135804336351929, 
     -0.5650036663011180322305141381076653, 
     -0.2245437989161187688779317294631657, 
     -0.3167088545245427646354993582346312, 
     -0.5648168693903931904010878260331929, 
     -0.5549827672775567477688079143971306, 
     -0.5120728253301845715218718158194338, 
     -0.3400477140537724278182036849058156, 
     -0.4464373961278773726299668315862868, 
     -0.5348456380201123672101389788730482, 
      0.8327133265910420576580234570528976, 
     -0.5648988533934863981864107760207050, 
      0.1127558109594723670980720926142653E+01, 
     -0.5393815319272739058039649892539370, 
     -0.5758062076985843575874824269557718 };
    static double ws[25] = { 
      0.9008428931929272978284555982082135E-02, 
      0.7624848013076671474981932002250592E-02, 
      0.2204184225038700283242061056342968E-02, 
      0.1185627435111220996327302583177423E-01, 
      0.8306372211706409316335445471098323E-02, 
      0.2235547623030763175546447504392690E-02, 
      0.1252247261761103193500179756735477E-01, 
      0.7561849278446502447261825323255299E-02, 
      0.1432466758007433215998170125431021E-01, 
      0.2084707600211569027917629357552928E-01, 
      0.1400325214573775947713571136961477E-01, 
      0.3349765842027190344172740201222808E-02, 
      0.2357591414900995147064953536612743E-01, 
      0.1047023089196828804219791222697079E-01, 
      0.4295322980586838056654945868959874E-02, 
      0.7178707806144740205045671854952271E-02, 
      0.6939081726476119924677787292839119E-02, 
      0.1808296563518180275001928371748139E-01, 
      0.8985018370003820265092102638348177E-02, 
      0.9626213990805757770654565513598145E-02, 
      0.6042636212818408481919357354311341E-02, 
      0.2228091765228852260396032167229650E-02, 
      0.5307136158214250141007056133594507E-03, 
      0.5556521038676869804672708399129014E-02, 
      0.1989511820786002256062338480298376E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule26 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE26 returns the rule of degree 26.
//
//  Discussion:
//
//    Order 26 (130 pts)
//    1/6 data for 26-th order quadrature with 27 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 27;
  static double xs[27] = { 
     -0.83506202914393021027297979613591, 
     -0.90755758065010263994368282125220, 
      0.00000000000000000000000000000000, 
     -0.82303054094670403650886493052722, 
     -0.92810181341106024720340709798891, 
      0.00000000000000000000000000000000, 
     -0.69330289802874450653761692118936, 
     -0.68258007506466821571212258497636, 
      0.00000000000000000000000000000000, 
     -0.37884097121369810169657547549398, 
     -0.19241793339756862919367432495484, 
     -0.19560294198083243085696241887367, 
      0.00000000000000000000000000000000, 
     -0.54994973953459751935852183658111, 
      0.00000000000000000000000000000000, 
     -0.16853064985401657599005512286299, 
     -0.19069066073587361426189189987171, 
     -0.70452513950898761906235841668280, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.38125784859592138782794589237690, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.54412206780215007658051815811596, 
      0.00000000000000000000000000000000, 
     -0.36472969593802653513797826614735, 
     -0.54666639745238929451134867455648 };
  static double ys[27] = { 
     -0.56904567277510180660027224756771, 
     -0.52685198805474673472457890788666, 
      0.92351636427246864597020638848735, 
     -0.53195431906859035819325095761412, 
     -0.56748083945007940043218005705387, 
      0.11327377305988868881940987951170E+01, 
     -0.50508092946854773292053923643775, 
     -0.40406610454040422598363913049004, 
     -0.55571014035627985057841942377748, 
     -0.36843999644558691382643352539671, 
     -0.52618905199929260038296168881663, 
     -0.42600632803845996222493394451907, 
     -0.19209885829653007293204006471099, 
     -0.44517417854698789674258323696018, 
      0.20849122488063464188861120733336, 
     -0.57378652101945731348551401873832, 
     -0.28207244982265046632145543744881, 
     -0.56349454153453596561203489840544, 
     -0.47976154988972600729665936044323, 
      0.62053111439688169923700812019236, 
     -0.48872837635331578543021077756047, 
      0.42016147390516239342109813663991, 
     -0.35559725708670166364713905486076, 
     -0.53788618565697879751607695617128, 
      0.00000000000000000000000000000000, 
     -0.56094203589043577139583078831749, 
     -0.57654658422239277276250600074448 };
  static double ws[27] = { 
      0.18405643148504851636845096658812E-02, 
      0.15866124696197450917004000667221E-02, 
      0.32334788927109915113560216754175E-02, 
      0.43503414953892562281646670251570E-02, 
      0.14288712197897272562630508508321E-02, 
      0.34675465145466701985841712380871E-03, 
      0.84276087852092015548095244572855E-02, 
      0.60726432878990522834955233622811E-02, 
      0.34890169558965309347553212683802E-02, 
      0.18924451034788801438603304001786E-01, 
      0.10870464987444917886855259943024E-01, 
      0.18067009059088050381787430564450E-01, 
      0.12810709081456309849270481375969E-01, 
      0.13684071075521919559699025257164E-01, 
      0.12855178620950701539792686120650E-01, 
      0.24441435248010210095439152942177E-02, 
      0.23161806576119345400252126809257E-01, 
      0.39044779268459560333810079419447E-02, 
      0.75861820209001453910764490554566E-02, 
      0.87224512475867079263341905985161E-02, 
      0.13301723810749135614147082883679E-01, 
      0.11150006112175140167331504007601E-01, 
      0.10799966961615883203916038065002E-01, 
      0.82509126176461880950437831298571E-02, 
      0.44936607076337887632843933697535E-02, 
      0.60428306726093991162278907047231E-02, 
      0.14997307146623417159024790140175E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule27 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE27 returns the rule of degree 29.
//
//  Discussion:
//
//    Order 27 (141 pts)
//    1/6 data for 27-th order quadrature with 29 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 29;
  static double xs[29] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.39518500208573395738563214277924, 
      -0.18067167650541349528945895328146, 
       0.00000000000000000000000000000000, 
      -0.22010619980742642366080422865661, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.15697790099606730648797798168280, 
      -0.11950430106482160834133342082078, 
      -0.53186222626517032602651279993926, 
      -0.33699460349041562232377307811535, 
       0.00000000000000000000000000000000, 
      -0.65833539230075168484436817404932, 
       0.00000000000000000000000000000000, 
      -0.34054463966780382785672888075568, 
      -0.48855088636612117128392523236728, 
      -0.11439817909230903220469204061727, 
       0.00000000000000000000000000000000, 
      -0.32048121511116152782316941977578, 
      -0.82912386550847366110475798016042, 
      -0.92517649077195119501691703892095, 
      -0.80474142977810041392358745728250, 
       0.00000000000000000000000000000000, 
      -0.68622113488411601370378338983970, 
      -0.53322822093527312322581065671644, 
      -0.70419862931679748024077449725206 };
  static double ys[29] = { 
      -0.16413151730008831921378959939565, 
      -0.39260212186133936931290882404374, 
      -0.28685548406622891421303264493955, 
       0.87651707827225350568041079993505, 
      -0.52412330045467293432896313035733, 
      -0.35365133466807523825318937197025, 
       0.34617781266034228041185039673984, 
      -0.52879484695937722926887109757979, 
       0.11143276407223406542419307680993, 
       0.54924497463551990458204709806308, 
      -0.21512488834550352002334718330313, 
      -0.46296751685726527586026447389785, 
      -0.50628322547494511263210318063682, 
      -0.56817102280382668923519285087070, 
      -0.52734615007802219988019679750528, 
      -0.46810290449456450390041493630006, 
       0.10418694827155926471229736250245E+01, 
      -0.31933904292428472681005197948760, 
      -0.41333011108690526666610282415552, 
      -0.56768418377653115453110237898598, 
       0.71278197104667799412847982912227, 
      -0.44731325970043847901381172777520, 
      -0.56525617363432326597137998165819, 
      -0.56679589738285794930315523889142, 
      -0.51666886935781835232132779079591, 
       0.11317016084912388514794935172145E+01, 
      -0.54383161967864243507848420517516, 
      -0.56465005418958368460329453488175, 
      -0.57650099541315185890062370844174 };
  static double ws[29] = { 
       0.62908896936458875949825510480729E-02, 
       0.61922333914763423137964574348267E-02, 
       0.79294953183311154387033613859802E-02, 
       0.34300980930711509534751164464939E-02, 
       0.72802514257052478410465398503156E-02, 
       0.16526515457441278585611916710974E-01, 
       0.88646240578297716926902784284285E-02, 
       0.84164942766149210573921235710888E-02, 
       0.10362744260858240831238545743647E-01, 
       0.74242647178535630111397120532648E-02, 
       0.18050472608152491107489365033227E-01, 
       0.12979477412320056876744738018169E-01, 
       0.84957485240956053915702977121675E-02, 
       0.38532361892010405972617343674236E-02, 
       0.46834056015487726028700708666340E-02, 
       0.93844438279733984004320414434716E-02, 
       0.18275921896011569137278968805752E-02, 
       0.16250438567218251268133507322890E-01, 
       0.14073701448288708561589196775091E-01, 
       0.42673273430936345003913524443316E-02, 
       0.64114157223758639239090891131632E-02, 
       0.14385493205073417350192096372167E-01, 
       0.26520511807255348074799784674523E-02, 
       0.15750426454989104155554651654759E-02, 
       0.56948602402767961466213505376975E-02, 
       0.37866239802690712287683217368769E-03, 
       0.60834035556216461895312643871598E-02, 
       0.44670891390107651705739274528023E-02, 
       0.11141963344849334695096764261554E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule28 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE28 returns the rule of degree 28.
//
//  Discussion:
//
//    Order 28 (150 pts)
//    1/6 data for 28-th order quadrature with 30 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 30;
  static double xs[30] = { 
      0.00000000000000000000000000000000, 
     -0.88746481351385107345127069412778, 
     -0.52414142110826588794486451384092, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.33813652586499035853666003249694, 
      0.00000000000000000000000000000000, 
     -0.15580617966892316052928987225377, 
     -0.71154528411036221771832562317249, 
     -0.56152923713272450305409461113685, 
     -0.15254044569676028218764760621125, 
      0.00000000000000000000000000000000, 
     -0.32755930579423371309923889584162, 
     -0.62891280297171951404385772231420, 
     -0.42562113948455467103037050486861, 
     -0.94140501712803148246323405238881, 
     -0.40187196592587673308438468971732, 
     -0.53523254211119856990013167468684, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.34633469607556379390839833888773, 
     -0.70700191257838956126529249352041, 
     -0.22662440063670192933975301253849, 
     -0.10878662096604677283105449627330, 
     -0.20631360388374033702271668616836, 
     -0.79729933059714082413378391932539, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.85821826344224735007618585707670 };
  static double ys[30] = { 
      0.10167280551560398976220492653467, 
     -0.54006894012803566318386100608171, 
     -0.49235910875172712775847678059386, 
      0.11380585570620829080112594647327E+01, 
     -0.43282761808513505745779087778619, 
     -0.26963680688236239861351493942337, 
     -0.18337097274071869661622880573497, 
     -0.24858562418623951123781187881298, 
     -0.56962812500573696889109864019987, 
     -0.42548755921389883340768034668690, 
     -0.46791873168797857582730644925957, 
      0.26004754891078197570252274618257, 
     -0.57016318648548524036109080252581, 
     -0.53786850890931757394897545264776, 
     -0.53803125336295947108113127503277, 
     -0.56669876546503991682732553682320, 
     -0.38330233101247314078613426299384, 
     -0.56906848596694110780690954139249, 
      0.78786657352369156249728778058495, 
     -0.33331193876680684460752457704886, 
     -0.46918587486573210609633873244633, 
     -0.49038194407928164918552182401052, 
     -0.53278806818989164900646131251018, 
     -0.56757013000065183447611373867027, 
     -0.37281400062453798641014279266735, 
     -0.54575369430987430575525395306512, 
     -0.52483841201603278034315141768996, 
      0.60516346416700599387343286263650, 
      0.94394739374781944379137982654772, 
     -0.57527136550578091340619030932143 };
  static double ws[30] = { 
      0.94510343300931478888807464349302E-02, 
      0.27868387832627458914989081520856E-02, 
      0.77591410167297521729327740749191E-02, 
      0.20473848133013403509762806085984E-03, 
      0.58247494675789234908196566576895E-02, 
      0.16526278580754572911455183578579E-01, 
      0.93510917286426097352060797088729E-02, 
      0.17248616209516818765814273146624E-01, 
      0.23411580330424100035578129796391E-02, 
      0.10544291927832240582270995892779E-01, 
      0.11140166083431255028726297685594E-01, 
      0.86155130881102002892400222873137E-02, 
      0.31279976876179532211383324626933E-02, 
      0.56287726620134695271354834891188E-02, 
      0.67531408157501848324491455305183E-02, 
      0.12483494047479294769456871723105E-02, 
      0.13484476650919357880017855769673E-01, 
      0.30592737461003165508862947478784E-02, 
      0.51392317566586590413074635030704E-02, 
      0.89148192557630731697759322343474E-02, 
      0.11746113414800820107755108695172E-01, 
      0.77883730444682187225229982545828E-02, 
      0.82393269507997219649806927417388E-02, 
      0.41327061926490009393858952877356E-02, 
      0.16828820578478445853343577421317E-01, 
      0.42179069123934048296232732334268E-02, 
      0.50551130485581726523795360821379E-02, 
      0.78587185794036915888389747222587E-02, 
      0.33745795976141645086914914308322E-02, 
      0.95433079635401447385836219413628E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule29 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE29 returns the rule of degree 29.
//
//  Discussion:
//
//    Order 29 (159 pts)
//    1/6 data for 29-th order quadrature with 32 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 32;
  static double xs[32] = { 
      -0.87938703907162049393738648870781, 
      -0.14324628011414286412542273880746, 
      -0.35189538626807214171525397778246, 
       0.00000000000000000000000000000000, 
      -0.61554248268722133288024470174652, 
      -0.93007010487368394982733169286299, 
      -0.19742931116675328671024827830067, 
      -0.82766874673327517194764945031097, 
      -0.72491944440609982649676161027856, 
      -0.97606406728761785676875297522146, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.48236271621499169294961905466403, 
      -0.17405270877030766542159522724513, 
      -0.18672840783491580658483096931098, 
      -0.49956457122986497525117031568907, 
      -0.64840539960610806885812016702526, 
       0.00000000000000000000000000000000, 
      -0.76445927822389460704185764718864, 
       0.00000000000000000000000000000000, 
      -0.34535027682270612143829088272072, 
      -0.51325936500558699339036907718902, 
      -0.35532987585108812184220678363790, 
       0.00000000000000000000000000000000, 
      -0.32573169227336622410380909993533, 
       0.00000000000000000000000000000000, 
      -0.18040262369805029614059002911543, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.62723509783699348733500748068739 };
  static double ys[32] = { 
      -0.57262394724341608181908512435579, 
      -0.30511050720283478498316524226896, 
      -0.57371128821577701886180134696590, 
      -0.57359111146432687487737878042109, 
      -0.45928362498519066569504746311567, 
      -0.55859049859349923530121259720141, 
      -0.19814832945992239613993592077499, 
      -0.54048028336522896166963698554618, 
      -0.50660081473204674367166958460713, 
      -0.57457293178871935440458908351462, 
      -0.35003739955028665529791783308077, 
       0.10123351916554935691973791840677E+01, 
      -0.40146746347211072112690881875070, 
      -0.50542328918299929386212475640451, 
      -0.41479896168463444735644949768926, 
      -0.56230562425237497147426436905503, 
      -0.54688501268470743019578735701145, 
       0.43276338426682519796086695407376, 
      -0.56778324409762155852431972748668, 
       0.59786046049173946734694450915043, 
      -0.53602539551087235345778481154042, 
      -0.50755670245009570927755668717151, 
      -0.45978142418549716875276682521809, 
      -0.53717224317544601306267373599653, 
      -0.33566750808701799356359087057475, 
       0.10720328289321576370831152259048, 
      -0.56337853290640048165947671108094, 
       0.75831443067019621851290191888360, 
      -0.45528374978288560170739457877104, 
       0.89932044998481280953329700019771, 
      -0.19844398974763577810899514361613, 
      -0.57713645303871158784974277557031 };
  static double ws[32] = { 
       0.10123938451463759970089359052465E-02, 
       0.13755888140419727618158412742168E-01, 
       0.17804164437108562499941904599297E-02, 
       0.99788818279604048507240176161182E-03, 
       0.84743804866337155687860067687316E-02, 
       0.16743020381843061998817864135626E-02, 
       0.17826547012452998087762815394247E-01, 
       0.36343788555051997590886515586824E-02, 
       0.59390298859384334961673129484654E-02, 
       0.35467172304180440218981679487014E-03, 
       0.73509973972781686569789083069102E-02, 
       0.18823092677622066092420624317723E-02, 
       0.12353105880035723559433384108907E-01, 
       0.10295539688423531976513712860782E-01, 
       0.13990388400987746999801020264183E-01, 
       0.40156561227163462310099803253767E-02, 
       0.50344656666784366169434955781421E-02, 
       0.82512598969804302641269300260534E-02, 
       0.22921998323995127511965450098240E-02, 
       0.69563841616376041067990591845970E-02, 
       0.74610255726590545841053915581726E-02, 
       0.86623962103173048554261254624063E-02, 
       0.12080143560322868032076843683152E-01, 
       0.40366629789560296898640227002180E-02, 
       0.16475096833731066604740732689647E-01, 
       0.10732740719960036973478529643535E-01, 
       0.47922734876003179932483846983130E-02, 
       0.53780564637593511319450846071907E-02, 
       0.67864123981058792365600157845442E-02, 
       0.36908289770711881028317444421851E-02, 
       0.10471484548953909836949668724134E-01, 
       0.90634414524923745915451079527243E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule30 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE30 returns the rule of degree 30.
//
//  Discussion:
//
//    Order 30 (171 pts)
//    1/6 data for 30-th order quadrature with 34 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 34;
  static double xs[34] = { 
      -0.43405710781709278714048879994740, 
       0.00000000000000000000000000000000, 
      -0.13719609648588667050981701694674, 
      -0.23002603213365866154430597150794, 
       0.00000000000000000000000000000000, 
      -0.35671356772393844272156480720059, 
       0.00000000000000000000000000000000, 
      -0.49401409080195115810253005052280, 
      -0.37062552685345952318116499093046, 
      -0.51430581999606514802110718905564, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.54607380091083807654188300887860, 
      -0.84275872174416357830668461312136, 
       0.00000000000000000000000000000000, 
      -0.57724024570193215242188127722161, 
      -0.39821033011703771207712265483060, 
      -0.30605985322806922570625655304903, 
      -0.68776182141386689567302952570524, 
      -0.19684758015531914348305332007873, 
      -0.96321844715106605579144337343699E-01, 
      -0.28286197991210025994286433769979, 
       0.00000000000000000000000000000000, 
      -0.10333424594611078538852107076057, 
      -0.64379410387887037331359309442254, 
       0.00000000000000000000000000000000, 
      -0.89873358431890446012484882993312, 
      -0.78478618004244445633333163168023, 
      -0.72349108151471496801951044440950, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.20239462350406048608674549751121, 
       0.00000000000000000000000000000000, 
      -0.94049820614395764193267866092540 };
  static double ys[34] = { 
      -0.49449740552349402097090260436227, 
       0.11432041379660233658011039218632E+01, 
      -0.43937592123719222447451521408216, 
      -0.47742236692181668651784560731673, 
       0.90399516562083738727059247641795, 
      -0.44352953384065130754596083984676, 
       0.99134074534282262769308436316314, 
      -0.53793159129250810586076585776398, 
      -0.36290163137317815520540870562751, 
      -0.37663681425566764191708385069014, 
      -0.46660354209881582130470470379699, 
       0.11107528506358136609304671750972E+01, 
      -0.46213867523774957532303180078394, 
      -0.56966207432086504530851133019367, 
       0.73349519835786388112994533410054, 
      -0.56927271273635792416671462689250, 
      -0.56920325339091635438895039029322, 
      -0.53373365021893670540255969319794, 
      -0.46342190137736602703775671709799, 
      -0.35889692861698105392115484048066, 
      -0.23982388127083061470496201687854, 
      -0.24652609751239687685958293412409, 
       0.52281528171007302347717037989158, 
      -0.52966101087647732119998666908062, 
      -0.52874213045823877050213253058208, 
      -0.10029941117116486849674700598946, 
      -0.54980647624872101250977315077756, 
      -0.53007527494669219480580008827439, 
      -0.56749280180954759083517414571717, 
      -0.35822245431797926788870429663662, 
       0.21098212229083163080456906246771, 
      -0.56840882187053609024013730511891, 
      -0.56806418709095800771552792081267, 
      -0.57642585867286344650886243834268 };
  static double ws[34] = { 
       0.55469340852267197383357734501536E-02, 
       0.11300463022066637667549395901160E-03, 
       0.77976368535667454572847114089957E-02, 
       0.78230337119887396150629609702565E-02, 
       0.25018123919483940113671864156512E-02, 
       0.90508836547186347263240388076065E-02, 
       0.18717426158053290329737464178677E-02, 
       0.53860315402766333154068593578253E-02, 
       0.11677070021275261169640320523198E-01, 
       0.99157266127682177040571302280702E-02, 
       0.51222447690487278622934837026924E-02, 
       0.55706132346425325140873159634701E-03, 
       0.89545764756588636881737594179016E-02, 
       0.15881019953772277632781829119128E-02, 
       0.48608613509509536172532545379762E-02, 
       0.25781213160755392530469248853581E-02, 
       0.29933391441706776190700124360794E-02, 
       0.70603223137458408522017573340119E-02, 
       0.77254463220316942248431581666276E-02, 
       0.14875615552650235020592865135842E-01, 
       0.18841554220226669798776264579558E-01, 
       0.17011603553515205522964322169417E-01, 
       0.70757738603786203845878772981467E-02, 
       0.82670460199627466892622746876773E-02, 
       0.61880970884828196145627428599983E-02, 
       0.10262805248776557218482365083790E-01, 
       0.24090825302303188401346360674076E-02, 
       0.50298776686404106640379168657878E-02, 
       0.25267049096575867286260646250422E-02, 
       0.81629221043416840383932440429541E-02, 
       0.98143263368157515388390143296277E-02, 
       0.34778149711369936288944433476055E-02, 
       0.18367906196371196197648328999357E-02, 
       0.44170301264357154992013311250048E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule31 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE31 returns the rule of degree 31.
//
//  Discussion:
//
//    Order 31 (181 pts)
//    1/6 data for 31-th order quadrature with 37 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 37;
  static double xs[37] = { 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.61514906219173326038112958639798, 
      0.00000000000000000000000000000000, 
     -0.60227258275696561313011448174463, 
     -0.50981188019097632901132686993753, 
      0.00000000000000000000000000000000, 
     -0.93815432009933315159290706133886, 
      0.00000000000000000000000000000000, 
     -0.20647073303633426254589462709978, 
     -0.14400331667261834628330435968823, 
     -0.48161033312898040755998550507708, 
     -0.10037070511088677024119003378018, 
     -0.69003405314433981556933467756463, 
     -0.45131857503676764396920554079559, 
     -0.84660740995916577326254742794424, 
     -0.14911203462776176458914390755779, 
     -0.55167853197329218224541056875306, 
     -0.88689641508961216382151190186048, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.27328678036217323411592305906149, 
     -0.33518301010596487973294137974615, 
     -0.41995122085501584151906684983796, 
     -0.18328764721779652677839208920087, 
     -0.36991861294507902457366355095885, 
     -0.19826213773515198256349732520883, 
      0.00000000000000000000000000000000, 
     -0.32037552704774453797884509286080, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.78037306571931027109013054924730, 
      0.00000000000000000000000000000000, 
     -0.63930924222334969496336881416464, 
      0.00000000000000000000000000000000, 
     -0.73284741377468351220017427799278 };
  static double ys[37] = { 
      0.00000000000000000000000000000000, 
     -0.52938101092158436217989944700065, 
     -0.39259296230588589577692024633618, 
     -0.29879521023518401246203229051833, 
     -0.57170178188865743017272606769990, 
     -0.55459143708784673167882967475209, 
      0.88032014667050929273459147232799, 
     -0.57328995909056364047175904022464, 
      0.10932773428047454311567544236664E+01, 
     -0.57395302409498325748592469605497, 
     -0.31297598217411644944584755018888, 
     -0.38089192422910876315233159687303, 
     -0.48291646123243918459436362644262, 
     -0.47237945269405575937245782085877, 
     -0.51224429424952928367019627443114, 
     -0.56975292247630686483425195330683, 
     -0.54743183795864664399668537119941, 
     -0.46896195562367563745199360999814, 
     -0.54546676966840707030987521015721, 
      0.39801431704076458186329655275394, 
      0.55731348796377460751706837305857, 
      0.97770163996668454891174645987393, 
     -0.50303973149300230512177279576695, 
     -0.55450610597908806787792991531332, 
     -0.57527720939459095316883310890960, 
     -0.19780119278968994222147469745240, 
     -0.43682917497899069302689014620258, 
     -0.41232373405740865931876786653227, 
     -0.56955126226549368079464495424716, 
     -0.32293156387181891591095890622273, 
     -0.17956217803484541051959756627326, 
      0.12415978610574965002707225626963, 
     -0.52966581912912838443677982459313, 
      0.11370931699099784490463357476566E+01, 
     -0.53253326710134984412140457620381, 
     -0.40532624526153890211957544149763, 
     -0.56730901867998811775513569725020 };
  static double ws[37] = { 
      0.18300274150131503530218918260118E-02, 
      0.26610940186702783315047214273191E-02, 
      0.76415116557529348094479438274031E-02, 
      0.58190756304798839792080800346536E-02, 
      0.16896438572357084969785073782623E-02, 
      0.37199019857877275953974147517983E-02, 
      0.27793705473843310754776111206512E-02, 
      0.66528731990096257950775283571449E-03, 
      0.64939444386637490706930367932606E-03, 
      0.17144436271854393525611545784889E-02, 
      0.13414188931651306592837821888099E-01, 
      0.10353383822703655875877567178813E-01, 
      0.93715119938232228413834490764269E-02, 
      0.65331420598530105269402018636481E-02, 
      0.71547540812742549665898738315076E-02, 
      0.15168844552172759195105122436401E-02, 
      0.54989950731761515927595152639621E-02, 
      0.83734196392879596672486354136481E-02, 
      0.23857385517888400677977593831317E-02, 
      0.73943641501760668515497495971065E-02, 
      0.62903046982657749327353939866577E-02, 
      0.21445996231349418275012329739199E-02, 
      0.88044815101158239329257671015790E-02, 
      0.47901323321320787218271749114584E-02, 
      0.12231201498467946837374439661990E-02, 
      0.16971546653439500088953398556642E-01, 
      0.11154908721966859054415465858694E-01, 
      0.12663626772734683766437904485010E-01, 
      0.15316868409913918465476865568287E-02, 
      0.15266022971056427587974399198100E-01, 
      0.90410304859438967096315948145332E-02, 
      0.81327060196892968093589471350676E-02, 
      0.48254994785913298058391441904342E-02, 
      0.22662885264281957658336382036483E-03, 
      0.58503938090777195674071616162139E-02, 
      0.68659431853565733752147231623839E-02, 
      0.23969034602009614667762140991332E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule32 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE32 returns the rule of degree 32.
//
//  Discussion:
//
//    Order 32 (193 pts)
//    1/6 data for 32-th order quadrature with 39 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 39;
  static double xs[39] = { 
      0.00000000000000000000000000000000, 
     -0.42746868613914916208067002398662, 
     -0.91710704502811275874385926715291E-01, 
      0.00000000000000000000000000000000, 
     -0.48363645746700960769905958976308, 
     -0.74677040622149170601230763257760, 
     -0.33090034155997483681589019846520, 
     -0.46023597931812041521733845612288, 
     -0.55152849427095291026768125922720, 
     -0.20577133469764029385252414280604, 
     -0.20092303305428395370178245682808, 
      0.00000000000000000000000000000000, 
     -0.25879477956175550900580265426342, 
      0.00000000000000000000000000000000, 
     -0.72262398109805332027279578058226, 
     -0.32997225542915476024532745217312, 
     -0.40104907184189418029309583833701, 
      0.00000000000000000000000000000000, 
     -0.69151287469732085458587108553859, 
     -0.81882407794345790149465089749751, 
      0.00000000000000000000000000000000, 
     -0.15979469018140228096803155997403, 
      0.00000000000000000000000000000000, 
     -0.59877821961446415276093759190966, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.42015173485047987871656798146767, 
     -0.28990563943088228305314863137265, 
     -0.84085334890089953008631910023907, 
     -0.12723619924623591749401870286216, 
     -0.12551628141173309089569225449278, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.55077303550427213809118028350091, 
     -0.61718975394462310684988723266656, 
     -0.90920085094674871327008780126851, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.97068830122109422411832138489096 };
  static double ys[39] = { 
      0.11496517991990804340276319918613E+01, 
     -0.30218292328705052033455030789710, 
     -0.55177404065234470341871573116633, 
     -0.57267682745711810079098827470171, 
     -0.48826894371900250566656288607256, 
     -0.54347180509529604896284843807779, 
     -0.55549427045965430240427224834972, 
     -0.54086954371044273649778573751227, 
     -0.37375893906115265347255148781801, 
     -0.52660317145013860310284351646008, 
     -0.57138712716648846304698735307588, 
     -0.51211470939471719438254705852599, 
     -0.42082435075243723269873387979224, 
     -0.41206328125638506041335006800986, 
     -0.48953938794088206965231133840212, 
     -0.49430799576545804368134930045435, 
     -0.57513903374312331695098163467832, 
     -0.14894052019732994524801638469381, 
     -0.56996729278416821373273663004497, 
     -0.57067752268839624045481729672027, 
      0.17429645803883591285042159190117, 
     -0.22764895370246828181434158209954, 
     -0.27832827065225098989964142553377, 
     -0.46007242792992039354306838595294, 
      0.10498449848071351612896202347394E+01, 
      0.11103247160622636896577670049089E+01, 
     -0.41181398695180711307526515295114, 
     -0.31690029177033716996969306425727, 
     -0.53517509064720514075419136628404, 
     -0.35324390902713134915617189643158, 
     -0.46783716452744596679002275357204, 
      0.00000000000000000000000000000000, 
      0.93988934616034742698090141489668, 
     -0.56822819321283957200082177635155, 
     -0.53204175653167234191507273763322, 
     -0.56772912901378600226993106149821, 
      0.36045825890517508775557027645869, 
      0.80150480143647657259482305879194, 
     -0.57627103176029262025065680355071 };
  static double ws[39] = { 
      0.35089101065112855884820534496463E-04, 
      0.10988190324830709992823329341192E-01, 
      0.40017234188803658556418944416512E-02, 
      0.93087707933178618554865620359439E-03, 
      0.69568246166267704238682724201103E-02, 
      0.32813277551329795738555421297416E-02, 
      0.39753082327573825644781696272502E-02, 
      0.47960047753629003693425154224810E-02, 
      0.98456677416312611477210982996463E-02, 
      0.60120289449877842202024432440538E-02, 
      0.22178527178401095681968551311850E-02, 
      0.36718564604057522445250408427746E-02, 
      0.10803923534388294394797611307667E-01, 
      0.55134227454329365414727793248372E-02, 
      0.58771612493962009794350932546429E-02, 
      0.81633895836356217774190986558429E-02, 
      0.11698483763495028290383571359353E-02, 
      0.83316003706962053904317473673468E-02, 
      0.19546717234154615848859112275635E-02, 
      0.14644039054791717956440739265680E-02, 
      0.84773151751406951264947448004187E-02, 
      0.15757826762525756648322149683424E-01, 
      0.71597492934303095837773594696045E-02, 
      0.77415633062502329911450279233184E-02, 
      0.12234738959380913656810149160128E-02, 
      0.52757528161206686084277006574669E-03, 
      0.11119171725536887499648558701843E-01, 
      0.13439031765939374980048747865869E-01, 
      0.32769705348614632852685205970683E-02, 
      0.13014005319508587872203236381061E-01, 
      0.97509412664434843717560320628702E-02, 
      0.28580896374081882633069151128591E-02, 
      0.26165845261107117861016768586572E-02, 
      0.25605447042338306730739080320977E-02, 
      0.56075763180908099390771878498063E-02, 
      0.13889507817898641123545037024783E-02, 
      0.82488872017764297885488574625852E-02, 
      0.43623899484098009158061336724647E-02, 
      0.22384872276251377786582863606787E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule33 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE33 returns the rule of degree 33.
//
//  Discussion:
//
//    Order 33 (204 pts)
//    1/6 data for 33-th order quadrature with 41 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 41;
  static double xs[41] = { 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.39676550645769949223549079823286, 
     -0.31598316073300286874362023003082, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.64019940087608054570302313946691, 
     -0.13023310833724125410060567348685, 
      0.00000000000000000000000000000000, 
     -0.77591825174786530058425772428980, 
     -0.19029241130244378942621074383698, 
     -0.66275085646437554778774954725515, 
      0.00000000000000000000000000000000, 
     -0.13890838462723977815572667811625, 
     -0.52593106350448959206487843135769, 
     -0.74112180754338096384725945089585, 
     -0.49295570914510750117532576951889, 
     -0.51670323488631757686566331872787, 
     -0.64505375945992157587556971416084, 
     -0.61872928859268166187915712963621, 
     -0.87492200353713827567393201652485, 
     -0.33145619841996054218301251037654, 
     -0.97409390068179509856666677924027E-01, 
      0.00000000000000000000000000000000, 
     -0.86499638747373789749985796834781, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.24972283611883747878887622813179, 
     -0.47779523112848852807050836313123, 
     -0.16616445416618183551509983345254, 
     -0.36279240961074976308711123134325, 
     -0.77375327602359261325494586474574, 
      0.00000000000000000000000000000000, 
     -0.30711355053101709037576993455102, 
      0.00000000000000000000000000000000, 
     -0.18486459332759578217311403693356, 
      0.00000000000000000000000000000000, 
     -0.42184939757610775040328718080002, 
      0.00000000000000000000000000000000, 
     -0.94593887427662142045742644100072 };
  static double ys[41] = { 
      0.85441752281491986952742561507422, 
     -0.26359952191974260301537245052622, 
     -0.42678545506855360887308799127774, 
     -0.37973706727153893912099920904950, 
     -0.15588528269354452279063922983879, 
     -0.57381590818073294683247861206343, 
     -0.41858125784281633932442886389431, 
     -0.20333227220053604130569869278831, 
      0.25520959894249213064365318384422, 
     -0.53856910378372629185609575858799, 
     -0.40545819052596609782356418100840, 
     -0.54194078514648208041737749935800, 
     -0.50291356104093599630234744015635, 
     -0.47583110208346099057285412944826, 
     -0.54370533785038485941253752308528, 
     -0.48924051171923119683779249416411, 
     -0.57075764711609390637957695741143, 
     -0.41461766926988511673184437368405, 
     -0.57099887763461860406270665329696, 
     -0.49001038870971575141844877415585, 
     -0.57022781943877866932429672746875, 
     -0.56983893491575439638788166530905, 
     -0.32834341163271889420656088222861, 
      0.10732337258984243455437606728791E+01, 
     -0.53947756141543200561909733382606, 
      0.94602784584493370468772946259542E-01, 
      0.64724758119173276973587947838225, 
      0.96632545246884377179105403347845, 
     -0.28583342001808206268947748302965, 
     -0.49332159111642935520258261160453, 
     -0.56914512902820977690651047487837, 
     -0.53901095406565886298444238871685, 
     -0.56974068564230494042674463877735, 
      0.39886122381719090545979996295666, 
     -0.48180999095667246889253860020090, 
     -0.54960747258649002726722665090028, 
     -0.53508842266119758145398714640325, 
     -0.41592300096438369677257895846249, 
     -0.31644289654629315370197407113192, 
      0.11383269060703847163981510995473E+01, 
     -0.56955190911504415394155156282000 };
  static double ws[41] = { 
      0.21072740283763039964528368286828E-02, 
      0.44974723243006456842650776075088E-02, 
      0.65946145346910002385096276412626E-02, 
      0.82474035280720032563304370435810E-02, 
      0.58790090916007076226614232725982E-02, 
      0.68507656125260760567769180226135E-03, 
      0.67199953862190751868955790397041E-02, 
      0.12685962295354290443758741902940E-01, 
      0.64085284625279664423720130246453E-02, 
      0.31865284621268815383673784463267E-02, 
      0.10130277996161918616988936043147E-01, 
      0.39327511834173456227281919359377E-02, 
      0.35957285155409189611031640004502E-02, 
      0.81581506936299909146481456622838E-02, 
      0.43654139395527354256349689343937E-02, 
      0.48812258789986847763530046372840E-02, 
      0.20302369941307448210755634209178E-02, 
      0.87021082793685719718873832596530E-02, 
      0.17835148063472755026619148855589E-02, 
      0.64535378321272648894894135877724E-02, 
      0.12052160414072912400513700553737E-02, 
      0.23976450419859318170667605866734E-02, 
      0.12163656006474426058122295662701E-01, 
      0.90899196307377903030326256967572E-03, 
      0.25803358156466235366918919356997E-02, 
      0.79548382691120458192515870842605E-02, 
      0.48546328081496349968870545571536E-02, 
      0.21118713588645547377857629192177E-02, 
      0.13437007473057861497945927056313E-01, 
      0.71907262853635647670060927729590E-02, 
      0.26467590924742060643044789761419E-02, 
      0.56174443807205437983099746393445E-02, 
      0.16949406736417696510033206105949E-02, 
      0.65319228663443823255328913298067E-02, 
      0.86859349355210683062499905616676E-02, 
      0.25415771481263707182156665796150E-02, 
      0.65111877505010114950629655024188E-02, 
      0.61419843433928965634632131426334E-02, 
      0.12083509278790420262130603712492E-01, 
      0.19312249945614427939614114959174E-03, 
      0.84755399951394965389373925159078E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule34 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE34 returns the rule of degree 34.
//
//  Discussion:
//
//    Order 34 (214 pts)
//    1/6 data for 34-th order quadrature with 44 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 44;
  static double xs[44] = { 
     -0.20865073605333605259133081748132, 
     -0.89097175315861241014243393897698, 
      0.00000000000000000000000000000000, 
     -0.87999414063426440899680777462237, 
      0.00000000000000000000000000000000, 
     -0.68124048521147079022287669759749, 
     -0.78844430308137329865810016089767, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.28883149528014893210917834177554, 
     -0.80036260521166327118670792806244, 
      0.00000000000000000000000000000000, 
     -0.14656285308642525770840545989704, 
     -0.14791223778178358225012515887081, 
     -0.14888251010183166048943112172279, 
     -0.14412341550112265052276598150328, 
     -0.56021943491268762662833110324255, 
     -0.56775489079636195195819292013723, 
      0.00000000000000000000000000000000, 
     -0.43438042322624036290148813406327, 
     -0.42883411933471152581183344420389, 
     -0.69612085434081448660482017170361, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.29429125011188233571351981822902, 
      0.00000000000000000000000000000000, 
     -0.69057074957270704503616084061837, 
     -0.14635945730951804290474707357531, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.43610043846910337710191570745380, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.29361459879318225343484932066898, 
     -0.29217524245073863070988162069564, 
      0.00000000000000000000000000000000, 
     -0.29404250943783623341704803057639, 
     -0.95251437348711362491282542183012, 
      0.00000000000000000000000000000000, 
     -0.57112662404759029584585818705102, 
      0.00000000000000000000000000000000, 
     -0.57210738184377731997454761477038, 
     -0.43576534554125907389165005866616, 
     -0.80626144515788714899642746572657 };
  static double ys[44] = { 
     -0.54080960852670874887477002688154, 
     -0.57148188680340885268957968075793, 
      0.99541215623761383488320268143123, 
     -0.54581251858802002963651717370634, 
      0.10825129994552011239408759106716E+01, 
     -0.47254254003026177485543685475103, 
     -0.51513519234555459140932034419392, 
      0.89225266841244135802626470038980, 
     -0.14988699059565575877251167626283, 
     -0.29318823199080331034093749074232, 
     -0.55884631397612334218448641893041, 
      0.15927070953202400003601747376533, 
     -0.52342688512318759090959018821675, 
     -0.34665387857625309237224070668568, 
     -0.44866353740446774321091903079561, 
     -0.22199147096348917763877176210378, 
     -0.42108532754196956310264613688775, 
     -0.50128574606027466276226199416617, 
     -0.48871256531163417297535112276726, 
     -0.45790141180639478630175258882333, 
     -0.36052903462331485421264904911911, 
     -0.56926107479496737145288159890053, 
      0.76958700429572619540848652285881, 
     -0.54916967076805156376085658979603, 
     -0.49260808560575517970684223356176, 
      0.63201013425192545991254006499216, 
     -0.53464092774804677680478041519000, 
     -0.56702734467084155713753902149103, 
      0.11402969514608451589323215865362E+01, 
      0.32180957687924744604856694624021, 
     -0.52795439225981246678322066497151, 
     -0.28467385501300634121071117733737, 
     -0.57576434340793074190688904900431, 
     -0.55039498029219558132681857598373, 
     -0.57574687817495954263227255109469, 
     -0.39917165954426298474560583795016, 
     -0.40538150617520055089614596703473, 
     -0.57046622308344295569063522362823, 
      0.00000000000000000000000000000000, 
     -0.55348323362954892556316690015126, 
      0.48132579880630916313439867481111, 
     -0.57601288500270791096736728894872, 
     -0.56789016986039254576904292632732, 
     -0.57734211176549375189279322347485 };
  static double ws[44] = { 
      0.13382626571418695846413420959349E-03, 
      0.86077553205023877867589308006139E-03, 
      0.15128050789281141606424158029908E-02, 
      0.20001846280203786077172851326386E-02, 
      0.72363050803554618757272999467355E-03, 
      0.61455888057448998470789535942896E-02, 
      0.38852026646137484310391372656345E-02, 
      0.26304196612070312222228001589138E-02, 
      0.77253587933225472077027128073975E-02, 
      0.12977714128349407406048789350468E-01, 
      0.22719033752163729514493246145023E-02, 
      0.77393098675853113647899856957649E-02, 
      0.66279889144916813462020297846060E-02, 
      0.12710485190182604950478151438930E-01, 
      0.10017843547790636952167111558575E-01, 
      0.14563495893338662830933135825511E-01, 
      0.85595344302267414192383009908361E-02, 
      0.63601741070384869392179065287274E-02, 
      0.42668601490321023633055488964966E-02, 
      0.86775530450313316070466948384888E-02, 
      0.10904023948277298869126932531947E-01, 
      0.18417526470910829459939716211994E-02, 
      0.38386152359715808879730672139045E-02, 
      0.24540640916463922042739572558992E-02, 
      0.80218185952419318245374428395049E-02, 
      0.50585349116444063800445633458038E-02, 
      0.41366919289449465339947189819144E-02, 
      0.29365977875069631100768219443565E-02, 
      0.14934897179533907877154395103507E-03, 
      0.71158143270123038007615755398533E-02, 
      0.58377502625954269481520792337924E-02, 
      0.70064127307723735616871240238792E-02, 
      0.50071757656760638993478483434402E-03, 
      0.45908969684050017035075831240866E-02, 
      0.96316333582374746961054859980204E-03, 
      0.58308259716404371642003962574066E-02, 
      0.10858660084749162408673215897113E-01, 
      0.65536218237774572732745631345446E-03, 
      0.26542942460412844806471662947295E-02, 
      0.36926643355893676151865055391569E-02, 
      0.61856406685725288566265475154246E-02, 
      0.76023090931745106769467860453136E-03, 
      0.25829647679686151643244490553125E-02, 
      0.37816775394238441141631154528165E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule35 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE35 returns the rule of degree 35.
//
//  Discussion:
//
//    Order 35 (228 pts)
//    1/6 data for 35-th order quadrature with 46 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 46;
  static double xs[46] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.88784082040153075383913617100652, 
      -0.10034207598327724033582632199321, 
      -0.69311601738832964374229025488167, 
       0.00000000000000000000000000000000, 
      -0.20170550651440610611345276285483, 
       0.00000000000000000000000000000000, 
      -0.64267456637991915772112509507131, 
       0.00000000000000000000000000000000, 
      -0.59680262574112570330708343507959, 
      -0.88583797294857920950727599039044, 
      -0.21481054147146732337946003185254, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.97750969393388806071547935914909E-01, 
      -0.55225246543256484131744271704192, 
      -0.73330302740947963883782960844096, 
      -0.35377855031176987704663393999516, 
      -0.51564990450718314801776121220598, 
      -0.21973987181432191570664654404377, 
      -0.60482829086183809094609158655982, 
      -0.10540651029785986602106161987810, 
       0.00000000000000000000000000000000, 
      -0.44461399827537925801209489543884, 
      -0.23876011572551469292227793760715, 
      -0.10680624500393084714250643212956, 
      -0.49845820778025588899721360580822, 
      -0.47118728299279689577235187734666, 
       0.00000000000000000000000000000000, 
      -0.33861070795850472330848428220997, 
      -0.70498953851611867618296732663118, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.11540856290209161674967217397369, 
       0.00000000000000000000000000000000, 
      -0.35516859452256584748215994566211, 
      -0.95143357699664369331063333170408, 
      -0.80625540528872521995344375173798, 
      -0.39094238101462525709590836324919, 
      -0.26704107522244421375617609033818, 
      -0.80012284035409657436386927229493, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[46] = { 
       0.11049848391267878160942425437432, 
      -0.54839524322647351656910274265446, 
      -0.84731241225870783138119631392495E-01, 
      -0.57376800669530175521786194508562, 
      -0.56192636681814703241337147158792, 
      -0.57368528854276496320410931239880, 
       0.24567191532977353520899699695488, 
      -0.57170157791299617481071719582110, 
      -0.24069657869901295964844222588206, 
      -0.49911718050236027863216490838426, 
      -0.48412771111544869326449479743582, 
      -0.56244420234421444672958208328214, 
      -0.55314597935499411028041061203470, 
      -0.25264039633717009571608632345015, 
       0.93691922814514350972936873734162, 
      -0.37503784513189034216209663360014, 
       0.10867100617799430553354204728332E+01, 
      -0.17426334519355714993668083479551, 
      -0.52711753152174214154822817555781, 
      -0.48456866722910133562449136216106, 
      -0.57129468294377403989555329497135, 
      -0.57531330331362051925464151910002, 
      -0.37997627488267737489770886024621, 
      -0.43023229931186750876181610790903, 
      -0.31309961497003354510619289311518, 
       0.10151427309957352682750877998298E+01, 
      -0.55376202677277905332098019475239, 
      -0.47996330941134033797610668708564, 
      -0.43574907178404019355497831707649, 
      -0.46711834658382799515239504110706, 
      -0.37748161580172233726943066542578, 
       0.81339592685777538075733642892247, 
      -0.32068265475981446715948171913088, 
      -0.54513781159869833529380365664930, 
       0.67133897977061262700017196648536, 
       0.38727728656145265547624996623580, 
      -0.52024278749878884072085026914167, 
       0.52856305553369069237166742739323, 
      -0.42946854793782251649705803709131, 
      -0.57105763188292920761880945578841, 
      -0.53228217266152343324657869863566, 
      -0.50900886837618163733840811909782, 
      -0.54310638162348937349259853010943, 
      -0.56864030976986894926046466151582, 
       0.11398838577722204060015136163662E+01, 
      -0.57704966306400061720563422243113 };
  static double ws[46] = { 
       0.37823533023071117552498633546329E-02, 
       0.15110951850926768597510849257431E-02, 
       0.52612367337077241481150311006990E-02, 
       0.59776000384333233272193666149629E-03, 
       0.25949680545647278336391650188551E-02, 
       0.10300467854022705106865695429020E-02, 
       0.57842676116151487639100767425388E-02, 
       0.16644297264251806023973860743196E-02, 
       0.54751091383596566484521677047910E-02, 
       0.46816125223955569111350042656706E-02, 
       0.34712691077662295851692973049338E-02, 
       0.24253041782253153086402273893358E-02, 
       0.17490359570497127164220966999176E-02, 
       0.11875479166746748231460750630995E-01, 
       0.19531898113644349044152767302173E-02, 
       0.50044915460995835332657889685798E-02, 
       0.66667113158108803350566280388053E-03, 
       0.11325480833009634041680955354389E-01, 
       0.47516431054402191637098512902079E-02, 
       0.46372667587850735462627390787042E-02, 
       0.18136506984725309737491530658395E-02, 
       0.91471037630436898650123223221922E-03, 
       0.10699314011695035989362484260043E-01, 
       0.72247940752599472272602792262149E-02, 
       0.10839168886619908505604562425067E-01, 
       0.13587382545922657958478769707436E-02, 
       0.36583876907265655319114947846825E-02, 
       0.83330935899037325736589164743249E-02, 
       0.91410693800006639894202414819246E-02, 
       0.73251904599782639661305917852163E-02, 
       0.96568380444916204234246412306618E-02, 
       0.32925800341072999195193696158169E-02, 
       0.11464882195010925904138963396459E-01, 
       0.36017648083679220604863221387225E-02, 
       0.44214031660381166110271695499832E-02, 
       0.60288751331023346137831096345316E-02, 
       0.61877293600659226814496412622813E-02, 
       0.55090243819186250387827772801666E-02, 
       0.94428428348687329331168226963113E-02, 
       0.62282804392723232782667478677435E-03, 
       0.33358029363864306610289766382805E-02, 
       0.66727022641263705281781993426708E-02, 
       0.53354589201841413526093533076666E-02, 
       0.16509864052584545067812812415452E-02, 
       0.15814849311712982200162693233010E-03, 
       0.41297372110944178234379022954739E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule36 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE36 returns the rule of degree 36.
//
//  Discussion:
//
//    Order 36 (243 pts)
//    1/6 data for 36-th order quadrature with 46 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 46;
  static double xs[46] = { 
      -0.45770306387539169487915756826196, 
      -0.36361676298891185295386221029935, 
      -0.74962911892600870983032844542448, 
      -0.26737758237748029380439006185063, 
      -0.37317732959222391549799088251508, 
      -0.14975882871475764975924141602265, 
      -0.98352138144106149192683092192778, 
       0.00000000000000000000000000000000, 
      -0.53788829941643881420320023927424, 
      -0.25563046668462188936725827478982, 
       0.00000000000000000000000000000000, 
      -0.83136119811518354383269937195326, 
      -0.13003286298009251819823947174789, 
       0.00000000000000000000000000000000, 
      -0.60110699315720762998855042089700, 
      -0.63205241025983805117322394914783, 
      -0.52716585532238106513416234938719, 
      -0.28230973261381914542579826867683, 
      -0.71047632569890861227787925936018, 
      -0.95056420733325419196436293042879, 
      -0.28616499791046129448997529539718, 
       0.00000000000000000000000000000000, 
      -0.13633349630546001365985579457732, 
      -0.80918371229769335783766611724880E-01, 
      -0.44770259704283524854791791199685, 
      -0.48011047223406774553422252507270, 
      -0.41323708163607651211902687274295, 
      -0.65560184610300882914185718548851, 
      -0.57620923725166014923400127915066, 
      -0.24925193834014909515080159065743, 
      -0.78377703132887513374849375170181, 
      -0.98306369371972070224522934512975E-01, 
       0.00000000000000000000000000000000, 
      -0.15401622396033279817412492899376, 
      -0.35730985892987794548052370157585, 
      -0.72653937769578913812744606839955, 
       0.00000000000000000000000000000000, 
      -0.83027114581134724322172615184608, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.88144405800714427701736717808599, 
       0.00000000000000000000000000000000, 
      -0.91784551568339890614029102074855, 
      -0.17476902776510310686941360192574, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[46] = { 
      -0.52141842729563126012047058298966, 
      -0.41429464277319715623346981701418, 
      -0.46477792693460978350876606581729, 
      -0.44967445714422062537104312954445, 
      -0.49063516603887538827137270218634, 
      -0.46321446291835242797318542339315, 
      -0.57521271319697780857478481745947, 
      -0.27374719021510399523990036642840, 
      -0.39427972611017856706920948711432, 
      -0.36279824244610613498688160674974, 
      -0.36620767514629492498613572227248, 
      -0.57248771387722136932451719517558, 
      -0.38562994925042077844843752485788, 
       0.34045539811756807080574481636674, 
      -0.57257600517555581330832753808181, 
      -0.45255306558366174094217160197888, 
      -0.55235918086652448153563150959453, 
      -0.57260000432656212732742831820917, 
      -0.50975814717158154135559916631247, 
      -0.56440377229715042680745590073970, 
      -0.27325102827656031685861861579114, 
       0.63223872333035094330098051711355, 
      -0.28077682324477788489838234773238, 
      -0.51793028571538908930758218474956, 
      -0.57210278655469173709707059163532, 
      -0.45179917425831580211503767552840, 
      -0.34557238584596557842323790994865, 
      -0.54673762175281423328617405554397, 
      -0.50424349278017526914470934564358, 
      -0.51773945365247080611701750038930, 
      -0.54712080331353357214461326692730, 
      -0.57522359685526201507706056294147, 
      -0.45359100459111562262893708000554, 
      -0.16813042003647684568253959866018, 
      -0.55046376829560201599703664333391, 
      -0.57068525840467165716710741013583, 
       0.49249815594813255797249024870368, 
      -0.51140681795255970991370425926398, 
       0.10530289859103386181163034523952E+01, 
       0.87254187575228190842967598185962E-01, 
      -0.55627077065397527838016283877844, 
       0.76956773175760798113534871535368, 
      -0.57586320130293542669591032679490, 
      -0.55598397134820992601771756543425, 
      -0.16513031625331098943557133545280, 
      -0.55823188110938432057255915241886 };
  static double ws[46] = { 
       0.40294016520725152813048610101848E-02, 
       0.64793903053765145181715901682133E-02, 
       0.43742626283010051684766029946910E-02, 
       0.62417922487631623395376249125259E-02, 
       0.54140669757067043817444691914728E-02, 
       0.68900982432856468561811786486874E-02, 
       0.17153699268662009594629984868070E-03, 
       0.50474692507570505960794202203366E-02, 
       0.68800457554916895631729667842012E-02, 
       0.87270677947457016043671196486654E-02, 
       0.45248570442564970883554024651669E-02, 
       0.91498024223602954229350467916617E-03, 
       0.91050902425828265328119448284342E-02, 
       0.53090410657138158843888963628538E-02, 
       0.13393775372397616562830163994317E-02, 
       0.59631198989805405320442107067512E-02, 
       0.31313092298346458015049957705066E-02, 
       0.16015029731880952140633711164508E-02, 
       0.40969742526751036868296166462955E-02, 
       0.85305987209231664440378714676407E-03, 
       0.10982854456089695449388211470504E-01, 
       0.40607862976056182266227222376222E-02, 
       0.11718063546352663453334801471048E-01, 
       0.62695317676519732186303155045592E-02, 
       0.16229081033223297129772742579429E-02, 
       0.69502345508257539200837869664982E-02, 
       0.94809421287883652563970965623795E-02, 
       0.33560498861862835325040698854943E-02, 
       0.53868503871236664135869203894291E-02, 
       0.58820503274730197860770659304124E-02, 
       0.27653642263897575237831549083688E-02, 
       0.10819301104675578264822260694178E-02, 
       0.43833382315646555911862314735351E-02, 
       0.13441518749900172007698067515263E-01, 
       0.40011563187151066613912579689780E-02, 
       0.14603290876836360969200891041503E-02, 
       0.53256416015449335120955671226641E-02, 
       0.35098393196898440174721922929817E-02, 
       0.10314228067656089530374958038011E-02, 
       0.72438530919178505538177420621044E-02, 
       0.17577599858575697519932889346354E-02, 
       0.35688958985136608389688752740980E-02, 
       0.35030971685040415436169439133341E-03, 
       0.38004098068192230521610341307517E-02, 
       0.69135644753295898893485091923938E-02, 
       0.19056197400002277482559131629861E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule37 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE37 returns the rule of degree 37.
//
//  Discussion:
//
//    Order 37 (252 pts)
//    1/6 data for 37-th order quadrature with 49 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 49;
  static double xs[49] = { 
       0.00000000000000000000000000000000, 
      -0.33412271632399888890129019512936, 
       0.00000000000000000000000000000000, 
      -0.75613518652324771803726882164662, 
      -0.84800448035565502181663502981330, 
      -0.70551641007683071137220723820188, 
      -0.12346882279099691408072860937378, 
      -0.14115293561514993591833853717433, 
       0.00000000000000000000000000000000, 
      -0.31129235394249030634514266109515, 
      -0.25517618371949366971568295794296, 
      -0.22738491727535507077740696347565, 
      -0.82431411808054765778463063897013, 
      -0.96025084028127376878602592691577, 
      -0.88743516782176053029087010855756, 
      -0.14743801672844793582653341492985, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.63852651231947468903625801359977, 
      -0.91244418004082918428945065210440, 
       0.00000000000000000000000000000000, 
      -0.53650720489459350535466965366810, 
      -0.85277881576583532724521330388158E-01, 
      -0.16435656319638607339573924940553, 
      -0.81381272571201343203517786852045, 
       0.00000000000000000000000000000000, 
      -0.41910547578514426267307171420310, 
      -0.48677494045468852012238429335006, 
      -0.35909034523353810535799934366341, 
       0.00000000000000000000000000000000, 
      -0.25248858363296016573959873583852, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.39948855862139271795827453837350, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.74822550833094322230914938071543, 
      -0.61214302766947216204021932607674, 
      -0.82456232718866158786784556422156E-01, 
       0.00000000000000000000000000000000, 
      -0.28849633742385887838307907599387, 
      -0.69338845905777882473767890474803, 
      -0.17221861996121408821061372532618, 
      -0.49166762890981157665024971850791, 
      -0.33772735494586915349190113708600, 
      -0.56199209754321689406198698034801, 
      -0.62908155899883917171575081362878, 
      -0.47085424810000757697255495637523, 
       0.00000000000000000000000000000000 };
  static double ys[49] = { 
       0.80030932175661049330529129464111, 
      -0.51875656535329760898448830295380, 
       0.10333868827205359544088392857986E+01, 
      -0.51993973121116010442502084667746, 
      -0.57382818530981477957536290227481, 
      -0.48120343056361501578369820373486, 
      -0.16501678669005198586553208123293, 
      -0.44697785411625922396526760726341, 
       0.25665474437476898507539344982761, 
      -0.48447626446390567150790800247712, 
      -0.41922309411584432631191867018353, 
      -0.34332292689716552850094371729457, 
      -0.50673434341218750366276527969451, 
      -0.57217095597866797982914123136412, 
      -0.54747072920426972017807350995914, 
      -0.26364744637852948251388994918932, 
       0.67166192797723746626700329929959, 
      -0.44283386491499269166138431840481, 
      -0.42124005904839257204244609933248, 
      -0.57039245293560439673032555344414, 
      -0.15872646837427968615407557344963, 
      -0.43711422090147355697061182546047, 
      -0.54811072534421336521000509462082, 
      -0.50551055464016532277335287814759, 
      -0.55356126813574196138633993389223, 
       0.10946289447927698298356343395322E+01, 
      -0.54631433608883669102312593267134, 
      -0.35483129429226142896168273650094, 
      -0.33950304358022085403206973553427, 
       0.88453895937762153767287377553499, 
      -0.54714013952455212484601864326434, 
      -0.50545261693006535308892323009918, 
      -0.27272204187887408822459840414627, 
      -0.42805065874848260516272058041530, 
      -0.57184433744987283372049694766113, 
       0.83856840257338103751495755746255E-01, 
      -0.57144178004035639723831793597820, 
      -0.50017785048732220313797140417802, 
      -0.36587111104193859470979519681048, 
       0.11419818578211332577350031418259E+01, 
      -0.23836996263996572403671622331600, 
      -0.54563400792324411093065715450446, 
      -0.57169987252525305855247428576888, 
      -0.57133062792022124429581677109025, 
      -0.57122090496927329315283644781684, 
      -0.54470987727952378493452147796822, 
      -0.57104889149424033456087085228856, 
      -0.49810325155934541668515841885833, 
       0.49872349530628412097743604970645 };
  static double ws[49] = { 
       0.99626129730910408577949008400168E-03, 
       0.27626921675168826283774018443100E-02, 
       0.72692273468131503092373352604759E-03, 
       0.25303170257571451860803668676728E-02, 
       0.63229856016814473011800273541530E-03, 
       0.39910569052688676630224931092243E-02, 
       0.96322224231485228364344323861314E-02, 
       0.70914264027495165913238721401754E-02, 
       0.49912714680319196696789614840510E-02, 
       0.60134263305403734057625688074224E-02, 
       0.75284170176253568400370133216087E-02, 
       0.88334951566927715848362109377594E-02, 
       0.30494454960037220959042145422116E-02, 
       0.40623159427664301328639018693128E-03, 
       0.15854780081932780104600463914237E-02, 
       0.10744094290784490004353067282619E-01, 
       0.34547338677747299508167133011287E-02, 
       0.38698068479336796185125547322801E-02, 
       0.57440366524074741575884680707207E-02, 
       0.72357295271785506116849988932337E-03, 
       0.56884381303256696164588543099345E-02, 
       0.69047009212170979910711955177161E-02, 
       0.42884132950362761878776351004925E-02, 
       0.62924591180266519167258563983064E-02, 
       0.21172628941575284289686456117100E-02, 
       0.57748488832703718923029245549870E-03, 
       0.39252782695269835705405080071355E-02, 
       0.80200877092618206400095291824248E-02, 
       0.93274649768244852461773810704269E-02, 
       0.23598744244779275367141637367377E-02, 
       0.41453075810464728064443485645257E-02, 
       0.32896867854980720978623838536080E-02, 
       0.59742143446146166209795010947358E-02, 
       0.82279235776610154350252369267537E-02, 
       0.93156858510830936728471952983330E-03, 
       0.67163446268129400718735416067323E-02, 
       0.12493560653618925032728437197660E-02, 
       0.53054261077983319226673244613556E-02, 
       0.10755243143420365190907771759900E-01, 
       0.11403474399353258205391049010235E-03, 
       0.11790258338603034214469882053681E-01, 
       0.32403041710010828670599208753256E-02, 
       0.18705266785451222748868057200548E-02, 
       0.17176825471147128641632325800221E-02, 
       0.19029239546031546635251984688755E-02, 
       0.38621442944772685484213090942755E-02, 
       0.15746286043986147827363630712349E-02, 
       0.64138925190696845281721401600690E-02, 
       0.54555303295239163064914865711407E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule38 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE38 returns the rule of degree 38.
//
//  Discussion:
//
//    Order 38 (267 pts)
//    1/6 data for 38-th order quadrature with 51 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 51;
  static double xs[51] = { 
       0.00000000000000000000000000000000, 
      -0.38017185669450718257184066048831, 
       0.00000000000000000000000000000000, 
      -0.64388710735877083783416282124044, 
       0.00000000000000000000000000000000, 
      -0.91532058656497198671471150555246, 
      -0.78472586424382273509683032899814, 
      -0.99643590876239824517695137609747E-01, 
      -0.41959606976158875740459488513640, 
      -0.90709796569683108888235443082745, 
       0.00000000000000000000000000000000, 
      -0.31006871991021548748342323634250, 
       0.00000000000000000000000000000000, 
      -0.55804669150858522261961383013874, 
      -0.23391252232540839712384994231650, 
      -0.84191118240738060068532375730213, 
      -0.53335107942756502512226568928935, 
       0.00000000000000000000000000000000, 
      -0.47291121570664026204475560732776, 
       0.00000000000000000000000000000000, 
      -0.72785722160360793350295551572105, 
      -0.49855221184504938895139013840901, 
      -0.54925671022190333222562093502386, 
      -0.44686452030684963518394665288121, 
      -0.61953974615995630087673220975416, 
      -0.62648548681833042070303279649082, 
      -0.67197877643540341416126929672985, 
      -0.41335885239677796031890142990207, 
      -0.25346463551838318570541525347334, 
      -0.69573906885921915728856348088559, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.38636910084523127223009550231297, 
       0.00000000000000000000000000000000, 
      -0.33658056387534557504382625690100, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.87843541547485157244218200006356E-01, 
      -0.85039120745944430746903967678568E-01, 
      -0.84482245043808986751701243019350, 
      -0.23361755931365234413765013776635, 
      -0.84712737676102112301556527859528E-01, 
      -0.30228200141771618391567732367974, 
      -0.17281546621195622411401122869585, 
       0.00000000000000000000000000000000, 
      -0.16293712483088177437982114203995, 
      -0.96366288222453288706533712520860, 
      -0.25373681319593324282811397911950, 
      -0.15001122172665224948595080689653, 
      -0.75738968673027804907685927008861, 
      -0.76664557469367476317397564232357 };
  static double ys[51] = { 
      -0.43557198951312314492372904354947, 
      -0.48676901082212536688741381364887, 
      -0.20761979727194592905202199533091, 
      -0.40573908618816893429385257585425, 
      -0.46713152507187031292455455867277, 
      -0.57377461458362246990957759635438, 
      -0.51373173879402534769265540295112, 
      -0.18600062141874083923383175875705, 
      -0.52071016051529641903708216132188, 
      -0.55471321549486532778867447040816, 
       0.10940631576403911595307511796686E+01, 
      -0.44966761533984648209489215134754, 
       0.65074263823985079600667000386443, 
      -0.57266700781635460274026638478223, 
      -0.19775702200253896797934274544035, 
      -0.54145730613515149557950920114789, 
      -0.40306833625953531748818108675870, 
       0.11441448170815456331207001620839E+01, 
      -0.45699037409328199099206927354961, 
       0.95123251126820600706088723471076, 
      -0.45891721551063656195753394995739, 
      -0.55254653716155075726305443495177, 
      -0.50985596144990674077481289770119, 
      -0.31858800667671763885666636445179, 
      -0.46956784465504062728813836459427, 
      -0.54683284777191952223714461151061, 
      -0.57051260429903260207066265388305, 
      -0.57227331670065037577673308824224, 
      -0.57253866498558213253629629823534, 
      -0.51790249949511235759270582432134, 
       0.10252000964578561221300395045227E+01, 
       0.16568297477093707545432698117652, 
      -0.38903470904529446210903397210864, 
      -0.30137742399876299186434574375181, 
      -0.55145765476506148518201936114546, 
       0.43524327738107859551153511292157, 
      -0.55031736573002371942681152441060, 
      -0.51270307120416046582718931909697, 
      -0.57221773164738784945037090035811, 
      -0.57007049747206373797913896032469, 
      -0.37942157065303255897465810627733, 
      -0.38182236985847584746876547919066, 
      -0.29988148386350240360102081846111, 
      -0.55110081638940231919507254561012, 
      -0.80478911639902275408961840267122E-01, 
      -0.45603365215790654433852111390607, 
      -0.57138683817501246632052994868585, 
      -0.51210804711263968662217366453516, 
      -0.28983256020905742573818686292564, 
      -0.55650931396143459003034201519312, 
      -0.57683373408531432268827479475092 };
  static double ws[51] = { 
       0.20820056888817244409479850097856E-02, 
       0.32704508067578457571759580719843E-02, 
       0.38571880274891361301769303800511E-02, 
       0.48648533235668256195606333219063E-02, 
       0.25862281661138416465531703526730E-02, 
       0.45760472370487254239987376847120E-03, 
       0.27293096561559624215019181897416E-02, 
       0.97008356269688197848677107128175E-02, 
       0.41132552385450103084160052046070E-02, 
       0.11925949669960168440297068927338E-02, 
       0.46027699181420719864281441155474E-03, 
       0.65746158861474390224655188929927E-02, 
       0.34233865355919380271157248680394E-02, 
       0.12182478611229097763044222506603E-02, 
       0.10899944052496863003448798018155E-01, 
       0.20225984312508308593119096478526E-02, 
       0.68735133627418967754698853779575E-02, 
       0.81627294555617067479425317996855E-04, 
       0.62097541317418902117276820215705E-02, 
       0.16492922508404918879311251694587E-02, 
       0.47406525818327189115185327580062E-02, 
       0.31194776872854070583692733919320E-02, 
       0.47516718437350612556088524941901E-02, 
       0.86077222408115445210495576178103E-02, 
       0.53685573833170132885235134226360E-02, 
       0.30607129436300396172490108601521E-02, 
       0.13236335711461000683695750999937E-02, 
       0.15044854999557370818299509578742E-02, 
       0.15641030619522655360987944091216E-02, 
       0.38234831557293443707212103289394E-02, 
       0.11236312245926095551123994880102E-02, 
       0.59128635891188981067020381172796E-02, 
       0.83264254588952817413009549759947E-02, 
       0.52233031594939664397142395297434E-02, 
       0.35872048025524414964720174624147E-02, 
       0.50988347474248667740715631552740E-02, 
       0.19977329387543947701265023377701E-02, 
       0.59507492016878607729260830919901E-02, 
       0.16966487646916691851186225113499E-02, 
       0.10748164815034162115762746127887E-02, 
       0.91437918154557240969683908824785E-02, 
       0.93423065980415739850720220419123E-02, 
       0.10123136480463870914740930355634E-01, 
       0.38498047327371845578787166272780E-02, 
       0.65645729720023034651464119787433E-02, 
       0.79710615214023377342492350020114E-02, 
       0.43918131172356037933829206838480E-03, 
       0.59611163393367252996223811895425E-02, 
       0.11058361133847283290469370750187E-01, 
       0.23833729931365671861610978427443E-02, 
       0.38466956567350313890347038963597E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule39 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE39 returns the rule of degree 39.
//
//  Discussion:
//
//    Order 39 (282 pts)
//    1/6 data for 39-th order quadrature with 54 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 54;
  static double xs[54] = { 
     -0.51907510039262467181240723005295, 
     -0.27085312542591789026981871557739, 
     -0.62535084249163554801261088188681E-01, 
     -0.44171982481354024257677278787610, 
      0.00000000000000000000000000000000, 
     -0.15015342269982087499746548552200, 
      0.00000000000000000000000000000000, 
     -0.73242413817247197416912819947887E-01, 
     -0.61576635195673085768179561363787, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.48003862478268162304965734842567, 
     -0.36508484113357361074797983209786, 
      0.00000000000000000000000000000000, 
     -0.56713187575703966770973366858629, 
     -0.14031613163580384457468777152879, 
     -0.25349195822799364362686756031925, 
     -0.23248267893550702450813771192708, 
     -0.52163635975054715145448041789043, 
      0.00000000000000000000000000000000, 
     -0.90690106293315290819021535664027, 
     -0.12932342290720507011192674546810, 
     -0.13784410422442617726158008200667, 
     -0.66294103217865996758685557932486, 
     -0.64122202646501331223433885902278, 
     -0.81530008183876068918941172976357, 
     -0.90238671492070698682084691210471, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.36521604921152885859291119415780, 
     -0.81466971749275774599927221958440E-01, 
     -0.50100475560027347803268966859046, 
      0.00000000000000000000000000000000, 
     -0.42538501367251799249997046103522, 
     -0.40434782435372039643598964415201, 
      0.00000000000000000000000000000000, 
     -0.16442016253358942481700772352343, 
     -0.27207964379670182797824599426077, 
     -0.31259592531293370325932513075940, 
      0.00000000000000000000000000000000, 
     -0.23061575503872384396183209443032, 
     -0.73174496781264133588838797582992, 
     -0.55942199595539976634905924759128, 
     -0.28187070033874590044167993901641, 
     -0.63181925527793262592959928916708, 
     -0.95845151882990183009015464856753, 
     -0.83806929801903041855622721765482, 
      0.00000000000000000000000000000000, 
     -0.73014202636963239486058389001507, 
     -0.38481648423499752372487549629974, 
     -0.74786883044034466526785200148550, 
      0.00000000000000000000000000000000, 
     -0.82569754738333487928005814615627, 
      0.00000000000000000000000000000000 };
  static double ys[54] = { 
     -0.56977344334973953415009221378062, 
     -0.11615238463710973900986144045138, 
     -0.16841403059063582703799040319519, 
     -0.57273896147926731246264811881990, 
     -0.48177391174488216182747312872034, 
     -0.21805432181592896862251556772349, 
     -0.43015355484213642389086647596671, 
     -0.55672301202474317753018700126839, 
     -0.57361430082863865627675251023354, 
     -0.71062437062144605538166336071450E-01, 
     -0.57425388147699645590520477695232, 
     -0.39157484051927385758800981076725, 
     -0.38260125928736646468021583934304, 
      0.76210815817498206981551181976317, 
     -0.37974002964326184444303529558949, 
     -0.47005865572418340661540287763308, 
     -0.36610435880584296669017948490237, 
     -0.51041551846753395228797796020113, 
     -0.54750015447026086336476242420687, 
      0.11000415018640989305561291676418E+01, 
     -0.55212750394609514626346834644682, 
     -0.31033539412300214897853090795097, 
     -0.40217042416953081560570497781101, 
     -0.45041551128090115382074323796237, 
     -0.55393733579619755005613977870791, 
     -0.50261462137952603540673793837763, 
     -0.57241680160442647300899412288108, 
      0.58655605653182625196365551350355, 
     -0.26028645488038103161700300407634, 
     -0.50737982602401753884163334078742, 
     -0.52365228029452559568466482181594, 
     -0.50930584971577307589941102218584, 
     -0.35752411415512963812671499048932, 
     -0.45339755648003318420528024155758, 
     -0.30993099181146801236771669846609, 
      0.10264156710616244871686286975762E+01, 
     -0.57254647239454006820290004274071, 
     -0.27970449294975765316042730395238, 
     -0.57171692937562625460547538450204, 
      0.42124024273454953036405961333276, 
     -0.54949681464841745643044481184475, 
     -0.50173031881501531356011446242568, 
     -0.45905866008218877559480919998558, 
     -0.44792709733285852771474443606843, 
     -0.51472732382392613411706398642860, 
     -0.57244321215608923333186877139555, 
     -0.54626015892447713295423347103218, 
      0.16065133195264123594739464348022, 
     -0.57182062792860608241992084882251, 
     -0.54913090515973984664975730120048, 
     -0.54631988756965195283384974383402, 
      0.87522865731065508967027909558189, 
     -0.57117161215160068889104816203223, 
      0.11420800543675486817879702813617E+01 };
  static double ws[54] = { 
      0.10297002991753660501136103351341E-02, 
      0.71409200454957424138919555702907E-02, 
      0.84670476850581833366989493797193E-02, 
      0.10061525828127805509276236358357E-02, 
      0.25008518613742033230286581854602E-02, 
      0.82379503874691016163177015741387E-02, 
      0.31645234824978289216154005178556E-02, 
      0.28454670499457543125480683962546E-02, 
      0.90921496761479735112870080703174E-03, 
      0.53162773859582871544697296720480E-02, 
      0.55174868571100333634113037140252E-03, 
      0.57658467066484687646578171852548E-02, 
      0.64457871132908107545457606397354E-02, 
      0.24406168935034987856598180879293E-02, 
      0.58315103312031023476998033502050E-02, 
      0.61853878639532303634906561313174E-02, 
      0.74949719610539269221700456157038E-02, 
      0.49000462324377305053951530096155E-02, 
      0.29889444058930598950088528689227E-02, 
      0.43996809067738625074189199968923E-03, 
      0.12505146482203086295392343808630E-02, 
      0.97825105755546869735655412908223E-02, 
      0.81418004583935088155390428216798E-02, 
      0.50503317417279426125202432815639E-02, 
      0.24951091420867912577241746709087E-02, 
      0.30478657588270378768815691124488E-02, 
      0.64260046142092336339976252566979E-03, 
      0.38189427195873096307259859473584E-02, 
      0.50663628845451396008230346642253E-02, 
      0.51764441731398308561678437452717E-02, 
      0.49181464598226297253857270629541E-02, 
      0.46465662734804209629119572618546E-02, 
      0.44591898920056385785598751636219E-02, 
      0.64136254490796000834798012697594E-02, 
      0.79881644631115754408349534866859E-02, 
      0.10476275740484627938272086913711E-02, 
      0.14500024170436655618156678953819E-02, 
      0.96194255594278304779191550064611E-02, 
      0.15445701615508816200489059359044E-02, 
      0.46888846948789891123342520895612E-02, 
      0.37215235791213064892823990494283E-02, 
      0.37005626494781645355545034117062E-02, 
      0.57612742538281776709493292272587E-02, 
      0.76066250692246841971988592588158E-02, 
      0.43715876705429943949816550672880E-02, 
      0.42564655054254956754287557583451E-03, 
      0.20530208586099517375295895636913E-02, 
      0.57438263789397392552513701846135E-02, 
      0.11429571862723729149581721669044E-02, 
      0.36214122455131093638488339463188E-02, 
      0.26987232674585178812827363148936E-02, 
      0.24397941805574716783727001505373E-02, 
      0.10319996026161764514307891440984E-02, 
      0.11509579298275706789740693353329E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule40 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE40 returns the rule of degree 40.
//
//  Discussion:
//
//    Order 40 (295 pts)
//    1/6 data for 40-th order quadrature with 58 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 58;
  static double xs[58] = { 
     -0.66411637423956578739874309909381, 
     -0.24207197863738398306862447542249, 
     -0.67514797177612497189134399646036, 
     -0.87225015944191233374717005469166, 
     -0.76868567879223402156192530043099, 
     -0.30999131154608929623930725445334, 
     -0.85323240841065682042462033708838, 
     -0.53727030292133569609030595551578, 
     -0.76757602207091515344003545816282, 
     -0.41709706260058898440960904281784, 
     -0.64250324642862691442102353901370, 
     -0.73395730528035636744565078095249, 
     -0.43310798704792516669286002220314, 
     -0.84693242965372475993997325395966, 
     -0.52988575538367881737400733664267, 
     -0.81395203621001148192221956701526, 
     -0.13262238321158207908085298608017, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.12565893500209059730381606002675, 
     -0.25582836243545168466615122195358, 
     -0.41982097163871183619680245728246, 
     -0.39584244571866434364968812226225, 
     -0.13831881077136084719064254716209, 
     -0.55558835329547161475407470013715, 
     -0.26904113625079335044301504651780, 
     -0.38617376855258266995224135587201, 
     -0.90712290003268365546639318057557, 
     -0.13211254083412060410003056730975, 
     -0.13532200006090754449420198166522, 
     -0.51005162143835566342423604409501, 
     -0.51567412876653561433770558683539, 
     -0.72835797984670651201858245642167, 
     -0.62260740206958077067832393367003, 
     -0.63528538940605775324683142950171, 
     -0.13545478402323854515885777980547, 
      0.00000000000000000000000000000000, 
     -0.91834711208585899367680729069800, 
     -0.26841456178680462127694701451016, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.29450956974500233083745593494225, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.15766005001997588398000791219776, 
      0.00000000000000000000000000000000, 
     -0.39500420341927114646355151929679, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.26583541850389618832328064074218, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.96488163632640580772795528150790 };
  static double ys[58] = { 
     -0.57470378380535839956317699309199, 
     -0.55141831198597603831507286345131, 
     -0.56556494724168897227639164041706, 
     -0.51618467078933129623213777000412, 
     -0.57460672418387407371421627702399, 
     -0.54008229992950442561714306153720, 
     -0.57384241303948358177304914166735, 
     -0.55074388383078986548196008116722, 
     -0.55974233776493979903089310270512, 
     -0.52536474950509016897082596544373, 
     -0.54449319151536051271909078020569, 
     -0.53103387228266569983647483538750, 
     -0.57587141183196674121843247917653, 
     -0.55611797625038840229584838026927, 
     -0.51296985906932015988541254967364, 
     -0.52158827793142738546379723447015, 
     -0.56187932687266342953378211324032, 
      0.11003011687008237336733434466617E+01, 
      0.13663584323344421793456675059387, 
     -0.12737340589454779203738364883373, 
     -0.18929257322289204287244570649377, 
     -0.25237212852983271455023246046020, 
     -0.56044801444994193030345087522119, 
     -0.47566799110898290470517971819636, 
     -0.52410356380692977799609815127570, 
     -0.57170940642241625229835788381429, 
     -0.49908301034879271805709825139645, 
     -0.31612494675971333985536967359669, 
     -0.54847061860786860395005750661379, 
     -0.29657575362028796765551235884180, 
     -0.46558820109446506663286885353751, 
     -0.37923882528380845099383801241958, 
     -0.45684803758863608291625036678218, 
     -0.48029591121869609066590186543256, 
     -0.43603282267888446054028374651391, 
     -0.50005009951296356150064101928238, 
     -0.38891721146901775531494918295119, 
      0.27933482596036362512959496310960, 
     -0.57179705694461496519432143619911, 
     -0.43363067211489722444853507200367, 
      0.42523996391594593545924154391851, 
      0.11441090108726184667750585003180E+01, 
     -0.24234470017351667426024672589549, 
     -0.49681561132122740420073482664742, 
     -0.57125217727233967706758863119402, 
      0.93003191476975548336483377359380, 
      0.56962631638217381356419889136285, 
      0.00000000000000000000000000000000, 
     -0.57631399232172626259918827974601, 
      0.82111436271392253740638163140260, 
     -0.40469300239711789646100476636550, 
     -0.54513955979531860983983136957607, 
     -0.57203141564320759086953238463017, 
     -0.35075660264700114997055999049565, 
     -0.34273351742727767767145394888763, 
      0.70353037724633104142734855780744, 
     -0.42794958420782980247198112865939, 
     -0.57216175902539330683834475050253 };
  static double ws[58] = { 
      0.49688685986076335500813888089767E-03, 
      0.20421322651151844017762046611226E-02, 
      0.11376335398163363442669464305672E-02, 
      0.13123691749460771142584854517446E-02, 
      0.53837565018711308131679874030211E-03, 
      0.26896172405115405725156808433702E-02, 
      0.54307048276333835762068445369337E-03, 
      0.25198826411461542126409010399677E-02, 
      0.15319142782715693899177741951315E-02, 
      0.38411906228772124799394024016653E-02, 
      0.25562607929590933659755695536622E-02, 
      0.26896816019905783851160754457047E-02, 
      0.60036514084449199244017297704017E-03, 
      0.14715716460982722722064995019804E-02, 
      0.39962224319702302486375553437309E-02, 
      0.25534279525811038857142356722120E-02, 
      0.25201810346020069847466083850254E-02, 
      0.41167230724256272754400665111648E-03, 
      0.56555181082450911185900550967659E-02, 
      0.56586789117358307508299653903984E-02, 
      0.10834097702542315480967100797076E-01, 
      0.99225669815767172945689578107284E-02, 
      0.24415471132838240081093498781218E-02, 
      0.57970709422402784901464121587120E-02, 
      0.50105635706238936498550178535846E-02, 
      0.12758608174913493052074050132044E-02, 
      0.55592023555177532864401800823748E-02, 
      0.87115223434972838093595378352456E-02, 
      0.13179261092062112504683572450985E-02, 
      0.99446646389628687350930319624097E-02, 
      0.69734880107002948832233705466873E-02, 
      0.73103828906969970566980130701137E-02, 
      0.59302091834307146929707082953472E-02, 
      0.40971162551422433954212787599193E-02, 
      0.57773519974100858251570872273170E-02, 
      0.42621746164839825812632309788527E-02, 
      0.86522628980919465385085382556473E-02, 
      0.52900487244341074555301014519952E-02, 
      0.60794404513723632450482344169488E-03, 
      0.73967700759539224320693784348360E-02, 
      0.47456468096285949486134939736546E-02, 
      0.80914786877822581165281163168211E-04, 
      0.53082561926097671621877127564791E-02, 
      0.30789395107446736330440294150100E-02, 
      0.16224159778341893367110818294581E-02, 
      0.17880830179715579824918122606522E-02, 
      0.41061127743554988075477642960788E-02, 
      0.19270579632685484848497884675987E-02, 
      0.59093070603955750104585123796686E-03, 
      0.25603623606255944146797327158473E-02, 
      0.74860701764472599200155271263176E-02, 
      0.19931737634037231520593140572046E-02, 
      0.79295305880649344360913599297743E-03, 
      0.88807475399574394813610532572560E-02, 
      0.47542421811652855221590656198943E-02, 
      0.33799358355588133678392461897413E-02, 
      0.40037958142737030132660946973599E-02, 
      0.36660639965830984726685636110253E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule41 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE41 returns the rule of degree 41.
//
//  Discussion:
//
//    Order 41 (309 pts)
//    1/6 data for 41-th order quadrature with 58 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 58;
  static double xs[58] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.90707509366054513441839735030000, 
      -0.83678266843431511046582515644399, 
      -0.94247957093734042414919187055867, 
      -0.47818348033753518695149771130636, 
      -0.47786178727052685018858128857025, 
       0.00000000000000000000000000000000, 
      -0.82977481051684242521921168766497, 
      -0.80506207604730995380052139962066, 
      -0.75892338634994084017742272212251, 
      -0.66697937351822451848579519358950E-01, 
      -0.93039119591469919672730326618751, 
      -0.10179246544689757723568411908004, 
       0.00000000000000000000000000000000, 
      -0.11964378950220348999333477645940, 
      -0.89370315180977679131356546281919, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.12537506708828195458087022987466, 
       0.00000000000000000000000000000000, 
      -0.86060405291382569981156454363751, 
      -0.37246528778512515239684348047261, 
      -0.25224760836480194469992557046715, 
      -0.77250795122178423033233897226693, 
      -0.97281903597154664921569844697325, 
      -0.20361834007684187876306722478542, 
      -0.14229399114449109919839299325896, 
       0.00000000000000000000000000000000, 
      -0.71543818078714833155329003431939E-01, 
       0.00000000000000000000000000000000, 
      -0.34072584339074050558489837485303, 
      -0.25420500317361023538457045471745, 
      -0.13311178834618634999839119772768, 
      -0.42835120414075877047497032365690, 
      -0.28308678232916800143508359259067, 
      -0.38982610581001483073643842554906, 
      -0.71358109961060030283458818362239, 
      -0.67177893900431715796469854575323, 
      -0.75059812077502214633644715294448, 
      -0.59464715111368400987629742768832, 
      -0.55914730411574729626305044953771, 
      -0.47513641482741620402395711593316, 
      -0.55098939175943078591819138664469, 
       0.00000000000000000000000000000000, 
      -0.66906775613675318713614198789292, 
      -0.42086181836543388245829813956203, 
      -0.21245173597369898203729302395941, 
      -0.25800420807868548956655134212985, 
      -0.29370336187542435095775000461405, 
      -0.60047596226514389492291181185375, 
      -0.64378659978871008420405373260754, 
      -0.52124100620529880265092320687752, 
      -0.35101185488933557658235806095652, 
      -0.14878346974678461727537723668515, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[58] = { 
       0.11484303236078189653199134730205E+01, 
       0.64867318675503157654602022450247E-01, 
      -0.13336487938297912910317293715113, 
      -0.50914331749408594839577823579785, 
      -0.54016766573740110842549247829795, 
      -0.55659144758280954666377007954550, 
      -0.41036116782557410608240092285442, 
      -0.35568548495723222120621918256555, 
      -0.51447322281830402728376303287297, 
      -0.50391286234357388108253695826115, 
      -0.55710779590286589800031532601371, 
      -0.46781272099145912035205586389519, 
      -0.55164499869306851582964960390351, 
      -0.57416322708744445609317793776209, 
      -0.13038223693758752754332823314062, 
       0.66182123823402367862263609405308, 
      -0.51461673052711459891571855445469, 
      -0.55908758428750220767836487535657, 
       0.54975195106341269371593215963727, 
      -0.31822881181210445991186635948848, 
      -0.31576402673138936793850567085049, 
       0.42731758372108878724118641820065, 
      -0.57305121360034517332191987171046, 
      -0.30895554368016513106704832983371, 
      -0.31414513744457014202355344166677, 
      -0.57323223918730283416506419525826, 
      -0.57206883390058007347834073977517, 
      -0.55123879656195327272974224654403, 
      -0.57235407491813576252917143617978, 
       0.24454257138245059331980778970388, 
      -0.39589172778651022907750111971765, 
      -0.57251313359481391919877047590178, 
      -0.55166089418268212499172347118887, 
      -0.51370572700173861398292804694009, 
      -0.22603324747774107959414273668825, 
      -0.46374623851729501975734765797479, 
      -0.57244266277976037743059136383490, 
      -0.51490548014106772939857400939256, 
      -0.55283101652348038522394903514786, 
      -0.46689252269607195573705009252807, 
      -0.51923400712748281983594166768317, 
      -0.40487430456712229688846342248602, 
      -0.46537465321601328935261456909651, 
      -0.55154491728129470377025145643426, 
      -0.57240622147657917252995100547841, 
      -0.22943850805370712324762812120392, 
      -0.57244722857070476735714053106836, 
      -0.57243826525037158906984485750360, 
      -0.39438908961583798768936375721498, 
      -0.22157162536728802662765970391773, 
      -0.46047385175478172230586878054650, 
      -0.55148536852211283444706530727333, 
      -0.51606114317066720323109583010886, 
      -0.51455997052960895786677902884130, 
      -0.39271417796087884909400497182759, 
      -0.46204815328845467442262509683753, 
      -0.46220719577895767618412315479348, 
       0.79756730028191280539809974543073 };
  static double ws[58] = { 
       0.34355666699665572908084305656243E-04, 
       0.42399098861468427455877995758264E-02, 
       0.37154110284138247834042735653067E-02, 
       0.12458404025090231132925565144087E-02, 
       0.14185150154704695672229773526671E-02, 
       0.72851462067084252604660686856037E-03, 
       0.49283298823596391034884407259565E-02, 
       0.53613423093296115706064365551043E-02, 
       0.19577303889005218175716306314385E-02, 
       0.23541334060017086291762687738959E-02, 
       0.13574312941869499732389430016423E-02, 
       0.33325721812546689467184010795628E-02, 
       0.29644284948173195699449497133337E-02, 
       0.37091902076730237290376100062468E-03, 
       0.80403871323958185654768572258748E-02, 
       0.27933379921942981118737786316620E-02, 
       0.43548510260672555779566602140770E-02, 
       0.98750346722284897815375163954063E-03, 
       0.33140399984357008523996740964860E-02, 
       0.39548685266133002402731823562045E-02, 
       0.81907961822989578237987369917385E-02, 
       0.39990415768416250914981967446877E-02, 
       0.64543016375148533251238361162029E-03, 
       0.77525238376409612296793051209554E-02, 
       0.82098277784420498504581064200822E-02, 
       0.77505992315050355696141616800802E-03, 
       0.31460954824696309149912473051966E-03, 
       0.30771570801706524162859379333622E-02, 
       0.13758671651214991382998623408107E-02, 
       0.44530286851553646989910342775063E-02, 
       0.78835322339101614146632253638963E-02, 
       0.67359298176593092050484809442566E-03, 
       0.30115058169438009830398484639718E-02, 
       0.47516681892344821461215335464754E-02, 
       0.92897161081570613163976545109188E-02, 
       0.57624992565964513866295347402457E-02, 
       0.13386293617916449653435311528237E-02, 
       0.45476733851440924400633719363253E-02, 
       0.21881299772745267054828351965957E-02, 
       0.44103886097359621239681297088827E-02, 
       0.31408957167973836083480622860722E-02, 
       0.54951077871372566664971370622657E-02, 
       0.50946823920529338224773265834011E-02, 
       0.28852806434905299950945859499645E-02, 
       0.12012233088773826615649104421917E-02, 
       0.47576797651663080029937010231751E-02, 
       0.10565898230048048225468444873170E-02, 
       0.12889097209731337671344489509095E-02, 
       0.78373349308822302058968405245672E-02, 
       0.88868099593492497317748964708364E-02, 
       0.64917648662630346100035827874498E-02, 
       0.26315125248814009986883523206499E-02, 
       0.37615631205534999461032753334140E-02, 
       0.42869570500049250611205071236900E-02, 
       0.77283572873216484561057801262313E-02, 
       0.67055003387097703603196736320573E-02, 
       0.33782766604571404544711257767029E-02, 
       0.26121233276609877149517818702563E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule42 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE42 returns the rule of degree 42.
//
//  Discussion:
//
//    Order 42 (324 pts)
//    1/6 data for 42-th order quadrature with 62 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 62;
  static double xs[62] = { 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.33732186702983072527377866348931, 
       0.00000000000000000000000000000000, 
      -0.96033929602597311557761977477877, 
      -0.86369982899641926289375115271667, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.91964505433817813657264837251525E-01, 
      -0.66309355251219477483888932315586, 
      -0.51372916997472296488281957020097, 
       0.00000000000000000000000000000000, 
      -0.93438263780905817601029667101954, 
       0.00000000000000000000000000000000, 
      -0.90377238041680626119957295786590, 
      -0.97636560454756210194148963029962, 
      -0.80976763015790175270655294347859, 
      -0.11503943729553541704016281376984, 
      -0.45642869772220040454862543963416, 
      -0.26081632436550845472752474357332, 
      -0.23737442933070842287301049364327, 
      -0.48002697917900107156245025485608, 
      -0.55650162501787387537981069378810, 
      -0.12478176574417714212774268424213, 
       0.00000000000000000000000000000000, 
      -0.57250293774060918335092055220851, 
      -0.35082904042609762405060413462916, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
      -0.48329385573259627109084340701262, 
      -0.60723425879311620189016063272093, 
      -0.74479369220077260375364623234338, 
      -0.88899121520171718881192636679144E-01, 
      -0.64093927012580804843536900337703, 
      -0.72840614316212545054184733409861, 
      -0.34632334976132410062989430769368, 
      -0.29328359092418193020774881585849, 
      -0.52075640841362497935400311743224, 
      -0.60555027935120159340828726711903, 
      -0.40905454399991287433613498814895, 
       0.00000000000000000000000000000000, 
      -0.68384877689179067855591079553230, 
      -0.71845211882128619528987877814782, 
      -0.21667354415688649871182254929595, 
      -0.78401417474182939603497549800273, 
      -0.73829866676407222514323372659105E-01, 
      -0.13454529284893910979363670735151, 
      -0.81224627731204864480024020247250, 
      -0.21779874120415795023777409358666, 
       0.00000000000000000000000000000000, 
      -0.42173742894104082835679644928203, 
      -0.86724517518683599270071051238882, 
      -0.35974948772649811066870719069685, 
      -0.38462016485447150927117414922040, 
      -0.14709579616219243877526898806227, 
      -0.14236698354849056304093448822305, 
       0.00000000000000000000000000000000, 
      -0.21806859696894763488778637237838, 
      -0.28058519641950796678974392628735, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000, 
       0.00000000000000000000000000000000 };
  static double ys[62] = { 
      -0.23915669727967103118380139768186, 
       0.11499706797515916998220806401861E+01, 
      -0.57463021225838519900335995927075, 
      -0.56879992979770302521884041359830, 
      -0.54196210567272946906489712244628, 
      -0.53941219280079277870643403803930, 
      -0.50772285124871792168226530442572, 
       0.92445442054483368113851672487764, 
      -0.57476336247113177371961348404936, 
      -0.48132277490581309659947493263723, 
      -0.37176697642006607510778844628504, 
       0.10049252920684353614918795242733E+01, 
      -0.57493742331160895521040050048537, 
      -0.54376585322681251262386416118604, 
      -0.55951105488338182442194601355457, 
      -0.57287342016015178837509330570583, 
      -0.50897875397137854857739836711645, 
      -0.55605671624505328368102486115945, 
      -0.57286458799044539871933065364668, 
      -0.51230069551399328623142182353309, 
      -0.54774164573700082169258700885382, 
      -0.42881250819630884111523888550059, 
      -0.47691618600212933305894604508553, 
      -0.51872351401291148614495359408056, 
       0.10611738733054864308352869822745E+01, 
      -0.57297391133539792503817569661015, 
      -0.55889259681253444754878375605606, 
       0.67490581305685926714955055515407, 
      -0.46308066751251361359457711328632, 
      -0.55408939280498580695400859722506, 
      -0.42433782075470535918840983649065, 
      -0.52032748775510677492310958477448, 
      -0.23994545417022376699165559718980, 
      -0.52285819108589385499701905311092, 
      -0.46066493316460819615558728243469, 
      -0.25420978730877827119837106988602, 
      -0.46670083015487351222116455572352, 
      -0.52093525817617790107842375316544, 
      -0.55435992266257116866227670316197, 
      -0.34138644830088341585231904393856, 
       0.78549511301950868330086537396862, 
      -0.57282030241995550160692247503933, 
      -0.55400012643472741059711642382680, 
      -0.24499725128675631298290375529787, 
      -0.57293824798152066163974973440469, 
      -0.40234854190615793941968454639826, 
      -0.14476929626336181753272157091150, 
      -0.55359930548468954074315532763614, 
      -0.57155725996820959476202397056243, 
       0.74817466669102061182393603775266E-01, 
      -0.48091898059119731232834971500976, 
      -0.57249401685345627316687746770457, 
      -0.41046630772117739922581443734300, 
      -0.52677593302205091340352395053193, 
      -0.46723591737438723999494490231575, 
      -0.32805865585485880306160348707314, 
       0.30647862202829743938620359376295, 
      -0.40551576066406470095805909037152, 
      -0.33418952446785743728949888770460, 
       0.54905382739186092457188995702651, 
      -0.32644308263379735102228951317930, 
      -0.14180931586841959084250934615779 };
  static double ws[62] = { 
       0.23336887717611856261356883070280E-02, 
       0.22846739164921267903695204331260E-04, 
       0.70111433879772498743791597802441E-03, 
       0.66736115582873727562875088212462E-03, 
       0.57785786514264065438532593664687E-03, 
       0.13648017145261238754376858784678E-02, 
       0.18715796130141111886024261596391E-02, 
       0.12311750904528985305433912829374E-02, 
       0.73641178059201489942713406653212E-03, 
       0.33400912254975143418927736967276E-02, 
       0.53408360185582866807409529152528E-02, 
       0.87691656602387380470263190564516E-03, 
       0.29827840345921013265121144680296E-03, 
       0.14428656899060768897862657143696E-02, 
       0.89736946555330205057357879483966E-03, 
       0.24622832560853176787327086769666E-03, 
       0.22197294779296558696560399072253E-02, 
       0.25038963454662955321871600907872E-02, 
       0.10178082444714138863786539439139E-02, 
       0.41750972622189718039985200715331E-02, 
       0.28539390170943817441698264694008E-02, 
       0.53618401327742042860547051137744E-02, 
       0.44379187212093563379576658435829E-02, 
       0.43923515022689111826190793996648E-02, 
       0.62671836400915488143041776852213E-03, 
       0.97955449417192964190566775890366E-03, 
       0.23018208943548581359392383454856E-02, 
       0.28689426702522553567921914933795E-02, 
       0.29825739568587811806411249714095E-02, 
       0.25047104675084795836734599265446E-02, 
       0.52092111037958286914269299684586E-02, 
       0.28471024101248386546695067189237E-02, 
       0.82474954978823005019031886199947E-02, 
       0.31999957236382169770094914684218E-02, 
       0.37199859798352839350442391912941E-02, 
       0.85220480529686853126175155374059E-02, 
       0.55913173485803385482484428528564E-02, 
       0.38298884403760618569504568998379E-02, 
       0.23089127874100000149196844935704E-02, 
       0.70117980296802485289704436077025E-02, 
       0.24321512738535918412817989755737E-02, 
       0.93669661545711090880311728905074E-03, 
       0.20642879082447440946641719056348E-02, 
       0.90936551597467743554493293539138E-02, 
       0.79986620899126243521988101233587E-03, 
       0.77120700714742987371071118791067E-02, 
       0.10219559553814186449346549825429E-01, 
       0.17241106961996584893369591793848E-02, 
       0.13912182187991904328859295409100E-02, 
       0.53283478006725665670848896066726E-02, 
       0.53185152093575276335966863060532E-02, 
       0.68994902614553096607967297139181E-03, 
       0.68499915239657437963246941954227E-02, 
       0.39920622096982370641881855622790E-02, 
       0.64073752422181504921770295563876E-02, 
       0.86898690635024727222144915241917E-02, 
       0.47974852921883757762301535128063E-02, 
       0.74399462871562674711612756662098E-02, 
       0.81490368631924462501235690710386E-02, 
       0.39847607162807112945480365181103E-02, 
       0.44476356973374278433219749730315E-02, 
       0.52129984983515280965046257082396E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule43 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE43 returns the rule of degree 43.
//
//  Discussion:
//
//    Order 43 (339 pts)
//    1/6 data for 43-th order quadrature with 65 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 65;
  static double xs[65] = { 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.39506733155766901012421952821161, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.24646837298065380013273617416435, 
     -0.79575263644167424708927638836510, 
     -0.48253886188633001772041546826783, 
     -0.62546966261193450084215465475134E-01, 
     -0.30130814355531906520968875282228, 
      0.00000000000000000000000000000000, 
     -0.35231165502281585353535246732493, 
      0.00000000000000000000000000000000, 
     -0.12800312652182834294317879260650, 
      0.00000000000000000000000000000000, 
     -0.75118628320269714215689196084969, 
     -0.67034428298913296399375404066474, 
     -0.77580876066209624496610974380034, 
      0.00000000000000000000000000000000, 
     -0.12785487625900558619560024121639, 
      0.00000000000000000000000000000000, 
     -0.71268785841843536059482927597042, 
     -0.27040014117930207740167640614761, 
      0.00000000000000000000000000000000, 
     -0.18478688196683115844840486556670, 
      0.00000000000000000000000000000000, 
     -0.80938151535440202130076995993744, 
     -0.31942848667068561771894552587106, 
     -0.85468012516227660038415955640997E-01, 
     -0.25825613481914215491217537746903, 
     -0.58279946491371434111606534307889, 
     -0.87969623175882104649066373698787, 
     -0.46561039262024648141549373583729, 
     -0.93340399559638336764338113699127, 
     -0.51396050819491791778081081792183, 
     -0.61691723863061563788005694821551, 
     -0.86363081428402156734207066959522, 
     -0.97099189034945842829889245820946, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.62361278877871758607587901074940, 
     -0.91989119779265366986023199875079, 
      0.00000000000000000000000000000000, 
     -0.57692441054824302159307937523806, 
     -0.50723010453494962638105205800859, 
     -0.84394330147763359879321357561682E-01, 
     -0.72274543590291663139776270288902, 
     -0.47250099165753816338998461789401, 
     -0.35080029727614026576846631364065, 
     -0.12275092202821139885819853243053, 
     -0.38479753550946700083527634968480, 
     -0.39605917796945196543812828983515, 
     -0.84481754164315782828022282277629, 
     -0.20507268425580293655253183641047, 
     -0.52353562229447874684448906255323, 
     -0.13693189087059081697918132267328, 
     -0.62320689311818688506948006480502, 
     -0.14379384782311729469209427032452, 
     -0.27449236039881046467163004113221, 
     -0.68227667429915713110204049156621, 
      0.00000000000000000000000000000000, 
     -0.40313963669589460865080761368403, 
     -0.22073052548745090606319136542233, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000 };
  static double ys[65] = { 
     -0.27174781657834587367794828403112, 
     -0.49274259095938795471287794986388, 
     -0.32541690542417216472121438315947, 
     -0.57450914029685444604595657642940, 
      0.12846740408215895203363603738439, 
     -0.52120501981219468628421175620930, 
     -0.55585048147520898296581708536554, 
     -0.33857300476573982851491415786097, 
     -0.34294121393179207633049958993479, 
     -0.32875457338449147616580244617275, 
     -0.68781371078380984368245670302187E-01, 
     -0.51594631133528140591593205811015, 
     -0.56039688938398587880112568739761, 
     -0.51600357027661430283719712442229, 
      0.10406975917564771325612019382002E+01, 
     -0.48071245484494572753863498860590, 
     -0.46654018005672172190079466787693, 
     -0.52493394350877366855931608373780, 
     -0.53080573778005055886715862589101, 
     -0.55199717395419102378360697714430, 
      0.93744923795883541057323646131539, 
     -0.55378528333018398654846926976093, 
     -0.57306164103038438002759473015017, 
      0.65329413919843916560086065083107, 
     -0.32979522157781538320615791050048, 
      0.11084041345444688623805682120094E+01, 
     -0.57315735914007292751468298272925, 
     -0.23933996521344827119141006030721, 
     -0.46640357835051818210779901438991, 
     -0.55450865650449884106657600267936, 
     -0.47018294735939914621510849234674, 
     -0.57282533435858176236473269144498, 
     -0.51813075309225748766670475143277, 
     -0.57275889464282222126170595591666, 
     -0.57286153985489888900353792433004, 
     -0.55262810787321224192869259072191, 
     -0.55335026706806535430401211392438, 
     -0.57303083539166774608654680411579, 
      0.11458153208366305893213396871360E+01, 
      0.50423840889193913255432466884584, 
     -0.40338147343700279281658672539043, 
     -0.55362560310096266386701313905319, 
     -0.17054110158363053238654111959974, 
     -0.51821251917615874954294201331871, 
     -0.55343995679985670491884892559641, 
     -0.25432761855273552279988165105579, 
     -0.57294781282937958180985481791538, 
     -0.46817056136138252363928334992228, 
     -0.46671221117152839234190202161524, 
     -0.15827564830985659745747052529915, 
     -0.55206358961559880212009617361087, 
     -0.57246915152619833238307709646113, 
     -0.51822758166184266566614267250055, 
     -0.24244875862513978976057242564588, 
     -0.40871265880657287984484480492628, 
     -0.57252078775710633720596892661685, 
     -0.57257680220856426874713292476538, 
     -0.40790874923298185895067251523869, 
     -0.40488358656967401446800162607396, 
     -0.51826247883638312778785904937525, 
     -0.41416427390059670796562799333653, 
     -0.40328359141014275624720285726773, 
     -0.47052810769116907040147328262294, 
      0.83299965359609554712354240867349, 
      0.28006707539158895931997168962991 };
  static double ws[65] = { 
      0.20267313231941893676520332064545E-02, 
      0.16367671769150146246172346306318E-02, 
      0.51493216424315589396152049973413E-02, 
      0.40049166791821914724573044327252E-03, 
      0.38679957027277816063713365008715E-02, 
      0.34975453538986470037287749440995E-02, 
      0.13692247906521867994012332715043E-02, 
      0.50916564180668567506293917701109E-02, 
      0.71541887598019862130827731596776E-02, 
      0.67554356473948161976278050806714E-02, 
      0.45197199185406262459677900374973E-02, 
      0.35861506794982540830122357172132E-02, 
      0.10622739063192419924108232636689E-02, 
      0.40560451238345220270219217309857E-02, 
      0.71306365950803849055449621912648E-03, 
      0.28901178272685770916650253073937E-02, 
      0.35628382297785204095556205282981E-02, 
      0.23337289887646671489220572776197E-02, 
      0.17876266669191564461852559716746E-02, 
      0.27892998079816587686734710640336E-02, 
      0.13320835956879952672184157241896E-02, 
      0.18191641344329572250280779466083E-02, 
      0.10907352970205207235715895372753E-02, 
      0.27936513125358662464080516863599E-02, 
      0.75818102244409073466682226454199E-02, 
      0.31066065630647406501495742839881E-03, 
      0.64045511782145175218433235311442E-03, 
      0.79769491206495963710875582745134E-02, 
      0.52105480536968519556750471994035E-02, 
      0.25683586574526020251013582443614E-02, 
      0.41796721048381089538295214085205E-02, 
      0.54981489589321546105773135980457E-03, 
      0.36900063605190763288766065453592E-02, 
      0.40460366093556559458970343028029E-03, 
      0.99764710822395591487080703901083E-03, 
      0.21433057408916717136795809729450E-02, 
      0.13181554039184077697417192202970E-02, 
      0.24966385865177930712178241334390E-03, 
      0.56415208115924918653995444365333E-04, 
      0.35108086150369173319630529016857E-02, 
      0.51624236567181411900078094419075E-02, 
      0.96766243407032355999668290740009E-03, 
      0.45731471950392280996000615522500E-02, 
      0.34583136168620142190249880837043E-02, 
      0.24237265422044902192234994737096E-02, 
      0.77977816580077882796552949686740E-02, 
      0.80588323551714052775101742970255E-03, 
      0.50339349950361747376652268826760E-02, 
      0.54331525883336925722484960144977E-02, 
      0.88072094639810414255439357569340E-02, 
      0.26943758842246860082040659620724E-02, 
      0.11555336183140089850951172772477E-02, 
      0.20764819925413973070582238653112E-02, 
      0.83010524908151710705626972584161E-02, 
      0.56231102351940614439232706639500E-02, 
      0.12763777281762985620080832782995E-02, 
      0.97128666683373741323713272921117E-03, 
      0.69250253173029413525682166964403E-02, 
      0.69755502078755805652513910381880E-02, 
      0.32320158906108914651035701130350E-02, 
      0.35091833149609820457735281284170E-02, 
      0.66535574606795091616463991758331E-02, 
      0.58357235029408977357372091365790E-02, 
      0.21654064214107489596481597275945E-02, 
      0.48130202892800976037200791733793E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule44 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE44 returns the rule of degree 44.
//
//  Discussion:
//
//    Order 44 (354 pts)
//    1/6 data for 44-th order quadrature with 67 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 67;
  static double xs[67] = { 
      0.0000000000000000000000000000000000, 
     -0.9497329452118762389747014470132091, 
     -0.8001862513364488071091260917428001, 
     -0.6910725775767816500425682030963931, 
     -0.7699895351147678691136080969564274E-01, 
     -0.9389616539760584684425697161489036, 
     -0.3318014560907163140854925533609147, 
      0.0000000000000000000000000000000000, 
      0.0000000000000000000000000000000000, 
     -0.7536007546194699584282605909703006, 
     -0.6418838401483368830531210405780724, 
     -0.8569413055073822333863570460181954, 
     -0.4713099774418540641295127244304687, 
     -0.4587710901095560036744622028451049, 
     -0.3454684993219427932117970013334666, 
     -0.5784975768218464125847982734150031, 
     -0.9757061286693509772479542529554229, 
     -0.9062217459170026624185611454216153, 
     -0.6998855908781156060583940844630905, 
      0.0000000000000000000000000000000000, 
     -0.7828939494539899851949840039195377, 
     -0.5976618502172688281903349218410519, 
      0.0000000000000000000000000000000000, 
     -0.2089756393930555645731581562882102, 
     -0.7011575058845152681611414003803953, 
     -0.5064990519200558712677355544916721, 
     -0.6042610388710118243964657680731453, 
     -0.2426529067325121256109910016105379, 
     -0.4856256157498482835842484912328102, 
      0.0000000000000000000000000000000000, 
     -0.7004374206319495265033410739625405E-01, 
     -0.3757615711643657550148681484974144, 
     -0.2919871922611893373479530841081256, 
     -0.9994102280758133535657317453111649E-01, 
      0.0000000000000000000000000000000000, 
     -0.7857254577834318544379545866687309, 
     -0.3721620253239425678261695592905290, 
      0.0000000000000000000000000000000000, 
     -0.2824529226357365089257330994701188, 
      0.0000000000000000000000000000000000, 
      0.0000000000000000000000000000000000, 
     -0.3953463775800016380774481750278100, 
     -0.8361108176791535739269469385088765, 
     -0.8814365377132475290846046360400148, 
     -0.1406663437184717854538742687564680, 
     -0.4650987866606519008128443728604466, 
     -0.6964523670523666770018160868777266, 
     -0.4113816477062986076261766146150044, 
     -0.1610536446585760688667173579011407, 
     -0.5028509942935229229810541094112460, 
     -0.5448555197999273032450518681162338, 
     -0.6058936424345374593927006219197785, 
     -0.1081927332768446316060675202726035, 
     -0.2295037521073542228709062704046674, 
     -0.5898196614388808672207378844102958E-01, 
      0.0000000000000000000000000000000000, 
      0.0000000000000000000000000000000000, 
     -0.2324653629444445842929660552851472, 
      0.0000000000000000000000000000000000, 
     -0.1819171108132164261481020014379600, 
      0.0000000000000000000000000000000000, 
     -0.2500032323171110666706981553600442, 
     -0.3375819432490318985102882332004112, 
     -0.1161841402458460505963725617476231, 
      0.0000000000000000000000000000000000, 
      0.0000000000000000000000000000000000, 
     -0.1381839808516535777015623012008872 };
  static double ys[67] = { 
      0.1148768930925886262151310465623686E+01, 
     -0.5587190843946788319800702097639321, 
     -0.5759982547513931308236141599939908, 
     -0.5739929781026895530093619240854897, 
     -0.5742565796637987766332210625225684, 
     -0.5743118620654462339294955029963182, 
     -0.5764542892494308227983752102621023, 
      0.1053661869059013126390801025098173E+01, 
      0.9830543452221882010316337450645482, 
     -0.5642804373589460086506076815201640, 
     -0.5574445948250216759499998540402298, 
     -0.5353255763542691825073068820035573, 
     -0.5603801665927884443776809238750006, 
     -0.5740608434639325204622941212626694, 
     -0.5640573337469584267853102545754350, 
     -0.5726612537821012774116991313991636, 
     -0.5726626528471781321664603205600438, 
     -0.5592282822642103798482359755790968, 
     -0.5400889021600139741287112668257262, 
     -0.5348248905002649083989599110468093, 
     -0.5378170138107226759506693144865885, 
     -0.4039918494468322341133479531414074, 
     -0.5637079596777493822434797158628362, 
     -0.5718777342221263057688087793323300, 
     -0.5018890152406374639656521891734000, 
     -0.4685508460352765972480529422854412, 
     -0.5172633713942475966596304508297206, 
     -0.5497203595501590605031193106499199, 
     -0.5140166345967331496952227323818916, 
      0.8893638105980905452621258664733128, 
     -0.3645379523625655382701633893722455, 
     -0.3659340829551636496622119107012566, 
     -0.4252279065273847102094159009996268, 
     -0.4635872615991260396686282785197555, 
     -0.4940450202599493985592231391958893, 
     -0.4981606513359721550788104276334523, 
     -0.4858920787064337532175409440511338, 
      0.7767824952808555984008258545514951, 
     -0.5197375697152384646052054170576424, 
      0.3536374193876506397331972717972702, 
      0.4760546324939761283315701222537493, 
     -0.5373555685767024547552789232104194, 
     -0.5627800449469781889100343115674311, 
     -0.5741358518948284514699798237494199, 
     -0.3143189178536211088756277347160222, 
     -0.3232690739017036811465453364851697, 
     -0.4508096294389419700119239124729038, 
     -0.4330349046242016037349211869278261, 
     -0.4081806944417965784076695553371667, 
     -0.4004602966965545264378386469600689, 
     -0.5480091160789313012122431861765678, 
     -0.4654469444625449906357368312464820, 
     -0.2408240697048009236946195102972083, 
     -0.4755237162067393099887086042912068, 
     -0.1554487149960662937738036182401154, 
     -0.3082651195390385495689250052901405, 
     -0.4210560663521613009285800897713914, 
     -0.2549706395457159319661921698363317, 
      0.6562159018607634589554760136404371, 
     -0.1645038009839102097914264096582779, 
     -0.2390476568695849024278482393008327, 
     -0.3500539054481347352752143235181112, 
     -0.2855707593455760613475971357331510, 
     -0.5546162891229904937801707067639082, 
      0.1352301011214314104583760897705539, 
     -0.6333228349610746330155009819900319E-01, 
     -0.5177444417481954592613791928463002 };
  static double ws[67] = { 
      0.2915373135979705462354669962222584E-04, 
      0.5421904422449839736435313083383586E-03, 
      0.3596051767025996319384007236006714E-03, 
      0.7097961495022424493012412105851562E-03, 
      0.8830579864206354364504096263197627E-03, 
      0.2973426233859937360534817748447138E-03, 
      0.4238631002729904678977975696568931E-03, 
      0.7389606753338469971149288564795965E-03, 
      0.1148374312030633712597584302352535E-02, 
      0.1273653397793799775219820208340736E-02, 
      0.1772310668755885176871822523120592E-02, 
      0.1760259894656385175896251017077239E-02, 
      0.1591696809621534312060249900305618E-02, 
      0.7861912161306376969359324550937249E-03, 
      0.1932261265768957224441734639436218E-02, 
      0.1054246990932389772147217163517965E-02, 
      0.2442199620811447342051688108995863E-03, 
      0.9494447813509459504714159753943138E-03, 
      0.2016444809029372508564748979562772E-02, 
      0.1796212502278787821852696540715106E-02, 
      0.1985990719999202087365783713106505E-02, 
      0.4476918545924685844030658743575517E-02, 
      0.7926095313482192607115919198160405E-03, 
      0.1363437100111953041990336948238837E-02, 
      0.3128949768581953030780660470515785E-02, 
      0.4468586167117141199592656650169745E-02, 
      0.3558780448911060819562296812185439E-02, 
      0.2654525745999437126651160840433368E-02, 
      0.3674322205637929576474475378034387E-02, 
      0.1556779548166703629723862770843712E-02, 
      0.5659722156959053747811933867465999E-02, 
      0.6559124062355508037637436719411665E-02, 
      0.6124082302973937834781741237000250E-02, 
      0.5432590387902965678105625419649959E-02, 
      0.2310495765495803855211312734514998E-02, 
      0.2708245744720212530179525275012375E-02, 
      0.5124555916440916550164562572827943E-02, 
      0.1936910902435187766318955933882164E-02, 
      0.3646527070098702883931381699630589E-02, 
      0.4424142630723912210455942532277462E-02, 
      0.3276931385664301994721671259073124E-02, 
      0.3264053064814834336425695074922234E-02, 
      0.1135643135915063583134474479188575E-02, 
      0.4223718148997970920874463333812100E-03, 
      0.6555375990726206244644356106443562E-02, 
      0.6606136053977267350584211383354891E-02, 
      0.3644515347971302216660166325983147E-02, 
      0.4871554113050847936724114713691689E-02, 
      0.6679124914763547321632436773214387E-02, 
      0.5616578389429630050906029716034673E-02, 
      0.2275791117266683412259386884315207E-02, 
      0.4243996397034961679236635195585269E-02, 
      0.6959330566880718539247596084458503E-02, 
      0.5290070422157258189378060153617274E-02, 
      0.7976505080967100010251007180513698E-02, 
      0.3238093010218998159411476040299402E-02, 
      0.3304316395913585112942109463963793E-02, 
      0.7632093580324891025483614928517932E-02, 
      0.2926022905168081777270491914230688E-02, 
      0.8413855313504570550721436602221883E-02, 
      0.3057643153861432817013016432152174E-02, 
      0.7441970831021575269373809133451851E-02, 
      0.6989734119088538264792355031494924E-02, 
      0.2596917279271723848812657753581774E-02, 
      0.4228198393093499219382464343859585E-02, 
      0.4122988694755634966269440274723658E-02, 
      0.4679274136115308847400441526701599E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule45 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE45 returns the rule of degree 45.
//
//  Discussion:
//
//    Order 45 (370 pts)
//    1/6 data for 45-th order quadrature with 70 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 70;
  static double xs[70] = { 
     -0.96749796278768431752695073746986, 
     -0.98943904955500689467050133968185, 
     -0.85724783085411714050106137680522, 
     -0.86009422327718343131869480185428, 
     -0.54678180057409276098388610849729E-01, 
     -0.49732663296837674567709737943128, 
     -0.96239582473679000067551742403253, 
     -0.27961070235673144034368245271631, 
     -0.40402890296751249677645751975093, 
     -0.29442839410535855309352484266892, 
      0.00000000000000000000000000000000, 
     -0.43275119028613663914670241562698, 
     -0.44044389245586822931301266425545, 
     -0.61279883033232538372205093469042, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.40788451897774191002745407172850, 
     -0.69268856593560823389350549597623, 
      0.00000000000000000000000000000000, 
     -0.76676405679321178282204182266112, 
     -0.88142494556619689489685668224494, 
     -0.30788563677170870686746228513001, 
     -0.51690806387591814992569974157337, 
     -0.36522335396393722816446107730499, 
      0.00000000000000000000000000000000, 
     -0.54085071215301248177693480228527, 
     -0.83755436630162203436831978893993, 
     -0.39400450916939910389752569327788, 
     -0.33507526831923579205690213103223, 
     -0.93207383259146302186264561603832, 
      0.00000000000000000000000000000000, 
     -0.24725397282679600089234012793255, 
     -0.59548295462739307737882050980714, 
     -0.21190291318947065883045380978784, 
     -0.48065853549210622637194679384859, 
     -0.58563257009522274955669451569327, 
      0.00000000000000000000000000000000, 
     -0.36489263643259686426301607599478, 
     -0.10662122389563936276094546991520, 
     -0.63098646216816865329341790460037, 
     -0.10951950034472544896736583276610, 
     -0.14147353942031966219611357193348, 
     -0.32195104149719150226124582347371, 
     -0.79909003824949698793707745791143, 
     -0.10582426030688274991436338221365, 
     -0.11204595369159801121778549049946, 
     -0.61200208537203979332083793516425, 
      0.00000000000000000000000000000000, 
     -0.69310856751575954320818550365882, 
     -0.52182537044957879662933711116161, 
      0.00000000000000000000000000000000, 
     -0.21330305926164315152747136574011, 
     -0.76353151836479030982601065330425, 
     -0.67884482751176626602735151717562, 
      0.00000000000000000000000000000000, 
     -0.71309872195622120220433316560504, 
     -0.23702202213867019459513762948482, 
     -0.49629255197942964513516871846921, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.21703474679150985775348937929931, 
     -0.12791164301131482078856717349660, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.91918133084949909274827769251893, 
     -0.78297196007811101735741996204489, 
     -0.14461283051109006329463020198732, 
      0.00000000000000000000000000000000, 
     -0.26734836171405202882645924901152, 
      0.00000000000000000000000000000000 };
  static double ys[70] = { 
     -0.56471941069694844856954469658100, 
     -0.57501200548013563442214533934183, 
     -0.56546982796685102989871197087213, 
     -0.57522024853914446403759310165459, 
     -0.56214268954647688248601400572909, 
     -0.57443518202158385272046704794250, 
     -0.57430537135662979697492111738520, 
     -0.57479443137000535377443032432917, 
     -0.31119792517448315773562579509789, 
     -0.18309182975606309835985477370170, 
      0.70503710102613623076996776216914, 
     -0.37954480542530077301394984603928, 
     -0.43813658252872489137658216590349, 
     -0.52431770764693138795353125680755, 
      0.43988432838417361094660056896797, 
      0.54818973150740130019484360782211, 
     -0.53068517830945825690351866312342, 
     -0.52054382897965220831059074194035, 
     -0.57563407311774753488214306675980, 
     -0.52066524714991149578624581683571, 
     -0.55109220346610688251784303249678, 
     -0.28136936804129377678526074056198, 
     -0.52770400117212017756160700268800, 
     -0.55536072100117397125270524523176, 
      0.96341139552317968299064589061772, 
     -0.42510034407525727010496874963780, 
     -0.52414106878540363881135726636148, 
     -0.57247532448088443486712017171738, 
     -0.43906434358617035923094688491037, 
     -0.55427727760435318351185124343632, 
      0.11840495907177995488840441336373, 
     -0.56071879514252223550399956253037, 
     -0.57301592049564293598812570618238, 
     -0.22794222552637819257304990263196, 
     -0.48610984701993063250622405741056, 
     -0.48075527463104237978284559370263, 
     -0.11656336099177678980294986355914, 
     -0.49380185173253747396240691815876, 
     -0.17257418086954546445487547770085, 
     -0.41685471639662622825619311565566, 
     -0.27691779802002702321949855249422, 
     -0.54643161443323592717327912197685, 
     -0.36959650358697796007635862331191, 
     -0.55312820158517175524163271806254, 
     -0.37061655711808016877061539071484, 
     -0.44989790809755280647980798667131, 
     -0.55532421580891103975772368666732, 
      0.23698504908978765362861518397758, 
     -0.57312340633371483370089991608323, 
     -0.35237310557002597123601067611590, 
     -0.22563680509593199790438697496259, 
     -0.32814170713377593953818022888866, 
     -0.47540208897837744967251438098628, 
     -0.47291221034981050011352430482453, 
     -0.41322179544695816391477350068005, 
     -0.55444311333444255226437372754811, 
     -0.47815353791871963047978697256372, 
     -0.55841263605812562308032038707870, 
     -0.32398951274847599369007130866498, 
      0.00000000000000000000000000000000, 
     -0.41120300953009715133424788540273, 
     -0.50863906936384494919608030957095, 
      0.82709200615759949176722321241391, 
     -0.48378305160084771830753762022321, 
     -0.57250910447612774904686920854490, 
     -0.57273156232313427034131360506698, 
     -0.57295167407019810131965208299179, 
      0.10366261458457849101902218371190E+01, 
     -0.52758606695132708413591860408544, 
     -0.53231347309540434406510538720317 };
  static double ws[70] = { 
      0.25572163619669275347030002593120E-03, 
      0.78130681960434454688570462017890E-04, 
      0.63262191626397677725423132811818E-03, 
      0.30286144471206534618188349487863E-03, 
      0.17036087426692587176471947445535E-02, 
      0.60607958666190953424951755070468E-03, 
      0.20408229164689422092414467361640E-03, 
      0.69282349702782943622699287812265E-03, 
      0.53664270813767778699249041260574E-02, 
      0.35893775549966510288707779886923E-02, 
      0.21073089951061989088486797402699E-02, 
      0.50245578513041966694332486164703E-02, 
      0.42688044135583548301459101168732E-02, 
      0.24807543972727024276437592523865E-02, 
      0.30515471493979938480916234139891E-02, 
      0.26894370339966747219062415099150E-02, 
      0.28699920637932601309100432329669E-02, 
      0.24144446935100309739631679872824E-02, 
      0.31790657648418738011813407939271E-03, 
      0.21416996150124135443929791578456E-02, 
      0.10904105148236371090982266813776E-02, 
      0.67218349746074238515473176347085E-02, 
      0.29468389036919417335879077316037E-02, 
      0.21585882255773739347679137982008E-02, 
      0.11138001100223078504096300255199E-02, 
      0.48321913046622933674773381791937E-02, 
      0.18596190838124657102588736539343E-02, 
      0.10060622517439881677085687475045E-02, 
      0.52037982694409060390307613806173E-02, 
      0.84864867119809259804133052600579E-03, 
      0.44968639304153501776835069862992E-02, 
      0.20461751558888743165348782978731E-02, 
      0.82449377452441123159564369966813E-03, 
      0.78856598384631329009983207242156E-02, 
      0.40251042064986514657905454665344E-02, 
      0.39434945748613745276711516127143E-02, 
      0.45355232247853231691788919927919E-02, 
      0.42553589471708845699257172105916E-02, 
      0.87471074213636353080888215779404E-02, 
      0.42047313560550802314674632374969E-02, 
      0.81096371219963978174715168622098E-02, 
      0.28740657651753097232044695511438E-02, 
      0.63674671313288725555458273689796E-02, 
      0.15410483366597558629474357722504E-02, 
      0.71166660878827645405484492115954E-02, 
      0.59662756163675577400108631767509E-02, 
      0.20094402952755521669759725025506E-02, 
      0.41552779846196640707863828648633E-02, 
      0.78774367939671502747274922841800E-03, 
      0.55555924933787222053248810519126E-02, 
      0.43156090349809314605947306977033E-02, 
      0.72575720454752208315574044802253E-02, 
      0.30758989792726489058477276759038E-02, 
      0.36691998090203889068724774157613E-02, 
      0.33012791980584823111731500154720E-02, 
      0.18827425749868378780107906069778E-02, 
      0.52275015728844655248822057721928E-02, 
      0.21360163256866199017108940749397E-02, 
      0.38709086959575945186011536418464E-02, 
      0.15225813039724090290900702612314E-02, 
      0.63577625651425152836180330975545E-02, 
      0.46584793041288563398185969929494E-02, 
      0.20643564927291033331571316680505E-02, 
      0.26933883371835105597852952390420E-02, 
      0.48308600971180510860882325582472E-03, 
      0.75143865994468267960236104870498E-03, 
      0.11530302612959225729996460632643E-02, 
      0.82341136968804111621935856491347E-03, 
      0.40628665291981280528823925473427E-02, 
      0.20108332814602822754605793773821E-02 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule46 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE46 returns the rule of degree 46.
//
//  Discussion:
//
//    Order 46 (385 pts)
//    1/6 data for 46-th order quadrature with 73 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 73;
  static double xs[73] = { 
      0.00000000000000000000000000000000, 
     -0.49937321552503064883212174284696E-01, 
     -0.58248078099624925941940739070475E-01, 
     -0.26752517246704927261140630108639, 
     -0.43867813853252091527068916646400, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.23976542950787972353232922301146, 
     -0.87549594799702810906321628660430, 
     -0.41617075719872941659009240564594, 
     -0.54326121142564273632982220217311, 
     -0.49334938461621505385336504601822, 
     -0.98141250190288778713761374334280, 
     -0.65605474425672036542123830767580E-01, 
     -0.16097654175693748178184943252463, 
     -0.32922143422069145595634871480200, 
     -0.23167087272607281619362898673378, 
     -0.32291103943378255878422004279321, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.79135174960542667310417788753063E-01, 
     -0.79501562539103118019172171437469, 
     -0.43520622836381734705078227452941, 
     -0.64094592211243916763929215200559, 
     -0.31732449637355273069551804224371, 
     -0.19862767910431307864204482880642, 
     -0.12403617126406631952361500698963, 
     -0.43416110602837889188093965553450, 
     -0.38864284496138298848368820422816, 
     -0.49742313896156865617690862094750, 
     -0.65506909060071668258633513435152E-01, 
     -0.54046107587222447376898336935916, 
     -0.19374430311448627449139164957895, 
     -0.21797347584325374123660718055068, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.19584532512471844401439133027531, 
     -0.63168655681204850681563549874707, 
     -0.90233030436033463513181070921112, 
     -0.60642348255765936124722001014165, 
      0.00000000000000000000000000000000, 
     -0.80972189677200294143541375522637, 
     -0.58209787586990173647256990544253, 
     -0.83686548692344391617162664381025, 
     -0.10243921312001814379582985121977, 
     -0.52704295544169988915202539886442, 
     -0.71820502549667139042476315500771, 
     -0.78211800995999025958688636663884, 
     -0.85578878348643964344785103905108, 
     -0.26080043841003711556582180762586, 
     -0.13061737230200327923828027514913, 
      0.00000000000000000000000000000000, 
     -0.44521414457403926773317497715994, 
     -0.69780908932146560679288082039321, 
     -0.66899206817769503267088297222047, 
      0.00000000000000000000000000000000, 
     -0.72891023006189965013712326310293, 
     -0.37024530768151522800672063932221, 
      0.00000000000000000000000000000000, 
     -0.32275953620672120686976828661180, 
     -0.66470762672695841262017667178171, 
     -0.32306760913108771241386793073162, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.13683702205917657289237473998807, 
      0.00000000000000000000000000000000, 
     -0.55947816484895364063783653823836, 
     -0.75690025149534103382015622011351, 
      0.00000000000000000000000000000000, 
     -0.92282620861069337654871187945177, 
      0.00000000000000000000000000000000, 
     -0.94945587611159447764356569726866, 
      0.00000000000000000000000000000000 };
  static double ys[73] = { 
     -0.50516777744773643918214212704283, 
     -0.31040114421146323467610402783791, 
     -0.38139174632823850428680564625081, 
     -0.45087881559146180751303614188317, 
     -0.38716524171664208079585405722511, 
      0.98605740376479278980387642403069, 
     -0.53179537383682808117522612094495, 
     -0.39923330820610651661243817919637, 
     -0.55924165141919532170843072337147, 
     -0.32183421694439171166523920516920, 
     -0.48706522111890310455837522078968, 
     -0.43864037323932895072146723439064, 
     -0.57360768783134360345778307830250, 
     -0.55766758615272760217178705034940, 
     -0.35256727255257676283714405416718, 
     -0.49258841625772331573291383409050, 
     -0.30296300817765420709023711921593, 
     -0.27337520447762240864886970930914, 
      0.47471537268435794928244865456278, 
      0.58035693653140731668646213878644, 
     -0.48533267753756863676915233617052, 
     -0.48917333384362543356099193333710, 
     -0.55808744003610720897481815391397, 
     -0.55694713731397587210816279802739, 
     -0.55689799448943849551272945603815, 
     -0.48859325854023263980525452506683, 
     -0.26117953849550246736508673038521, 
     -0.48398948921467227547423443087059, 
     -0.52999026102693285956849576671356, 
     -0.52583495263659649090093095454021, 
     -0.57355483546767908414270209478254, 
     -0.55557204112323420607295734852969, 
     -0.55739761697092742823342363228578, 
     -0.21246252441995259506958905076888, 
      0.36139661177510802837475693850211, 
      0.11582303928093664478281478578489, 
     -0.57359024620476593895960903232245, 
     -0.48075999227146968059407118917937, 
     -0.57371658547584443139058969065517, 
     -0.52698535358924523864722032814143, 
     -0.11344215054773373022279305101957, 
     -0.55728798211338842288776835769559, 
     -0.42586101506138965898586810863037, 
     -0.57375019535379531179673731067538, 
     -0.16775873288237485180210443333935, 
     -0.36463881029295166822714518212865, 
     -0.48207876525200853027431004337685, 
     -0.52747840671084502211591497419342, 
     -0.53180888262584484912774717875574, 
     -0.52832053993193271018974971786141, 
     -0.52863298211428052722594126999808, 
      0.10481124042981156776747972952527E+01, 
     -0.57371140175475285675686432495456, 
     -0.52448444053919330486337551923650, 
     -0.42468204214393876878868764552758, 
      0.86891668006615384049718425264833, 
     -0.55555967179393677106465759845972, 
     -0.43207267925178557569834574602145, 
     -0.44163717527386104241934308216581, 
     -0.36120487686156226967665052805572, 
     -0.57348736410677225487636160697833, 
     -0.57339310967499044449513445911666, 
      0.23111663125846296538268627177402, 
      0.70977328433019335226869776419213, 
     -0.43325524227150080894750839759120, 
      0.00000000000000000000000000000000, 
     -0.57318116904317379016261023974720, 
     -0.57316379219729504712466298101902, 
     -0.22231995200693748608468315007358, 
     -0.55502955244943758240255763984666, 
      0.11100839542250144489340222725729E+01, 
     -0.57286876316098503074069033064773, 
      0.11529930056431959330783131719851E+01 };
  static double ws[73] = { 
      0.99124219640185234470989336544560E-03, 
      0.58614096060962842645487250845942E-02, 
      0.56854522710881067383564039723551E-02, 
      0.40625620962929522664929908767625E-02, 
      0.45144310543124329374817531802773E-02, 
      0.82088088596127442161702412322479E-03, 
      0.14407828064975009824715695382818E-02, 
      0.49783947422097924090273775604626E-02, 
      0.88233634287350469355364288087551E-03, 
      0.54934439723781169400547752605367E-02, 
      0.33669998728502865263960539461844E-02, 
      0.41903624852692709261378072222874E-02, 
      0.17334944288966619414269024913119E-03, 
      0.21658266572394104793711645431275E-02, 
      0.60479958760001231392858259471875E-02, 
      0.37756322991798540494330453652444E-02, 
      0.62207318637145172822591743503502E-02, 
      0.63286469162156006735700584173533E-02, 
      0.29826294052344452485495959644700E-02, 
      0.27059360522727772191971339645741E-02, 
      0.42019177882590680248910223133272E-02, 
      0.23015042374454809666901890875641E-02, 
      0.18828680455513843606219753257332E-02, 
      0.16690125284755693166888299832115E-02, 
      0.20781252759372409612923667460460E-02, 
      0.42317638788900425084313484715466E-02, 
      0.71037803619188976533554128557295E-02, 
      0.40307467270170115306414668368486E-02, 
      0.30522630792848607523818935976199E-02, 
      0.30218379230350043130768695564237E-02, 
      0.96574870847494658906685789105858E-03, 
      0.19000953811516758513791103731380E-02, 
      0.21829276451288062623559064112356E-02, 
      0.73314920756006735071788427464883E-02, 
      0.35114723476610777601133441893946E-02, 
      0.42611988684563371970090214166312E-02, 
      0.94946505946997602521228000171027E-03, 
      0.34509142883519454104321786989457E-02, 
      0.39953322686081774873355011018233E-03, 
      0.27638699271507528771931724478937E-02, 
      0.42958741400998447897262423018398E-02, 
      0.12717205008773558943510641599495E-02, 
      0.42175860264852262389840957601069E-02, 
      0.51628931464429952176360586413594E-03, 
      0.82354400118366882096718549841164E-02, 
      0.49763923773181868425612339373364E-02, 
      0.30353430462986957072339039682070E-02, 
      0.21010204804081172094616080287391E-02, 
      0.16435606239301717159637463528701E-02, 
      0.33633961349195807382580363482430E-02, 
      0.35248211978801784050294560000961E-02, 
      0.63942094224216359477928434782097E-03, 
      0.84411673200603475386814798360338E-03, 
      0.25146454557737862996260311126641E-02, 
      0.39858155043769138905211326067054E-02, 
      0.16047534667867884438074054540413E-02, 
      0.15714180013667609459726436385097E-02, 
      0.52508538616460060672169443423580E-02, 
      0.27559932699033631806510178138074E-02, 
      0.62831477470730771167890815345462E-02, 
      0.74766335440469584911791085907174E-03, 
      0.95979934038783073189952725536941E-03, 
      0.39391900179981746067668114449212E-02, 
      0.24830895862025367863318765158069E-02, 
      0.57901415438434356432649769437734E-02, 
      0.14553688489280240366387824569707E-02, 
      0.89356029388787655872574198511159E-03, 
      0.69941978699045079847532710585071E-03, 
      0.42651637563773055644827276665484E-02, 
      0.85219038226410641706285284477815E-03, 
      0.29568676410262965948870235894266E-03, 
      0.34392691859187999137468851025595E-03, 
      0.92751764638845732676762722810718E-05 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule47 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE47 returns the rule of degree 47.
//
//  Discussion:
//
//    Order 47 (399 pts)
//    1/6 data for 47-th order quadrature with 75 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 75;
  static double xs[75] = { 
      0.00000000000000000000000000000000, 
     -0.96663742233248238790417222682341, 
     -0.55421586636589606362103864964496, 
     -0.61264177937586419170392244079506, 
      0.00000000000000000000000000000000, 
     -0.79684085415639357802005234937672, 
     -0.97128586885482393995603476389312, 
     -0.93648904794849557432223757784636, 
      0.00000000000000000000000000000000, 
     -0.53945273449791007382233425306043E-01, 
     -0.10288337054182960580026068100176, 
      0.00000000000000000000000000000000, 
     -0.73248213764906407079372614104364, 
     -0.68426811117296314009035410342951, 
      0.00000000000000000000000000000000, 
     -0.11259199413113422069683485156813, 
      0.00000000000000000000000000000000, 
     -0.87984784872866267520426458088785, 
     -0.23123809819906032963749193832378, 
     -0.46889550912339936686025371549248, 
     -0.11739233105524143722624127853417, 
     -0.10992961573914339919329881126213, 
     -0.16204920039724447196768977188244, 
      0.00000000000000000000000000000000, 
     -0.56792554236958167052318250024853, 
     -0.36310842599378440119880946169444, 
     -0.46903656396230065793902775073785, 
     -0.33803039308384563532443360051617, 
     -0.12603730801931571541242657552702, 
     -0.44018371642510357125560685343620, 
      0.00000000000000000000000000000000, 
     -0.20847045004866279897685055766704, 
     -0.22797318474342093841199473155967, 
     -0.81969334154814590110268106658620, 
     -0.24826155062103396453153632431776, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.12427095901943726054227825242898, 
     -0.26070189438266727648706156956545, 
     -0.56248468735254147016623038597020, 
     -0.66071472756340151327417396690268, 
     -0.82669686777638652733059075709640, 
      0.00000000000000000000000000000000, 
     -0.88966644367593888841875654650478, 
      0.00000000000000000000000000000000, 
     -0.34452655299507786176572793372419, 
      0.00000000000000000000000000000000, 
     -0.65494542725425338003061612036810, 
     -0.93012524112428256602193934913638, 
      0.00000000000000000000000000000000, 
     -0.66785089424544850131055902834174E-01, 
     -0.41597899021645575618793010361463, 
     -0.24627705576691568442688105699070, 
     -0.74973726200095632101264662967133, 
     -0.74550793444160419487461731848021, 
     -0.45718177804191062450002601345813, 
     -0.31237691698758109357216488302742, 
     -0.36244609664691797249270623595478, 
     -0.32616646260370706311727811093441, 
     -0.53784292001756868772192703221020, 
     -0.55851533886826749020734059929719, 
     -0.66003657270757138896159641690619, 
     -0.44695335875345503292611264014711, 
     -0.19870749346215743992312176932942, 
      0.00000000000000000000000000000000, 
     -0.80239000604949848720107906069622, 
     -0.13449553992786908006612921402555, 
     -0.71709284775736678626288283834417, 
      0.00000000000000000000000000000000, 
     -0.51015724874008600140209964145657, 
     -0.87380123018010155280547358512783, 
     -0.39184908288647526514367672290932, 
     -0.61920375973213862119233932838043, 
     -0.26594843220537421747193002664485, 
      0.00000000000000000000000000000000 };
  static double ys[75] = { 
     -0.54885817731744273760020596669005E-01, 
     -0.56535389802881173311914125214849, 
     -0.43595283427447394280235040784487, 
     -0.42805378386933825639851755340423, 
      0.10543613068074012401800062516168, 
     -0.49023992523710045116805254072657, 
     -0.57553352291636312128177900316235, 
     -0.55999425828931551781156236993222, 
      0.97965289344360951855945358784080, 
     -0.22457611472173417062464887127478, 
     -0.30138835771850783413933825956514, 
     -0.30545083028762013017508143191825, 
     -0.48491296481354098368779641623556, 
     -0.42961315180326213617646939344682, 
     -0.48945570323521847726831699400261, 
     -0.48818031348608414902160574354381, 
      0.87905212764608203798450778699722, 
     -0.52962597960343054791384488347290, 
     -0.37175232717159904337497914590613, 
     -0.43001929575977665531448764589702, 
     -0.37274742279795921566656912482601, 
     -0.13964656636214571305124861250780, 
     -0.22098784556599346718036712899788, 
      0.10726698183580133144561578611120E+01, 
     -0.52773439710278615552647179725219, 
     -0.52763898191178919130394070654322, 
     -0.52697835910181444235797110038292, 
     -0.36648653111336403999140387539267, 
     -0.52894478654176106374465306810932, 
     -0.36579746792427023484428727362546, 
      0.24344631435472381136691814975209, 
     -0.29950682881818228902768609651588, 
     -0.48721710760636867571802595633338, 
     -0.52883882159690679649594382117023, 
     -0.52822203489770240968743803141947, 
     -0.37667277468535763242146591083137, 
     -0.52940747517054265256495756528320, 
     -0.43564261027324440083824439259671, 
     -0.20992358276179358318845033026211, 
     -0.48651251086220905025316938031502, 
     -0.52681305116132787774755531194190, 
     -0.55712806420844685964294934833455, 
      0.59300989346090695220519494616955, 
     -0.55768637399333409575145201865462, 
     -0.13980368807481589507662723468431, 
     -0.48599981614691996341237894437611, 
      0.41840134446650681930079926553608, 
     -0.48392971086401977797578179121298, 
     -0.57383958132479530247972005058149, 
     -0.43820483878434062239610080574033, 
     -0.55764972588734990620416241805372, 
     -0.29204731777806559769068818139152, 
     -0.43448156854177003457395184650007, 
     -0.55679770960375493548218080656374, 
     -0.52714057557693106272153869645905, 
     -0.48449632102066625308448025105653, 
     -0.29161853684292105771088446811658, 
     -0.43162684421663934029018768245950, 
     -0.55705385397930168896131546768592, 
     -0.36980700529962546759999533399443, 
     -0.55689739217360638744554116629651, 
     -0.55674192040916277835962510290896, 
     -0.55682962368161313519342300149308, 
     -0.55738068167530312022049700934039, 
     -0.57360394426682101478744320464048, 
     -0.57345366439244921550225256351043, 
     -0.57357771856416956297048773459624, 
     -0.57343681973260213887282500888328, 
      0.72177614769930334697640951850176, 
     -0.57345497171610556737079624481070, 
     -0.57352108080599339578631706071895, 
     -0.57345327364010964490936936951959, 
     -0.57342883815744037296971547838128, 
     -0.57351635407928127465513028410519, 
      0.11456762187114260229690797261936E+01 };
  static double ws[75] = { 
      0.28392039488336948310179896994156E-02, 
      0.28003686921136784757704775202544E-03, 
      0.29631179312114912348720642679995E-02, 
      0.29809999165091017111077408004804E-02, 
      0.30922285699322994610269077138085E-02, 
      0.18460959404365729590181431759235E-02, 
      0.12872508972618665537965952389878E-03, 
      0.53151358172545237088027219688402E-03, 
      0.81365776381310603141873461966954E-03, 
      0.67825918391732126919840300376994E-02, 
      0.59347976694683184123957244470317E-02, 
      0.29716933960007650620912395217221E-02, 
      0.25197407852253980023714633624274E-02, 
      0.34300378782082533220752669591608E-02, 
      0.19375775064770987319833820229901E-02, 
      0.40484506484434265904918087580818E-02, 
      0.13706927206297658573031048360427E-02, 
      0.13592574447542674909853690088402E-02, 
      0.57306616583541857081486048497868E-02, 
      0.44642387582970339954371597118078E-02, 
      0.59659697355700592641664212659135E-02, 
      0.70669480375351196705523756968934E-02, 
      0.65951783851403357899976425662021E-02, 
      0.43628370961428043795773932668710E-03, 
      0.25708114420998187417114840014581E-02, 
      0.29876897524925467219456452120523E-02, 
      0.28146519351937002714666355261641E-02, 
      0.56303412225285581963359647583548E-02, 
      0.32994230091903712639399659211836E-02, 
      0.52297242356955663867894660142373E-02, 
      0.36838585883225044871577564370140E-02, 
      0.62387352632335361474645491103618E-02, 
      0.41599811873456112930679098763696E-02, 
      0.17194415272340685412270563416839E-02, 
      0.31880794971043660654134886269745E-02, 
      0.29952623965446603243032973139989E-02, 
      0.16498839143299521006632241989400E-02, 
      0.54403322566293624509259954994073E-02, 
      0.67187331124535844522324030229043E-02, 
      0.35502719665898356713324365688902E-02, 
      0.24754883633746269136218392735926E-02, 
      0.12031798815691044211245980202696E-02, 
      0.26745453014501257900256883065918E-02, 
      0.92265422131087006168788486251726E-03, 
      0.36225684350627313248590370760841E-02, 
      0.41984151428468319683013713735675E-02, 
      0.32537703172163361922509390986613E-02, 
      0.32030949415312354270600565161129E-02, 
      0.33013381207853754356962435323648E-03, 
      0.26754401539508943172771916930207E-02, 
      0.22515566984635737794357729457255E-02, 
      0.59591990991400095363194619176182E-02, 
      0.52750078810945658157831528545914E-02, 
      0.14733980267684772821583690219304E-02, 
      0.21789567303724062164154906204501E-02, 
      0.40463448470928386904872275539629E-02, 
      0.63567301660162872416841128607952E-02, 
      0.51309743375347267914123716886908E-02, 
      0.21660630347248746919200853426513E-02, 
      0.48851415391794237971659363471398E-02, 
      0.18669504929149491480935586903022E-02, 
      0.16964959923743123037357734611633E-02, 
      0.20550364827452665795472373860314E-02, 
      0.22311032805628682468829744578561E-02, 
      0.49275841826369605327759276116059E-03, 
      0.59592130393246716586641585985888E-03, 
      0.98092238248484605161450626802125E-03, 
      0.70060956245054444684115888791512E-03, 
      0.23900824377240423856080438652344E-02, 
      0.86547617602673961987264272734724E-03, 
      0.47650966908745127742574114884819E-03, 
      0.92916397453096297560102228965669E-03, 
      0.79211481322285497661392660284067E-03, 
      0.96420760486050977248262575972495E-03, 
      0.58732212176584084295881188275225E-04 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule48 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE48 returns the rule of degree 48.
//
//  Discussion:
//
//    Order 48 (423 pts)
//    1/6 data for 48-th order quadrature with 78 nodes
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 78;
  static double xs[78] = { 
      0.00000000000000000000000000000000, 
     -0.45164772650050731682068089766497E-01, 
     -0.27473470409127158336658796332646, 
     -0.32601858877976093934124855013368, 
      0.00000000000000000000000000000000, 
     -0.99490990154620226158720849494482, 
     -0.32465958030201586136906707429503, 
     -0.84012618143751376567667347453966E-01, 
     -0.68050562893254251474678136205112, 
     -0.55412485218420527438896911621850, 
      0.00000000000000000000000000000000, 
     -0.33038679441849822366871793787534, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.72722633108802342397155155606768, 
     -0.23810041415405527504557089888405, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.76837647852811595992367941824138, 
     -0.39481835688065754917935872570598, 
     -0.35808200103666623461734551126576, 
     -0.70743720676693189266156017090859, 
     -0.97010789493490351143718625991234E-01, 
     -0.22510606888237807083444019159841, 
     -0.89142246128908740222364354150092, 
     -0.11231501558393124954330372677768, 
     -0.24247749198133974832677108412529, 
     -0.83724014027395502676007507010710, 
     -0.13318347049091707247482863126316, 
     -0.72689161196677377776022254195996, 
     -0.64034267680893964392473256901576, 
     -0.68321717220184284728388581404162, 
     -0.83839805257724055363465708705911, 
      0.00000000000000000000000000000000, 
     -0.59225145333921036546705814041242, 
     -0.36426477146050005896562735592106, 
     -0.64466682188081663940279587210100, 
     -0.18326940651609773666251309353216, 
     -0.64024696057757473624745250904559E-01, 
     -0.46015923358724641317739926182026, 
     -0.62850907517187595403082243802386, 
     -0.93449067868256757675651398163129, 
     -0.49718634910646430159295459732375, 
      0.00000000000000000000000000000000, 
     -0.76477834359061656278550903563318, 
     -0.86560964360183956199073793066554, 
     -0.80347270468749417960796064245399, 
     -0.12053965097836346999960081953107, 
     -0.61721699934178919899314762692642E-01, 
     -0.94600280954375241238529550780606, 
     -0.20973158081832926767896749904045, 
     -0.53588808819149533597940929507170, 
     -0.89909074000354327737143587267527, 
     -0.14226265688360420003675348239732, 
     -0.58106813824728181662821714714662, 
     -0.43230361497566010923448640099201, 
     -0.54615078926623396814736475110416, 
     -0.79506240679562089185948571554395, 
     -0.25184574558858296991088593038170, 
      0.00000000000000000000000000000000, 
     -0.97712556750745254153359066235566, 
      0.00000000000000000000000000000000, 
     -0.26607245123060807138658601473649, 
     -0.12724847472110119640623202948370, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.30727380337735402963322838421622, 
     -0.58589180790116215572933101719513, 
     -0.17860222460692283965662443246084, 
     -0.43751879768235861761962323575309, 
     -0.48056445088740205085548935531431, 
     -0.48038610507302961545832372554540, 
     -0.41788380954447150618204025259598, 
      0.00000000000000000000000000000000, 
     -0.36944533995870546586326908358313, 
      0.00000000000000000000000000000000, 
     -0.72247744117483068359801481921580E-01, 
     -0.23245589269716346191888401881117 };
  static double ys[78] = { 
     -0.49852320663764975836658715881651E-01, 
     -0.20374579426251387086718546231968, 
     -0.38976656878879817546660173662947, 
     -0.51725511617300734515148746244054, 
     -0.13890333151925469882109904912617, 
     -0.57319686228940163286673996687428, 
     -0.54414549714157621787319018880544, 
     -0.98677996812038879657787374589079E-01, 
     -0.57450768223436031703493865457993, 
     -0.48499178896207101508248204799370, 
     -0.54041966741780855485936015475384, 
     -0.43371331113028682857825644993170, 
      0.98700017550217162316302585294308, 
      0.87701971609267538163224358348922, 
     -0.53696054670467352362246242651342, 
     -0.33709214277689455129828382991245, 
     -0.50491682572683646092981172306082, 
     -0.46179620632104391793324280621199, 
     -0.55909334619049486411241922583373, 
     -0.55942830548805619506092139545342, 
     -0.48419731762996786439303530812544, 
     -0.45652178474401761799466199325239, 
     -0.55049949314164685016946873588734, 
     -0.50872207349408828698853958022725, 
     -0.55689353511683902690020296113291, 
     -0.51622984675368820169968858727220, 
     -0.46526044009365978341053361548518, 
     -0.55929784568894692081450974794614, 
     -0.18399301337776744727308882527111, 
     -0.50260381043912900882119219273873, 
     -0.49252317057731876434255148922154, 
     -0.56097451023259417695151253611004, 
     -0.57406923516192628548625680508365, 
      0.11141547030632238086479160754135E+01, 
     -0.55651660365661076171303549417789, 
     -0.57405163423686233358717892168751, 
     -0.53196149015755170375584796642517, 
     -0.41372183327215987269752027182733, 
     -0.41051333470792190295409628981832, 
     -0.48470232197639096590316447145498, 
     -0.43960983180755192818725222543803, 
     -0.55824444146228311378010225156953, 
     -0.55624309454314381245877434541105, 
      0.22685299762337478108310204018125, 
     -0.57369006779153222884968212905112, 
     -0.53020559443910391226024877057732, 
     -0.53187309456333599368352110014477, 
     -0.46914022619776744572721261875474, 
     -0.27636104453730623282084492698214, 
     -0.57380107877191676765658650529341, 
     -0.54463175727474219404497175590015, 
     -0.43425004501953381185987824757432, 
     -0.57340126819590793499251738305266, 
     -0.56908448605073339254550804081805, 
     -0.37608188392552828091032473401365, 
     -0.42921616619109235112346077631242, 
     -0.52577535892830112703001792903022, 
     -0.49030858978345230828266827255943, 
     -0.20270732080067755549950587581872, 
      0.41636790443295163120820765376968, 
     -0.57333501569761494464853713844014, 
     -0.34742627962761774860725212273245, 
     -0.56621285723554951319184951691942, 
     -0.34679505484408231261353167027046, 
      0.60078752240377683267001545331127, 
      0.10574851563475822963976796065719E+01, 
     -0.28592841602633954760583222174030, 
     -0.57323206949747451195434362599767, 
     -0.26683655423740594471342520417422, 
     -0.52805029804391126423728142010436, 
     -0.57349391506508209544925398288834, 
     -0.36887306566949226352897864054757, 
     -0.29271463211375329134525196252490, 
      0.78041430345998157811281878081512, 
     -0.36277089606217184970553665226938, 
     -0.56756638958566119435478459622114, 
     -0.57707312306547335610792852443827, 
     -0.57716232573911703125929908613260 };
  static double ws[78] = { 
      0.28724228383926440358132669094431E-02, 
      0.48957050616893594998704796075844E-02, 
      0.39911404555441563247326888994291E-02, 
      0.24308030008482638440637786996581E-02, 
      0.26118640649911296596022610648512E-02, 
      0.31286906542800058837727462589604E-04, 
      0.20332915152781863944904896111764E-02, 
      0.66223356940445958324045825592685E-02, 
      0.51160504059193445246422419150521E-03, 
      0.30177173547108975416398128210504E-02, 
      0.11990187952425307052673011943937E-02, 
      0.42890959944533322728503697326550E-02, 
      0.79717712193663934454694268352422E-03, 
      0.12706395119960644037665264677177E-02, 
      0.17101518366924235948036417308251E-02, 
      0.53225083732062422101624101122955E-02, 
      0.17372337371711722394712009961994E-02, 
      0.21611199292316263486466123355414E-02, 
      0.11377184875874711922363652238389E-02, 
      0.16661121086858607539082516065793E-02, 
      0.36094495784172349233648222679572E-02, 
      0.27840329081828913602050362052188E-02, 
      0.21361914923668763800444433760586E-02, 
      0.33445043203334672620049150278461E-02, 
      0.81943436438923681400473901191879E-03, 
      0.34768574704553087820389981511291E-02, 
      0.42973989829883201728563607299665E-02, 
      0.10074738530685819860190968486265E-02, 
      0.67963853132839489825520666772776E-02, 
      0.24418950258761731688970926384011E-02, 
      0.30003987104228644513285637413249E-02, 
      0.13513377123312761519271287419019E-02, 
      0.43804618210474632627816499919792E-03, 
      0.21250973370504179328882832139026E-03, 
      0.16000160750664913858277795560230E-02, 
      0.77557764662015942773935106392613E-03, 
      0.22824221568478339384080338386456E-02, 
      0.53299353429123467214358800295839E-02, 
      0.56252060523420948976026938340142E-02, 
      0.37291069649491200558250666737351E-02, 
      0.36110603880506100367942684638801E-02, 
      0.62253124627782549500409494910883E-03, 
      0.18073877946168386060665209809962E-02, 
      0.35738294519148201788176936416809E-02, 
      0.55709543451833858845707308928538E-03, 
      0.14035418542487990820244339057528E-02, 
      0.18193421280631652518829218910697E-02, 
      0.48148722968830507266694307528103E-02, 
      0.68567519968462246935420491727088E-02, 
      0.27283604140411874997242863942056E-03, 
      0.25772891742822183465897792887364E-02, 
      0.42231795512038431505510095354369E-02, 
      0.41379361426485253173206620037588E-03, 
      0.14194049785849893147037237210493E-02, 
      0.47384283794563909094585996951638E-02, 
      0.48183988814673540748933494366651E-02, 
      0.28449733997703946018961571694460E-02, 
      0.23703233612168300918965933073201E-02, 
      0.71693912212143004097243510739388E-02, 
      0.33905755163797439812159959659759E-02, 
      0.18396767112155352948407062614303E-03, 
      0.32855821621882024772752035184439E-02, 
      0.15742892940985067189187325607995E-02, 
      0.65179431778550066437367125243887E-02, 
      0.27923790445661769025995748763738E-02, 
      0.58407401612327712992571037001792E-03, 
      0.66530250807598734354239162260138E-02, 
      0.79700265904184871610937261087522E-03, 
      0.70512484634479567571833106295500E-02, 
      0.31132194356885414154457811749791E-02, 
      0.83725477449697097193875164230539E-03, 
      0.54179864271945087535990988436145E-02, 
      0.62428841652845444519001502321154E-02, 
      0.20885623095275804800293379629396E-02, 
      0.60569245588953023785142088907429E-02, 
      0.85221050803831089796247049975189E-03, 
      0.30928913568865192645879463322352E-03, 
      0.31569151523254203691074928546368E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule49 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE49 returns the rule of degree 49.
//
//  Discussion:
//
//    Order 49 (435 pts)
//    1/6 data for 49-th order quadrature with 82 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 82;
  static double xs[82] = { 
     -0.61057213632871754777556178617093, 
     -0.21740746844089577030217736555995, 
     -0.74657702398083537596708720880628, 
      0.00000000000000000000000000000000, 
     -0.18652159062065316394708803878930, 
     -0.88761825157466116267033225250300E-01, 
      0.00000000000000000000000000000000, 
     -0.48639838036980400557648362218933, 
     -0.37022598751812238568894342780029, 
      0.00000000000000000000000000000000, 
     -0.55660589033253986016354022616185, 
     -0.23093959523489028823297605424361, 
     -0.12830890012394010891659450444936, 
     -0.63417236502579470761226709395580, 
     -0.36128311243234453339020706484578, 
     -0.24458985250678919382029730465663, 
      0.00000000000000000000000000000000, 
     -0.76440468282859730372500804391987, 
     -0.43881396034292203335647632860058, 
     -0.29031421623299546923298899864555, 
     -0.11022106138867966101881171838643, 
     -0.82697827385445359866410833664403, 
      0.00000000000000000000000000000000, 
     -0.22649962118473459738073399292893, 
     -0.34409570678022498852315838982214, 
     -0.75185639971866986636651769262636, 
     -0.53461908075309836850101213437195, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.67842852868361945529379240804657, 
     -0.98128049170705739598259613001700E-01, 
     -0.81198101881640100102568599136862, 
      0.00000000000000000000000000000000, 
     -0.38816954030832809949204934487173, 
     -0.51214253977457617025282281965566, 
      0.00000000000000000000000000000000, 
     -0.89378361429922217289759901811815, 
     -0.20438937977062137351398598793386, 
     -0.29077462427214531822794640367537, 
     -0.45530327016209634617378737911148, 
      0.00000000000000000000000000000000, 
     -0.50049811961058966626537547990924, 
     -0.59653316717651885600358396085261, 
     -0.18241222915476293063943966605377, 
     -0.17504275018488792290118428567192, 
     -0.87027029250770403203069802318185, 
     -0.25906638561252602666842654816923, 
     -0.25389417516222969869856982682409, 
     -0.63924734492284199961203241117988E-01, 
     -0.32336677856616207857128015137548, 
     -0.45933867955078441295151086540347, 
     -0.73910550014397140051390907061818, 
     -0.11606977913379114755850687691385, 
     -0.65577583679224655444244929369890, 
     -0.83552616064107205337213906126893, 
     -0.67748842104198500405815554991182, 
     -0.37231932331276101094567014835110, 
     -0.93676060044306703638748111096616, 
     -0.66619173394165069756523303173151, 
      0.00000000000000000000000000000000, 
     -0.80515347766473982936691157759254, 
     -0.94025356560109772303978498971830, 
     -0.55800138185704783040837528977761, 
     -0.33515583932564482010856277294034, 
     -0.72304850966873315163769622529180, 
     -0.46826261407944706264024812254397, 
      0.00000000000000000000000000000000, 
     -0.11247883641677583766316101770443, 
     -0.41172923406071229009367773338534, 
     -0.89091288563351724918827721227521, 
      0.00000000000000000000000000000000, 
     -0.70999559860544022018934170422900E-01, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.13798316432331004899182819510688, 
     -0.57133871576000905289179002614761, 
     -0.97448997895894366957466437859893, 
      0.00000000000000000000000000000000, 
     -0.58709342896187276256324804582086 };
  static double ys[82] = { 
     -0.50805411191391851523789864972562, 
     -0.56234652242839114914654701128374, 
     -0.57509668203848429431038797517268, 
      0.53679828476119178289324025007759, 
     -0.15190921808712561549658969393441, 
     -0.14305301576398435281760117117876, 
      0.90237311257497758569938319452643, 
     -0.57473532714413321696267840291545, 
     -0.38858551748024275977250958169595, 
     -0.33135871726650782696587314579878, 
     -0.40073527072101794974371578083120, 
     -0.54310366079605492812547381587275, 
     -0.57500475146778529513914433479050, 
     -0.42554839571165244196875402232376, 
     -0.44153050231697094128484412680053, 
     -0.46866793891867074878457016677705, 
     -0.16570704844176947394320757260228, 
     -0.56504278301868705651681263115196, 
     -0.50652589028067622329647801023762, 
     -0.52106748079988378853417843443633, 
     -0.56223112097793733954156361308556, 
     -0.57462270717388530892162623588585, 
      0.45932243496595672349082629601422, 
     -0.22828624984036516968491608086478, 
     -0.48829139755087173802506438507515, 
     -0.54439150942275613023935368560183, 
     -0.50670327829951912727787293959844, 
      0.99438807073001743653229293189603, 
      0.12536668432105287250159280886041, 
     -0.52796151262664728760570958076544, 
     -0.31111906291100543384916079312981, 
     -0.53473260414787373087982349525345, 
      0.11171423268477136498973199720050E+01, 
     -0.53658558406526073583077834124086, 
     -0.33881791091529569834478815622842, 
      0.79840597283762599147168675735209, 
     -0.55834047448625304405740354168855, 
     -0.29941559804457568276715427694835, 
     -0.34707538263885206082450977546688, 
     -0.56218624883920679775484228539619, 
     -0.52252588465952329020710730711399, 
     -0.54142101949433652578422629700355, 
     -0.54312444348452881138890792354339, 
     -0.50297724458337249430337647089853, 
     -0.37315461111231110867297450418990, 
     -0.53339994853211595567498621491810, 
     -0.41535262741815023278680392917739, 
     -0.57390850498838587799366720382219, 
     -0.38909332692463915420277441179215, 
     -0.27441054007776771321877791597838, 
     -0.46136190485467458900510307928511, 
     -0.50637232744855604509324845495701, 
     -0.53497319519919178836708838160351, 
     -0.47916361999904640273373623874426, 
     -0.56063167758821432142013730704571, 
     -0.55812525411705885484285303911831, 
     -0.57393107655513157172595640909445, 
     -0.55819574918551185878858176169832, 
     -0.57308501713419226805469571310331, 
     -0.55400636586589158164274508913012, 
     -0.49699174948588019440752603403060, 
     -0.57375002851113489476058869278243, 
     -0.45895888627009702143807907926916, 
     -0.55885836769979743482431579044375, 
     -0.45754555590751383576718554173853, 
     -0.40298677651135379575927050935698, 
     -0.25079797310498966721630314170512, 
     -0.22883057799614009961932750289238, 
     -0.32888195428957532835127033424983, 
     -0.57366938316409890583305452094740, 
      0.70290670754518056713385541693329, 
     -0.48828388229492647796479685307693, 
     -0.44375651708543360206049930519612, 
      0.10591815381653453969290932401725E+01, 
      0.35157474889658538587578359215421, 
     -0.63699436116132156076951185465648E-01, 
      0.11470364861459781695524661983521E+01, 
     -0.44198593821724404560352142656255, 
     -0.56609726595635302503750770819831, 
     -0.57372975226826569608034581860957, 
     -0.57315719828247563816289202130035, 
     -0.57710801847006901827791738420396 };
  static double ws[82] = { 
      0.19977337356156323426514692597692E-02, 
      0.12687431336688509410843998391028E-02, 
      0.37027790592768656352817137656994E-03, 
      0.21253790361811728001071706475289E-02, 
      0.61674141211813126477419858184194E-02, 
      0.59458894649162451804820476952354E-02, 
      0.10391008318227642432772429159470E-02, 
      0.56119170151050531358685160135421E-03, 
      0.42730773131837205407382425934036E-02, 
      0.25866812740239429576418384733441E-02, 
      0.35683371351746158116269879423198E-02, 
      0.21488306719133579150142490513559E-02, 
      0.63043916628928240061101412753203E-03, 
      0.34221501813368796618088746749049E-02, 
      0.40453583847800237927752889872429E-02, 
      0.39321920602677497806251470735232E-02, 
      0.32874602491734720571848038384576E-02, 
      0.91150147690764261925763826535926E-03, 
      0.29798409559495009891833524489293E-02, 
      0.28537247963181996687266646898085E-02, 
      0.16586934600140952309818245115017E-02, 
      0.39735960174484489743453114770676E-03, 
      0.26834215475038735697529994323634E-02, 
      0.60034335233209257724964806936636E-02, 
      0.34785715040962553678932832030749E-02, 
      0.14890850096270047389519178852968E-02, 
      0.27754331998416555843944854095898E-02, 
      0.77107891245567947444663501212492E-03, 
      0.33454955156733941459045360333934E-02, 
      0.20541529314760399047675461654284E-02, 
      0.60945428068954117127542229222828E-02, 
      0.15330127901593561414098039723028E-02, 
      0.19148445202446191652573211720709E-03, 
      0.25108686300732809342950959448632E-02, 
      0.47651564692312408698835096361491E-02, 
      0.16240239076313059324667223643076E-02, 
      0.79583774233079393748085634291181E-03, 
      0.58849276406227123448983694298859E-02, 
      0.56877909758060890782274443197947E-02, 
      0.15575024842219775881243735068023E-02, 
      0.16105816942358188547189746892777E-02, 
      0.22502840437046037026372796523992E-02, 
      0.20084168279399749118935611766491E-02, 
      0.36230228738241653515806389685717E-02, 
      0.57423299789474669137183051012294E-02, 
      0.12717030464342037309770710927106E-02, 
      0.49768031950563028511667619620966E-02, 
      0.79434188752183035723354695167328E-03, 
      0.56001107097121521065256828080448E-02, 
      0.63817130650436342049602870849597E-02, 
      0.39817711545630631406587510631784E-02, 
      0.24087654021930613127612933755615E-02, 
      0.30418302126603186311999594057584E-02, 
      0.30238355103351015420463768411162E-02, 
      0.10183025304896036412292516963758E-02, 
      0.14596718890771209157014286627469E-02, 
      0.77875785285251059052497877539906E-03, 
      0.61372358148713085688282377971969E-03, 
      0.66600519529457007903113678794478E-03, 
      0.11077863336612002623045984707220E-02, 
      0.20349205397800009217113402131491E-02, 
      0.29455943798576609790138160596699E-03, 
      0.38905584979980386058934940581983E-02, 
      0.18807248411144958164960089336880E-02, 
      0.28241228846436807197405817934160E-02, 
      0.47178192432261121880126220071801E-02, 
      0.35306465870984426703979784223797E-02, 
      0.71762038572211377060639788033076E-02, 
      0.57665964607857117382384704479685E-02, 
      0.40426442422356089657524686097409E-03, 
      0.23026478604784709918110573713097E-02, 
      0.42616806207016132242718913366787E-02, 
      0.26058395691911076244778888792290E-02, 
      0.54745447246951839739615843265097E-03, 
      0.35299796865257136122691199548024E-02, 
      0.40316437640725431828744499609816E-02, 
      0.42516711921106531477513916666353E-04, 
      0.54029876938551540928891230894040E-02, 
      0.14060391641463104982472631429696E-02, 
      0.18730351088069459968645316638647E-03, 
      0.51889848051624996463295221022567E-03, 
      0.21130683065318872697498930840672E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule50 ( double x[], double y[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE50 returns the rule of degree 50.
//
//  Discussion:
//
//    Order 50 (453 pts)
//    1/6 data for 50-th order quadrature with 84 nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double X[*], Y[*], the coordinates of the nodes.
//
//    Output, double W[*], the weights.
//
{
  const int n = 84;
  static double xs[84] = { 
      0.00000000000000000000000000000000, 
     -0.11243137017724544706384538604138, 
     -0.92568822098610615625930541796974, 
     -0.29692106117972817002503164575225, 
     -0.70777454997646828139185546400833, 
     -0.90907207527032819051622115148967, 
     -0.17776222920358422699613435141006, 
     -0.63358007216059497217331539880521, 
      0.00000000000000000000000000000000, 
     -0.48639150346057012854861946724002, 
     -0.67960972311168179219672603013974, 
     -0.80987858028162034515744403376905, 
      0.00000000000000000000000000000000, 
     -0.10197291667485080764243827190062, 
     -0.63905313551409030298179847776955, 
     -0.54848856966069179497342801711918, 
     -0.28743291564930093738030476563064, 
     -0.25853575380565982943550456940571, 
     -0.85293311497853524347068807827847, 
     -0.69117532096345247573444072455055E-01, 
      0.00000000000000000000000000000000, 
     -0.56115415639613667473122134344096E-01, 
     -0.24397413927543615902606748562471, 
     -0.27325078883205130181299635420966, 
      0.00000000000000000000000000000000, 
     -0.16627205471955132096796330738419, 
     -0.93789271474578126741505762079105, 
     -0.11930506754016075016780500019251, 
     -0.90151108565688587759737165100797, 
     -0.32049535943729789391066639033515, 
     -0.59104204439569195659289399130103, 
     -0.21751462540323724688100842714558, 
     -0.15800907798262943177922632788479, 
     -0.34192996428318888182680259909761, 
     -0.50348314662449992878796511382845, 
     -0.72296664168698516928824798397709, 
     -0.83410412765541264189179663045556, 
     -0.45220969294953177503238310311635, 
     -0.89091508745388259553341479950777, 
     -0.43793566633133321281524272450307, 
     -0.54799008713988301043149151643902, 
     -0.20885455826291604847376889739060, 
     -0.19653261626463977858040946368613, 
     -0.34026425560719830519871252095558, 
     -0.30334121141262021200705945773916, 
     -0.49535294281864462436231396682347, 
     -0.54967499236570302852186466340436E-01, 
     -0.79641981144983823490340768373671, 
     -0.77252621434083133665996549756551, 
     -0.22122265719250289382440633768970, 
      0.00000000000000000000000000000000, 
     -0.63310729924221124505794362578563, 
      0.00000000000000000000000000000000, 
     -0.38407586571187361708454487843458, 
     -0.53170841252843484837013141905012, 
     -0.17511720551530423320716712387413, 
     -0.30063413976092480084347753701585, 
     -0.10868692205689496593149118515046, 
     -0.41113970829952388068679218228879, 
      0.00000000000000000000000000000000, 
     -0.69178376289359713732790198931403, 
      0.00000000000000000000000000000000, 
     -0.97953562628985148065013964889033, 
     -0.95090315064315556393142498648752, 
     -0.11195104828353529165299273055927, 
     -0.39218251627345475097146235231760, 
     -0.48703846544346851999956533815843, 
     -0.74320583688319031979489348519293, 
      0.00000000000000000000000000000000, 
     -0.69393567486204502934677097789401, 
     -0.59416809204083916290436093693814, 
      0.00000000000000000000000000000000, 
     -0.59168695264172826158307513888901, 
     -0.42226270556224588488343911175678, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.77777971180155417861382042517561, 
      0.00000000000000000000000000000000, 
      0.00000000000000000000000000000000, 
     -0.37473331227815313999612621965198, 
      0.00000000000000000000000000000000, 
     -0.77988920685466297897759157424588E-01, 
      0.00000000000000000000000000000000, 
     -0.84965569809328811423376218283034 };
  static double ys[84] = { 
      0.26665777797796383627973623253967, 
     -0.57575252834615965036824127302641, 
     -0.53848652830351311098524847745854, 
     -0.19146540832960751941458046699634, 
     -0.53730488391129554992245069202693, 
     -0.57586628305377182252547466076703, 
     -0.56902146444580641033160663668216, 
     -0.53665917816023464212621359008944, 
      0.11490520500034347586508846660531E+01, 
     -0.41558253110498620541881452613918, 
     -0.50495637092434172779930285170762, 
     -0.49186089726616651601947223505482, 
     -0.54384156132974709175980042825510, 
     -0.31266847173129537509553966714426, 
     -0.46456703619309813450876313794865, 
     -0.53671264515430239685282775620237, 
     -0.48722771965406562202415837416678, 
     -0.57476168616596147470350970491546, 
     -0.55648988758506322732625931002792, 
     -0.47916934061220515500684077992539, 
     -0.30918435458414123508211511701385, 
     -0.37840666740342688089437456064073, 
     -0.52319798442392801521139581532973, 
     -0.40025389355433590091307180243175, 
     -0.51006243855179098239070640229809, 
     -0.49360682183497845303635246898716, 
     -0.56006513490944252779041221120332, 
     -0.43200879236359204554404676070914, 
     -0.56569680324298668363410710315288, 
     -0.27024376159027951510165751092348, 
     -0.50341926055772996189589528959240, 
     -0.44849821522859359956466291790167, 
     -0.15547647645906741636976279609852, 
     -0.44503618911963343982519208289006, 
     -0.36218696833371919155331150775928, 
     -0.56041252984416499585686185156419, 
     -0.52895122351604286441157165271526, 
     -0.53438108768101891752999282805610, 
     -0.54075896308780913586275538622291, 
     -0.45931123138039248803011481724000, 
     -0.45819330231577922836160500737553, 
     -0.30653810287769073602767881586011, 
     -0.55000331539356022007925762637213, 
     -0.53673989933539527787145346842252, 
     -0.34191292599382455707842137133247, 
     -0.49993522428841323273059831641668, 
     -0.15116834147609688134917263547298, 
     -0.55870738754201772527722045828665, 
     -0.53093517049471556661138723835977, 
     -0.23111044337275539896270189389374, 
     -0.43735474019364034368275094309953, 
     -0.56061901718698357069687788950861, 
      0.76210269896315505375472989314401, 
     -0.50180464055963906570282323273135, 
     -0.56066404582769439179312259695008, 
     -0.37564082456563228243592800723225, 
     -0.56093583647028840530479402530989, 
     -0.52891066884917480700076516219948, 
     -0.32227834482125130031233025928199, 
      0.11168066855539393695052411427936E+01, 
     -0.57411035226075729901604357841345, 
      0.88226154891869896377217079959897, 
     -0.57376987737405061910707525396861, 
     -0.57387511563190780439997792525620, 
     -0.23684295020028737262235915750548, 
     -0.39119141236066998475482225138722, 
     -0.57412552079576272088571038657122, 
     -0.49038868751056796001858236250926, 
      0.56490646899011486213357544715935, 
     -0.44034718635067890842582839203872, 
     -0.57418257398341738868103015579407, 
      0.46314994945504571565754224755338, 
     -0.40847798164087445421263923221521, 
     -0.55956030132080839871226342430982, 
     -0.23325735799867068880579964798512, 
      0.66759979014975408805793640512288, 
     -0.57398131471203676458813819810279, 
     -0.61858570316929177810688167413239E-01, 
      0.10083409893504061073884069251960E+01, 
     -0.57386414467099150410234406203086, 
     -0.57279115920674012280479926283406, 
     -0.56023791168005362040639295690520, 
      0.12766415047174134768367989248219, 
     -0.57349445987280550794758489384318 };
  static double ws[84] = { 
      0.21904046082402179331777860407172E-02, 
      0.45710956467466776333595601564413E-03, 
      0.39488814170737652484817523330772E-03, 
      0.35934111950642942202851471876703E-02, 
      0.14910742851108025176218601254730E-02, 
      0.17967284713619579793465412642966E-03, 
      0.10125558476824481605944822924240E-02, 
      0.17069127169250074598843709053414E-02, 
      0.24079695740395151054526221911175E-04, 
      0.36927130629621416994224679037866E-02, 
      0.22489391525853586565409506926398E-02, 
      0.19524312574849497675486111229175E-02, 
      0.10866799154497416799396832774298E-02, 
      0.55988347396351114436180639136753E-02, 
      0.28924669217577866067492578922788E-02, 
      0.20380484585475350431644792897213E-02, 
      0.34261098541603004157722338638113E-02, 
      0.61351457019749324652202710478653E-03, 
      0.86542020701437739404493838669357E-03, 
      0.36671947134331494740326630077914E-02, 
      0.28026432284100394242484806225548E-02, 
      0.54349341984383407453878533427381E-02, 
      0.28351957640421590232032942807150E-02, 
      0.45746103628843216234268299714345E-02, 
      0.15952659260343677779327391072587E-02, 
      0.35621046832988029971593531686097E-02, 
      0.50771409433251890919140609918512E-03, 
      0.45268556007289308320192321203838E-02, 
      0.57647586887781341017507334442447E-03, 
      0.57986412207422415390596586269670E-02, 
      0.26551803851187666770994169697277E-02, 
      0.42618143181478615021232589520584E-02, 
      0.61410787717774668691311895702589E-02, 
      0.41577749997957744966853239679784E-02, 
      0.45714943165212042708880269013658E-02, 
      0.11784033012900141932339284269968E-02, 
      0.14219896538536742581747920401589E-02, 
      0.23860085586204575770728358003760E-02, 
      0.10070034152583604416622993805894E-02, 
      0.37609206282914322618831855234874E-02, 
      0.34952452756601121036043544579388E-02, 
      0.58926392977540359306346878698652E-02, 
      0.21761510975303425549110319628113E-02, 
      0.25597503150258418017093231009733E-02, 
      0.53360452230413602215172989866309E-02, 
      0.30890849879669851554048232262039E-02, 
      0.71328535685235871693729072620986E-02, 
      0.10810483517258989817564835528970E-02, 
      0.17544886347878652397548372765366E-02, 
      0.61494904649466669709635349908104E-02, 
      0.23575565603961950228417434213703E-02, 
      0.13783891229758136398436483499410E-02, 
      0.18290788738495703825316717211575E-02, 
      0.33652181243651173378348463861955E-02, 
      0.15184413550438083995525042414313E-02, 
      0.55587775921145309411995749989960E-02, 
      0.17601404353050383171865064267571E-02, 
      0.31341293249895626259540265265832E-02, 
      0.56365508134121191078584408821979E-02, 
      0.19737963916342285028504669269696E-03, 
      0.57998586649655338308131100175788E-03, 
      0.13607360815822191768414153771193E-02, 
      0.15219771844702045466850220251726E-03, 
      0.24008493579443378017998981782239E-03, 
      0.67581925200909458799817643884140E-02, 
      0.51393423392662772922503829588641E-02, 
      0.69517340368258314203719954545839E-03, 
      0.24993080895095462395459033935751E-02, 
      0.26503828818360718009873158711194E-02, 
      0.30717948054284111140681604233590E-02, 
      0.63523033610323033970102101013160E-03, 
      0.30076725635797333204424140098894E-02, 
      0.42759500648291081429579191318110E-02, 
      0.17369109423123418439135661109544E-02, 
      0.33989953518860807286401877670786E-02, 
      0.22934322638873487811005284243824E-02, 
      0.52212507063839734314031356362762E-03, 
      0.38540006507433082413851188895204E-02, 
      0.78895266863538848229357197934055E-03, 
      0.77068430988334312857080924995683E-03, 
      0.49308164137002449756423311284903E-03, 
      0.18435234238175567747584066472141E-02, 
      0.37992680605946786583823819678412E-02, 
      0.48761272444903304947626163241666E-03 };

  r8vec_copy ( n, xs, x );
  r8vec_copy ( n, ys, y );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

double *simplex_to_triangle ( double tvert1[], double tvert2[], 
  double tvert3[], double s[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_TO_TRIANGLE maps points from the simplex to a triangle.
//
//  Discussion:
//
//    The simplex has vertices:
//
//      (  0, 0 )
//      (  1, 0 )
//      (  0, 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TVERT1[2], TVERT2[2], TVERT3[2], the coordinates
//    of the vertices of the triangle.  These vertices will be taken
//    to be the images of (0,0), (1,0) and (0,1) respectively.
//
//    Input, double S[2], the coordinates of the point in the simplex.
//
//    Output, double SIMPLEX_TO_TRIANGLE[2], the coordinates of the point in
//    the triangle.
//
{
  int i;
  double *t;

  t = new double[2];

  for ( i = 0; i < 2; i++ )
  {
    t[i] = tvert1[i] * ( 1.0 - s[0] - s[1] ) 
         + tvert2[i] * s[0]
         + tvert3[i] * s[1];
  }

  return t;
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

double triangle_area ( double vert1[], double vert2[], double vert3[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA returns the area of a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double VERT1[2], VERT2[2], VERT3[2], the vertices of
//    the triangle.
//
//    Output, double TRIANGLE_AREA, the area of the triangle.
//
{
  double value;

  value = 0.5 * 
    ( 
        ( vert2[0] - vert1[0] ) * ( vert3[1] - vert1[1] ) 
      - ( vert3[0] - vert1[0] ) * ( vert2[1] - vert1[1] ) 
    );

  return value;
}
//****************************************************************************80

double *triangle_to_ref ( double tvert1[], double tvert2[], 
  double tvert3[], double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_TO_REF maps points from any triangle to the reference triangle.
//
//  Discussion:
//
//    The reference triangle has vertices:
//
//      ( -1, -1/sqrt(3) )
//      ( +1, -1/sqrt(3) )
//      (  0, +2/sqrt(3) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TVERT1[2], TVERT2[2], TVERT3[2], the coordinates
//    of the vertices of the triangle.  These vertices will be taken
//    to be the images of (0,0), (1,0) and (0,1) respectively.
//
//    Input, double T[2], the coordinates of a point in the triangle.
//
//    Output, double TRIANGLE_TO_REF[2], the coordinates of a point in the
//    reference triangle.
//
{
  double *r;
  double rvert1[2];
  double rvert2[2];
  double rvert3[2];
  double *s;

  s = triangle_to_simplex ( tvert1, tvert2, tvert3, t );

  rvert1[0] = - 1.0;
  rvert1[1] = - 1.0 / sqrt ( 3.0 );
  rvert2[0] = + 1.0;
  rvert2[1] = - 1.0 / sqrt ( 3.0 );
  rvert3[0] =   0.0;
  rvert3[1] =   2.0 / sqrt ( 3.0 );

  r = simplex_to_triangle ( rvert1, rvert2, rvert3, s );

  delete [] s;

  return r;
}
//****************************************************************************80

double *triangle_to_simplex ( double tvert1[], double tvert2[], 
  double tvert3[], double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_TO_SIMPLEX maps points from any triangle to the simplex.
//
//  Discussion:
//
//    The simplex has vertices:
//
//      (  0, 0 )
//      (  1, 0 )
//      (  0, 1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TVERT1[2], TVERT2[2], TVERT3[2], the coordinates
//    of the vertices of the triangle.  These vertices will be taken
//    to be the images of (0,0), (1,0) and (0,1) respectively.
//
//    Input, double T[2]), the coordinates of a point in the triangle.
//
//    Output, double TRIANGLE_TO_SIMPLEX[2], the coordinates of the
//    point in the simplex.
//
{
  double *s;

  s = new double[2];

  s[0] = ( ( tvert3[1] - tvert1[1] ) * ( t[0]      - tvert1[0] ) 
         - ( tvert3[0] - tvert1[0] ) * ( t[1]      - tvert1[1] ) ) 
       / ( ( tvert3[1] - tvert1[1] ) * ( tvert2[0] - tvert1[0] ) 
         - ( tvert3[0] - tvert1[0] ) * ( tvert2[1] - tvert1[1] ) );

  s[1] = ( ( tvert2[0] - tvert1[0] ) * ( t[1]      - tvert1[1] ) 
         - ( tvert2[1] - tvert1[1] ) * ( t[0]      - tvert1[0] ) ) 
       / ( ( tvert3[1] - tvert1[1] ) * ( tvert2[0] - tvert1[0] ) 
         - ( tvert3[0] - tvert1[0] ) * ( tvert2[1] - tvert1[1] ) );

  return s;
}
//****************************************************************************80

void trianmap ( int numnodes, double vert1[], double vert2[], double vert3[], 
  double rnodes[], double whts[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANMAP maps rules from the reference triangle to the user triangle.
//
//  Discussion:
//
//    This routine maps the quadrature nodes on the reference
//    triangle into a user-defined triangle specified by its vertices.
//
//    The weights are rescaled accordingly, so that
//    the resulting quadrature will integrate correctly constants, i.e.
//    the sum of the weights is the area of the user-defined triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int NUMNODES, the number of nodes.
//
//    Input, double VERT1[2], VERT2[2], VERT3[2], the vertices of
//    the triangle on which the quadrature rule is to be constructed.
//
//    Input/output, double RNODES[2*NUMNODES], the nodes.
//
//    Input/output, double WHTS[NUMNODES], the weights.
//
{
  double area;
  int j;
  double scale;
  double u;
  double v;
  double x;
  double y;

  area = fabs ( triangle_area ( vert1, vert2, vert3 ) );

  scale = r8vec_sum ( numnodes, whts );

  scale = area / scale;

  for ( j = 0; j < numnodes; j++ )
  {
    x = rnodes[0+j*2];
    y = rnodes[1+j*2];
    triasimp ( x, y, u, v );
    x = ( vert2[0] - vert1[0] ) * u 
      + ( vert3[0] - vert1[0] ) * v + vert1[0];
    y = ( vert2[1] - vert1[1] ) * u 
      + ( vert3[1] - vert1[1] ) * v + vert1[1];
    rnodes[0+j*2] = x;
    rnodes[1+j*2] = y;
    whts[j] = whts[j] * scale;
  }

  return;
}
//****************************************************************************80

void triasimp ( double x, double y, double &uout, double &vout )

//****************************************************************************80
//
//  Purpose:
//
//    TRIASIMP maps a point from the reference triangle to the simplex.
//
//  Discussion:
//
//    Map the reference triangle with vertices
//      (-1,-1/sqrt(3)), (1,-1/sqrt(3)), (0,2/sqrt(3))
//    to the simplex with vertices
//      (0,0), (1,0), (0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double X, Y, the coordinates of a point in the
//    reference triangle.
//
//    Output, double &UOUT, &VOUT, the coordinates of the corresponding
//    point in the simplex.
//
{
  double scale;

  scale = 1.0 / sqrt ( 3.0 );

  uout = 0.5 * ( x + 1.0 ) - 0.5 * scale * ( y + scale );
  vout = 1.0 * scale * ( y + scale );

  return;
}
//****************************************************************************80

void triasymq ( int n, double vert1[], double vert2[], double vert3[], 
  double rnodes[], double weights[], int numnodes )

//****************************************************************************80
//
//  Purpose:
//
//    TRIASYMQ returns a symmetric quadrature formula for a user triangle.
//
//  Discussion:
//
//    This routine constructs (or rather retrieves) a symmetric
//    quadrature formula on the user-defined triangle in the plane
//    specified by its vertices (vert1,vert2,vert3).
//
//    The total number of nodes and weights to be created is numnodes.
//
//      n       1     2     3     4     5     6     7     8     9    10
//     -----------------------------------------------------------------
//    nodes     1     3     6     6     7    12    15    16    19    25
//
//
//      n      11    12    13    14    15    16    17    18    19    20
//     -----------------------------------------------------------------
//    nodes    28    33    37    42    49    55    60    67    73    79
//
//
//      n      21    22    23    24    25    26    27    28    29    30
//     -----------------------------------------------------------------
//    nodes    87    96   103   112   120   130   141   150   159   171
//
//
//      n      31    32    33    34    35    36    37    38    39    40
//     -----------------------------------------------------------------
//    nodes   181   193   204   214   228   243   252   267   282   295
//
//
//      n      41    42    43    44    45    46    47    48    49    50
//     -----------------------------------------------------------------
//    nodes   309   324   339   354   370   385   399   423   435   453
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    30 June 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the degree of the quadrature (must not exceed 50).
//    Note that the total number of nodes to be created is numnodes.
//
//    Input, double VERT1[2], VERT2[2], VERT3[2], the vertices of
//    the triangle on which the quadrature rule is to be constructed.
//
//    Output, double RNODES[2*NUMNODES], the nodes in the plane
//    (all inside the user-supplied triangle)
//
//    Output, double WEIGHTS[NUMNODES], the quadrature weights.
//
//    Input, int NUMNODES, the number of nodes.
//
{
  int itype;

  itype = 0;

  quaequad ( itype, n, rnodes, weights, numnodes );

  trianmap ( numnodes, vert1, vert2, vert3, rnodes, weights );

  return;
}
//****************************************************************************80

void triasymq_gnuplot ( double vert1[], double vert2[], double vert3[], 
  int numnodes, double rnodes[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    TRIASYMQ_GNUPLOT: set up a GNUPLOT plot of the triangle quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double VERT1[2], VERT2[2], VERT3[2], the vertices
//    of the triangle.
//
//    Input, double NUMNODES, the number of nodes.
//
//    Input, double RNODES[2*NUMNODES], the coordinates of
//    the nodes.
//
//    Input, string HEADER, a string to be used to identify
//    the files created.
//
{
  string command_filename;
  ofstream command_unit;
  int i;
  int j;
  int l;
  int n;
  string node_filename;
  ofstream node_unit;
  string plot_filename;
  string vertex_filename;
  ofstream vertex_unit;
//
//  Create the vertex file.
//
  vertex_filename = header + "_vertices.txt";
  vertex_unit.open ( vertex_filename.c_str ( ) );
  vertex_unit << vert1[0] << "  "
                << vert1[1] << "\n";
  vertex_unit << vert2[0] << "  "
                << vert2[1] << "\n";
  vertex_unit << vert3[0] << "  "
                << vert3[1] << "\n";
  vertex_unit << vert1[0] << "  "
                << vert1[1] << "\n";
  vertex_unit.close ( );
  cout << "\n";
  cout << "  Created vertex file '" << vertex_filename << "'\n";
//
//  Create node file.
//
  node_filename = header + "_nodes.txt";
  node_unit.open ( node_filename.c_str ( ) );
  for ( j = 0; j < numnodes; j++ )
  {
    node_unit << rnodes[0+j*2] << "  "
              << rnodes[1+j*2] << "\n";
  }
  node_unit.close ( );
  cout << "  Created node file '" << node_filename << "'\n";
//
//  Create graphics command file.
//
  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  plot_filename = header + ".png";
  command_unit << "set output '" << plot_filename << "'\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set title '" << header << "'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "set size ratio -1\n";
  command_unit << "set style data lines\n";
  command_unit << "set timestamp\n";
  command_unit << "plot '" << vertex_filename << "' with lines lw 3, \\\n";
  command_unit << "     '" << node_filename << "' with points pt 7 lt 0\n";
  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  return;
}
