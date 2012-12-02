# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "toms446.hpp"

//****************************************************************************80

double *binom ( double x[], double xx[], int npl, int m, int nt )

//****************************************************************************80
//
//  Purpose:
//
//    BINOM: binomial expansion series for the (-1/M) power of a Chebyshev series.
//
//  Discussion:
//
//    This routine uses a certain number of terms of the binomial expansion 
//    series to estimate the (-1/M) power of a given Chebyshev series. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double X[NPL], the given Chebyshev series.
//
//    Input, double XX[NPL], an initial estimate for
//    the Chebyshev series for the input function raised to the (-1/M) power.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Input, int M, defines the exponent, (-1/M).
//    0 < M.
//
//    Input, int NT, the number of terms of the binomial
//    series to be used.
//
//    Output, double BINOM[NPL], the estimated Chebyshev series
//    for the input function raised to the (-1/M) power.
//
{
  double alfa;
  double coef;
  double dkm2;
  double dkmm;
  double dm;
  int k;
  double *w2;
  double *w3;
  double *ww;
  double *xa;

  dm = ( double ) ( m );
  alfa = -1.0 / dm;

  ww = r8vec_copy_new ( npl, x );

  w2 = new double[npl];

  for ( k = 1; k <= m; k++ )
  {
    mltply ( ww, xx, npl, w2 );
    r8vec_copy ( npl, w2, ww );
  }

  ww[0] = ww[0] - 2.0;

  xa = r8vec_zero_new ( npl );
  xa[0] = 2.0;

  w3 = r8vec_copy_new ( npl, xa );

  for ( k = 2; k <= nt; k++ )
  {
    dkmm = ( double ) ( k - 1 );
    dkm2 = ( double ) ( k - 2 );
    coef = ( alfa - dkm2 ) / dkmm;

    mltply ( w3, ww, npl, w2 );

    r8vec_copy ( npl, w2, w3 );
    r8vec_scale ( coef, npl, w3 );
    r8vec_add ( npl, w3, xa );
  }
  mltply ( xa, xx, npl, w2 );

  r8vec_copy ( npl, w2, xa );

  delete [] ww;
  delete [] w2;
  delete [] w3;

  return xa;
}
//****************************************************************************80

double *cheby ( int nf, int npl, double *functn ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY carries out the Chebyshev analysis of one or more functions.
//
//  Discussion:
//
//    This routine carries out the simultaneous Chebyshev analysis of 
//    NF functions.
//
//    The output is a matrix containing one Chebyshev series per column.
//
//    An example of a routine to compute the function values is:
//
//      double *functn ( double a )
//      {
//        double *val;
//        val = new double[2];
//        val[0] = sin(a);
//        val[1] = cos(a);
//        return val;
//      }
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, int NF, the number of functions to be analyzed.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Input, int NPLMAX, the leading dimension of X.
//
//    Input, external FUNCTN, the name of a routine which computes
//    the function values at any given point.
//
//    Output, double CHEBY[NPL*NF], the Chebyshev series.
//
{
  double enn;
  double enw;
  double fac;
  double fk;
  double *fxj;
  double *gc;
  int j;
  int k;
  int l;
  int lm;
  int n;
  int n2;
  double pen;
  double *x;
  double xj;

  x = r8vec_zero_new ( npl * nf );

  n = npl - 1;
  enn = ( double ) ( n );
  pen = 3.1415926535897932 / enn;

  gc = new double[2*n];
  for ( k = 1; k <= 2 * n; k++ )
  {
    fk = ( double ) ( k - 1 );
    gc[k-1] = cos ( fk * pen );
  }

  for ( j = 0; j < npl; j++ )
  {
    xj = gc[j];
    fxj = functn ( xj );

    if ( j == 0 || j == npl - 1 )
    {
      r8vec_scale ( 0.5, nf, fxj );
    }
    for ( l = 0; l < npl; l++ )
    {
      lm = ( l * j ) % ( 2 * n );
      for ( k = 0; k < nf; k++ )
      {
        x[l+k*npl] = x[l+k*npl] + fxj[k] * gc[lm];
      }
    }
  }

  r8vec_scale ( 2.0 / enn, npl * nf, x );

  delete [] fxj;
  delete [] gc;

  return x;
}
//****************************************************************************80

double *dfrnt ( double xx[], int npl )

//****************************************************************************80
//
//  Purpose:
//
//    DFRNT determines the derivative of a Chebyshev series.
//
//  Discussion:
//
//    This routine computes the Chebyshev series of the derivative of a 
//    function whose Chebyshev series is given.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double XX[NPL], the given Chebyshev series.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Output, double DFRNT[NPL], the Chebyshev series for the
//    derivative.
//
{
  double dl;
  double dn;
  int k;
  int l;
  double *x2;
  double xn;
  double xxl;
  double xxn;

  x2 = new double[npl];
  dn = ( double ) ( npl - 1 );
  xxn = xx[npl-2];
  x2[npl-2] = 2.0 * xx[npl-1] * dn;
  x2[npl-1] = 0.0;

  for ( k = 3; k <= npl; k++ )
  {
    l = npl - k + 1;
    dl = ( double ) ( l );
    xxl = xx[l-1];
    x2[l-1] = x2[l+1] + 2.0 * xxn * dl;
    xxn = xxl;
  }
  return x2;
}
//****************************************************************************80

double echeb ( double x, double coef[], int npl )

//****************************************************************************80
//
//  Purpose:
//
//    ECHEB evaluates a Chebyshev series at a point.
//
//  Discussion:
//
//    This routine evaluates a Chebyshev series at a point in [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NPL], the Chebyshev series.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Output, double ECHEB, the value of the Chebyshev series at X.
//
{
  double br;
  double brp2;
  double brpp;
  double fx;
  int j;
  int k;

  br = 0.0;
  brpp = 0.0;

  for ( k = 1; k <= npl; k++ )
  {
    j = npl - k + 1;
    brp2 = brpp;
    brpp = br;
    br = 2.0 * x * brpp - brp2 + coef[j-1];
  }
  fx = 0.5 * ( br - brp2 );
  return fx;
}
//****************************************************************************80

double edcheb ( double x, double coef[], int npl )

//****************************************************************************80
//
//  Purpose:
//
//    EDCHEB evaluates the derivative of a Chebyshev series at a point.
//
//  Discussion:
//
//    This routine evaluates the derivative of a Chebyshev series 
//    at a point in [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NPL], the Chebyshev series.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Output, double EDCHEB, the value of the derivative of the
//    Chebyshev series at X.
//
{
  double bf;
  double bj;
  double bjp2;
  double bjpl;
  double dj;
  double fx;
  int j;
  int k;
  int n;
  double xj;
  double xjp2;
  double xjpl;

  xjp2 = 0.0;
  xjpl = 0.0;
  bjp2 = 0.0;
  bjpl = 0.0;
  n = npl - 1;

  for ( k = 1; k <= n; k++ )
  {
    j = npl - k;
    dj = ( double ) ( j );
    xj = 2.0 * coef[j] * dj + xjp2;
    bj = 2.0 * x * bjpl - bjp2 + xj;
    bf = bjp2;
    bjp2 = bjpl;
    bjpl = bj;
    xjp2 = xjpl;
    xjpl = xj;
  }
  fx = 0.5 * ( bj - bf );
  return fx;
}
//****************************************************************************80

double *invert ( double x[], double xx[], int npl, int net )

//****************************************************************************80
//
//  Purpose:
//
//    INVERT computes the inverse Chebyshev series.
//
//  Discussion:
//
//    This routine computes the inverse of a Chebyshev series, starting with
//    an initial approximation XX. 
//
//    The routine uses the Euler method and computes all powers EPS^K 
//    up to K=2^(NET+1), where EPS = 1 - X * ( XX inverse ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double X[NPL], the Chebyshev series.
//
//    Input, double XX[NPL], an initial approximation for the
//    inverse Chebyshev series.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Input, int NET, the number of iterations to take.
//
//    Output, double INVERT[NPL], the estimated Chebyshev
//    series of the inverse function.
//
{
  int k;
  double s;
  double *w2;
  double *ww;
  double *xnvse;

  ww = mltply_new ( x, xx, npl );

  s = - 1.0;
  r8vec_scale ( s, npl, ww );
  ww[0] = ww[0] + 2.0;

  w2 = mltply_new ( ww, ww, npl );
  ww[0] = 2.0 * ww[0];

  xnvse = new double[npl];

  for ( k = 1; k <= net; k++ )
  {
    mltply ( ww, w2, npl, xnvse );

    r8vec_add ( npl, xnvse, ww );

    mltply ( w2, w2, npl, xnvse );

    r8vec_copy ( npl, xnvse, w2 );
  }
  mltply ( ww, xx, npl, xnvse );

  delete [] w2;
  delete [] ww;

  return xnvse;
}
//****************************************************************************80

void mltply ( double xx[], double x2[], int npl, double x3[] )

//****************************************************************************80
//
//  Purpose:
//
//    MLTPLY_NEW multiplies two Chebyshev series.
//
//  Discussion:
//
//    This routine multiplies two given Chebyshev series, XX and X2,
//    to produce an output Chebyshev series, X3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double XX[NPL], the first Chebyshev series.
//
//    Input, double X2[NPL], the second Chebyshev series.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Output, double X3[NPL], the Chebyshev series of the
//    product.
//
{
  double ex;
  int k;
  int l;
  int m;
  int mm;

  for ( k = 1; k <= npl; k++ )
  {
    ex = 0.0;
    mm = npl - k + 1;
    for ( m = 1; m <= mm; m++ )
    {
      l = m + k - 1;
      ex = ex + xx[m-1] * x2[l-1] + xx[l-1] * x2[m-1];
    }
    x3[k-1] = 0.5 * ex;
  }

  x3[0] = x3[0] - 0.5 * xx[0] * x2[0];

  for ( k = 3; k <= npl; k++ )
  {
    ex = 0.0;
    mm = k - 1;
    for ( m = 2; m <= mm; m++ )
    {
      l = k - m + 1;
      ex = ex + xx[m-1] * x2[l-1];
    }
    x3[k-1] = 0.5 * ex + x3[k-1];
  }
  return;
}
//****************************************************************************80

double *mltply_new ( double xx[], double x2[], int npl )

//****************************************************************************80
//
//  Purpose:
//
//    MLTPLY_NEW multiplies two Chebyshev series.
//
//  Discussion:
//
//    This routine multiplies two given Chebyshev series, XX and X2,
//    to produce an output Chebyshev series, X3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double XX[NPL], the first Chebyshev series.
//
//    Input, double X2[NPL], the second Chebyshev series.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Output, double MLTPLY_NEW[NPL], the Chebyshev series of the
//    product.
//
{
  double *x3;

  x3 = r8vec_zero_new ( npl );

  mltply ( xx, x2, npl, x3 );

  return x3;
}
//****************************************************************************80

double *ntgrt ( double xx[], int npl )

//****************************************************************************80
//
//  Purpose:
//
//    NTGRT determines the integral of a Chebyshev series.
//
//  Discussion:
//
//    This routine computes the Chebyshev series for the integral of a 
//    function whose Chebyshev series is given.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double XX[NPL], the Chebyshev series.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Output, double NTGRT[NPL], the Chebyshev series for the
//    integral of the function.
//
{
  double dk;
  int k;
  int n;
  double term;
  double *x2;
  double xpr;

  x2 = new double[npl];

  xpr = xx[0];
  x2[0] = 0.0;
  n = npl - 1;

  for ( k = 2; k <= n; k++ )
  {
    dk = ( double ) ( k - 1 );
    term = ( xpr - xx[k] ) / ( 2.0 * dk );
    xpr = xx[k-1];
    x2[k-1] = term;
  }

  dk = ( double ) n;
  x2[npl-1] = xpr / ( 2.0 * dk );

  return x2;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

void r8vec_add ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ADD adds one R8VEC to another.
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
//    22 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be added.
//
//    Input/output, double A2[N], the vector to be increased.
//    On output, A2 = A2 + A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a2[i] + a1[i];
  }
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

double *r8vec_copy_new ( int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY_NEW copies an R8VEC to a "new" R8VEC.
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
//    03 July 2008
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
//    Output, double R8VEC_COPY_NEW[N], the copy of A1.
//
{
  double *a2;
  int i;

  a2 = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

void r8vec_scale ( double s, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SCALE multiples an R8VEC by a scale factor.
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
//    22 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double S, the scale factor.
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double A[N], the vector to be scaled.
//    On output, A[] = S * A[].
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = s * a[i];
  }
  return;
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

void xalfa2 ( double x[], double xx[], int npl, int m, int maxet, double epsln, 
  int &net )

//****************************************************************************80
//
//  Purpose:
//
//    XALFA2 computes a Chebyshev series raised to the (-1/M) power.
//
//  Discussion:
//
//    This routine estimates the Chebyshev series for a function raised
//    to the (-1/M) power, given the Chebyshev series for the function,
//    and a starting estimate for the desired series.
//
//    The convergence is quadratic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double X[NPL], the Chebyshev series of the function.
//
//    Input/output, double XX[NPL].  On input, an initial
//    approximation to the Chebyshev series of the function raised to the
//    (-1/M) power.  On output, an improved approximation.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Input, int M, determines the exponent (-1/M).
//
//    Input, int MAXET, the maximum number of iterations.
//
//    Input, double EPSLN, the required precision.
//
//    Output, int &NET, the actual number of iterations.
//
{
  double dalfa;;
  double dm;
  int i;
  int jx;
  int k;
  double s;
  double t;
  double tdmm;
  double *w2;
  double *ww;

  dm = ( double ) ( m );
  dalfa = 1.0 / dm;
  tdmm = 2.0 * ( dm + 1.0 );
 
  ww = new double[npl];
  w2 = new double[npl];

  for ( jx = 1; jx <= maxet; jx++ )
  {
    r8vec_copy ( npl, x, ww );
 
    for ( k = 1; k <= m; k++ )
    {
      mltply ( ww, xx, npl, w2 );
      r8vec_copy ( npl, w2, ww );
    }

    t = - 2.0;
    for ( i = 0; i < npl; i++ )
    {
      t = t + r8_abs ( ww[i] );
    }

    s = - 1.0;
    r8vec_scale ( s, npl, ww );
    ww[0] = ww[0] + tdmm;

    mltply ( ww, xx, npl, w2 );

    r8vec_copy ( npl, w2, xx );
    r8vec_scale ( dalfa, npl, xx );

    net = jx;

    if ( r8_abs ( t ) < epsln )
    {
      break;
    }
  }

  delete [] ww;
  delete [] w2;

  return;
}
//****************************************************************************80

void xalfa3 ( double x[], double xx[], int npl, int m, int maxet, double epsln, 
  int &net )

//****************************************************************************80
//
//  Purpose:
//
//    XALFA3 computes a Chebyshev series raised to the (-1/M) power.
//
//  Discussion:
//
//    This routine estimates the Chebyshev series for a function raised
//    to the (-1/M) power, given the Chebyshev series for the function,
//    and a starting estimate for the desired series.
//
//    The convergence is of order three.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Roger Broucke.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    October 1973, Volume 16, Number 4, pages 254-256.
//
//  Parameters:
//
//    Input, double X[NPL], the Chebyshev series of the function.
//
//    Input/output, double XX[NPL].  On input, an initial
//    approximation to the Chebyshev series of the function raised to the
//    (-1/M) power.  On output, an improved approximation.
//
//    Input, int NPL, the number of terms in the 
//    Chebyshev series.
//
//    Input, int M, determines the exponent (-1/M).
//
//    Input, int MAXET, the maximum number of iterations.
//
//    Input, double EPSLN, the required precision.
//
//    Output, int &NET, the actual number of iterations.
//
{
  double dalfa;
  double dm;
  int i;
  int jx;
  int k;
  double p5dml;
  double s;
  double t;
  double tdmm;
  double *w2;
  double *ww;

  dm = ( double ) ( m );
  dalfa = 1.0 / dm;
  tdmm = 2.0 * ( dm + 1.0 );
  p5dml = 0.5 * ( dm + 1.0 );

  ww = new double[npl];
  w2 = new double[npl];

  for ( jx = 1; jx <= maxet; jx++ )
  {
    r8vec_copy ( npl, x, ww );

    for ( k = 1; k <= m; k++ )
    {
      mltply ( ww, xx, npl, w2 );
      r8vec_copy ( npl, w2, ww );
    }

    t = - 2.0;
    for ( i = 0; i < npl; i++ )
    {
      t = t + r8_abs ( ww[i] );
    }

    ww[0] = ww[0] - 2.0;
    r8vec_scale ( dalfa, npl, ww );

    mltply ( ww, ww, npl, w2 );

    s = - 1.0;
    r8vec_scale ( s, npl, ww );

    r8vec_scale ( p5dml, npl, w2 );

    ww[0] = ww[0] + 2.0;

    r8vec_add ( npl, ww, w2 );

    mltply ( w2, xx, npl, ww );

    r8vec_copy ( npl, ww, xx );

    net = jx;

    if ( r8_abs ( t ) < epsln )
    {
      break;
    }
  }

  delete [] ww;
  delete [] w2;

  return;
}
