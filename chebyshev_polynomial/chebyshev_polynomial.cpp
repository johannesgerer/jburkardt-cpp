# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

# include "chebyshev_polynomial.hpp"

//****************************************************************************80

void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DA, the multiplier of DX.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Input/output, double DY[*], the second vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
//    Input, int INCY, the increment between successive entries of DY.
//
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( da == 0.0 )
  {
    return;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dy[iy] = dy[iy] + da * dx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 4;

    for ( i = 0; i < m; i++ )
    {
      dy[i] = dy[i] + da * dx[i];
    }

    for ( i = m; i < n; i = i + 4 )
    {
      dy[i  ] = dy[i  ] + da * dx[i  ];
      dy[i+1] = dy[i+1] + da * dx[i+1];
      dy[i+2] = dy[i+2] + da * dx[i+2];
      dy[i+3] = dy[i+3] + da * dx[i+3];
    }

  }

  return;
}
//****************************************************************************80

double ddot ( int n, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DDOT forms the dot product of two vectors.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries in DX.
//
//    Input, double DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries in DY.
//
//    Output, double DDOT, the sum of the product of the corresponding
//    entries of DX and DY.
//
{
  double dtemp;
  int i;
  int ix;
  int iy;
  int m;

  dtemp = 0.0;

  if ( n <= 0 )
  {
    return dtemp;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dtemp = dtemp + dx[ix] * dy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      dtemp = dtemp + dx[i] * dy[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      dtemp = dtemp + dx[i  ] * dy[i  ]
                    + dx[i+1] * dy[i+1]
                    + dx[i+2] * dy[i+2]
                    + dx[i+3] * dy[i+3]
                    + dx[i+4] * dy[i+4];
    }

  }

  return dtemp;
}
//****************************************************************************80

double dnrm2 ( int n, double x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DNRM2 returns the euclidean norm of a vector.
//
//  Discussion:
//
//     DNRM2 ( X ) = sqrt ( X' * X )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[*], the vector whose norm is to be computed.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, double DNRM2, the Euclidean norm of X.
//
{
  double absxi;
  int i;
  int ix;
  double norm;
  double scale;
  double ssq;
  double value;

  if ( n < 1 || incx < 1 )
  {
    norm = 0.0;
  }
  else if ( n == 1 )
  {
    norm = r8_abs ( x[0] );
  }
  else
  {
    scale = 0.0;
    ssq = 1.0;
    ix = 0;

    for ( i = 0; i < n; i++ )
    {
      if ( x[ix] != 0.0 )
      {
        absxi = r8_abs ( x[ix] );
        if ( scale < absxi )
        {
          ssq = 1.0 + ssq * ( scale / absxi ) * ( scale / absxi );
          scale = absxi;
        }
        else
        {
          ssq = ssq + ( absxi / scale ) * ( absxi / scale );
        }
      }
      ix = ix + incx;
    }

    norm  = scale * sqrt ( ssq );
  }

  return norm;
}
//****************************************************************************80

void drot ( int n, double x[], int incx, double y[], int incy, double c,
  double s )

//****************************************************************************80
//
//  Purpose:
//
//    DROT applies a plane rotation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to be rotated.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to be rotated.
//
//    Input, int INCY, the increment between successive elements of Y.
//
//    Input, double C, S, parameters (presumably the cosine and
//    sine of some angle) that define a plane rotation.
//
{
  int i;
  int ix;
  int iy;
  double stemp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      stemp = c * x[i] + s * y[i];
      y[i]  = c * y[i] - s * x[i];
      x[i]  = stemp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      stemp = c * x[ix] + s * y[iy];
      y[iy] = c * y[iy] - s * x[ix];
      x[ix] = stemp;
      ix = ix + incx;
      iy = iy + incy;
    }

  }

  return;
}
//****************************************************************************80

void drotg ( double *sa, double *sb, double *c, double *s )

//****************************************************************************80
//
//  Purpose:
//
//    DROTG constructs a Givens plane rotation.
//
//  Discussion:
//
//    Given values A and B, this routine computes
//
//    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
//          = sign ( B ) if abs ( A ) <= abs ( B );
//
//    R     = SIGMA * ( A * A + B * B );
//
//    C = A / R if R is not 0
//      = 1     if R is 0;
//
//    S = B / R if R is not 0,
//        0     if R is 0.
//
//    The computed numbers then satisfy the equation
//
//    (  C  S ) ( A ) = ( R )
//    ( -S  C ) ( B ) = ( 0 )
//
//    The routine also computes
//
//    Z = S     if abs ( A ) > abs ( B ),
//      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
//      = 1     if C is 0.
//
//    The single value Z encodes C and S, and hence the rotation:
//
//    If Z = 1, set C = 0 and S = 1;
//    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
//    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input/output, double *SA, *SB,  On input, SA and SB are the values
//    A and B.  On output, SA is overwritten with R, and SB is
//    overwritten with Z.
//
//    Output, double *C, *S, the cosine and sine of the Givens rotation.
//
{
  double r;
  double roe;
  double scale;
  double z;

  if ( r8_abs ( *sb ) < r8_abs ( *sa ) )
  {
    roe = *sa;
  }
  else
  {
    roe = *sb;
  }

  scale = r8_abs ( *sa ) + r8_abs ( *sb );

  if ( scale == 0.0 )
  {
    *c = 1.0;
    *s = 0.0;
    r = 0.0;
  }
  else
  {
    r = scale * sqrt ( ( *sa / scale ) * ( *sa / scale )
                     + ( *sb / scale ) * ( *sb / scale ) );
    r = r8_sign ( roe ) * r;
    *c = *sa / r;
    *s = *sb / r;
  }

  if ( 0.0 < r8_abs ( *c ) && r8_abs ( *c ) <= *s )
  {
    z = 1.0 / *c;
  }
  else
  {
    z = *s;
  }

  *sa = r;
  *sb = z;

  return;
}
//****************************************************************************80

void dscal ( int n, double sa, double x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DSCAL scales a vector by a constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double SA, the multiplier.
//
//    Input/output, double X[*], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of X.
//
{
  int i;
  int ix;
  int m;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 )
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for ( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }

  }

  return;
}
//****************************************************************************80

int dsvdc ( double a[], int lda, int m, int n, double s[], double e[], 
  double u[], int ldu, double v[], int ldv, double work[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DSVDC computes the singular value decomposition of a real rectangular matrix.
//
//  Discussion:
//
//    This routine reduces an M by N matrix A to diagonal form by orthogonal
//    transformations U and V.  The diagonal elements S(I) are the singular
//    values of A.  The columns of U are the corresponding left singular
//    vectors, and the columns of V the right singular vectors.
//
//    The form of the singular value decomposition is then
//
//      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the M by N matrix whose
//    singular value decomposition is to be computed.  On output, the matrix
//    has been destroyed.  Depending on the user's requests, the matrix may 
//    contain other useful information.
//
//    Input, int LDA, the leading dimension of the array A.
//    LDA must be at least M.
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix A.
//
//    Output, double S[MM], where MM = min(M+1,N).  The first
//    min(M,N) entries of S contain the singular values of A arranged in
//    descending order of magnitude.
//
//    Output, double E[MM], where MM = min(M+1,N), ordinarily contains zeros.
//    However see the discussion of INFO for exceptions.
//
//    Output, double U[LDU*K].  If JOBA = 1 then K = M;
//    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of left singular
//    vectors.  U is not referenced if JOBA = 0.  If M <= N or if JOBA = 2, then
//    U may be identified with A in the subroutine call.
//
//    Input, int LDU, the leading dimension of the array U.
//    LDU must be at least M.
//
//    Output, double V[LDV*N], the N by N matrix of right singular vectors.
//    V is not referenced if JOB is 0.  If N <= M, then V may be identified
//    with A in the subroutine call.
//
//    Input, int LDV, the leading dimension of the array V.
//    LDV must be at least N.
//
//    Workspace, double WORK[M].
//
//    Input, int JOB, controls the computation of the singular
//    vectors.  It has the decimal expansion AB with the following meaning:
//      A =  0, do not compute the left singular vectors.
//      A =  1, return the M left singular vectors in U.
//      A >= 2, return the first min(M,N) singular vectors in U.
//      B =  0, do not compute the right singular vectors.
//      B =  1, return the right singular vectors in V.
//
//    Output, int *DSVDC, status indicator INFO.
//    The singular values (and their corresponding singular vectors)
//    S(*INFO+1), S(*INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
//    Thus if *INFO is 0, all the singular values and their vectors are
//    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
//    matrix with the elements of S on its diagonal and the elements of E on
//    its superdiagonal.  Thus the singular values of A and B are the same.
//
{
  double b;
  double c;
  double cs;
  double el;
  double emm1;
  double f;
  double g;
  int i;
  int info;
  int iter;
  int j;
  int jobu;
  int k;
  int kase;
  int kk;
  int l;
  int ll;
  int lls;
  int ls;
  int lu;
  int maxit = 30;
  int mm;
  int mm1;
  int mn;
  int mp1;
  int nct;
  int nctp1;
  int ncu;
  int nrt;
  int nrtp1;
  double scale;
  double shift;
  double sl;
  double sm;
  double smm1;
  double sn;
  double t;
  double t1;
  double test;
  bool wantu;
  bool wantv;
  double ztest;
//
//  Determine what is to be computed.
//
  info = 0;
  wantu = false;
  wantv = false;
  jobu = ( job % 100 ) / 10;

  if ( 1 < jobu )
  {
    ncu = i4_min ( m, n );
  }
  else
  {
    ncu = m;
  }

  if ( jobu != 0 )
  {
    wantu = true;
  }

  if ( ( job % 10 ) != 0 )
  {
    wantv = true;
  }
//
//  Reduce A to bidiagonal form, storing the diagonal elements
//  in S and the super-diagonal elements in E.
//
  nct = i4_min ( m-1, n );
  nrt = i4_max ( 0, i4_min ( m, n-2 ) );
  lu = i4_max ( nct, nrt );

  for ( l = 1; l <= lu; l++ )
  {
//
//  Compute the transformation for the L-th column and
//  place the L-th diagonal in S(L).
//
    if ( l <= nct )
    {
      s[l-1] = dnrm2 ( m-l+1, a+l-1+(l-1)*lda, 1 );

      if ( s[l-1] != 0.0 )
      {
        if ( a[l-1+(l-1)*lda] != 0.0 )
        {
          s[l-1] = r8_sign ( a[l-1+(l-1)*lda] ) * r8_abs ( s[l-1] );
        }
        dscal ( m-l+1, 1.0 / s[l-1], a+l-1+(l-1)*lda, 1 );
        a[l-1+(l-1)*lda] = 1.0 + a[l-1+(l-1)*lda];
      }
      s[l-1] = -s[l-1];
    }

    for ( j = l+1; j <= n; j++ )
    {
//
//  Apply the transformation.
//
      if ( l <= nct && s[l-1] != 0.0 )
      {
        t = - ddot ( m-l+1, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 ) 
          / a[l-1+(l-1)*lda];
        daxpy ( m-l+1, t, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 );
      }
//
//  Place the L-th row of A into E for the
//  subsequent calculation of the row transformation.
//
      e[j-1] = a[l-1+(j-1)*lda];
    }
//
//  Place the transformation in U for subsequent back multiplication.
//
    if ( wantu && l <= nct )
    {
      for ( i = l; i <= m; i++ )
      {
        u[i-1+(l-1)*ldu] = a[i-1+(l-1)*lda];
      }
    }

    if ( l <= nrt )
    {
//
//  Compute the L-th row transformation and place the
//  L-th superdiagonal in E(L).
//
      e[l-1] = dnrm2 ( n-l, e+l, 1 );

      if ( e[l-1] != 0.0 )
      {
        if ( e[l] != 0.0 )
        {
          e[l-1] = r8_sign ( e[l] ) * r8_abs ( e[l-1] );
        }
        dscal ( n-l, 1.0 / e[l-1], e+l, 1 );
        e[l] = 1.0 + e[l];
      }

      e[l-1] = -e[l-1];
//
//  Apply the transformation.
//
      if ( l+1 <= m && e[l-1] != 0.0 )
      {
        for ( j = l+1; j <= m; j++ )
        {
          work[j-1] = 0.0;
        }

        for ( j = l+1; j <= n; j++ )
        {
          daxpy ( m-l, e[j-1], a+l+(j-1)*lda, 1, work+l, 1 );
        }

        for ( j = l+1; j <= n; j++ )
        {
          daxpy ( m-l, -e[j-1]/e[l], work+l, 1, a+l+(j-1)*lda, 1 );
        }
      }
//
//  Place the transformation in V for subsequent back multiplication.
//
      if ( wantv )
      {
        for ( j = l+1; j <= n; j++ )
        {
          v[j-1+(l-1)*ldv] = e[j-1];
        }
      }
    }
  }
//
//  Set up the final bidiagonal matrix of order MN.
//
  mn = i4_min ( m + 1, n );
  nctp1 = nct + 1;
  nrtp1 = nrt + 1;

  if ( nct < n )
  {
    s[nctp1-1] = a[nctp1-1+(nctp1-1)*lda];
  }

  if ( m < mn )
  {
    s[mn-1] = 0.0;
  }

  if ( nrtp1 < mn )
  {
    e[nrtp1-1] = a[nrtp1-1+(mn-1)*lda];
  }

  e[mn-1] = 0.0;
//
//  If required, generate U.
//
  if ( wantu )
  {
    for ( i = 1; i <= m; i++ )
    {
      for ( j = nctp1; j <= ncu; j++ )
      {
        u[(i-1)+(j-1)*ldu] = 0.0;
      }
    }

    for ( j = nctp1; j <= ncu; j++ )
    {
      u[j-1+(j-1)*ldu] = 1.0;
    }

    for ( ll = 1; ll <= nct; ll++ )
    {
      l = nct - ll + 1;

      if ( s[l-1] != 0.0 )
      {
        for ( j = l+1; j <= ncu; j++ )
        {
          t = - ddot ( m-l+1, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 ) 
            / u[l-1+(l-1)*ldu];
          daxpy ( m-l+1, t, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 );
        }

        dscal ( m-l+1, -1.0, u+(l-1)+(l-1)*ldu, 1 );
        u[l-1+(l-1)*ldu] = 1.0 + u[l-1+(l-1)*ldu];
        for ( i = 1; i <= l-1; i++ )
        {
          u[i-1+(l-1)*ldu] = 0.0;
        }
      }
      else
      {
        for ( i = 1; i <= m; i++ )
        {
          u[i-1+(l-1)*ldu] = 0.0;
        }
        u[l-1+(l-1)*ldu] = 1.0;
      }
    }
  }
//
//  If it is required, generate V.
//
  if ( wantv )
  {
    for ( ll = 1; ll <= n; ll++ )
    {
      l = n - ll + 1;

      if ( l <= nrt && e[l-1] != 0.0 )
      {
        for ( j = l+1; j <= n; j++ )
        {
          t = - ddot ( n-l, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 ) 
            / v[l+(l-1)*ldv];
          daxpy ( n-l, t, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 );
        }

      }
      for ( i = 1; i <= n; i++ )
      {
        v[i-1+(l-1)*ldv] = 0.0;
      }
      v[l-1+(l-1)*ldv] = 1.0;
    }
  }
//
//  Main iteration loop for the singular values.
//
  mm = mn;
  iter = 0;

  while ( 0 < mn )
  {
//
//  If too many iterations have been performed, set flag and return.
//
    if ( maxit <= iter )
    {
      info = mn;
      return info;
    }
//
//  This section of the program inspects for
//  negligible elements in the S and E arrays.
//
//  On completion the variables KASE and L are set as follows:
//
//  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
//  KASE = 2     if S(L) is negligible and L < MN
//  KASE = 3     if E(L-1) is negligible, L < MN, and
//               S(L), ..., S(MN) are not negligible (QR step).
//  KASE = 4     if E(MN-1) is negligible (convergence).
//
    for ( ll = 1; ll <= mn; ll++ )
    {
      l = mn - ll;

      if ( l == 0 )
      {
        break;
      }

      test = r8_abs ( s[l-1] ) + r8_abs ( s[l] );
      ztest = test + r8_abs ( e[l-1] );

      if ( ztest == test )
      {
        e[l-1] = 0.0;
        break;
      }
    }

    if ( l == mn - 1 )
    {
      kase = 4;
    }
    else
    {
      mp1 = mn + 1;

      for ( lls = l+1; lls <= mn+1; lls++ )
      {
        ls = mn - lls + l + 1;

        if ( ls == l )
        {
          break;
        }

        test = 0.0;
        if ( ls != mn )
        {
          test = test + r8_abs ( e[ls-1] );
        }

        if ( ls != l + 1 )
        {
          test = test + r8_abs ( e[ls-2] );
        }

        ztest = test + r8_abs ( s[ls-1] );

        if ( ztest == test )
        {
          s[ls-1] = 0.0;
          break;
        }

      }

      if ( ls == l )
      {
        kase = 3;
      }
      else if ( ls == mn )
      {
        kase = 1;
      }
      else
      {
        kase = 2;
        l = ls;
      }
    }

    l = l + 1;
//
//  Deflate negligible S(MN).
//
    if ( kase == 1 )
    {
      mm1 = mn - 1;
      f = e[mn-2];
      e[mn-2] = 0.0;

      for ( kk = 1; kk <= mm1; kk++ )
      {
        k = mm1 - kk + l;
        t1 = s[k-1];
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = t1;

        if ( k != l )
        {
          f = -sn * e[k-2];
          e[k-2] = cs * e[k-2];
        }

        if ( wantv )
        {
          drot ( n, v+0+(k-1)*ldv, 1, v+0+(mn-1)*ldv, 1, cs, sn );
        }
      }
    }
//
//  Split at negligible S(L).
//
    else if ( kase == 2 )
    {
      f = e[l-2];
      e[l-2] = 0.0;

      for ( k = l; k <= mn; k++ )
      {
        t1 = s[k-1];
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = t1;
        f = - sn * e[k-1];
        e[k-1] = cs * e[k-1];
        if ( wantu )
        {
          drot ( m, u+0+(k-1)*ldu, 1, u+0+(l-2)*ldu, 1, cs, sn );
        }
      }
    }
//
//  Perform one QR step.
//
    else if ( kase == 3 )
    {
//
//  Calculate the shift.
//
      scale = r8_max ( r8_abs ( s[mn-1] ), 
              r8_max ( r8_abs ( s[mn-2] ), 
              r8_max ( r8_abs ( e[mn-2] ), 
              r8_max ( r8_abs ( s[l-1] ), r8_abs ( e[l-1] ) ) ) ) );

      sm = s[mn-1] / scale;
      smm1 = s[mn-2] / scale;
      emm1 = e[mn-2] / scale;
      sl = s[l-1] / scale;
      el = e[l-1] / scale;
      b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0;
      c = ( sm * emm1 ) * ( sm * emm1 );
      shift = 0.0;

      if ( b != 0.0 || c != 0.0 )
      {
        shift = sqrt ( b * b + c );
        if ( b < 0.0 )
        {
          shift = -shift;
        }
        shift = c / ( b + shift );
      }

      f = ( sl + sm ) * ( sl - sm ) - shift;
      g = sl * el;
//
//  Chase zeros.
//
      mm1 = mn - 1;

      for ( k = l; k <= mm1; k++ )
      {
        drotg ( &f, &g, &cs, &sn );

        if ( k != l )
        {
          e[k-2] = f;
        }

        f = cs * s[k-1] + sn * e[k-1];
        e[k-1] = cs * e[k-1] - sn * s[k-1];
        g = sn * s[k];
        s[k] = cs * s[k];

        if ( wantv )
        {
          drot ( n, v+0+(k-1)*ldv, 1, v+0+k*ldv, 1, cs, sn );
        }

        drotg ( &f, &g, &cs, &sn );
        s[k-1] = f;
        f = cs * e[k-1] + sn * s[k];
        s[k] = -sn * e[k-1] + cs * s[k];
        g = sn * e[k];
        e[k] = cs * e[k];

        if ( wantu && k < m )
        {
          drot ( m, u+0+(k-1)*ldu, 1, u+0+k*ldu, 1, cs, sn );
        }
      }
      e[mn-2] = f;
      iter = iter + 1;
    }
//
//  Convergence.
//
    else if ( kase == 4 )
    {
//
//  Make the singular value nonnegative.
//
      if ( s[l-1] < 0.0 )
      {
        s[l-1] = -s[l-1];
        if ( wantv )
        {
          dscal ( n, -1.0, v+0+(l-1)*ldv, 1 );
        }
      }
//
//  Order the singular value.
//
      for ( ; ; )
      {
        if ( l == mm )
        {
          break;
        }

        if ( s[l] <= s[l-1] )
        {
          break;
        }

        t = s[l-1];
        s[l-1] = s[l];
        s[l] = t;

        if ( wantv && l < n )
        {
          dswap ( n, v+0+(l-1)*ldv, 1, v+0+l*ldv, 1 );
        }

        if ( wantu && l < m )
        {
          dswap ( m, u+0+(l-1)*ldu, 1, u+0+l*ldu, 1 );
        }

        l = l + 1;
      }
      iter = 0;
      mn = mn - 1;
    }
  }

  return info;
}
//****************************************************************************80

void dswap ( int n, double x[], int incx, double y[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DSWAP interchanges two vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to swap.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to swap.
//
//    Input, int INCY, the increment between successive elements of Y.
//
{
  int i;
  int ix;
  int iy;
  int m;
  double temp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    m = n % 3;

    for ( i = 0; i < m; i++ )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;
    }

    for ( i = m; i < n; i = i + 3 )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;

      temp = x[i+1];
      x[i+1] = y[i+1];
      y[i+1] = temp;

      temp = x[i+2];
      x[i+2] = y[i+2];
      y[i+2] = temp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      temp = x[ix];
      x[ix] = y[iy];
      y[iy] = temp;
      ix = ix + incx;
      iy = iy + incy;
    }

  }

  return;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

string i4_to_string ( int i4, string format )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  char i4_char[80];
  string i4_string;

  sprintf ( i4_char, format.c_str ( ), i4 );

  i4_string = string ( i4_char );

  return i4_string;
}
//****************************************************************************80

int i4_uniform ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2006
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
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

void imtqlx ( int n, double d[], double e[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to 
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine. 
//
//    It has been modified to produce the product Q' * Z, where Z is an input 
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
//    The changes consist (essentialy) of applying the orthogonal transformations
//    directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, double E(N), the subdiagonal entries of the 
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, double Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( r8_abs ( e[m-1] ) <= prec * ( r8_abs ( d[m-1] ) + r8_abs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        cout << "\n";
        cout << "IMTQLX - Fatal error!\n";
        cout << "  Iteration limit exceeded\n";
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r =  sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + r8_abs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( r8_abs ( g ) <= r8_abs ( f ) )
        {
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
//
//  Sorting.
//
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
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
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( - x + 0.5 );
  }
  else
  {
    value =   ( int ) (  x + 0.5 );
  }
  return value;
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

double r8_add ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ADD adds two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the numbers to be added.
//
//    Output, double R8_ADD, the sum of X and Y.
//
{
  double value;

  value = x + y;

  return value;
}
//****************************************************************************80

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  double value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
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
//    11 August 2010
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
  double one;
  double temp;
  double test;
  double value;

  one = ( double ) ( 1 );

  value = one;
  temp = value / 2.0;
  test = r8_add ( one, temp );

  while ( one < test )
  {
    value = temp;
    temp = value / 2.0;
    test = r8_add ( one, temp );
  }
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

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

double *r8mat_copy_new ( int m, int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
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
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
//
{
  double *a2;
  int i;
  int j;

  a2 = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
}
//****************************************************************************80

double *r8mat_mtv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MTV_NEW multiplies a transposed matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
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
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8MAT_MTV_NEW[N], the product A'*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[n];

  for ( j = 0; j < n; j++ )
  {
    y[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      y[j] = y[j] + a[i+j*m] * x[i];
    }
  }

  return y;
}
//****************************************************************************80

double *r8mat_mv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV_NEW multiplies a matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
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
//    Output, double R8MAT_MV_NEW[M], the product A*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

bool r8vec_in_ab ( int n, double x[], double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_IN_AB is TRUE if the entries of an R8VEC are in the range [A,B].
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
//    15 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, double X[N], the vector
//
//    Input, double A, B, the limits of the range.
//
//    Output, bool R8VEC_IN_AB, is TRUE if every entry is
//    between A and B.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < a || b < x[i] )
    {
      return false;
    }
  }
  return true;
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

double *r8vec_uniform_new ( int n, double b, double c, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_NEW returns a scaled pseudorandom R8VEC.
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
//    30 January 2005
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
//    Input, double B, C, the lower and upper limits of the pseudorandom values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = b + ( c - b ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int *seed )

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
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void r8vec2_print ( int n, double a1[], double a2[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_PRINT prints an R8VEC2.
//
//  Discussion:
//
//    An R8VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A1[N], double A2[N], the vectors to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    cout << setw(6)  << i
         << ": " << setw(14) << a1[i]
         << "  " << setw(14) << a2[i] << "\n";
  }

  return;
}
//****************************************************************************80

double *svd_solve ( int m, int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    SVD_SOLVE solves a linear system in the least squares sense.
//
//  Discussion:
//
//    The vector X returned by this routine should always minimize the 
//    Euclidean norm of the residual ||A*x-b||.
//
//    If the matrix A does not have full column rank, then there are multiple
//    vectors that attain the minimum residual.  In that case, the vector
//    X returned by this routine is the unique such minimizer that has the 
//    the minimum possible Euclidean norm, that is, ||A*x-b|| and ||x||
//    are both minimized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input, double A[M*N], the matrix.
//
//    Input, double B[M], the right hand side.
//
//    Output, double SVD_SOLVE[N], the least squares solution.
//
{
  double *a_copy;
  double *e;
  int i;
  int info;
  int lda;
  int ldu;
  int ldv;
  int job;
  double *sdiag;
  double smax;
  double stol;
  double *sub;
  double *u;
  double *ub;
  double *v;
  double *work;
  double *x;
//
//  Get the SVD.
//
  a_copy = r8mat_copy_new ( m, n, a );
  lda = m;
  sdiag = new double[i4_max ( m + 1, n )];
  e = new double[i4_max ( m + 1, n )];
  u = new double[m*m];
  ldu = m;
  v = new double[n*n];
  ldv = n;
  work = new double[m];
  job = 11;

  info = dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job );

  delete [] work;

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "SVD_SOLVE - Failure!\n";
    cerr << "  The SVD could not be calculated.\n";
    cerr << "  LINPACK routine DSVDC returned a nonzero\n";
    cerr << "  value of the error flag, INFO = " << info << "\n";
    exit ( 1 );
  }

  ub = r8mat_mtv_new ( m, m, u, b );
//
//  For singular problems, there may be tiny but nonzero singular values
//  that should be ignored.  This is a reasonable attempt to avoid such 
//  problems, although in general, the user might wish to control the tolerance.
//
  smax = r8vec_max ( n, sdiag );
  if ( smax <= r8_epsilon ( ) )
  {
    smax = 1.0;
  }

  stol = r8_epsilon ( ) * smax;

  sub = new double[n];

  for ( i = 0; i < n; i++ )
  {
    sub[i] = 0.0;
    if ( i < m )
    {
      if ( stol <= sdiag[i] )
      {
        sub[i] = ub[i] / sdiag[i];
      }
    }
  }

  x = r8mat_mv_new ( n, n, v, sub );

  delete [] a_copy;
  delete [] e;
  delete [] sdiag;
  delete [] sub;
  delete [] u;
  delete [] ub;
  delete [] v;

  return x;
}
//****************************************************************************80

double t_double_product_integral ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    T_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) T(i,x)*T(j,x)/sqrt(1-x^2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the polynomial indices.
//    0 <= I, J.
//
//    Output, double T_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( i < 0 )
  {
    cout << "\n";
    cout << "T_DOUBLE_PRODUCT_INTEGRAL - Fatal error!";
    cout << "  0 <= I, is required.\n";
    exit ( 1 );
  }
  if ( j < 0 )
  {
    cout << "\n";
    cout << "T_DOUBLE_PRODUCT_INTEGRAL - Fatal error!";
    cout << "  0 <= J is required.\n";
    exit ( 1 );
  }
  if ( i != j )
  {
    value = 0.0;
  }
  else if ( i == 0 )
  {
    value = pi;
  }
  else if ( 0 < i )
  {
    value = pi / 2.0;
  }
  return value;
}
//****************************************************************************80

double t_integral ( int e )

//****************************************************************************80
//
//  Purpose:
//
//    T_INTEGRAL: integral ( -1 <= x <= +1 ) x^e dx / sqrt ( 1 - x^2 ).
//
//  Discussion:
//
//    Set 
//      x = cos ( theta ), 
//      dx = - sin ( theta ) d theta = - sqrt ( 1 - x^2 ) d theta
//    to transform the integral to
//      integral ( 0 <= theta <= pi ) - ( cos ( theta ) )^e d theta
//    which becomes
//      0 if E is odd,
//      (1/2^e) * choose ( e, e/2 ) * pi if E is even.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int E, the exponent of X.
//    0 <= E.
//
//    Output, double T_INTEGRAL, the value of the integral.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( ( e % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = r8_choose ( e, e / 2 ) * pi / pow ( 2.0, e );
  }

  return value;
}
//****************************************************************************80

double *t_polynomial ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL evaluates Chebyshev polynomials T(n,x).
//
//  Discussion:
//
//    Chebyshev polynomials are useful as a basis for representing the
//    approximation of functions since they are well conditioned, in the sense
//    that in the interval [-1,1] they each have maximum absolute value 1.
//    Hence an error in the value of a coefficient of the approximation, of
//    size epsilon, is exactly reflected in an error of size epsilon between
//    the computed approximation and the theoretical approximation.
//
//    Typical usage is as follows, where we assume for the moment
//    that the interval of approximation is [-1,1].  The value
//    of N is chosen, the highest polynomial to be used in the
//    approximation.  Then the function to be approximated is
//    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
//    Chebyshev polynomial.  Let these values be denoted by F(XJ).
//
//    The coefficients of the approximation are now defined by
//
//      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I,XJ)
//
//    except that C(0) is given a value which is half that assigned
//    to it by the above formula,
//
//    and the representation is
//
//    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
//
//    Now note that, again because of the fact that the Chebyshev polynomials
//    have maximum absolute value 1, if the higher order terms of the
//    coefficients C are small, then we have the option of truncating
//    the approximation by dropping these terms, and we will have an
//    exact value for maximum perturbation to the approximation that
//    this will cause.
//
//    It should be noted that typically the error in approximation
//    is dominated by the first neglected basis function (some multiple of
//    T(N+1,X) in the example above).  If this term were the exact error,
//    then we would have found the minimax polynomial, the approximating
//    polynomial of smallest maximum deviation from the original function.
//    The minimax polynomial is hard to compute, and another important
//    feature of the Chebyshev approximation is that it tends to behave
//    like the minimax polynomial while being easy to compute.
//
//    To evaluate a sum like
//
//      sum ( 0 <= J <= N ) C(J) T(J,X),
//
//    Clenshaw's recurrence formula is recommended instead of computing the
//    polynomial values, forming the products and summing.
//
//    Assuming that the coefficients C(J) have been computed
//    for J = 0 to N, then the coefficients of the representation of the
//    indefinite integral of the function may be computed by
//
//      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1,
//
//    with
//
//      C(N+1)=0
//      B(0) arbitrary.
//
//    Also, the coefficients of the representation of the derivative of the
//    function may be computed by:
//
//      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0,
//
//    with
//
//      D(N+1) = D(N)=0.
//
//    Some of the above may have to adjusted because of the irregularity of C(0).
//
//    The formula is:
//
//      T(N,X) = COS(N*ARCCOS(X))
//
//  Differential equation:
//
//    (1-X*X) Y'' - X Y' + N N Y = 0
//
//  First terms:
//
//    T(0,X) =  1
//    T(1,X) =  1 X
//    T(2,X) =  2 X^2 -   1
//    T(3,X) =  4 X^3 -   3 X
//    T(4,X) =  8 X^4 -   8 X^2 +  1
//    T(5,X) = 16 X^5 -  20 X^3 +  5 X
//    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
//    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
//
//  Inequality:
//
//    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
//
//  Orthogonality:
//
//    For integration over [-1,1] with weight
//
//      W(X) = 1 / sqrt(1-X*X),
//
//    if we write the inner product of T(I,X) and T(J,X) as
//
//      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
//
//    then the result is:
//
//      0 if I /= J
//      PI/2 if I == J /= 0
//      PI if I == J == 0
//
//    A discrete orthogonality relation is also satisfied at each of
//    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
//                              = 0 if I /= J
//                              = N/2 if I == J /= 0
//                              = N if I == J == 0
//
//  Recursion:
//
//    T(0,X) = 1,
//    T(1,X) = X,
//    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
//
//    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
//
//  Special values:
//
//    T(N,1) = 1
//    T(N,-1) = (-1)^N
//    T(2N,0) = (-1)^N
//    T(2N+1,0) = 0
//    T(N,X) = (-1)**N * T(N,-X)
//
//  Zeroes:
//
//    M-th zero of T(N,X) is cos((2*M-1)*PI/(2*N)), M = 1 to N
//
//  Extrema:
//
//    M-th extremum of T(N,X) is cos(PI*M/N), M = 0 to N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest polynomial to compute.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double T_POLYNOMIAL[M*(N+1)], the values of the Chebyshev polynomials.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i] = 1.0;
  }
  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = 2.0 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
    }
  }

  return v;
}
//****************************************************************************80

double *t_polynomial_ab ( double a, double b, int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_AB: Chebyshev polynomials T(n,x) in [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the domain of definition.
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest polynomial to compute.
//
//    Input, double X[M], the evaluation points.
//    It must be the case that A <= X(*) <= B.
//
//    Output, double T_POLYNOMIAL_AB[M*(N+1)], the values.
//
{
  int i;
  double *v;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = ( ( b - x[i]     )   
           - (     x[i] - a ) ) 
           / ( b        - a );
  }

  v = t_polynomial ( m, n, x );
 
  return v;
}
//****************************************************************************80

double *t_polynomial_coefficients ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_COEFFICIENTS: coefficients of the Chebyshev polynomial T(n,x).
//
//  First terms:
//
//    N/K     0     1      2      3       4     5      6    7      8    9   10
//
//     0      1
//     1      0     1
//     2     -1     0      2
//     3      0    -3      0      4
//     4      1     0     -8      0       8
//     5      0     5      0    -20       0    16
//     6     -1     0     18      0     -48     0     32
//     7      0    -7      0     56       0  -112      0    64
//
//  Recursion:
//
//    T(0,X) = 1,
//    T(1,X) = X,
//    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Output, double T_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients 
//    of the Chebyshev T polynomials.
//
{
  double *c;
  int i;
  int j;

  if ( n < 0 )
  {
    return NULL;
  }

  c = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n == 0 )
  {
    return c;
  }

  c[1+1*(n+1)] = 1.0;

  for ( i = 2; i <= n; i++ )
  {
    c[i+0*(n+1)]     =                    - c[i-2+0*(n+1)];
    for ( j = 1; j <= i - 2; j++ )
    {
      c[i+j*(n+1)] = 2.0 * c[i-1+(j-1)*(n+1)] - c[i-2+j*(n+1)];
    }
    c[i+(i-1)*(n+1)] = 2.0 * c[i-1+(i-2)*(n+1)];
    c[i+ i   *(n+1)] = 2.0 * c[i-1+(i-1)*(n+1)];
  }

  return c;
}
//****************************************************************************80

double t_polynomial_value ( int n, double x )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_VALUE: returns the single value T(n,x).
//
//  Discussion:
//
//    In cases where calling T_POLYNOMIAL is inconvenient, because it returns
//    a vector of values for multiple arguments X, this simpler interface
//    may be appropriate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double X, the argument of the polynomial.
//
//    Output, double T_POLYNOMIAL_VALUE, the value of T(n,x).
//
{
  int m;
  double *v_vec;
  double value;
  double x_vec[1];

  m = 1;
  x_vec[0] = x;

  v_vec = t_polynomial ( m, n, x_vec );

  value = v_vec[n];

  delete [] v_vec;

  return value;
}
//****************************************************************************80

void t_polynomial_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_VALUES returns values of the Chebyshev polynomial T(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ChebyshevT[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.8000000000000000E+00,
      0.2800000000000000E+00,
     -0.3520000000000000E+00,
     -0.8432000000000000E+00,
     -0.9971200000000000E+00,
     -0.7521920000000000E+00,
     -0.2063872000000000E+00,
      0.4219724800000000E+00,
      0.8815431680000000E+00,
      0.9884965888000000E+00,
      0.7000513740800000E+00,
      0.1315856097280000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *t_polynomial_zeros ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial T(n,x).
//
//  Discussion:
//
//    The I-th zero is cos((2*I-1)*PI/(2*N)), I = 1 to N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double T_POLYNOMIAL_ZEROS[N], the zeroes.
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;
  double *z;

  z = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    angle = ( double) ( 2 * i - 1 ) * pi / ( double ) ( 2 * n );
    z[i-1] = cos ( angle );
  }

  return z;
}
//****************************************************************************80

double *t_project_coefficients ( int n, double f ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    T_PROJECT_COEFFICIENTS: function projected onto Chebyshev polynomials T(n,x).
//
//  Discussion:
//
//    It is assumed that the interval of definition is -1 <= x <= +1.
//
//    Over this interval, f(x) will be well approximated by
//
//      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,x)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//
//    Input, double F ( double X ), evaluates the function.
//
//    Output, double T_PROJECT_COEFFICIENTS[N+1], the projection coefficients 
//    of f(x) onto T(0,x) through T(n,x).
//
{
  double *c;
  double *d;
  double fac;
  int j;
  int k;
  double pi = 3.141592653589793;
  double total;
  double y;

  d = new double[n+1];

  for ( k = 0; k <= n; k++ )
  {
    y = cos ( pi * ( ( double ) ( k ) + 0.5 ) / ( double ) ( n + 1 ) );
    d[k] = f ( y );
  }

  fac = 2.0 / ( double ) ( n + 1 );

  c = new double[n+1];

  for ( j = 0; j <= n; j++ )
  {
    total = 0.0;
    for ( k = 0; k <= n; k++ )
    {
      total = total + d[k] * cos ( ( pi * ( double ) ( j ) ) 
        * ( ( ( double ) ( k ) + 0.5 ) / ( double ) ( n + 1 ) ) );
    }
    c[j] = fac * total;
  }

  c[0] = c[0] / 2.0;

  delete [] d;

  return c;
}
//****************************************************************************80

double *t_project_coefficients_ab ( int n, double f ( double x ), double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    T_PROJECT_COEFFICIENTS_AB: function projected onto T(n,x) over [a,b]
//
//  Discussion:
//
//    It is assumed that the interval of definition is a <= x <= b.
//
//    Over this interval, f(x) will be well approximated by
//
//      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,(2x-a-b)/(b-a))
//
//    where x* = ( x - b - a ) / ( b - a )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//
//    Input, double F ( double X ), evaluates the function.
//
//    Input, double A, B, the interval of definition.
//
//    Output, double T_PROJECT_COEFFICIENTS_AB[N+1], the projection coefficients 
//    of f(x) onto T(0,x) through T(n,x).
//
{
  double *c;
  double *d;
  double fac;
  int j;
  int k;
  double pi = 3.141592653589793;
  double t;
  double total;
  double y;

  d = new double[n+1];

  for ( k = 0; k <= n; k++ )
  {
    t = cos ( pi * ( ( double ) ( k ) + 0.5 ) / ( double ) ( n + 1 ) );

    y = ( ( 1.0 + t ) * b   
        + ( 1.0 - t ) * a ) 
        /   2.0;

    d[k] = f ( y );
  }

  fac = 2.0 / ( double ) ( n + 1 );

  c = new double[n+1];

  for ( j = 0; j <= n; j++ )
  {
    total = 0.0;
    for ( k = 0; k <= n; k++ )
    {
      total = total + d[k] * cos ( ( pi * ( double ) ( j ) ) 
        * ( ( ( double ) ( k ) + 0.5 ) / ( double ) ( n + 1 ) ) );
    }
    c[j] = fac * total;
  }

  c[0] = c[0] / 2.0;

  delete [] d;

  return c;
}
//****************************************************************************80

double *t_project_coefficients_data ( double a, double b, int m, int n, 
  double x[], double d[] )

//****************************************************************************80
//
//  Purpose:
//
//    T_PROJECT_COEFFICIENTS_DATA: project data onto Chebyshev polynomials T(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the domain of definition.
//
//    Input, int M, the number of data values.
//
//    Input, int N, the desired order of the Chebyshev 
//    expansion.
//
//    Input, double X[M], the data abscissas.  These need not 
//    be sorted.  It must be the case that A <= X() <= B.
//
//    Input, double D[M], the data values.
//
//    Output, double T_PROJECT_COEFFICIENTS_DATA[N+1], the approximate 
//    Chebshev coefficients.
//
{
  double *c;
  double *v;

  if ( !r8vec_in_ab ( m, x, a, b ) )
  {
    cerr << "\n";
    cerr << "T_PROJECT_COEFFICIENTS_DATA- Fatal error!\n";
    cerr << "  Some X not in [A,B].\n";
    exit ( 1 );
  }
//
//  Compute the M by N+1 Chebyshev Vandermonde matrix V.
//
  v = t_polynomial_ab ( a, b, m, n, x );
//
//  Compute the least-squares solution C.
//
  c = svd_solve ( m, n + 1, v, d );

  return c;
}
//****************************************************************************80

double *t_project_value ( int m, int n, double x[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    T_PROJECT_VALUE evaluates an expansion in Chebyshev polynomials T(n,x).
//
//  Discussion:
//
//    The projection is assumed to be based on the interval [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//
//    Input, double X[M], the evaluation points.
//
//    Input, double C[N+1], the expansion coefficients.
//
//    Output, double T_PROJECT_VALUE[M], the value of the Chebyshev function.
//
{
  double *b0;
  double *b1;
  double *b2;
  int i;
  int j;
  double *v;

  b0 = new double[m];
  b1 = new double[m];
  b2 = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b0[i] = 0.0;
  }
  for ( i = 0; i < m; i++ )
  {
    b1[i] = 0.0;
  }

  for ( j = n; 0 <= j; j-- )
  {
    for ( i = 0; i < m; i++ )
    {
      b2[i] = b1[i];
      b1[i] = b0[i];
      b0[i] = c[j] + 2.0 * x[i] * b1[i] - b2[i];
    }
  }

  v = new double[m];

  for ( i = 0; i < m; i++ )
  {
    v[i] = 0.5 * ( c[0] + b0[i] - b2[i] );
  }
  return v;
}
//****************************************************************************80

double *t_project_value_ab ( int m, int n, double x[], double c[], double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    T_PROJECT_VALUE_AB evaluates an expansion in Chebyshev polynomials T(n,x).
//
//  Discussion:
//
//    The projection is assumed to be based on the interval [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//
//    Input, double X[M], the evaluation points.
//
//    Input, double C[N+1], the expansion coefficients.
//
//    Input, double A, B, the interval of definition.
//
//    Output, double T_PROJECT_VALUE_AB[M], the value of the Chebyshev function.
//
{
  double *b0;
  double *b1;
  double *b2;
  int i;
  int j;
  double *v;

  b0 = new double[m];
  b1 = new double[m];
  b2 = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b0[i] = 0.0;
  }
  for ( i = 0; i < m; i++ )
  {
    b1[i] = 0.0;
  }

  for ( j = n; 0 <= j; j-- )
  {
    for ( i = 0; i < m; i++ )
    {
      b2[i] = b1[i];
      b1[i] = b0[i];
      b0[i] = c[j] + 2.0 / ( b - a ) * ( 2.0 * x[i] - a - b ) * b1[i] - b2[i];
    }
  }

  v = new double[m];

  for ( i = 0; i < m; i++ )
  {
    v[i] = 0.5 * ( c[0] + b0[i] - b2[i] );
  }
  return v;
}
//****************************************************************************80

void t_quadrature_rule ( int n, double t[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    T_QUADRATURE_RULE: quadrature rule for T(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Output, double T[N], W[N], the points and weights of the rule.
//
{
  double *bj;
  int i;
  double pi = 3.141592653589793;

  for ( i = 0; i < n; i++ )
  {
    t[i] = 0.0;
  }

  bj = new double[n];
  bj[0] = sqrt ( 0.5 );
  for ( i = 1; i < n; i++ )
  {
    bj[i] = 0.5;
  }

  w[0] = sqrt ( pi );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
  imtqlx ( n, t, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }
  return;
}
//****************************************************************************80

double t_triple_product_integral ( int i, int j, int k )

//****************************************************************************80
//
//  Purpose:
//
//    T_TRIPLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) T(i,x)*T(j,x)*T(k,x)/sqrt(1-x^2) dx
//
//  Discussion:
//
//    
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Mason, David Handscomb,
//    Chebyshev Polynomials,
//    CRC Press, 2002,
//    ISBN: 0-8493-035509,
//    LC: QA404.5.M37.
//
//  Parameters:
//
//    Input, int I, J, K, the polynomial indices.
//    0 <= I, J.
//
//    Output, double T_TRIPLE_PRODUCT_INTEGRAL, the integral.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( i < 0 )
  {
    cerr << "\n";
    cerr << "T_TRIPLE_PRODUCT_INTEGRAL - Fatal error!\n";
    cerr << "  0 <= I is required.\n";
    exit ( 1 );
  }

  if ( j < 0 )
  {
    cerr << "\n";
    cerr << "T_TRIPLE_PRODUCT_INTEGRAL - Fatal error!\n";
    cerr << "  0 <= J is required.\n";
    exit ( 1 );
  }

  if ( k < 0 )
  {
    cerr << "\n";
    cerr << "T_TRIPLE_PRODUCT_INTEGRAL - Fatal error!\n";
    cerr << "  0 <= K is required.\n";
    exit ( 1 );
  }

  value = 0.5 * ( 
      t_double_product_integral (       i + j,   k ) + 
    + t_double_product_integral ( abs ( i - j ), k ) );

  return value;
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

double u_double_product_integral ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    U_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) U(i,x)*U(j,x)*sqrt(1-x^2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the polynomial indices.
//    0 <= I, J.
//
//    Output, double U_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( i < 0 )
  {
    cout << "\n";
    cout << "U_DOUBLE_PRODUCT_INTEGRAL - Fatal error!";
    cout << "  0 <= I, is required.\n";
    exit ( 1 );
  }
  if ( j < 0 )
  {
    cout << "\n";
    cout << "U_DOUBLE_PRODUCT_INTEGRAL - Fatal error!";
    cout << "  0 <= J is required.\n";
    exit ( 1 );
  }
  if ( i != j )
  {
    value = 0.0;
  }
  else
  {
    value = pi / 2.0;
  }
  return value;
}
//****************************************************************************80

double u_integral ( int e )

//****************************************************************************80
//
//  Purpose:
//
//    U_INTEGRAL: integral ( -1 <= x <= +1 ) x^e sqrt ( 1 - x^2 ) dx.
//
//  Discussion:
//
//     E    U_INTEGRAL
//    --    -------------- 
//     0         pi /    2 
//     2         pi /    8
//     4         pi /   16
//     6     5 * pi /  128
//     8     7 * pi /  256
//    10    21 * pi / 1024
//    12    33 * pi / 2048
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int E, the exponent of X.
//    0 <= E.
//
//    Output, double U_INTEGRAL, the value of the integral.
//
{
  double arg1;
  double arg2;
  double pi = 3.141592653589793;
  double value;

  if ( ( e % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    arg1 = 0.5 * ( double ) ( 1 + e );
    arg2 = 2.0 + 0.5 * ( double ) ( e );
    value = 0.5 * sqrt ( pi ) * gamma ( arg1 ) / gamma ( arg2 );
  }

  return value;
}
//****************************************************************************80

double *u_polynomial ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL evaluates Chebyshev polynomials U(n,x).
//
//  Differential equation:
//
//    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
//
//  First terms:
//
//    U(0,X) =   1
//    U(1,X) =   2 X
//    U(2,X) =   4 X^2 -   1
//    U(3,X) =   8 X^3 -   4 X
//    U(4,X) =  16 X^4 -  12 X^2 +  1
//    U(5,X) =  32 X^5 -  32 X^3 +  6 X
//    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
//    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
//
//  Recursion:
//
//    U(0,X) = 1,
//    U(1,X) = 2 * X,
//    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
//
//  Norm:
//
//    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = PI/2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest polynomial to compute.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double U_POLYNOMIAL[M*(N+1)], the values of the N+1 Chebyshev polynomials.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = 2.0 * x[i];
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 2; j <= n; j++ )
    {
      v[i+j*m] = 2.0 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
    }
  }
  return v;
}
//****************************************************************************80

double *u_polynomial_coefficients ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_COEFFICIENTS evaluates coefficients of Chebyshev polynomials U(n,x).
//
//  First terms:
//
//    N/K     0     1      2      3       4     5      6    7      8    9   10
//
//     0      1
//     1      0     2
//     2     -1     0      4
//     3      0    -4      0      8
//     4      1     0    -12      0      16
//     5      0     6      0    -32       0    32
//     6     -1     0     24      0     -80     0     64
//     7      0    -8      0     80       0  -192      0   128
//
//  Recursion:
//
//    U(0,X) = 1,
//    U(1,X) = 2*X,
//    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Output, double U_POLYNOMIAL_COEFFICIENTS[(N+1)*((N+1)], the coefficients 
//    of the Chebyshev U polynomials.
//
{
  double *c;
  int i;
  int j;

  if ( n < 0 )
  {
    return NULL;
  }

  c = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n == 0 )
  {
    return c;
  }

  c[1+1*(n+1)] = 2.0;

  for ( i = 2; i <= n; i++ )
  {
    c[i+0*(n+1)]     =                    - c[i-2+0*(n+1)];
    for ( j = 1; j <= i-2; j++ )
    {
      c[i+j*(n+1)] = 2.0 * c[i-1+(j-1)*(n+1)] - c[i-2+j*(n+1)];
    }

    c[i+(i-1)*(n+1)] = 2.0 * c[i-1+(i-2)*(n+1)];
    c[i+ i   *(n+1)] = 2.0 * c[i-1+(i-1)*(n+1)];
  }

  return c;
}
//****************************************************************************80

void u_polynomial_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_VALUES returns values of Chebyshev polynomials U(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ChebyshevU[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.1600000000000000E+01,
      0.1560000000000000E+01,
      0.8960000000000000E+00,
     -0.1264000000000000E+00,
     -0.1098240000000000E+01,
     -0.1630784000000000E+01,
     -0.1511014400000000E+01,
     -0.7868390400000000E+00,
      0.2520719360000000E+00,
      0.1190154137600000E+01,
      0.1652174684160000E+01,
      0.1453325357056000E+01 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *u_polynomial_zeros ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_ZEROS returns zeroes of Chebyshev polynomials U(n,x).
//
//  Discussion:
//
//    The I-th zero is cos((I-1)*PI/(N-1)), I = 1 to N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double U_POLYNOMIAL_ZEROS[N], the zeroes.
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;
  double *z;

  z = new double[n];

  for ( i = 1; i <= n; i++ )
  {
    angle = ( double) ( i ) * pi / ( double ) ( n + 1 );
    z[i-1] = cos ( angle );
  }

  return z;
}
//****************************************************************************80

void u_quadrature_rule ( int n, double t[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    U_QUADRATURE_RULE: quadrature rule for U(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Output, double T[N], W[N], the points and weights of the rule.
//
{
  double *bj;
  int i;
  double pi = 3.141592653589793;

  for ( i = 0; i < n; i++ )
  {
    t[i] = 0.0;
  }

  bj = new double[n];
  for ( i = 0; i < n; i++ )
  {
    bj[i] = 0.5;
  }

  w[0] = sqrt ( pi / 2.0 );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
  imtqlx ( n, t, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }
  return;
}
//****************************************************************************80

double v_double_product_integral ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    V_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) V(i,x)*V(j,x)*sqrt(1+x)/sqrt(1-x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the polynomial indices.
//    0 <= I, J.
//
//    Output, double V_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( i < 0 )
  {
    cout << "\n";
    cout << "V_DOUBLE_PRODUCT_INTEGRAL - Fatal error!";
    cout << "  0 <= I, is required.\n";
    exit ( 1 );
  }
  if ( j < 0 )
  {
    cout << "\n";
    cout << "V_DOUBLE_PRODUCT_INTEGRAL - Fatal error!";
    cout << "  0 <= J is required.\n";
    exit ( 1 );
  }
  if ( i != j )
  {
    value = 0.0;
  }
  else
  {
    value = pi;
  }
  return value;
}
//****************************************************************************80

double *v_polynomial ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL evaluates Chebyshev polynomials V(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest polynomial to compute.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double V_POLYNOMIAL[M*(N+1)], the values of the N+1 Chebyshev polynomials.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = 2.0 * x[i] - 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 2; j <= n; j++ )
    {
      v[i+j*m] = 2.0 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
    }
  }
  return v;
}
//****************************************************************************80

void v_polynomial_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL_VALUES returns values of Chebyshev polynomials V(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      u = Sqrt[(x+1)/2],
//      ChebyshevT[2*n+1,u] / u
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
     1.0000000000000000E+00, 
     0.6000000000000000E+00, 
    -0.0400000000000000E+00, 
    -0.6640000000000000E+00, 
    -1.0224000000000000E+00, 
    -0.9718400000000000E+00, 
    -0.5325440000000000E+00, 
     0.1197696000000000E+00, 
     0.7241753600000000E+00, 
     1.0389109760000000E+00, 
     0.9380822016000000E+00, 
     0.4620205465600000E+00, 
    -0.1988493271040000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *v_polynomial_zeros ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial V(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double V_POLYNOMIAL_ZEROS[N], the zeroes.
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;
  double *z;

  z = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = ( double) ( 2 * n - 2 * i - 1 ) * pi / ( double ) ( 2 * n + 1 );
    z[i] = cos ( angle );
  }

  return z;
}
//****************************************************************************80

double w_double_product_integral ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    W_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) W(i,x)*W(j,x)*sqrt(1-x)/sqrt(1+x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the polynomial indices.
//    0 <= I, J.
//
//    Output, double W_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( i < 0 )
  {
    cout << "\n";
    cout << "W_DOUBLE_PRODUCT_INTEGRAL - Fatal error!";
    cout << "  0 <= I, is required.\n";
    exit ( 1 );
  }
  if ( j < 0 )
  {
    cout << "\n";
    cout << "W_DOUBLE_PRODUCT_INTEGRAL - Fatal error!";
    cout << "  0 <= J is required.\n";
    exit ( 1 );
  }
  if ( i != j )
  {
    value = 0.0;
  }
  else
  {
    value = pi;
  }
  return value;
}
//****************************************************************************80

double *w_polynomial ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL evaluates Chebyshev polynomials W(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest polynomial to compute.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double W_POLYNOMIAL[M*(N+1)], the values of the N+1 Chebyshev polynomials.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = 2.0 * x[i] + 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 2; j <= n; j++ )
    {
      v[i+j*m] = 2.0 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
    }
  }
  return v;
}
//****************************************************************************80

void w_polynomial_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL_VALUES returns values of Chebyshev polynomials W(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      u = Sqrt[(x+1)/2],
//      ChebyshevU[2*n,u]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
     1.000000000000000E+00, 
     2.600000000000000E+00, 
     3.160000000000000E+00, 
     2.456000000000000E+00, 
     0.769600000000000E+00, 
    -1.224640000000000E+00, 
    -2.729024000000000E+00, 
    -3.141798400000000E+00, 
    -2.297853440000000E+00, 
    -0.534767104000000E+00, 
     1.442226073600000E+00, 
     2.842328821760000E+00, 
     3.105500041216000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *w_polynomial_zeros ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial W(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double W_POLYNOMIAL_ZEROS[N], the zeroes.
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;
  double *z;

  z = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = ( double) ( 2 * ( n - i ) * pi / ( double ) ( 2 * n + 1 );
    z[i] = cos ( angle );
  }

  return z;
}
