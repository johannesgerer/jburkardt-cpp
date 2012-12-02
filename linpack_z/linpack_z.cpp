# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "blas1_z.hpp"
# include "linpack_z.hpp"

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
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
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

int zchdc ( complex <double> a[], int lda, int p, int ipvt[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZCHDC: Cholesky decomposition of a Hermitian positive definite matrix.
//
//  Discussion:
//
//    A pivoting option allows the user to estimate the condition of a
//    Hermitian positive definite matrix or determine the rank of a
//    Hermitian positive semidefinite matrix.
//
//    For Hermitian positive definite matrices, INFO = P is the normal return.
//
//    For pivoting with Hermitian positive semidefinite matrices, INFO will
//    in general be less than P.  However, INFO may be greater than
//    the rank of A, since rounding error can cause an otherwise zero
//    element to be positive.  Indefinite systems will always cause
//    INFO to be less than P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*P].  On input, A contains the matrix
//    whose decomposition is to be computed.  Only the upper half of A
//    need be stored.  The lower part of the array A is not referenced.
//    On output, A contains in its upper half the Cholesky factor
//    of the matrix A as it has been permuted by pivoting.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int P, the order of the matrix.
//
//    Input/output, int IPVT[P].  IPVT is not referenced if JOB == 0.
//    On input, IPVT contains integers that control the selection of the
//    pivot elements, if pivoting has been requested.  Each diagonal element
//    A(K,K) is placed in one of three classes according to the input
//    value of IPVT(K):
//      IPVT(K) >  0, X(K) is an initial element.
//      IPVT(K) == 0, X(K) is a free element.
//      IPVT(K) <  0, X(K) is a final element.
//    Before the decomposition is computed, initial elements are moved by
//    symmetric row and column interchanges to the beginning of the array A
//    and final elements to the end.  Both initial and final elements
//    are frozen in place during the computation and only free elements
//    are moved.  At the K-th stage of the reduction, if A(K,K) is occupied
//    by a free element, it is interchanged with the largest free element
//    A(L,L) with K <= L.
//    On output, IPVT(K) contains the index of the diagonal element
//    of A that was moved into the J-th position, if pivoting was requested.
//
//    Input, int JOB, specifies whether column pivoting is to be done.
//    0, no pivoting is done.
//    nonzero, pivoting is done.
//
//    Output, int ZCHDC, contains the index of the last positive
//    diagonal element of the Cholesky factor.
//
{
  int i_temp;
  int info;
  int j;
  int k;
  int kb;
  int l;
  double maxdia;
  int maxl;
  bool negk;
  int pl;
  int plp1;
  int pu;
  bool swapk;
  complex <double> temp;
  complex <double> *work;

  pl = 1;
  pu = 0;
  info = p;

  work = new complex <double> [p];

  if ( job != 0 )
  {
//
//  Pivoting has been requested.  Rearrange the elements according to IPVT.
//
    for ( k = 1; k <= p; k++ )
    {
      swapk = ( 0 < ipvt[k-1] );
      negk = ( ipvt[k-1] < 0 );

      if ( negk )
      {
        ipvt[k-1] = -k;
      }
      else
      {
        ipvt[k-1] = k;
      }

      if ( swapk )
      {
        if ( k != pl )
        {
          zswap ( pl-1, a+0+(k-1)*lda, 1, a+0+(pl-1)*lda, 1 );

          temp               = a[k-1+(k-1)*lda];
          a[k-1+(k-1)*lda]   = a[pl-1+(pl-1)*lda];
          a[pl-1+(pl-1)*lda] = temp;

          a[pl-1+(k-1)*lda] = conj ( a[pl-1+(k-1)*lda] );
          plp1 = pl + 1;

          for ( j = plp1; j <= p; j++ )
          {
            if ( j < k )
            {
              temp              = conj ( a[pl-1+(j-1)*lda] );
              a[pl-1+(j-1)*lda] = conj ( a[j-1+(k-1)*lda] );
              a[j-1+(k-1)*lda]  = temp;
            }
            else if ( j != k )
            {
              temp              = a[pl-1+(j-1)*lda];
              a[pl-1+(j-1)*lda] = a[k-1+(j-1)*lda];
              a[k-1+(j-1)*lda]  = temp;
            }
          }
          ipvt[k-1] = ipvt[pl-1];
          ipvt[pl-1] = k;
        }
        pl = pl + 1;
      }
    }

    pu = p;

    for ( kb = pl; kb <= p; kb++ )
    {
      k = p - kb + pl;

      if ( ipvt[k-1] < 0 )
      {
        ipvt[k-1] = -ipvt[k-1];

        if ( pu != k )
        {
          zswap ( k-1, a+0+(k-1)*lda, 1, a+0+(pu-1)*lda, 1 );

          temp               = a[k-1+(k-1)*lda];
          a[k-1+(k-1)*lda]   = a[pu-1+(pu-1)*lda];
          a[pu-1+(pu-1)*lda] = temp;

          a[k-1+(pu-1)*lda] = conj ( a[k-1+(pu-1)*lda] );

          for ( j = k + 1; j <= p; j++ )
          {
            if ( j < pu )
            {
              temp              = conj ( a[k-1+(j-1)*lda] );
              a[k-1+(j-1)*lda]  = conj ( a[j-1+(pu-1)*lda] );
              a[j-1+(pu-1)*lda] = temp;
            }
            else if ( j != pu )
            {
              temp              = a[k-1+(j-1)*lda];
              a[k-1+(j-1)*lda]  = a[pu-1+(j-1)*lda];
              a[pu-1+(j-1)*lda] = temp;
            }
          }
          i_temp     = ipvt[k-1];
          ipvt[k-1]  = ipvt[pu-1];
          ipvt[pu-1] = i_temp;
        }
        pu = pu - 1;
      }
    }
  }

  for ( k = 1; k <= p; k++ )
  {
//
//  Reduction loop.
//
    maxdia = real ( a[k-1+(k-1)*lda] );
    maxl = k;
//
//  Determine the pivot element.
//
    if ( pl <= k && k < pu )
    {
      for ( l = k + 1; l <= pu; l++ )
      {
        if ( maxdia < real ( a[l-1+(l-1)*lda] ) )
        {
          maxdia = real ( a[l-1+(l-1)*lda] );
          maxl = l;
        }
      }
    }
//
//  Quit if the pivot element is not positive.
//
    if ( maxdia <= 0.0 )
    {
      info = k - 1;
      delete [] work;
      return info;
    }
//
//  Start the pivoting and update IPVT.
//
    if ( k != maxl )
    {
      zswap ( k-1, a+0+(k-1)*lda, 1, a+0+(maxl-1)*lda, 1 );
      a[maxl-1+(maxl-1)*lda] = a[k-1+(k-1)*lda];
      a[k-1+(k-1)*lda] = complex <double> ( maxdia, 0.0 );

      i_temp       = ipvt[maxl-1];
      ipvt[maxl-1] = ipvt[k-1];
      ipvt[k-1]    = i_temp;

      a[k-1+(maxl-1)*lda] = conj ( a[k-1+(maxl-1)*lda] );
    }
//
//  Reduction step.  Pivoting is contained across the rows.
//
    work[k-1] = complex <double> ( sqrt ( real ( a[k-1+(k-1)*lda] ) ), 0.0 );
    a[k-1+(k-1)*lda] = work[k-1];

    for ( j = k + 1; j <= p; j++ )
    {
      if ( k != maxl )
      {
        if ( j < maxl )
        {
          temp                = conj ( a[k-1+(j-1)*lda] );
          a[k-1+(j-1)*lda]    = conj ( a[j-1+(maxl-1)*lda] );
          a[j-1+(maxl-1)*lda] = temp;
        }
        else if ( j != maxl )
        {
          temp                = a[k-1+(j-1)*lda];
          a[k-1+(j-1)*lda]    = a[maxl-1+(j-1)*lda];
          a[maxl-1+(j-1)*lda] = temp;
        }
      }
      a[k-1+(j-1)*lda] = a[k-1+(j-1)*lda] / work[k-1];
      work[j-1] = conj ( a[k-1+(j-1)*lda] );
      temp = -a[k-1+(j-1)*lda];
      zaxpy ( j-k, temp, work+k, 1, a+k+(j-1)*lda, 1 );
    }
  }
  delete [] work;
  return info;
}
//****************************************************************************80

int zcdhd ( complex <double> r[], int ldr, int p, complex <double> x[],
  complex <double> z[], int ldz, int nz, complex <double> y[], double rho[],
  double c[], complex <double> s[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZCHDD downdates an augmented Cholesky decomposition.
//
//  Discussion:
//
//    zcdhD downdates an augmented Cholesky decomposition or the
//    triangular factor of an augmented QR decomposition.
//    Specifically, given an upper triangular matrix R of order P,  a
//    row vector X, a column vector Z, and a scalar Y, ZCHDD
//    determines a unitary matrix U and a scalar ZETA such that
//
//          ( R   Z  )     ( RR  ZZ )
//      U * (        )  =  (        ),
//          ( 0 ZETA )     (  X   Y )
//
//    where RR is upper triangular.  If R and Z have been obtained
//    from the factorization of a least squares problem, then
//    RR and ZZ are the factors corresponding to the problem
//    with the observation (X,Y) removed.  In this case, if RHO
//    is the norm of the residual vector, then the norm of
//    the residual vector of the downdated problem is
//      sqrt ( RHO**2 - ZETA**2 ).
//    zcdhD will simultaneously downdate several triplets (Z,Y,RHO)
//    along with R.
//
//    For a less terse description of what ZCHDD does and how
//    it may be applied, see the LINPACK guide.
//
//    The matrix U is determined as the product U(1)*...*U(P)
//    where U(I) is a rotation in the (P+1,I)-plane of the
//    form
//
//      ( C(I)  -conj ( S(I) ) )
//      (                       ).
//      ( S(I)           C(I)   )
//
//    The rotations are chosen so that C(I) is real.
//
//    The user is warned that a given downdating problem may
//    be impossible to accomplish or may produce
//    inaccurate results.  For example, this can happen
//    if X is near a vector whose removal will reduce the
//    rank of R.  Beware.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> R[LDR*P]; on input, the upper triangular matrix
//    that is to be downdated.  On output, the downdated matrix.  The
//    part of R below the diagonal is not referenced.
//
//    Input, int LDR, the leading dimension of R.  P <= LDR.
//
//    Input, int P, the order of the matrix.
//
//    Input, complex <double> X(P), the row vector that is to
//    be removed from R.
//
//    Input/output, complex <double> Z[LDZ*NZ]; on input, an array of NZ
//    P-vectors which are to be downdated along with R.  On output,
//    the downdated vectors.
//
//    Input, int LDZ, the leading dimension of Z.  P <= LDZ.
//
//    Input, int NZ, the number of vectors to be downdated.
//    NZ may be zero, in which case Z, Y, and R are not referenced.
//
//    Input, complex <double> Y[NZ], the scalars for the downdating
//    of the vectors Z.
//
//    Input/output, double RHO[NZ].  On input, the norms of the residual
//    vectors that are to be downdated.  On output, the downdated norms.
//
//    Output, double C[P], the cosines of the transforming rotations.
//
//    Output, complex <double> S[P], the sines of the transforming rotations.
//
//    Output, int ZCHDD:
//     0, if the entire downdating was successful.
//    -1, if R could not be downdated.  In this case, all quantities
//        are left unaltered.
//     1, if some RHO could not be downdated.  The offending RHO's are
//        set to -1.
//
{
  double a;
  double alpha;
  double azeta;
  complex <double> b;
  int i;
  int ii;
  int info;
  int j;
  double norm;
  double scale;
  complex <double> t;
  complex <double> xx;
  complex <double> zeta;
//
//  Solve the system hermitian(R) * A = X, placing the result in S.
//
  info = 0;
  s[0] = conj ( x[0] ) / conj ( r[0+0*ldr] );

  for ( j = 2; j <= p; j++ )
  {
    s[j-1] = conj ( x[j-1] ) - zdotc ( j-1, r+0+(j-1)*ldr, 1, s, 1 );
    s[j-1] = s[j-1] / conj ( r[j-1+(j-1)*ldr] );
  }

  norm = dznrm2 ( p, s, 1 );

  if ( 1.0 <= norm )
  {
    info = -1;
    return info;
  }

  alpha = sqrt ( 1.0 - norm * norm );
//
//  Determine the transformations.
//
  for ( ii = 1; ii <= p; ii++ )
  {
    i = p - ii + 1;
    scale = alpha + abs ( s[i-1] );
    a = alpha / scale;
    b = s[i-1] / scale;
    norm = sqrt ( a * a + real ( b ) * real ( b ) + imag ( b ) * imag ( b ) );
    c[i-1] = a / norm;
    s[i-1] = conj ( b ) / norm;
    alpha = scale * norm;
  }
//
//  Apply the transformations to R.
//
  for ( j = 1; j <= p; j++ )
  {
    xx = complex <double> ( 0.0, 0.0 );
    for ( ii = 1; ii <= j; ii++ )
    {
      i = j - ii + 1;
      t = c[i-1] * xx + s[i-1] * r[i-1+(j-1)*ldr];
      r[i-1+(j-1)*ldr] = c[i-1] * r[i-1+(j-1)*ldr] - conj ( s[i-1] ) * xx;
      xx = t;
    }
  }
//
//  If required, downdate Z and RHO.
//
  for ( j = 1; j <= nz; j++ )
  {
    zeta = y[j-1];

    for ( i = 1; i <= p; i++ )
    {
      z[i-1+(j-1)*ldz] = ( z[i-1+(j-1)*ldz]
        - conj ( s[i-1] ) * zeta ) / c[i-1];
      zeta = c[i-1] * zeta - s[i-1] * z[i-1+(j-1)*ldz];
    }

    azeta = abs ( zeta );

    if ( rho[j-1] < azeta )
    {
      info = 1;
      rho[j-1] = -1.0;
    }
    else
    {
      rho[j-1] = rho[j-1]
        * sqrt ( 1.0 - ( azeta / rho[j-1] ) * ( azeta / rho[j-1] ) );
    }
  }

  return info;
}
//****************************************************************************80

void zchex ( complex <double> r[], int ldr, int p, int k, int l,
  complex <double> z[], int ldz, int nz, double c[], complex <double> s[],
  int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZCHEX updates a Cholesky factorization.
//
//  Discussion:
//
//    ZCHEX updates a Cholesky factorization
//
//      A = hermitian(R) * R
//
//    of a positive definite matrix A of order P under diagonal
//    permutations of the form
//
//      E' * A * E
//
//    where E is a permutation matrix.  Specifically, given
//    an upper triangular matrix R and a permutation matrix
//    E (which is specified by K, L, and JOB), ZCHEX determines
//    a unitary matrix U such that
//
//      U * R * E = RR,
//
//    where RR is upper triangular.  At the user's option, the
//    transformation U will be multiplied into the array Z.
//
//    If A = hermitian(X)*X, so that R is the triangular part of the
//    QR factorization of X, then RR is the triangular part of the
//    QR factorization of X * E, that is, X with its columns permuted.
//
//    For a less terse description of what ZCHEX does and how
//    it may be applied, see the LINPACK guide.
//
//    The matrix Q is determined as the product U(L-K)*...*U(1)
//    of plane rotations of the form
//
//      (    C(I)       S(I) )
//      (                    ) ,
//      ( -conj(S(i))  C(I) )
//
//    where C(I) is real, the rows these rotations operate on
//    are described below.
//
//    There are two types of permutations, which are determined
//    by the value of job.
//
//    JOB = 1, right circular shift:
//    The columns are rearranged in the following order.
//
//      1, ..., K-1, L, K, K+1, ..., L-1, L+1, ..., P.
//
//    U is the product of L-K rotations U(I), where U(I)
//    acts in the (L-I,L-I+1)-plane.
//
//    JOB = 2, left circular shift:
//    The columns are rearranged in the following order
//
//      1, ..., K-1, K+1, K+2, ..., L, L, L+1, ..., P.
//
//    U is the product of L-K rotations U(I), where U(I)
//    acts in the (K+I-1,K+I)-plane.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> R[LDR*P]; On input, the upper triangular factor
//    that is to be updated.  On output, the updated factor.  Elements
//    below the diagonal are not referenced.
//
//    Input, int LDR, the leading dimension of R, which is at least P.
//
//    Input, int P, the order of the matrix.
//
//    Input, int K, the first column to be permuted.
//
//    Input, int L, the last column to be permuted.
//    L must be strictly greater than K.
//
//    Input/output, complex <double> Z[LDZ*NZ]; on input, an array of NZ P-vectors into
//    which the transformation U is multiplied.  On output, the updated
//    matrix.  Z is not referenced if NZ = 0.
//
//    Input, int LDZ, the leading dimension of Z, which must
//    be at least P.
//
//    Input, int NZ, the number of columns of the matrix Z.
//
//    Output, double C[P], the cosines of the transforming rotations.
//
//    Output, complex <double> S[P], the sines of the transforming rotations.
//
//    Input, int JOB, determines the type of permutation.
//    1, right circular shift.
//    2, left circular shift.
//
{
  int i;
  int ii;
  int il;
  int iu;
  int j;
  int jj;
  complex <double> t;

  if ( job == 1 )
  {
//
//  Right circular shift.
//
//  Reorder the columns.
//
    for ( i = 1; i <= l; i++ )
    {
      ii = l - i + 1;
      s[i-1] = r[ii-1+(l-1)*ldr];
    }

    for ( jj = k; jj <= l - 1; jj++ )
    {
      j = l - 1 - jj + k;
      for ( i = 1; i <= j; i++ )
      {
        r[i-1+j*ldr] = r[i-1+(j-1)*ldr];
      }
      r[j+j*ldr] = complex <double> ( 0.0, 0.0 );
    }
    for ( i = 1; i <= k-1; i++ )
    {
      ii = l - i + 1;
      r[i-1+(k-1)*ldr] = s[ii-1];
    }
//
//  Calculate the rotations.
//
    t = s[0];
    for ( i = 1; i <= l - k; i++ )
    {
      zrotg ( s+i, t, c+i-1, s+i-1 );
      t = s[i];
    }

    r[k-1+(k-1)*ldr] = t;
    for ( j = k+1; j <= p; j++ )
    {
      il = i4_max ( 1, l-j+1 );
      for ( ii = il; ii <= l - k; ii++ )
      {
        i = l - ii;
        t = c[ii-1] * r[i-1+(j-1)*ldr] + s[ii-1] * r[i+(j-1)*ldr];
        r[i+(j-1)*ldr] = c[ii-1] * r[i+(j-1)*ldr]
          - conj ( s[ii-1] ) * r[i-1+(j-1)*ldr];
        r[i-1+(j-1)*ldr] = t;
      }
    }
//
//  If required, apply the transformations to Z.
//
    for ( j = 1; j <= nz; j++ )
    {
      for ( ii = 1; ii <= l - k; ii++ )
      {
        i = l - ii;
        t = c[ii-1] * z[i-1+(j-1)*ldz] + s[ii-1] * z[i+(j-1)*ldz];
        z[i+(j-1)*ldz] = c[ii-1] * z[i+(j-1)*ldz]
          - conj ( s[ii-1] ) * z[i-1+(j-1)*ldz];
        z[i-1+(j-1)*ldz] = t;
      }
    }
  }
  else
  {
//
//  Left circular shift.
//
//  Reorder the columns.
//
    for ( i = 1; i <= k; i++ )
    {
      ii = l - k + i;
      s[ii-1] = r[i-1+(k-1)*ldr];
    }

    for ( j = k; j <= l - 1; j++ )
    {
      for ( i = 1; i <= j; i++ )
      {
        r[i-1+(j-1)*ldr] = r[i-1+j*ldr];
      }
      jj = j - k + 1;
      s[jj-1] = r[j+j*ldr];
    }

    for ( i = 1; i <= k; i++ )
    {
      ii = l - k + i;
      r[i-1+(l-1)*ldr] = s[ii-1];
    }

    for ( i = k + 1; i <= l; i++ )
    {
      r[i-1+(l-1)*ldr] = complex <double> ( 0.0, 0.0 );
    }
//
//  Reduction loop.
//
    for ( j = k; j <= p; j++ )
    {
//
//  Apply the rotations.
//
      if ( j != k )
      {
        iu = i4_min ( j - 1, l - 1 );
        for ( i = k; i <= iu; i++ )
        {
          ii = i - k + 1;
          t = c[ii-1] * r[i-1+(j-1)*ldr] + s[ii-1] * r[i+(j-1)*ldr];
          r[i+(j-1)*ldr] = c[ii-1] * r[i+(j-1)*ldr]
            - conj ( s[ii-1] ) * r[i-1+(j-1)*ldr];
          r[i-1+(j-1)*ldr] = t;
        }
      }

      if ( j < l )
      {
        jj = j - k + 1;
        t = s[jj-1];
        zrotg ( r+j-1+(j-1)*ldr, t, c+jj-1, s+jj-1 );
      }
    }
//
//  Apply the rotations to Z.
//
    for ( j = 1; j <= nz; j++ )
    {
      for ( i = k; i <= l - 1; i++ )
      {
        ii = i - k + 1;
        t = c[ii-1] * z[i-1+(j-1)*ldz] + s[ii-1] * z[i+(j-1)*ldz];
        z[i+(j-1)*ldz] = c[ii-1] * z[i+(j-1)*ldz] - conj ( s[ii-1] )
          * z[i-1+(j-1)*ldz];
        z[i-1+(j-1)*ldz] = t;
      }
    }
  }
  return;
}
//****************************************************************************80

void zchud ( complex <double> r[], int ldr, int p, complex <double> x[],
  complex <double> z[], int ldz, int nz, complex <double> y[], double rho[],
  double c[], complex <double> s[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZCHUD updates an augmented Cholesky decomposition.
//
//  Discussion:
//
//    ZCHUD updates an augmented Cholesky decomposition of the
//    triangular part of an augmented QR decomposition.  Specifically,
//    given an upper triangular matrix R of order P, a row vector
//    X, a column vector Z, and a scalar Y, ZCHUD determines a
//    unitary matrix U and a scalar ZETA such that
//
//           ( R  Z )     ( RR   ZZ  )
//      U  * (      )  =  (          ),
//           ( X  Y )     (  0  ZETA )
//
//    where RR is upper triangular.  If R and Z have been
//    obtained from the factorization of a least squares
//    problem, then RR and ZZ are the factors corresponding to
//    the problem with the observation (X,Y) appended.  In this
//    case, if RHO is the norm of the residual vector, then the
//    norm of the residual vector of the updated problem is
//    sqrt ( RHO**2 + ZETA**2 ).  ZCHUD will simultaneously update
//    several triplets (Z,Y,RHO).
//
//    For a less terse description of what ZCHUD does and how
//    it may be applied see the LINPACK guide.
//
//    The matrix U is determined as the product U(P)*...*U(1),
//    where U(I) is a rotation in the (I,P+1) plane of the
//    form
//
//      (          C(I)    S(I) )
//      (                       ).
//      ( -conjg ( S(I) )  C(I) )
//
//    The rotations are chosen so that C(I) is real.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> R[LDR*P], the upper triangular matrix
//    that is to be updated.  The part of R below the diagonal is
//    not referenced.
//
//    Input, int LDR, the leading dimension of R.
//    P <= LDR.
//
//    Input, int P, the order of the matrix.
//
//    Input, complex <double> X[P], the row to be added to R.
//
//    Input/output, complex <double> Z[LDZ*NZ], NZ P-vectors to
//    be updated with R.
//
//    Input, int LDZ, the leading dimension of Z.
//    P <= LDZ.
//
//    Input, int NZ, the number of vectors to be updated.
//    NZ may be zero, in which case Z, Y, and RHO are not referenced.
//
//    Input, complex <double> Y[NZ], the scalars for updating the vectors Z.
//
//    Input/output, double RHO[NZ]; on input, the norms of the residual
//    vectors that are to be updated.  If RHO(J) is negative, it is
//    left unaltered.  On output, the updated values.
//
//    Output, double C[P]. the cosines of the transforming rotations.
//
//    Output, complex <double> S[P], the sines of the transforming rotations.
//
{
  double azeta;
  int i;
  int j;
  double scale;
  complex <double> t;
  complex <double> xj;
  complex <double> zeta;
//
//  Update R.
//
  for ( j = 1; j <= p; j++ )
  {
    xj = x[j-1];
//
//  Apply the previous rotations.
//
    for ( i = 1; i <= j - 1; i++ )
    {
      t = c[i-1] * r[i-1+(j-1)*ldr] + s[i-1] * xj;
      xj = c[i-1] * xj - conj ( s[i-1] ) * r[i-1+(j-1)*ldr];
      r[i-1+(j-1)*ldr] = t;
    }
//
//  Compute the next rotation.
//
    zrotg ( r+j-1+(j-1)*ldr, xj, c+j-1, s+j-1 );
  }
//
//  If required, update Z and RHO.
//
  for ( j = 1; j <= nz; j++ )
  {
    zeta = y[j-1];

    for ( i = 1; i <= p; i++ )
    {
      t = c[i-1] * z[i-1+(j-1)*ldz] + s[i-1] * zeta;
      zeta = c[i-1] * zeta - conj ( s[i-1] ) * z[i-1+(j-1)*ldz];
      z[i-1+(j-1)*ldz] = t;
    }

    azeta = abs ( zeta );

    if ( azeta != 0.0 && 0.0 <= rho[j-1] )
    {
      scale = azeta + rho[j-1];
      rho[j-1] = scale * sqrt ( pow ( azeta / scale, 2 )
                            + pow ( rho[j-1] / scale, 2 ) );
    }
  }
  return;
}
//****************************************************************************80

double zgbco ( complex <double> abd[], int lda, int n, int ml, int mu,
  int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZGBCO factors a complex band matrix and estimates its condition.
//
//  Discussion:
//
//    If RCOND is not needed, ZGBFA is slightly faster.
//
//    To solve A*X = B, follow ZGBCO by ZGBSL.
//
//    To compute inverse(A)*C, follow ZGBCO by ZGBSL.
//
//    To compute determinant(A), follow ZGBCO by ZGBDI.
//
//  Band storage:
//
//    If A is a band matrix, the following program segment
//    will set up the input.
//
//      ml = (band width below the diagonal)
//      mu = (band width above the diagonal)
//      m = ml + mu + 1
//      do j = 1, n
//        i1 = max ( 1, j - mu )
//        i2 = min ( n, j + ml )
//        do i = i1, i2
//          k = i - j + m
//          abd(k,j) = a(i,j)
//        }
//      }
//
//    This uses rows ML+1 through 2*ML+MU+1 of ABD.
//    In addition, the first ML rows in ABD are used for
//    elements generated during the triangularization.
//    The total number of rows needed in ABD is 2*ML+MU+1.
//    The ML+MU by ML+MU upper left triangle and the
//    ML by ML lower right triangle are not referenced.
//
//  Example:
//
//    If the original matrix A is
//
//      11 12 13  0  0  0
//      21 22 23 24  0  0
//       0 32 33 34 35  0
//       0  0 43 44 45 46
//       0  0  0 54 55 56
//       0  0  0  0 65 66
//
//     Then N = 6, ML = 1, MU = 2, 5 <= LDA and ABD should contain
//
//       *  *  *  +  +  +
//       *  * 13 24 35 46
//       * 12 23 34 45 56
//      11 22 33 44 55 66
//      21 32 43 54 65  *
//
//    * = not used,
//    + = used for pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> ABD[LDA*N], on input, contains the matrix in
//    band storage.  The columns of the matrix are stored in the columns
//    of ABD and the diagonals of the matrix are stored in rows ML+1
//    through 2*ML+MU+1 of ABD.  On output, an upper triangular matrix
//    in band storage and the multipliers which were used to obtain it.
//    The factorization can be written A = L*U where L is a product of
//    permutation and unit lower triangular matrices and U is upper triangular.
//
//    Input, int LDA, the leading dimension of ABD.
//    LDA must be at least 2*ML+MU+1.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the number of diagonals below the main diagonal.
//    0 <= ML < N.
//
//    Input, int MU, the number of diagonals above the main diagonal.
//    0 <= MU < N.
//    More efficient if ML <= MU.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, double ZGBCO, an estimate of the reciprocal condition RCOND of A.
//    For the system A*X = B, relative perturbations in A and B of size
//    epsilon may cause relative perturbations in X of size (EPSILON/RCOND).
//    If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate
//    underflows.
//
//  Local Parameters:
//
//    Workspace, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
//
{
  double anorm;
  complex <double> ek;
  int info;
  int is;
  int j;
  int ju;
  int k;
  int l;
  int la;
  int lm;
  int lz;
  int m;
  int mm;
  double rcond;
  double s;
  double sm;
  complex <double> t;
  complex <double> wk;
  complex <double> wkm;
  double ynorm;
  complex <double> *z;

  z = new complex <double> [n];
//
//  Compute 1-norm of A.
//
  anorm = 0.0;
  l = ml + 1;
  is = l + mu;

  for ( j = 1; j <= n; j++ )
  {
    anorm = r8_max ( anorm, dzasum ( l, abd+is-1+(j-1)*lda, 1 ) );

    if ( ml + 1 < is )
    {
      is = is - 1;
    }

    if ( j <= mu )
    {
      l = l + 1;
    }

    if ( n - ml <= j )
    {
      l = l - 1;
    }
  }
//
//  Factor
//
  info = zgbfa ( abd, lda, n, ml, mu, ipvt );
//
//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and hermitian(A)*Y = E.
//
//  Hermitian(A) is the conjugate transpose of A.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where hermitian(U)*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve hermitian(U) * W = E.
//
  ek = complex <double> ( 1.0, 0.0 );

  for ( j = 1; j <= n; j++ )
  {
    z[j-1] = complex <double> ( 0.0, 0.0 );
  }

  m = ml + mu + 1;
  ju = 0;

  for ( k = 1; k <= n; k++ )
  {
    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, -z[k-1] );
    }

    if ( zabs1 ( abd[m-1+(k-1)*lda] ) < zabs1 ( ek - z[k-1] ) )
    {
      s = zabs1 ( abd[m-1+(k-1)*lda] ) / zabs1 ( ek - z[k-1] );
      zdscal ( n, s, z, 1 );
      ek = complex <double> ( s, 0.0 ) * ek;
    }

    wk = ek - z[k-1];
    wkm = -ek - z[k-1];
    s = zabs1 ( wk );
    sm = zabs1 ( wkm );

    if ( zabs1 ( abd[m-1+(k-1)*lda] ) != 0.0 )
    {
      wk = wk / conj ( abd[m-1+(k-1)*lda] );
      wkm = wkm / conj ( abd[m-1+(k-1)*lda] );
    }
    else
    {
      wk = complex <double> ( 1.0, 0.0 );
      wkm = complex <double> ( 1.0, 0.0 );
    }

    ju = i4_min ( i4_max ( ju, mu + ipvt[k-1] ), n );
    mm = m;

    if ( k+1 <= ju )
    {
      for ( j = k+1; j <= ju; j++ )
      {
        mm = mm - 1;
        sm = sm + zabs1 ( z[j-1] + wkm * conj ( abd[mm-1+(j-1)*lda] ) );
        z[j-1] = z[j-1] + wk * conj ( abd[mm-1+(j-1)*lda] );
        s = s + zabs1 ( z[j-1] );
      }

      if ( s < sm )
      {
        t = wkm - wk;
        wk = wkm;
        mm = m;
        for ( j = k+1; j <= ju; j++ )
        {
          mm = mm - 1;
          z[j-1] = z[j-1] + t * conj ( abd[mm-1+(j-1)*lda] );
        }
      }
    }
    z[k-1] = wk;
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
//
//  Solve hermitian(L) * Y = W.
//
  for ( k = n; 1 <= k; k-- )
  {
    lm = i4_min ( ml, n - k );

    if ( k < n )
    {
      z[k-1] = z[k-1] + zdotc ( lm, abd+m+(k-1)*lda, 1, z+k, 1 );
    }

    if ( 1.0 < zabs1 ( z[k-1] ) )
    {
      s = 1.0 / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
    }

    l = ipvt[k-1];

    t      = z[l-1];
    z[l-1] = z[k-1];
    z[k-1] = t;
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = 1.0;
//
//  Solve L * V = Y.
//
  for ( k = 1; k <= n; k++ )
  {
    l = ipvt[k-1];

    t      = z[l-1];
    z[l-1] = z[k-1];
    z[k-1] = t;

    lm = i4_min ( ml, n - k );

    if ( k < n )
    {
      zaxpy ( lm, t, abd+m+(k-1)*lda, 1, z+k, 1 );
    }

    if ( 1.0 < zabs1 ( z[k-1] ) )
    {
      s = 1.0 / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;
//
//  Solve U * Z = W.
//
  for ( k = n; 1 <= k; k-- )
  {
    if ( zabs1 ( abd[m-1+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
    {
      s = zabs1 ( abd[m-1+(k-1)*lda] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }

    if ( zabs1 ( abd[m-1+(k-1)*lda] ) != 0.0 )
    {
      z[k-1] = z[k-1] / abd[m-1+(k-1)*lda];
    }
    else
    {
      z[k-1] = complex <double> ( 1.0, 0.0 );
    }

    lm = i4_min ( k, m ) - 1;
    la = m - lm;
    lz = k - lm;
    t = -z[k-1];
    zaxpy ( lm, t, abd+la-1+(k-1)*lda, 1, z+lz-1, 1 );
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }

  delete [] z;

  return rcond;
}
//****************************************************************************80

void zgbdi ( complex <double> abd[], int lda, int n, int ml, int mu, int ipvt[],
  complex <double> det[2] )

//****************************************************************************80
//
//  Purpose:
//
//    ZGBDI computes the determinant of a band matrix factored by ZGBCO or ZGBFA.
//
//  Discussion:
//
//    If the inverse is needed, use ZGBSL N times.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> ABD[LDA*N], the output from ZGBCO or ZGBFA.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the number of diagonals below the main diagonal.
//
//    Input, int MU, the number of diagonals above the main diagonal.
//
//    Input, int IPVT[N], the pivot vector from ZGBCO or ZGBFA.
//
//    Output, complex <double> DET[2], determinant of original matrix.
//    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= zabs1 ( DET(1) ) < 10.0
//    or DET(1) = 0.0.  Also, DET(2) is strictly real.
//
{
  int i;
  int m;

  m = ml + mu + 1;
  det[0] = complex <double> ( 1.0, 0.0 );
  det[1] = complex <double> ( 0.0, 0.0 );

  for ( i = 1; i <= n; i++ )
  {
    if ( ipvt[i-1] != i )
    {
      det[0] = -det[0];
    }

    det[0] = det[0] * abd[m-1+(i-1)*lda];

    if ( zabs1 ( det[0] ) == 0.0 )
    {
      break;
    }

    while ( zabs1 ( det[0] ) < 1.0 )
    {
      det[0] = det[0] * complex <double> ( 10.0, 0.0 );
      det[1] = det[1] - complex <double> ( 1.0, 0.0 );
    }

    while ( 10.0 <= zabs1 ( det[0] ) )
    {
      det[0] = det[0] / complex <double> ( 10.0, 0.0 );
      det[1] = det[1] + complex <double> ( 1.0, 0.0 );
    }

  }

  return;
}
//****************************************************************************80

int zgbfa ( complex <double> abd[], int lda, int n, int ml, int mu, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZGBFA factors a complex band matrix by elimination.
//
//  Discussion:
//
//    ZGBFA is usually called by ZGBCO, but it can be called
//    directly with a saving in time if RCOND is not needed.
//
//  Band storage:
//
//    If A is a band matrix, the following program segment
//    will set up the input.
//
//      ml = (band width below the diagonal)
//      mu = (band width above the diagonal)
//      m = ml + mu + 1
//      do j = 1, n
//        i1 = max ( 1, j - mu )
//        i2 = min ( n, j + ml )
//        do i = i1, i2
//          k = i - j + m
//          abd(k,j) = a(i,j)
//        end do
//      end do
//
//    This uses rows ML+1 through 2*ML+MU+1 of ABD.
//    In addition, the first ML rows in ABD are used for
//    elements generated during the triangularization.
//    The total number of rows needed in ABD is 2*ML+MU+1.
//    The ML+MU by ML+MU upper left triangle and the
//    ML by ML lower right triangle are not referenced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> ABD[LDA*N], on input, contains the matrix in
//    band storage.  The columns of the matrix are stored in the columns
//    of ABD and the diagonals of the matrix are stored in rows ML+1
//    through 2*ML+MU+1 of ABD.  On output, an upper triangular matrix
//    in band storage and the multipliers which were used to obtain it.
//    The factorization can be written A = L*U where L is a product of
//    permutation and unit lower triangular matrices and U is upper triangular.
//
//    Input, int LDA, the leading dimension of ABD.
//    LDA must be at least 2*ML+MU+1.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the number of diagonals below the main diagonal.
//    0 <= ML < N.
//
//    Input, int MU, the number of diagonals above the main diagonal.
//    0 <= MU < N.  More efficient if ML <= MU.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int ZGBFA.
//    0, normal value.
//    K, if U(K,K) == 0.0.  This is not an error condition for this
//    subroutine, but it does indicate that ZGBSL will divide by zero if
//    called.  Use RCOND in ZGBCO for a reliable indication of singularity.
//
{
  int i;
  int i0;
  int info;
  int j;
  int j0;
  int j1;
  int ju;
  int jz;
  int k;
  int l;
  int lm;
  int m;
  int mm;
  complex <double> t;

  m = ml + mu + 1;
  info = 0;
//
//  Zero initial fill-in columns.
//
  j0 = mu + 2;
  j1 = i4_min ( n, m ) - 1;

  for ( jz = j0; jz <= j1; jz++ )
  {
    i0 = m + 1 - jz;
    for ( i = i0; i <= ml; i++ )
    {
      abd[i-1+(jz-1)*lda] = complex <double> ( 0.0, 0.0 );
    }
  }

  jz = j1;
  ju = 0;
//
//  Gaussian elimination with partial pivoting.
//
  for ( k = 1; k <= n-1; k++ )
  {
//
//  Zero next fill-in column
//
    jz = jz + 1;
    if ( jz <= n )
    {
      for ( i = 1; i <= ml; i++ )
      {
        abd[i-1+(jz-1)*lda] = complex <double> ( 0.0, 0.0 );
      }
    }
//
//  Find L = pivot index.
//
    lm = i4_min ( ml, n - k );
    l = izamax ( lm+1, abd+m-1+(k-1)*lda, 1 ) + m - 1;
    ipvt[k-1] = l + k - m;
//
//  Zero pivot implies this column already triangularized.
//
    if ( zabs1 ( abd[l-1+(k-1)*lda] ) == 0.0 )
    {
      info = k;
      continue;
    }
//
//  Interchange if necessary.
//
    if ( l != m )
    {
      t                  = abd[l-1+(k-1)*lda];
      abd[l-1+(k-1)*lda] = abd[m-1+(k-1)*lda];
      abd[m-1+(k-1)*lda] = t;
    }
//
//  Compute multipliers.
//
    t = - complex <double> ( 1.0, 0.0 ) / abd[m-1+(k-1)*lda];
    zscal ( lm, t, abd+m+(k-1)*lda, 1 );
//
//  Row elimination with column indexing.
//
    ju = i4_min ( i4_max ( ju, mu + ipvt[k-1] ), n );
    mm = m;

    for ( j = k+1; j <= ju; j++ )
    {
      l = l - 1;
      mm = mm - 1;
      t = abd[l-1+(j-1)*lda];
      if ( l != mm )
      {
        abd[l-1+(j-1)*lda] = abd[mm-1+(j-1)*lda];
        abd[mm-1+(j-1)*lda] = t;
      }
      zaxpy ( lm, t, abd+m+(k-1)*lda, 1, abd+mm+(j-1)*lda, 1 );
    }
  }

  ipvt[n-1] = n;

  if ( zabs1 ( abd[m-1+(n-1)*lda] ) == 0.0 )
  {
    info = n;
  }

  return info;
}
//****************************************************************************80

void zgbsl ( complex <double> abd[], int lda, int n, int ml, int mu,
  int ipvt[], complex <double> b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZGBSL solves a complex band system factored by ZGBCO or ZGBFA.
//
//  Discussion:
//
//    ZGBSL can solve A * X = B or hermitan ( A ) * X = B.
//
//    A division by zero will occur if the input factor contains a
//    zero on the diagonal.  Technically this indicates singularity
//    but it is often caused by improper arguments or improper
//    setting of LDA.  It will not occur if the subroutines are
//    called correctly and if ZGBCO has set 0.0 < RCOND
//    or ZGBFA has set INFO = 0.
//
//    To compute inverse ( A ) * C where C is a matrix with P columns:
//
//      call zgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
//
//      if ( rcond is not too small ) then
//        do j = 1, p
//          call zgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
//        end do
//      end if
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> ABD[LDA*N], the output from ZGBCO or ZGBFA.
//
//    Input, int LDA, the leading dimension of ABD.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, the number of diagonals below the main diagonal.
//
//    Input, int MU, the number of diagonals above the main diagonal.
//
//    Input, int IPVT[N], the pivot vector from ZGBCO or ZGBFA.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
//    Input, int JOB.
//    0, to solve A*x = b,
//    nonzero, to solve hermitian(A)*x = b, where hermitian(A) is the
//    conjugate transpose.
//
{
  int k;
  int l;
  int la;
  int lb;
  int lm;
  int m;
  complex <double> t;

  m = mu + ml + 1;
//
//  JOB = 0, solve A * X = B.
//
  if ( job == 0 )
  {
//
//  First solve L * Y = B.
//
    if ( ml != 0 )
    {
      for ( k = 1; k <= n-1; k++ )
      {
        lm = i4_min ( ml, n - k );
        l = ipvt[k-1];
        t = b[l-1];

        if ( l != k )
        {
          b[l-1] = b[k-1];
          b[k-1] = t;
        }
        zaxpy ( lm, t, abd+m+(k-1)*lda, 1, b+k, 1 );
      }
    }
//
//  Now solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] / abd[m-1+(k-1)*lda];
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      t = -b[k-1];
      zaxpy ( lm, t, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
    }
  }
//
//  JOB = nonzero, solve hermitian(A) * X = B.
//
  else
  {
//
//  First solve hermitian ( U ) * Y = B.
//
    for ( k = 1; k <= n; k++ )
    {
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      t = zdotc ( lm, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
      b[k-1] = ( b[k-1] - t ) / conj ( abd[m-1+(k-1)*lda] );
    }
//
//  Now solve hermitian ( L ) * X = Y.
//
    if ( ml != 0 )
    {
      for ( k = n-1; 1 <= k; k-- )
      {
        lm = i4_min ( ml, n - k );
        b[k-1] = b[k-1] + zdotc ( lm, abd+m+(k-1)*lda, 1, b+k, 1 );
        l = ipvt[k-1];

        if ( l != k )
        {
          t      = b[l-1];
          b[l-1] = b[k-1];
          b[k-1] = t;
        }
      }
    }
  }
  return;
}
//****************************************************************************80

double zgeco ( complex <double> a[], int lda, int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZGECO factors a complex matrix and estimates its condition.
//
//  Discussion:
//
//    If RCOND is not needed, ZGEFA is slightly faster.
//
//    To solve A*X = B, follow ZGECO by ZGESL.
//
//    To compute inverse(A)*C, follow ZGECO by ZGESL.
//
//    To compute determinant(A), follow ZGECO by ZGEDI.
//
//    To compute inverse(A), follow ZGECO by ZGEDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N], on input, the matrix to be
//    factored.  On output, an upper triangular matrix and the multipliers
//    used to obtain it.  The factorization can be written A = L*U where
//    L is a product of permutation and unit lower triangular matrices
//    and U is upper triangular.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, double SGECO, an estimate of the reciprocal condition of A.
//    For the system A*X = B, relative perturbations in A and B of size
//    EPSILON may cause relative perturbations in X of size (EPSILON/RCOND).
//    If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate
//    underflows.
//
//  Local Parameters:
//
//    Workspace, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is
//    an approximate null vector in the sense that
//      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
//
{
  double anorm;
  complex <double> ek;
  int i;
  int info;
  int j;
  int k;
  int l;
  double rcond;
  double s;
  double sm;
  complex <double> t;
  complex <double> wk;
  complex <double> wkm;
  double ynorm;
  complex <double> *z;

  z = new complex <double> [n];
//
//  Compute the 1-norm of A.
//
  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    anorm = r8_max ( anorm, dzasum ( n, a+0+j*lda, 1 ) );
  }
//
//  Factor.
//
  info = zgefa ( a, lda, n, ipvt );
//
//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and hermitian(A)*Y = E.
//
//  Hermitian(A) is the conjugate transpose of A.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where hermitian(U)*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve hermitian(U)*W = E.
//
  ek = complex <double> ( 1.0, 0.0 );
  for ( i = 0; i < n; i++ )
  {
    z[i] = complex <double> ( 0.0, 0.0 );
  }

  for ( k = 1; k <= n; k++ )
  {
    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, -z[k-1] );
    }

    if ( zabs1 ( a[k-1+(k-1)*lda] ) < zabs1 ( ek - z[k-1] ) )
    {
      s = zabs1 ( a[k-1+(k-1)*lda] ) / zabs1 ( ek - z[k-1] );
      zdscal ( n, s, z, 1 );
      ek = complex <double> ( s, 0.0 ) * ek;
    }

    wk = ek - z[k-1];
    wkm = -ek - z[k-1];
    s = zabs1 ( wk );
    sm = zabs1 ( wkm );

    if ( zabs1 ( a[k-1+(k-1)*lda] ) != 0.0 )
    {
      wk = wk / conj ( a[k-1+(k-1)*lda] );
      wkm = wkm / conj ( a[k-1+(k-1)*lda] );
    }
    else
    {
      wk = complex <double> ( 1.0, 0.0 );
      wkm = complex <double> ( 1.0, 0.0 );
    }

    for ( j = k+1; j <= n; j++ )
    {
      sm = sm + zabs1 ( z[j-1] + wkm * conj ( a[k-1+(j-1)*lda] ) );
      z[j-1] = z[j-1] + wk * conj ( a[k-1+(j-1)*lda] );
      s = s + zabs1 ( z[j-1] );
    }

    if ( s < sm )
    {
      t = wkm - wk;
      wk = wkm;
      for ( j = k+1; j <= n; j++ )
      {
        z[j-1] = z[j-1] + t * conj ( a[k-1+(j-1)*lda] );
      }
    }
    z[k-1] = wk;
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
//
//  Solve hermitian(L) * Y = W.
//
  for ( k = n; 1 <= k; k-- )
  {
    if ( k < n )
    {
      z[k-1] = z[k-1] + zdotc ( n-k, a+k+(k-1)*lda, 1, z+k, 1 );
    }

    if ( 1.0 < zabs1 ( z[k-1] ) )
    {
      s = 1.0 / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
    }
    l = ipvt[k-1];

    t      = z[l-1];
    z[l-1] = z[k-1];
    z[k-1] = t;
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );

  ynorm = 1.0;
//
//  Solve L * V = Y.
//
  for ( k = 1; k <= n; k++ )
  {
    l = ipvt[k-1];

    t      = z[l-1];
    z[l-1] = z[k-1];
    z[k-1] = t;

    if ( k < n )
    {
      zaxpy ( n-k, t, a+k+(k-1)*lda, 1, z+k, 1 );
    }

    if ( 1.0 < zabs1 ( z[k-1] ) )
    {
      s = 1.0 / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;
//
//  Solve U * Z = V.
//
  for ( k = n; 1 <= k; k-- )
  {
    if ( zabs1 ( a[k-1+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
    {
      s = zabs1 ( a[k-1+(k-1)*lda] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }

    if ( zabs1 ( a[k-1+(k-1)*lda] ) != 0.0 )
    {
      z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
    }
    else
    {
      z[k-1] = complex <double> ( 1.0, 0.0 );
    }

    t = -z[k-1];
    zaxpy ( k-1, t, a+0+(k-1)*lda, 1, z, 1 );
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }

  delete [] z;

  return rcond;
}
//****************************************************************************80

void zgedi ( complex <double> a[], int lda, int n, int ipvt[],
  complex <double> det[2], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZGEDI computes the determinant and inverse of a matrix.
//
//  Discussion:
//
//    The matrix must have been factored by ZGECO or ZGEFA.
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal and the inverse is requested.
//    It will not occur if the subroutines are called correctly
//    and if ZGECO has set 0.0 < RCOND or ZGEFA has set
//    INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the factor information
//    from ZGECO or ZGEFA.  On output, the inverse matrix, if it
//    was requested,
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZGECO or ZGEFA.
//
//    Output, complex <double> DET[2], the determinant of the original matrix,
//    if requested.  Otherwise not referenced.
//    Determinant = DET(1) * 10.0**DET(2) with
//    1.0 <= zabs1 ( DET(1) ) < 10.0 or DET(1) == 0.0.
//    Also, DET(2) is strictly real.
//
//    Input, int JOB.
//    11, both determinant and inverse.
//    01, inverse only.
//    10, determinant only.
//
{
  int i;
  int j;
  int k;
  int l;
  complex <double> t;
  complex <double> *work;
//
//  Compute the determinant.
//
  if ( job / 10 != 0 )
  {
    det[0] = complex <double> ( 1.0, 0.0 );
    det[1] = complex <double> ( 0.0, 0.0 );

    for ( i = 1; i <= n; i++ )
    {
      if ( ipvt[i-1] != i )
      {
        det[0] = -det[0];
      }

      det[0] = a[i-1+(i-1)*lda] * det[0];

      if ( zabs1 ( det[0] ) == 0.0 )
      {
        break;
      }

      while ( zabs1 ( det[0] ) < 1.0 )
      {
        det[0] = det[0] * complex <double> ( 10.0, 0.0 );
        det[1] = det[1] - complex <double> ( 1.0, 0.0 );
      }

      while ( 10.0 <= zabs1 ( det[0] ) )
      {
        det[0] = det[0] / complex <double> ( 10.0, 0.0 );
        det[1] = det[1] + complex <double> ( 1.0, 0.0 );
      }
    }
  }
//
//  Compute inverse(U).
//
  if ( ( job % 10 ) != 0 )
  {
    work = new complex <double>[n];

    for ( k = 1; k <= n; k++ )
    {
      a[k-1+(k-1)*lda] = complex <double> ( 1.0, 0.0 ) / a[k-1+(k-1)*lda];
      t = -a[k-1+(k-1)*lda];
      zscal ( k-1, t, a+0+(k-1)*lda, 1 );

      for ( j = k+1; j <= n; j++ )
      {
        t = a[k-1+(j-1)*lda];
        a[k-1+(j-1)*lda] = complex <double> ( 0.0, 0.0 );
        zaxpy ( k, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
      }
    }
//
//  Form inverse(U) * inverse(L).
//
    for ( k = n-1; 1 <= k; k-- )
    {
      for ( i = k+1; i <= n; i++ )
      {
        work[i-1] = a[i-1+(k-1)*lda];
        a[i-1+(k-1)*lda] = complex <double> ( 0.0, 0.0 );
      }

      for ( j = k+1; j <= n; j++ )
      {
        t = work[j-1];
        zaxpy ( n, t, a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
      }

      l = ipvt[k-1];

      if ( l != k )
      {
        zswap( n, a+0+(k-1)*lda, 1, a+0+(l-1)*lda, 1 );
      }
    }

    delete [] work;
  }

  return;
}
//****************************************************************************80

int zgefa ( complex <double> a[], int lda, int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZGEFA factors a complex matrix by Gaussian elimination.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the matrix to be factored.
//    On output, an upper triangular matrix and the multipliers which were
//    used to obtain it.  The factorization can be written A = L*U where
//    L is a product of permutation and unit lower triangular matrices and
//    U is upper triangular.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int ZGEFA,
//    0, normal value.
//    K, if U(K,K) == 0.0.  This is not an error condition for this
//    subroutine, but it does indicate that ZGESL or ZGEDI will divide by zero
//    if called.  Use RCOND in ZGECO for a reliable indication of singularity.
//
{
  int info;
  int j;
  int k;
  int l;
  complex <double> t;
//
//  Gaussian elimination with partial pivoting.
//
  info = 0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Find L = pivot index.
//
    l = izamax ( n-k+1, a+(k-1)+(k-1)*lda, 1 ) + k - 1;
    ipvt[k-1] = l;
//
//  Zero pivot implies this column already triangularized.
//
    if ( zabs1 ( a[l-1+(k-1)*lda] ) == 0.0 )
    {
      info = k;
      continue;
    }
//
//  Interchange if necessary.
//
    if ( l != k )
    {
      t                = a[l-1+(k-1)*lda];
      a[l-1+(k-1)*lda] = a[k-1+(k-1)*lda];
      a[k-1+(k-1)*lda] = t;
    }
//
//  Compute multipliers
//
    t = - complex <double> ( 1.0, 0.0 ) / a[k-1+(k-1)*lda];
    zscal ( n-k, t, a+k+(k-1)*lda, 1 );
//
//  Row elimination with column indexing
//
    for ( j = k+1; j <= n; j++ )
    {
      t = a[l-1+(j-1)*lda];
      if ( l != k )
      {
        a[l-1+(j-1)*lda] = a[k-1+(j-1)*lda];
        a[k-1+(j-1)*lda] = t;
      }
      zaxpy ( n-k, t, a+k+(k-1)*lda, 1, a+k+(j-1)*lda, 1 );
    }

  }

  ipvt[n-1] = n;

  if ( zabs1 ( a[n-1+(n-1)*lda] ) == 0.0 )
  {
    info = n;
  }

  return info;
}
//****************************************************************************80

void zgesl ( complex <double> a[], int lda, int n, int ipvt[],
  complex <double> b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZGESL solves a complex system factored by ZGECO or ZGEFA.
//
//  Discussion:
//
//    A division by zero will occur if the input factor contains a
//    zero on the diagonal.  Technically this indicates singularity
//    but it is often caused by improper arguments or improper
//    setting of LDA.  It will not occur if the subroutines are
//    called correctly and if ZGECO has set 0.0 < RCOND
//    or ZGEFA has set INFO == 0.
//
//    To compute inverse(A) * C where C is a matrix with P columns:
//
//      call zgeco(a,lda,n,ipvt,rcond,z)
//
//      if (rcond is not too small) then
//        do j = 1, p
//          call zgesl ( a, lda, n, ipvt, c(1,j), 0 )
//        end do
//      end if
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> A[LDA*N], the factored matrix information,
//    as output from ZGECO or ZGEFA.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZGECO or ZGEFA.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
//    Input, int JOB.
//    0, to solve A*X = B,
//    nonzero, to solve hermitian(A)*X = B where hermitian(A) is the
//    conjugate transpose.
//
{
  int k;
  int l;
  complex <double> t;
//
//  JOB = 0, solve A * X = B.
//
//  First solve L * Y = B.
//
  if ( job == 0 )
  {
    for ( k = 1; k <= n-1; k++ )
    {
      l = ipvt[k-1];
      t = b[l-1];
      if ( l != k )
      {
        b[l-1] = b[k-1];
        b[k-1] = t;
      }
      zaxpy ( n-k, t, a+k+(k-1)*lda, 1, b+k, 1 );
    }
//
//  Now solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
      t = -b[k-1];
      zaxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
    }
  }
//
//  JOB nonzero, solve hermitian(A) * X = B.
//
//  First solve hermitian(U) * Y = B.
//
  else
  {
    for ( k = 1; k <= n; k++ )
    {
      t = zdotc ( k-1, a+0+(k-1)*lda, 1, b, 1 );
      b[k-1] = ( b[k-1] - t ) / conj ( a[k-1+(k-1)*lda] );
    }
//
//  Now solve hermitian(L) * X = Y.
//
    for ( k = n-1; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] + zdotc ( n-k, a+k+(k-1)*lda, 1, b+k, 1 );
      l = ipvt[k-1];
      if ( l != k )
      {
        t      = b[l-1];
        b[l-1] = b[k-1];
        b[k-1] = t;
      }
    }
  }

  return;
}
//****************************************************************************80

int zgtsl ( int n, complex <double> c[], complex <double> d[],
  complex <double> e[], complex <double> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZGTSL solves a complex general tridiagonal system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, complex <double> C[N]; on input, the subdiagonal
//    of the tridiagonal matrix in entries C(2:N).  On output, C has
//    been overwritten.
//
//    Input/output, complex <double> D[N]; on input, the diagonal of
//    the tridiagonal matrix.  On output, D has been overwritten.
//
//    Input/output, complex <double> E[N]; on input, the superdiagonal
//    of the tridiagonal matrix in entries E(1:N-1).  On output, E
//    has been overwritten.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
//    Output, int ZGTSL.
//    0, normal value.
//    K, if the K-th element of the diagonal becomes exactly zero.  The
//    subroutine returns when this is detected.
//
{
  int info;
  int k;
  complex <double> t;

  info = 0;
  c[0] = d[0];

  if ( 1 <= n-1 )
  {
    d[0] = e[0];
    e[0] = complex <double> ( 0.0, 0.0 );
    e[n-1] = complex <double> ( 0.0, 0.0 );

    for ( k = 1; k <= n-1; k++ )
    {
      if ( zabs1 ( c[k-1] ) <= zabs1 ( c[k] ) )
      {
        t      = c[k];
        c[k]   = c[k-1];
        c[k-1] = t;

        t      = d[k];
        d[k]   = d[k-1];
        d[k-1] = t;

        t      = e[k];
        e[k]   = e[k-1];
        e[k-1] = t;

        t      = b[k];
        b[k]   = b[k-1];
        b[k-1] = t;
      }

      if ( zabs1 ( c[k-1] ) == 0.0 )
      {
        info = k;
        return info;
      }

      t = -c[k] / c[k-1];
      c[k] = d[k] + t * d[k-1];
      d[k] = e[k] + t * e[k-1];
      e[k] = complex <double> ( 0.0, 0.0 );
      b[k] = b[k] + t * b[k-1];
    }
  }

  if ( zabs1 ( c[n-1] ) == 0.0 )
  {
    info = n;
    return info;
  }
//
//  Back solve.
//
  b[n-1] = b[n-1] / c[n-1];

  if ( 1 < n )
  {
    b[n-2] = ( b[n-2] - d[n-2] * b[n-1] ) / c[n-2];

    for ( k = n-2; 1 <= k; k-- )
    {
      b[k-1] = ( b[k-1] - d[k-1] * b[k] - e[k-1] * b[k+1] ) / c[k-1];
    }
  }

  return info;
}
//****************************************************************************80

double zhico ( complex <double> a[], int lda, int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZHICO factors a complex hermitian matrix and estimates its condition.
//
//  Discussion:
//
//    If RCOND is not needed, ZHIFA is slightly faster.
//
//    To solve A*X = B, follow ZHICO by ZHISL.
//
//    To compute inverse(A)*C, follow ZHICO by ZHISL.
//
//    To compute inverse(A), follow ZHICO by ZHIDI.
//
//    To compute determinant(A), follow ZHICO by ZHIDI.
//
//    To compute inertia(A), follow ZHICO by ZHIDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the hermitian matrix
//    to be factored.  On output, a block diagonal matrix and the multipliers
//    which were used to obtain it.  The factorization can be written
//    A = U*D*hermitian(U) where U is a product of permutation and unit
//    upper triangular matrices, hermitian(U) is the conjugate transpose
//    of U, and D is block diagonal with 1 by 1 and 2 by 2 blocks.
//    Only the diagonal and upper triangle are used.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, double ZHICO, an estimate of RCOND, the reciprocal condition of
//    the matrix.  For the system A*X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    (EPSILON/RCOND).  If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
//  Local Parameter:
//
//    Workspace, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
{
  complex <double> ak;
  complex <double> akm1;
  double anorm;
  complex <double> bk;
  complex <double> bkm1;
  complex <double> denom;
  complex <double> ek;
  int i;
  int info;
  int j;
  int k;
  int kp;
  int kps;
  int ks;
  double rcond;
  double s;
  complex <double> t;
  double ynorm;
  complex <double> *z;
//
//  Find norm of A using only upper half.
//
  z = new complex <double> [n];

  for ( j = 1; j <= n; j++ )
  {
    z[j-1] = complex <double> ( dzasum ( j, a+0+(j-1)*lda, 1 ), 0.0 );
    for ( i = 1; i <= j-1; i++ )
    {
      z[i-1] =
        complex <double> ( real ( z[i-1] ) + zabs1 ( a[i-1+(j-1)*lda] ), 0.0 );
    }
  }

  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    anorm = r8_max ( anorm, real ( z[j] ) );
  }
//
//  Factor.
//
  info = zhifa ( a, lda, n, ipvt );
//
//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where U*D*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve U*D*W = E.
//
  ek = complex <double> ( 1.0, 0.0 );
  for ( i = 0; i < n; i++ )
  {
    z[i] = complex <double> ( 0.0, 0.0 );
  }

  k = n;

  while ( 0 < k )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    kp = abs ( ipvt[k-1] );
    kps = k + 1 - ks;

    if ( kp != kps )
    {
      t        = z[kps-1];
      z[kps-1] = z[kp-1];
      z[kp-1]  = t;
    }

    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, z[k-1] );
    }

    z[k-1] = z[k-1] + ek;
    zaxpy ( k-ks, z[k-1], a+0+(k-1)*lda, 1, z, 1 );

    if ( ks != 1 )
    {
      if ( zabs1 ( z[k-2] ) != 0.0 )
      {
        ek = zsign1 ( ek, z[k-2] );
      }
      z[k-2] = z[k-2] + ek;
      zaxpy ( k-ks, z[k-2], a+0+(k-2)*lda, 1, z, 1 );
    }

    if ( ks != 2 )
    {
      if ( zabs1 ( a[k-1+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
      {
        s = zabs1 ( a[k-1+(k-1)*lda] ) / zabs1 ( z[k-1] );
        zdscal ( n, s, z, 1 );
        ek = complex <double> ( s, 0.0 ) * ek;
      }

      if ( zabs1 ( a[k-1+(k-1)*lda] ) != 0.0 )
      {
        z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
      }
      else
      {
        z[k-1] = complex <double> ( 1.0, 0.0 );
      }
    }
    else
    {
      ak = a[k-1+(k-1)*lda] / conj ( a[k-2+(k-1)*lda] );
      akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
      bk = z[k-1] / conj ( a[k-2+(k-1)*lda] );
      bkm1 = z[k-2] / a[k-2+(k-1)*lda];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      z[k-1] = ( akm1 * bk - bkm1 ) / denom;
      z[k-2] = ( ak * bkm1 - bk ) / denom;
    }
    k = k - ks;
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
//
//  Solve hermitian(U) * Y = W.
//
  k = 1;

  while ( k <= n )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    if ( k != 1 )
    {
      z[k-1] = z[k-1] + zdotc ( k-1, a+0+(k-1)*lda, 1, z, 1 );

      if ( ks == 2 )
      {
        z[k] = z[k] + zdotc ( k-1, a+0+k*lda, 1, z, 1 );
      }

      kp = abs ( ipvt[k-1] );

      if ( kp != k )
      {
        t       = z[k-1];
        z[k-1]  = z[kp-1];
        z[kp-1] = t;
      }
    }
    k = k + ks;
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = 1.0;
//
//  Solve U*D*V = Y.
//
  k = n;

  while ( 0 < k )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    if ( k != ks )
    {
      kp = abs ( ipvt[k-1] );
      kps = k + 1 - ks;

      if ( kp != kps )
      {
        t        = z[kps-1];
        z[kps-1] = z[kp-1];
        z[kp-1]  = t;
      }

      zaxpy ( k-ks, z[k-1], a+0+(k-1)*lda, 1, z, 1 );

      if ( ks == 2 )
      {
        zaxpy ( k-ks, z[k-2], a+0+(k-2)*lda, 1, z, 1 );
      }
    }

    if ( ks != 2 )
    {
      if ( zabs1 ( a[k-1+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
      {
        s = zabs1 ( a[k-1+(k-1)*lda] ) / zabs1 ( z[k-1] );
        zdscal ( n, s, z, 1 );
        ynorm = s * ynorm;
      }

      if ( zabs1 ( a[k-1+(k-1)*lda] ) != 0.0 )
      {
        z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
      }
      else
      {
        z[k-1] = complex <double> ( 1.0, 0.0 );
      }
    }
    else
    {
      ak = a[k-1+(k-1)*lda] / conj ( a[k-2+(k-1)*lda] );
      akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
      bk = z[k-1] / conj ( a[k-2+(k-1)*lda] );
      bkm1 = z[k-2] / a[k-2+(k-1)*lda];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      z[k-1] = ( akm1 * bk - bkm1 ) / denom;
      z[k-2] = ( ak * bkm1 - bk ) / denom;
    }
    k = k - ks;
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;
//
//  Solve hermitian(U) * Z = V.
//
  k = 1;

  while ( k <= n )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    if ( k != 1 )
    {
      z[k-1] = z[k-1] + zdotc ( k-1, a+0+(k-1)*lda, 1, z, 1 );

      if ( ks == 2 )
      {
        z[k] = z[k] + zdotc ( k-1, a+0+k*lda, 1, z, 1 );
      }

      kp = abs ( ipvt[k-1] );

      if ( kp != k )
      {
        t       = z[k-1];
        z[k-1]  = z[kp-1];
        z[kp-1] = t;
      }
    }
    k = k + ks;
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }

  delete [] z;

  return rcond;
}
//****************************************************************************80

void zhidi ( complex <double> a[], int lda, int n, int ipvt[], double det[2],
  int inert[3], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZHIDI computes the determinant and inverse of a matrix factored by ZHIFA.
//
//  Discussion:
//
//    ZHIDI computes the determinant, inertia (number of positive, zero,
//    and negative eigenvalues) and inverse of a complex hermitian matrix
//    using the factors from ZHIFA.
//
//    A division by zero may occur if the inverse is requested
//    and ZHICO has set RCOND == 0.0 or ZHIFA has set INFO /= 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the factored matrix
//    from ZHIFA.  On output, if the inverse was requested, A contains
//    the inverse matrix.  The strict lower triangle of A is never
//    referenced.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZHIFA.
//
//    Output, double DET[2], the determinant of the original matrix.
//    Determinant = det[0] * 10.0**det[1] with 1.0 <= abs ( det[0] ) < 10.0
//    or det[0] = 0.0.
//
//    Output, int INERT[3], the inertia of the original matrix.
//    INERT(1) = number of positive eigenvalues.
//    INERT(2) = number of negative eigenvalues.
//    INERT(3) = number of zero eigenvalues.
//
//    Input, int JOB, has the decimal expansion ABC where:
//    if C /= 0, the inverse is computed,
//    if B /= 0, the determinant is computed,
//    if A /= 0, the inertia is computed.
//    For example, JOB = 111 gives all three.
//
{
  double ak;
  complex <double> akkp1;
  double akp1;
  double d;
  int i;
  int j;
  int jb;
  int k;
  int km1;
  int ks;
  int kstep;
  bool nodet;
  bool noert;
  bool noinv;
  double t;
  complex <double> t2;
  complex <double> *work;

  noinv = ( job %   10 )       == 0;
  nodet = ( job %  100 ) /  10 == 0;
  noert = ( job % 1000 ) / 100 == 0;

  if ( !nodet || !noert )
  {
    if ( !noert )
    {
      for ( i = 0; i < 3; i++ )
      {
        inert[i] = 0;
      }
    }

    if ( !nodet )
    {
      det[0] = 1.0;
      det[1] = 0.0;
    }

    t = 0.0;

    for ( k = 0; k < n; k++ )
    {
      d = real ( a[k+k*lda] );
//
//  Check if 1 by 1.
//
      if ( ipvt[k] <= 0 )
      {
//
//  2 by 2 block
//  Use DET = ( D / T * C - T ) * T, T = abs ( S )
//  to avoid underflow/overflow troubles.
//  Take two passes through scaling.  Use T for flag.
//
        if ( t == 0.0 )
        {
          t = abs ( a[k+(k+1)*lda] );
          d = ( d / t ) * real ( a[k+1+(k+1)*lda] ) - t;
        }
        else
        {
          d = t;
          t = 0.0;
        }
      }

      if ( !noert )
      {
        if ( 0.0 < d )
        {
          inert[0] = inert[0] + 1;
        }
        else if ( d < 0.0 )
        {
          inert[1] = inert[1] + 1;
        }
        else if ( d == 0.0 )
        {
          inert[2] = inert[2] + 1;
        }
      }

      if ( !nodet )
      {
        det[0] = det[0] * d;

        if ( det[0] != 0.0 )
        {
          while ( fabs ( det[0] ) < 1.0 )
          {
            det[0] = det[0] * 10.0;
            det[1] = det[1] -  1.0;
          }

          while ( 10.0 <= fabs ( det[0] ) )
          {
            det[0] = det[0] / 10.0;
            det[1] = det[1] + 1.0;
          }
        }
      }
    }
  }
//
//  Compute inverse(A).
//
  if ( !noinv )
  {
    work = new complex <double> [n];

    k = 1;

    while ( k <= n )
    {
      km1 = k - 1;

      if ( 0 <= ipvt[k-1] )
      {
//
//  1 by 1
//
        a[k-1+(k-1)*lda] =
          complex <double> ( 1.0 / real ( a[k-1+(k-1)*lda] ), 0.0 );

        if ( 1 <= km1 )
        {
          for ( i = 1; i <= km1; i++ )
          {
            work[i-1] = a[i-1+(k-1)*lda];
          }

          for ( j = 1; j <= km1; j++ )
          {
            a[j-1+(k-1)*lda] = zdotc ( j, a+0+(j-1)*lda, 1, work, 1 );
            zaxpy ( j-1, work[j-1], a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
          }

          a[k-1+(k-1)*lda] = a[k-1+(k-1)*lda] + complex <double> (
            real ( zdotc ( km1, work, 1, a+0+(k-1)*lda, 1 ) ), 0.0 );
        }
        kstep = 1;
      }
      else
      {
//
//  2 by 2
//
        t = abs ( a[k-1+k*lda] );
        ak = real ( a[k-1+(k-1)*lda] ) / t;
        akp1 = real ( a[k+k*lda] ) / t;
        akkp1 = a[k-1+k*lda] / t;
        d = t * ( ak * akp1 - 1.0 );
        a[k-1+(k-1)*lda] = complex <double> ( akp1 / d, 0.0 );
        a[k+k*lda] = complex <double> ( ak / d, 0.0 );
        a[k-1+k*lda] = -akkp1 / d;

        if ( 1 <= km1 )
        {
          for ( i = 1; i <= km1; i++ )
          {
            work[i-1] = a[i-1+k*lda];
          }

          for ( j = 1; j <= km1; j++ )
          {
            a[j-1+k*lda] = zdotc ( j, a+0+(j-1)*lda, 1, work, 1 );
            zaxpy ( j-1, work[j-1], a+0+(j-1)*lda, 1, a+0+k*lda, 1 );
          }

          a[k+k*lda] = a[k+k*lda] + complex <double> (
            real ( zdotc ( km1, work, 1, a+0+k*lda, 1 ) ), 0.0 );

          a[k-1+k*lda] = a[k-1+k*lda]
            + zdotc ( km1, a+0+(k-1)*lda, 1, a+0+k*lda, 1 );

          for ( i = 1; i <= km1; i++ )
          {
            work[i-1] = a[i-1+(k-1)*lda];
          }

          for ( j = 1; j <= km1; j++ )
          {
            a[j-1+(k-1)*lda] = zdotc ( j, a+0+(j-1)*lda, 1, work, 1 );
            zaxpy ( j-1, work[j-1], a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
          }

          a[k-1+(k-1)*lda] = a[k-1+(k-1)*lda] + complex <double> (
            real ( zdotc ( km1, work, 1, a+0+(k-1)*lda, 1 ) ), 0.0 );
        }
        kstep = 2;
      }
//
//  Swap
//
      ks = abs ( ipvt[k-1] );

      if ( ks != k )
      {
        zswap ( ks, a+0+(ks-1)*lda, 1, a+0+(k-1)*lda, 1 );

        for ( j = k; ks <= j; j-- )
        {
          t2                = conj ( a[j-1+(k-1)*lda] );
          a[j-1+(k-1)*lda]  = conj ( a[ks-1+(j-1)*lda] );
          a[ks-1+(j-1)*lda] = t2;
        }

        if ( kstep != 1 )
        {
          t2            = a[ks-1+k*lda];
          a[ks-1+k*lda] = a[k-1+k*lda];
          a[k-1+k*lda]  = t2;
        }
      }
      k = k + kstep;
    }
    delete [] work;
  }
  return;
}
//****************************************************************************80

int zhifa ( complex <double> a[], int lda, int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZHIFA factors a complex hermitian matrix.
//
//  Discussion:
//
//    ZHIFA performs the factoring by elimination with symmetric pivoting.
//
//    To solve A*X = B, follow ZHIFA by ZHISL.
//
//    To compute inverse(A)*C, follow ZHIFA by ZHISL.
//
//    To compute determinant(A), follow ZHIFA by ZHIDI.
//
//    To compute inertia(A), follow ZHIFA by ZHIDI.
//
//    To compute inverse(A), follow ZHIFA by ZHIDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the hermitian matrix to be
//    factored.  On output, a block diagonal matrix and the multipliers which
//    were used to obtain it.  The factorization can be written
//    A = U*D*hermitian(U) where U is a product of permutation and unit upper
//    triangular matrices, hermitian(U) is the conjugate transpose of U, and
//    D is block diagonal with 1 by 1 and 2 by 2 blocks.  Only the diagonal
//    and upper triangle are used.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int ZHIFA.
//    0, normal value.
//    K, if the K-th pivot block is singular.  This is not an error condition
//    for this subroutine, but it does indicate that ZHISL or ZHIDI may
//    divide by zero if called.
//
{
  double absakk;
  complex <double> ak;
  complex <double> akm1;
  double alpha;
  complex <double> bk;
  complex <double> bkm1;
  double colmax;
  complex <double> denom;
  int imax;
  int info;
  int j;
  int jj;
  int jmax;
  int k;
  int km1;
  int km2;
  int kstep;
  complex <double> mulk;
  complex <double> mulkm1;
  double rowmax;
  bool swap;
  complex <double> t;
//
//  Initialize.
//
//  ALPHA is used in choosing pivot block size.
//
  alpha = ( 1.0 + sqrt ( 17.0 ) ) / 8.0;

  info = 0;
//
//  Main loop on K, which goes from N to 1.
//
  k = n;

  for ( ; ; )
  {
//
//  Leave the loop if K = 0 or K = 1.
//
    if ( k == 0 )
    {
      break;
    }

    if ( k == 1 )
    {
      ipvt[0] = 1;
      if ( zabs1 ( a[0+0*lda] ) == 0.0 )
      {
        info = 1;
      }
      break;
    }
//
//  This section of code determines the kind of
//  elimination to be performed.  When it is completed,
//  KSTEP will be set to the size of the pivot block, and
//  SWAP will be set to .true. if an interchange is
//  required.
//
    km1 = k - 1;
    absakk = zabs1 ( a[k-1+(k-1)*lda] );
//
//  Determine the largest off-diagonal element in column K.
//
    imax = izamax ( k-1, a+0+(k-1)*lda, 1 );
    colmax = zabs1 ( a[imax-1+(k-1)*lda] );

    if ( alpha * colmax <= absakk )
    {
      kstep = 1;
      swap = false;
    }
    else
    {
//
//  Determine the largest off-diagonal element in row IMAX.
//
      rowmax = 0.0;
      for ( j = imax + 1; j <= k; j++ )
      {
        rowmax = r8_max ( rowmax, zabs1 ( a[imax-1+(j-1)*lda] ) );
      }

      if ( imax != 1 )
      {
        jmax = izamax ( imax-1, a+0+(imax-1)*lda, 1 );
        rowmax = max ( rowmax, zabs1 ( a[jmax-1+(imax-1)*lda] ) );
      }

      if ( alpha * rowmax <= zabs1 ( a[imax-1+(imax-1)*lda] )  )
      {
        kstep = 1;
        swap = true;
      }
      else if ( alpha * colmax * ( colmax / rowmax ) <= absakk )
      {
        kstep = 1;
        swap = false;
      }
      else
      {
        kstep = 2;
        swap = ( imax != km1 );
      }
    }
//
//  Column K is zero.  Set INFO and iterate the loop.
//
    if ( r8_max ( absakk, colmax ) == 0.0 )
    {
      ipvt[k-1] = k;
      info = k;
      k = k - kstep;
      continue;
    }

    if ( kstep != 2 )
    {
//
//  1 x 1 pivot block.
//
      if ( swap )
      {
        zswap ( imax, a+0+(imax-1)*lda, 1, a+0+(k-1)*lda, 1 );

        for ( jj = imax; jj <= k; jj++ )
        {
          j = k + imax - jj;
          t                   = conj ( a[j-1+(k-1)*lda] );
          a[j-1+(k-1)*lda]    = conj ( a[imax-1+(j-1)*lda] );
          a[imax-1+(j-1)*lda] = t;
        }
      }
//
//  Perform the elimination.
//
      for ( jj = 1; jj <= km1; jj++ )
      {
        j = k - jj;
        mulk = -a[j-1+(k-1)*lda] / a[k-1+(k-1)*lda];
        t = conj ( mulk );
        zaxpy ( j, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
        a[j-1+(j-1)*lda] = complex <double> ( real ( a[j-1+(j-1)*lda] ), 0.0 );
        a[j-1+(k-1)*lda] = mulk;
      }
//
//  Set the pivot array.
//
      ipvt[k-1] = k;

      if ( swap )
      {
        ipvt[k-1] = imax;
      }
    }
    else
    {
//
//  2 x 2 pivot block.
//
      if ( swap )
      {
        zswap ( imax, a+0+(imax-1)*lda, 1, a+0+(k-2)*lda, 1 );

        for ( jj = imax; jj <= km1; jj++ )
        {
          j = km1 + imax - jj;

          t                   = conj ( a[j-1+(k-2)*lda] );
          a[j-1+(k-2)*lda]    = conj ( a[imax-1+(j-1)*lda] );
          a[imax-1+(j-1)*lda] = t;
        }
        t                   = a[k-2+(k-1)*lda];
        a[k-2+(k-1)*lda]    = a[imax-1+(k-1)*lda];
        a[imax-1+(k-1)*lda] = t;
      }
//
//  Perform the elimination.
//
      km2 = k - 2;

      if ( 0 < k - 2 )
      {
        ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
        akm1 = a[k-2+(k-2)*lda] / conj ( a[k-2+(k-1)*lda] );
        denom = complex <double> ( 1.0, 0.0 ) - ak * akm1;

        for ( jj = 1; jj <= k-2; jj++ )
        {
          j = km1 - jj;
          bk = a[j-1+(k-1)*lda] / a[k-2+(k-1)*lda];
          bkm1 = a[j-1+(k-2)*lda] / conj ( a[k-2+(k-1)*lda] );
          mulk = ( akm1 * bk - bkm1 ) / denom;
          mulkm1 = ( ak * bkm1 - bk ) / denom;
          t = conj ( mulk );
          zaxpy ( j, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
          t = conj ( mulkm1 );
          zaxpy ( j, t, a+0+(k-2)*lda, 1, a+0+(j-1)*lda, 1 );
          a[j-1+(k-1)*lda] = mulk;
          a[j-1+(k-2)*lda] = mulkm1;
          a[j-1+(j-1)*lda] = complex <double> ( real ( a[j-1+(j-1)*lda] ), 0.0 );
        }
      }
//
//  Set the pivot array.
//
      if ( swap )
      {
        ipvt[k-1] = -imax;
      }
      else
      {
        ipvt[k-1] = 1 - k;
      }

      ipvt[k-2] = ipvt[k-1];
    }
    k = k - kstep;
  }
  return info;
}
//*****************************************************************************

void zhisl ( complex <double> a[], int lda, int n, int ipvt[],
  complex <double> b[] )

//*****************************************************************************
//
//  Purpose:
//
//    ZHISL solves a complex hermitian system factored by ZHIFA.
//
//  Discussion:
//
//    A division by zero may occur if ZHICO has set RCOND == 0.0
//    or ZHIFA has set INFO != 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> A[LDA*N], the output from ZHIFA.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZHIFA.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  complex <double> ak;
  complex <double> akm1;
  complex <double> bk;
  complex <double> bkm1;
  complex <double> denom;
  int k;
  int kp;
  complex <double> t;
//
//  Loop backward applying the transformations and D inverse to B.
//
  k = n;

  while ( 0 < k )
  {
//
//  1 x 1 pivot block.
//
    if ( 0 <= ipvt[k-1] )
    {
      if ( k != 1 )
      {
        kp = ipvt[k-1];

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
        zaxpy ( k-1, b[k-1], a+0+(k-1)*lda, 1, b, 1 );
      }
//
//  Apply D inverse.
//
      b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
      k = k - 1;
    }
//
//  2 x 2 pivot block.
//
    else
    {
      if ( k != 2 )
      {
        kp = abs ( ipvt[k-1] );

        if ( kp != k - 1 )
        {
          t       = b[k-2];
          b[k-2]  = b[kp-1];
          b[kp-1] = t;
        }
        zaxpy ( k-2, b[k-1], a+0+(k-1)*lda, 1, b, 1 );
        zaxpy ( k-2, b[k-2], a+0+(k-2)*lda, 1, b, 1 );
      }
//
//  Apply D inverse.
//
      ak = a[k-1+(k-1)*lda] / conj ( a[k-2+(k-1)*lda] );
      akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
      bk = b[k-1] / conj ( a[k-2+(k-1)*lda] );
      bkm1 = b[k-2] / a[k-2+(k-1)*lda];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      b[k-1] = ( akm1 * bk - bkm1 ) / denom;
      b[k-2] = ( ak * bkm1 - bk ) / denom;
      k = k - 2;
    }
  }
//
//  Loop forward applying the transformations.
//
  k = 1;
  while ( k <= n )
  {
//
//  1 x 1 pivot block.
//
    if ( 0 <= ipvt[k-1] )
    {
      if ( k != 1 )
      {
        b[k-1] = b[k-1] + zdotc ( k-1, a+0+(k-1)*lda, 1, b, 1 );
        kp = ipvt[k-1];

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
      }
      k = k + 1;
    }
    else
    {
//
//  2 x 2 pivot block.
//
      if ( k != 1 )
      {
        b[k-1] = b[k-1] + zdotc ( k-1, a+0+(k-1)*lda, 1, b, 1 );
        b[k]   = b[k]   + zdotc ( k-1, a+0+k*lda, 1, b, 1 );
        kp = abs ( ipvt[k-1] );

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
      }
      k = k + 2;
    }
  }
  return;
}
//****************************************************************************80

double zhpco ( complex <double> ap[], int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZHPCO factors a complex hermitian packed matrix and estimates its condition.
//
//  Discussion:
//
//    If RCOND is not needed, ZHPFA is slightly faster.
//
//    To solve A*X = B, follow ZHPCO by ZHPSL.
//
//    To compute inverse(A)*C, follow ZHPCO by ZHPSL.
//
//    To compute inverse(A), follow ZHPCO by ZHPDI.
//
//    To compute determinant(A), follow ZHPCO by ZHPDI.
//
//    To compute inertia(A), follow ZHPCO by ZHPDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> AP[N*(N+1)/2]; on input, the packed form of a
//    hermitian matrix A.  The columns of the upper triangle are stored
//    sequentially in a one-dimensional array of length N*(N+1)/2.  On
//    output, a block diagonal matrix and the multipliers which were used
//    to obtain it stored in packed form.  The factorization can be written
//    A = U*D*hermitian(U) where U is a product of permutation and unit
//    upper triangular matrices, hermitian(U) is the conjugate transpose
//    of U, and D is block diagonal with 1 by 1 and 2 by 2 blocks.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, double ZHPCO, an estimate of RCOND, the reciprocal condition of
//    the matrix.  For the system A*X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    (EPSILON/RCOND).  If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
//  Local Parameters:
//
//    Workspace, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
{
  complex <double> ak;
  complex <double> akm1;
  double anorm;
  complex <double> bk;
  complex <double> bkm1;
  complex <double> denom;
  complex <double> ek;
  int i;
  int ij;
  int ik;
  int ikm1;
  int ikp1;
  int info;
  int j;
  int j1;
  int k;
  int kk;
  int km1k;
  int km1km1;
  int kp;
  int kps;
  int ks;
  double rcond;
  double s;
  complex <double> t;
  double ynorm;
  complex <double> *z;

  z = new complex <double> [n];
//
//  Find norm of A using only upper half.
//
  j1 = 1;

  for ( j = 1; j <= n; j++ )
  {
    z[j-1] = complex <double> ( dzasum ( j, ap+j1-1, 1 ), 0.0 );
    ij = j1;
    j1 = j1 + j;

    for ( i = 1; i <= j-1; i++ )
    {
      z[i-1] = complex <double> ( real ( z[i-1] ) + zabs1 ( ap[ij-1] ), 0.0 );
      ij = ij + 1;
    }
  }

  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    anorm = r8_max ( anorm, real ( z[j] ) );
  }
//
//  Factor.
//
  info = zhpfa ( ap, n, ipvt );
//
//  RCOND = 1/(norm(A) * (estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where U*D*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve U*D*W = E.
//
  ek = complex <double> ( 1.0, 0.0 );
  for ( i = 0; i < n; i++ )
  {
    z[i] = complex <double> ( 0.0, 0.0 );
  }

  k = n;
  ik = ( n * ( n - 1 ) ) / 2;

  while ( 0 < k )
  {
    kk = ik + k;
    ikm1 = ik - ( k - 1 );

    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    kp = abs ( ipvt[k-1] );
    kps = k + 1 - ks;

    if ( kp != kps )
    {
      t        = z[kps-1];
      z[kps-1] = z[kp-1];
      z[kp-1]  = t;
    }

    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, z[k-1] );
    }

    z[k-1] = z[k-1] + ek;
    zaxpy ( k-ks, z[k-1], ap+ik, 1, z, 1 );

    if ( ks != 1 )
    {
      if ( zabs1 ( z[k-2] ) != 0.0 )
      {
        ek = zsign1 ( ek, z[k-2] );
      }
      z[k-2] = z[k-2] + ek;
      zaxpy ( k-ks, z[k-2], ap+ikm1, 1, z, 1 );
    }

    if ( ks != 2 )
    {
      if ( zabs1 ( ap[kk-1] ) < zabs1 ( z[k-1] ) )
      {
        s = zabs1 ( ap[kk-1] ) / zabs1 ( z[k-1] );
        zdscal ( n, s, z, 1 );
        ek = complex <double> ( s, 0.0 ) * ek;
      }

      if ( zabs1 ( ap[kk-1] ) != 0.0 )
      {
        z[k-1] = z[k-1] / ap[kk-1];
      }
      else
      {
        z[k-1] = complex <double> ( 1.0, 0.0 );
      }
    }
    else
    {
      km1k = ik + k - 1;
      km1km1 = ikm1 + k - 1;
      ak = ap[kk-1] / conj ( ap[km1k-1] );
      akm1 = ap[km1km1-1] / ap[km1k-1];
      bk = z[k-1] / conj ( ap[km1k-1] );
      bkm1 = z[k-2] / ap[km1k-1];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      z[k-1] = ( akm1 * bk - bkm1 ) / denom;
      z[k-2] = ( ak * bkm1 - bk ) / denom;
    }

    k = k - ks;
    ik = ik - k;

    if ( ks == 2 )
    {
      ik = ik - ( k + 1 );
    }
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
//
//  Solve hermitian(U) * Y = W.
//
  k = 1;
  ik = 0;

  while ( k <= n )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    if ( k != 1 )
    {
      z[k-1] = z[k-1] + zdotc ( k-1, ap+ik, 1, z, 1 );
      ikp1 = ik + k;

      if ( ks == 2 )
      {
        z[k] = z[k] + zdotc ( k-1, ap+ikp1, 1, z, 1 );
      }

      kp = abs ( ipvt[k-1] );

      if ( kp != k )
      {
        t       = z[k-1];
        z[k-1]  = z[kp-1];
        z[kp-1] = t;
      }
    }

    ik = ik + k;
    if ( ks == 2 )
    {
      ik = ik + ( k + 1 );
    }
    k = k + ks;
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = 1.0;
//
//  Solve U*D*V = Y.
//
  k = n;
  ik = ( n * ( n - 1 ) ) / 2;

  while ( 0 < k )
  {
    kk = ik + k;
    ikm1 = ik - ( k - 1 );

    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    if ( k != ks )
    {
      kp = abs ( ipvt[k-1] );
      kps = k + 1 - ks;

      if ( kp != kps )
      {
        t        = z[kps-1];
        z[kps-1] = z[kp-1];
        z[kp-1]  = t;
      }

      zaxpy ( k-ks, z[k-1], ap+ik, 1, z, 1 );

      if ( ks == 2 )
      {
        zaxpy ( k-ks, z[k-2], ap+ikm1, 1, z, 1 );
      }
    }

    if ( ks != 2 )
    {
      if ( zabs1 ( ap[kk-1] ) < zabs1 ( z[k-1] ) )
      {
        s = zabs1 ( ap[kk-1] ) / zabs1 ( z[k-1] );
        zdscal ( n, s, z, 1 );
        ynorm = s * ynorm;
      }

      if ( zabs1 ( ap[kk-1] ) != 0.0 )
      {
        z[k-1] = z[k-1] / ap[kk-1];
      }
      else
      {
        z[k-1] = complex <double> ( 1.0, 0.0 );
      }
    }
    else
    {
      km1k = ik + k - 1;
      km1km1 = ikm1 + k - 1;
      ak = ap[kk-1] / conj ( ap[km1k-1] );
      akm1 = ap[km1km1-1] / ap[km1k-1];
      bk = z[k-1] / conj ( ap[km1k-1] );
      bkm1 = z[k-2] / ap[km1k-1];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      z[k-1] = ( akm1 * bk - bkm1 ) / denom;
      z[k-2] = ( ak * bkm1 - bk ) / denom;
    }

    k = k - ks;
    ik = ik - k;

    if ( ks == 2 )
    {
      ik = ik - ( k + 1 );
    }
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;
//
//  Solve hermitian(U) * Z = V.
//
  k = 1;
  ik = 0;

  while ( k <= n )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    if ( k != 1 )
    {
      z[k-1] = z[k-1] + zdotc ( k-1, ap+ik, 1, z, 1 );
      ikp1 = ik + k;

      if ( ks == 2 )
      {
        z[k] = z[k] + zdotc ( k-1, ap+ikp1, 1, z, 1 );
      }

      kp = abs ( ipvt[k-1] );

      if ( kp != k )
      {
        t       = z[k-1];
        z[k-1]  = z[kp-1];
        z[kp-1] = t;
      }
    }

    ik = ik + k;

    if ( ks == 2 )
    {
      ik = ik + ( k + 1 );
    }
    k = k + ks;
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }

  delete [] z;

  return rcond;
}
//****************************************************************************80

void zhpdi ( complex <double> ap[], int n, int ipvt[], double det[2],
  int inert[3], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZHPDI: determinant, inertia and inverse of a complex hermitian matrix.
//
//  Discussion:
//
//    The routine uses the factors from ZHPFA.
//
//    The matrix is stored in packed form.
//
//    A division by zero will occur if the inverse is requested and ZHPCO has
//    set RCOND == 0.0 or ZHPFA has set INFO != 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> AP[N*(N+1)/2]; on input, the factored matrix
//    from ZHPFA.  If the inverse was requested, then on output, AP contains
//    the upper triangle of the inverse of the original matrix, stored in packed
//    form.  The columns of the upper triangle are stored sequentially in a
//    one-dimensional array.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZHPFA.
//
//    Output, double DET[2], if requested, the determinant of the original
//    matrix.  Determinant = DET(1) * 10.0**DET(2) with
//    1.0 <= abs ( DET(1) ) < 10.0 or DET(1) = 0.0.
//
//    Output, int INERT[3], if requested, the inertia of the original matrix.
//    INERT(1) = number of positive eigenvalues.
//    INERT(2) = number of negative eigenvalues.
//    INERT(3) = number of zero eigenvalues.
//
//    Input, int JOB, has the decimal expansion ABC where:
//    if C != 0, the inverse is computed,
//    if B != 0, the determinant is computed,
//    if A != 0, the inertia is computed.
//    For example, JOB = 111 gives all three.
//
{
  double ak;
  complex <double> akkp1;
  double akp1;
  double d;
  int ij;
  int ik;
  int ikp1;
  int iks;
  int j;
  int jb;
  int jk;
  int jkp1;
  int k;
  int kk;
  int kkp1;
  int km1;
  int ks;
  int ksj;
  int kskp1;
  int kstep;
  bool nodet;
  bool noert;
  bool noinv;
  double t;
  complex <double> t2;
  complex <double> *work;

  noinv = ( job % 10 ) == 0;
  nodet = ( job % 100 ) / 10 == 0;
  noert = ( job % 1000 ) / 100 == 0;

  if ( !nodet || !noert )
  {
    if ( !noert )
    {
      inert[0] = 0;
      inert[1] = 0;
      inert[2] = 0;
    }

    if ( !nodet )
    {
      det[0] = 1.0;
      det[1] = 0.0;
    }

    t = 0.0;
    ik = 0;

    for ( k = 1; k <= n; k++ )
    {
      kk = ik + k;
      d = real ( ap[kk-1] );
//
//  Check if 1 by 1
//
      if ( ipvt[k-1] <= 0 )
      {
//
//  2 by 2 block
//  Use DET (D  S; S  C)  =  ( D / T * C - T ) * T, T = abs ( S )
//  to avoid underflow/overflow troubles.
//  Take two passes through scaling.  Use T for flag.
//
        if ( t == 0.0 )
        {
          ikp1 = ik + k;
          kkp1 = ikp1 + k;
          t = abs ( ap[kkp1-1] );
          d = ( d / t ) * real ( ap[kkp1] ) - t;
        }
        else
        {
          d = t;
          t = 0.0;
        }
      }

      if ( !noert )
      {
        if ( 0.0 < d )
        {
          inert[0] = inert[0] + 1;
        }
        else if ( d < 0.0 )
        {
          inert[1] = inert[1] + 1;
        }
        else if ( d == 0.0 )
        {
          inert[2] = inert[2] + 1;
        }
      }

      if ( !nodet )
      {
        det[0] = det[0] * d;

        if ( det[0] != 0.0 )
        {
          while ( fabs ( det[0] ) < 1.0 )
          {
            det[0] = det[0] * 10.0;
            det[1] = det[1] - 1.0;
          }
          while ( 10.0 <= fabs ( det[0] ) )
          {
            det[0] = det[0] / 10.0;
            det[1] = det[1] + 1.0;
          }
        }
      }
      ik = ik + k;
    }
  }
//
//  Compute inverse(A).
//
  if ( !noinv )
  {
    work = new complex <double> [n];

    k = 1;
    ik = 0;

    while ( k <= n )
    {
      km1 = k - 1;
      kk = ik + k;
      ikp1 = ik + k;
      kkp1 = ikp1 + k;
//
//  1 by 1
//
      if ( 0 <= ipvt[k-1] )
      {
        ap[kk-1] = complex <double> ( 1.0 / real ( ap[kk-1] ), 0.0 );

        if ( 1 <= km1 )
        {
          for ( j = 1; j <= km1; j++ )
          {
            work[j-1] = ap[ik+j-1];
          }
          ij = 0;
          for ( j = 1; j <= km1; j++ )
          {
            jk = ik + j;
            ap[jk-1] = zdotc ( j, ap+ij, 1, work, 1 );
            zaxpy ( j-1, work[j-1], ap+ij, 1, ap+ik, 1 );
            ij = ij + j;
          }
          ap[kk-1] = ap[kk-1] + complex <double>
            ( real ( zdotc ( km1, work, 1, ap+ik, 1) ), 0.0 );
        }
        kstep = 1;
      }
//
//  2 by 2
//
      else
      {
        t = abs ( ap[kkp1-1] );
        ak = real ( ap[kk-1] ) / t;
        akp1 = real ( ap[kkp1] ) / t;
        akkp1 = ap[kkp1-1] / t;
        d = t * ( ak * akp1 - 1.0 );
        ap[kk-1] = complex <double> ( akp1 / d, 0.0 );
        ap[kkp1] = complex <double> ( ak / d, 0.0 );
        ap[kkp1-1] = -akkp1 / d;

        if ( 1 <= km1 )
        {
          for ( j = 1; j <= km1; j++ )
          {
            work[j-1] = ap[ikp1+j-1];
          }
          ij = 0;
          for ( j = 1; j <= km1; j++ )
          {
            jkp1 = ikp1 + j;
            ap[jkp1-1] = zdotc ( j, ap+ij, 1, work, 1 );
            zaxpy ( j-1, work[j-1], ap+ij, 1, ap+ikp1, 1 );
            ij = ij + j;
          }

          ap[kkp1] = ap[kkp1] + complex <double>
            ( real ( zdotc ( km1, work, 1, ap+ikp1, 1 ) ), 0.0 );

          ap[kkp1-1] = ap[kkp1-1] + zdotc ( km1, ap+ik, 1, ap+ikp1, 1 );
          for ( j = 1; j <= km1; j++ )
          {
            work[j-1] = ap[ik+j-1];
          }
          ij = 0;

          for ( j = 1; j <= km1; j++ )
          {
            jk = ik + j;
            ap[jk-1] = zdotc ( j, ap+ij, 1, work, 1 );
            zaxpy ( j-1, work[j-1], ap+ij, 1, ap+ik, 1 );
            ij = ij + j;
          }
          ap[kk-1] = ap[kk-1] + complex <double>
            ( real ( zdotc ( km1, work, 1, ap+ik, 1 ) ), 0.0 );
        }
        kstep = 2;
      }
//
//  Swap
//
      ks = abs ( ipvt[k-1] );

      if ( ks != k )
      {
        iks = ( ks * ( ks - 1 ) ) / 2;

        zswap ( ks, ap+iks, 1, ap+ik, 1 );
        ksj = ik + ks;

        for ( jb = ks; jb <= k; jb++ )
        {
          j = k + ks - jb;
          jk = ik + j;

          t2        = conj ( ap[jk-1] );
          ap[jk-1]  = conj ( ap[ksj-1] );
          ap[ksj-1] = t2;

          ksj = ksj - ( j - 1 );
        }

        if ( kstep != 1 )
        {
          kskp1 = ikp1 + ks;

          t2          = ap[kskp1-1];
          ap[kskp1-1] = ap[kkp1-1];
          ap[kkp1-1]  = t2;
        }
      }

      ik = ik + k;

      if ( kstep == 2 )
      {
        ik = ik + k + 1;
      }
      k = k + kstep;
    }
    delete [] work;
  }

  return;
}
//****************************************************************************80

int zhpfa ( complex <double> ap[], int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZHPFA factors a complex hermitian packed matrix.
//
//  Discussion:
//
//    To solve A*X = B, follow ZHPFA by ZHPSL.
//
//    To compute inverse(A)*C, follow ZHPFA by ZHPSL.
//
//    To compute determinant(A), follow ZHPFA by ZHPDI.
//
//    To compute inertia(A), follow ZHPFA by ZHPDI.
//
//    To compute inverse(A), follow ZHPFA by ZHPDI.
//
//  Packed storage:
//
//    The following program segment will pack the upper
//    triangle of a hermitian matrix.
//
//      k = 0
//      do j = 1, n
//        do i = 1, j
//          k = k + 1
//          ap(k) = a(i,j)
//        }
//      }
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> AP[N*(N+1)/2]; on input, the packed form
//    of a hermitian matrix.  The columns of the upper triangle are
//    stored sequentially in a one-dimensional array.  On output, a
//    block diagonal matrix and the multipliers which were used to
//    obtain it stored in packed form.  The factorization can be
//    written A = U*D*hermitian(U) where U is a product of permutation
//    and unit upper triangular matrices , hermitian(U) is the
//    conjugate transpose of U, and D is block diagonal with 1 by 1
//    and 2 by 2 blocks.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int ZHPFA.
//    0, normal value.
//    K, if the K-th pivot block is singular.  This is not an error condition
//    for this subroutine, but it does indicate that ZHPSL or ZHPDI may divide
//    by zero if called.
//
{
  double absakk;
  complex <double> ak;
  complex <double> akm1;
  double alpha;
  complex <double> bk;
  complex <double> bkm1;
  double colmax;
  complex <double> denom;
  int ij;
  int ijj;
  int ik;
  int ikm1;
  int im;
  int imax;
  int imim;
  int imj;
  int imk;
  int info;
  int j;
  int jj;
  int jk;
  int jkm1;
  int jmax;
  int jmim;
  int k;
  int kk;
  int km1;
  int km1k;
  int km1km1;
  int km2;
  int kstep;
  complex <double> mulk;
  complex <double> mulkm1;
  double rowmax;
  bool swap;
  complex <double> t;
//
//  Initialize.
//
//  ALPHA is used in choosing pivot block size.
//
  alpha = ( 1.0 + sqrt ( 17.0 ) ) / 8.0;

  info = 0;
//
//  Main loop on K, which goes from N to 1.
//
  k = n;
  ik = ( n * ( n - 1 ) ) / 2;

  for ( ; ; )
  {
//
//  Leave the loop if K = 0 or K = 1.
//
    if ( k == 0 )
    {
      break;
    }

    if ( k == 1 )
    {
      ipvt[0] = 1;
      if ( zabs1 ( ap[0] ) == 0.0 )
      {
        info = 1;
      }
      break;
    }
//
//  This section of code determines the kind of
//  elimination to be performed.  When it is completed,
//  KSTEP will be set to the size of the pivot block, and
//  SWAP will be set to .true. if an interchange is
//  required.
//
    km1 = k - 1;
    kk = ik + k;
    absakk = zabs1 ( ap[kk-1] );
//
//  Determine the largest off-diagonal element in column K.
//
    imax = izamax ( k-1, ap+ik, 1 );
    imk = ik + imax;
    colmax = zabs1 ( ap[imk-1] );

    if ( alpha * colmax <= absakk )
    {
      kstep = 1;
      swap = false;
    }
//
//  Determine the largest off-diagonal element in row IMAX.
//
    else
    {
      rowmax = 0.0;
      im = imax * ( imax - 1 ) / 2;
      imj = im + 2 * imax;

      for ( j = imax + 1; j <= k; j++ )
      {
        rowmax = r8_max ( rowmax, zabs1 ( ap[imj-1] ) );
        imj = imj + j;
      }

      if ( imax != 1 )
      {
        jmax = izamax ( imax-1, ap+im, 1 );
        jmim = jmax + im;
        rowmax = max ( rowmax, zabs1 ( ap[jmim-1] ) );
      }

      imim = imax + im;

      if ( alpha * rowmax <= zabs1 ( ap[imim-1] ) )
      {
        kstep = 1;
        swap = true;
      }
      else if ( alpha * colmax * ( colmax / rowmax ) <= absakk )
      {
        kstep = 1;
        swap = false;
      }
      else
      {
        kstep = 2;
        swap = ( imax != km1 );
      }
    }
//
//  Column K is zero.  Set INFO and iterate the loop.
//
    if ( r8_max ( absakk, colmax ) == 0.0 )
    {
      ipvt[k-1] = k;
      info = k;
      ik = ik - ( k - 1 );
      if ( kstep == 2 )
      {
        ik = ik - ( k - 2 );
      }
      k = k - kstep;
      continue;
    }

    if ( kstep != 2 )
    {
//
//  1 x 1 pivot block.
//
      if ( swap )
      {
        zswap ( imax, ap+im, 1, ap+ik, 1 );
        imj = ik + imax;

        for ( jj = imax; jj <= k; jj++ )
        {
          j = k + imax - jj;
          jk = ik + j;

          t         = conj ( ap[jk-1] );
          ap[jk-1]  = conj ( ap[imj-1] );
          ap[imj-1] = t;

          imj = imj - ( j - 1 );
        }
      }
//
//  Perform the elimination.
//
      ij = ik - ( k - 1 );
      for ( jj = 1; jj <= km1; jj++ )
      {
        j = k - jj;
        jk = ik + j;
        mulk = -ap[jk-1] / ap[kk-1];
        t = conj ( mulk );
        zaxpy ( j, t, ap+ik, 1, ap+ij, 1 );
        ijj = ij + j;
        ap[ijj-1] = complex <double> ( real ( ap[ijj-1] ), 0.0 );
        ap[jk-1] = mulk;
        ij = ij - ( j - 1 );
      }
//
//  Set the pivot array.
//
      if ( swap )
      {
        ipvt[k-1] = imax;
      }
      else
      {
        ipvt[k-1] = k;
      }
    }
//
//  2 x 2 pivot block.
//
    else
    {
      km1k = ik + k - 1;
      ikm1 = ik - ( k - 1 );

      if ( swap )
      {
        zswap ( imax, ap+im, 1, ap+ikm1, 1 );
        imj = ikm1 + imax;

        for ( jj = imax; jj <= km1; jj++ )
        {
          j = km1 + imax - jj;
          jkm1 = ikm1 + j;

          t          = conj ( ap[jkm1-1] );
          ap[jkm1-1] = conj ( ap[imj-1] );
          ap[imj-1]   = t;

          imj = imj - ( j - 1 );
        }
        t          = ap[km1k-1];
        ap[km1k-1] = ap[imk-1];
        ap[imk-1]  = t;
      }
//
//  Perform the elimination.
//
      km2 = k - 2;

      if ( km2 != 0 )
      {
        ak = ap[kk-1] / ap[km1k-1];
        km1km1 = ikm1 + k - 1;
        akm1 = ap[km1km1-1] / conj ( ap[km1k-1] );
        denom = complex <double> ( 1.0, 0.0 ) - ak * akm1;
        ij = ik - ( k - 1 ) - ( k - 2 );

        for ( jj = 1; jj <= km2; jj++ )
        {
          j = km1 - jj;
          jk = ik + j;
          bk = ap[jk-1] / ap[km1k-1];
          jkm1 = ikm1 + j;
          bkm1 = ap[jkm1-1] / conj ( ap[km1k-1] );
          mulk = ( akm1 * bk - bkm1 ) / denom;
          mulkm1 = ( ak * bkm1 - bk ) / denom;
          t = conj ( mulk );
          zaxpy ( j, t, ap+ik, 1, ap+ij, 1 );
          t = conj ( mulkm1 );
          zaxpy ( j, t, ap+ikm1, 1, ap+ij, 1 );
          ap[jk-1] = mulk;
          ap[jkm1-1] = mulkm1;
          ijj = ij + j;
          ap[ijj-1] = complex <double> ( real ( ap[ijj-1] ), 0.0 );
          ij = ij - ( j - 1 );
        }
      }
//
//  Set the pivot array.
//
      if ( swap )
      {
        ipvt[k-1] = -imax;
      }
      else
      {
        ipvt[k-1] = 1 - k;
      }
      ipvt[k-2] = ipvt[k-1];
    }

    ik = ik - ( k - 1 );
    if ( kstep == 2 )
    {
      ik = ik - ( k - 2 );
    }
    k = k - kstep;
  }

  return info;
}
//****************************************************************************80

void zhpsl ( complex <double> ap[], int n, int ipvt[], complex <double> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZHPSL solves a complex hermitian system factored by ZHPFA.
//
//  Discussion:
//
//    A division by zero may occur if ZHPCO set RCOND to 0.0
//    or ZHPFA set INFO nonzero.
//
//    To compute
//
//      inverse ( A ) * C
//
//    where C is a matrix with P columns
//
//      call zhpfa(ap,n,ipvt,info)
//
//      if ( info == 0 )
//        do j = 1, p
//          call zhpsl(ap,n,ipvt,c(1,j))
//        end do
//      }
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> AP[N*(N+1)/2], the output from ZHPFA.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZHPFA.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  complex <double> ak;
  complex <double> akm1;
  complex <double> bk;
  complex <double> bkm1;
  complex <double> denom;
  int ik;
  int ikm1;
  int ikp1;
  int k;
  int kk;
  int km1k;
  int km1km1;
  int kp;
  complex <double> t;
//
//  Loop backward applying the transformations and inverse ( D ) to B.
//
  k = n;
  ik = ( n * ( n - 1 ) ) / 2;

  while ( 0 < k )
  {
    kk = ik + k;
//
//  1 x 1 pivot block.
//
    if ( 0 <= ipvt[k-1] )
    {
      if ( k != 1 )
      {
        kp = ipvt[k-1];

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
        zaxpy ( k-1, b[k-1], ap+ik, 1, b, 1 );
      }
//
//  Apply D inverse.
//
      b[k-1] = b[k-1] / ap[kk-1];
      k = k - 1;
      ik = ik - k;
    }
    else
    {
//
//  2 x 2 pivot block.
//
      ikm1 = ik - ( k - 1 );

      if ( k != 2 )
      {
        kp = abs ( ipvt[k-1] );

        if ( kp != k - 1 )
        {
          t       = b[k-2];
          b[k-2]  = b[kp-1];
          b[kp-1] = t;
        }
        zaxpy ( k-2, b[k-1], ap+ik, 1, b, 1 );
        zaxpy ( k-2, b[k-2], ap+ikm1, 1, b, 1 );
      }
//
//  Apply D inverse.
//
      km1k = ik + k - 1;
      kk = ik + k;
      ak = ap[kk-1] / conj ( ap[km1k-1] );
      km1km1 = ikm1 + k - 1;
      akm1 = ap[km1km1-1] / ap[km1k-1];
      bk = b[k-1] / conj ( ap[km1k-1] );
      bkm1 = b[k-2] / ap[km1k-1];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      b[k-1] = ( akm1 * bk - bkm1 ) / denom;
      b[k-2] = ( ak * bkm1 - bk ) / denom;
      k = k - 2;
      ik = ik - ( k + 1 ) - k;
    }
  }
//
//  Loop forward applying the transformations.
//
  k = 1;
  ik = 0;

  while ( k <= n )
  {
//
//  1 x 1 pivot block.
//
    if ( 0 <= ipvt[k-1] )
    {
      if ( k != 1 )
      {
        b[k-1] = b[k-1] + zdotc ( k-1, ap+ik, 1, b, 1 );
        kp = ipvt[k-1];

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
      }
      ik = ik + k;
      k = k + 1;
    }
//
//  2 x 2 pivot block.
//
    else
    {
      if ( k != 1 )
      {
        b[k-1] = b[k-1] + zdotc ( k-1, ap+ik, 1, b, 1 );
        ikp1 = ik + k;
        b[k] = b[k] + zdotc ( k-1, ap+ikp1, 1, b, 1 );
        kp = abs ( ipvt[k-1] );

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
      }
      ik = ik + k + k + 1;
      k = k + 2;
    }
  }
  return;
}
//****************************************************************************80

double zpbco ( complex <double> abd[], int lda, int n, int m, int *info )

//****************************************************************************80
//
//  Purpose:
//
//    ZPBCO factors a complex <double> hermitian positive definite band matrix.
//
//  Discussion:
//
//    The routine also estimates the condition number of the matrix.
//
//    If RCOND is not needed, ZPBFA is slightly faster.
//
//    To solve A*X = B, follow ZPBCO by ZPBSL.
//
//    To compute inverse(A)*C, follow ZPBCO by ZPBSL.
//
//    To compute determinant(A), follow ZPBCO by ZPBDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> ABD[LDA*N]; on input, the matrix to be factored.
//    The columns of the upper triangle are stored in the columns of ABD,
//    and the diagonals of the upper triangle are stored in the rows of ABD.
//    On output, an upper triangular matrix R, stored in band form, so that
//    A = hermitian(R) * R.  If INFO != 0, the factorization is not complete.
//
//    Input, int LDA, the leading dimension of ABD.
//    LDA must be at least M+1.
//
//    Input, int N, the order of the matrix.
//
//    Input, int M, the number of diagonals above the main diagonal.
//    0 <= M < N.
//
//    Output, double ZPBCO, an estimate of RCOND, the reciprocal condition of
//    the matrix.  For the system A*X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    (EPSILON/RCOND).  If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
//    Output, int *INFO.
//    0, for normal return.
//    K, signals an error condition.  The leading minor of order K is not
//    positive definite.
//
//  Local Parameter:
//
//    Workspace, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If A is singular to working precision, then Z is
//    an approximate null vector in the sense that
//    norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
//    If INFO != 0, Z is unchanged.
//
{
  double anorm;
  complex <double> ek;
  int i;
  int j;
  int j2;
  int k;
  int kp1;
  int l;
  int la;
  int lb;
  int lm;
  int mu;
  double rcond;
  double s;
  double sm;
  complex <double> t;
  complex <double> wk;
  complex <double> wkm;
  double ynorm;
  complex <double> *z;
//
//  Find the norm of A.
//
  z = new complex <double> [n];

  for ( j = 1; j <= n; j++ )
  {
    l = i4_min ( j, m + 1 );
    mu = i4_max ( m + 2 - j, 1 );
    z[j-1] = complex <double> ( dzasum ( l, abd+mu-1+(j-1)*lda, 1 ), 0.0 );
    k = j - l;

    for ( i = mu; i <= m; i++ )
    {
      k = k + 1;
      z[k-1] = complex <double> ( real ( z[k-1] )
        + zabs1 ( abd[i-1+(j-1)*lda] ), 0.0 );
    }
  }

  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    anorm = r8_max ( anorm, real ( z[j] ) );
  }
//
//  Factor.
//
  *info = zpbfa ( abd, lda, n, m );

  if ( *info != 0 )
  {
    rcond = 0.0;
    delete [] z;
    return rcond;
  }
//
//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where hermitian(R)*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve hermitian(R)*W = E.
//
  ek = complex <double> ( 1.0, 0.0 );

  for ( i = 0; i < n; i++ )
  {
    z[i] = complex <double> ( 0.0, 0.0 );
  }

  for ( k = 1; k <= n; k++ )
  {
    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, -z[k-1] );
    }

    if ( real ( abd[m+(k-1)*lda] ) < zabs1 ( ek - z[k-1] ) )
    {
      s = real ( abd[m+(k-1)*lda] ) / zabs1 ( ek - z[k-1] );
      zdscal ( n, s, z, 1 );
      ek = complex <double> ( s, 0.0 ) * ek;
    }

    wk = ek - z[k-1];
    wkm = - ek - z[k-1];
    s = zabs1 ( wk );
    sm = zabs1 ( wkm );
    wk = wk / abd[m+(k-1)*lda];
    wkm = wkm / abd[m+(k-1)*lda];
    j2 = i4_min ( k + m, n );
    i = m + 1;

    if ( k+1 <= j2 )
    {
      for ( j = k + 1; j <= j2; j++ )
      {
        i = i - 1;
        sm = sm + zabs1 ( z[j-1] + wkm * conj ( abd[i-1+(j-1)*lda] ) );
        z[j-1] = z[j-1] + wk * conj ( abd[i-1+(j-1)*lda] );
        s = s + zabs1 ( z[j-1] );
      }

      if ( s < sm )
      {
        t = wkm - wk;
        wk = wkm;
        i = m + 1;
        for ( j = k + 1; j <= j2; j++ )
        {
          i = i - 1;
          z[j-1] = z[j-1] + t * conj ( abd[i-1+(j-1)*lda] );
        }
      }
    }
    z[k-1] = wk;
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
//
//  Solve R * Y = W.
//
  for ( k = n; 1 <= k; k-- )
  {
    if ( real ( abd[m+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
    {
      s = real ( abd[m+(k-1)*lda] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
    }
    z[k-1] = z[k-1] / abd[m+(k-1)*lda];
    lm = i4_min ( k - 1, m );
    la = m + 1 - lm;
    lb = k - lm;
    t = -z[k-1];
    zaxpy ( lm, t, abd+la-1+(k-1)*lda, 1, z+lb-1, 1 );
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = 1.0;
//
//  Solve hermitian(R)*V = Y.
//
  for ( k = 1; k <= n; k++ )
  {
    lm = i4_min ( k - 1, m );
    la = m + 1 - lm;
    lb = k - lm;
    z[k-1] = z[k-1] - zdotc ( lm, abd+la-1+(k-1)*lda, 1, z+lb-1, 1 );

    if ( real ( abd[m+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
    {
      s = real ( abd[m+(k-1)*lda] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }
    z[k-1] = z[k-1] / abd[m+(k-1)*lda];
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;
//
//  Solve R * Z = W.
//
  for ( k = n; 1 <= k; k-- )
  {
    if ( real ( abd[m+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
    {
      s = real ( abd[m+(k-1)*lda] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }
    z[k-1] = z[k-1] / abd[m+(k-1)*lda];
    lm = i4_min ( k - 1, m );
    la = m + 1 - lm;
    lb = k - lm;
    t = -z[k-1];
    zaxpy ( lm, t, abd+la-1+(k-1)*lda, 1, z+lb-1, 1 );
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }
  delete [] z;

  return rcond;
}
//****************************************************************************80

void zpbdi ( complex <double> abd[], int lda, int n, int m, double det[2] )

//****************************************************************************80
//
//  Purpose:
//
//    ZPBDI gets the determinant of a hermitian positive definite band matrix.
//
//  Discussion:
//
//    ZPBDI uses the factors computed by ZPBCO or ZPBFA.
//
//    If the inverse is needed, use ZPBSL N times.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> ABD[LDA*N], the output from ZPBCO or ZPBFA.
//
//    Input, int LDA, the leading dimension of the array ABD.
//
//    Input, int N, the order of the matrix.
//
//    Input, int M, the number of diagonals above the main diagonal.
//
//    Output, double DET[2], the determinant of the original matrix in the
//    form determinant = DET(1) * 10.0**DET(2) with 1.0 <= DET(1) < 10.0
//    or DET(1) == 0.0.
//
{
  int i;

  det[0] = 1.0;
  det[1] = 0.0;

  for ( i = 1; i <= n; i++ )
  {
    det[0] = det[0] * real ( abd[m+(i-1)*lda] ) * real ( abd[m+(i-1)*lda] );

    if ( det[0] == 0.0 )
    {
      break;
    }

    while ( det[0] < 1.0 )
    {
      det[0] = det[0] * 10.0;
      det[1] = det[1] - 1.0;
    }

    while ( 10.0 <= det[0] )
    {
      det[0] = det[0] / 10.0;
      det[1] = det[1] + 1.0;
    }

  }
  return;
}
//****************************************************************************80

int zpbfa ( complex <double> abd[], int lda, int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    ZPBFA factors a complex hermitian positive definite band matrix.
//
//  Discussion:
//
//    ZPBFA is usually called by ZPBCO, but it can be called
//    directly with a saving in time if RCOND is not needed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> ABD[LDA*N]; on input, the matrix to be factored.
//    The columns of the upper triangle are stored in the columns of ABD
//    and the diagonals of the upper triangle are stored in the rows of ABD.
//    On output, an upper triangular matrix R, stored in band form, so that
//    A = hermitian(R)*R.
//
//    Input, int LDA, the leading dimension of ABD.
//    LDA must be at least M+1.
//
//    Input, int N, the order of the matrix.
//
//    Input, int M, the number of diagonals above the main diagonal.
//    0 <= M < N.
//
//    Output, int ZSPFA.
//    0, for normal return.
//    K, if the leading minor of order K is not positive definite.
//
{
  int ik;
  int info;
  int j;
  int jk;
  int k;
  int mu;
  double s;
  complex <double> t;

  info = 0;

  for ( j = 1; j <= n; j++ )
  {
    s = 0.0;
    ik = m + 1;
    jk = i4_max ( j - m, 1 );
    mu = i4_max ( m + 2 - j, 1 );

    for ( k = mu; k <= m; k++ )
    {
      t = abd[k-1+(j-1)*lda]
        - zdotc ( k-mu, abd+ik-1+(jk-1)*lda, 1, abd+mu-1+(j-1)*lda, 1 );
      t = t / abd[m+(jk-1)*lda];
      abd[k-1+(j-1)*lda] = t;
      s = s + real ( t * conj ( t ) );
      ik = ik - 1;
      jk = jk + 1;
    }

    s = real ( abd[m+(j-1)*lda] ) - s;

    if ( s <= 0.0 || imag ( abd[m+(j-1)*lda] ) != 0.0 )
    {
      info = j;
      break;
    }
    abd[m+(j-1)*lda] = complex <double> ( sqrt ( s ), 0.0 );
  }

  return info;
}
//****************************************************************************80

void zpbsl ( complex <double> abd[], int lda, int n, int m,
  complex <double> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZPBSL solves a complex hermitian positive definite band system.
//
//  Discussion:
//
//    The system matrix must have been factored by ZPBCO or ZPBFA.
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal.  Technically this indicates
//    singularity but it is usually caused by improper subroutine
//    arguments.  It will not occur if the subroutines are called
//    correctly and INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> ABD[LDA*N], the output from ZPBCO or ZPBFA.
//
//    Input, int LDA, the leading dimension of ABD.
//
//    Input, int N, the order of the matrix.
//
//    Input, int M, the number of diagonals above the main diagonal.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  int k;
  int la;
  int lb;
  int lm;
  complex <double> t;
//
//  Solve hermitian(R) * Y = B.
//
  for ( k = 1; k <= n; k++ )
  {
    lm = i4_min ( k - 1, m );
    la = m + 1 - lm;
    lb = k - lm;
    t = zdotc ( lm, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
    b[k-1] = ( b[k-1] - t ) / abd[m+(k-1)*lda];
  }
//
//  Solve R * X = Y.
//
  for ( k = n; 1 <= k; k-- )
  {
    lm = i4_min ( k - 1, m );
    la = m + 1 - lm;
    lb = k - lm;
    b[k-1] = b[k-1] / abd[m+(k-1)*lda];
    t = -b[k-1];
    zaxpy ( lm, t, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
  }

  return;
}
//****************************************************************************80

double zpoco ( complex <double> a[], int lda, int n, int *info )

//****************************************************************************80
//
//  Purpose:
//
//    ZPOCO factors a complex hermitian positive definite matrix.
//
//  Discussion:
//
//    The routine also estimates the condition of the matrix.
//
//    If RCOND is not needed, ZPOFA is slightly faster.
//
//    To solve A*X = B, follow ZPOCO by ZPOSL.
//
//    To compute inverse(A)*C, follow ZPOCO by ZPOSL.
//
//    To compute determinant(A), follow ZPOCO by ZPODI.
//
//    To compute inverse(A), follow ZPOCO by ZPODI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the hermitian matrix to be
//    factored.  On output, an upper triangular matrix R so that
//      A = hermitian(R)*R
//    where hermitian(R) is the conjugate transpose.  The strict lower
//    triangle is unaltered.  If INFO /= 0, the factorization is not complete.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int *INFO.
//    0, for normal return.
//    K, signals an error condition.  The leading minor of order K is not
//    positive definite.
//
//    Output, double ZPOCO, an estimate of RCOND, the reciprocal condition of
//    the matrix.  For the system A*X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    (EPSILON/RCOND).  If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
//  Local Parameters:
//
//    Workspace, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
{
  double anorm;
  complex <double> ek;
  int i;
  int j;
  int k;
  int kp1;
  double rcond;
  double s;
  double sm;
  complex <double> t;
  complex <double> wk;
  complex <double> wkm;
  double ynorm;
  complex <double> *z;
//
//  Find norm of A using only upper half.
//
  z = new complex <double> [n];

  for ( j = 1; j <= n; j++ )
  {
    z[j-1] = complex <double> ( dzasum ( j, a+0+(j-1)*lda, 1 ), 0.0 );
    for ( i = 1; i < j; i++ )
    {
      z[i-1] =
       complex <double> ( real ( z[i-1] ) + zabs1 ( a[i-1+(j-1)*lda] ), 0.0 );
    }
  }
  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    anorm = r8_max ( anorm, real ( z[j] ) );
  }
//
//  Factor.
//
  *info = zpofa ( a, lda, n );

  if ( *info != 0 )
  {
    rcond = 0.0;
    return rcond;
  }
//
//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where hermitian(R)*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve hermitian(R)*W = E.
//
  ek = complex <double> ( 1.0, 0.0 );
  for ( j = 0; j < n; j++ )
  {
    z[j] = complex <double> ( 0.0, 0.0 );
  }

  for ( k = 1; k <= n; k++ )
  {
    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, -z[k-1] );
    }

    if ( real ( a[k-1+(k-1)*lda] ) < zabs1 ( ek - z[k-1] ) )
    {
      s = real ( a[k-1+(k-1)*lda] ) / zabs1 ( ek - z[k-1] );
      zdscal ( n, s, z, 1 );
      ek = complex <double> ( s, 0.0 ) * ek;
    }

    wk = ek - z[k-1];
    wkm = -ek - z[k-1];
    s = zabs1 ( wk );
    sm = zabs1 ( wkm );
    wk = wk / a[k-1+(k-1)*lda];
    wkm = wkm / a[k-1+(k-1)*lda];
    kp1 = k + 1;

    if ( kp1 <= n )
    {
      for ( j = kp1; j <= n; j++ )
      {
        sm = sm + zabs1 ( z[j-1] + wkm * conj ( a[k-1+(j-1)*lda] ) );
        z[j-1] = z[j-1] + wk * conj ( a[k-1+(j-1)*lda] );
        s = s + zabs1 ( z[j-1] );
      }

      if ( s < sm )
      {
        t = wkm - wk;
        wk = wkm;
        for ( j = kp1; j <= n; j++ )
        {
          z[j-1] = z[j-1] + t * conj ( a[k-1+(j-1)*lda] );
        }
      }
    }
    z[k-1] = wk;
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
//
//  Solve R * Y = W.
//
  for ( k = n; 1 <= k; k-- )
  {
    if ( real ( a[k-1+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
    {
      s = real ( a[k-1+(k-1)*lda] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
    }

    z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
    t = -z[k-1];
    zaxpy ( k-1, t, a+0+(k-1)*lda, 1, z, 1 );
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = 1.0;
//
//  Solve hermitian(R) * V = Y.
//
  for ( k = 1; k <= n; k++ )
  {
    z[k-1] = z[k-1] - zdotc ( k-1, a+0+(k-1)*lda, 1, z, 1 );

    if ( real ( a[k-1+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
    {
      s = real ( a[k-1+(k-1)*lda] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }
    z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;
//
//  Solve R * Z = V.
//
  for ( k = n; 1 <= k; k-- )
  {
    if ( real ( a[k-1+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
    {
      s = real ( a[k-1+(k-1)*lda] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }
    z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
    t = -z[k-1];
    zaxpy ( k-1, t, a+0+(k-1)*lda, 1, z, 1 );
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }

  delete [] z;

  return rcond;
}
//****************************************************************************80

void zpodi ( complex <double> a[], int lda, int n, double det[2], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZPODI: determinant, inverse of a complex hermitian positive definite matrix.
//
//  Discussion:
//
//    The matrix is assumed to have been factored by ZPOCO, ZPOFA or ZQRDC.
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal and the inverse is requested.
//    It will not occur if the subroutines are called correctly
//    and if ZPOCO or ZPOFA has set INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the output A from ZPOCO or
//    ZPOFA, or the output X from ZQRDC.  On output, if ZPOCO or ZPOFA was
//    used to factor A, then ZPODI produces the upper half of inverse(A).
//    If ZQRDC was used to decompose X, then ZPODI produces the upper half
//    of inverse(hermitian(X)*X) where hermitian(X) is the conjugate transpose.
//    Elements of A below the diagonal are unchanged.
//    If the units digit of JOB is zero, A is unchanged.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Output, double DET[2], if requested, the determinant of A or of
//    hermitian(X)*X.  Determinant = DET(1) * 10.0**DET(2) with
//    1.0 <= abs ( DET(1) ) < 10.0 or DET(1) = 0.0.
//
//    Input, int JOB.
//    11, both determinant and inverse.
//    01, inverse only.
//    10, determinant only.
//
{
  int i;
  int j;
  int k;
  complex <double> t;
//
//  Compute determinant
//
  if ( ( job / 10 ) != 0 )
  {
    det[0] = 1.0;
    det[1] = 0.0;

    for ( i = 0; i < n; i++ )
    {
      det[0] = det[0] * real ( a[i+i*lda] ) * real ( a[i+i*lda] );

      if ( det[0] == 0.0 )
      {
        break;
      }
      while ( det[0] < 1.0 )
      {
        det[0] = det[0] * 10.0;
        det[1] = det[1] - 1.0;
      }
      while ( 10.0 <= det[0] )
      {
        det[0] = det[0] / 10.0;
        det[1] = det[1] + 1.0;
      }
    }
  }
//
//  Compute inverse(R).
//
  if ( ( job % 10 ) != 0 )
  {
    for ( k = 1; k <= n; k++ )
    {
      a[k-1+(k-1)*lda] = complex <double> ( 1.0, 0.0 ) / a[k-1+(k-1)*lda];
      t = -a[k-1+(k-1)*lda];
      zscal ( k-1, t, a+0+(k-1)*lda, 1 );

      for ( j = k+1; j <= n; j++ )
      {
        t = a[k-1+(j-1)*lda];
        a[k-1+(j-1)*lda] = complex <double> ( 0.0, 0.0 );
        zaxpy ( k, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
      }
    }
//
//  Form inverse(R) * hermitian(inverse(R)).
//
    for ( j = 1; j <= n; j++ )
    {
      for ( k = 1; k <= j-1; k++ )
      {
        t = conj ( a[k-1+(j-1)*lda] );
        zaxpy ( k, t, a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
      }

      t = conj ( a[j-1+(j-1)*lda] );
      zscal ( j, t, a+0+(j-1)*lda, 1 );
    }
  }
  return;
}
//****************************************************************************80

int zpofa ( complex <double> a[], int lda, int n )

//****************************************************************************80
//
//  Purpose:
//
//    ZPOFA factors a complex hermitian positive definite matrix.
//
//  Discussion:
//
//    ZPOFA is usually called by ZPOCO, but it can be called
//    directly with a saving in time if RCOND is not needed.
//    (time for ZPOCO) = (1 + 18/N) * (time for ZPOFA).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; On input, the hermitian matrix to be
//    factored.  On output, an upper triangular matrix R so that
//      A = hermitian(R)*R
//    where hermitian(R) is the conjugate transpose.  The strict lower
//    triangle is unaltered.  If INFO /= 0, the factorization is not
//    complete.  Only the diagonal and upper triangle are used.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int ZPOFA.
//    0, for normal return.
//    K, signals an error condition.  The leading minor of order K is
//    not positive definite.
//
{
  int info;
  int j;
  int k;
  double s;
  complex <double> t;

  info = 0;

  for ( j = 1; j <= n; j++ )
  {
    s = 0.0;
    for ( k = 1; k <= j-1; k++ )
    {
      t = a[k-1+(j-1)*lda] - zdotc ( k-1, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
      t = t / a[k-1+(k-1)*lda];
      a[k-1+(j-1)*lda] = t;
      s = s + real ( t * conj ( t ) );
    }

    s = real ( a[j-1+(j-1)*lda] ) - s;

    if ( s <= 0.0 || imag ( a[j-1+(j-1)*lda] ) != 0.0 )
    {
      info = j;
      break;
    }
    a[j-1+(j-1)*lda] = complex <double> ( sqrt ( s ), 0.0 );
  }
  return info;
}
//****************************************************************************80

void zposl ( complex <double> a[], int lda, int n, complex <double> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZPOSL solves a complex hermitian positive definite system.
//
//  Discussion:
//
//    ZPOSL uses the factors computed by ZPOCO or ZPOFA.
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal.  Technically this indicates
//    singularity but it is usually caused by improper subroutine
//    arguments.  It will not occur if the subroutines are called
//    correctly and INFO == 0.
//
//    To compute inverse(A) * C where C is a matrix with  p  columns
//
//      call zpoco(a,lda,n,rcond,z,info)
//
//      if (rcond is too small .or. info /= 0) then
//        error
//      end if
//
//      do j = 1, p
//        call zposl(a,lda,n,c(1,j))
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> A[LDA*N], the output from ZPOCO or ZPOFA.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  int k;
  complex <double> t;
//
//  Solve hermitian(R) * Y = B.
//
  for ( k = 1; k <= n; k++ )
  {
    t = zdotc ( k-1, a+0+(k-1)*lda, 1, b, 1 );
    b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
  }
//
//  Solve R * X = Y.
//
  for ( k = n; 1 <= k; k-- )
  {
    b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
    t = -b[k-1];
    zaxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
  }
  return;
}
//****************************************************************************80

double zppco ( complex <double> ap[], int n, int *info )

//****************************************************************************80
//
//  Purpose:
//
//    ZPPCO factors a complex <double> hermitian positive definite matrix.
//
//  Discussion:
//
//    The routine also estimates the condition of the matrix.
//
//    The matrix is stored in packed form.
//
//    If RCOND is not needed, ZPPFA is slightly faster.
//
//    To solve A*X = B, follow ZPPCO by ZPPSL.
//
//    To compute inverse(A)*C, follow ZPPCO by ZPPSL.
//
//    To compute determinant(A), follow ZPPCO by ZPPDI.
//
//    To compute inverse(A), follow ZPPCO by ZPPDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> AP[N*(N+1)/2]; on input, the packed form of a
//    hermitian matrix.  The columns of the upper triangle are stored
//    sequentially in a one-dimensional array.  On output, an upper
//    triangular matrix R, stored in packed form, so that A = hermitian(R) * R.
//    If INFO != 0 , the factorization is not complete.
//
//    Input, int N, the order of the matrix.
//
//    Output, double ZPPCO, an estimate of RCOND, the reciprocal condition of
//    the matrix.  For the system A*X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    (EPSILON/RCOND).  If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
//    Output, int *INFO.
//    0, for normal return.
//    K, signals an error condition.  The leading minor of order K is not
//    positive definite.
//
//  Local Parameters:
//
//    Local, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
{
  double anorm;
  complex <double> ek;
  int i;
  int ij;
  int j;
  int j1;
  int k;
  int kj;
  int kk;
  double rcond;
  double s;
  double sm;
  complex <double> t;
  complex <double> wk;
  complex <double> wkm;
  double ynorm;
  complex <double> *z;
//
//  Find norm of A.
//
  z = new complex <double> [n];

  j1 = 1;

  for ( j = 1; j <= n; j++ )
  {
    z[j-1] = complex <double> ( dzasum ( j, ap+j1-1, 1 ), 0.0 );
    ij = j1;
    j1 = j1 + j;

    for ( i = 1; i <= j-1; i++ )
    {
      z[i-1] = complex <double> ( real ( z[i-1] ) + zabs1 ( ap[ij-1] ), 0.0 );
      ij = ij + 1;
    }
  }
  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    anorm = r8_max ( anorm, real ( z[j] ) );
  }
//
//  Factor.
//
  *info = zppfa ( ap, n );

  if ( *info != 0 )
  {
    delete [] z;
    rcond = 0.0;
    return rcond;
  }
//
//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where hermitian(R)*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve hermitian(R)*W = E.
//
  ek = complex <double> ( 1.0, 0.0 );
  for ( j = 0; j < n; j++ )
  {
    z[j] = complex <double> ( 0.0, 0.0 );
  }

  kk = 0;

  for ( k = 1; k <= n; k++ )
  {
    kk = kk + k;

    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, -z[k-1] );
    }

    if ( real ( ap[kk-1] ) < zabs1 ( ek - z[k-1] ) )
    {
      s = real ( ap[kk-1] ) / zabs1 ( ek - z[k-1] );
      zdscal ( n, s, z, 1 );
      ek = complex <double> ( s, 0.0 ) * ek;
    }

    wk = ek - z[k-1];
    wkm = -ek - z[k-1];
    s = zabs1 ( wk );
    sm = zabs1 ( wkm );
    wk = wk / ap[kk-1];
    wkm = wkm / ap[kk-1];
    kj = kk + k;

    if ( k+1 <= n )
    {
      for ( j = k + 1; j <= n; j++ )
      {
        sm = sm + zabs1 ( z[j-1] + wkm * conj ( ap[kj-1] ) );
        z[j-1] = z[j-1] + wk * conj ( ap[kj-1] );
        s = s + zabs1 ( z[j-1] );
        kj = kj + j;
      }

      if ( s < sm )
      {
        t = wkm - wk;
        wk = wkm;
        kj = kk + k;
        for ( j = k + 1; j <= n; j++ )
        {
          z[j-1] = z[j-1] + t * conj ( ap[kj-1] );
          kj = kj + j;
        }
      }
    }
    z[k-1] = wk;
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
//
//  Solve R * Y = W.
//
  for ( k = n; 1 <= k; k-- )
  {
    if ( real ( ap[kk-1] ) < zabs1 ( z[k-1] ) )
    {
      s = real ( ap[kk-1] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
    }
    z[k-1] = z[k-1] / ap[kk-1];
    kk = kk - k;
    t = -z[k-1];
    zaxpy ( k-1, t, ap+kk, 1, z, 1 );
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = 1.0;
//
//  Solve hermitian(R) * V = Y.
//
  for ( k = 1; k <= n; k++ )
  {
    z[k-1] = z[k-1] - zdotc ( k-1, ap+kk, 1, z, 1 );
    kk = kk + k;

    if ( real ( ap[kk-1] ) < zabs1 ( z[k-1] ) )
    {
      s = real ( ap[kk-1] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }
    z[k-1] = z[k-1] / ap[kk-1];
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;
//
//  Solve R * Z = V.
//
  for ( k = n; 1 <= k; k-- )
  {
    if ( real ( ap[kk-1] ) < zabs1 ( z[k-1] ) )
    {
      s = real ( ap[kk-1] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }
    z[k-1] = z[k-1] / ap[kk-1];
    kk = kk - k;
    t = -z[k-1];
    zaxpy ( k-1, t, ap+kk, 1, z, 1 );
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }
  delete [] z;
  return rcond;
}
//****************************************************************************80

void zppdi ( complex <double> ap[], int n, double det[2], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZPPDI: determinant, inverse of a complex hermitian positive definite matrix.
//
//  Discussion:
//
//    The matrix is assumed to have been factored by ZPPCO or ZPPFA.
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal and the inverse is requested.
//    It will not occur if the subroutines are called correctly
//    and if ZPOCO or ZPOFA has set INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[(N*(N+1))/2]; on input, the output from ZPPCO
//    or ZPPFA.  On output, the upper triangular half of the inverse.
//    The strict lower triangle is unaltered.
//
//    Input, int N, the order of the matrix.
//
//    Output, double DET[2], the determinant of original matrix if requested.
//    Otherwise not referenced.  Determinant = DET(1) * 10.0**DET(2)
//    with 1.0 <= DET(1) < 10.0 or DET(1) == 0.0.
//
//    Input, int JOB.
//    11, both determinant and inverse.
//    01, inverse only.
//    10, determinant only.
//
{
  int i;
  int ii;
  int j;
  int j1;
  int jj;
  int k;
  int k1;
  int kj;
  int kk;
  int kp1;
  complex <double> t;
//
//  Compute determinant.
//
  if ( ( job / 10 ) != 0 )
  {
    det[0] = 1.0;
    det[1] = 0.0;
    ii = 0;

    for ( i = 1; i <= n; i++ )
    {
      ii = ii + i;
      det[0] = det[0] * real ( ap[ii-1] ) * real ( ap[ii-1] );

      if ( det[0] == 0.0 )
      {
        break;
      }

      while ( det[0] < 1.0 )
      {
        det[0] = det[0] * 10.0;
        det[1] = det[1] - 1.0;
      }

      while ( 10.0 <= det[0] )
      {
        det[0] = det[0] / 10.0;
        det[1] = det[1] + 1.0;
      }
    }
  }
//
//  Compute inverse ( R ).
//
  if ( ( job % 10 ) != 0 )
  {
    kk = 0;

    for ( k = 1; k <= n; k++ )
    {
      k1 = kk + 1;
      kk = kk + k;
      ap[kk-1] = complex <double> ( 1.0, 0.0 ) / ap[kk-1];
      t = -ap[kk-1];
      zscal ( k-1, t, ap+k1-1, 1 );
      kp1 = k + 1;
      j1 = kk + 1;
      kj = kk + k;

      for ( j = kp1; j <= n; j++ )
      {
        t = ap[kj-1];
        ap[kj-1] = complex <double> ( 0.0, 0.0 );
        zaxpy ( k, t, ap+k1-1, 1, ap+j1-1, 1 );
        j1 = j1 + j;
        kj = kj + j;
      }
    }
//
//  Form inverse ( R ) * hermitian ( inverse ( R ) ).
//
    jj = 0;
    for ( j = 1; j <= n; j++ )
    {
      j1 = jj + 1;
      jj = jj + j;
      k1 = 1;
      kj = j1;

      for ( k = 1; k <= j-1; k++ )
      {
        t = conj ( ap[kj-1] );
        zaxpy ( k, t, ap+j1-1, 1, ap+k1-1, 1 );
        k1 = k1 + k;
        kj = kj + 1;
      }
      t = conj ( ap[jj-1] );
      zscal ( j, t, ap+j1-1, 1 );
    }
  }
  return;
}
//****************************************************************************80

int zppfa ( complex <double> ap[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    ZPPFA factors a complex hermitian positive definite packed matrix.
//
//  Discussion:
//
//    The following program segment will pack the upper triangle of a
//    hermitian matrix.
//
//      k = 0
//      do j = 1, n
//        do i = 1, j
//          k = k + 1
//          ap(k) = a(i,j)
//        end do
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> AP[N*(N+1)/2]; on input, the packed form
//    of a hermitian matrix A.  The columns of the upper triangle are
//    stored sequentially in a one-dimensional array.  On output, an
//    upper triangular matrix R, stored in packed form, so that
//      A = hermitian(R) * R.
//
//    Input, int N, the order of the matrix.
//
//    Output, int ZPPFA.
//    0, for normal return.
//    K, if the leading minor of order K is not positive definite.
//
{
  int info;
  int j;
  int jj;
  int k;
  int kj;
  int kk;
  double s;
  complex <double> t;

  info = 0;
  jj = 0;

  for ( j = 1; j <= n; j++ )
  {
    s = 0.0;
    kj = jj;
    kk = 0;

    for ( k = 1; k <= j-1; k++ )
    {
      kj = kj + 1;
      t = ap[kj-1] - zdotc ( k-1, ap+kk, 1, ap+jj, 1 );
      kk = kk + k;
      t = t / ap[kk-1];
      ap[kj-1] = t;
      s = s + real ( t * conj ( t ) );
    }

    jj = jj + j;
    s = real ( ap[jj-1] ) - s;

    if ( s <= 0.0 || imag ( ap[jj-1] ) != 0.0 )
    {
      info = j;
      break;
    }
    ap[jj-1] = complex <double> ( sqrt ( s ), 0.0 );
  }
  return info;
}
//****************************************************************************80

void zppsl ( complex <double> ap[], int n, complex <double> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZPPSL solves a complex hermitian positive definite linear system.
//
//  Discussion:
//
//    The matrix is assumed to have been factored by ZPPCO or ZPPFA.
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal.  Technically this indicates
//    singularity but it is usually caused by improper subroutine
//    arguments.  It will not occur if the subroutines are called
//    correctly and INFO == 0.
//
//    To compute inverse(A) * C where C is a matrix with P columns:
//
//      call zppco(ap,n,rcond,z,info)
//
//      if (rcond is too small .or. info /= 0) then
//        error
//      end if
//
//      do j = 1, p
//        call zppsl(ap,n,c(1,j))
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> AP[N*(N+1)/2], the output from ZPPCO or ZPPFA.
//
//    Input, int N, the order of the matrix.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  int k;
  int kk;
  complex <double> t;

  kk = 0;
  for ( k = 1; k <= n; k++ )
  {
    t = zdotc ( k-1, ap+kk, 1, b, 1 );
    kk = kk + k;
    b[k-1] = ( b[k-1] - t ) / ap[kk-1];
  }

  for ( k = n; 1 <= k; k-- )
  {
    b[k-1] = b[k-1] / ap[kk-1];
    kk = kk - k;
    t = -b[k-1];
    zaxpy ( k-1, t, ap+kk, 1, b, 1 );
  }

  return;
}
//****************************************************************************80

void zptsl ( int n, complex <double> d[], complex <double> e[],
  complex <double> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZPTSL solves a Hermitian positive definite tridiagonal linear system.
//
//  Discussion;
//
//    The system does not have to be factored first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, complex <double> D[N].  On input, the diagonal of the
//    matrix.  On output, this has been overwritten by other information.
//
//    Input/output, complex <double> E[N].  On input, the superdiagonal
//    entries of the matrix in locations E(1:N-1).  On output, this has
//    been overwritten by other information.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  int k;
  int kbm1;
  int ke;
  int kf;
  int kp1;
  int nm1;
  int nm1d2;
  complex <double> t1;
  complex <double> t2;
//
//  Check for 1 x 1 case.
//
  if ( n == 1 )
  {
    b[0] = b[0] / d[0];
    return;
  }

  nm1 = n - 1;
  nm1d2 = ( n - 1 ) / 2;

  if ( n != 2 )
  {
    kbm1 = n - 1;
//
//  Zero top half of subdiagonal and bottom half of superdiagonal.
//
    for ( k = 1; k <= nm1d2; k++ )
    {
      t1 = conj ( e[k-1] ) / d[k-1];
      d[k] = d[k] - t1 * e[k-1];
      b[k] = b[k] - t1 * b[k-1];
      t2 = e[kbm1-1] / d[kbm1];
      d[kbm1-1] = d[kbm1-1] - t2 * conj ( e[kbm1-1] );
      b[kbm1-1] = b[kbm1-1] - t2 * b[kbm1];
      kbm1 = kbm1 - 1;
    }
  }

  kp1 = nm1d2 + 1;
//
//  Clean up for possible 2 x 2 block at center.
//
  if ( ( n % 2 ) == 0 )
  {
    t1 = conj ( e[kp1-1] ) / d[kp1-1];
    d[kp1] = d[kp1] - t1 * e[kp1-1];
    b[kp1] = b[kp1] - t1 * b[kp1-1];
    kp1 = kp1 + 1;
  }
//
//  Back solve starting at the center, going towards the top and bottom.
//
  b[kp1-1] = b[kp1-1] / d[kp1-1];

  if ( n != 2 )
  {
    k = kp1 - 1;
    ke = kp1 + nm1d2 - 1;

    for ( kf = kp1; kf <= ke; kf++ )
    {
      b[k-1] = ( b[k-1] - e[k-1] * b[k] ) / d[k-1];
      b[kf] = ( b[kf] - conj ( e[kf-1] ) * b[kf-1] ) / d[kf];
      k = k - 1;
    }
  }

  if ( ( n % 2 ) == 0 )
  {
    b[0] = ( b[0] - e[0] * b[1] ) / d[0];
  }

  return;
}
//****************************************************************************80

void zqrdc ( complex <double> x[], int ldx, int n, int p,
  complex <double> qraux[], int ipvt[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZQRDC computes the QR factorization of an N by P complex <double> matrix.
//
//  Discussion:
//
//    ZQRDC uses Householder transformations to compute the QR factorization
//    of an N by P matrix X.  Column pivoting based on the 2-norms of the
//    reduced columns may be performed at the user's option.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> X[LDX*P]; on input, the matrix whose decomposition
//    is to be computed.  On output, the upper triangle contains the upper
//    triangular matrix R of the QR factorization.  Below its diagonal, X
//    contains information from which the unitary part of the decomposition
//    can be recovered.  If pivoting has been requested, the decomposition is
//    not that of the original matrix X, but that of X with its columns
//    permuted as described by IPVT.
//
//    Input, int LDX, the leading dimension of X.  N <= LDX.
//
//    Input, int N, the number of rows of the matrix.
//
//    Input, int P, the number of columns in the matrix X.
//
//    Output, complex <double> QRAUX[P], further information required to recover
//    the unitary part of the decomposition.
//
//    Input/output, int IPVT[P]; on input, ints that control the
//    selection of the pivot columns.  The K-th column X(K) of X is placed
//    in one of three classes according to the value of IPVT(K):
//      IPVT(K) > 0, then X(K) is an initial column.
//      IPVT(K) == 0, then X(K) is a free column.
//      IPVT(K) < 0, then X(K) is a final column.
//    Before the decomposition is computed, initial columns are moved to the
//    beginning of the array X and final columns to the end.  Both initial
//    and final columns are frozen in place during the computation and only
//    free columns are moved.  At the K-th stage of the reduction, if X(K)
//    is occupied by a free column it is interchanged with the free column
//    of largest reduced norm.
//    On output, IPVT(K) contains the index of the column of the
//    original matrix that has been interchanged into
//    the K-th column, if pivoting was requested.
//    IPVT is not referenced if JOB == 0.
//
//    Input, int JOB, initiates column pivoting.
//    0, no pivoting is done.
//    nonzero, pivoting is done.
//
{
  int itemp;
  int j;
  int jj;
  int l;
  int lp1;
  int lup;
  int maxj;
  double maxnrm;
  bool negj;
  complex <double> nrmxl;
  int pl;
  int pu;
  bool swapj;
  complex <double> t;
  double tt;
  complex <double> *work;

  pl = 1;
  pu = 0;
  work = new complex <double> [p];

  if ( job != 0 )
  {
//
//  Pivoting has been requested.  Rearrange the columns according to IPVT.
//
    for ( j = 1; j <= p; j++ )
    {
      swapj = ( 0 < ipvt[j-1] );
      negj = ( ipvt[j-1] < 0 );

      if ( negj )
      {
        ipvt[j-1] = -j;
      }
      else
      {
        ipvt[j-1] = j;
      }

      if ( swapj )
      {
        if ( j != pl )
        {
          zswap ( n, x+0+(pl-1)*ldx, 1, x+0+(j-1)*ldx, 1 );
        }
        ipvt[j-1] = ipvt[pl-1];
        ipvt[pl-1] = j;
        pl = pl + 1;
      }
    }
    pu = p;

    for ( jj = 1; jj <= p; jj++ )
    {
      j = p - jj + 1;

      if ( ipvt[j-1] < 0 )
      {
        ipvt[j-1] = -ipvt[j-1];

        if ( j != pu )
        {
          zswap ( n, x+0+(pu-1)*ldx, 1, x+0+(j-1)*ldx, 1 );

          itemp      = ipvt[pu-1];
          ipvt[pu-1] = ipvt[j-1];
          ipvt[j-1]  = itemp;
        }
        pu = pu - 1;
      }
    }
  }
//
//  Compute the norms of the free columns.
//
  for ( j = pl; j <= pu; j++ )
  {
    qraux[j-1] = complex <double> ( dznrm2 ( n, x+0+(j-1)*ldx, 1 ), 0.0 );
    work[j-1] = qraux[j-1];
  }
//
//  Perform the Householder reduction of X.
//
  lup = i4_min ( n, p );

  for ( l = 1; l <= lup; l++ )
  {
//
//  Locate the column of largest norm and bring it
//  into the pivot position.
//
    if ( pl <= l && l < pu )
    {
      maxnrm = 0.0;
      maxj = l;

      for ( j = l; j <= pu; j++ )
      {
        if ( maxnrm < real ( qraux[j-1] ) )
        {
          maxnrm = real ( qraux[j-1] );
          maxj = j;
        }
      }

      if ( maxj != l )
      {
        zswap ( n, x+0+(l-1)*ldx, 1, x+0+(maxj-1)*ldx, 1 );
        qraux[maxj-1] = qraux[l-1];
        work[maxj-1] = work[l-1];

        itemp        = ipvt[maxj-1];
        ipvt[maxj-1] = ipvt[l-1];
        ipvt[l-1]    = itemp;
      }
    }
    qraux[l-1] = complex <double> ( 0.0, 0.0 );

    if ( l != n )
    {
//
//  Compute the Householder transformation for column L.
//
      nrmxl = complex <double> ( dznrm2 ( n-l+1, x+l-1+(l-1)*ldx, 1 ), 0.0 );

      if ( zabs1 ( nrmxl ) != 0.0 )
      {
        if ( zabs1 ( x[l-1+(l-1)*ldx] ) != 0.0 )
        {
          nrmxl = zsign2 ( nrmxl, x[l-1+(l-1)*ldx] );
        }

        t = complex <double> ( 1.0, 0.0 ) / nrmxl;
        zscal ( n-l+1, t, x+l-1+(l-1)*ldx, 1 );
        x[l-1+(l-1)*ldx] = complex <double> ( 1.0, 0.0 ) + x[l-1+(l-1)*ldx];
//
//  Apply the transformation to the remaining columns,
//  updating the norms.
//
        lp1 = l + 1;

        for ( j = l+1; j <= p; j++ )
        {
          t = -zdotc ( n-l+1, x+l-1+(l-1)*ldx, 1, x+l-1+(j-1)*ldx, 1 )
            / x[l-1+(l-1)*ldx];
          zaxpy ( n-l+1, t, x+l-1+(l-1)*ldx, 1, x+l-1+(j-1)*ldx, 1 );

          if ( j < pl || pu < j )
          {
            continue;
          }

          if ( zabs1 ( qraux[j-1] ) == 0.0 )
          {
            continue;
          }

          tt = 1.0 - pow ( abs ( x[l-1+(j-1)*ldx] ) / real ( qraux[j-1] ), 2 );
          tt = r8_max ( tt, 0.0 );
          t = complex <double> ( tt, 0.0 );
          tt = 1.0 + 0.05 * tt
            * pow ( real ( qraux[j-1] ) / real ( work[j-1] ), 2 );

          if ( tt != 1.0 )
          {
            qraux[j-1] = qraux[j-1] * sqrt ( t );
          }
          else
          {
            qraux[j-1] =
              complex <double> ( dznrm2 ( n-l, x+l+(j-1)*ldx, 1 ), 0.0 );
            work[j-1] = qraux[j-1];
          }
        }
//
//  Save the transformation.
//
        qraux[l-1] = x[l-1+(l-1)*ldx];
        x[l-1+(l-1)*ldx] = -nrmxl;
      }
    }
  }
  delete [] work;

  return;
}
//****************************************************************************80

int zqrsl ( complex <double> x[], int ldx, int n, int k,
  complex <double> qraux[], complex <double> y[], complex <double> qy[],
  complex <double> qty[], complex <double> b[], complex <double> rsd[],
  complex <double> xb[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZQRSL solves, transforms or projects systems factored by ZQRDC.
//
//  Discussion:
//
//    The routine applies the output of ZQRDC to compute coordinate
//    transformations, projections, and least squares solutions.
//
//    For K <= min ( N, P ), let XK be the matrix
//
//      XK = ( X(IPVT(1)), X(IPVT(2)), ... ,X(IPVT(k)) )
//
//    formed from columnns IPVT(1), ... ,IPVT(K) of the original
//    N by P matrix X that was input to ZQRDC (if no pivoting was
//    done, XK consists of the first K columns of X in their
//    original order).  ZQRDC produces a factored unitary matrix Q
//    and an upper triangular matrix R such that
//
//      XK = Q * ( R )
//               ( 0 )
//
//    This information is contained in coded form in the arrays
//    X and QRAUX.
//
//    The parameters QY, QTY, B, RSD, and XB are not referenced
//    if their computation is not requested and in this case
//    can be replaced by dummy variables in the calling program.
//
//    To save storage, the user may in some cases use the same
//    array for different parameters in the calling sequence.  A
//    frequently occuring example is when one wishes to compute
//    any of B, RSD, or XB and does not need Y or QTY.  In this
//    case one may identify Y, QTY, and one of B, RSD, or XB, while
//    providing separate arrays for anything else that is to be
//    computed.  Thus the calling sequence
//
//      zqrsl ( x, ldx, n, k, qraux, y, dum, y, b, y, dum, 110, info )
//
//    will result in the computation of B and RSD, with RSD
//    overwriting Y.  More generally, each item in the following
//    list contains groups of permissible identifications for
//    a single callinng sequence.
//
//    1. ( Y, QTY, B )   ( RSD )      ( XB )  ( QY )
//    2. ( Y, QTY, RSD ) ( B )        ( XB )  ( QY )
//    3. ( Y, QTY, XB )  ( B )        ( RSD ) ( QY )
//    4. ( Y, QY )       ( QTY, B )   ( RSD ) ( XB )
//    5. ( Y, QY )       ( QTY, RSD ) ( B )   ( XB )
//    6. ( Y, QY )       ( QTY, XB )  ( B )   ( RSD )
//
//    In any group the value returned in the array allocated to
//    the group corresponds to the last member of the group.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> X[LDX*P], the output of ZQRDC.
//
//    Input, int LDX, the leading dimension of X.
//
//    Input, int N, the number of rows of the matrix XK, which
//    must have the same value as N in ZQRDC.
//
//    Input, int K, the number of columns of the matrix XK.  K must not
//    be greater than min ( N, P), where P is the same as in the calling
//    sequence to ZQRDC.
//
//    Input, complex <double> QRAUX[P], the auxiliary output from ZQRDC.
//
//    Input, complex <double> Y[N], a vector that is to be manipulated by ZQRSL.
//
//    Output, complex <double> QY[N], contains Q*Y, if it has been requested.
//
//    Output, complex <double> QTY[N], contains hermitian(Q)*Y, if it has
//    been requested.  Here hermitian(Q) is the conjugate transpose
//    of the matrix Q.
//
//    Output, complex <double> B[K], the solution of the least squares problem
//      minimize norm2 ( Y - XK * B ),
//    if it has been requested.  If pivoting was requested in ZQRDC,
//    the J-th component of B will be associated with column IPVT(J)
//    of the original matrix X that was input into ZQRDC.
//
//    Output, complex <double> RSD[N], the least squares residual Y - XK*B,
//    if it has been requested.  RSD is also the orthogonal projection
//    of Y onto the orthogonal complement of the column space of XK.
//
//    Output, complex <double> XB[N], the least squares approximation XK*N,
//    if its computation has been requested.  XB is also the orthogonal
//    projection of Y onto the column space of X.
//
//    Input, int JOB, specifies what is to be computed.  JOB has
//    the decimal expansion ABCDE, meaning:
//    if A != 0, compute QY.
//    if B, D, D, or E != 0, compute QTY.
//    if C != 0, compute B.
//    if D != 0, compute RSD.
//    if E != 0, compute XB.
//    A request to compute B, RSD, or XB automatically triggers the
//    computation of QTY, for which an array must be provided in the
//    calling sequence.
//
//    Output, int ZQRSL, the value of INFO, which is zero unless
//    the computation of B has been requested and R is exactly singular.
//    In this case, INFO is the index of the first zero diagonal element
//    of R and B is left unaltered.
//
{
  bool cb;
  bool cqty;
  bool cqy;
  bool cr;
  bool cxb;
  int i;
  int info;
  int j;
  int jj;
  int ju;
  int kp1;
  complex <double> t;
  complex <double> temp;

  info = 0;
//
//  Determine what is to be computed.
//
  cqy =  (   job / 10000         != 0 );
  cqty = ( ( job %  10000 )      != 0 );
  cb =   ( ( job %   1000 ) /100 != 0 );
  cr =   ( ( job %    100 ) / 10 != 0 );
  cxb =  ( ( job %     10 )      != 0 );

  ju = i4_min ( k, n - 1 );
//
//  Special action when N=1.
//
  if ( ju == 0 )
  {
    if ( cqy )
    {
      qy[0] = y[0];
    }
    if ( cqty )
    {
      qty[0] = y[0];
    }
    if ( cxb )
    {
      xb[0] = y[0];
    }
    if ( cb )
    {
      if ( zabs1 ( x[0+0*ldx] ) == 0.0 )
      {
        info = 1;
      }
      else
      {
        b[0] = y[0] / x[0+0*ldx];
      }
    }
    if ( cr )
    {
      rsd[0] = complex <double> ( 0.0, 0.0 );
    }
    return info;
  }
//
//  Set up to compute QY or QTY.
//
  if ( cqy )
  {
    for ( i = 0; i < n; i++ )
    {
      qy[i] = y[i];
    }
  }

  if ( cqty )
  {
    for ( i = 0; i < n; i++ )
    {
      qty[i] = y[i];
    }
  }
//
//  Compute QY.
//
  if ( cqy )
  {
    for ( jj = 1; jj <= ju; jj++ )
    {
      j = ju - jj + 1;

      if ( zabs1 ( qraux[j-1] ) != 0.0 )
      {
        temp = x[j-1+(j-1)*ldx];
        x[j-1+(j-1)*ldx] = qraux[j-1];
        t = -zdotc ( n-j+1, x+j-1+(j-1)*ldx, 1, qy+j-1, 1 ) / x[j-1+(j-1)*ldx];
        zaxpy ( n-j+1, t, x+j-1+(j-1)*ldx, 1, qy+j-1, 1 );
        x[j-1+(j-1)*ldx] = temp;
      }
    }
  }
//
//  Compute hermitian ( A ) * Y.
//
  if ( cqty )
  {
    for ( j = 1; j <= ju; j++ )
    {
      if ( zabs1 ( qraux[j-1] ) != 0.0 )
      {
        temp = x[j-1+(j-1)*ldx];
        x[j-1+(j-1)*ldx] = qraux[j-1];
        t = -zdotc ( n-j+1, x+j-1+(j-1)*ldx, 1, qty+j-1, 1 ) / x[j-1+(j-1)*ldx];
        zaxpy ( n-j+1, t, x+j-1+(j-1)*ldx, 1, qty+j-1, 1 );
        x[j-1+(j-1)*ldx] = temp;
      }
    }
  }
//
//  Set up to compute B, RSD, or XB.
//
  if ( cb )
  {
    for ( i = 0; i < k; i++ )
    {
      b[i] = qty[i];
    }
  }

  kp1 = k + 1;

  if ( cxb )
  {
    for ( i = 0; i < k; i++ )
    {
      xb[i] = qty[i];
    }
  }

  if ( cr && k < n )
  {
    for ( i = k; i < n; i++ )
    {
      rsd[i] = qty[i];
    }
  }

  if ( cxb )
  {
    for ( i = k; i < n; i++ )
    {
      xb[i] = complex <double> ( 0.0, 0.0 );
    }
  }

  if ( cr )
  {
    for ( i = 0; i < k; i++ )
    {
      rsd[i] = complex <double> ( 0.0, 0.0 );
    }
  }
//
//  Compute B.
//
  if ( cb )
  {
    for ( jj = 1; jj <= k; jj++ )
    {
      j = k - jj + 1;

      if ( zabs1 ( x[j-1+(j-1)*ldx] ) == 0.0 )
      {
        info = j;
        break;
      }

      b[j-1] = b[j-1] / x[j-1+(j-1)*ldx];

      if ( j != 1 )
      {
        t = -b[j-1];
        zaxpy ( j-1, t, x+0+(j-1)*ldx, 1, b, 1 );
      }
    }
  }

  if ( cr || cxb )
  {
//
//  Compute RSD or XB as required.
//
    for ( jj = 1; jj <= ju; jj++ )
    {
      j = ju - jj + 1;

      if ( zabs1 ( qraux[j-1] ) != 0.0 )
      {
        temp = x[j-1+(j-1)*ldx];
        x[j-1+(j-1)*ldx] = qraux[j-1];

        if ( cr )
        {
          t = -zdotc ( n-j+1, x+j-1+(j-1)*ldx, 1, rsd+j-1, 1 )
            / x[j-1+(j-1)*ldx];
          zaxpy ( n-j+1, t,x+j-1+(j-1)*ldx, 1, rsd+j-1, 1 );
        }

        if ( cxb )
        {
          t = -zdotc ( n-j+1, x+j-1+(j-1)*ldx, 1, xb+j-1, 1 )
            / x[j-1+(j-1)*ldx];
          zaxpy ( n-j+1, t, x+j-1+(j-1)*ldx, 1, xb+j-1, 1 );
        }
        x[j-1+(j-1)*ldx] = temp;
      }
    }
  }
  return info;
}
//****************************************************************************80

double zsico ( complex <double> a[], int lda, int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZSICO factors a complex symmetric matrix.
//
//  Discussion:
//
//    The factorization is done by symmetric pivoting.
//
//    The routine also estimates the condition of the matrix.
//
//    If RCOND is not needed, ZSIFA is slightly faster.
//
//    To solve A*X = B, follow ZSICO by ZSISL.
//
//    To compute inverse(A)*C, follow ZSICO by ZSISL.
//
//    To compute inverse(A), follow ZSICO by ZSIDI.
//
//    To compute determinant(A), follow ZSICO by ZSIDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the symmetric matrix to be
//    factored.  On output, a block diagonal matrix and the multipliers which
//    were used to obtain it.  The factorization can be written A = U*D*U'
//    where U is a product of permutation and unit upper triangular matrices,
//    U' is the transpose of U, and D is block diagonal with 1 by 1 and
//    2 by 2 blocks.  Only the diagonal and upper triangle are used.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, double ZSICO, an estimate of RCOND, the reciprocal condition of
//    the matrix.  For the system A*X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    (EPSILON/RCOND).  If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
//  Local Parameters:
//
//    Workspace, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
{
  complex <double> ak;
  complex <double> akm1;
  double anorm;
  complex <double> bk;
  complex <double> bkm1;
  complex <double> denom;
  complex <double> ek;
  int i;
  int info;
  int j;
  int k;
  int kp;
  int kps;
  int ks;
  double rcond;
  double s;
  double ynorm;
  complex <double> t;
  complex <double> *z;

  z = new complex <double> [n];
//
//  Find norm of A using only upper half.
//
  for ( j = 1; j <= n; j++ )
  {
    z[j-1] = complex <double> ( dzasum ( j, a+0+(j-1)*lda, 1 ), 0.0 );
    for ( i = 1; i <= j-1; i++ )
    {
      z[i-1] =
        complex <double> ( real ( z[i-1] ) + zabs1 ( a[i-1+(j-1)*lda] ), 0.0 );
    }
  }
  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    anorm = r8_max ( anorm, real ( z[j] ) );
  }
//
//  Factor.
//
  info = zsifa ( a, lda, n, ipvt );
//
//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where U*D*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve U*D*W = E.
//
  ek = complex <double> ( 1.0, 0.0 );
  for ( j = 0; j < n; j++ )
  {
    z[j] = complex <double> ( 0.0, 0.0 );
  }

  k = n;

  while ( 0 < k )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    kp = abs ( ipvt[k-1] );
    kps = k + 1 - ks;

    if ( kp != kps )
    {
      t        = z[kps-1];
      z[kps-1] = z[kp-1];
      z[kp-1]  = t;
    }

    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, z[k-1] );
    }

    z[k-1] = z[k-1] + ek;
    zaxpy ( k-ks, z[k-1], a+0+(k-1)*lda, 1, z, 1 );

    if ( ks != 1 )
    {
      if ( zabs1 ( z[k-2] ) != 0.0 )
      {
        ek = zsign1 ( ek, z[k-2] );
      }
      z[k-2] = z[k-2] + ek;
      zaxpy ( k-ks, z[k-2], a+0+(k-2)*lda, 1, z, 1 );
    }
    if ( ks != 2 )
    {
      if ( zabs1 ( a[k-1+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
      {
        s = zabs1 ( a[k-1+(k-1)*lda] ) / zabs1 ( z[k-1] );
        zdscal ( n, s, z, 1 );
        ek = complex <double> ( s, 0.0 ) * ek;
      }

      if ( zabs1 ( a[k-1+(k-1)*lda] ) != 0.0 )
      {
        z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
      }
      else
      {
        z[k-1] = complex <double> ( 1.0, 0.0 );
      }
    }
    else
    {
      ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
      akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
      bk = z[k-1] / a[k-2+(k-1)*lda];
      bkm1 = z[k-2] / a[k-2+(k-1)*lda];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      z[k-1] = ( akm1 * bk - bkm1 ) / denom;
      z[k-2] = ( ak * bkm1 - bk ) / denom;
    }
    k = k - ks;
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
//
//  Solve U' * Y = W.
//
  k = 1;
  while ( k <= n )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }
    if ( k != 1 )
    {
      z[k-1] = z[k-1] + zdotu ( k-1, a+0+(k-1)*lda, 1, z, 1 );

      if ( ks == 2 )
      {
        z[k] = z[k] + zdotu ( k-1, a+0+k*lda, 1, z, 1 );
      }
      kp = abs ( ipvt[k-1] );

      if ( kp != k )
      {
        t       = z[k-1];
        z[k-1]  = z[kp-1];
        z[kp-1] = t;
      }
    }
    k = k + ks;
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = 1.0;
//
//  Solve U*D*V = Y.
//
  k = n;

  while ( 0 < k )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    if ( k != ks )
    {
      kp = abs ( ipvt[k-1] );
      kps = k + 1 - ks;

      if ( kp != kps )
      {
        t        = z[kps-1];
        z[kps-1] = z[kp-1];
        z[kp-1]  = t;
      }
      zaxpy ( k-ks, z[k-1], a+0+(k-1)*lda, 1, z, 1 );

      if ( ks == 2 )
      {
        zaxpy ( k-ks, z[k-2], a+0+(k-2)*lda, 1, z, 1 );
      }
    }

    if ( ks != 2 )
    {
      if ( zabs1 ( a[k-1+(k-1)*lda] ) < zabs1 ( z[k-1] ) )
      {
        s = zabs1 ( a[k-1+(k-1)*lda] ) / zabs1 ( z[k-1] );
        zdscal ( n, s, z, 1 );
        ynorm = s * ynorm;
      }

      if ( zabs1 ( a[k-1+(k-1)*lda] ) != 0.0 )
      {
        z[k-1] = z[k-1] / a[k-1+(k-1)*lda];
      }
      else
      {
        z[k-1] = complex <double> ( 1.0, 0.0 );
      }
    }
    else
    {
      ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
      akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
      bk = z[k-1] / a[k-2+(k-1)*lda];
      bkm1 = z[k-2] / a[k-2+(k-1)*lda];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      z[k-1] = ( akm1 * bk - bkm1 ) / denom;
      z[k-2] = ( ak * bkm1 - bk ) / denom;
    }
    k = k - ks;
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;
//
//  Solve U' * Z = V.
//
  k = 1;

  while ( k <= n )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }
    if ( k != 1 )
    {
      z[k-1] = z[k-1] + zdotu ( k-1, a+0+(k-1)*lda, 1, z, 1 );

      if ( ks == 2 )
      {
        z[k] = z[k] + zdotu ( k-1, a+0+k*lda, 1, z, 1 );
      }
      kp = abs ( ipvt[k-1] );

      if ( kp != k )
      {
        t       = z[k-1];
        z[k-1]  = z[kp-1];
        z[kp-1] = t;
      }
    }
    k = k + ks;
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }

  delete [] z;
  return rcond;
}
//****************************************************************************80

void zsidi ( complex <double> a[], int lda, int n, int ipvt[],
  complex <double> det[2], int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZSIDI computes the determinant and inverse of a matrix factored by ZSIFA.
//
//  Discussion:
//
//    It is assumed the complex symmetric matrix has already been factored
//    by ZSIFA.
//
//    A division by zero may occur if the inverse is requested
//    and ZSICO set RCOND == 0.0 or ZSIFA set INFO nonzero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the output from ZSIFA.
//    If the inverse was requested, then on output, A contains the upper triangle
//    of the inverse of the original matrix.  The strict lower triangle
//    is never referenced.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZSIFA.
//
//    Output, complex <double> DET[2], if requested, the determinant of the matrix.
//    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= abs ( DET(1) ) < 10.0
//    or DET(1) = 0.0.  Also, DET(2) is strictly real.
//
//    Input, int JOB, has the decimal expansion AB where
//    if B != 0, the inverse is computed,
//    if A != 0, the determinant is computed,
//    For example, JOB = 11 gives both.
//
{
  complex <double> ak;
  complex <double> akkp1;
  complex <double> akp1;
  complex <double> d;
  int i;
  int j;
  int jb;
  int k;
  int km1;
  int ks;
  int kstep;
  bool nodet;
  bool noinv;
  complex <double> t;
  complex <double> *work;

  noinv = ( job % 10 ) == 0;
  nodet = ( job % 100 ) / 10 == 0;

  if ( !nodet )
  {
    det[0] = complex <double> ( 1.0, 0.0 );
    det[1] = complex <double> ( 0.0, 0.0 );
    t = complex <double> ( 0.0, 0.0 );

    for ( k = 1; k <= n; k++ )
    {
      d = a[k-1+(k-1)*lda];
//
//   2 by 2 block.
//   Use det ( D  T ) = ( D / T * C - T ) * T
//           ( T  C )
//   to avoid underflow/overflow troubles.
//   Take two passes through scaling.  Use T for flag.
//
      if ( ipvt[k-1] <= 0 )
      {
        if ( zabs1 ( t ) == 0.0 )
        {
          t = a[k-1+k*lda];
          d = ( d / t ) * a[k+k*lda] - t;
        }
        else
        {
          d = t;
          t = complex <double> ( 0.0, 0.0 );
        }
      }
      det[0] = det[0] * d;

      if ( zabs1 ( det[0] ) != 0.0 )
      {
        while ( zabs1 ( det[0] ) < 1.0 )
        {
          det[0] = det[0] * complex <double> ( 10.0, 0.0 );
          det[1] = det[1] - complex <double> ( 1.0, 0.0 );
        }
        while ( 10.0 <= zabs1 ( det[0] ) )
        {
          det[0] = det[0] / complex <double> ( 10.0, 0.0 );
          det[1] = det[1] + complex <double> ( 1.0, 0.0 );
        }
      }
    }
  }
//
//  Compute inverse ( A ).
//
  if ( !noinv )
  {
    work = new complex <double> [n];

    k = 1;

    while ( k <= n )
    {
      km1 = k - 1;
//
//  1 by 1
//
      if ( 0 <= ipvt[k-1] )
      {
        a[k-1+(k-1)*lda] = complex <double> ( 1.0, 0.0 ) / a[k-1+(k-1)*lda];

        if ( 1 <= km1 )
        {
          for ( i = 1; i <= km1; i++ )
          {
            work[i-1] = a[i-1+(k-1)*lda];
          }

          for ( j = 1; j <= km1; j++ )
          {
            a[j-1+(k-1)*lda] = zdotu ( j, a+0+(j-1)*lda, 1, work, 1 );
            zaxpy ( j-1, work[j-1], a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
          }
          a[k-1+(k-1)*lda] = a[k-1+(k-1)*lda]
            + zdotu ( km1, work, 1, a+0+(k-1)*lda, 1 );
        }
        kstep = 1;
      }
//
//  2 by 2
//
      else
      {
        t = a[k-1+k*lda];
        ak = a[k-1+(k-1)*lda] / t;
        akp1 = a[k+k*lda] / t;
        akkp1 = a[k-1+k*lda] / t;
        d = t * ( ak * akp1 - complex <double> ( 1.0, 0.0 ) );
        a[k-1+(k-1)*lda] = akp1 / d;
        a[k+k*lda] = ak / d;
        a[k-1+k*lda] = -akkp1 / d;

        if ( 1 <= km1 )
        {
          for ( i = 1; i <= km1; i++ )
          {
            work[i-1] = a[i-1+k*lda];
          }

          for ( j = 1; j <= km1; j++ )
          {
            a[j-1+k*lda] = zdotu ( j, a+0+(j-1)*lda, 1, work, 1 );
            zaxpy ( j-1, work[j-1], a+0+(j-1)*lda, 1, a+0+k*lda, 1 );
          }

          a[k+k*lda] = a[k+k*lda] + zdotu ( km1, work, 1, a+0+k*lda, 1 );
          a[k-1+k*lda] = a[k-1+k*lda]
            + zdotu ( km1, a+0+(k-1)*lda, 1, a+0+k*lda, 1 );

          for ( i = 1; i <= km1; i++ )
          {
            work[i-1] = a[i-1+(k-1)*lda];
          }

          for ( j = 1; j <= km1; j++ )
          {
            a[j-1+(k-1)*lda] = zdotu ( j, a+0+(j-1)*lda, 1, work, 1 );
            zaxpy ( j-1, work[j-1], a+0+(j-1)*lda, 1, a+0+(k-1)*lda, 1 );
          }
          a[k-1+(k-1)*lda] = a[k-1+(k-1)*lda]
            + zdotu ( km1, work, 1, a+0+(k-1)*lda, 1 );
        }
        kstep = 2;
      }
//
//  Swap.
//
      ks = abs ( ipvt[k-1] );

      if ( ks != k )
      {
        zswap ( ks, a+0+(ks-1)*lda, 1, a+0+(k-1)*lda, 1 );

        for ( jb = ks; jb <= k; jb++ )
        {
          j = k + ks - jb;

          t                 = a[j-1+(k-1)*lda];
          a[j-1+(k-1)*lda]  = a[ks-1+(j-1)*lda];
          a[ks-1+(j-1)*lda] = t;
        }

        if ( kstep != 1 )
        {
          t             = a[ks-1+k*lda];
          a[ks-1+k*lda] = a[k-1+k*lda];
          a[k-1+k*lda]  = t;
        }
      }
      k = k + kstep;
    }
    delete [] work;
  }
  return;
}
//****************************************************************************80

int zsifa ( complex <double> a[], int lda, int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZSIFA factors a complex symmetric matrix.
//
//  Discussion:
//
//    The factorization is accomplished by elimination with symmetric pivoting.
//
//    To solve A*X = B, follow ZSIFA by ZSISL.
//
//    To compute inverse(A)*C, follow ZSIFA by ZSISL.
//
//    To compute determinant(A), follow ZSIFA by ZSIDI.
//
//    To compute inverse(A), follow ZSIFA by ZSIDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> A[LDA*N]; on input, the symmetric matrix to be
//    factored.  On output, a block diagonal matrix and the multipliers which
//    were used to obtain it.  The factorization can be written A = U*D*U'
//    where U is a product of permutation and unit upper triangular matrices,
//    U' is the transpose of U, and D is block diagonal with 1 by 1 and 2 by 2
//    blocks.  Only the diagonal and upper triangle are used.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int ZSIFA.
//    0, normal value.
//    K, if the K-th pivot block is singular.  This is not an error condition
//    for this subroutine, but it does indicate that ZSISL or ZSIDI may
//    divide by zero if called.
//
{
  double absakk;
  complex <double> ak;
  complex <double> akm1;
  double alpha;
  complex <double> bk;
  complex <double> bkm1;
  double colmax;
  complex <double> denom;
  int imax;
  int info;
  int j;
  int jj;
  int jmax;
  int k;
  int km1;
  int km2;
  int kstep;
  complex <double> mulk;
  complex <double> mulkm1;
  double rowmax;
  bool swap;
  complex <double> t;
//
//  Initialize.
//
//  ALPHA is used in choosing pivot block size.
//
  alpha = ( 1.0 + sqrt ( 17.0 ) ) / 8.0;

  info = 0;
//
//  Main loop on K, which goes from N to 1.
//
  k = n;

  for ( ; ; )
  {
//
//  Leave the loop if K = 0 or K = 1.
//
    if ( k == 0 )
    {
      break;
    }

    if ( k == 1 )
    {
      ipvt[0] = 1;
      if ( zabs1 ( a[0+0*lda] ) == 0.0 )
      {
        info = 1;
      }
      break;
    }
//
//  This section of code determines the kind of
//  elimination to be performed.  When it is completed,
//  KSTEP will be set to the size of the pivot block, and
//  SWAP will be set to TRUE if an interchange is
//  required.
//
    km1 = k - 1;
    absakk = zabs1 ( a[k-1+(k-1)*lda] );
//
//  Determine the largest off-diagonal element in column K.
//
    imax = izamax ( k-1, a+0+(k-1)*lda, 1 );
    colmax = zabs1 ( a[imax-1+(k-1)*lda] );

    if ( alpha * colmax < absakk )
    {
      kstep = 1;
      swap = false;
    }
//
//  Determine the largest off-diagonal element in row IMAX.
//
    else
    {
      rowmax = 0.0;

      for ( j = imax + 1; j <= k; j++ )
      {
        rowmax = r8_max ( rowmax, zabs1 ( a[imax-1+(j-1)*lda] ) );
      }

      if ( imax != 1 )
      {
        jmax = izamax ( imax-1, a+0+(imax-1)*lda, 1 );
        rowmax = r8_max ( rowmax, zabs1 ( a[jmax-1+(imax-1)*lda] ) );
      }

      if ( alpha * rowmax <= zabs1 ( a[imax-1+(imax-1)*lda] ) )
      {
        kstep = 1;
        swap = true;
      }
      else if ( alpha * colmax * ( colmax / rowmax ) <= absakk )
      {
        kstep = 1;
        swap = false;
      }
      else
      {
        kstep = 2;
        swap = ( imax != km1 );
      }
    }
//
//  Column K is zero.  Set INFO and iterate the loop.
//
    if ( r8_max ( absakk, colmax ) == 0.0 )
    {
      ipvt[k-1] = k;
      info = k;
      k = k - kstep;
      continue;
    }

    if ( kstep != 2 )
    {
//
//  1 x 1 pivot block.
//
//  Perform an interchange.
//
      if ( swap )
      {
        zswap ( imax, a+0+(imax-1)*lda, 1, a+0+(k-1)*lda, 1 );

        for ( jj = imax; jj <= k; jj++ )
        {
          j = k + imax - jj;

          t                   = a[j-1+(k-1)*lda];
          a[j-1+(k-1)*lda]    = a[imax-1+(j-1)*lda];
          a[imax-1+(j-1)*lda] = t;
        }
      }
//
//  Perform the elimination.
//
      for ( jj = 1; jj <= km1; jj++ )
      {
        j = k - jj;
        mulk = -a[j-1+(k-1)*lda] / a[k-1+(k-1)*lda];
        t = mulk;
        zaxpy ( j, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
        a[j-1+(k-1)*lda] = mulk;
      }
//
//  Set the pivot array.
//
      if ( swap )
      {
        ipvt[k-1] = imax;
      }
      else
      {
        ipvt[k-1] = k;
      }
    }
//
//  2 x 2 pivot block.
//
    else
    {
      if ( swap )
      {
        zswap ( imax, a+0+(imax-1)*lda, 1, a+0+(k-2)*lda, 1 );

        for ( jj = imax; jj <= km1; jj++ )
        {
          j = km1 + imax - jj;

          t                   = a[j-1+(k-2)*lda];
          a[j-1+(k-2)*lda]    = a[imax-1+(j-1)*lda];
          a[imax-1+(j-1)*lda] = t;
        }
        t                   = a[k-2+(k-1)*lda];
        a[k-2+(k-1)*lda]    = a[imax-1+(k-1)*lda];
        a[imax-1+(k-1)*lda] = t;
      }
//
//  Perform the elimination.
//
      km2 = k - 2;

      if ( km2 != 0 )
      {
        ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
        akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
        denom = complex <double> ( 1.0, 0.0 ) - ak * akm1;

        for ( jj = 1; jj <= km2; jj++ )
        {
          j = km1 - jj;
          bk = a[j-1+(k-1)*lda] / a[k-2+(k-1)*lda];
          bkm1 = a[j-1+(k-2)*lda] / a[k-2+(k-1)*lda];
          mulk = ( akm1 * bk - bkm1 ) / denom;
          mulkm1 = ( ak * bkm1 - bk ) / denom;
          t = mulk;
          zaxpy ( j, t, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
          t = mulkm1;
          zaxpy ( j, t, a+0+(k-2)*lda, 1, a+0+(j-1)*lda, 1 );
          a[j-1+(k-1)*lda] = mulk;
          a[j-1+(k-2)*lda] = mulkm1;
        }
      }
//
//  Set the pivot array.
//
      if ( swap )
      {
        ipvt[k-1] = -imax;
      }
      else
      {
        ipvt[k-1] = 1 - k;
      }
      ipvt[k-2] = ipvt[k-1];
    }
    k = k - kstep;
  }
  return info;
}
//****************************************************************************80

void zsisl ( complex <double> a[], int lda, int n, int ipvt[],
  complex <double> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZSISL solves a complex symmetric system that was factored by ZSIFA.
//
//  Discussion:
//
//    A division by zero may occur if ZSICO has set RCOND == 0.0
//    or ZSIFA has set INFO != 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> A[LDA*N], the output from ZSICO or ZSIFA.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZSICO or ZSIFA.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  complex <double> ak;
  complex <double> akm1;
  complex <double> bk;
  complex <double> bkm1;
  complex <double> denom;
  int k;
  int kp;
  complex <double> t;
//
//  Loop backward applying the transformations and D inverse to B.
//
  k = n;

  while ( 0 < k )
  {
//
//  1 x 1 pivot block.
//
    if ( 0 <= ipvt[k-1] )
    {
      if ( k != 1 )
      {
        kp = ipvt[k-1];

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
        zaxpy ( k-1, b[k-1], a+0+(k-1)*lda, 1, b, 1 );
      }

      b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
      k = k - 1;
    }
//
//  2 x 2 pivot block.
//
    else
    {
      if ( k != 2 )
      {
        kp = abs ( ipvt[k-1] );

        if ( kp != k - 1 )
        {
          t       = b[k-2];
          b[k-2]  = b[kp-1];
          b[kp-1] = t;
        }
        zaxpy ( k-2, b[k-1], a+0+(k-1)*lda, 1, b, 1 );
        zaxpy ( k-2, b[k-2], a+0+(k-2)*lda, 1, b, 1 );
      }
      ak = a[k-1+(k-1)*lda] / a[k-2+(k-1)*lda];
      akm1 = a[k-2+(k-2)*lda] / a[k-2+(k-1)*lda];
      bk = b[k-1] / a[k-2+(k-1)*lda];
      bkm1 = b[k-2] / a[k-2+(k-1)*lda];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      b[k-1] = ( akm1 * bk - bkm1 ) / denom;
      b[k-2] = ( ak * bkm1 - bk ) / denom;
      k = k - 2;
    }
  }
//
//  Loop forward applying the transformations.
//
  k = 1;

  while ( k <= n )
  {
    if ( 0 <= ipvt[k-1] )
    {
//
//  1 x 1 pivot block.
//
      if ( k != 1 )
      {
        b[k-1] = b[k-1] + zdotu ( k-1, a+0+(k-1)*lda, 1, b, 1 );
        kp = ipvt[k-1];

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
      }
      k = k + 1;
    }
//
//  2 x 2 pivot block.
//
    else
    {
      if ( k != 1 )
      {
        b[k-1] = b[k-1] + zdotu ( k-1, a+0+(k-1)*lda, 1, b, 1 );
        b[k] = b[k] + zdotu ( k-1, a+0+k*lda, 1, b, 1 );
        kp = abs ( ipvt[k-1] );

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
      }
      k = k + 2;
    }
  }
  return;
}
//****************************************************************************80

double zspco ( complex <double> ap[], int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZSPCO factors a complex <double> symmetric matrix stored in packed form.
//
//  Discussion:
//
//    The routine also estimates the condition of the matrix.
//
//    If RCOND is not needed, ZSPFA is slightly faster.
//
//    To solve A*X = B, follow ZSPCO by ZSPSL.
//
//    To compute inverse(A)*C, follow ZSPCO by ZSPSL.
//
//    To compute inverse(A), follow ZSPCO by ZSPDI.
//
//    To compute determinant(A), follow ZSPCO by ZSPDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> AP[N*(N+1)/2]; on input, the packed form of a
//    symmetric matrix.  The columns of the upper triangle are stored
//    sequentially in a one-dimensional array.  On output, a block diagonal
//    matrix and the multipliers which were used to obtain it, stored in packed
//    form.  The factorization can be written A = U*D*U' where U is a product
//    of permutation and unit upper triangular matrices, U' is the transpose
//    of U, and D is block diagonal with 1 by 1 and 2 by 2 blocks.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, double ZSPCO, an estimate of RCOND, the reciprocal condition of
//    the matrix.  For the system A*X = B, relative perturbations in A and B
//    of size EPSILON may cause relative perturbations in X of size
//    (EPSILON/RCOND).  If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then A may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate underflows.
//
//  Local Parameters:
//
//    Local, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If A is close to a singular matrix, then Z is an
//    approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
{
  complex <double> ak;
  complex <double> akm1;
  double anorm;
  complex <double> bk;
  complex <double> bkm1;
  complex <double> denom;
  complex <double> ek;
  int i;
  int ij;
  int ik;
  int ikm1;
  int ikp1;
  int info;
  int j;
  int j1;
  int k;
  int kk;
  int km1k;
  int km1km1;
  int kp;
  int kps;
  int ks;
  double rcond;
  double s;
  complex <double> t;
  double ynorm;
  complex <double> *z;

  z = new complex <double> [n];
//
//  Find norm of A using only upper half.
//
  j1 = 1;

  for ( j = 1; j <= n; j++ )
  {
    z[j-1] = complex <double> ( dzasum ( j, ap+j1-1, 1 ), 0.0 );
    ij = j1;
    j1 = j1 + j;

    for ( i = 1; i <= j-1; i++ )
    {
      z[i-1] = complex <double> ( real ( z[i-1] ) + zabs1 ( ap[ij-1] ), 0.0 );
      ij = ij + 1;
    }
  }

  anorm = 0.0;
  for ( j = 0; j < n; j++ )
  {
    anorm = r8_max ( anorm, real ( z[j] ) );
  }
//
//  Factor.
//
  info = zspfa ( ap, n, ipvt );
//
//  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
//
//  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of W where U*D*W = E.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve U*D*W = E.
//
  ek = complex <double> ( 1.0, 0.0 );
  for ( j = 0; j < n; j++ )
  {
    z[j] = complex <double> ( 0.0, 0.0 );
  }

  k = n;
  ik = ( n * ( n - 1 ) ) / 2;

  while ( 0 < k )
  {
    kk = ik + k;
    ikm1 = ik - ( k - 1 );

    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    kp = abs ( ipvt[k-1] );
    kps = k + 1 - ks;

    if ( kp != kps )
    {
      t       = z[kps-1];
      z[kps-1]  = z[kp-1];
      z[kp-1] = t;
    }

    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, z[k-1] );
    }
    z[k-1] = z[k-1] + ek;
    zaxpy ( k-ks, z[k-1], ap+ik, 1, z, 1 );

    if ( ks != 1 )
    {
      if ( zabs1 ( z[k-2] ) != 0.0 )
      {
        ek = zsign1 ( ek, z[k-2] );
      }
      z[k-2] = z[k-2] + ek;
      zaxpy ( k-ks, z[k-2], ap+ikm1, 1, z, 1 );
    }

    if ( ks != 2 )
    {
      if ( zabs1 ( ap[kk-1] ) < zabs1 ( z[k-1] ) )
      {
        s = zabs1 ( ap[kk-1] ) / zabs1 ( z[k-1] );
        zdscal ( n, s, z, 1 );
        ek = complex <double> ( s, 0.0 ) * ek;
      }
      if ( zabs1 ( ap[kk-1] ) != 0.0 )
      {
        z[k-1] = z[k-1] / ap[kk-1];
      }
      else
      {
        z[k-1] = complex <double> ( 1.0, 0.0 );
      }
    }
    else
    {
      km1k = ik + k - 1;
      km1km1 = ikm1 + k - 1;
      ak = ap[kk-1] / ap[km1k-1];
      akm1 = ap[km1km1-1] / ap[km1k-1];
      bk = z[k-1] / ap[km1k-1];
      bkm1 = z[k-2] / ap[km1k-1];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      z[k-1] = ( akm1 * bk - bkm1 ) / denom;
      z[k-2] = ( ak * bkm1 - bk ) / denom;
    }
    k = k - ks;
    ik = ik - k;
    if ( ks == 2 )
    {
      ik = ik - ( k + 1 );
    }
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
//
//  Solve trans(U) * Y = W.
//
  k = 1;
  ik = 0;

  while ( k <= n )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }
    if ( k != 1 )
    {
      z[k-1] = z[k-1] + zdotu ( k-1, ap+ik, 1, z, 1 );
      ikp1 = ik + k;

      if ( ks == 2 )
      {
        z[k] = z[k] + zdotu ( k-1, ap+ikp1, 1, z, 1 );
      }

      kp = abs ( ipvt[k-1] );

      if ( kp != k )
      {
        t       = z[k-1];
        z[k-1]  = z[kp-1];
        z[kp-1] = t;
      }
    }
    ik = ik + k;
    if ( ks == 2 )
    {
      ik = ik + ( k + 1 );
    }
    k = k + ks;
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = 1.0;
//
//  Solve U*D*V = Y.
//
  k = n;
  ik = ( n * ( n - 1 ) ) / 2;

  while ( 0 < k )
  {
    kk = ik + k;
    ikm1 = ik - ( k - 1 );

    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    if ( k != ks )
    {
      kp = abs ( ipvt[k-1] );
      kps = k + 1 - ks;

      if ( kp != kps )
      {
        t        = z[kps-1];
        z[kps-1] = z[kp-1];
        z[kp-1]  = t;
      }

      zaxpy ( k-ks, z[k-1], ap+ik, 1, z, 1 );

      if ( ks == 2 )
      {
        zaxpy ( k-ks, z[k-2], ap+ikm1, 1, z, 1 );
      }
    }

    if ( ks != 2 )
    {
      if ( zabs1 ( ap[kk-1] ) < zabs1 ( z[k-1] ) )
      {
        s = zabs1 ( ap[kk-1] ) / zabs1 ( z[k-1] );
        zdscal ( n, s, z, 1 );
        ynorm = s * ynorm;
      }

      if ( zabs1 ( ap[kk-1] ) != 0.0 )
      {
        z[k-1] = z[k-1] / ap[kk-1];
      }
      else
      {
        z[k-1] = complex <double> ( 1.0, 0.0 );
      }
    }
    else
    {
      km1k = ik + k - 1;
      km1km1 = ikm1 + k - 1;
      ak = ap[kk-1] / ap[km1k-1];
      akm1 = ap[km1km1-1] / ap[km1k-1];
      bk = z[k-1] / ap[km1k-1];
      bkm1 = z[k-2] / ap[km1k-1];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      z[k-1] = ( akm1 * bk - bkm1 ) / denom;
      z[k-2] = ( ak * bkm1 - bk ) / denom;
    }

    k = k - ks;
    ik = ik - k;

    if ( ks == 2 )
    {
      ik = ik - ( k + 1 );
    }
  }
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;
//
//  Solve U' * Z = V.
//
  k = 1;
  ik = 0;

  while ( k <= n )
  {
    if ( ipvt[k-1] < 0 )
    {
      ks = 2;
    }
    else
    {
      ks = 1;
    }

    if ( k != 1 )
    {
      z[k-1] = z[k-1] + zdotu ( k-1, ap+ik, 1, z, 1 );
      ikp1 = ik + k;

      if ( ks == 2 )
      {
        z[k] = z[k] + zdotu ( k-1, ap+ikp1, 1, z, 1 );
      }

      kp = abs ( ipvt[k-1] );

      if ( kp != k )
      {
        t       = z[k-1];
        z[k-1]  = z[kp-1];
        z[kp-1] = t;
      }
    }

    ik = ik + k;

    if ( ks == 2 )
    {
      ik = ik + ( k + 1 );
    }
    k = k + ks;
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( anorm != 0.0 )
  {
    rcond = ynorm / anorm;
  }
  else
  {
    rcond = 0.0;
  }

  delete [] z;

  return rcond;
}
//****************************************************************************80

void zspdi ( complex <double> ap[], int n, int ipvt[], complex <double> det[2],
  int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZSPDI sets the determinant and inverse of a complex symmetric packed matrix.
//
//  Discussion:
//
//    ZSPDI uses the factors from ZSPFA.
//
//    The matrix is stored in packed form.
//
//    A division by zero will occur if the inverse is requested and ZSPCO has
//    set RCOND to 0.0 or ZSPFA has set INFO nonzero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> AP[N*(N+1)/2]; on input, the matrix factors
//    from ZSPFA.  On output, if the inverse was requested, the upper
//    triangle of the inverse of the original matrix, stored in packed
//    form.  The columns of the upper triangle are stored sequentially
//    in a one-dimensional array.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZSPFA.
//
//    Output, complex <double> DET[2], the determinant of the original matrix.
//    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= abs ( DET(1) ) < 10.0
//    or DET(1) = 0.0.  Also, DET(2) is strictly real.
//
//    Input, int JOB, has the decimal expansion AB where
//    if B != 0, the inverse is computed,
//    if A != 0, the determinant is computed,
//    For example, JOB = 11 gives both.
//
{
  complex <double> ak;
  complex <double> akkp1;
  complex <double> akp1;
  complex <double> d;
  int i;
  int ij;
  int ik;
  int ikp1;
  int iks;
  int j;
  int jb;
  int jk;
  int jkp1;
  int k;
  int kk;
  int kkp1;
  int km1;
  int ks;
  int ksj;
  int kskp1;
  int kstep;
  bool nodet;
  bool noinv;
  complex <double> t;
  complex <double> *work;

  noinv = ( job % 10 ) == 0;
  nodet = ( job % 100 ) / 10 == 0;

  if ( !nodet )
  {
    det[0] = complex <double> ( 1.0, 0.0 );
    det[1] = complex <double> ( 0.0, 0.0 );
    t = complex <double> ( 0.0, 0.0 );
    ik = 0;

    for ( k = 1; k <= n; k++ )
    {
      kk = ik + k;
      d = ap[kk-1];
//
//  2 by 2 block
//  Use det (D  T)  =  ( D / T * C - T ) * T
//          (T  C)
//  to avoid underflow/overflow troubles.
//  Take two passes through scaling.  Use T for flag.
//
      if ( ipvt[k-1] <= 0 )
      {
        if ( zabs1 ( t ) == 0.0 )
        {
          ikp1 = ik + k;
          kkp1 = ikp1 + k;
          t = ap[kkp1-1];
          d = ( d / t ) * ap[kkp1] - t;
        }
        else
        {
          d = t;
          t = complex <double> ( 0.0, 0.0 );
        }
      }

      if ( !nodet )
      {
        det[0] = det[0] * d;

        if ( zabs1 ( det[0] ) != 0.0 )
        {
          while ( zabs1 ( det[0] ) < 1.0 )
          {
            det[0] = det[0] * complex <double> ( 10.0, 0.0 );
            det[1] = det[1] - complex <double> ( 1.0, 0.0 );
          }

          while ( 10.0 <= zabs1 ( det[0] ) )
          {
            det[0] = det[0] / complex <double> ( 10.0, 0.0 );
            det[1] = det[1] + complex <double> ( 1.0, 0.0 );
          }
        }
      }
      ik = ik + k;
    }
  }
//
//  Compute inverse ( A ).
//
  if ( !noinv )
  {
    work = new complex <double> [n];
    k = 1;
    ik = 0;

    while ( k <= n )
    {
      km1 = k - 1;
      kk = ik + k;
      ikp1 = ik + k;

      if ( 0 <= ipvt[k-1] )
      {
//
//  1 by 1
//
        ap[kk-1] = complex <double> ( 1.0, 0.0 ) / ap[kk-1];

        if ( 1 <= km1 )
        {
          for ( i = 1; i <= km1; i++ )
          {
            work[i-1] = ap[ik+i-1];
          }
          ij = 0;

          for ( j = 1; j <= km1; j++ )
          {
            jk = ik + j;
            ap[jk-1] = zdotu ( j, ap+ij, 1, work, 1 );
            zaxpy ( j-1, work[j-1], ap+ij, 1, ap+ik, 1 );
            ij = ij + j;
          }
          ap[kk-1] = ap[kk-1] + zdotu ( km1, work, 1, ap+ik, 1 );
        }
        kstep = 1;
      }
//
//  2 by 2
//
      else
      {
        kkp1 = ikp1 + k;
        t = ap[kkp1-1];
        ak = ap[kk-1] / t;
        akp1 = ap[kkp1] / t;
        akkp1 = ap[kkp1-1] / t;
        d = t * ( ak * akp1 - complex <double> ( 1.0, 0.0 ) );
        ap[kk-1] = akp1 / d;
        ap[kkp1] = ak / d;
        ap[kkp1-1] = -akkp1 / d;

        if ( 1 <= km1 )
        {
          for ( i = 1; i <= km1; i++ )
          {
            work[i-1] = ap[ikp1-1];
          }
          ij = 0;

          for ( j = 1; j <= km1; j++ )
          {
            jkp1 = ikp1 + j;
            ap[jkp1-1] = zdotu ( j, ap+ij, 1, work, 1 );
            zaxpy ( j-1, work[j-1], ap+ij, 1, ap+ikp1, 1 );
            ij = ij + j;
          }

          ap[kkp1] = ap[kkp1] + zdotu ( km1, work, 1, ap+ikp1, 1 );
          ap[kkp1-1] = ap[kkp1-1] + zdotu ( km1, ap+ik, 1, ap+ikp1, 1 );

          for ( i = 1; i <= km1; i++ )
          {
            work[i-1] = ap[ik+i-1];
          }
          ij = 0;

          for ( j = 1; j <= km1; j++ )
          {
            jk = ik + j;
            ap[jk-1] = zdotu ( j, ap+ij, 1, work, 1 );
            zaxpy ( j-1, work[j-1], ap+ij, 1, ap+ik, 1 );
            ij = ij + j;
          }
          ap[kk-1] = ap[kk-1] + zdotu ( km1, work, 1, ap+ik, 1 );
        }
        kstep = 2;
      }
//
//  Swap.
//
      ks = abs ( ipvt[k-1] );

      if ( ks != k )
      {
        iks = ( ks * ( ks - 1 ) ) / 2;
        zswap ( ks, ap+iks, 1, ap+ik, 1 );
        ksj = ik + ks;

        for ( jb = ks; jb <= k; jb++ )
        {
          j = k + ks - jb;
          jk = ik + j;

          t         = ap[jk-1];
          ap[jk-1]  = ap[ksj-1];
          ap[ksj-1] = t;

          ksj = ksj - ( j - 1 );
        }

        if ( kstep != 1 )
        {
          kskp1 = ikp1 + ks;

          t           = ap[kskp1-1];
          ap[kskp1-1] = ap[kkp1-1];
          ap[kkp1-1]  = t;
        }
      }
      ik = ik + k;

      if ( kstep == 2 )
      {
        ik = ik + k + 1;
      }
      k = k + kstep;
    }
    delete [] work;
  }
  return;
}
//****************************************************************************80

int zspfa ( complex <double> ap[], int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZSPFA factors a complex symmetric matrix stored in packed form.
//
//  Discussion:
//
//    The factorization is done by elimination with symmetric pivoting.
//
//    To solve A*X = B, follow ZSPFA by ZSPSL.
//
//    To compute inverse(A)*C, follow ZSPFA by ZSPSL.
//
//    To compute determinant(A), follow ZSPFA by ZSPDI.
//
//    To compute inverse(A), follow ZSPFA by ZSPDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> AP[N*(N+1)/2]; On input, the packed form of a
//    symmetric matrix A.  The columns of the upper triangle are stored
//    sequentially in a one-dimensional array.  On output, a block diagonal
//    matrix and the multipliers which were used to obtain it stored in
//    packed form.  The factorization can be written A = U*D*U' where U
//    is a product of permutation and unit upper triangular matrices,
//    U' is the transpose of U, and D is block diagonal with 1 by 1 and
//    2 by 2 blocks.
//
//    Input, int N, the order of the matrix.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int ZSPFA.
//    0, normal value.
//    K, if the K-th pivot block is singular.  This is not an error condition
//    for this subroutine, but it does indicate that ZSPSL or ZSPDI may
//    divide by zero if called.
//
{
  double absakk;
  complex <double> ak;
  complex <double> akm1;
  double alpha;
  complex <double> bk;
  complex <double> bkm1;
  double colmax;
  complex <double> denom;
  int ij;
  int ijj;
  int ik;
  int ikm1;
  int im;
  int imax;
  int imim;
  int imj;
  int imk;
  int info;
  int j;
  int jj;
  int jk;
  int jkm1;
  int jmax;
  int jmim;
  int k;
  int kk;
  int km1;
  int km1k;
  int km1km1;
  int km2;
  int kstep;
  complex <double> mulk;
  complex <double> mulkm1;
  double rowmax;
  bool swap;
  complex <double> t;
//
//  Initialize.
//
//  ALPHA is used in choosing pivot block size.
//
  alpha = ( 1.0 + sqrt ( 17.0 ) ) / 8.0;

  info = 0;
//
//  Main loop on K, which goes from N to 1.
//
  k = n;
  ik = ( n * ( n - 1 ) ) / 2;

  for ( ; ; )
  {
//
//  Leave the loop if K = 0 or K = 1.
//
    if ( k == 0 )
    {
      break;
    }

    if ( k == 1 )
    {
      ipvt[0] = 1;
      if ( zabs1 ( ap[0] ) == 0.0 )
      {
        info = 1;
      }
      break;
    }
//
//  This section of code determines the kind of
//  elimination to be performed.  When it is completed,
//  KSTEP will be set to the size of the pivot block, and
//  SWAP will be set to .true. if an interchange is
//  required.
//
    km1 = k - 1;
    kk = ik + k;
    absakk = zabs1 ( ap[kk-1] );
//
//  Determine the largest off-diagonal element in column K.
//
    imax = izamax ( k-1, ap+ik, 1 );
    imk = ik + imax;
    colmax = zabs1 ( ap[imk-1] );

    if ( alpha * colmax <= absakk )
    {
      kstep = 1;
      swap = false;
    }
//
//  Determine the largest off-diagonal element in row IMAX.
//
    else
    {
      rowmax = 0.0;
      im = ( imax * ( imax - 1 ) ) / 2;
      imj = im + 2 * imax;

      for ( j = imax + 1; j <= k; j++ )
      {
        rowmax = r8_max ( rowmax, zabs1 ( ap[imj-1] ) );
        imj = imj + j;
      }

      if ( imax != 1 )
      {
        jmax = izamax ( imax-1, ap+im, 1 );
        jmim = jmax + im;
        rowmax = r8_max ( rowmax, zabs1 ( ap[jmim-1] ) );
      }

      imim = imax + im;

      if ( alpha * rowmax <= zabs1 ( ap[imim-1] ) )
      {
        kstep = 1;
        swap = true;
      }
      else if ( alpha * colmax * ( colmax / rowmax ) <= absakk )
      {
        kstep = 1;
        swap = false;
      }
      else
      {
        kstep = 2;
        swap = ( imax != km1 );
      }
    }
//
//  Column K is zero.  Set INFO and iterate the loop.
//
    if ( r8_max ( absakk, colmax ) == 0.0 )
    {
      ipvt[k-1] = k;
      info = k;
      ik = ik - ( k - 1 );
      if ( kstep == 2 )
      {
        ik = ik - ( k - 2 );
      }
      k = k - kstep;
      continue;
    }

    if ( kstep != 2 )
    {
//
//  1 x 1 pivot block.
//
      if ( swap )
      {
        zswap ( imax, ap+im, 1, ap+ik, 1 );
        imj = ik + imax;

        for ( jj = imax; jj <= k; jj++ )
        {
          j = k + imax - jj;
          jk = ik + j;

          t         = ap[jk-1];
          ap[jk-1]  = ap[imj-1];
          ap[imj-1] = t;

          imj = imj - ( j - 1 );
        }
      }
//
//  Perform the elimination.
//
      ij = ik - ( k - 1 );

      for ( jj = 1; jj <= km1; jj++ )
      {
        j = k - jj;
        jk = ik + j;
        mulk = -ap[jk-1] / ap[kk-1];
        t = mulk;
        zaxpy ( j, t, ap+ik, 1, ap+ij, 1 );
        ijj = ij + j;
        ap[jk-1] = mulk;
        ij = ij - ( j - 1 );
      }
//
//  Set the pivot array.
//
      if ( swap )
      {
        ipvt[k-1] = imax;
      }
      else
      {
        ipvt[k-1] = k;
      }
    }
//
//  2 x 2 pivot block.
//
    else
    {
      km1k = ik + k - 1;
      ikm1 = ik - ( k - 1 );

      if ( swap )
      {
        zswap ( imax, ap+im, 1, ap+ikm1, 1 );
        imj = ikm1 + imax;

        for ( jj = imax; jj <= km1; jj++ )
        {
          j = km1 + imax - jj;
          jkm1 = ikm1 + j;

          t          = ap[jkm1-1];
          ap[jkm1-1] = ap[imj-1];
          ap[imj-1]  = t;

          imj = imj - ( j - 1 );
        }
        t          = ap[km1k-1];
        ap[km1k-1] = ap[imk-1];
        ap[imk-1]  = t;
      }
//
//  Perform the elimination.
//
      km2 = k - 2;

      if ( km2 != 0 )
      {
        ak = ap[kk-1] / ap[km1k-1];
        km1km1 = ikm1 + k - 1;
        akm1 = ap[km1km1-1] / ap[km1k-1];
        denom = complex <double> ( 1.0, 0.0 ) - ak * akm1;
        ij = ik - ( k - 1 ) - ( k - 2 );

        for ( jj = 1; jj <= km2; jj++ )
        {
          j = km1 - jj;
          jk = ik + j;
          bk = ap[jk-1] / ap[km1k-1];
          jkm1 = ikm1 + j;
          bkm1 = ap[jkm1-1] / ap[km1k-1];
          mulk = ( akm1 * bk - bkm1 ) / denom;
          mulkm1 = ( ak * bkm1 - bk ) / denom;
          t = mulk;
          zaxpy ( j, t, ap+ik, 1, ap+ij, 1 );
          t = mulkm1;
          zaxpy ( j, t, ap+ikm1, 1, ap+ij, 1 );
          ap[jk-1] = mulk;
          ap[jkm1-1] = mulkm1;
          ijj = ij + j;
          ij = ij - ( j - 1 );
        }
      }
//
//  Set the pivot array.
//
      if ( swap )
      {
        ipvt[k-1] = -imax;
      }
      else
      {
        ipvt[k-1] = 1 - k;
      }

      ipvt[k-2] = ipvt[k-1];
    }
    ik = ik - ( k - 1 );

    if ( kstep == 2 )
    {
      ik = ik - ( k - 2 );
    }
    k = k - kstep;
  }
  return info;
}
//****************************************************************************80

void zspsl ( complex <double> ap[], int n, int ipvt[], complex <double> b[] )

//****************************************************************************80
//
//  Purpose:
//
//    ZSPSL solves a complex symmetric system factored by ZSPFA.
//
//  Discussion:
//
//    A division by zero may occur if ZSPCO has set RCOND == 0.0
//    or ZSPFA has set INFO != 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> AP[N*(N+1)/2], the output from ZSPFA.
//
//    Input, int N, the order of the matrix.
//
//    Input, int IPVT[N], the pivot vector from ZSPFA.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  complex <double> ak;
  complex <double> akm1;
  complex <double> bk;
  complex <double> bkm1;
  complex <double> denom;
  int ik;
  int ikm1;
  int ikp1;
  int k;
  int kk;
  int km1k;
  int km1km1;
  int kp;
  complex <double> t;
//
//  Loop backward applying the transformations and d inverse to b.
//
  k = n;
  ik = ( n * ( n - 1 ) ) / 2;

  while ( 0 < k )
  {
    kk = ik + k;
    if ( 0 <= ipvt[k-1] )
    {
//
//  1 x 1 pivot block.
//
      if ( k != 1 )
      {
        kp = ipvt[k-1];
        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
        zaxpy ( k-1, b[k-1], ap+ik, 1, b, 1 );
      }
//
//  Apply D inverse.
//
      b[k-1] = b[k-1] / ap[kk-1];
      k = k - 1;
      ik = ik - k;
    }
//
//  2 x 2 pivot block.
//
    else
    {
      ikm1 = ik - ( k - 1 );

      if ( k != 2 )
      {
        kp = abs ( ipvt[k-1] );

        if ( kp != k - 1 )
        {
          t       = b[k-2];
          b[k-2]  = b[kp-1];
          b[kp-2] = t;
        }
        zaxpy ( k-2, b[k-1], ap+ik, 1, b, 1 );
        zaxpy ( k-2, b[k-2], ap+ikm1, 1, b, 1 );
      }
//
//  Apply D inverse.
//
      km1k = ik + k - 1;
      kk = ik + k;
      ak = ap[kk-1] / ap[km1k-1];
      km1km1 = ikm1 + k - 1;
      akm1 = ap[km1km1-1] / ap[km1k-1];
      bk = b[k-1] / ap[km1k-1];
      bkm1 = b[k-2] / ap[km1k-1];
      denom = ak * akm1 - complex <double> ( 1.0, 0.0 );
      b[k-1] = ( akm1 * bk - bkm1 ) / denom;
      b[k-2] = ( ak * bkm1 - bk ) / denom;
      k = k - 2;
      ik = ik - ( k + 1 ) - k;
    }
  }
//
//  Loop forward applying the transformations.
//
  k = 1;
  ik = 0;

  while ( k <= n )
  {
//
//  1 x 1 pivot block.
//
    if ( 0 <= ipvt[k-1] )
    {
      if ( k != 1 )
      {
        b[k-1] = b[k-1] + zdotu ( k-1, ap+ik, 1, b, 1 );
        kp = ipvt[k-1];
        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
      }
      ik = ik + k;
      k = k + 1;
    }
//
//  2 x 2 pivot block.
//
    else
    {
      if ( k != 1 )
      {
        b[k-1] = b[k-1] + zdotu ( k-1, ap+ik, 1, b, 1 );
        ikp1 = ik + k;
        b[k] = b[k] + zdotu ( k-1, ap+ikp1, 1, b, 1 );
        kp = abs ( ipvt[k-1] );

        if ( kp != k )
        {
          t       = b[k-1];
          b[k-1]  = b[kp-1];
          b[kp-1] = t;
        }
      }
      ik = ik + k + k + 1;
      k = k + 2;
    }
  }
  return;
}
//****************************************************************************80

int zsvdc ( complex <double> x[], int ldx, int n, int p,
  complex <double> s[], complex <double> e[], complex <double> u[], int ldu,
  complex <double> v[], int ldv, int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZSVDC applies the singular value decompostion to an N by P matrix.
//
//  Discussion:
//
//    The routine reduces a complex <double> N by P matrix X, by unitary
//    transformations U and V, to diagonal form.
//
//    The diagonal elements, S(I), are the singular values of Z.  The
//    columns of U are the corresponding left singular vectors,
//    and the columns of V the right singular vectors.
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
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> X[LDX*P]; on input, the matrix whose singular
//    value decomposition is to be computed.  X is destroyed on output.
//
//    Input, int LDX, the leading dimension of X.  N <= LDX.
//
//    Input, int N, the number of rows of the matrix.
//
//    Input, int P, the number of columns of the matrix X.
//
//    Output, complex <double> S[MM], where MM = min ( N + 1, P ), the first min ( N, P )
//    entries of S contain the singular values of X arranged in descending
//    order of magnitude.
//
//    Output, complex <double> E[MM], where MM = min ( N + 1, P ),
//    ordinarily contains zeros on output.  However, see the discussion
//    of INFO for exceptions.
//
//    Output, complex <double> U[LDU*K].  If JOBA == 1 then K == n; if JOBA >= 2,
//    then K == min ( N, P ).  U contains the matrix of left singular vectors.
//    U is not referenced if JOBA == 0.  If N <= P or if JOBA > 2,
//    then U may be identified with X in the subroutine call.
//
//    Input, int LDU, the leading dimension of U.  N <= LDU.
//
//    Output, complex <double> V[LDV*P], if requested, the matrix of right singular
//    vectors.  If P <= N, V may be identified with X in the subroutine call.
//
//    Input, int LDV, the leading dimension of V.  P <= LDV.
//
//    Input, int JOB, controls the computation of the singular vectors.
//    It has the decimal expansion AB meaning:
//    A =  0, do not compute the left singular vectors.
//    A =  1, return the N left singular vectors in U.
//    A >= 2, returns the first min ( N, P ) left singular vectors in U.
//    B =  0, do not compute the right singular vectors.
//    B =  1, return the right singular vectors in V.
//
//    Output, int ZSVDC, the value of INFO.  The singular values and their
//    corresponding singular vectors are correct for entries,
//    S(INFO+1), S(INFO+2), ..., S(M).  Here M = min ( N, P ).  Thus if
//    INFO == 0, all the singular values and their vectors are correct.
//    In any event, the matrix
//      B = hermitian(U)*X*V
//    is the bidiagonal matrix with the elements of S on its diagonal
//    and the elements of E on its super-diagonal.  Hermitian(U)
//    is the conjugate-transpose of U.  Thus the singular values of X
//    and B are the same.
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
  int lp1;
  int ls;
  int lu;
  int m;
  int maxit = 30;
  int mm;
  int mm1;
  int mp1;
  int nct;
  int nctp1;
  int ncu;
  int nrt;
  int nrtp1;
  complex <double> r;
  double scale;
  double shift;
  double sl;
  double sm;
  double smm1;
  double sn;
  complex <double> t;
  double t1;
  double test;
  bool wantu;
  bool wantv;
  complex <double> *work;
  double ztest;

  work = new complex <double> [n];
//
//  Determine what is to be computed.
//
  wantu = false;
  wantv = false;
  jobu = ( job % 100 ) / 10;

  if ( 1 < jobu )
  {
    ncu = i4_min ( n, p );
  }
  else
  {
    ncu = n;
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
//  Reduce X to bidiagonal form, storing the diagonal elements
//  in S and the super-diagonal elements in E.
//
  info = 0;
  nct = i4_min ( n - 1, p );
  nrt = i4_max ( 0, i4_min ( p - 2, n ) );
  lu = i4_max ( nct, nrt );

  for ( l = 1; l <= lu; l++ )
  {
    lp1 = l + 1;
//
//  Compute the transformation for the L-th column and
//  place the L-th diagonal in S(L).
//
    if ( l <= nct )
    {
      s[l-1] = complex <double> ( dznrm2 ( n-l+1, x+l-1+(l-1)*ldx, 1 ), 0.0 );

      if ( zabs1 ( s[l-1] ) != 0.0 )
      {
        if ( zabs1 ( x[l-1+(l-1)*ldx] ) != 0.0 )
        {
          s[l-1] = zsign2 ( s[l-1], x[l-1+(l-1)*ldx] );
        }
        t = complex <double> ( 1.0, 0.0 ) / s[l-1];
        zscal ( n-l+1, t, x+l-1+(l-1)*ldx, 1 );
        x[l-1+(l-1)*ldx] = complex <double> ( 1.0, 0.0 ) + x[l-1+(l-1)*ldx];
      }
      s[l-1] = -s[l-1];
    }

    for ( j = lp1; j <= p; j++ )
    {
      if ( l <= nct )
      {
        if ( zabs1 ( s[l-1] ) != 0.0 )
        {
          t = -zdotc ( n-l+1, x+l-1+(l-1)*ldx, 1, x+l-1+(j-1)*ldx, 1 )
            / x[l-1+(l-1)*ldx];
          zaxpy ( n-l+1, t, x+l-1+(l-1)*ldx, 1, x+l-1+(j-1)*ldx, 1 );
        }
      }
//
//  Place the L-th row of X into E for the
//  subsequent calculation of the row transformation.
//
      e[j-1] = conj ( x[l-1+(j-1)*ldx] );
    }
//
//  Place the transformation in U for subsequent back multiplication.
//
    if ( wantu && l <= nct )
    {
      for ( i = l; i <= n; i++ )
      {
        u[i-1+(l-1)*ldu] = x[i-1+(l-1)*ldx];
      }
    }

    if ( l <= nrt )
    {
//
//  Compute the L-th row transformation and place the
//  L-th super-diagonal in E(L).
//
      e[l-1] = complex <double> ( dznrm2 ( p-l, e+lp1-1, 1 ), 0.0 );

      if ( zabs1 ( e[l-1] ) != 0.0 )
      {
        if ( zabs1 ( e[lp1-1] ) != 0.0 )
        {
          e[l-1] = zsign2 ( e[l-1], e[lp1-1] );
        }
        t = complex <double> ( 1.0, 0.0 ) / e[l-1];
        zscal ( p-l, t, e+lp1-1, 1 );
        e[lp1-1] = complex <double> ( 1.0, 0.0 ) + e[lp1-1];
      }

      e[l-1] = -conj ( e[l-1] );
//
//  Apply the transformation.
//
      if ( lp1 <= n && zabs1 ( e[l-1] ) != 0.0 )
      {
        for ( j = lp1; j <= n; j++ )
        {
          work[j-1] = complex <double> ( 0.0, 0.0 );
        }
        for ( j = lp1; j <= p; j++ )
        {
          zaxpy ( n-l, e[j-1], x+lp1-1+(j-1)*ldx, 1, work+lp1-1, 1 );
        }
        for ( j = lp1; j <= p; j++ )
        {
          zaxpy ( n-l, conj ( -e[j-1] / e[lp1-1] ), work+lp1-1,
            1, x+lp1-1+(j-1)*ldx, 1 );
        }
      }
//
//  Place the transformation in V for subsequent back multiplication.
//
      if ( wantv )
      {
        for ( i = lp1; i <= p; i++ )
        {
          v[i-1+(l-1)*ldv] = e[i-1];
        }
      }
    }
  }
//
//  Set up the final bidiagonal matrix of order M.
//
  m = i4_min ( p, n + 1 );
  nctp1 = nct + 1;
  nrtp1 = nrt + 1;

  if ( nct < p )
  {
    s[nctp1-1] = x[nctp1-1+(nctp1-1)*ldx];
  }

  if ( n < m )
  {
    s[m-1] = complex <double> ( 0.0, 0.0 );
  }

  if ( nrtp1 < m )
  {
    e[nrtp1-1] = x[nrtp1-1+(m-1)*ldx];
  }

  e[m-1] = complex <double> ( 0.0, 0.0 );
//
//  If required, generate U.
//
  if ( wantu )
  {
    for ( j = nctp1; j <= ncu; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        u[i-1+(j-1)*ldu] = complex <double> ( 0.0, 0.0 );
      }
      u[j-1+(j-1)*ldu] = complex <double> ( 1.0, 0.0 );
    }
    for ( ll = 1; ll <= nct; ll++ )
    {
      l = nct - ll + 1;

      if ( zabs1 ( s[l-1] ) != 0.0 )
      {
        lp1 = l + 1;

        for ( j = l+1; j <= ncu; j++ )
        {
          t = -zdotc ( n-l+1, u+l-1+(l-1)*ldu, 1, u+l-1+(j-1)*ldu, 1 )
            / u[l-1+(l-1)*ldu];
          zaxpy ( n-l+1, t, u+l-1+(l-1)*ldu, 1, u+l-1+(j-1)*ldu, 1 );
        }

        zscal ( n-l+1, complex <double> ( -1.0, 0.0 ), u+l-1+(l-1)*ldu, 1 );
        u[l-1+(l-1)*ldu] = complex <double> ( 1.0, 0.0 ) + u[l-1+(l-1)*ldu];
        for ( i = 1; i <= l-1; i++ )
        {
          u[i-1+(l-1)*ldu] = complex <double> ( 0.0, 0.0 );
        }
      }
      else
      {
        for ( i = 1; i <= n; i++ )
        {
          u[i-1+(l-1)*ldu] = complex <double> ( 0.0, 0.0 );
        }
        u[l-1+(l-1)*ldu] = complex <double> ( 1.0, 0.0 );
      }
    }
  }
//
//  If it is required, generate V.
//
  if ( wantv )
  {
    for ( ll = 1; ll <= p; ll++ )
    {
      l = p - ll + 1;
      lp1 = l + 1;

      if ( l <= nrt )
      {
        if ( zabs1 ( e[l-1] ) != 0.0 )
        {
          for ( j = lp1; j <= p; j++ )
          {
            t = -zdotc ( p-l, v+lp1-1+(l-1)*ldv, 1, v+lp1-1+(j-1)*ldv, 1 )
              / v[lp1-1+(l-1)*ldv];
            zaxpy ( p-l, t, v+lp1-1+(l-1)*ldv, 1, v+lp1-1+(j-1)*ldv, 1 );
          }
        }
      }
      for ( i = 1; i <= p; i++ )
      {
        v[i-1+(l-1)*ldv] = complex <double> ( 0.0, 0.0 );
      }
      v[l-1+(l-1)*ldv] = complex <double> ( 1.0, 0.0 );
    }
  }
//
//  Transform S and E so that they are real.
//
  for ( i = 1; i <= m; i++ )
  {
    if ( zabs1 ( s[i-1] ) != 0.0 )
    {
      t = complex <double> ( abs ( s[i-1] ), 0.0 );
      r = s[i-1] / t;
      s[i-1] = t;

      if ( i < m )
      {
        e[i-1] = e[i-1] / r;
      }

      if ( wantu )
      {
        zscal ( n, r, u+0+(i-1)*ldu, 1 );
      }
    }

    if ( i == m )
    {
      break;
    }

    if ( zabs1 ( e[i-1] ) != 0.0 )
    {
      t = complex <double> ( abs ( e[i-1] ), 0.0 );
      r = t / e[i-1];
      e[i-1] = t;
      s[i] = s[i] * r;

      if ( wantv )
      {
        zscal ( p, r, v+0+i*ldv, 1 );
      }
    }
  }
//
//  Main iteration loop for the singular values.
//
  mm = m;
  iter = 0;

  for ( ; ; )
  {
//
//  Quit if all the singular values have been found.
//
    if ( m == 0 )
    {
      break;
    }
//
//  If too many iterations have been performed, set flag and return.
//
    if ( maxit <= iter )
    {
      info = m;
      break;
    }
//
//  This section of the program inspects for negligible elements in S and E.
//
//  On completion, the variables KASE and L are set as follows.
//
//  KASE = 1     if S(M) and E(L-1) are negligible and L < M
//  KASE = 2     if S(L) is negligible and L < M
//  KASE = 3     if E(L-1) is negligible, L < M, and
//               S(L), ..., S(M) are not negligible (QR step).
//  KASE = 4     if E(M-1) is negligible (convergence).
//
    for ( ll = 1; ll <= m; ll++ )
    {
      l = m - ll;

      if ( l == 0 )
      {
        break;
      }

      test = abs ( s[l-1] ) + abs ( s[l] );
      ztest = test + abs ( e[l-1] );

      if ( ztest == test )
      {
        e[l-1] = complex <double> ( 0.0, 0.0 );
        break;
      }
    }

    if ( l == m - 1 )
    {
      kase = 4;
    }
    else
    {
      lp1 = l + 1;
      mp1 = m + 1;

      for ( lls = lp1; lls <= mp1; lls++ )
      {
        ls = m - lls + lp1;

        if ( ls == l )
        {
          break;
        }

        test = 0.0;

        if ( ls != m )
        {
          test = test + abs ( e[ls-1] );
        }

        if ( ls != l + 1 )
        {
          test = test + abs ( e[ls-2] );
        }

        ztest = test + abs ( s[ls-1] );

        if ( ztest == test)
        {
          s[ls-1] = complex <double> ( 0.0, 0.0 );
          break;
        }
      }
      if ( ls == l )
      {
        kase = 3;
      }
      else if ( ls == m )
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
//  Deflate negligible S(M).
//
    if ( kase == 1 )
    {
      mm1 = m - 1;
      f = real ( e[m-2] );
      e[m-2] = complex <double> ( 0.0, 0.0 );

      for ( kk = 1; kk <= mm1; kk++ )
      {
        k = mm1 - kk + l;
        t1 = real ( s[k-1] );
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = complex <double> ( t1, 0.0 );

        if ( k != l )
        {
          f = -sn * real ( e[k-2] );
          e[k-2] = cs * e[k-2];
        }

        if ( wantv )
        {
          zdrot ( p, v+0+(k-1)*ldv, 1, v+0+(m-1)*ldv, 1, cs, sn );
        }
      }
    }
//
//  Split at negligible S(L).
//
    else if ( kase == 2 )
    {
      f = real ( e[l-2] );
      e[l-2] = complex <double> ( 0.0, 0.0 );

      for ( k = l; k <= m; k++ )
      {
        t1 = real ( s[k-1] );
        drotg ( &t1, &f, &cs, &sn );
        s[k-1] = complex <double> ( t1, 0.0 );
        f = -sn * real ( e[k-1] );
        e[k-1] = cs * e[k-1];

        if ( wantu )
        {
          zdrot ( n, u+0+(k-1)*ldu, 1, u+0+(l-2)*ldu, 1, cs, sn );
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
      scale = r8_max ( abs ( s[m-1] ),
              r8_max ( abs ( s[m-2] ),
              r8_max ( abs ( e[m-2] ),
              r8_max ( abs ( s[l-1] ), abs ( e[l-1] ) ) ) ) );

      sm = real ( s[m-1] ) / scale;
      smm1 = real ( s[m-2] ) / scale;
      emm1 = real ( e[m-2] ) / scale;
      sl = real ( s[l-1] ) / scale;
      el = real ( e[l-1] ) / scale;
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

      f = ( sl + sm ) * ( sl - sm ) + shift;
      g = sl * el;
//
//  Chase zeros.
//
      mm1 = m - 1;

      for ( k = l; k <= mm1; k++ )
      {
        drotg ( &f, &g, &cs, &sn );

        if ( k != l )
        {
         e[k-2] = complex <double> ( f, 0.0 );
        }

        f = cs * real ( s[k-1] ) + sn * real ( e[k-1] );
        e[k-1] = cs * e[k-1] - sn * s[k-1];
        g = sn * real ( s[k] );
        s[k] = cs * s[k];

        if ( wantv )
        {
          zdrot ( p, v+0+(k-1)*ldv, 1, v+0+k*ldv, 1, cs, sn );
        }

        drotg ( &f, &g, &cs, &sn );
        s[k-1] = complex <double> ( f, 0.0 );
        f = cs * real ( e[k-1] ) + sn * real ( s[k] );
        s[k] = -sn * e[k-1] + cs * s[k];
        g = sn * real ( e[k] );
        e[k] = cs * e[k];

        if ( wantu && k < n )
        {
          zdrot ( n, u+0+(k-1)*ldu, 1, u+0+k*ldu, 1, cs, sn );
        }
      }
      e[m-2] = complex <double> ( f, 0.0 );
      iter = iter + 1;
    }
//
//  Convergence.
//
    else if ( kase == 4 )
    {
//
//  Make the singular value positive.
//
      if ( real ( s[l-1] ) < 0.0 )
      {
        s[l-1] = -s[l-1];
        if ( wantv )
        {
          zscal ( p, complex <double> ( -1.0, 0.0 ), v+0+(l-1)*ldv, 1 );
        }
      }
//
//  Order the singular values.
//
      while ( l != mm )
      {
        if ( real ( s[l] ) <= real ( s[l-1] ) )
        {
          break;
        }

        t      = s[l-1];
        s[l-1] = s[l];
        s[l]   = t;

        if ( wantv && l < p )
        {
          zswap ( p, v+0+(l-1)*ldv, 1, v+0+l*ldv, 1 );
        }

        if ( wantu && l < n )
        {
          zswap ( n, u+0+(l-1)*ldu, 1, u+0+l*ldu, 1 );
        }
        l = l + 1;
      }
      iter = 0;
      m = m - 1;
    }
  }
  delete [] work;

  return info;
}
//****************************************************************************80

double ztrco ( complex <double> t[], int ldt, int n, int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZTRCO estimates the condition of a complex triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> T[LDT*N], the triangular matrix.  The zero
//    elements of the matrix are not referenced, and the corresponding
//    elements of the array can be used to store other information.
//
//    Input, int LDT, the leading dimension of T.
//
//    Input, int N, the order of the matrix.
//
//    Input, int JOB, indicates if matrix is upper or lower triangular.
//    0, lower triangular.
//    nonzero, upper triangular.
//
//    Output, double ZTRCO, an estimate of RCOND, the reciprocal condition of T.
//    For the system T*X = B, relative perturbations in T and B of size
//    EPSILON may cause relative perturbations in X of size (EPSILON/RCOND).
//    If RCOND is so small that the logical expression
//      1.0 + RCOND == 1.0
//    is true, then T may be singular to working precision.  In particular,
//    RCOND is zero if exact singularity is detected or the estimate
//    underflows.
//
//  Local Parameters:
//
//    Workspace, complex <double> Z[N], a work vector whose contents are usually
//    unimportant.  If T is close to a singular matrix, then Z is
//    an approximate null vector in the sense that
//      norm(A*Z) = RCOND * norm(A) * norm(Z).
//
{
  complex <double> ek;
  int i;
  int i1;
  int j;
  int j1;
  int j2;
  int k;
  int kk;
  int l;
  bool lower;
  double rcond;
  double s;
  double sm;
  double tnorm;
  complex <double> w;
  complex <double> wk;
  complex <double> wkm;
  double ynorm;
  complex <double> *z;

  lower = ( job == 0 );
//
//  Compute 1-norm of T
//
  tnorm = 0.0;

  for ( j = 1; j <= n; j++ )
  {
    if ( lower )
    {
      l = n + 1 - j;
      i1 = j;
    }
    else
    {
      l = j;
      i1 = 1;
    }
    tnorm = r8_max ( tnorm, dzasum ( l, t+i1-1+(j-1)*ldt, 1 ) );
  }
//
//  RCOND = 1/(norm(T)*(estimate of norm(inverse(T)))).
//
//  Estimate = norm(Z)/norm(Y) where T*Z = Y and hermitian(T)*Y = E.
//
//  Hermitian(T) is the conjugate transpose of T.
//
//  The components of E are chosen to cause maximum local
//  growth in the elements of Y.
//
//  The vectors are frequently rescaled to avoid overflow.
//
//  Solve hermitian(T)*Y = E.
//
  ek = complex <double> ( 1.0, 0.0 );
  z = new complex <double> [n];
  for ( i = 0; i < n; i++ )
  {
    z[i] = complex <double> ( 0.0, 0.0 );
  }

  for ( kk = 1; kk <= n; kk++ )
  {
    if ( lower )
    {
      k = n + 1 - kk;
    }
    else
    {
      k = kk;
    }

    if ( zabs1 ( z[k-1] ) != 0.0 )
    {
      ek = zsign1 ( ek, -z[k-1] );
    }

    if ( zabs1 ( t[k-1+(k-1)*ldt] ) < zabs1 ( ek - z[k-1] ) )
    {
      s = zabs1 ( t[k-1+(k-1)*ldt] ) / zabs1 ( ek - z[k-1] );
      zdscal ( n, s, z, 1 );
      ek = complex <double> ( s, 0.0 ) * ek;
    }

    wk = ek - z[k-1];
    wkm = - ek - z[k-1];
    s = zabs1 ( wk );
    sm = zabs1 ( wkm );

    if ( zabs1 ( t[k-1+(k-1)*ldt] ) != 0.0 )
    {
      wk = wk / conj ( t[k-1+(k-1)*ldt] );
      wkm = wkm / conj ( t[k-1+(k-1)*ldt] );
    }
    else
    {
      wk = complex <double> ( 1.0, 0.0 );
      wkm = complex <double> ( 1.0, 0.0 );
    }

    if ( kk != n )
    {
      if (lower)
      {
        j1 = 1;
        j2 = k - 1;
      }
      else
      {
        j1 = k + 1;
        j2 = n;
      }

      for ( j = j1; j <= j2; j++ )
      {
        sm = sm + zabs1 ( z[j-1] + wkm * conj ( t[k-1+(j-1)*ldt] ) );
        z[j-1] = z[j-1] + wk * conj ( t[k-1+(j-1)*ldt] );
        s = s + zabs1 ( z[j-1] );
      }

      if ( s < sm )
      {
        w = wkm - wk;
        wk = wkm;
        for ( j = j1; j <= j2; j++ )
        {
          z[j-1] = z[j-1] + w * conj ( t[k-1+(j-1)*ldt] );
        }
      }
    }
    z[k-1] = wk;
  }

  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = 1.0;
//
//  Solve T*Z = Y.
//
  for ( kk = 1; kk <= n; kk++ )
  {
    if ( lower )
    {
      k = kk;
    }
    else
    {
      k = n + 1 - kk;
    }

    if ( zabs1 ( t[k-1+(k-1)*ldt] ) < zabs1 ( z[k-1] ) )
    {
      s = zabs1 ( t[k-1+(k-1)*ldt] ) / zabs1 ( z[k-1] );
      zdscal ( n, s, z, 1 );
      ynorm = s * ynorm;
    }

    if ( zabs1 ( t[k-1+(k-1)*ldt] ) != 0.0 )
    {
      z[k-1] = z[k-1] / t[k-1+(k-1)*ldt];
    }
    else
    {
      z[k-1] = complex <double> ( 1.0, 0.0 );
    }

    if ( lower )
    {
      i1 = k + 1;
    }
    else
    {
      i1 = 1;
    }

    if ( kk < n )
    {
      w = -z[k-1];
      zaxpy ( n-kk, w, t+i1-1+(k-1)*ldt, 1, z+i1-1, 1 );
    }
  }
//
//  Make ZNORM = 1.
//
  s = 1.0 / dzasum ( n, z, 1 );
  zdscal ( n, s, z, 1 );
  ynorm = s * ynorm;

  if ( tnorm != 0.0 )
  {
    rcond = ynorm / tnorm;
  }
  else
  {
    rcond = 0.0;
  }

  delete [] z;

  return rcond;
}
//****************************************************************************80

int ztrdi ( complex <double> t[], int ldt, int n, complex <double> det[2],
  int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZTRDI computes the determinant and inverse of a complex triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input/output, complex <double> T[LDT*N], the triangular matrix.  The zero
//    elements of the matrix are not referenced, and the corresponding
//    elements of the array can be used to store other information.
//    On output, if an inverse was requested, then T has been overwritten
//    by its inverse.
//
//    Input, int LDT, the leading dimension of T.
//
//    Input, int N, the order of the matrix.
//
//    Input, int JOB.
//    010, no determinant,    inverse, matrix is lower triangular.
//    011, no determinant,    inverse, matrix is upper triangular.
//    100,    determinant, no inverse.
//    110,    determinant,    inverse, matrix is lower triangular.
//    111,    determinant,    inverse, matrix is upper triangular.
//
//    Output, complex <double> DET[2], the determinant of the original matrix,
//    if requested.  Otherwise not referenced.
//    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= zabs1 ( DET(1) ) < 10.0
//    or DET(1) == 0.0.  Also, DET(2) is strictly real.
//
//    Output, int ZTRDI.
//    0, an inverse was requested and the matrix is nonsingular.
//    K, an inverse was requested, but the K-th diagonal element
//    of T is zero.
//
{
  int i;
  int info;
  int j;
  int k;
  complex <double> temp;

  if ( ( job / 100 ) != 0 )
  {
    det[0] = complex <double> ( 1.0, 0.0 );
    det[1] = complex <double> ( 0.0, 0.0 );

    for ( i = 0; i < n; i++ )
    {
      det[0] = det[0] * t[i+i*ldt];

      if ( zabs1 ( det[0] ) == 0.0 )
      {
        break;
      }

      while ( zabs1 ( det[0] ) < 1.0 )
      {
        det[0] = det[0] * complex <double> ( 10.0, 0.0 );
        det[1] = det[1] - complex <double> (  1.0, 0.0 );
      }

      while ( 10.0 <= zabs1 ( det[0] ) )
      {
        det[0] = det[0] / complex <double> ( 10.0, 0.0 );
        det[1] = det[1] + complex <double> (  1.0, 0.0 );
      }
    }
  }
//
//  Compute inverse of upper triangular matrix.
//
  if ( ( ( job / 10 ) % 10 ) != 0 )
  {
    if ( ( job % 10 ) != 0 )
    {
      info = 0;

      for ( k = 0; k < n; k++ )
      {
        if ( zabs1 ( t[k+k*ldt] ) == 0.0 )
        {
          info = k + 1;
          break;
        }

        t[k+k*ldt] = complex <double> ( 1.0, 0.0 ) / t[k+k*ldt];
        temp = -t[k+k*ldt];
        zscal ( k, temp, t+0+k*ldt, 1 );

        for ( j = k+1; j < n; j++ )
        {
          temp = t[k+j*ldt];
          t[k+j*ldt] = complex <double> ( 0.0, 0.0 );
          zaxpy ( k+1, temp, t+0+k*ldt, 1, t+0+j*ldt, 1 );
        }
      }
    }
//
//  Compute inverse of lower triangular matrix.
//
    else
    {
      info = 0;

      for ( k = n-1; 0 <= k; k-- )
      {
        if ( zabs1 ( t[k+k*ldt] ) == 0.0 )
        {
          info = k + 1;
          break;
        }

        t[k+k*ldt] = complex <double> ( 1.0, 0.0 ) / t[k+k*ldt];

        if ( k != n - 1 )
        {
          temp = -t[k+k*ldt];
          zscal ( n-k-1, temp, t+k+1+k*ldt, 1 );
        }

        for ( j = 0; j < k; j++ )
        {
          temp = t[k+j*ldt];
          t[k+j*ldt] = complex <double> ( 0.0, 0.0 );
          zaxpy ( n-k, temp, t+k+k*ldt, 1, t+k+j*ldt, 1 );
        }
      }
    }
  }

  return info;
}
//****************************************************************************80

int ztrsl ( complex <double> t[], int ldt, int n, complex <double> b[],
  int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZTRSL solves triangular systems T*X=B or Hermitian(T)*X=B.
//
//  Discussion:
//
//    Hermitian ( T ) denotes the conjugate transpose of the matrix T.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//
//  Parameters:
//
//    Input, complex <double> T[LDT*N], the matrix of the system.  The zero
//    elements of the matrix are not referenced, and the corresponding
//    elements of the array can be used to store other information.
//
//    Input, int LDT, the leading dimension of T.
//
//    Input, int N, the order of the matrix.
//
//    Input/output, complex <double> B[N].  On input, the right hand side.
//    On output, the solution.
//
//    Input, int JOB, specifies what kind of system is to be solved.
//    00, solve T*X=B, T lower triangular,
//    01, solve T*X=B, T upper triangular,
//    10, solve hermitian(T)*X=B, T lower triangular,
//    11, solve hermitian(T)*X=B, T upper triangular.
//
//    Output, int ZTRSL.
//    0, the system is nonsingular.
//    K, the index of the first zero diagonal element of T.
//
{
  int kase;
  int i;
  int info;
  int j;
  int jj;
  complex <double> temp;
//
//  Check for zero diagonal elements.
//
  for ( i = 0; i < n; i++ )
  {
    if ( zabs1 ( t[i+i*ldt] ) == 0.0 )
    {
      info = i + 1;
      return info;
    }
  }

  info = 0;
//
//  Determine the task and go to it.
//
  kase = 1;

  if ( ( job % 10 ) != 0 )
  {
    kase = 2;
  }

  if ( ( job % 100 ) / 10 != 0 )
  {
    kase = kase + 2;
  }
//
//  Solve T * X = B for T lower triangular.
//
  if ( kase == 1 )
  {
    b[0] = b[0] / t[0+0*ldt];

    for ( j = 2; j <= n; j++ )
    {
      temp = -b[j-2];
      zaxpy ( n-j+1, temp, t+j-1+(j-2)*ldt, 1, b+j-1, 1 );
      b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
    }
  }
//
//  Solve T * X = B for T upper triangular.
//
  else if ( kase == 2 )
  {
    b[n-1] = b[n-1] / t[n-1+(n-1)*ldt];

    for ( jj = 2; jj <= n; jj++ )
    {
      j = n - jj + 1;
      temp = -b[j];
      zaxpy ( j, temp, t+0+j*ldt, 1, b, 1 );
      b[j-1] = b[j-1] / t[j-1+(j-1)*ldt];
    }
  }
//
//  Solve hermitian(T) * X = B for T lower triangular.
//
  else if ( kase == 3 )
  {
    b[n-1] = b[n-1] / conj ( t[n-1+(n-1)*ldt] );

    for ( jj = 2; jj <= n; jj++ )
    {
      j = n - jj + 1;
      b[j-1] = b[j-1] - zdotc ( jj-1, t+j+(j-1)*ldt, 1, b+j, 1 );
      b[j-1] = b[j-1] / conj ( t[j-1+(j-1)*ldt] );
    }
  }
//
//  Solve hermitian(T) * X = B for T upper triangular.
//
  else if ( kase == 4 )
  {
    b[0] = b[0] / conj ( t[0+0*ldt] );

    for ( j = 2; j <= n; j++ )
    {
      b[j-1] = b[j-1] - zdotc ( j-1, t+0+(j-1)*ldt, 1, b, 1 );
      b[j-1] = b[j-1] / conj ( t[j-1+(j-1)*ldt] );
    }
  }

  return info;
}
