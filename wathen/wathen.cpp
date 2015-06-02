# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "wathen.hpp"

//****************************************************************************80

void bandwidth ( int m, int n, double a[], int &b, int &l, int &d, int &u )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH returns the bandwidth of a matrix.
//
//  Discussion:
//
//    If the nonzeros of a matrix only occur in entries that are "close"
//    to the main diagonal, we say the matrix is banded.
//
//    Roughly speaking, the bandwidth B of a matrix is the number of 
//    diagonals containing nonzeros.  More precisely, it is the minimum number
//    of contiguous diagonals that contain all the nonzeros.  It is presumed
//    that the main diagonal is nonzero.
//
//    We can also measure U and L, the upper and lower "half-bandwidths" which
//    count the number of contiguous diagonals above or below the main
//    diagonal.
//
//    We may write
//      B = L + D + U
//    where D is presumably 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], the matrix.
//
//    Output, int &B, the total bandwidth.
//
//    Output, int &L, &D, &U, the lower, diagonal, and upper 
//    bandwidths.
//
{
  int i;
  int j;

  l = 0;
  d = 0;
  u = 0;

  for ( i = 0; i < n; i++ )
  {
    j = 0;
    while ( l < i - j )
    {
      if ( a[i+j*m] != 0.0 )
      {
        l = i - j;
        break;
      }
      j = j + 1;
    }

    if ( a[i+i*m] != 0.0 )
    {
      d = 1;
    }

    j = n - 1;
    while ( u < j - i )
    {
      if ( a[i+j*m] != 0.0 )
      {
        u = j - i;
        break;
      }
      j = j - 1;
    }
  }

  b = l + d + u;

  return;
}
//****************************************************************************80

void cg_gb ( int n, int ml, int mu, double a[], double b[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    CG_GB uses the conjugate gradient method for a general banded (GB) matrix.
//
//  Discussion:
//
//    The linear system has the form A*x=b, where A is a positive-definite
//    symmetric matrix.
//
//    The method is designed to reach the solution to the linear system
//      A * x = b
//    after N computational steps.  However, roundoff may introduce
//    unacceptably large errors for some problems.  In such a case,
//    calling the routine a second time, using the current solution estimate
//    as the new starting guess, should result in improved results.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Frank Beckman,
//    The Solution of Linear Equations by the Conjugate Gradient Method,
//    in Mathematical Methods for Digital Computers,
//    edited by John Ralston, Herbert Wilf,
//    Wiley, 1967,
//    ISBN: 0471706892,
//    LC: QA76.5.R3.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//
//    Input, double A[(2*ML+MU+1)*N], the band matrix.
//
//    Input, double B[N], the right hand side vector.
//
//    Input/output, double X[N].
//    On input, an estimate for the solution, which may be 0.
//    On output, the approximate solution vector.  
//
{
  double alpha;
  double *ap;
  double beta;
  int i;
  int it;
  double *p;
  double pap;
  double pr;
  double *r;
  double rap;
//
//  Initialize
//    AP = A * x,
//    R  = b - A * x,
//    P  = b - A * x.
//
  ap = mv_gb ( n, n, ml, mu, a, x );

  r = new double[n];
  for ( i = 0; i < n; i++ )
  {
    r[i] = b[i] - ap[i];
  }

  p = new double[n];
  for ( i = 0; i < n; i++ )
  {
    p[i] = b[i] - ap[i];
  }
//
//  Do the N steps of the conjugate gradient method.
//
  for ( it = 1; it <= n; it++ )
  {
//
//  Compute the matrix*vector product AP = A*P.
//
    delete [] ap;
    ap = mv_gb ( n, n, ml, mu, a, p );
//
//  Compute the dot products
//    PAP = P*AP,
//    PR  = P*R
//  Set
//    ALPHA = PR / PAP.
//
    pap = r8vec_dot_product ( n, p, ap );
    pr = r8vec_dot_product ( n, p, r );

    if ( pap == 0.0 )
    {
      break;
    }
    alpha = pr / pap;
//
//  Set
//    X = X + ALPHA * P
//    R = R - ALPHA * AP.
//
    for ( i = 0; i < n; i++ )
    {
      x[i] = x[i] + alpha * p[i];
    }
    for ( i = 0; i < n; i++ )
    {
      r[i] = r[i] - alpha * ap[i];
    }
//
//  Compute the vector dot product
//    RAP = R*AP
//  Set
//    BETA = - RAP / PAP.
//

    rap = r8vec_dot_product ( n, r, ap );

    beta = - rap / pap;
//
//  Update the perturbation vector
//    P = R + BETA * P.
//
    for ( i = 0; i < n; i++ )
    {
      p[i] = r[i] + beta * p[i];
    }
  }

  delete [] ap;
  delete [] p;
  delete [] r;

  return;
}
//****************************************************************************80

void cg_ge ( int n, double a[], double b[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    CG_GE uses the conjugate gradient method for a general storage (GE) matrix.
//
//  Discussion:
//
//    The linear system has the form A*x=b, where A is a positive-definite
//    symmetric matrix, stored as a full storage matrix.
//
//    The method is designed to reach the solution to the linear system
//      A * x = b
//    after N computational steps.  However, roundoff may introduce
//    unacceptably large errors for some problems.  In such a case,
//    calling the routine a second time, using the current solution estimate
//    as the new starting guess, should result in improved results.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Frank Beckman,
//    The Solution of Linear Equations by the Conjugate Gradient Method,
//    in Mathematical Methods for Digital Computers,
//    edited by John Ralston, Herbert Wilf,
//    Wiley, 1967,
//    ISBN: 0471706892,
//    LC: QA76.5.R3.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Input, double B[N], the right hand side vector.
//
//    Input/output, double X[N].
//    On input, an estimate for the solution, which may be 0.
//    On output,  the approximate solution vector.  
//
{
  double alpha;
  double *ap;
  double beta;
  int i;
  int it;
  double *p;
  double pap;
  double pr;
  double *r;
  double rap;
//
//  Initialize
//    AP = A * x,
//    R  = b - A * x,
//    P  = b - A * x.
//
  ap = mv_ge ( n, n, a, x );
 
  r = new double[n];
  for ( i = 0; i < n; i++ )
  {
    r[i] = b[i] - ap[i];
  }

  p = new double[n];
  for ( i = 0; i < n; i++ )
  {
    p[i] = b[i] - ap[i];
  }
//
//  Do the N steps of the conjugate gradient method.
//
  for ( it = 1; it <= n; it++ )
  {
//
//  Compute the matrix*vector product AP = A*P.
//

    delete [] ap;
    ap = mv_ge ( n, n, a, p );
//
//  Compute the dot products
//    PAP = P*AP,
//    PR  = P*R
//  Set
//    ALPHA = PR / PAP.
//
    pap = r8vec_dot_product ( n, p, ap );
    pr = r8vec_dot_product ( n, p, r );

    if ( pap == 0.0 )
    {
      break;
    }

    alpha = pr / pap;
//
//  Set
//    X = X + ALPHA * P
//    R = R - ALPHA * AP.
//
    for ( i = 0; i < n; i++ )
    {
      x[i] = x[i] + alpha * p[i];
    }
    for ( i = 0; i < n; i++ )
    {
      r[i] = r[i] - alpha * ap[i];
    }
//
//  Compute the vector dot product
//    RAP = R*AP
//  Set
//    BETA = - RAP / PAP.
//
    rap = r8vec_dot_product ( n, r, ap );

    beta = - rap / pap;
//
//  Update the perturbation vector
//    P = R + BETA * P.
//
    for ( i = 0; i < n; i++ )
    {
      p[i] = r[i] + beta * p[i];
    }
  }

  delete [] ap;
  delete [] p;
  delete [] r;

  return;
}
//****************************************************************************80

void cg_st ( int n, int nz_num, int row[], int col[], double a[], double b[], 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    CG_ST uses the conjugate gradient method for a sparse triplet (ST) matrix.
//
//  Discussion:
//
//    The linear system has the form A*x=b, where A is a positive-definite
//    symmetric matrix, stored as a full storage matrix.
//
//    The method is designed to reach the solution to the linear system
//      A * x = b
//    after N computational steps.  However, roundoff may introduce
//    unacceptably large errors for some problems.  In such a case,
//    calling the routine a second time, using the current solution estimate
//    as the new starting guess, should result in improved results.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Frank Beckman,
//    The Solution of Linear Equations by the Conjugate Gradient Method,
//    in Mathematical Methods for Digital Computers,
//    edited by John Ralston, Herbert Wilf,
//    Wiley, 1967,
//    ISBN: 0471706892,
//    LC: QA76.5.R3.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int NZ_NUM, the number of nonzeros.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column 
//    indices of the nonzero entries.
//
//    Input, double A[NZ_NUM], the nonzero entries.
//
//    Input, double B[N], the right hand side vector.
//
//    Input/output, double X[N].
//    On input, an estimate for the solution, which may be 0.
//    On output, the approximate solution vector.  
//
{
  double alpha;
  double *ap;
  double beta;
  int i;
  int it;
  double *p;
  double pap;
  double pr;
  double *r;
  double rap;
//
//  Initialize
//    AP = A * x,
//    R  = b - A * x,
//    P  = b - A * x.
//
  ap = mv_st ( n, n, nz_num, row, col, a, x );

  r = new double[n];
  for ( i = 0; i < n; i++ )
  {
    r[i] = b[i] - ap[i];
  }

  p = new double[n];
  for ( i = 0; i < n; i++ )
  {
    p[i] = b[i] - ap[i];
  }
//
//  Do the N steps of the conjugate gradient method.
//
  for ( it = 1; it <= n; it++ )
  {
//
//  Compute the matrix*vector product AP = A*P.
//
    delete [] ap;
    ap = mv_st ( n, n, nz_num, row, col, a, p );
//
//  Compute the dot products
//    PAP = P*AP,
//    PR  = P*R
//  Set
//    ALPHA = PR / PAP.
//
    pap = r8vec_dot_product ( n, p, ap );
    pr = r8vec_dot_product ( n, p, r );

    if ( pap == 0.0 )
    {
      break;
    }

    alpha = pr / pap;
//
//  Set
//    X = X + ALPHA * P
//    R = R - ALPHA * AP.
//
    for ( i = 0; i < n; i++ )
    {
      x[i] = x[i] + alpha * p[i];
    }
    for ( i = 0; i < n; i++ )
    {
      r[i] = r[i] - alpha * ap[i];
    }
//
//  Compute the vector dot product
//    RAP = R*AP
//  Set
//    BETA = - RAP / PAP.
//
    rap = r8vec_dot_product ( n, r, ap );

    beta = - rap / pap;
//
//  Update the perturbation vector
//    P = R + BETA * P.
//
    for ( i = 0; i < n; i++ )
    {
      p[i] = r[i] + beta * p[i];
    }
  }

  delete [] ap;
  delete [] p;
  delete [] r;

  return;
}
//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
//
//    CPU_TIME returns the current reading on the CPU clock.
//
//  Discussion:
//
//    The CPU time measurements available through this routine are often
//    not very accurate.  In some cases, the accuracy is no better than
//    a hundredth of a second.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
//
{
  double value;

  value = ( double ) clock ( ) 
        / ( double ) CLOCKS_PER_SEC;

  return value;
}
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
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
//    David Kincaid, Fred Krogh.
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
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
//    David Kincaid, Fred Krogh.
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

int dgbfa ( double abd[], int lda, int n, int ml, int mu, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGBFA factors a real band matrix by elimination.
//
//  Discussion:
//
//    DGBFA is usually called by DGBCO, but it can be called
//    directly with a saving in time if RCOND is not needed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
//    Pete Stewart.
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
//    Input/output, double ABD[LDA*N].  On input, the matrix in band
//    storage.  The columns of the matrix are stored in the columns of ABD
//    and the diagonals of the matrix are stored in rows ML+1 through
//    2*ML+MU+1 of ABD.  On output, an upper triangular matrix in band storage
//    and the multipliers which were used to obtain it.  The factorization
//    can be written A = L*U where L is a product of permutation and unit lower
//    triangular matrices and U is upper triangular.
//
//    Input, int LDA, the leading dimension of the array ABD.
//    2*ML + MU + 1 <= LDA is required.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, MU, the number of diagonals below and above the
//    main diagonal.  0 <= ML < N, 0 <= MU < N.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, integer DGBFA, error flag.
//    0, normal value.
//    K, if U(K,K) == 0.0D+00.  This is not an error condition for this
//      subroutine, but it does indicate that DGBSL will divide by zero if
//      called.  Use RCOND in DGBCO for a reliable indication of singularity.
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
  double t;

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
      abd[i-1+(jz-1)*lda] = 0.0;
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
//  Zero out the next fill-in column.
//
    jz = jz + 1;
    if ( jz <= n )
    {
      for ( i = 1; i <= ml; i++ )
      {
        abd[i-1+(jz-1)*lda] = 0.0;
      }
    }
//
//  Find L = pivot index.
//
    lm = i4_min ( ml, n-k );
    l = idamax ( lm+1, abd+m-1+(k-1)*lda, 1 ) + m - 1;
    ipvt[k-1] = l + k - m;
//
//  Zero pivot implies this column already triangularized.
//
    if ( abd[l-1+(k-1)*lda] == 0.0 )
    {
      info = k;
    }
//
//  Interchange if necessary.
//
    else
    {
      if ( l != m )
      {
        t = abd[l-1+(k-1)*lda];
        abd[l-1+(k-1)*lda] = abd[m-1+(k-1)*lda];
        abd[m-1+(k-1)*lda] = t;
      }
//
//  Compute multipliers.
//
      t = -1.0 / abd[m-1+(k-1)*lda];
      dscal ( lm, t, abd+m+(k-1)*lda, 1 );
//
//  Row elimination with column indexing.
//
      ju = i4_min ( i4_max ( ju, mu+ipvt[k-1] ), n );
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
        daxpy ( lm, t, abd+m+(k-1)*lda, 1, abd+mm+(j-1)*lda, 1 );
      }

    }

  }

  ipvt[n-1] = n;

  if ( abd[m-1+(n-1)*lda] == 0.0 )
  {
    info = n;
  }

  return info;
}
//****************************************************************************80

void dgbsl ( double abd[], int lda, int n, int ml, int mu, int ipvt[], 
  double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DGBSL solves a real banded system factored by DGBCO or DGBFA.
//
//  Discussion:
//
//    DGBSL can solve either A * X = B  or  A' * X = B.
//
//    A division by zero will occur if the input factor contains a
//    zero on the diagonal.  Technically this indicates singularity
//    but it is often caused by improper arguments or improper
//    setting of LDA.  It will not occur if the subroutines are
//    called correctly and if DGBCO has set 0.0 < RCOND
//    or DGBFA has set INFO == 0.
//
//    To compute inverse(A) * C  where C is a matrix with P columns:
//
//      call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
//
//      if ( rcond is too small ) then
//        exit
//      end if
//
//      do j = 1, p
//        call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
//      end do
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
//    Pete Stewart.
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
//    Input, double ABD[LDA*N], the output from DGBCO or DGBFA.
//
//    Input, integer LDA, the leading dimension of the array ABD.
//
//    Input, int N, the order of the matrix.
//
//    Input, int ML, MU, the number of diagonals below and above the
//    main diagonal.  0 <= ML < N, 0 <= MU < N.
//
//    Input, int IPVT[N], the pivot vector from DGBCO or DGBFA.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
//    Input, int JOB, job choice.
//    0, solve A*X=B.
//    nonzero, solve A'*X=B.
//
{
  int k;
  int l;
  int la;
  int lb;
  int lm;
  int m;
  double t;

  m = mu + ml + 1;
//
//  JOB = 0, Solve A * x = b.
//
//  First solve L * y = b.
//
  if ( job == 0 )
  {
    if ( 0 < ml )
    {
      for ( k = 1; k <= n-1; k++ )
      {
        lm = i4_min ( ml, n-k );
        l = ipvt[k-1];
        t = b[l-1];
        if ( l != k )
        {
          b[l-1] = b[k-1];
          b[k-1] = t;
        }
        daxpy ( lm, t, abd+m+(k-1)*lda, 1, b+k, 1 );
      }
    }
//
//  Now solve U * x = y.
//
    for ( k = n; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] / abd[m-1+(k-1)*lda];
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      t = -b[k-1];
      daxpy ( lm, t, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
    }
  }
//
//  JOB nonzero, solve A' * x = b.
//
//  First solve U' * y = b.
//
  else
  {
    for ( k = 1; k <= n; k++ )
    {
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      t = ddot ( lm, abd+la-1+(k-1)*lda, 1, b+lb-1, 1 );
      b[k-1] = ( b[k-1] - t ) / abd[m-1+(k-1)*lda];
    }
//
//  Now solve L' * x = y.
//
    if ( 0 < ml )
    {
      for ( k = n-1; 1 <= k; k-- )
      {
        lm = i4_min ( ml, n-k );
        b[k-1] = b[k-1] + ddot ( lm, abd+m+(k-1)*lda, 1, b+k, 1 );
        l = ipvt[k-1];
        if ( l != k )
        {
          t = b[l-1];
          b[l-1] = b[k-1];
          b[k-1] = t;
        }
      }
    }
  }

  return;
}
//****************************************************************************80

int dgefa ( double a[], int lda, int n, int ipvt[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGEFA factors a real general matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
//    Pete Stewart.
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
//    Input/output, double A[LDA*N].
//    On intput, the matrix to be factored.
//    On output, an upper triangular matrix and the multipliers used to obtain
//    it.  The factorization can be written A=L*U, where L is a product of
//    permutation and unit lower triangular matrices, and U is upper triangular.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix A.
//
//    Output, int IPVT[N], the pivot indices.
//
//    Output, int DGEFA, singularity indicator.
//    0, normal value.
//    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
//    but it does indicate that DGESL or DGEDI will divide by zero if called.
//    Use RCOND in DGECO for a reliable indication of singularity.
//
{
  int info;
  int j;
  int k;
  int l;
  double t;
//
//  Gaussian elimination with partial pivoting.
//
  info = 0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Find L = pivot index.
//
    l = idamax ( n-k+1, a+(k-1)+(k-1)*lda, 1 ) + k - 1;
    ipvt[k-1] = l;
//
//  Zero pivot implies this column already triangularized.
//
    if ( a[l-1+(k-1)*lda] == 0.0 )
    {
      info = k;
      continue;
    }
//
//  Interchange if necessary.
//
    if ( l != k )
    {
      t = a[l-1+(k-1)*lda];
      a[l-1+(k-1)*lda] = a[k-1+(k-1)*lda];
      a[k-1+(k-1)*lda] = t;
    }
//
//  Compute multipliers.
//
    t = -1.0 / a[k-1+(k-1)*lda];

    dscal ( n-k, t, a+k+(k-1)*lda, 1 );
//
//  Row elimination with column indexing.
//
    for ( j = k+1; j <= n; j++ )
    {
      t = a[l-1+(j-1)*lda];
      if ( l != k )
      {
        a[l-1+(j-1)*lda] = a[k-1+(j-1)*lda];
        a[k-1+(j-1)*lda] = t;
      }
      daxpy ( n-k, t, a+k+(k-1)*lda, 1, a+k+(j-1)*lda, 1 );
    }

  }

  ipvt[n-1] = n;

  if ( a[n-1+(n-1)*lda] == 0.0 )
  {
    info = n;
  }

  return info;
}
//****************************************************************************80

void dgesl ( double a[], int lda, int n, int ipvt[], double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DGESL solves a real general linear system A * X = B.
//
//  Discussion:
//
//    DGESL can solve either of the systems A * X = B or A' * X = B.
//
//    The system matrix must have been factored by DGECO or DGEFA.
//
//    A division by zero will occur if the input factor contains a
//    zero on the diagonal.  Technically this indicates singularity
//    but it is often caused by improper arguments or improper
//    setting of LDA.  It will not occur if the subroutines are
//    called correctly and if DGECO has set 0.0 < RCOND
//    or DGEFA has set INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
//    Pete Stewart.
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
//    Input, double A[LDA*N], the output from DGECO or DGEFA.
//
//    Input, int LDA, the leading dimension of A.
//
//    Input, int N, the order of the matrix A.
//
//    Input, int IPVT[N], the pivot vector from DGECO or DGEFA.
//
//    Input/output, double B[N].
//    On input, the right hand side vector.
//    On output, the solution vector.
//
//    Input, int JOB.
//    0, solve A * X = B;
//    nonzero, solve A' * X = B.
//
{
  int k;
  int l;
  double t;
//
//  Solve A * X = B.
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

      daxpy ( n-k, t, a+k+(k-1)*lda, 1, b+k, 1 );

    }

    for ( k = n; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
      t = -b[k-1];
      daxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
    for ( k = 1; k <= n; k++ )
    {
      t = ddot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
      b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
    }

    for ( k = n-1; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] + ddot ( n-k, a+k+(k-1)*lda, 1, b+k, 1 );
      l = ipvt[k-1];

      if ( l != k )
      {
        t = b[l-1];
        b[l-1] = b[k-1];
        b[k-1] = t;
      }
    }
  }
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
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
//    David Kincaid, Fred Krogh.
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

int idamax ( int n, double dx[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    IDAMAX finds the index of the vector element of maximum absolute value.
//
//  Discussion:
//
//    WARNING: This index is a 1-based index, not a 0-based index!
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
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
//    David Kincaid, Fred Krogh.
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
//    Input, double X[*], the vector to be examined.
//
//    Input, int INCX, the increment between successive entries of SX.
//
//    Output, int IDAMAX, the index of the element of maximum
//    absolute value.
//
{
  double dmax;
  int i;
  int ix;
  int value;

  value = 0;

  if ( n < 1 || incx <= 0 )
  {
    return value;
  }

  value = 1;

  if ( n == 1 )
  {
    return value;
  }

  if ( incx == 1 )
  {
    dmax = fabs ( dx[0] );

    for ( i = 1; i < n; i++ )
    {
      if ( dmax < fabs ( dx[i] ) )
      {
        value = i + 1;
        dmax = fabs ( dx[i] );
      }
    }
  }
  else
  {
    ix = 0;
    dmax = fabs ( dx[0] );
    ix = ix + incx;

    for ( i = 1; i < n; i++ )
    {
      if ( dmax < fabs ( dx[ix] ) )
      {
        value = i + 1;
        dmax = fabs ( dx[ix] );
      }
      ix = ix + incx;
    }
  }

  return value;
}
//****************************************************************************80

double *mv_gb ( int m, int n, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MV_GB multiplies a banded matrix by an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2014
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
//    Input, int ML, MU, the lower and upper bandwidths.
//
//    Input, double A[(2*ML+MU+1)*N], the matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double MV_GB[M], the product A * x.
//
{
  double *b;
  int i;
  int j;
  int jhi;
  int jlo;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
  }

  for ( i = 0; i < m; i++ )
  {
    jlo = i4_max ( 0, i - ml );
    jhi = i4_min ( n - 1, i + mu );
    for ( j = jlo; j <= jhi; j++ )
    {
      b[i] = b[i] + a[i-j+ml+mu+j*(2*ml+mu+1)] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

double *mv_ge ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MV_GE multiplies a GE matrix by an R8VEC.
//
//  Discussion:
//
//    The GE storage format is used for a general M by N matrix.  A storage 
//    space is made for each entry.  The two dimensional logical
//    array can be thought of as a vector of M*N entries, starting with
//    the M entries in the column 1, then the M entries in column 2
//    and so on.  Considered as a vector, the entry A(I,J) is then stored
//    in vector location I+(J-1)*M.
//
//    GE storage is used by LINPACK and LAPACK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2014
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
//    Input, double A(M,N), the matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double MV_GE[M], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

double *mv_st ( int m, int n, int nz_num, int row[], int col[], double a[], 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MV_ST multiplies a sparse triple matrix times a vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int NZ_NUM, the number of nonzero values.
//
//    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and 
//    column indices.
//
//    Input, double A[NZ_NUM], the nonzero values in the matrix.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Output, double MV_ST[M], the product A*X.
//
{
  double *b;
  int i;
  int k;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
  }

  for ( k = 0; k < nz_num; k++ )
  {
    b[row[k]] = b[row[k]] + a[k] * x[col[k]];
  }

  return b;
}
//****************************************************************************80

int nonzeros ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    NONZEROS counts the nonzeros in a matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], the matrix.
//
//    Output, int NONZEROS, the number of nonzero entries.
//
{
  int i;
  int j;
  int nnz;

  nnz = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] != 0.0 )
      {
        nnz = nnz + 1;
      }
    }
  }

  return nnz;
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

double *r8mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
//    in column-major order.
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
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
//    03 October 2005
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
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }
      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

double r8vec_diff_norm_li ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIFF_NORM_LI returns the L-oo norm of the difference of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L-oo norm is defined as:
//
//      R8VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], B[N], the vectors.
//
//    Output, double R8VEC_DIFF_NORM_LI, the L-oo norm of A - B.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = r8_max ( value, fabs ( a[i] - b[i] ) );
  }
  return value;
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

void wathen_bandwidth ( int nx, int ny, int &l, int &d, int &u )

//****************************************************************************80
//
//  Purpose:
//
//    WATHEN_BANDWIDTH returns the bandwidth of the WATHEN matrix.
//
//  Discussion:
//
//    The bandwidth measures the minimal number of contiguous diagonals,
//    including the central diagonal, which contain all the nonzero elements
//    of a matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
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
//    Output, int &L, &D, &U, the lower, diagonal, and upper 
//    bandwidths of the matrix,
//
{
  l = 3 * nx + 4;
  d = 1;
  u = 3 * nx + 4;

  return;
}
//****************************************************************************80

double *wathen_gb ( int nx, int ny, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    WATHEN_GB returns the Wathen matrix, using general banded (GB) storage.
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
//    The matrix is the consistent mass matrix for a regular NX by NY grid
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
//    The matrix is symmetric positive definite for any positive values of the
//    density RHO(X,Y).
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
//    Input, int NX, NY, values which determine the size
//    of the matrix.
//
//    Input, int N, the number of rows and columns.
//
//    Input/output, int &SEED, the random number seed.
//
//    Output, double WATHEN_GB[(9*NX+13)*N], the matrix.
//
{
  double *a;
  const double em[8*8] = {
     6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, 
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, 
     2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, 
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, 
     3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, 
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, 
     2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, 
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 };
  int i;
  int ii;
  int j;
  int jj;
  int kcol;
  int krow;
  int lda;
  int ml;
  int mu;
  int node[8];
  double rho;

  ml = 3 * nx + 4;
  mu = 3 * nx + 4;
  lda = 2 * ml + mu + 1;
  a = new double[lda*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < lda; i++ )
    {
      a[i+j*lda] = 0.0;
    }
  }

  for ( j = 0; j < nx; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      node[0] = 3 * ( j + 1 ) * nx + 2 * ( j + 1 ) + 2 * ( i + 1 );
      node[1] = node[0] - 1;
      node[2] = node[0] - 2;
      node[3] = ( 3 * ( j + 1 ) - 1 ) * nx + 2 * ( j + 1 ) + ( i + 1 ) - 2;
      node[4] = ( 3 * ( j + 1 ) - 3 ) * nx + 2 * ( j + 1 ) + 2 * ( i + 1 ) - 4;
      node[5] = node[4] + 1;
      node[6] = node[4] + 2;
      node[7] = node[3] + 1;

      rho = 100.0 * r8_uniform_01 ( seed );

      for ( krow = 0; krow < 8; krow++ )
      {
        ii = node[krow];
        for ( kcol = 0; kcol < 8; kcol++ )
        {
          jj = node[kcol];
          a[ii-jj+ml+mu+jj*lda] = a[ii-jj+ml+mu+jj*lda] 
            + rho[i+j*nx] * em[krow+kcol*8];
        }
      }
    }
  }

  return a;
}
//****************************************************************************80

double *wathen_ge ( int nx, int ny, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    WATHEN_GE returns the Wathen matrix as a general storage (GE) matrix.
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
//    The matrix is the consistent mass matrix for a regular NX by NY grid
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
//    The matrix is symmetric positive definite for any positive values of the
//    density RHO(X,Y).
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
//    Input, int NX, NY, values which determine the size 
//    of the matrix.
//
//    Input, int N, the number of rows and columns.
//
//    Input/output, int &SEED, the random number seed.
//
//    Output, double WATHEN_GE[N*N], the matrix.
//
{
  double *a;
  const double em[8*8] = {
     6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, 
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, 
     2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, 
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, 
     3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, 
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, 
     2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, 
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 };
  int i;
  int ii;
  int j;
  int jj;
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

  for ( j = 0; j < nx; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      node[0] = 3 * ( j + 1 ) * nx + 2 * ( j + 1 ) + 2 * ( i + 1 );
      node[1] = node[0] - 1;
      node[2] = node[0] - 2;
      node[3] = ( 3 * ( j + 1 ) - 1 ) * nx + 2 * ( j + 1 ) + ( i + 1 ) - 2;
      node[4] = ( 3 * ( j + 1 ) - 3 ) * nx + 2 * ( j + 1 ) + 2 * ( i + 1 ) - 4;
      node[5] = node[4] + 1;
      node[6] = node[4] + 2;
      node[7] = node[3] + 1;

      rho = 100.0 * r8_uniform_01 ( seed );

      for ( krow = 0; krow < 8; krow++ )
      {
        ii = node[krow];
        for ( kcol = 0; kcol < 8; kcol++ )
        {
          jj = node[kcol];
          a[ii+jj*n] = a[ii+jj*n] + rho * em[krow+kcol*8];
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
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
//    Output, int WATHEN_ORDER, the order of the matrix,
//    as determined by NX and NY.
//
{
  int n;

  n = 3 * nx * ny + 2 * nx + 2 * ny + 1;

  return n;
}
//****************************************************************************80

double *wathen_st ( int nx, int ny, int nz_num, int &seed, int row[], 
  int col[] )

//****************************************************************************80
//
//  Purpose:
//
//    WATHEN_ST: Wathen matrix stored in sparse triplet (ST) format.
//
//  Discussion:
//
//    When dealing with sparse matrices in MATLAB, it can be much more efficient
//    to work first with a triple of I, J, and X vectors, and only once
//    they are complete, convert to MATLAB's sparse format.
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
//
//      N = 3*NX*NY + 2*NX + 2*NY + 1,
//
//    The matrix is the consistent mass matrix for a regular NX by NY grid
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
//    The matrix is symmetric positive definite for any positive values of the
//    density RHO(X,Y).
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
//    John Burkardt.
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
//    Input, int NX, NY, values which determine the size of 
//    the matrix.
//
//    Input, int NZ_NUM, the number of values used to 
//    describe the matrix.
//
//    Input/output, int &SEED, the random number seed.
//
//    Output, int ROW[NZ_NUM], COL[NZ_NUM], the row and 
//    column indices of the nonzero entries.
//
//    Output, double WATHEN_ST[NZ_NUM], the nonzero entries of the matrix.
//
{
  double *a;
  const double em[8*8] = {
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
  int k;
  int kcol;
  int krow;
  int node[8];
  double rho;

  a = new double[nz_num];
 
  for ( k = 0; k < nz_num; k++ )
  {
    row[k] = 0;
    col[k] = 0;
    a[k] = 0.0;
  }

  k = 0;

  for ( j = 0; j < nx; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      node[0] = 3 * ( j + 1 ) * nx + 2 * ( j + 1 ) + 2 * ( i + 1 );
      node[1] = node[0] - 1;
      node[2] = node[0] - 2;
      node[3] = ( 3 * ( j + 1 ) - 1 ) * nx + 2 * ( j + 1 ) + ( i + 1 ) - 2;
      node[4] = ( 3 * ( j + 1 ) - 3 ) * nx + 2 * ( j + 1 ) + 2 * ( i + 1 ) - 4;
      node[5] = node[4] + 1;
      node[6] = node[4] + 2;
      node[7] = node[3] + 1;

      rho = 100.0 * r8_uniform_01 ( seed );

      for ( krow = 0; krow < 8; krow++ )
      {
        for ( kcol = 0; kcol < 8; kcol++ )
        {
          row[k] = node[krow];
          col[k] = node[kcol];
          a[k] = rho * em[krow+kcol*8];
          k = k + 1;
        }
      }
    }
  }

  return a;
}
//****************************************************************************80

int wathen_st_size ( int nx, int ny )

//****************************************************************************80
//
//  Purpose:
//
//    WATHEN_ST_SIZE: Size of Wathen matrix stored in sparse triplet format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2014
//
//  Author:
//
//    John Burkardt.
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
//    Input, integer NX, NY, values which determine the size of the matrix.
//
//    Output, integer NZ_NUM, the number of items of data used to describe
//    the matrix.
//
{
  int nz_num;

  nz_num = nx * ny * 64;

  return nz_num;
}
