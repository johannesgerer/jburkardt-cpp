# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa007.hpp"

//****************************************************************************80

void cholesky ( double a[], int n, int nn, double u[], int *nullty, 
  int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    CHOLESKY computes the Cholesky factorization of a PDS matrix.
//
//  Discussion:
//
//    For a positive definite symmetric matrix A, the Cholesky factor U
//    is an upper triangular matrix such that A = U' * U.
//
//    This routine was originally named "CHOL", but that conflicted with
//    a built in MATLAB routine name.
//
//    The missing initialization "II = 0" has been added to the code.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Michael Healy.
//    Modifications by AJ Miller.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Michael Healy,
//    Algorithm AS 6:
//    Triangular decomposition of a symmetric matrix,
//    Applied Statistics,
//    Volume 17, Number 2, 1968, pages 195-197.
//
//  Parameters:
//
//    Input, double A((N*(N+1))/2), a positive definite matrix 
//    stored by rows in lower triangular form as a one dimensional array, 
//    in the sequence
//    A(1,1),
//    A(2,1), A(2,2),
//    A(3,1), A(3,2), A(3,3), and so on.
//
//    Input, int N, the order of A.
//
//    Input, int NN, the dimension of the array used to store A, 
//    which should be at least (N*(N+1))/2.
//
//    Output, double U((N*(N+1))/2), an upper triangular matrix,
//    stored by columns, which is the Cholesky factor of A.  The program is
//    written in such a way that A and U can share storage.
//
//    Output, int NULLTY, the rank deficiency of A.  If NULLTY is zero,
//    the matrix is judged to have full rank.
//
//    Output, int IFAULT, an error indicator.
//    0, no error was detected;
//    1, if N < 1;
//    2, if A is not positive semi-definite.
//    3, if NN < (N*(N+1))/2.
//
//  Local Parameters:
//
//    Local, double ETA, should be set equal to the smallest positive
//    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
//    accuracy being used.
//
{
  double eta = 1.0E-09;
  int i;
  int icol;
  int ii;
  int irow;
  int j;
  int k;
  int kk;
  int l;
  int m;
  double w;
  double x;

  *ifault = 0;
  *nullty = 0;

  if ( n <= 0 )
  {
    *ifault = 1;
    return;
  }

  if ( nn < ( n * ( n + 1 ) ) / 2 )
  {
    *ifault = 3;
    return;
  }

  j = 1;
  k = 0;
  ii = 0;
//
//  Factorize column by column, ICOL = column number.
//
  for ( icol = 1; icol <= n; icol++ )
  {
    ii = ii + icol;
    x = eta * eta * a[ii-1];
    l = 0;
    kk = 0;
//
//  IROW = row number within column ICOL.
//
    for ( irow = 1; irow <= icol; irow++ )
    {
      kk = kk + irow;
      k = k + 1;
      w = a[k-1];
      m = j;

      for ( i = 1; i < irow; i++ )
      {
        l = l + 1;
        w = w - u[l-1] * u[m-1];
        m = m + 1;
      }

      l = l + 1;

      if ( irow == icol ) 
      {
        break;
      }

      if ( u[l-1] != 0.0 )
      {
        u[k-1] = w / u[l-1];
      }
      else
      {
        u[k-1] = 0.0;

        if ( r8_abs ( x * a[k-1] ) < w * w )
        {
          *ifault = 2;
          return;
        }
      }
    }
//
//  End of row, estimate relative accuracy of diagonal element.
//
    if ( r8_abs ( w ) <= r8_abs ( eta * a[k-1] ) )
    {
      u[k-1] = 0.0;
      *nullty = *nullty + 1;
    }
    else
    {
      if ( w < 0.0 )
      {
        *ifault = 2;
        return;
      }
      u[k-1] = sqrt ( w );
    }
    j = j + icol;
  }
  return;
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
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

void syminv ( double a[], int n, double c[], double w[], int *nullty, 
  int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    SYMINV computes the inverse of a symmetric matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Michael Healy.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Michael Healy,
//    Algorithm AS 7:
//    Inversion of a Positive Semi-Definite Symmetric Matrix,
//    Applied Statistics,
//    Volume 17, Number 2, 1968, pages 198-199.
//
//  Parameters:
//
//    Input, double A((N*(N+1))/2), a positive definite matrix stored 
//    by rows in lower triangular form as a one dimensional array, in the sequence
//    A(1,1),
//    A(2,1), A(2,2),
//    A(3,1), A(3,2), A(3,3), and so on.
//
//    Input, int N, the order of A.
//
//    Output, double C((N*(N+1))/2), the inverse of A, or generalized
//    inverse if A is singular, stored using the same storage scheme employed
//    for A.  The program is written in such a way that A and U can share storage.
//
//    Workspace, double W(N).
//
//    Output, int *NULLTY, the rank deficiency of A.  If NULLTY is zero,
//    the matrix is judged to have full rank.
//
//    Output, int *IFAULT, error indicator.
//    0, no error detected.
//    1, N < 1.
//    2, A is not positive semi-definite.
//
{
  int i;
  int icol;
  int irow;
  int j;
  int jcol;
  int k;
  int l;
  int mdiag;
  int ndiag;
  int nn;
  int nrow;
  double x;

  *ifault = 0;

  if ( n <= 0 )
  {
    *ifault = 1;
    return;
  }

  nrow = n;
//
//  Compute the Cholesky factorization of A.
//  The result is stored in C.
//
  nn = ( n * ( n + 1 ) ) / 2;

  cholesky ( a, n, nn, c, nullty, ifault );

  if ( *ifault != 0 )
  {
    return;
  }
//
//  Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse
//  of C, row by row starting with the last row.
//  IROW = the row number, 
//  NDIAG = location of last element in the row.
//
  irow = nrow;
  ndiag = nn;
//
//  Special case, zero diagonal element.
//
  for ( ; ; )
  {
    if ( c[ndiag-1] == 0.0 )
    {
      l = ndiag;
      for ( j = irow; j <= nrow; j++ )
      {
        c[l-1] = 0.0;
        l = l + j;
      }
    }
    else
    {
      l = ndiag;
      for ( i = irow; i <= nrow; i++ )
      {
        w[i-1] = c[l-1];
        l = l + i;
      }

      icol = nrow;
      jcol = nn;
      mdiag = nn;

      for ( ; ; )
      {
        l = jcol;
 
        if ( icol == irow )
        {
          x = 1.0 / w[irow-1];
        }
        else
        {
          x = 0.0;
        }

        k = nrow;

        while ( irow < k )
        {
          x = x - w[k-1] * c[l-1];
          k = k - 1;
          l = l - 1;

          if ( mdiag < l )
          {
            l = l - k + 1;
          }
        }

        c[l-1] = x / w[irow-1];

        if ( icol <= irow )
        {
          break;
        }
        mdiag = mdiag - icol;
        icol = icol - 1;
        jcol = jcol - 1;
      }
    }

    ndiag = ndiag - irow;
    irow = irow - 1;

    if ( irow <= 0 )
    {
      break;
    }
  }

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
//    24 September 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
