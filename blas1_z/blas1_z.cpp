# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "blas0.hpp"
# include "blas1_z.hpp"

//****************************************************************************80

double dzasum ( int n, complex <double> x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DZASUM takes the sum of the absolute values of a complex <double> vector.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, complex <double> X[], the vector.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, double DZASUM, the sum of the absolute values.
//
{
  int i;
  int ix;
  double value;

  value = 0.0;

  if ( n <= 0 || incx <= 0 )
  {
    return value;
  }

  if ( incx == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      value = value + fabs ( real ( x[i] ) )
                    + fabs ( imag ( x[i] ) );
    }
  }
  else
  {
    ix = 0;
    for ( i = 0; i < n; i++ )
    {
      value = value + fabs ( real ( x[ix] ) )
                    + fabs ( imag ( x[ix] ) );
      ix = ix + incx;
    }
  }
  return value;
}
//****************************************************************************80

double dznrm2 ( int n, complex <double> x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DZNRM2 returns the euclidean norm of a complex <double> vector.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    DZNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
//            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, complex <double> X[], the vector.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, double DZNRM2, the norm of the vector.
//
{
  int i;
  int ix;
  double scale;
  double ssq;
  double temp;
  double value;

  if ( n < 1 || incx < 1 )
  {
    value  = 0.0;
  }
  else
  {
    scale = 0.0;
    ssq = 1.0;
    ix = 0;

    for ( i = 0; i < n; i++ )
    {
      if ( real ( x[ix] ) != 0.0 )
      {
        temp = fabs ( real ( x[ix] ) );
        if ( scale < temp )
        {
          ssq = 1.0 + ssq * pow ( scale / temp, 2 );
          scale = temp;
        }
        else
        {
          ssq = ssq + pow ( temp / scale, 2 );
        }
      }

      if ( imag ( x[ix] ) != 0.0 )
      {
        temp = fabs ( imag ( x[ix] ) );
        if ( scale < temp )
        {
          ssq = 1.0 + ssq * pow ( scale / temp, 2 );
          scale = temp;
        }
        else
        {
          ssq = ssq + pow ( temp / scale, 2 );
        }
      }
      ix = ix + incx;
    }
    value  = scale * sqrt ( ssq );
  }
  return value;
}
//****************************************************************************80

int izamax ( int n, complex <double> x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    IZAMAX indexes the complex <double> vector element of maximum absolute value.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    WARNING: This index is a 1-based index, not a 0-based index!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, complex <double> X[], the vector.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, int IZAMAX, the index of the element of maximum
//    absolute value.
//
{
  int i;
  int ix;
  double smax;
  int value;

  value = 0;

  if ( n < 1 || incx  <=  0 )
  {
    return value;
  }

  value = 1;

  if ( n == 1 )
  {
    return value;
  }

  if ( incx != 1 )
  {
    ix = 0;
    smax = zabs1 ( x[0] );
    ix = ix + incx;

    for ( i = 1; i < n; i++ )
    {
      if ( smax < zabs1 ( x[ix] ) )
      {
        value = i + 1;
        smax = zabs1 ( x[ix] );
      }
      ix = ix + incx;
    }
  }
  else
  {
    smax = zabs1 ( x[0] );
    for ( i = 1; i < n; i++ )
    {
      if ( smax < zabs1 ( x[i] ) )
      {
        value = i + 1;
        smax = zabs1 ( x[i] );
      }
    }
  }

  return value;
}
//****************************************************************************80

void zaxpy ( int n, complex <double> ca, complex <double> cx[],
  int incx, complex <double> cy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    ZAXPY adds a multiple of one complex <double> vector to another.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Input, int N, the number of elements in CX and CY.
//
//    Input, complex <double> CA, the multiplier of CX.
//
//    Input, complex <double> CX[], the first vector.
//
//    Input, int INCX, the increment between successive entries of CX.
//
//    Input/output, complex <double> CY[], the second vector.
//    On output, CY(*) has been replaced by CY(*) + CA * CX(*).
//
//    Input, int INCY, the increment between successive entries of CY.
//
{
  int i;
  int ix;
  int iy;

  if ( n <= 0 )
  {
    return;
  }

  if ( zabs1 ( ca ) == 0.0 )
  {
    return;
  }

  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      cy[iy] = cy[iy] + ca * cx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      cy[i] = cy[i] + ca * cx[i];
    }

  }

  return;
}
//****************************************************************************80

void zcopy ( int n, complex <double> cx[], int incx, complex <double> cy[],
  int incy )

//****************************************************************************80
//
//  Purpose:
//
//    ZCOPY copies a complex <double> vector.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Input, int N, the number of elements in CX and CY.
//
//    Input, complex <double> CX[], the first vector.
//
//    Input, int INCX, the increment between successive entries of CX.
//
//    Output, complex <double> CY[], the second vector.
//
//    Input, int INCY, the increment between successive entries of CY.
//
{
  int i;
  int ix;
  int iy;

  if ( n <= 0 )
  {
    return;
  }

  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      cy[iy] = cx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      cy[i] = cx[i];
    }
  }
  return;
}
//****************************************************************************80

complex <double> zdotc ( int n, complex <double> cx[], int incx,
  complex <double> cy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    ZDOTC forms the conjugated dot product of two complex <double> vectors.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Input, int N, the number of entries in the vectors.
//
//    Input, complex <double> CX[], the first vector.
//
//    Input, int INCX, the increment between successive entries in CX.
//
//    Input, complex <double> CY[], the second vector.
//
//    Input, int INCY, the increment between successive entries in CY.
//
//    Output, complex <double> ZDOTC, the conjugated dot product of
//    the corresponding entries of CX and CY.
//
{
  int i;
  int ix;
  int iy;
  complex <double> value;

  value = complex <double> ( 0.0, 0.0 );

  if ( n <= 0 )
  {
    return value;
  }

  if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      value = value + conj ( cx[i] ) * cy[i];
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
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      value = value + conj ( cx[ix] ) * cy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  return value;
}
//****************************************************************************80

complex <double> zdotu ( int n, complex <double> cx[], int incx,
  complex <double> cy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    ZDOTU forms the unconjugated dot product of two complex <double> vectors.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Input, int N, the number of entries in the vectors.
//
//    Input, complex <double> CX[], the first vector.
//
//    Input, int INCX, the increment between successive entries in CX.
//
//    Input, complex <double> CY[], the second vector.
//
//    Input, int INCY, the increment between successive entries in CY.
//
//    Output, complex <double> ZDOTU, the unconjugated dot product of
//    the corresponding entries of CX and CY.
//
{
  int i;
  int ix;
  int iy;
  complex <double> value;

  value = complex <double> ( 0.0, 0.0 );

  if ( n <= 0 )
  {
    return value;
  }

  if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      value = value + cx[i] * cy[i];
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
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      value = value + cx[ix] * cy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  return value;
}
//****************************************************************************80

void zdrot ( int n, complex <double> cx[], int incx, complex <double> cy[],
  int incy, double c, double s )

//****************************************************************************80
//
//  Purpose:
//
//    ZDROT applies a plane rotation to complex <double> vectors.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The cosine C and sine S are real and the vectors CX and CY are complex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, complex <double> CX[], one of the vectors to be rotated.
//
//    Input, int INCX, the increment between successive entries of CX.
//
//    Input/output, complex <double> CY[], one of the vectors to be rotated.
//
//    Input, int INCY, the increment between successive elements of CY.
//
//    Input, double C, S, parameters (presumably the cosine and sine of
//    some angle) that define a plane rotation.
//
{
  complex <double> ctemp;
  int i;
  int ix;
  int iy;

  if ( n <= 0 )
  {
    return;
  }

  if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      ctemp = c * cx[i] + s * cy[i];
      cy[i] = c * cy[i] - s * cx[i];
      cx[i] = ctemp;
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
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      ctemp  = c * cx[ix] + s * cy[iy];
      cy[iy] = c * cy[iy] - s * cx[ix];
      cx[ix] = ctemp;
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  return;
}
//****************************************************************************80

void zdscal ( int n, double sa, complex <double> cx[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    ZDSCAL scales a complex <double> vector by a double.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The scaling constant is double precision real.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Input, int N, the number of entries in the vector.
//
//    Input, double SA, the multiplier.
//
//    Input/output, complex <double> CX[], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of
//    the vector CX.
//
{
  int i;

  if ( n <= 0 || incx <= 0 )
  {
    return;
  }

  if ( incx == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      cx[i] = sa * cx[i];
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      cx[i*incx] = sa * cx[i*incx];
    }
  }
  return;
}
//****************************************************************************80

void zrotg ( complex <double> *ca, complex <double> cb, double *c,
  complex <double> *s )

//****************************************************************************80
//
//  Purpose:
//
//    ZROTG determines a complex <double> Givens rotation.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    Given values A and B, this routine computes:
//
//    If A = 0:
//
//      R = B
//      C = 0
//      S = (1,0).
//
//    If A /= 0:
//
//      ALPHA = A / abs ( A )
//      NORM  = sqrt ( ( abs ( A ) )^2 + ( abs ( B ) )^2 )
//      R     = ALPHA * NORM
//      C     = abs ( A ) / NORM
//      S     = ALPHA * conj ( B ) / NORM
//
//    In either case, the computed numbers satisfy the equation:
//
//    (         C    S ) * ( A ) = ( R )
//    ( -conj ( S )  C )   ( B ) = ( 0 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2007
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
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input/output, complex <double> *CA, on input, the value A.  On output,
//    the value R.
//
//    Input, complex <double> CB, the value B.
//
//    Output, double *C, the cosine of the Givens rotation.
//
//    Output, complex <double> *S, the sine of the Givens rotation.
//
{
  complex <double> alpha;
  double norm;
  double scale;

  if ( zabs2 ( *ca ) == 0.0 )
  {
    *c = 0.0;
    *s = complex <double> ( 1.0, 0.0 );
    *ca = cb;
  }
  else
  {
    scale = zabs2 ( *ca ) + zabs2 ( cb );
    norm = scale * sqrt ( pow ( zabs2 ( *ca / scale ), 2 )
                        + pow ( zabs2 (  cb / scale ), 2 ) );
    alpha = *ca / zabs2 ( *ca );
    *c = zabs2 ( *ca ) / norm;
    *s = alpha * conj ( cb ) / norm;
    *ca = alpha * norm;
  }

  return;
}
//****************************************************************************80

void zscal ( int n, complex <double> ca, complex <double> cx[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    ZSCAL scales a complex <double> vector by a complex <double>.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    C++ version by John Burkardt.
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
//    Input, int N, the number of entries in the vector.
//
//    Input, complex <double> CA, the multiplier.
//
//    Input/output, complex <double> CX[], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of CX.
//
{
  int i;

  if ( n <= 0 || incx <= 0 )
  {
    return;
  }

  if ( incx == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      cx[i] = ca * cx[i];
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      cx[i*incx] = ca * cx[i*incx];
    }
  }
  return;
}
//****************************************************************************80

void zswap ( int n, complex <double> cx[], int incx, complex <double> cy[],
  int incy )

//****************************************************************************80
//
//  Purpose:
//
//    ZSWAP interchanges two complex <double> vectors.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
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
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, complex <double> CX[], one of the vectors to swap.
//
//    Input, int INCX, the increment between successive entries of CX.
//
//    Input/output, complex <double> CY[], one of the vectors to swap.
//
//    Input, int INCY, the increment between successive elements of CY.
//
{
  complex <double> ctemp;
  int i;
  int ix;
  int iy;

  if ( n <= 0 )
  {
    return;
  }

  if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      ctemp = cx[i];
      cx[i] = cy[i];
      cy[i] = ctemp;
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
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      ctemp = cx[ix];
      cx[ix] = cy[iy];
      cy[iy] = ctemp;
      ix = ix + incx;
      iy = iy + incy;
    }
  }

  return;
}
