# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "asa314.hpp"

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

void i4mat_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT prints an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
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
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT_SOME prints some of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
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
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 10

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
//  Print the columns of the matrix, in strips of INCX.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << "  " << setw(6) << j - 1;
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to INCX) entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ":";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << "  " << setw(6) << a[i-1+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void invmod ( int mat[], int imat[], int rmod[], int cmod[], int nrow, 
  int &ifault )

//****************************************************************************80
//
//  Purpose:
//
//    INVMOD inverts a matrix using modulo arithmetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2013
//
//  Author:
//
//    Original FORTRAN77 version by Roger Payne.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Payne,
//    Inversion of matrices with contents subject to modulo arithmetic,
//    Applied Statistics,
//    Volume 46, Number 2, 1997, pages 295-298.
//
//  Parameters:
//
//    Input/output, int MAT[NROW*NROW].
//    On input, the matrix to be inverted.
//    On output, the product of the input matrix and IMAT.
//
//    Output, int IMAT[NROW*NROW], the inverse matrix.  
//    If IFAULT = -1 on output, then IMAT is only a left inverse.
//
//    Input, int RMOD[NROW], the modulus for values in each row.
//
//    Input, int CMOD[NROW], the modulus for values 
//    in each column.
//
//    Input, int NROW, the order of the matrix.
//
//    Output, int &IFAULT, an error flag.
//    0, no error was detected.
//    -1, only a left inverse could be formed.
//    1, the matrix contains elements that are negative, or too large.
//    2, the matrix contains nonzero elements in mixed modulus positions.
//    3, the matrix cannot be inverted.
//
{
  int all_zero;
  int *csort;
  int i;
  int ir;
  int j;
  int k;
  int kir;
  int kjr;
  int n;
  int *rsort;
//
//  Check that elements in 'mixed-moduli' positions are all zero.
//
  n = 0;
  for ( i = 1; i <= nrow; i++ )
  {
    for ( j = 1; j <= nrow; j++ )
    {
      n = n + 1;

      if ( ( rmod[i-1] != cmod[j-1] ) && ( 0 < mat[n-1] ) )
      {
        ifault = 2;
        return;
      }

      if ( ( rmod[i-1] < mat[n-1] ) ||( mat[n-1] < 0 ) )
      {
        ifault = 1;
        return;
      }
    }
  }

  n = 0;
  for ( i = 1; i <= nrow; i++ )
  {
    for ( j = 1; j <= nrow; j++ )
    {
      n = n + 1;
      imat[n-1] = 0;
    }
  }
//
//  Sort rows and columns into ascending order of moduli
//
  rsort = new int[nrow];
  csort = new int[nrow];

  msort ( mat, imat, rmod, cmod, rsort, csort, nrow );
//
//  Complete initialization of inverse matrix 
//
  for ( n = 1; n <= nrow * nrow; n = n + nrow + 1 )
  {
    imat[n-1] = 1;
  }
//
//  Invert the matrix.
//
  for ( ir = 1; ir <= nrow; ir++ )
  {
    kir = ( ir - 1 ) * nrow;

    if ( mat[kir+ir-1] == 0 )
    {
//
//  Find a row JR below IR such that K(JR,IR)>0
//
      all_zero = 1;

      for ( kjr = kir + nrow + ir; kjr <= nrow * nrow; kjr = kjr + nrow )
      {
        if ( 0 < mat[kjr-1] )
        {
          all_zero = 0;
          break;
        }
      }
//
//  Column IR contains all zeros in rows IR or below:
//  look for a row above with zeros to left of column IR 
//  and K(JR,IR)>0
//
      if ( all_zero )
      {
        for ( kjr = ir; kjr <= kir; kjr = kjr + nrow )
        {
          if ( 0 < mat[kjr-1] )
          {
            for ( i = kjr - ir + 1; i < kjr; i++ )
            {
              if ( 0 < mat[i-1] )
              {
                ifault = 3;
                return;
              }
            }
            all_zero = 0;
            break;
          }
        }
      }
//
//  Column IR contains all zeros
//
      if ( all_zero )
      {
        continue;
      }
//
//  Switch row JR with row IR
//
      kjr = kjr - ir;

      for ( i = 1; i <= nrow; i++ )
      {
        k = mat[kir+i-1];
        mat[kir+i-1] = mat[kjr+i-1];
        mat[kjr+i-1] = k;

        k = imat[kir+i-1];
        imat[kir+i-1] = imat[kjr+i-1];
        imat[kjr+i-1] = k;
      }
    }
//
//  Find a multiplier N such that N*MAT(IR,IR)=1 mod(P{IR})
//
    k = mat[kir+ir-1];
    for ( n = 1; n < rmod[ir-1]; n++ )
    {
      if ( ( n * k ) % rmod[ir-1] == 1 )
      {
        break;
      }
    }
//
//  Multiply row IR by N.
//
    if ( 1 < n )
    {
      for ( i = kir + 1; i <= ir * nrow; i++ )
      {
        mat[i-1] = mat[i-1] * n;
        imat[i-1] = imat[i-1] * n;
      }
    }
//
//  Subtract MAT(JR,IR) * row IR from each row JR
//
    for ( kjr = 0; kjr < nrow * nrow; kjr = kjr + nrow )
    {
      n = rmod[ir-1] - mat[kjr+ir-1];
      if ( ( kjr != kir ) && ( n != 0 ) )
      {
        for ( i = 1; i <= nrow; i++ )
        {
          mat[kjr+i-1]  = (  mat[kjr+i-1] + n *  mat[kir+i-1] ) % cmod[i-1];
          imat[kjr+i-1] = ( imat[kjr+i-1] + n * imat[kir+i-1] ) % cmod[i-1];
        }
      }
    }

  }
//
//  Check inversion was possible - that result has
//  non-zero elements only on diagonal.
//
  ifault = 0;
//
//  If we encounter a zero diagonal element, then only a left inverse
//  will be formed.
//
  for ( n = 1; n <= nrow * nrow; n = n + nrow + 1 )
  {
    if ( mat[n-1] == 0 )
    {
      ifault = -1;
    }
    mat[n-1] = - mat[n-1];
  }

  for ( n = 1; n <= nrow * nrow; n++ )
  {
    if ( 0 < mat[n-1] )
    {
      ifault = 3;
      return;
    }
  }

  for ( n = 1; n <= nrow * nrow; n = n + nrow + 1 )
  {
    mat[n-1] = - mat[n-1];
  }
//
//  Unsort the rows and columns back into their original order.
//
  musort ( mat, imat, rmod, cmod, rsort, csort, nrow );

  delete [] csort;
  delete [] rsort;

  return;
}
//****************************************************************************80

void msort ( int mat[], int imat[], int rmod[], int cmod[], int rsort[], 
  int csort[], int nrow )

//****************************************************************************80
//
//  Purpose:
//
//    MSORT sorts matrix rows and columns in ascending order of moduli.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2013
//
//  Author:
//
//    Original FORTRAN77 version by Roger Payne.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Payne,
//    Inversion of matrices with contents subject to modulo arithmetic,
//    Applied Statistics,
//    Volume 46, Number 2, 1997, pages 295-298.
//
//  Parameters:
//
//    Input/output, int MAT[NROW*NROW].
//    On output, the matrix has been sorted.
//
//    Ignoreput, int IMAT[NROW*NROW].  
//    This quantity is ignored.
//
//    Input/output, int RMOD[NROW], the modulus for values in 
//    each row.  On output, these have been rearranged according to the sorting.
//
//    Input/output, int CMOD[NROW], the modulus for values in 
//    each column.  On output, these have been rearranged according to the 
//    sorting.
//
//    Output, int RSORT[NROW], the sorted row indices.
//
//    Output, int CSORT[NROW], the sorted column indices.
//
//    Input, int NROW, the order of the matrix.
//
{
  int i;
  int irc;
  int j;
  int jrc;
  int kirc;
  int kjrc;
  int p;
//
//  Initialize row and column addresses.
//
  for ( i = 1; i <= nrow; i++ )
  {
    rsort[i-1] = i;
    csort[i-1] = i;
  }
//
//  Sort the rows.
//
  for ( irc = 1; irc <= nrow; irc++ )
  {
//
//  Find the next row.
//
    jrc = irc;
    p = rmod[irc-1];

    for ( i = irc + 1; i <= nrow; i++ )
    {
      if ( rmod[i-1] < p )
      {
        p = rmod[i-1];
        jrc = i;
      }
    }

    if ( irc != jrc )
    {
      i = rmod[irc-1];
      rmod[irc-1] = rmod[jrc-1];
      rmod[jrc-1] = i;

      i = rsort[irc-1];
      rsort[irc-1] = rsort[jrc-1];
      rsort[jrc-1] = i;
//
//  Switch the rows.
//
      kirc = ( irc - 1 ) * nrow;
      kjrc = ( jrc - 1 ) * nrow;

      for ( j = 1; j <= nrow; j++ )
      {
        i = mat[kirc+j-1];
        mat[kirc+j-1] = mat[kjrc+j-1];
        mat[kjrc+j-1] = i;
      }
    }
  }
//
//  Sort the columns.
//
  for ( irc = 1; irc <= nrow; irc++ )
  {
//
//  Find the next column.
//
    jrc = irc;
    p = cmod[irc-1];

    for ( i = irc + 1; i <= nrow; i++ )
    {
      if ( cmod[i-1] < p )
      {
        p = cmod[i-1];
        jrc = i;
      }
    }

    if ( irc != jrc )
    {
      i = cmod[irc-1];
      cmod[irc-1] = cmod[jrc-1];
      cmod[jrc-1] = i;

      i = csort[irc-1];
      csort[irc-1] = csort[jrc-1];
      csort[jrc-1] = i;
//
//  Switch the columns.
//
      for ( j = 0; j < nrow * nrow; j = j + nrow )
      {
        i = mat[irc+j-1];
        mat[irc+j-1] = mat[jrc+j-1];
        mat[jrc+j-1] = i;
      }
    }
  }
  return;
}
//****************************************************************************80

void musort ( int mat[], int imat[], int rmod[], int cmod[], int rsort[], 
  int csort[], int nrow )

//****************************************************************************80
//
//  Purpose:
//
//    MUSORT unsorts the inverse matrix rows and columns into the original order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2013
//
//  Author:
//
//    Original FORTRAN77 version by Roger Payne.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Payne,
//    Inversion of matrices with contents subject to modulo arithmetic,
//    Applied Statistics,
//    Volume 46, Number 2, 1997, pages 295-298.
//
//  Parameters:
//
//    Input/output, int MAT[NROW*NROW].
//    On output, the matrix has been "unsorted".
//
//    Input/output, int IMAT[NROW*NROW].
//    On output, the matrix has been "unsorted".
//
//    Input/output, int RMOD[NROW], the modulus for values in 
//    each row.  On output, these have been restored to their original ordering.
//
//    Input/output, int CMOD[NROW], the modulus for values in
//    each column.  On output, these have been restored to their original 
//    ordering.
//
//    Input/output, int RSORT[NROW], the sorted row indices.
//
//    Input/output, int CSORT[NROW], the sorted column indices.
//
//    Input, int NROW, the order of the matrix.
//
{
  int i;
  int irc;
  int j;
  int jrc;
  int kirc;
  int kjrc;
//
//  Sort rows of inverse (= columns of original).
//
  for ( irc = 1; irc <= nrow; irc++ )
  {
//
//  Find next row.
//
    if ( csort[irc-1] != irc )
    {
      for ( jrc = irc + 1; jrc <= nrow; jrc++ )
      {
        if ( csort[jrc-1] == irc )
        {
          break;
        }
      }

      i = cmod[irc-1];
      cmod[irc-1] = cmod[jrc-1];
      cmod[jrc-1] = i;

      i = csort[irc-1];
      csort[irc-1] = csort[jrc-1];
      csort[jrc-1] = i;
//
//  Switch rows.
//
      kirc = ( irc - 1 ) * nrow;
      kjrc = ( jrc - 1 ) * nrow;

      for ( j = 1; j <= nrow; j++ )
      {
        i = imat[kirc+j-1];
        imat[kirc+j-1] = imat[kjrc+j-1];
        imat[kjrc+j-1] = i;
      }
    }
  }
//
//  Sort the columns of the inverse (= rows of original).
//
  for ( irc = 1; irc <= nrow; irc++ )
  {
//
//  Find the next column.
//
    if ( rsort[irc-1] != irc )
    {
      for ( jrc = irc + 1; jrc <= nrow; jrc++ )
      {
        if ( rsort[jrc-1] == irc )
        {
          break;
        }
      }

      i = rmod[irc-1];
      rmod[irc-1] = rmod[jrc-1];
      rmod[jrc-1] = i;

      i = rsort[irc-1];
      rsort[irc-1] = rsort[jrc-1];
      rsort[jrc-1] = i;
//
//  Switch the columns of IMAT.
//
      for ( j = 0; j < nrow * nrow; j = j + nrow )
      {
        i = imat[irc+j-1];
        imat[irc+j-1] = imat[jrc+j-1];
        imat[jrc+j-1] = i;
      }
//
//  Switch the diagonal elements of MAT (others are zero).
//
      kirc = ( irc - 1 ) * nrow + irc;
      kjrc = ( jrc - 1 ) * nrow + jrc;

      i = mat[kirc-1];
      mat[kirc-1] = mat[kjrc-1];
      mat[kjrc-1] = i;
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
