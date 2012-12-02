# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "cell.hpp"

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

int i4vec_max ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MAX returns the value of the maximum element in an I4VEC.
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
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MAX, the value of the maximum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }

  return value;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
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
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
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
    cout << "  " << setw(8) << i
         << ": " << setw(8) << a[i]  << "\n";
  }
  return;
}
//****************************************************************************80

double r8cvv_iget ( int mn, double a[], int m, int roff[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_IGET gets item J from row I in an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, int I, the row of the item.
//    0 <= I < M.
//
//    Input, int J, the column of the item.
//    0 <= J < NR[I].
//
//    Output, double R8CVV_IGET, the value of item A(I,J).
//
{
  double aij;
  int k;

  k = roff[i] + j;
  aij = a[k];

  return aij;
}
//****************************************************************************80

void r8cvv_iinc ( int mn, double a[], int m, int roff[], int i, int j, 
  double daij )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_IINC increments item J from row I in an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input/output, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, int I, the row of the item.
//    0 <= I < M.
//
//    Input, int J, the column of the item.
//    0 <= J < NR(I).
//
//    Input, double DAIJ, the increment to the value of item A(I,J).
//
{
  int k;

  k = roff[i] + j;
  a[k] = a[k] + daij;

  return;
}
//****************************************************************************80

void r8cvv_iset ( int mn, double a[], int m, int roff[], int i, int j, 
  double aij )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_ISET sets item J from row I in an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input/output, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, int I, the row of the item.
//    0 <= I < M.
//
//    Input, int J, the column of the item.
//    0 <= J < NR[I].
//
//    Input, double AIJ, the new value of item A(I,J).
//
{
  int k;

  k = roff[i] + j;
  a[k] = aij;

  return;
}
//****************************************************************************80

double *r8cvv_nget_new ( int mn, double a[], int m, int roff[], int nn, 
  int in[], int jn[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_NGET_NEW gets N items JN(*) from row IN(*) in an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, int NN, the number of items.
//
//    Input, int IN[NN], the rows of the items.
//    0 <= IN(*) < M.
//
//    Input, int JN[NN], the columns of the items.
//    0 <= JN(*) < NR(IN(*)).
//
//    Output, double R8CVV_NGET[NN], the value of items A(IN(*),JN(*)).
//
{
  int i;
  int k;
  double *vn;

  vn = new double[nn];

  for ( i = 0; i < nn; i++ )
  {
    k = roff[in[i]] + jn[i];
    vn[i] = a[k];
  }
  return vn;
}
//****************************************************************************80

void r8cvv_ninc ( int mn, double a[], int m, int roff[], int nn, int in[], 
  int jn[], double dvn[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_NINC increments items JN(*) from row IN(*) in an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input/output, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, int NN, the number of items.
//
//    Input, int IN[NN], the rows of the items.
//    0 <= IN(*) < M.
//
//    Input, int JN[NN], the columns of the items.
//    0 <= JN(*) < NR(IN(*)).
//
//    Input, double DVN[NN], the increments of items A(IN(*),JN(*)).
//
{
  int i;
  int k;

  for ( i = 0; i < nn; i++ )
  {
    k = roff[in[i]] + jn[i];
    a[k] = a[k] + dvn[i];
  }
  return;
}
//****************************************************************************80

void r8cvv_nset ( int mn, double a[], int m, int roff[], int nn, int in[], 
  int jn[], double vn[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_NSET sets items JN(*) from row IN(*) in an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input/output, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, int NN, the number of items.
//
//    Input, int IN[NN], the rows of the items.
//    0 <= IN(*) < M.
//
//    Input, int JN[NN], the columns of the items.
//    0 <= JN(*) < NR(IN(*)).
//
//    Input, double VN[NN], the new value of items A(IN(*),JN(*)).
//
{
  int i;
  int k;

  for ( i = 0; i < nn; i++ )
  {
    k = roff[in[i]] + jn[i];
    a[k] = vn[i];
  }
  return;
}
//****************************************************************************80

int *r8cvv_offset ( int m, int nr[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_OFFSET determines the row offsets of an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in the array.
//
//    Input, int NR[M], the row sizes.
//
//    Output, int R8CVV_OFFSET[M+1], the row offsets.
//
{
  int i;
  int *roff;

  roff = new int[m+1];

  roff[0] = 0;
  for ( i = 0; i < m; i++ )
  {
    roff[i+1] = roff[i] + nr[i];
  }

  return roff;
}
//****************************************************************************80

void r8cvv_print ( int mn, double a[], int m, int roff[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_PRINT prints an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int j;
  int k;
  int k1;
  int k2;
  int khi;
  int klo;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  for ( i = 0; i < m; i++ )
  {
    k1 = roff[i];
    k2 = roff[i+1];

    for ( klo = k1; klo < k2; klo = klo + 5 )
    {
      khi = i4_min ( klo + 5, k2 );
      if ( klo == k1 )
      {
        cout << setw(5) << i;
      }
      else
      {
        cout << "     ";
      }
      cout << "  ";
      for ( k = klo; k < khi; k++ )
      {
        cout << setw(14) << a[k];
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

double *r8cvv_rget_new ( int mn, double a[], int m, int roff[], int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_RGET_NEW gets row I from an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, int I, the row.
//    0 <= I < M.
//
//    Output, double R8CVV_GET[NR[I]], the value of A(I,*).
//
{
  double *ai;
  int j;
  int k1;
  int k2;
  int nv;

  k1 = roff[i];
  k2 = roff[i+1];
  nv = k2 - k1;
  ai = new double[nv];
  for ( j = 0; j < nv; j++ )
  {
    ai[j] = a[k1+j];
  }

  return ai;
}
//****************************************************************************80

void r8cvv_rinc ( int mn, double a[], int m, int roff[], int i, double dai[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_RINC increments row I in an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input/output, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, int I, the row.
//    0 <= I < M.
//
//    Input, double DAI[NR[I]], the increment for A(I,*).
//
{
  int j;
  int k1;
  int k2;
  int nv;

  k1 = roff[i];
  k2 = roff[i+1];
  nv = k2 - k1;
  for ( j = 0; j < nv; j++ )
  {
    a[k1+j] = a[k1+j] + dai[j];
  }

  return;
}
//****************************************************************************80

void r8cvv_rset ( int mn, double a[], int m, int roff[], int i, double ai[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_RSET sets row I from an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MN, the size of the cell array.
//
//    Input/output, double A[MN], the cell array.
//
//    Input, int M, the number of rows in the array.
//
//    Input, int ROFF[M+1], the row offsets.
//
//    Input, int I, the row.
//    0 <= I < M.
//
//    Input, double AI[NR[I]], the new value of A(I,*).
//
{
  int j;
  int k1;
  int k2;
  int nv;

  k1 = roff[i];
  k2 = roff[i+1];
  nv = k2 - k1;
  for ( j = 0; j < nv; j++ )
  {
    a[k1+j] = ai[j];
  }

  return;
}
//****************************************************************************80

int r8cvv_size ( int m, int nr[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8CVV_SIZE determines the size of an R8CVV.
//
//  Discussion:
//
//    An R8CVV is a "vector of vectors" of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in the array.
//
//    Input, int NR[M], the size of each row.
//
//    Output, int R8CVV_SIZE, the size of the cell array.
//
{
  int i;
  int mn;

  mn = 0;
  for ( i = 0; i < m; i++ )
  {
    mn = mn + nr[i];
  }
  return mn;
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

void r8vec_transpose_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Example:
//
//    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
//    TITLE = 'My vector:  '
//
//    My vector:
//
//        1.0    2.1    3.2    4.3    5.4
//        6.5    7.6    8.7    9.8   10.9
//       11.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2010
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
  int ihi;
  int ilo;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= 0 )
  {
    cout << "  (Empty)\n";
    return;
  }

  for ( ilo = 0; ilo < n; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 5, n );
    for ( i = ilo; i < ihi; i++ )
    {
      cout << "  " << setw(12) << a[i];
    }
    cout << "\n";
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
