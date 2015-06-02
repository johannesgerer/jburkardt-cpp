# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "pdflib.hpp"
# include "problem1c_covariance.hpp"

//****************************************************************************80

Covariance::Covariance ( int par_num )

//****************************************************************************80
{
  int i;
  int j;

  order = par_num;

  Covariance::array_set ( );

  Covariance::factor_set ( );

  Covariance::det_set ( );

  Covariance::inv_set ( );

  Covariance::mean_set ( );

  return;
}
//****************************************************************************80

double *Covariance::array_get ( )

//****************************************************************************80
{
  double *array2;

  array2 = r8mat_copy_new ( order, order, array );

  return array2;
}
//****************************************************************************80

void Covariance::array_set ( )

//****************************************************************************80
{
  int i;
  int j;

  array = new double[order*order];

  for ( j = 0; j < order; j++ )
  {
    for ( i = 0; i < order; i++ )
    {
      array[i+j*order] = 0.5;
    }
  }

  for ( i = 0; i < order; i++ )
  {
    array[i+i*order] = ( double ) ( i + 1 );
  }
  return;
}
//****************************************************************************80

double Covariance::det_get ( )

//****************************************************************************80
{
  double det2;

  det2 = det;

  return det2;
}
//****************************************************************************80

void Covariance::det_set ( )

//****************************************************************************80
{
  det = r8mat_podet ( order, factor );

  return;
}
//****************************************************************************80

double *Covariance::factor_get ( )

//****************************************************************************80
{
  double *factor2;

  factor2 = r8mat_copy_new ( order, order, factor );

  return factor2;
}
//****************************************************************************80

void Covariance::factor_set ( )

//****************************************************************************80
{
  factor = r8mat_pofac ( order, array );

  return;
}
//****************************************************************************80

double *Covariance::inv_get ( )

//****************************************************************************80
{
  double *inv2;

  inv2 = r8mat_copy_new ( order, order, inv );

  return inv2;
}
//****************************************************************************80

void Covariance::inv_set ( )

//****************************************************************************80
{
  inv = r8mat_poinv ( order, factor );

  return;
}
//****************************************************************************80

double *Covariance::mean_get ( )

//****************************************************************************80
{
  double *mean2;

  mean2 = r8vec_copy_new ( order, mean );

  return mean2;
}
//****************************************************************************80

void Covariance::mean_set ( )

//****************************************************************************80
{
  int i;

  mean = new double[order];

  for ( i = 0; i < order; i++ )
  {
    mean[i] = 0.0;
  }
  return;
}
//****************************************************************************80

Covariance::~Covariance ( )

//****************************************************************************80
{
  return;
}
//****************************************************************************80

void Covariance::print ( string title )

//****************************************************************************80
{
  r8mat_print ( order, order, array, title );

  return;
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

double *r8vec_copy_new ( int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY_NEW copies an R8VEC to a new R8VEC.
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

