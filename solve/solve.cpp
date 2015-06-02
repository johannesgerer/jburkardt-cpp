# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "solve.hpp"

//****************************************************************************80

double **r8rmat_copy_new ( int m, int n, double **a )

//****************************************************************************80
//
//  Purpose:
//
//    R8RMAT_COPY_NEW makes a new copy of an R8RMAT .
//
//  Discussion:
//
//    An R8RMAT is a matrix stored in row major form, using M pointers
//    to the beginnings of rows.
//
//    A declaration of the form
//      double **a;
//    is necesary.  Then an assignment of the form:
//      a = r8rmat_new ( m, n );
//    allows the user to assign entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17.0;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double **A, the array to copy.
//
//    Output, double **R8RMAT_COPY_NEW, the copied array.
//
{
  double **b;
  int i;
  int j;

  b = r8rmat_new ( m, n );

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i][j] = a[i][j];
    }
  }
  return b;
}
//****************************************************************************80

void r8rmat_delete ( int m, int n, double **a )

//****************************************************************************80
//
//  Purpose:
//
//    R8RMAT_DELETE frees memory associated with an R8RMAT.
//
//  Discussion:
//
//    This function releases the memory associated with an R8RMAT.
// 
//    An R8RMAT is a row-major array that was created by a 
//    command like:
//
//      double **a;
//      a = r8rmat_new ( m, n );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input, double **A, the pointer to the array.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    delete [] a[i];
  }

  delete [] a;

  return;
}
//****************************************************************************80

double *r8rmat_fs_new ( int n, double **a, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8RMAT_FS_NEW factors and solves an R8RMAT system with one right hand side.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double **A, the coefficient matrix of the linear system.
//
//    Input, double B[N], the right hand side of the linear system.
//
//    Output, double R8RMAT_FS_NEW[N], the solution of the linear system.
//
{
  double **a2;
  int i;
  int j;
  int k;
  int p;
  double t;
  double *x;
//
//  Create a copy of the matrix.
//
  a2 = new double *[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = new double[n];
  }

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a2[i][j] = a[i][j];
    }
  }
//
//  Create X and set it to B.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( k = 0; k < n; k++ )
  {
//
//  Find the maximum element in column I.
//
    p = k;

    for ( i = k + 1; i < n; i++ )
    {
      if ( fabs ( a2[p][k] ) < fabs ( a2[i][k] ) )
      {
        p = i;
      }
    }

    if ( a2[p][k] == 0.0 )
    {
      cerr << "\n";
      cerr << "R8RMAT_FS_NEW - Fatal error!\n";
      cerr << "  Zero pivot on step " << k << "\n";
      exit ( 1 );
    }
//
//  Switch rows K and P.
//
    if ( k != p )
    {
      for ( j = 0; j < n; j++ )
      {
        t        = a2[k][j];
        a2[k][j] = a2[p][j];
        a2[p][j] = t;
      }
      t    = x[k];
      x[k] = x[p];
      x[p] = t;
    }
//
//  Scale the pivot row.
//
    t = a2[k][k];
    a2[k][k] = 1.0;
    for ( j = k + 1; j < n; j++ )
    {
      a2[k][j] = a2[k][j] / t;
    }
    x[k] = x[k] / t;
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = k + 1; i < n; i++ )
    {
      if ( a2[i][k] != 0.0 )
      {
        t = - a2[i][k];
        a2[i][k] = 0.0;
        for ( j = k + 1; j < n; j++ )
        {
          a2[i][j] = a2[i][j] + t * a2[k][j];
        }
        x[i] = x[i] + t * x[k];
      }
    }
  }
//
//  Back solve.
//
  for ( j = n - 1; 1 <= j; j-- )
  {
    for ( i = 0; i < j; i++ )
    {
      x[i] = x[i] - a2[i][j] * x[j];
    }
  }
//
//  Delete the matrix copy.
//
  for ( i = 0; i < n; i++ )
  {
    delete [] a2[i];
  }

  delete [] a2;

  return x;
}
//****************************************************************************80

double **r8rmat_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8RMAT_NEW allocates a new R8RMAT.
//
//  Discussion:
//
//    An R8RMAT is a row-major array that was created by a 
//    command like:
//
//      double **a;
//      a = r8rmat_new ( m, n );
//
//    The user assigns entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17.0;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, double **R8RMAT_NEW, a new matrix.
//
{
  double **a;
  int i;

  a = new double *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8RMAT_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new double[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R8RMAT_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
}
//****************************************************************************80

double **r8rmat_zero ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8RMAT_ZERO allocates and zeroes a new R8RMAT.
//
//  Discussion:
//
//    An R8RMAT is a row-major array that was created by a 
//    command like:
//
//      double **a;
//      a = r8rmat_new ( m, n );
//
//    The user assigns entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17.0;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, double **R8RMAT_ZERO, a new matrix.
//
{
  double **a;
  int i;
  int j;

  a = new double *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "R8RMAT_ZERO - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new double[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "R8RMAT_ZERO - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i][j] = 0.0;
    }
  }
  return a;
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
