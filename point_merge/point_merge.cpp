# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "point_merge.hpp"

//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
//
//    CPU_TIME reports the elapsed CPU time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
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

int i4_uniform ( int a, int b, int *seed )

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
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
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

int point_radial_tol_unique_count ( int m, int n, double a[], double tol,
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_RADIAL_TOL_UNIQUE_COUNT counts the tolerably unique points.
//
//  Discussion:
//
//    The input data is an M x N array A, representing the M-dimensional
//    coordinates of N points.
//
//    The output is the number of tolerably unique points in the list.
//
//    This program performs the same task as POINT_TOL_UNIQUE_COUNT.
//    But that program is guaranteed to use N^2 comparisons.
//
//    It is hoped that this function, on the other hand, will tend
//    to use O(N) comparisons after an O(NLog(N)) sort.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, double A[M*N], the array of N columns of data.
//
//    Input, double TOL, a tolerance for equality.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, int POINT_RADIAL_TOL_UNIQUE_COUNT, the number of tolerably
//    unique points.
//
{
  double dist;
  int hi;
  int i;
  int *indx;
  int j;
  int k;
  double *r;
  bool *unique;
  int unique_num;
  double *w;
  double w_sum;
  double *z;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }
//
//  Assign a base point Z randomly in the convex hull.
//
  w = r8vec_uniform_01_new ( n, seed );
  w_sum = r8vec_sum ( n, w );
  for ( j = 0; j < n; j++ )
  {
    w[j] = w[j] / w_sum;
  }

  z = new double[m];
  for ( i = 0; i < m; i++ )
  {
    z[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      z[i] = z[i] + a[i+j*m] * w[j];
    }
  }
//
//  Compute the radial distance R of each point to Z.
//
  r = new double[n];

  for ( j = 0; j < n; j++ )
  {
    r[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r[j] = r[j] + pow ( a[i+j*m] - z[i], 2 );
    }
    r[j] = sqrt ( r[j] );
  }
//
//  Implicitly sort the R array.
//
  indx = r8vec_sort_heap_index_a ( n, r );
//
//  To determine if a point I is tolerably unique, we only have to check
//  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
//
  unique_num = 0;

  unique = new bool[n];
  for ( i = 0; i < n; i++ )
  {
    unique[i] = true;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( unique[indx[i]] )
    {
//
//  Point INDX(I) is unique, in that no earlier point is near it.
//
      unique_num = unique_num + 1;
//
//  Look for later points which are close to point INDX(I)
//  in terms of R.
//
      hi = i;

      while ( hi < n - 1 )
      {
        if ( r[indx[i]] + tol < r[indx[hi+1]] )
        {
          break;
        }
        hi = hi + 1;
      }
//
//  Points INDX(I+1) through INDX(HI) have an R value close to
//  point INDX(I).  Are they truly close to point INDEX(I)?
//
      for ( j = i + 1; j <= hi; j++ )
      {
        if ( unique[indx[j]] )
        {
          dist = 0.0;
          for ( k = 0; k < m; k++ )
          {
            dist = dist + pow ( a[k+indx[i]*m] - a[k+indx[j]*m], 2 );
          }
          dist = sqrt ( dist );

          if ( dist <= tol )
          {
            unique[indx[j]] = false;
          }
        }
      }
    }
  }

  delete [] indx;
  delete [] r;
  delete [] unique;
  delete [] w;
  delete [] z;

  return unique_num;
}
//****************************************************************************80

int point_radial_tol_unique_index ( int m, int n, double a[], double tol,
  int *seed, int undx[], int xdnu[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_RADIAL_TOL_UNIQUE_INDEX indexes the tolerably unique points.
//
//  Discussion:
//
//    The input data is an M x N array A, representing the M-dimensional
//    coordinates of N points.
//
//    The output is:
//    * the number of tolerably unique points in the list;
//    * the index, in the list of unique items, of the representatives
//      of each point;
//    * the index, in A, of the tolerably unique representatives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, double A[M*N], the array of N columns of data.
//
//    Input, double TOL, a tolerance for equality.
//
//    Input/output, int SEED, a seed for the random
//    number generator.
//
//    Output, int UNDX[UNIQUE_NUM], the index, in A, of the
//    tolerably unique points.
//
//    Output, int XDNU[N], the index, in UNDX, of the
//    tolerably unique point that "represents" this point.
//
//    Output, int POINT_RADIAL_TOL_UNIQUE_INDEX, the number of tolerably
//    unique points.
//
{
  double dist;
  int hi;
  int i;
  int *indx;
  int j;
  int k;
  double *r;
  bool *unique;
  int unique_num;
  double *w;
  double w_sum;
  double *z;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }
//
//  Assign a base point Z randomly in the convex hull.
//
  w = r8vec_uniform_01_new ( n, seed );
  w_sum = r8vec_sum ( n, w );
  for ( j = 0; j < n; j++ )
  {
    w[j] = w[j] / w_sum;
  }

  z = new double[m];
  for ( i = 0; i < m; i++ )
  {
    z[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      z[i] = z[i] + a[i+j*m] * w[j];
    }
  }
//
//  Compute the radial distance R of each point to Z.
//
  r = new double[n];

  for ( j = 0; j < n; j++ )
  {
    r[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r[j] = r[j] + pow ( a[i+j*m] - z[i], 2 );
    }
    r[j] = sqrt ( r[j] );
  }
//
//  Implicitly sort the R array.
//
  indx = r8vec_sort_heap_index_a ( n, r );
//
//  To determine if a point I is tolerably unique, we only have to check
//  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
//
  unique_num = 0;

  unique = new bool[n];
  for ( i = 0; i < n; i++ )
  {
    unique[i] = true;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( unique[indx[i]] )
    {
//
//  Point INDX(I) is unique, in that no earlier point is near it.
//
      xdnu[indx[i]] = unique_num;
      undx[unique_num] = indx[i];
      unique_num = unique_num + 1;
//
//  Look for later points which are close to point INDX(I)
//  in terms of R.
//
      hi = i;

      while ( hi < n - 1 )
      {
        if ( r[indx[i]] + tol < r[indx[hi+1]] )
        {
          break;
        }
        hi = hi + 1;
      }
//
//  Points INDX(I+1) through INDX(HI) have an R value close to
//  point INDX(I).  Are they truly close to point INDEX(I)?
//
      for ( j = i + 1; j <= hi; j++ )
      {
        if ( unique[indx[j]] )
        {
          dist = 0.0;
          for ( k = 0; k < m; k++ )
          {
            dist = dist + pow ( a[k+indx[i]*m] - a[k+indx[j]*m], 2 );
          }
          dist = sqrt ( dist );

          if ( dist <= tol )
          {
            unique[indx[j]] = false;
            xdnu[indx[j]] = xdnu[indx[i]];
          }
        }
      }
    }
  }

  delete [] indx;
  delete [] r;
  delete [] unique;
  delete [] w;
  delete [] z;

  return unique_num;
}
//****************************************************************************80

int point_radial_unique_count ( int m, int n, double a[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_RADIAL_UNIQUE_COUNT counts the unique points.
//
//  Discussion:
//
//    The input data is an M x N array A, representing the M-dimensional
//    coordinates of N points.
//
//    The output is the number of unique points in the list.
//
//    This program performs the same task as POINT_UNIQUE_COUNT, and
//    carries out more work.  Hence, it is not a substitute for
//    POINT_UNIQUE_COUNT.  Instead, it is intended to be a starting point
//    for a similar program which includes a tolerance.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, double A[M*N], the array of N columns of data.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, int POINT_RADIAL_UNIQUE_COUNT, the number of unique points.
//
{
  bool equal;
  int hi;
  int i;
  int *indx;
  int j;
  int j1;
  int j2;
  int lo;
  double *r;
  bool unique;
  int unique_index;
  int unique_num;
  double *w;
  double w_sum;
  double *z;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }
//
//  Assign a base point Z randomly in the convex hull.
//
  w = r8vec_uniform_01_new ( n, seed );
  w_sum = r8vec_sum ( n, w );
  for ( j = 0; j < n; j++ )
  {
    w[j] = w[j] / w_sum;
  }

  z = new double[m];
  for ( i = 0; i < m; i++ )
  {
    z[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      z[i] = z[i] + a[i+j*m] * w[j];
    }
  }
//
//  Compute the radial distance R of each point to Z.
//
  r = new double[n];

  for ( j = 0; j < n; j++ )
  {
    r[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r[j] = r[j] + pow ( a[i+j*m] - z[i], 2 );
    }
    r[j] = sqrt ( r[j] );
  }
//
//  Implicitly sort the R array.
//
  indx = r8vec_sort_heap_index_a ( n, r );
//
//  To determine if a point is unique, we only have to check
//  whether it is distinct from all points with the same
//  R value and lower ordering.
//
  unique_num = 0;
  hi = - 1;

  while ( hi < n - 1 )
  {
//
//  Advance LO.
//
    lo = hi + 1;
//
//  Extend HI.
//
    hi = lo;

    while ( hi < n - 1 )
    {
      if ( r[indx[hi+1]] == r[indx[lo]] )
      {
        hi = hi + 1;
      }
      else
      {
        break;
      }
    }
//
//  Points INDX(LO) through INDX(HI) have same R value.
//
//  Find the unique ones.
//
    unique_num = unique_num + 1;

    for ( j1 = lo + 1; j1 <= hi; j1++ )
    {
      for ( j2 = lo; j2 < j1; j2++ )
      {
        equal = true;
        for ( i = 0; i < m; i++ )
        {
          if ( a[i+indx[j2]*m] != a[i+indx[j1]*m] )
          {
            equal = false;
            break;
          }
        }
        if ( equal )
        {
          break;
        }
      }

      if ( !equal )
      {
        unique_num = unique_num + 1;
      }
    }
  }

  delete [] indx;
  delete [] r;
  delete [] w;
  delete [] z;

  return unique_num;
}
//****************************************************************************80

int point_tol_unique_count ( int m, int n, double a[], double tol )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_TOL_UNIQUE_COUNT counts the tolerably unique points.
//
//  Discussion:
//
//    The input data is an M x N array A, representing the M-dimensional
//    coordinates of N points.
//
//    This function uses a simple but expensive approach.  The first point
//    is accepted as unique.  Each subsequent point is accepted as unique
//    only if it is at least a tolerance away from all accepted unique points.
//    This means the expected amount of work is O(N^2).
//
//    The output is the number of unique points in the list.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, double A[M*N], the array of N columns of data.
//
//    Input, double TOL, a tolerance.
//
//    Output, int POINT_TOL_UNIQUE_COUNT, the number of unique points.
//
{
  double dist;
  int i;
  int j;
  int k;
  bool *unique;
  int unique_num;

  unique = new bool[n];

  for ( i = 0; i < n; i++ )
  {
    unique[i] = true;
  }
  unique_num = n;

  for ( i = 1; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      if ( unique[j] )
      {
        dist = 0.0;
        for ( k = 0; k < m; k++ )
        {
          dist = dist + pow ( a[k+i*m] - a[k+j*m], 2 );
        }
        dist = sqrt ( dist );
        if ( dist <= tol )
        {
          unique[i] = false;
          unique_num = unique_num - 1;
          break;
        }
      }
    }
  }
  delete [] unique;

  return unique_num;
}
//****************************************************************************80

int point_tol_unique_index ( int m, int n, double a[], double tol, int xdnu[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_TOL_UNIQUE_INDEX indexes the tolerably unique points.
//
//  Discussion:
//
//    This routine uses an algorithm that is O(N^2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the data values.
//
//    Input, int N, the number of data values.
//
//    Input, double A[M*N], the data values.
//
//    Input, double TOL, a tolerance for equality.
//
//    Output, int XDNU[N], the index, in A, of the tolerably unique
//    point that "represents" this point.
//
//    Output, int POINT_TOL_UNIQUE_INDEX, the number of tolerably
//    unique points.
//
{
  double dist;
  int i;
  int j;
  int k;
  bool *unique;
  int unique_num;

  unique = new bool[n];

  for ( i = 0; i < n; i++ )
  {
    unique[i] = true;
  }
  for ( i = 0; i < n; i++ )
  {
    xdnu[i] = i;
  }
  unique_num = n;

  i = 0;
  xdnu[0] = 0;

  for ( i = 1; i < n; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      if ( unique[j] )
      {
        dist = 0.0;
        for ( k = 0; k < m; k++ )
        {
          dist = dist + pow ( a[k+i*m] - a[k+j*m], 2 );
        }
        dist = sqrt ( dist );
        if ( dist <= tol )
        {
          unique[i] = false;
          unique_num = unique_num - 1;
          xdnu[i] = j;
          break;
        }
      }
    }
  }

  delete [] unique;

  return unique_num;
}
//****************************************************************************80

int point_unique_count ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_UNIQUE_COUNT counts the unique points.
//
//  Discussion:
//
//    The input data is an M x N array A, representing the M-dimensional
//    coordinates of N points.
//
//    The algorithm relies on the fact that, in a sorted list, points that
//    are exactly equal must occur consecutively.
//
//    The output is the number of unique points in the list.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, double A[M*N], the array of N columns of data.
//
//    Output, int POINT_UNIQUE_COUNT, the number of unique points.
//
{
  int i;
  int *indx;
  int j;
  int unique_index;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }
//
//  Implicitly sort the array.
//
  indx = r8col_sort_heap_index_a ( m, n, a );
//
//  Two points are considered equal only if they exactly match.
//  In that case, equal points can only occur as consecutive items
//  in the sorted list.   This makes counting easy.
//
  unique_num = 1;
  unique_index = indx[0];

  for ( j = 1; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+unique_index*m] != a[i+indx[j]*m] )
      {
        unique_num = unique_num + 1;
        unique_index = indx[j];
      }
    }
  }

  delete [] indx;

  return unique_num;
}
//****************************************************************************80

void point_unique_index ( int m, int n, double a[], int unique_num, int undx[],
  int xdnu[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINT_UNIQUE_INDEX indexes unique points.
//
//  Discussion:
//
//    An R8COL is an M by N array of R8's, regarded as an array of N columns,
//    each of length M.
//
//    The goal of this routine is to determine a vector UNDX,
//    which points to the unique elements of A, in sorted order,
//    and a vector XDNU, which identifies, for each entry of A, the index of
//    the unique sorted element of A.
//
//    This is all done with index vectors, so that the elements of
//    A are never moved.
//
//    The first step of the algorithm requires the indexed sorting
//    of A, which creates arrays INDX and XDNI.  (If all the entries
//    of A are unique, then these arrays are the same as UNDX and XDNU.)
//
//    We then use INDX to examine the entries of A in sorted order,
//    noting the unique entries, creating the entries of XDNU and
//    UNDX as we go.
//
//    Once this process has been completed, the vector A could be
//    replaced by a compressed vector XU, containing the unique entries
//    of A in sorted order, using the formula
//
//      XU(*) = A(UNDX(*)).
//
//    We could then, if we wished, reconstruct the entire vector A, or
//    any element of it, by index, as follows:
//
//      A(I) = XU(XDNU(I)).
//
//    We could then replace A by the combination of XU and XDNU.
//
//    Later, when we need the I-th entry of A, we can locate it as
//    the XDNU(I)-th entry of XU.
//
//    Here is an example of a vector A, the sort and inverse sort
//    index vectors, and the unique sort and inverse unique sort vectors
//    and the compressed unique sorted vector.
//
//      I     A  Indx  Xdni       XU  Undx  Xdnu
//    ----+-----+-----+-----+--------+-----+-----+
//      0 | 11.     0     0 |    11.     0     0
//      1 | 22.     2     4 |    22.     1     1
//      2 | 11.     5     1 |    33.     3     0
//      3 | 33.     8     7 |    55.     4     2
//      4 | 55.     1     8 |                  3
//      5 | 11.     6     2 |                  0
//      6 | 22.     7     5 |                  1
//      7 | 22.     3     6 |                  1
//      8 | 11.     4     3 |                  0
//
//    INDX(2) = 3 means that sorted item(2) is A(3).
//    XDNI(2) = 5 means that A(2) is sorted item(5).
//
//    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
//    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
//
//    XU(XDNU(I))) = A(I).
//    XU(I)        = A(UNDX(I)).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the data values.
//
//    Input, int N, the number of data values,
//
//    Input, double A[M*N], the data values.
//
//    Input, int UNIQUE_NUM, the number of unique values in A.
//    This value is only required for languages in which the size of
//    UNDX must be known in advance.
//
//    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
//
//    Output, int XDNU[N], the XDNU vector.
//
{
  double diff;
  int i;
  int *indx;
  int j;
  int k;
//
//  Implicitly sort the array.
//
  indx = r8col_sort_heap_index_a ( m, n, a );
//
//  Walk through the implicitly sorted array.
//
  i = 0;

  j = 0;
  undx[j] = indx[i];

  xdnu[indx[i]] = j;

  for ( i = 1; i < n; i++ )
  {
    diff = 0.0;
    for ( k = 0; k < m; k++ )
    {
      diff = r8_max ( diff, r8_abs ( a[k+indx[i]*m] - a[k+undx[j]*m] ) );
    }
    if ( 0.0 < diff )
    {
      j = j + 1;
      undx[j] = indx[i];
    }
    xdnu[indx[i]] = j;
  }
  delete [] indx;

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

double *r8col_duplicates ( int m, int n, int n_unique, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_DUPLICATES generates an R8COL with some duplicate columns.
//
//  Discussion:
//
//    An R8COL is an M by N array of R8's, regarded as an array of N columns,
//    each of length M.
//
//    This routine generates a random R8COL with a specified number of
//    duplicate columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in each column of A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int N_UNIQUE, the number of unique columns in A.
//    1 <= N_UNIQUE <= N.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, double R8COL_DUPLICATES[M*N], the array.
//
{
  double *a;
  int i;
  int j1;
  int j2;
  double temp;

  if ( n_unique < 1 || n < n_unique )
  {
    cout << "\n";
    cout << "R8COL_DUPLICATES - Fatal error!\n";
    cout << "  1 <= N_UNIQUE <= N is required.\n";
    exit ( 1 );
  }

  a = new double[m*n];

  r8mat_uniform_01 ( m, n_unique, seed, a );
//
//  Randomly copy unique columns.
//
  for ( j1 = n_unique; j1 < n; j1++ )
  {
    j2 = i4_uniform ( 0, n_unique - 1, seed );
    for ( i = 0; i < m; i++ )
    {
      a[i+j1*m] = a[i+j2*m];
    }
  }
//
//  Permute the columns.
//
  for ( j1 = 0; j1 < n; j1 ++ )
  {
    j2 = i4_uniform ( j1, n - 1, seed );
    for ( i = 0; i < m; i++ )
    {
      temp      = a[i+j1*m];
      a[i+j1*m] = a[i+j2*m];
      a[i+j2*m] = temp;
    }
  }
  return a;
}
//****************************************************************************80

int *r8col_sort_heap_index_a ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8COL.
//
//  Discussion:
//
//    An R8COL is an M by N array of R8's, regarded as an array of N columns,
//    each of length M.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2) is negative.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      A(*,INDX(*)) is sorted,
//
//    Note that the index vector is 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in each column of A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the array.
//
//    Output, int R8COL_SORT_HEAP_INDEX_A[N], contains the sort index.  The
//    I-th column of the sorted array is A(*,INDX(I)).
//
{
  double *column;
  int i;
  int *indx;
  int indxt;
  int ir;
  int isgn;
  int j;
  int k;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0];
    return indx;
  }

  column = new double[m];

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      for ( k = 0; k < m; k++ )
      {
        column[k] = a[k+indxt*m];
      }
    }
    else
    {
      indxt = indx[ir-1];
      for ( k = 0; k < m; k++ )
      {
        column[k] = a[k+indxt*m];
      }
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        isgn = r8vec_compare ( m, a+indx[j-1]*m, a+indx[j]*m );

        if ( isgn < 0 )
        {
          j = j + 1;
        }
      }

      isgn = r8vec_compare ( m, column, a+indx[j-1]*m );

      if ( isgn < 0 )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
  delete [] column;

  for ( i = 0; i < n; i++ )
  {
    indx[i] = indx[i];
  }

  return indx;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << ":";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_uniform_01 ( int m, int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
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
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
//****************************************************************************80

int r8vec_compare ( int n, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COMPARE compares two R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The lexicographic ordering is used.
//
//  Example:
//
//    Input:
//
//      A1 = ( 2.0, 6.0, 2.0 )
//      A2 = ( 2.0, 8.0, 12.0 )
//
//    Output:
//
//      ISGN = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A[N], B[N], the vectors to be compared.
//
//    Output, int R8VEC_COMPARE, the results of the comparison:
//    -1, A is lexicographically less than B,
//     0, A is equal to B,
//    +1, A is lexicographically greater than B.
//
{
  int isgn;
  int k;

  isgn = 0;

  for ( k = 0; k < n; k++ )
  {
    if ( a[k] < b[k] )
    {
      isgn = -1;
      return isgn;
    }
    else if ( b[k] < a[k] )
    {
      isgn = +1;
      return isgn;
    }
  }
  return isgn;
}
//****************************************************************************80

double r8vec_norm_l2 ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], the vector whose L2 norm is desired.
//
//    Output, double R8VEC_NORM_L2, the L2 norm of A.
//
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
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

int *r8vec_sort_heap_index_a ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(indx(*))
//
//    or explicitly, by the call
//
//      r8vec_permute ( n, indx, 0, a )
//
//    after which a(*) is sorted.
//
//    Note that the index vector is 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], an array to be index-sorted.
//
//    Output, int R8VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0];
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}
//****************************************************************************80

void r8vec_uniform_01 ( int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
//    Output, double R[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

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

  return;
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
