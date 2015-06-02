# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "simplex_grid.hpp"

//****************************************************************************80

void comp_next_grlex ( int kc, int xc[] )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT_GRLEX returns the next composition in grlex order.
//
//  Discussion:
//
//    Example:
//
//    KC = 3
//
//    #   XC(1) XC(2) XC(3)  Degree
//      +------------------------
//    1 |  0     0     0        0
//      |
//    2 |  0     0     1        1
//    3 |  0     1     0        1
//    4 |  1     0     0        1
//      |
//    5 |  0     0     2        2
//    6 |  0     1     1        2
//    7 |  0     2     0        2
//    8 |  1     0     1        2
//    9 |  1     1     0        2
//   10 |  2     0     0        2
//      |
//   11 |  0     0     3        3
//   12 |  0     1     2        3
//   13 |  0     2     1        3
//   14 |  0     3     0        3
//   15 |  1     0     2        3
//   16 |  1     1     1        3
//   17 |  1     2     0        3
//   18 |  2     0     1        3
//   19 |  2     1     0        3
//   20 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int KC, the number of parts of the composition.
//    1 <= KC.
//
//    Input/output, int XC[KC], the current composition.
//    Each entry of XC must be nonnegative.
//    On return, XC has been replaced by the next composition in the
//    grlex order.
//
{
  int i;
  int im1;
  int j;
  int t;
//
//  Ensure that 1 <= KC.
//
  if ( kc < 1 )
  {
    cerr << "\n";
    cerr << "COMP_NEXT_GRLEX - Fatal error!\n";
    cerr << "  KC < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 0 <= XC(I).
//
  for ( i = 0; i < kc; i++ )
  {
    if ( xc[i] < 0 )
    {
      cerr << "\n";
      cerr << "COMP_NEXT_GRLEX - Fatal error!\n";
      cerr << "  XC[I] < 0\n";
      exit ( 1 );
    }
  }
//
//  Find I, the index of the rightmost nonzero entry of X.
//
  i = 0;
  for ( j = kc; 1 <= j; j-- )
  {
    if ( 0 < xc[j-1] )
    {
      i = j;
      break;
    }
  }
//
//  set T = X(I)
//  set XC(I) to zero,
//  increase XC(I-1) by 1,
//  increment XC(KC) by T-1.
//
  if ( i == 0 )
  {
    xc[kc-1] = 1;
    return;
  }
  else if ( i == 1 )
  {
    t = xc[0] + 1;
    im1 = kc;
  }
  else if ( 1 < i )
  {
    t = xc[i-1];
    im1 = i - 1;
  }

  xc[i-1] = 0;
  xc[im1-1] = xc[im1-1] + 1;
  xc[kc-1] = xc[kc-1] + t - 1;

  return;
}
//****************************************************************************80

int *comp_random_new ( int n, int k, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_RANDOM_NEW selects a random composition of the integer N into K parts.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer to be decomposed.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int COMP_RANDOM_NEW[K], the parts of the composition.
//
{
  int *a;
  int i;
  int l;
  int m;

  a = ksub_random_new ( n+k-1, k-1, seed );

  a[k-1] = n + k;
  l = 0;

  for ( i = 0; i < k; i++ )
  {
    m = a[i];
    a[i] = a[i] - l - 1;
    l = m;
  }

  return a;
}
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
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
//    02 October 2012
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
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

void i4mat_transpose_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
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
//    31 January 2005
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
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
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
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i - 1 << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = jlo;
    if ( j2lo < 1 )
    {
      j2lo = 1;
    }
    j2hi = jhi;
    if ( n < jhi )
    {
      j2hi = n;
    }

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j - 1 << ":";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

int *ksub_random_new ( int n, int k, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_RANDOM_NEW selects a random subset of size K from a set of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, number of elements in desired subsets.  K must
//    be between 0 and N.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int KSUB_RANDOM_NEW[K].  A(I) is the I-th element of the
//    output set.  The elements of A are in order.
//
{
  int *a;
  int i;
  int ids;
  int ihi;
  int ip;
  int ir;
  int is;
  int ix;
  int l;
  int ll;
  int m;
  int m0;

  if ( k < 0 )
  {
    cerr << "\n";
    cerr << "KSUB_RANDOM - Fatal error!\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K is required!\n";
    exit ( 1 );
  }
  else if ( n < k )
  {
    cerr << "\n";
    cerr << "KSUB_RANDOM - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  K <= N is required!\n";
    exit ( 1 );
  }

  a = ( int * ) malloc ( k * sizeof ( int ) );

  if ( k == 0 )
  {
    return a;
  }

  for ( i = 1; i <= k; i++ )
  {
    a[i-1] = ( ( i - 1 ) * n ) / k;
  }

  for ( i = 1; i <= k; i++ )
  {
    for ( ; ; )
    {
      ix = i4_uniform_ab ( 1, n, seed );

      l = 1 + ( ix * k - 1 ) / n;

      if ( a[l-1] < ix )
      {
        break;
      }

    }

    a[l-1] = a[l-1] + 1;

  }

  ip = 0;
  is = k;

  for ( i = 1; i <= k; i++ )
  {
    m = a[i-1];
    a[i-1] = 0;

    if ( m != ( ( i - 1 ) * n ) / k )
    {
      ip = ip + 1;
      a[ip-1] = m;
    }

  }

  ihi = ip;

  for ( i = 1; i <= ihi; i++ )
  {
    ip = ihi + 1 - i;
    l = 1 + ( a[ip-1] * k - 1 ) / n;
    ids = a[ip-1] - ( ( l - 1 ) * n ) / k;
    a[ip-1] = 0;
    a[is-1] = l;
    is = is - ids;
  }

  for ( ll = 1; ll <= k; ll++ )
  {
    l = k + 1 - ll;

    if ( a[l-1] != 0 )
    {
      ir = l;
      m0 = 1 + ( ( a[l-1] - 1 ) * n ) / k;
      m = ( a[l-1] * n ) / k - m0 + 1;
    }

    ix = i4_uniform_ab ( m0, m0 + m - 1, seed );

    i = l + 1;

    while ( i <= ir )
    {
      if ( ix < a[i-1] )
      {
        break;
      }

      ix = ix + 1;
      a[i-2] = a[i-1];
      i = i + 1;
    }
    a[i-2] = ix;
    m = m - 1;
  }
  return a;
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
//    07 April 2014
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
  int i2lo_hi;
  int i2lo_lo;
  int inc;
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

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;

    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j - 1 << ":";
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

int *simplex_grid_index_all ( int m, int n, int ng )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_GRID_INDEX_ALL returns all the simplex grid indices.
//
//  Discussion:
//
//    The number of grid indices can be determined by calling 
//      ng = simplex_grid_size ( m, n )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of subintervals.
//
//    Input, int NG, the number of values in the grid.
//
//    Output, int SIMPLEX_GRID_INDEX_ALL[(M+1)*NG], the current, and then the next,
//    grid index.
//
{
  int *g;
  int *grid;
  int i;
  int k;

  g = new int[m+1];

  for ( i = 0; i < m; i++ )
  {
    g[i] = 0;
  }
  g[m] = n;

  grid = new int[(m+1)*ng];

  k = 0;
  for ( i = 0; i <= m; i++ )
  {
    grid[i+k*(m+1)] = g[i];
  }

  while ( k < ng )
  {
    comp_next_grlex ( m + 1, g );
    k = k + 1;
    for ( i = 0; i <= m; i++ )
    {
      grid[i+k*(m+1)] = g[i];
    }
  }

  delete [] g;

  return grid;
}
//****************************************************************************80

void simplex_grid_index_next ( int m, int n, int g[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_GRID_INDEX_NEXT returns the next simplex grid index.
//
//  Discussion:
//
//    The vector G has dimension M+1.  The first M entries may be regarded
//    as grid coordinates.  These coordinates must have a sum between 0 and N.
//    The M+1 entry contains the remainder, that is N minus the sum of the
//    first M coordinates.
//
//    Each time the function is called, it is given a current grid index, and
//    computes the next one.  The very first index is all zero except for a 
//    final value of N, and the very last index has all zero except for an'
//    intial value of N.
//
//    For example, here are the coordinates in order for M = 3, N = 3:
//
//     0  0 0 0 3
//     1  0 0 1 2
//     2  0 0 2 1
//     3  0 0 3 0
//     4  0 1 0 2
//     5  0 1 1 1
//     6  0 1 2 0
//     7  0 2 0 1
//     8  0 2 1 0
//     9  0 3 0 0
//    10  1 0 0 2
//    11  1 0 1 1
//    12  1 0 2 0
//    13  1 1 0 1
//    14  1 1 1 0
//    15  1 2 0 0
//    16  2 0 0 1
//    17  2 0 1 0
//    18  2 1 0 0
//    19  3 0 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of subintervals.
//
//    Input/output, int G[M+1], the current, and then the next,
//    grid index.
//
{
  comp_next_grlex ( m + 1, g );

  return;
}
//****************************************************************************80

int *simplex_grid_index_sample ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_GRID_INDEX_SAMPLE returns a random simplex grid index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of subintervals in
//    each dimension.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int SIMPLEX_GRID_INDEX_SAMPLE[M+1], a randomly selected index 
//    in the simplex grid.
//
{
  int *g;

  g = comp_random_new ( n, m + 1, seed );

  return g;
}
//****************************************************************************80

double *simplex_grid_index_to_point ( int m, int n, int ng, int g[],
  double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_GRID_INDEX_TO_POINT returns  points corresponding to simplex indices.
//
//  Discussion:
//
//    The M-dimensional simplex is defined by M+1 vertices.
//
//    Given a regular grid that uses N subintervals along the edge between
//    each pair of vertices, a simplex grid index G is a set of M+1 values
//    each between 0 and N, and summing to N. 
//
//    This function determines the coordinates X of the point corresponding
//    to the index G.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of subintervals.
//
//    Input, int NG, the number of grid indices to be converted.
//
//    Input, int G[(M+1)*NG], the grid indices of 1 
//    or more points.
//
//    Input, double V[M*(M+1)], the coordinates of the vertices 
//    of the simplex.
//
//    Output, double SIMPLEX_GRID_INDEX_TO_POINT[M*NG], the coordinates of one 
//    or more points.
//
{
  int i;
  int j;
  int k;
  double *x;
  
  x = new double[m*ng];

  for ( j = 0; j < ng; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = 0.0;
      for ( k = 0; k < m + 1; k++ )
      {
        x[i+j*m] = x[i+j*m] + v[i+k*m] * ( double ) ( g[k+j*(m+1)] );
      }
      x[i+j*m] = x[i+j*m] / ( double ) ( n );
    }
  }

  return x;
}
//****************************************************************************80

int simplex_grid_size ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_GRID_SIZE counts the grid points inside a simplex.
//
//  Discussion:
//
//    The size of a grid with parameters M, N is C(M+N,N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of subintervals.
//
//    Output, int SIMPLEX_GRID_SIZE, the number of grid points.
//
{
  int i;
  int ng;

  ng = 1;

  for ( i = 1; i <= m; i++ )
  {
    ng = ( ng * ( n + i ) ) / i;
  }

  return ng;
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
