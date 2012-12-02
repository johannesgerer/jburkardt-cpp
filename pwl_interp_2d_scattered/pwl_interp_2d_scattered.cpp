# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "pwl_interp_2d_scattered.hpp"
# include "r8lib.hpp"

//****************************************************************************80

int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 )

//****************************************************************************80
//
//  Purpose:
//
//    DIAEDG chooses a diagonal edge.
//
//  Discussion:
//
//    The routine determines whether 0--2 or 1--3 is the diagonal edge
//    that should be chosen, based on the circumcircle criterion, where
//    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
//    quadrilateral in counterclockwise order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double X0, Y0, X1, Y1, X2, Y2, X3, Y3, the coordinates of the
//    vertices of a quadrilateral, given in counter clockwise order.
//
//    Output, int DIAEDG, chooses a diagonal:
//    +1, if diagonal edge 02 is chosen;
//    -1, if diagonal edge 13 is chosen;
//     0, if the four vertices are cocircular.
//
{
  double ca;
  double cb;
  double dx10;
  double dx12;
  double dx30;
  double dx32;
  double dy10;
  double dy12;
  double dy30;
  double dy32;
  double s;
  double tol;
  double tola;
  double tolb;
  int value;

  tol = 100.0 * r8_epsilon ( );

  dx10 = x1 - x0;
  dy10 = y1 - y0;
  dx12 = x1 - x2;
  dy12 = y1 - y2;
  dx30 = x3 - x0;
  dy30 = y3 - y0;
  dx32 = x3 - x2;
  dy32 = y3 - y2;

  tola = tol * r8_max ( fabs ( dx10 ),
               r8_max ( fabs ( dy10 ),
               r8_max ( fabs ( dx30 ), fabs ( dy30 ) ) ) );

  tolb = tol * r8_max ( fabs ( dx12 ),
               r8_max ( fabs ( dy12 ),
               r8_max ( fabs ( dx32 ), fabs ( dy32 ) ) ) );

  ca = dx10 * dx30 + dy10 * dy30;
  cb = dx12 * dx32 + dy12 * dy32;

  if ( tola < ca && tolb < cb )
  {
    value = -1;
  }
  else if ( ca < -tola && cb < -tolb )
  {
    value = 1;
  }
  else
  {
    tola = r8_max ( tola, tolb );
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb
      + ( dx32 * dy12 - dx12 * dy32 ) * ca;

    if ( tola < s )
    {
      value = -1;
    }
    else if ( s < -tola )
    {
      value = 1;
    }
    else
    {
      value = 0;
    }

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
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

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
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

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

void i4vec_heap_d ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_HEAP_D reorders an I4VEC into a descending heap.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    A heap is an array A with the property that, for every index J,
//    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//
//            A(1)         A(2)
//
//      A(3)       A(4)  A(5) A(6)
//
//    A(7) A(8)  A(9) A(10)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the input array.
//
//    Input/output, int A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
//
{
  int i;
  int ifree;
  int key;
  int m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = ( n / 2 ) - 1; 0 <= i; i-- )
  {
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ; ; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }
      }
    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;
  }

  return;
}
//****************************************************************************80

int i4vec_min ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MIN returns the value of the minimum element in an I4VEC.
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
//    Output, int I4VEC_MIN, the value of the minimum element.  This
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
    if ( a[i] < value )
    {
      value = a[i];
    }
  }
  return value;
}
//****************************************************************************80

void i4vec_sort_heap_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
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
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int n1;
  int temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into descending heap form.
//
  i4vec_heap_d ( n, a );
//
//  2: Sort A.
//
//  The largest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    i4vec_heap_d ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }
  return;
}
//****************************************************************************80

int i4vec_sorted_unique ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
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
//    24 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements in A.
//
//    Input/output, int A[N].  On input, the sorted array.  On output, the 
//    unique elements in A.
//
//    Output, int I4VEC_SORTED_UNIQUE, the number of unique elements in A.
//
{
  int i;
  int unique_num;

  unique_num = 0;

  if ( n <= 0 )
  {
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] != a[unique_num-1] )
    {
      unique_num = unique_num + 1;
      a[unique_num-1] = a[i];
    }
  }

  return unique_num;
}
//****************************************************************************80

int lrline ( double xu, double yu, double xv1, double yv1, double xv2,
  double yv2, double dv )

//****************************************************************************80
//
//  Purpose:
//
//    LRLINE determines where a point lies in relation to a directed line.
//
//  Discussion:
//
//    LRLINE determines whether a point is to the left of, right of,
//    or on a directed line parallel to a line through given points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
//    directed line is parallel to and at signed distance DV to the left of
//    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
//    which the position relative to the directed line is to be determined.
//
//    Input, double DV, the signed distance, positive for left.
//
//    Output, int LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
//    to the right of, on, or left of the directed line.  LRLINE is 0 if
//    the line degenerates to a point.
//
{
  double dx;
  double dxu;
  double dy;
  double dyu;
  double t;
  double tol = 0.0000001;
  double tolabs;
  int value;

  dx = xv2 - xv1;
  dy = yv2 - yv1;
  dxu = xu - xv1;
  dyu = yu - yv1;

  tolabs = tol * r8_max ( fabs ( dx ),
                 r8_max ( fabs ( dy ),
                 r8_max ( fabs ( dxu ),
                 r8_max ( fabs ( dyu ), fabs ( dv ) ) ) ) );

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy );

  if ( tolabs < t )
  {
    value = 1;
  }
  else if ( -tolabs <= t )
  {
    value = 0;
  }
  else if ( t < -tolabs )
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int perm_check2 ( int n, int p[], int base )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK2 checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from 1
//    to N occurs among the N entries of the permutation.
//
//    Set the input quantity BASE to 0, if P is a 0-based permutation,
//    or to 1 if P is a 1-based permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the permutation, in standard index form.
//
//    Input, int BASE, the index base.
//
//    Output, int PERM_CHECK, is 1 if the array is NOT a permutation.
//
{
  int error;
  int ifind;
  int iseek;

  error = 0;

  for ( iseek = base; iseek < base + n; iseek++ )
  {
    error = 1;

    for ( ifind = 1; ifind <= n; ifind++ )
    {
      if ( p[ifind-1] == iseek )
      {
        error = 0;
        break;
      }
    }

    if ( error )
    {
      cerr << "\n";
      cerr << "PERM_CHECK2 - Fatal error!\n";
      cerr << "  Could not find occurrence of value " << iseek << "\n";
      return 1;
    }
  }

  return 0;
}
//****************************************************************************80

void perm_inverse ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE inverts a permutation "in place".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On output, P describes the inverse permutation
//
{
  int base;
  int error;
  int i;
  int i0;
  int i1;
  int i2;
  int is;
  int p_min;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE - Fatal error!\n";
    cerr << "  Input value of N = " << n << "\n";
    exit ( 1 );
  }
//
//  Find the least value, and shift data so it begins at 1.
//
  p_min = i4vec_min ( n, p );
  base = 1;

  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - p_min + base;
  }
//
//  Check the permutation.
//
  error = perm_check2 ( n, p, base );

  if ( error )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE - Fatal error!\n";
    cerr << "  The input array does not represent\n";
    cerr << "  a proper permutation.\n";
    exit ( 1 );
  }

  is = 1;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p[i-1];

    while ( i < i1 )
    {
      i2 = p[i1-1];
      p[i1-1] = -i2;
      i1 = i2;
    }

    is = - i4_sign ( p[i-1] );
    p[i-1] = abs ( p[i-1] ) * i4_sign ( is );

  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = -p[i-1];

    if ( 0 <= i1 )
    {
      i0 = i;

      for ( ; ; )
      {
        i2 = p[i1-1];
        p[i1-1] = i0;

        if ( i2 < 0 )
        {
          break;
        }

        i0 = i1;
        i1 = i2;
      }
    }
  }
//
//  Now we can restore the permutation.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + p_min - base;
  }

  return;
}
//****************************************************************************80

double *pwl_interp_2d_scattered_value ( int nd, double xyd[], double zd[], 
  int t_num, int t[], int t_neighbor[], int ni, double xyi[] )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_INTERP_2D_SCATTERED_VALUE evaluates a 2d interpolant of scattered data
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//
//    Input, double XYD[2*ND], the data point coordinates.
//
//    Input, double ZD[ND], the data values.
//
//    Input, int T_NUM, the number of triangles.
//
//    Input, int T[3*T_NUM], the triangle information.
//
//    Input, int T_NEIGHBOR[3*T_NUM], the triangle neighbors.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XYI[2*NI], the interpolation point coordinates.
//
//    Output, double PWL_INTERP_2D_SCATTERED_VALUE[NI], the interpolated values.
//
{
  double alpha;
  double beta;
  double gamma;
  int edge;
  int i;
  int j;
  int step_num;
  double *zi;

  zi = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    triangulation_search_delaunay ( nd, xyd, 3, t_num, t, t_neighbor, 
      xyi+2*i, j, alpha, beta, gamma, edge, step_num );

    if ( j == -1 )
    {
      zi[i] = -1.0;
    }

    zi[i] = alpha * zd[t[0+j*3]] 
          + beta  * zd[t[1+j*3]] 
          + gamma * zd[t[2+j*3]];
  }
  return zi;
}
//****************************************************************************80

int r8tris2 ( int node_num, double node_xy[], int &triangle_num,
  int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8TRIS2 constructs a Delaunay triangulation of 2D vertices.
//
//  Discussion:
//
//    The routine constructs the Delaunay triangulation of a set of 2D vertices
//    using an incremental approach and diagonal edge swaps.  Vertices are
//    first sorted in lexicographically increasing (X,Y) order, and
//    then are inserted one at a time from outside the convex hull.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input/output, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//    On output, the coordinates have been sorted into dictionary order.
//
//    Output, int &TRIANGLE_NUM, the number of triangles in the triangulation;
//    TRIANGLE_NUM is equal to 2*node_num - NB - 2, where NB is the number
//    of boundary vertices.
//
//    Output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up each
//    triangle.  The elements are indices of NODE_XY.  The vertices of the
//    triangles are in counterclockwise order.
//
//    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list.
//    Positive elements are indices of TIL; negative elements are used for links
//    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
//    where I, J = triangle, edge index; TRIANGLE_NEIGHBOR[I,J] refers to
//    the neighbor along edge from vertex J to J+1 (mod 3).
//
//    Output, int R8TRIS2, is 0 for no error.
{
  int base;
  double cmax;
  int e;
  int error;
  int i;
  int *indx;
  int j;
  int k;
  int l;
  int ledg;
  int lr;
  int ltri;
  int m;
  int m1;
  int m2;
  int n;
  int redg;
  int rtri;
  int *stack;
  int t;
  double tol;
  int top;

  stack = new int[node_num];

  tol = 100.0 * r8_epsilon ( );
//
//  Sort the vertices by increasing (x,y).
//
  base = 0;

  indx = r82vec_sort_heap_index_a ( node_num, base, node_xy );

  r82vec_permute ( node_num, indx, base, node_xy );
//
//  Make sure that the nodes are "reasonably" distinct.
//
  m1 = 1;

  for ( i = 2; i <= node_num; i++ )
  {
    m = m1;
    m1 = i;

    k = -1;

    for ( j = 0; j <= 1; j++ )
    {
      cmax = r8_max ( fabs ( node_xy[2*(m-1)+j] ),
                     fabs ( node_xy[2*(m1-1)+j] ) );

      if ( tol * ( cmax + 1.0 )
           < fabs ( node_xy[2*(m-1)+j] - node_xy[2*(m1-1)+j] ) )
      {
        k = j;
        break;
      }

    }

    if ( k == -1 )
    {
      cerr << "\n";
      cerr << "R8TRIS2 - Fatal error!\n";
      cerr << "  Fails for point number I = " << i << "\n";
      cerr << "  M =  " << m  << "\n";
      cerr << "  M1 = " << m1 << "\n";
      cerr << "  X,Y(M)  = " << node_xy[2*(m-1)+0] << "  "
                             << node_xy[2*(m-1)+1] << "\n";
      cerr << "  X,Y(M1) = " << node_xy[2*(m1-1)+0] << "  "
                             << node_xy[2*(m1-1)+1] << "\n";
      exit ( 1 );
    }

  }
//
//  Starting from nodes M1 and M2, search for a third point M that
//  makes a "healthy" triangle (M1,M2,M)
//
  m1 = 1;
  m2 = 2;
  j = 3;

  for ( ; ; )
  {
    if ( node_num < j )
    {
      cerr << "\n";
      cerr << "R8TRIS2 - Fatal error!\n";
      delete [] stack;
      return 225;
    }

    m = j;

    lr = lrline ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1],
      node_xy[2*(m1-1)+0], node_xy[2*(m1-1)+1],
      node_xy[2*(m2-1)+0], node_xy[2*(m2-1)+1], 0.0 );

    if ( lr != 0 )
    {
      break;
    }

    j = j + 1;

  }
//
//  Set up the triangle information for (M1,M2,M), and for any other
//  triangles you created because nodes were collinear with M1, M2.
//
  triangle_num = j - 2;

  if ( lr == -1 )
  {
    triangle_node[3*0+0] = m1;
    triangle_node[3*0+1] = m2;
    triangle_node[3*0+2] = m;
    triangle_neighbor[3*0+2] = -3;

    for ( i = 2; i <= triangle_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      triangle_node[3*(i-1)+0] = m1;
      triangle_node[3*(i-1)+1] = m2;
      triangle_node[3*(i-1)+2] = m;
      triangle_neighbor[3*(i-1)+0] = -3 * i;
      triangle_neighbor[3*(i-1)+1] = i;
      triangle_neighbor[3*(i-1)+2] = i - 1;

    }

    triangle_neighbor[3*(triangle_num-1)+0] = -3 * triangle_num - 1;
    triangle_neighbor[3*(triangle_num-1)+1] = -5;
    ledg = 2;
    ltri = triangle_num;
  }
  else
  {
    triangle_node[3*0+0] = m2;
    triangle_node[3*0+1] = m1;
    triangle_node[3*0+2] = m;
    triangle_neighbor[3*0+0] = -4;

    for ( i = 2; i <= triangle_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      triangle_node[3*(i-1)+0] = m2;
      triangle_node[3*(i-1)+1] = m1;
      triangle_node[3*(i-1)+2] = m;
      triangle_neighbor[3*(i-2)+2] = i;
      triangle_neighbor[3*(i-1)+0] = -3 * i - 3;
      triangle_neighbor[3*(i-1)+1] = i - 1;
    }

    triangle_neighbor[3*(triangle_num-1)+2] = -3 * triangle_num;
    triangle_neighbor[3*0+1] = -3 * triangle_num - 2;
    ledg = 2;
    ltri = 1;

  }
//
//  Insert the vertices one at a time from outside the convex hull,
//  determine visible boundary edges, and apply diagonal edge swaps until
//  Delaunay triangulation of vertices (so far) is obtained.
//
  top = 0;

  for ( i = j+1; i <= node_num; i++ )
  {
    m = i;
    m1 = triangle_node[3*(ltri-1)+ledg-1];

    if ( ledg <= 2 )
    {
      m2 = triangle_node[3*(ltri-1)+ledg];
    }
    else
    {
      m2 = triangle_node[3*(ltri-1)+0];
    }

    lr = lrline ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1],
      node_xy[2*(m1-1)+0], node_xy[2*(m1-1)+1],
      node_xy[2*(m2-1)+0], node_xy[2*(m2-1)+1], 0.0 );

    if ( 0 < lr )
    {
      rtri = ltri;
      redg = ledg;
      ltri = 0;
    }
    else
    {
      l = -triangle_neighbor[3*(ltri-1)+ledg-1];
      rtri = l / 3;
      redg = (l % 3) + 1;
    }

    vbedg ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1], node_num,
      node_xy, triangle_num, triangle_node, triangle_neighbor,
      ltri, ledg, rtri, redg );

    n = triangle_num + 1;
    l = -triangle_neighbor[3*(ltri-1)+ledg-1];

    for ( ; ; )
    {
      t = l / 3;
      e = ( l % 3 ) + 1;
      l = -triangle_neighbor[3*(t-1)+e-1];
      m2 = triangle_node[3*(t-1)+e-1];

      if ( e <= 2 )
      {
        m1 = triangle_node[3*(t-1)+e];
      }
      else
      {
        m1 = triangle_node[3*(t-1)+0];
      }

      triangle_num = triangle_num + 1;
      triangle_neighbor[3*(t-1)+e-1] = triangle_num;
      triangle_node[3*(triangle_num-1)+0] = m1;
      triangle_node[3*(triangle_num-1)+1] = m2;
      triangle_node[3*(triangle_num-1)+2] = m;
      triangle_neighbor[3*(triangle_num-1)+0] = t;
      triangle_neighbor[3*(triangle_num-1)+1] = triangle_num - 1;
      triangle_neighbor[3*(triangle_num-1)+2] = triangle_num + 1;
      top = top + 1;

      if ( node_num < top )
      {
        cerr << "\n";
        cerr << "R8TRIS2 - Fatal error!\n";
        cerr << "  Stack overflow.\n";
        delete [] stack;
        return 8;
      }

      stack[top-1] = triangle_num;

      if ( t == rtri && e == redg )
      {
        break;
      }

    }

    triangle_neighbor[3*(ltri-1)+ledg-1] = -3 * n - 1;
    triangle_neighbor[3*(n-1)+1] = -3 * triangle_num - 2;
    triangle_neighbor[3*(triangle_num-1)+2] = -l;
    ltri = n;
    ledg = 2;

    error = swapec ( m, top, ltri, ledg, node_num, node_xy, triangle_num,
      triangle_node, triangle_neighbor, stack );

    if ( error != 0 )
    {
      cerr << "\n";
      cerr << "R8TRIS2 - Fatal error!\n";
      cerr << "  Error return from SWAPEC.\n";
      delete [] stack;
      return error;
    }

  }
//
//  Now account for the sorting that we did.
//
  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < triangle_num; j++ )
    {
      triangle_node[i+j*3] = indx [ triangle_node[i+j*3] - 1 ];
    }
  }

  perm_inverse ( node_num, indx );

  r82vec_permute ( node_num, indx, base, node_xy );

  delete [] indx;
  delete [] stack;

  return 0;
}
//****************************************************************************80

int swapec ( int i, int &top, int &btri, int &bedg, int node_num,
  double node_xy[], int triangle_num, int triangle_node[],
  int triangle_neighbor[], int stack[] )

//****************************************************************************80
//
//  Purpose:
//
//    SWAPEC swaps diagonal edges until all triangles are Delaunay.
//
//  Discussion:
//
//    The routine swaps diagonal edges in a 2D triangulation, based on
//    the empty circumcircle criterion, until all triangles are Delaunay,
//    given that I is the index of the new vertex added to the triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int I, the index of the new vertex.
//
//    Input/output, int &TOP, the index of the top of the stack.
//    On output, TOP is zero.
//
//    Input/output, int &BTRI, &BEDG; on input, if positive, are the
//    triangle and edge indices of a boundary edge whose updated indices
//    must be recorded.  On output, these may be updated because of swaps.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input/output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the triangle incidence
//    list.  May be updated on output because of swaps.
//
//    Input/output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor
//    list; negative values are used for links of the counter-clockwise linked
//    list of boundary edges;  May be updated on output because of swaps.
//
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Workspace, int STACK[MAXST]; on input, entries 1 through TOP
//    contain the indices of initial triangles (involving vertex I)
//    put in stack; the edges opposite I should be in interior;  entries
//    TOP+1 through MAXST are used as a stack.
//
//    Output, int SWAPEC, is set to 8 for abnormal return.
//
{
  int a;
  int b;
  int c;
  int e;
  int ee;
  int em1;
  int ep1;
  int f;
  int fm1;
  int fp1;
  int l;
  int r;
  int s;
  int swap;
  int t;
  int tt;
  int u;
  double x;
  double y;
//
//  Determine whether triangles in stack are Delaunay, and swap
//  diagonal edge of convex quadrilateral if not.
//
  x = node_xy[2*(i-1)+0];
  y = node_xy[2*(i-1)+1];

  for ( ; ; )
  {
    if ( top <= 0 )
    {
      break;
    }

    t = stack[top-1];
    top = top - 1;

    if ( triangle_node[3*(t-1)+0] == i )
    {
      e = 2;
      b = triangle_node[3*(t-1)+2];
    }
    else if ( triangle_node[3*(t-1)+1] == i )
    {
      e = 3;
      b = triangle_node[3*(t-1)+0];
    }
    else
    {
      e = 1;
      b = triangle_node[3*(t-1)+1];
    }

    a = triangle_node[3*(t-1)+e-1];
    u = triangle_neighbor[3*(t-1)+e-1];

    if ( triangle_neighbor[3*(u-1)+0] == t )
    {
      f = 1;
      c = triangle_node[3*(u-1)+2];
    }
    else if ( triangle_neighbor[3*(u-1)+1] == t )
    {
      f = 2;
      c = triangle_node[3*(u-1)+0];
    }
    else
    {
      f = 3;
      c = triangle_node[3*(u-1)+1];
    }

    swap = diaedg ( x, y,
      node_xy[2*(a-1)+0], node_xy[2*(a-1)+1],
      node_xy[2*(c-1)+0], node_xy[2*(c-1)+1],
      node_xy[2*(b-1)+0], node_xy[2*(b-1)+1] );

    if ( swap == 1 )
    {
      em1 = i4_wrap ( e - 1, 1, 3 );
      ep1 = i4_wrap ( e + 1, 1, 3 );
      fm1 = i4_wrap ( f - 1, 1, 3 );
      fp1 = i4_wrap ( f + 1, 1, 3 );

      triangle_node[3*(t-1)+ep1-1] = c;
      triangle_node[3*(u-1)+fp1-1] = i;
      r = triangle_neighbor[3*(t-1)+ep1-1];
      s = triangle_neighbor[3*(u-1)+fp1-1];
      triangle_neighbor[3*(t-1)+ep1-1] = u;
      triangle_neighbor[3*(u-1)+fp1-1] = t;
      triangle_neighbor[3*(t-1)+e-1] = s;
      triangle_neighbor[3*(u-1)+f-1] = r;

      if ( 0 < triangle_neighbor[3*(u-1)+fm1-1] )
      {
        top = top + 1;
        stack[top-1] = u;
      }

      if ( 0 < s )
      {
        if ( triangle_neighbor[3*(s-1)+0] == u )
        {
          triangle_neighbor[3*(s-1)+0] = t;
        }
        else if ( triangle_neighbor[3*(s-1)+1] == u )
        {
          triangle_neighbor[3*(s-1)+1] = t;
        }
        else
        {
          triangle_neighbor[3*(s-1)+2] = t;
        }

        top = top + 1;

        if ( node_num < top )
        {
          return 8;
        }

        stack[top-1] = t;
      }
      else
      {
        if ( u == btri && fp1 == bedg )
        {
          btri = t;
          bedg = e;
        }

        l = - ( 3 * t + e - 1 );
        tt = t;
        ee = em1;

        while ( 0 < triangle_neighbor[3*(tt-1)+ee-1] )
        {
          tt = triangle_neighbor[3*(tt-1)+ee-1];

          if ( triangle_node[3*(tt-1)+0] == a )
          {
            ee = 3;
          }
          else if ( triangle_node[3*(tt-1)+1] == a )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }

        }

        triangle_neighbor[3*(tt-1)+ee-1] = l;

      }

      if ( 0 < r )
      {
        if ( triangle_neighbor[3*(r-1)+0] == t )
        {
          triangle_neighbor[3*(r-1)+0] = u;
        }
        else if ( triangle_neighbor[3*(r-1)+1] == t )
        {
          triangle_neighbor[3*(r-1)+1] = u;
        }
        else
        {
          triangle_neighbor[3*(r-1)+2] = u;
        }
      }
      else
      {
        if ( t == btri && ep1 == bedg )
        {
          btri = u;
          bedg = f;
        }

        l = - ( 3 * u + f - 1 );
        tt = u;
        ee = fm1;

        while ( 0 < triangle_neighbor[3*(tt-1)+ee-1] )
        {
          tt = triangle_neighbor[3*(tt-1)+ee-1];

          if ( triangle_node[3*(tt-1)+0] == b )
          {
            ee = 3;
          }
          else if ( triangle_node[3*(tt-1)+1] == b )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }

        }

        triangle_neighbor[3*(tt-1)+ee-1] = l;
      }
    }
  }

  return 0;
}
//****************************************************************************80

void triangulation_order3_print ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_PRINT prints information defining a triangulation.
//
//  Discussion:
//
//    Triangulations created by R8TRIS2 include extra information encoded
//    in the negative values of TRIANGLE_NEIGHBOR.
//
//    Because some of the nodes counted in NODE_NUM may not actually be
//    used in the triangulation, I needed to compute the true number
//    of vertices.  I added this calculation on 13 October 2001.
//
//    Ernest Fasse pointed out an error in the indexing of VERTEX_LIST,
//    which was corrected on 19 February 2004.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
//    the triangles.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
//    on each side.  If there is no triangle neighbor on a particular side,
//    the value of TRIANGLE_NEIGHBOR should be negative.  If the
//    triangulation data was created by R8TRIS2, then there is more
//    information encoded in the negative values.
//
{
# define DIM_NUM 2

  int boundary_num;
  int i;
  int j;
  int k;
  int n1;
  int n2;
  int s;
  int s1;
  int s2;
  bool skip;
  int t;
  int *vertex_list;
  int vertex_num;

  cout << "\n";
  cout << "TRIANGULATION_ORDER3_PRINT\n";
  cout << "  Information defining a triangulation.\n";
  cout << "\n";
  cout << "  The number of nodes is " << node_num << "\n";

  r8mat_transpose_print ( DIM_NUM, node_num, node_xy, "  Node coordinates" );

  cout << "\n";
  cout << "  The number of triangles is " << triangle_num << "\n";
  cout << "\n";
  cout << "  Sets of three nodes are used as vertices of\n";
  cout << "  the triangles.  For each triangle, the nodes\n";
  cout << "  are listed in counterclockwise order.\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_node, "  Triangle nodes" );

  cout << "\n";
  cout << "  On each side of a given triangle, there is either\n";
  cout << "  another triangle, or a piece of the convex hull.\n";
  cout << "  For each triangle, we list the indices of the three\n";
  cout << "  neighbors, or (if negative) the codes of the\n";
  cout << "  segments of the convex hull.\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_neighbor,
    "  Triangle neighbors" );
//
//  Determine VERTEX_NUM, the number of vertices.
//
  vertex_list = new int[3*triangle_num];

  k = 0;
  for ( t = 0; t < triangle_num; t++ )
  {
    for ( s = 0; s < 3; s++ )
    {
      vertex_list[k] = triangle_node[s+t*3];
      k = k + 1;
    }
  }

  i4vec_sort_heap_a ( 3*triangle_num, vertex_list );

  vertex_num = i4vec_sorted_unique ( 3*triangle_num, vertex_list );

  delete [] vertex_list;
//
//  Determine the number of boundary points.
//
  boundary_num = 2 * vertex_num - triangle_num - 2;

  cout << "\n";
  cout << "  The number of boundary points is " << boundary_num << "\n";
  cout << "\n";
  cout << "  The segments that make up the convex hull can be\n";
  cout << "  determined from the negative entries of the triangle\n";
  cout << "  neighbor list.\n";
  cout << "\n";
  cout << "     #   Tri  Side    N1    N2\n";
  cout << "\n";

  skip = false;

  k = 0;

  for ( i = 0; i < triangle_num; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      if ( triangle_neighbor[j+i*3] < 0 )
      {
        s = -triangle_neighbor[j+i*3];
        t = s / 3;

        if ( t < 1 || triangle_num < t )
        {
          cout << "\n";
          cout << "  Sorry, this data does not use the R8TRIS2\n";
          cout << "  convention for convex hull segments.\n";
          skip = true;
          break;
        }

        s1 = ( s % 3 ) + 1;
        s2 = i4_wrap ( s1+1, 1, 3 );
        k = k + 1;
        n1 = triangle_node[s1-1+(t-1)*3];
        n2 = triangle_node[s2-1+(t-1)*3];
        cout                  << "  "
             << setw(4) << k  << "  "
             << setw(4) << t  << "  "
             << setw(4) << s1 << "  "
             << setw(4) << n1 << "  "
             << setw(4) << n2 << "\n";
      }
    }

    if ( skip )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void triangulation_search_delaunay ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[],
  int triangle_neighbor[], double p[2], int &triangle_index, 
  double &alpha, double &beta, double &gamma, int &edge,
  int &step_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_SEARCH_DELAUNAY searches a triangulation for a point.
//
//  Discussion:
//
//    The algorithm "walks" from one triangle to its neighboring triangle,
//    and so on, until a triangle is found containing point P, or P is found
//    to be outside the convex hull.
//
//    The algorithm computes the barycentric coordinates of the point with
//    respect to the current triangle.  If all three quantities are positive,
//    the point is contained in the triangle.  If the I-th coordinate is
//    negative, then (X,Y) lies on the far side of edge I, which is opposite
//    from vertex I.  This gives a hint as to where to search next.
//
//    For a Delaunay triangulation, the search is guaranteed to terminate.
//    For other triangulations, a cycle may occur.
//
//    Note the surprising fact that, even for a Delaunay triangulation of
//    a set of nodes, the nearest point to (X,Y) need not be one of the
//    vertices of the triangle containing (X,Y).
//
//    The code can be called for triangulations of any order, but only
//    the first three nodes in each triangle are considered.  Thus, if
//    higher order triangles are used, and the extra nodes are intended
//    to give the triangle a polygonal shape, these will have no effect,
//    and the results obtained here might be misleading.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2012
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the nodes of each triangle.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list.
//
//    Input, double P[2], the coordinates of a point.
//
//    Output, int &TRIANGLE_INDEX, the index of the triangle where the search ended.
//    If a cycle occurred, then TRIANGLE_INDEX = -1.
//
//    Output, double &ALPHA, &BETA, &GAMMA, the barycentric coordinates
//    of the point with respect to triangle TRIANGLE_INDEX.
//
//    Output, int &EDGE, indicates the position of the point (X,Y) in
//    triangle TRIANGLE:
//    0, the interior or boundary of the triangle;
//    -1, outside the convex hull of the triangulation, past edge 1;
//    -2, outside the convex hull of the triangulation, past edge 2;
//    -3, outside the convex hull of the triangulation, past edge 3.
//
//    Output, int &STEP_NUM, the number of steps.
{
  int a;
  int b;
  int c;
  double det;
  double dxp;
  double dxa;
  double dxb;
  double dyp;
  double dya;
  double dyb;
  static int triangle_index_save = -1;

  step_num = - 1;
  edge = 0;

  if ( triangle_index_save < 0 || triangle_num <= triangle_index_save )
  {
    triangle_index = ( triangle_num + 1 ) / 2;
  }
  else
  {
    triangle_index = triangle_index_save;
  }

  for ( ; ; )
  {
    step_num = step_num + 1;

    if ( triangle_num < step_num )
    {
      cerr << "\n";
      cerr << "TRIANGULATION_SEARCH_DELAUNAY - Fatal error!\n";
      cerr << "  The algorithm seems to be cycling.\n";
      cerr << "  Current triangle is " << triangle_index << "\n";
      exit ( 1 );
    }
//
//  Get the vertices of triangle TRIANGLE.
//
    a = triangle_node[0+triangle_index*triangle_order];
    b = triangle_node[1+triangle_index*triangle_order];
    c = triangle_node[2+triangle_index*triangle_order];
//
//  Using vertex C as a base, compute the distances to vertices A and B,
//  and the point (X,Y).
//
    dxa = node_xy[0+a*2] - node_xy[0+c*2];
    dya = node_xy[1+a*2] - node_xy[1+c*2];

    dxb = node_xy[0+b*2] - node_xy[0+c*2];
    dyb = node_xy[1+b*2] - node_xy[1+c*2];

    dxp = p[0]           - node_xy[0+c*2];
    dyp = p[1]           - node_xy[1+c*2];

    det = dxa * dyb - dya * dxb;
//
//  Compute the barycentric coordinates of the point (X,Y) with respect
//  to this triangle.
//
    alpha = ( dxp * dyb - dyp * dxb ) / det;
    beta =  ( dxa * dyp - dya * dxp ) / det;
    gamma = 1.0 - alpha - beta;
//
//  If the barycentric coordinates are all positive, then the point
//  is inside the triangle and we're done.
//
    if ( 0.0 <= alpha &&
         0.0 <= beta  &&
         0.0 <= gamma )
    {
      break;
    }
//
//  At least one barycentric coordinate is negative.
//
//  If there is a negative barycentric coordinate for which there exists
//  an opposing triangle neighbor closer to the point, move to that triangle.
//
//  (Two coordinates could be negative, in which case we could go for the
//  most negative one, or the most negative one normalized by the actual
//  distance it represents).
//
    if ( alpha < 0.0 && 0 <= triangle_neighbor[1+triangle_index*3] )
    {
      triangle_index = triangle_neighbor[1+triangle_index*3];
      continue;
    }
    else if ( beta < 0.0 && 0 <= triangle_neighbor[2+triangle_index*3] )
    {
      triangle_index = triangle_neighbor[2+triangle_index*3];
      continue;
    }
    else if ( gamma < 0.0 && 0 <= triangle_neighbor[0+triangle_index*3] )
    {
      triangle_index = triangle_neighbor[0+triangle_index*3];
      continue;
    }
//
//  All negative barycentric coordinates correspond to vertices opposite
//  sides on the convex hull.
//
//  Note the edge and exit.
//
    if ( alpha < 0.0 )
    {
      edge = -2;
      break;
    }
    else if ( beta < 0.0 )
    {
      edge = -3;
      break;
    }
    else if ( gamma < 0.0 )
    {
      edge = -1;
      break;
    }
    else
    {
      cerr << "\n";
      cerr << "TRIANGULATION_ORDER3_SEARCH - Fatal error!\n";
      cerr << "  The algorithm seems to have reached a dead end\n";
      cerr << "  after " << step_num << " steps.\n";
      exit ( 1 );
    }
  }
  triangle_index_save = triangle_index;

  return;
}
//****************************************************************************80

void vbedg ( double x, double y, int node_num, double node_xy[],
  int triangle_num, int triangle_node[], int triangle_neighbor[], int &ltri,
  int &ledg, int &rtri, int &redg )

//****************************************************************************80
//
//  Purpose:
//
//    VBEDG determines which boundary edges are visible to a point.
//
//  Discussion:
//
//    The point (X,Y) is assumed to be outside the convex hull of the
//    region covered by the 2D triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2008
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double X, Y, the coordinates of a point outside the convex hull
//    of the current triangulation.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the triangle incidence list.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list;
//    negative values are used for links of a counter clockwise linked list
//    of boundary edges;
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Input/output, int &LTRI, &LEDG.  If LTRI != 0 then these values are
//    assumed to be already computed and are not changed, else they are updated.
//    On output, LTRI is the index of boundary triangle to the left of the
//    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
//    edge of triangle LTRI to the left of the leftmost boundary edge visible
//    from (X,Y).  1 <= LEDG <= 3.
//
//    Input/output, int &RTRI.  On input, the index of the boundary triangle
//    to begin the search at.  On output, the index of the rightmost boundary
//    triangle visible from (X,Y).
//
//    Input/output, int &REDG, the edge of triangle RTRI that is visible
//    from (X,Y).  1 <= REDG <= 3.
//
{
  int a;
  double ax;
  double ay;
  int b;
  double bx;
  double by;
  bool done;
  int e;
  int l;
  int lr;
  int t;
//
//  Find the rightmost visible boundary edge using links, then possibly
//  leftmost visible boundary edge using triangle neighbor information.
//
  if ( ltri == 0 )
  {
    done = false;
    ltri = rtri;
    ledg = redg;
  }
  else
  {
    done = true;
  }

  for ( ; ; )
  {
    l = -triangle_neighbor[3*(rtri-1)+redg-1];
    t = l / 3;
    e = 1 + l % 3;
    a = triangle_node[3*(t-1)+e-1];

    if ( e <= 2 )
    {
      b = triangle_node[3*(t-1)+e];
    }
    else
    {
      b = triangle_node[3*(t-1)+0];
    }

    ax = node_xy[2*(a-1)+0];
    ay = node_xy[2*(a-1)+1];

    bx = node_xy[2*(b-1)+0];
    by = node_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

    rtri = t;
    redg = e;

  }

  if ( done )
  {
    return;
  }

  t = ltri;
  e = ledg;

  for ( ; ; )
  {
    b = triangle_node[3*(t-1)+e-1];
    e = i4_wrap ( e-1, 1, 3 );

    while ( 0 < triangle_neighbor[3*(t-1)+e-1] )
    {
      t = triangle_neighbor[3*(t-1)+e-1];

      if ( triangle_node[3*(t-1)+0] == b )
      {
        e = 3;
      }
      else if ( triangle_node[3*(t-1)+1] == b )
      {
        e = 1;
      }
      else
      {
        e = 2;
      }

    }

    a = triangle_node[3*(t-1)+e-1];
    ax = node_xy[2*(a-1)+0];
    ay = node_xy[2*(a-1)+1];

    bx = node_xy[2*(b-1)+0];
    by = node_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

  }

  ltri = t;
  ledg = e;

  return;
}