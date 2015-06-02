# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "tet_mesh.hpp"

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

void i4_swap ( int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = *i;
  *i = *j;
  *j = k;
 
  return;
}
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4.
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
//    Output, int I4_UNIFORM_AB, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
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

int i4col_compare ( int m, int n, int a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_COMPARE compares columns I and J of an I4COL.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 4
//
//      A = (
//        1  2  3  4
//        5  6  7  8
//        9 10 11 12 )
//
//    Output:
//
//      I4COL_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], an array of N columns of vectors of length M.
//
//    Input, int I, J, the columns to be compared.
//    I and J must be between 1 and N.
//
//    Output, int I4COL_COMPARE, the results of the comparison:
//    -1, column I < column J,
//     0, column I = column J,
//    +1, column J < column I.
//
{
  int k;
//
//  Check.
//
  if ( i < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index I = " << i << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < i )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index I = " << i << ".\n";
    exit ( 1 );
  }

  if ( j < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index J = " << j << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < j )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index J = " << j << ".\n";
    exit ( 1 );
  }

  if ( i == j )
  {
    return 0;
  }

  k = 1;

  while ( k <= m )
  {
    if ( a[k-1+(i-1)*m] < a[k-1+(j-1)*m] )
    {
      return (-1);
    }
    else if ( a[k-1+(j-1)*m] < a[k-1+(i-1)*m] )
    {
      return 1;
    }
    k = k + 1;
  }

  return 0;
}
//****************************************************************************80

void i4col_sort_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT_A ascending sorts the columns of an I4COL.
//
//  Discussion:
//
//    In lexicographic order, the statement "X < Y", applied to two
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, X is less than Y if, at the first index where they
//    differ, the X value is less than the Y value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors;
//    On output, the columns of A have been sorted in ascending
//    lexicographic order.
//
{
  int i;
  int indx;
  int isgn;
  int j;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      i4col_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

void i4col_sort2_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A, and the length
//    of a vector of data.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors.
//    On output, the elements of each column of A have been sorted in ascending
//    order.
//
{
  int col;
  int i;
  int indx;
  int isgn;
  int j;
  int row;
  int temp;

  if ( m <= 1 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }
//
//  Initialize.
//
  for ( col = 0; col < n; col++ )
  {
    i = 0;
    indx = 0;
    isgn = 0;
    j = 0;
//
//  Call the external heap sorter.
//
    for ( ; ; )
    {
      sort_heap_external ( m, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
      if ( 0 < indx )
      {
        temp         = a[i-1+col*m];
        a[i-1+col*m] = a[j-1+col*m];
        a[j-1+col*m] = temp;
      }
//
//  Compare the I and J objects.
//
      else if ( indx < 0 )
      {
        if ( a[j-1+col*m] < a[i-1+col*m] )
        {
          isgn = +1;
        }
        else
        {
          isgn = -1;
        }
      }
      else if ( indx == 0 )
      {
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

int i4col_sorted_unique_count ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
//
//  Discussion:
//
//    The columns of the array may be ascending or descending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], a sorted array, containing
//    N columns of data.
//
//    Output, int I4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
//
{
  int i;
  int j1;
  int j2;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  unique_num = 1;
  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        unique_num = unique_num + 1;
        j1 = j2;
        break;
      }
    }
  }

  return unique_num;
}
//****************************************************************************80

void i4col_swap ( int m, int n, int a[], int icol1, int icol2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SWAP swaps two columns of an I4COL.
//
//  Discussion:
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based!  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], an array of data.
//
//    Input, int ICOL1, ICOL2, the two columns to swap.
//    These indices should be between 1 and N.
//
{
# define OFFSET 1

  int i;
  int t;
//
//  Check.
//
  if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL1 is out of range.\n";
    exit ( 1 );
  }

  if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL2 is out of range.\n";
    exit ( 1 );
  }

  if ( icol1 == icol2 )
  {
    return;
  }
  for ( i = 0; i < m; i++ )
  {
    t                     = a[i+(icol1-OFFSET)*m];
    a[i+(icol1-OFFSET)*m] = a[i+(icol2-OFFSET)*m];
    a[i+(icol2-OFFSET)*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

void i4i4_sort_a ( int i1, int i2, int *j1, int *j2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4I4_SORT_A ascending sorts a pair of I4's.
//
//  Discussion:
//
//    The program allows the reasonable call:
//
//      i4i4_sort_a ( i1, i2, &i1, &i2 );
//
//    and this will return the reasonable result.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, the values to sort.
//
//    Output, int J1, J2, the sorted values.
//
{
  int k1;
  int k2;
//
//  Copy arguments, so that the user can make "reasonable" calls like:
//
//    i4i4_sort_a ( i1, i2, &i1, &i2 );
//
  k1 = i1;
  k2 = i2;

  *j1 = i4_min ( k1, k2 );
  *j2 = i4_max ( k1, k2 );

  return;
}
//****************************************************************************80

void i4i4i4_sort_a ( int i1, int i2, int i3, int *j1, int *j2, int *j3 )

//****************************************************************************80
//
//  Purpose:
//
//    I4I4I4_SORT_A ascending sorts a triple of I4's.
//
//  Discussion:
//
//    The program allows the reasonable call:
//
//      i4i4i4_sort_a ( i1, i2, i3, &i1, &i2, &i3 );
//
//    and this will return the reasonable result.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, I3, the values to sort.
//
//    Output, int *J1, *J2, *J3, the sorted values.
//
{
  int k1;
  int k2;
  int k3;
//
//  Copy arguments, so that the user can make "reasonable" calls like:
//
//    i4i4i4_sort_a ( i1, i2, i3, &i1, &i2, &i3 );
//
  k1 = i1;
  k2 = i2;
  k3 = i3;

  *j1 = i4_min ( i4_min ( k1, k2 ), i4_min ( k2, k3 ) );
  *j2 = i4_min ( i4_max ( k1, k2 ), 
        i4_min ( i4_max ( k2, k3 ), i4_max ( k3, k1 ) ) );
  *j3 = i4_max ( i4_max ( k1, k2 ), i4_max ( k2, k3 ) );

  return;
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
//    Input, string TITLE, a title to be printed.
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
//    14 June 2005
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
//    Input, string TITLE, a title for the matrix.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
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
      cout << setw(6) << i << "  ";
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
      cout << setw(5) << j << "  ";
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
//    Input, string TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  " << setw(8) << i 
         << "  " << setw(8) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

void i4vec_zero ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO zeroes an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}
//****************************************************************************80

void mesh_base_one ( int node_num, int element_order, int element_num,
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_BASE_ONE ensures that the element definition is 1-based.
//
//  Discussion:
//
//    The ELEMENT_NODE array contains nodes indices that form elements.
//    The convention for node indexing might start at 0 or at 1.
//
//    If this function detects 0-based indexing, it converts to 1-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
//    definitions.
//
{
  int element;
  const int i4_huge = 2147483647;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = + i4_huge;
  node_max = - i4_huge;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      if ( node < node_min )
      {
        node_min = node;
      }
      if ( node_max < node )
      {
        node_max = node;
      }
    }
  }
  if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ONE:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  This will be converted to 1-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] + 1;
      }
    }
  }
  else if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ONE:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else
  {
    cout << "\n";
    cout << "MESH_BASE_ONE - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
    cout << "  NODE_MIN = " << node_min << "\n";
    cout << "  NODE_MAX = " << node_max << "\n";
    cout << "  NODE_NUM = " << node_num << "\n";
  }
  return;
}
//****************************************************************************80

void mesh_base_zero ( int node_num, int element_order, int element_num,
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_BASE_ZERO ensures that the element definition is zero-based.
//
//  Discussion:
//
//    The ELEMENT_NODE array contains nodes indices that form elements.
//    The convention for node indexing might start at 0 or at 1.
//    Since a C++ program will naturally assume a 0-based indexing, it is
//    necessary to check a given element definition and, if it is actually
//    1-based, to convert it.
//
//    This function attempts to detect 1-based node indexing and correct it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
//    definitions.
//
{
  int element;
  const int i4_huge = 2147483647;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = + i4_huge;
  node_max = - i4_huge;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      if ( node < node_min )
      {
        node_min = node;
      }
      if ( node_max < node )
      {
        node_max = node;
      }
    }
  }

  if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  This will be converted to 0-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] - 1;
      }
    }
  }
  else
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
    cout << "  NODE_MIN = " << node_min << "\n";
    cout << "  NODE_MAX = " << node_max << "\n";
    cout << "  NODE_NUM = " << node_num << "\n";
  }
  return;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

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
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal double precision number, 
//    and is usually defined in math.h, or sometimes in stdlib.h.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  return HUGE_VAL;
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

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

void r8_swap ( double *x, double *y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SWAP switches two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, double *X, *Y.  On output, the values of X and
//    Y have been interchanged.
//
{
  double z;

  z = *x;
  *x = *y;
  *y = z;
 
  return;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

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
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      r8_uniform_01 = seed / ( 2^31 - 1 )
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
//    11 August 2004
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
//    in Handbook of Simulation
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
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double r8mat_det_4d ( double a[4*4] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
//
//  Discussion:
//
//    The two dimensional array is stored as a one dimensional vector,
//    by COLUMNS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[4*4], the matrix whose determinant is desired.
//
//    Output, double R8MAT_DET_4D, the determinant of the matrix.
//
{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}
//****************************************************************************80

double *r8mat_mv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV_NEW multiplies a matrix times a vector.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8MAT_MV_NEW[M], the product A*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT, with an optional title.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
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
//    29 August 2003
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
//    Input, string TITLE, a title to be printed.
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
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2004
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
//    Input, string TITLE, a title for the matrix.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of 5.
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
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
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
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i;
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << "  " << setprecision(6) << setw(12) << a[i-1+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

int r8mat_solve ( int n, int rhs_num, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
//
//  Discussion: 							    
//
//    A R8MAT is a doubly dimensioned array of double precision values, which
//    may be stored as a vector in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*N]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int RHS_NUM, the number of right hand sides.  RHS_NUM
//    must be at least 0.
//
//    Input/output, double A[N*(N+RHS_NUM)], contains in rows and columns 1
//    to N the coefficient matrix, and in columns N+1 through
//    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
//    area has been destroyed, while the right hand sides have
//    been overwritten with the corresponding solutions.
//
//    Output, int R8MAT_SOLVE, singularity flag.
//    0, the matrix was not singular, the solutions were computed;
//    J, factorization failed on step J, and the solutions could not
//    be computed.
//
{
  double apivot;
  double factor;
  int i;
  int ipivot;
  int j;
  int k;
  double temp;

  for ( j = 0; j < n; j++ )
  {
//
//  Choose a pivot row.
//
    ipivot = j;
    apivot = a[j+j*n];

    for ( i = j; i < n; i++ )
    {
      if ( fabs ( apivot ) < fabs ( a[i+j*n] ) )
      {
        apivot = a[i+j*n];
        ipivot = i;
      }
    }

    if ( apivot == 0.0 )
    {
      return j;
    }
//
//  Interchange.
//
    for ( i = 0; i < n + rhs_num; i++ )
    {
      temp          = a[ipivot+i*n];
      a[ipivot+i*n] = a[j+i*n];
      a[j+i*n]      = temp;
    }
//
//  A(J,J) becomes 1.
//
    a[j+j*n] = 1.0;
    for ( k = j; k < n + rhs_num; k++ )
    {
      a[j+k*n] = a[j+k*n] / apivot;
    }
//
//  A(I,J) becomes 0.
//
    for ( i = 0; i < n; i++ )
    {
      if ( i != j )
      {
        factor = a[i+j*n];
        a[i+j*n] = 0.0;
        for ( k = j; k < n + rhs_num; k++ )
        {
          a[i+k*n] = a[i+k*n] - factor * a[j+k*n];
        }
      }
    }
  }

  return 0;
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
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
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
//    Input, string TITLE, an optional title.
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
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
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
//    Input, string TITLE, an optional title.
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

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

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
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setprecision ( 6 ) << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *r8mat_uniform_01_new ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW fills a double precision array with unit pseudorandom values.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
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
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number//
//
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

double *r8vec_cross_3d ( double v1[3], double v2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CROSS_3D computes the cross product of two vectors in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], the coordinates of the vectors.
//
//    Output, double R8VEC_CROSS_3D[3], the cross product vector.
//
{
  double *v3;

  v3 = new double[3];

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}
//****************************************************************************80

bool r8vec_is_nonnegative ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_IS_NONNEGATIVE is true if all entries in an R8VEC are nonnegative.
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
//    04 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[N], the vector to be checked.
//
//    Output, bool R8VEC_IS_NONNEGATIVE is true if all elements of X
//    are nonnegative.
//
{
  int i;

  for ( i = 0; i < n; i++ ) 
  {
    if ( x[i] < 0.0 ) 
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

bool r8vec_is_zero ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_IS_ZERO is true if the entries in an R8VEC are all zero.
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
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[N], the vector to be checked.
//
//    Output, bool R8VEC_IS_ZERO is true if all N elements of X
//    are zero.
//
{
  int i;

  for ( i = 0; i < n; i++ ) 
  {
    if ( x[i] != 0.0 ) 
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

double r8vec_length ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LENGTH returns the Euclidean length of a R8VEC
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the vector.
//
//    Output, double R8VEC_LENGTH, the Euclidean length of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    value = value + pow ( x[i], 2 );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

double r8vec_max ( int n, double dvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the maximum element in an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double DVEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double *r8vec_pointer;
  double value;

  value = - r8_huge ( );

  if ( n <= 0 ) 
  {
    return value;
  }

  for ( i = 0; i < n; i++ ) 
  {
    if ( value < dvec[i] )
    {
      value = dvec[i];
    }
  }
  return value;
}
//****************************************************************************80

double r8vec_mean ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MEAN returns the mean of a R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[N], the vector whose mean is desired.
//
//    Output, double R8VEC_MEAN, the mean, or average, of the vector entries.
//
{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ ) 
  {
    mean = mean + x[i];
  }

  mean = mean / ( double ) n;

  return mean;
}
//****************************************************************************80

double r8vec_min ( int n, double dvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the minimum element in an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double DVEC[N], the array to be checked.
//
//    Output, double R8VEC_MIN, the value of the minimum element.
//
{
  int i;
  double *r8vec_pointer;
  double value;

  value = r8_huge ( );

  if ( n <= 0 ) 
  {
    return value;
  }

  for ( i = 0; i < n; i++ ) 
  {
    if ( dvec[i] < value )
    {
      value = dvec[i];
    }
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
//    Input, string TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(14) << a[i]  << "\n";
  }

  return;
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

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a unit pseudorandom R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
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
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
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
  int k;
  double *r;

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double r8vec_variance ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_VARIANCE returns the variance of a double vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[N], the vector whose variance is desired.
//
//    Output, double R8VEC_VARIANCE, the variance of the vector entries.
//
{
  int i;
  double mean;
  double variance;

  mean = r8vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ ) 
  {
    variance = variance + ( x[i] - mean ) * ( x[i] - mean );
  }

  if ( 1 < n ) 
  {
    variance = variance / ( double ) ( n - 1 );
  }
  else
  {
    variance = 0.0;
  }

  return variance;
}
//****************************************************************************80

void r8vec_zero ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO zeroes a real vector.
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
//    Input, int N, the number of entries in the vector.
//
//    Output, double A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n ) 
  {
    if ( s[n-1] != ' ' && s[n-1] != '\n' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
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
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 ) 
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
//****************************************************************************80

int *tet_mesh_neighbor_tets ( int tetra_order, int tetra_num, 
  int tetra_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_NEIGHBOR_TETS determines tetrahedron neighbors.
//
//  Discussion:
//
//    A tet mesh of a set of nodes can be completely described by
//    the coordinates of the nodes, and the list of nodes that make up
//    each tetrahedron.  In the most common case, four nodes are used.
//    There is also a 10 node case, where nodes are also placed on
//    the midsides of the tetrahedral edges.
//
//    This routine can handle 4 or 10-node tetrahedral meshes.  The
//    10-node case is handled simply by ignoring the six midside nodes,
//    which are presumed to be listed after the vertices.
//
//    The tetrahedron adjacency information records which tetrahedron
//    is adjacent to a given tetrahedron on a particular face.
//
//    This routine creates a data structure recording this information.
//
//    The primary amount of work occurs in sorting a list of 4 * TETRA_NUM
//    data items.
//
//    The neighbor tetrahedrons are indexed by the face they share with
//    the tetrahedron.
//
//    Each face of the tetrahedron is indexed by the node which is NOT
//    part of the face.  That is:
//
//    * Neighbor 1 shares face 1 defined by nodes 2, 3, 4.
//    * Neighbor 2 shares face 2 defined by nodes 1, 3, 4;
//    * Neighbor 3 shares face 3 defined by nodes 1, 2, 4;
//    * Neighbor 4 shares face 4 defined by nodes 1, 2, 3.
//
//    For instance, if the (transposed) TETRA_NODE array was:
//
//    Row       1      2      3      4
//    Col
//
//      1       4      3      5      1
//      2       4      2      5      1
//      3       4      7      3      5
//      4       4      7      8      5
//      5       4      6      2      5
//      6       4      6      8      5
//
//    then the (transposed) TETRA_NEIGHBOR array should be:
//
//    Row       1      2      3      4
//    Col
//
//      1      -1      2     -1      3
//      2      -1      1     -1      5
//      3      -1      1      4     -1
//      4      -1      6      3     -1
//      5      -1      2      6     -1
//      6      -1      4      5     -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TETRA_ORDER, the order of the tetrahedrons.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
//
//    Output, int TET_MESH_NEIGHBORS[4*TETRA_NUM], the four tetrahedrons that
//    are direct neighbors of a given tetrahedron.  If there is no neighbor
//    sharing a given face, the index is set to -1.
//
{
  int a;
  int b;
  int c;
  int face;
  int face1;
  int face2;
  int *faces;
  int i;
  int j;
  int k;
  int l;
  int tetra;
  int *tetra_neighbor;
  int tetra1;
  int tetra2;

  faces = new int[5*(4*tetra_num)];
  tetra_neighbor = new int[4*tetra_num];
//
//  Step 1.
//  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
//  construct the four face relations:
//
//    (J,K,L,1,T)
//    (I,K,L,2,T)
//    (I,J,L,3,T)
//    (I,J,K,4,T)
//
//  In order to make matching easier, we reorder each triple of nodes
//  into ascending order.
//
  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    i = tetra_node[0+tetra*tetra_order];
    j = tetra_node[1+tetra*tetra_order];
    k = tetra_node[2+tetra*tetra_order];
    l = tetra_node[3+tetra*tetra_order];

    i4i4i4_sort_a ( j, k, l, &a, &b, &c );

    faces[0+0*5+tetra*5*4] = a;
    faces[1+0*5+tetra*5*4] = b;
    faces[2+0*5+tetra*5*4] = c;
    faces[3+0*5+tetra*5*4] = 0;
    faces[4+0*5+tetra*5*4] = tetra;

    i4i4i4_sort_a ( i, k, l, &a, &b, &c );

    faces[0+1*5+tetra*5*4] = a;
    faces[1+1*5+tetra*5*4] = b;
    faces[2+1*5+tetra*5*4] = c;
    faces[3+1*5+tetra*5*4] = 1;
    faces[4+1*5+tetra*5*4] = tetra;

    i4i4i4_sort_a ( i, j, l, &a, &b, &c );

    faces[0+2*5+tetra*5*4] = a;
    faces[1+2*5+tetra*5*4] = b;
    faces[2+2*5+tetra*5*4] = c;
    faces[3+2*5+tetra*5*4] = 2;
    faces[4+2*5+tetra*5*4] = tetra;

    i4i4i4_sort_a ( i, j, k, &a, &b, &c );

    faces[0+3*5+tetra*5*4] = a;
    faces[1+3*5+tetra*5*4] = b;
    faces[2+3*5+tetra*5*4] = c;
    faces[3+3*5+tetra*5*4] = 3;
    faces[4+3*5+tetra*5*4] = tetra;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:3; the routine we call here
//  sorts on rows 1 through 5 but that won't hurt us.
//
//  What we need is to find cases where two tetrahedrons share a face.
//  By sorting the columns of the FACES array, we will put shared faces
//  next to each other.
//
  i4col_sort_a ( 5, 4*tetra_num, faces );
//
//  Step 3. Neighboring tetrahedrons show up as consecutive columns with
//  identical first three entries.  Whenever you spot this happening,
//  make the appropriate entries in TETRA_NEIGHBOR.
//
  for ( j = 0; j < tetra_num; j++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      tetra_neighbor[i+j*4] = -1;
    }
  }

  face = 0;

  for ( ; ; )
  {
    if ( 4 * tetra_num - 1 <= face )
    {
      break;
    }

    if ( faces[0+face*5] == faces[0+(face+1)*5] &&
         faces[1+face*5] == faces[1+(face+1)*5] &&
         faces[2+face*5] == faces[2+(face+1)*5] )
    {
      face1 = faces[3+face*5];
      tetra1 = faces[4+face*5];
      face2 = faces[3+(face+1)*5];
      tetra2 = faces[4+(face+1)*5];
      tetra_neighbor[face1+tetra1*4] = tetra2;
      tetra_neighbor[face2+tetra2*4] = tetra1;
      face = face + 2;
    }
    else
    {
      face = face + 1;
    }
  }

  delete [] faces;

  return tetra_neighbor;
}
//****************************************************************************80

int *tet_mesh_node_order ( int tetra_order, int tetra_num, int tetra_node[], 
  int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_NODE_ORDER: determines the order of nodes.
//
//  Discussion:
//
//    The order of a node is the number of tetrahedrons that use that node
//    as a vertex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TETRA_ORDER, the order of the tetrahedrons.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, int TET_MESH_NODE_ORDER[NODE_NUM], the order of each node.
//
{
  int i;
  int node;
  int *node_order;
  int tetra;

  node_order = new int[node_num];

  i4vec_zero ( node_num, node_order );

  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    for ( i = 0; i < tetra_order; i++ )
    {
      node = tetra_node[i+tetra*tetra_order];
      if ( node < 0 || node_num <= node )
      {
        cout << "\n";
        cout << "TET_MESH_NODE_ORDER - Fatal error!\n";
        cout << "  Illegal entry in TETRA_NODE.\n";
        exit ( 1 );
      }
      else
      {
        node_order[node] = node_order[node] + 1;
      }
    }
  }

  return node_order;
}
//****************************************************************************80

void tet_mesh_order4_adj_count ( int node_num, int tetra_num, 
  int tetra_node[], int *adj_num, int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_ADJ_COUNT counts the number of nodal adjacencies.
//
//  Discussion:
//
//    Assuming that the tet mesh is to be used in a finite element
//    computation, we declare that two distinct nodes are "adjacent" if and
//    only if they are both included in some tetrahedron.
//
//    It is the purpose of this routine to determine the number of
//    such adjacency relationships.
//
//    The initial count gets only the (I,J) relationships, for which
//    node I is strictly less than node J.  This value is doubled
//    to account for symmetry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
//
//    Output, int *ADJ_NUM, the total number of adjacency relationships,
//
//    Output, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
//
{
  int i;
  int j;
  int k;
  int node;
  int *pair;
  int pair_num;
  int pair_unique_num;
  int tetra;
//
//  Each order 4 tetrahedron defines 6 adjacency pairs.
//
  pair = new int[2*6*tetra_num];

  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    pair[0+             tetra *2] = tetra_node[0+tetra*4];
    pair[1+             tetra *2] = tetra_node[1+tetra*4];

    pair[0+(  tetra_num+tetra)*2] = tetra_node[0+tetra*4];
    pair[1+(  tetra_num+tetra)*2] = tetra_node[2+tetra*4];

    pair[0+(2*tetra_num+tetra)*2] = tetra_node[0+tetra*4];
    pair[1+(2*tetra_num+tetra)*2] = tetra_node[3+tetra*4];

    pair[0+(3*tetra_num+tetra)*2] = tetra_node[1+tetra*4];
    pair[1+(3*tetra_num+tetra)*2] = tetra_node[2+tetra*4];

    pair[0+(4*tetra_num+tetra)*2] = tetra_node[1+tetra*4];
    pair[1+(4*tetra_num+tetra)*2] = tetra_node[3+tetra*4];

    pair[0+(5*tetra_num+tetra)*2] = tetra_node[2+tetra*4];
    pair[1+(5*tetra_num+tetra)*2] = tetra_node[3+tetra*4];
  }
  pair_num = 6 * tetra_num;
//
//  Force the nodes of each pair to be listed in ascending order.
//
  i4mat_transpose_print_some ( 2, pair_num, pair, 1, 1, 2, pair_num,
    "DEBUG: PAIR before first sort" );

  i4col_sort2_a ( 2, pair_num, pair );

  i4mat_transpose_print_some ( 2, pair_num, pair, 1, 1, 2, pair_num,
    "DEBUG: PAIR after first sort" );
//
//  Rearrange the columns in ascending order.
//
  i4col_sort_a ( 2, pair_num, pair );
//
//  Get the number of unique columns.
//
  pair_unique_num = i4col_sorted_unique_count ( 2, pair_num, pair );
//
//  The number of adjacencies is TWICE this value, plus the number of nodes.
//
  *adj_num = 2 * pair_unique_num;
//
//  Now set up the ADJ_ROW counts.
//
  for ( node = 0; node < node_num; node++ )
  {
    adj_row[node] = 0;
  }

  for ( k = 0; k < pair_num; k++ )
  {
    if ( 0 < k )
    {
      if ( pair[0+(k-1)*2] == pair[0+k*2] &&
           pair[1+(k-1)*2] == pair[1+k*2] )
      {
        continue;
      }
    }
    i = pair[0+k*2];
    j = pair[1+k*2];

    adj_row[i-1] = adj_row[i-1] + 1;
    adj_row[j-1] = adj_row[j-1] + 1;
  }
//
//  We used ADJ_ROW to count the number of entries in each row.
//  Convert it to pointers into the ADJ array.
//
  for ( node = node_num-1; 0 <= node; node-- )
  {
    adj_row[node] = adj_row[node+1];
  }

  adj_row[0] = 1;
  for ( node = 1; node <= node_num; node++ )
  {
    adj_row[node] = adj_row[node-1] + adj_row[i];
  }

  delete [] pair;

  return;
}
//****************************************************************************80

int *tet_mesh_order4_adj_set ( int node_num, int element_num, 
  int element_node[], int adj_num, int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_ADJ_SET sets the nodal adjacency matrix.
//
//  Discussion:
//
//    A compressed format is used for the nodal adjacency matrix.
//
//    It is assumed that we know ADJ_NUM, the number of adjacency entries
//    and the ADJ_ROW array, which keeps track of the list of slots
//    in ADJ where we can store adjacency information for each row.
//
//    We essentially repeat the work of TET_MESH_ORDER4_ADJ_COUNT, but
//    now we have a place to store the adjacency information.
//
//    A copy of the ADJ_ROW array is useful, as we can use it to keep track
//    of the next available entry in ADJ for adjacencies associated with
//    a given row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
//
//    Input, int ADJ_NUM, the total number of adjacency relationships,
//
//    Input, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
//
//    Output, int TET_MESH_ORDER4_ADJ_SET[ADJ_NUM], 
//    the adjacency information.
//
{
  int *adj;
  int *adj_row_copy;
  int i;
  int j;
  int k;
  int node;
  int *pair;
  int pair_num;
  int tetra;
//
//  Each order 4 tetrahedron defines 6 adjacency pairs.
//
  pair = new int[2*6*element_num];

  for ( tetra = 0; tetra < element_num; tetra++ )
  {
    pair[0+             tetra *2] = element_node[0+tetra*4];
    pair[1+             tetra *2] = element_node[1+tetra*4];

    pair[0+(  element_num+tetra)*2] = element_node[0+tetra*4];
    pair[1+(  element_num+tetra)*2] = element_node[2+tetra*4];

    pair[0+(2*element_num+tetra)*2] = element_node[0+tetra*4];
    pair[1+(2*element_num+tetra)*2] = element_node[3+tetra*4];

    pair[0+(3*element_num+tetra)*2] = element_node[1+tetra*4];
    pair[1+(3*element_num+tetra)*2] = element_node[2+tetra*4];

    pair[0+(4*element_num+tetra)*2] = element_node[1+tetra*4];
    pair[1+(4*element_num+tetra)*2] = element_node[3+tetra*4];

    pair[0+(5*element_num+tetra)*2] = element_node[2+tetra*4];
    pair[1+(5*element_num+tetra)*2] = element_node[3+tetra*4];
  }
  pair_num = 6 * element_num;
//
//  Force the nodes of each pair to be listed in ascending order.
//
  i4col_sort2_a ( 2, pair_num, pair );
//
//  Rearrange the columns in ascending order.
//
  i4col_sort_a ( 2, pair_num, pair );
//
//  Mark all entries of ADJ so we will know later if we missed one.
//
  adj = new int[adj_num];

  for ( i = 0; i < adj_num; i++ )
  {
    adj[i] = -1;
  }
//
//  Copy the ADJ_ROW array and use it to keep track of the next
//  free entry for each row.
//
  adj_row_copy = new int[node_num];

  for ( node = 0; node < node_num; node++ )
  {
    adj_row_copy[node] = adj_row[node];
  }
//
//  Now set up the ADJ_ROW counts.
//
  for ( k = 0; k < pair_num; k++ )
  {
    if ( 0 < k )
    {
      if ( pair[0+(k-1)*2] == pair[0+k*2] &&
           pair[1+(k-1)*2] == pair[1+k*2] )
      {
        continue;
      }
    }
    i = pair[0+k*2];
    j = pair[1+k*2];

    adj[adj_row_copy[i]] = j;
    adj_row_copy[i] = adj_row_copy[i] + 1;
    adj[adj_row_copy[j]] = i;
    adj_row_copy[j] = adj_row_copy[j] + 1;
  }
  delete [] adj_row_copy;
  delete [] pair;

  return adj;
}
//****************************************************************************80

int tet_mesh_order4_boundary_face_count ( int tetra_num, int tetra_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_BOUNDARY_FACE_COUNT counts the number of boundary faces.
//
//  Discussion:
//
//    This routine is given a tet mesh, an abstract list of 
//    quadruples of nodes.  It is assumed that the nodes forming each 
//    face of each tetrahedron are listed in a counterclockwise order, 
//    although the routine should work if the nodes are consistently 
//    listed in a clockwise order as well.
//
//    It is assumed that each face of the tet mesh is either
//    * an INTERIOR face, which is listed twice, once with positive
//      orientation and once with negative orientation, or;
//    * a BOUNDARY face, which will occur only once.
//
//    This routine should work even if the region has holes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
//
//    Output, int TET_MESH_ORDER4_BOUNDARY_FACE_COUNT, the number of 
//    boundary faces.
//
{
  int boundary_face_num;
  int *face;
  int face_num;
  int interior_face_num;
  int m;
  int tet;
  int unique_face_num;

  face = new int[3*4*tetra_num];

  m = 3;
  face_num = 4 * tetra_num;
//
//  Set up the face array:
//  (Omit node 1)
//  (Omit node 2)
//  (Omit node 3)
//  (Omit node 4)
//
  for ( tet = 0; tet < tetra_num; tet++ )
  {
    face[0+(            tet)*3] = tetra_node[1+tet*4];
    face[1+(            tet)*3] = tetra_node[2+tet*4];
    face[2+(            tet)*3] = tetra_node[3+tet*4];

    face[0+(  tetra_num+tet)*3] = tetra_node[0+tet*4];
    face[1+(  tetra_num+tet)*3] = tetra_node[2+tet*4];
    face[2+(  tetra_num+tet)*3] = tetra_node[3+tet*4];

    face[0+(2*tetra_num+tet)*3] = tetra_node[0+tet*4];
    face[1+(2*tetra_num+tet)*3] = tetra_node[1+tet*4];
    face[2+(2*tetra_num+tet)*3] = tetra_node[3+tet*4];
    
    face[0+(3*tetra_num+tet)*3] = tetra_node[0+tet*4];
    face[1+(3*tetra_num+tet)*3] = tetra_node[1+tet*4];
    face[2+(3*tetra_num+tet)*3] = tetra_node[2+tet*4];
  }
//
//  Force the nodes of each face to be listed in ascending order.
//
  i4col_sort2_a ( m, face_num, face );
//
//  Ascending sort the columns.
//
  i4col_sort_a ( m, face_num, face );
//
//  Get the number of unique columns.
//
  unique_face_num = i4col_sorted_unique_count ( m, face_num, face );
//
//  Determine the number of interior and boundary faces.
//
  interior_face_num = 4 * tetra_num - unique_face_num;

  boundary_face_num = 4 * tetra_num - 2 * interior_face_num;

  delete [] face;

  return boundary_face_num;
}
//****************************************************************************80

int tet_mesh_order4_edge_count ( int tetra_num, int tetra_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_EDGE_COUNT counts the number of edges.
//
//  Discussion:
//
//    This routine is given a tet mesh, an abstract list of
//    quadruples of nodes.  Each tetrahedron defines 6 edges; however,
//    assuming that tetrahedrons are touching each other, most edges
//    will be used more than once.  This routine determines the actual
//    number of "geometric" edges associated with the tet mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
//
//    Output, int TET_MESH_ORDER4_EDGE_COUNT, the number of edges.
//
{
  int *edge;
  int edge_num;
  int edge_num_raw;
  int m;
  int tet;

  edge = new int[2*6*tetra_num];

  m = 3;
  edge_num_raw = 6 * tetra_num;
//
//  Set up the raw edge array:
//
  for ( tet = 0; tet < tetra_num; tet++ )
  {
    edge[0+             tet *2] = tetra_node[0+tet*4];
    edge[1+             tet *2] = tetra_node[1+tet*4];

    edge[0+(  tetra_num+tet)*2] = tetra_node[0+tet*4];
    edge[1+(  tetra_num+tet)*2] = tetra_node[2+tet*4];

    edge[0+(2*tetra_num+tet)*2] = tetra_node[0+tet*4];
    edge[1+(2*tetra_num+tet)*2] = tetra_node[3+tet*4];

    edge[0+(3*tetra_num+tet)*2] = tetra_node[1+tet*4];
    edge[1+(3*tetra_num+tet)*2] = tetra_node[2+tet*4];

    edge[0+(4*tetra_num+tet)*2] = tetra_node[1+tet*4];
    edge[1+(4*tetra_num+tet)*2] = tetra_node[3+tet*4];

    edge[0+(5*tetra_num+tet)*2] = tetra_node[2+tet*4];
    edge[1+(5*tetra_num+tet)*2] = tetra_node[3+tet*4];
  }
//
//  Force the nodes of each face to be listed in ascending order.
//
  i4col_sort2_a ( m, edge_num_raw, edge );
//
//  Ascending sort the columns.
//
  i4col_sort_a ( m, edge_num_raw, edge );
//
//  Get the number of unique columns.
//
  edge_num = i4col_sorted_unique_count ( m, edge_num_raw, edge );

  delete [] edge;

  return edge_num;
}
//****************************************************************************80

void tet_mesh_order4_example_set ( int node_num, int tetra_num, 
  double node_xyz[], int tetra_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_EXAMPLE_SET sets an example linear tet mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Output, double NODE_XYZ[3*NODE_NUM], the node coordinates.
//
//    Output, int TETRA_NODE[4*TETRA_NUM], the nodes forming each tet.
//
{
  int i;
  int j;
  double node_xyz_save[3*63] = {
  0.0,  0.0,  0.0, 
  0.0,  0.0,  0.5, 
  0.0,  0.0,  1.0, 
  0.0,  0.5,  0.0, 
  0.0,  0.5,  0.5, 
  0.0,  0.5,  1.0, 
  0.0,  1.0,  0.0, 
  0.0,  1.0,  0.5, 
  0.0,  1.0,  1.0, 
  0.5,  0.0,  0.0, 
  0.5,  0.0,  0.5, 
  0.5,  0.0,  1.0, 
  0.5,  0.5,  0.0, 
  0.5,  0.5,  0.5, 
  0.5,  0.5,  1.0, 
  0.5,  1.0,  0.0, 
  0.5,  1.0,  0.5, 
  0.5,  1.0,  1.0, 
  1.0,  0.0,  0.0, 
  1.0,  0.0,  0.5, 
  1.0,  0.0,  1.0, 
  1.0,  0.5,  0.0, 
  1.0,  0.5,  0.5, 
  1.0,  0.5,  1.0, 
  1.0,  1.0,  0.0, 
  1.0,  1.0,  0.5, 
  1.0,  1.0,  1.0, 
  1.5,  0.0,  0.0, 
  1.5,  0.0,  0.5, 
  1.5,  0.0,  1.0, 
  1.5,  0.5,  0.0, 
  1.5,  0.5,  0.5, 
  1.5,  0.5,  1.0, 
  1.5,  1.0,  0.0, 
  1.5,  1.0,  0.5, 
  1.5,  1.0,  1.0, 
  2.0,  0.0,  0.0, 
  2.0,  0.0,  0.5, 
  2.0,  0.0,  1.0, 
  2.0,  0.5,  0.0, 
  2.0,  0.5,  0.5, 
  2.0,  0.5,  1.0, 
  2.0,  1.0,  0.0, 
  2.0,  1.0,  0.5, 
  2.0,  1.0,  1.0, 
  2.5,  0.0,  0.0, 
  2.5,  0.0,  0.5, 
  2.5,  0.0,  1.0, 
  2.5,  0.5,  0.0, 
  2.5,  0.5,  0.5, 
  2.5,  0.5,  1.0, 
  2.5,  1.0,  0.0, 
  2.5,  1.0,  0.5, 
  2.5,  1.0,  1.0, 
  3.0,  0.0,  0.0, 
  3.0,  0.0,  0.5, 
  3.0,  0.0,  1.0, 
  3.0,  0.5,  0.0, 
  3.0,  0.5,  0.5, 
  3.0,  0.5,  1.0, 
  3.0,  1.0,  0.0, 
  3.0,  1.0,  0.5, 
  3.0,  1.0,  1.0 };
  int tetra_node_save[4*144] = {
     1,   2,   4,  10, 
     2,   4,   5,  10, 
     2,   5,  10,  11, 
     2,   3,   5,  11, 
     4,   5,  10,  13, 
     3,   5,   6,  11, 
     5,  10,  11,  13, 
     4,   5,   7,  13, 
     5,   6,   8,  14, 
     5,   7,   8,  13, 
     6,   8,   9,  14, 
    11,  13,  14,  19, 
    12,  14,  15,  20, 
     3,   6,  11,  12, 
     5,   6,  11,  14, 
     6,   9,  14,  15, 
     6,  11,  12,  14, 
     6,  12,  14,  15, 
     7,   8,  13,  16, 
     5,   8,  13,  14, 
    10,  11,  13,  19, 
     8,   9,  14,  17, 
    11,  12,  14,  20, 
     5,  11,  13,  14, 
     8,  13,  14,  16, 
     9,  14,  15,  17, 
    13,  14,  16,  22, 
     8,  14,  16,  17, 
    14,  15,  17,  23, 
    14,  16,  17,  22, 
     9,  15,  17,  18, 
    15,  17,  18,  23, 
    14,  17,  22,  23, 
    13,  14,  19,  22, 
    11,  14,  19,  20, 
    14,  15,  20,  23, 
    15,  20,  21,  23, 
    21,  23,  24,  29, 
    20,  22,  23,  28, 
    14,  19,  20,  22, 
    15,  18,  23,  24, 
    12,  15,  20,  21, 
    15,  21,  23,  24, 
    16,  17,  22,  25, 
    19,  20,  22,  28, 
    17,  18,  23,  26, 
    20,  21,  23,  29, 
    14,  20,  22,  23, 
    17,  22,  23,  25, 
    18,  23,  24,  26, 
    22,  23,  25,  31, 
    17,  23,  25,  26, 
    23,  24,  26,  32, 
    23,  25,  26,  31, 
    18,  24,  26,  27, 
    24,  26,  27,  32, 
    23,  26,  31,  32, 
    22,  23,  28,  31, 
    20,  23,  28,  29, 
    23,  24,  29,  32, 
    24,  29,  30,  32, 
    30,  32,  33,  38, 
    29,  31,  32,  37, 
    23,  28,  29,  31, 
    24,  27,  32,  33, 
    21,  24,  29,  30, 
    24,  30,  32,  33, 
    25,  26,  31,  34, 
    28,  29,  31,  37, 
    26,  27,  32,  35, 
    29,  30,  32,  38, 
    23,  29,  31,  32, 
    26,  31,  32,  34, 
    27,  32,  33,  35, 
    31,  32,  34,  40, 
    26,  32,  34,  35, 
    32,  33,  35,  41, 
    32,  34,  35,  40, 
    27,  33,  35,  36, 
    33,  35,  36,  41, 
    32,  35,  40,  41, 
    31,  32,  37,  40, 
    29,  32,  37,  38, 
    32,  33,  38,  41, 
    33,  38,  39,  41, 
    39,  41,  42,  47, 
    38,  40,  41,  46, 
    32,  37,  38,  40, 
    33,  36,  41,  42, 
    30,  33,  38,  39, 
    33,  39,  41,  42, 
    34,  35,  40,  43, 
    37,  38,  40,  46, 
    35,  36,  41,  44, 
    38,  39,  41,  47, 
    32,  38,  40,  41, 
    35,  40,  41,  43, 
    36,  41,  42,  44, 
    40,  41,  43,  49, 
    35,  41,  43,  44, 
    41,  42,  44,  50, 
    41,  43,  44,  49, 
    36,  42,  44,  45, 
    42,  44,  45,  50, 
    41,  44,  49,  50, 
    40,  41,  46,  49, 
    38,  41,  46,  47, 
    41,  42,  47,  50, 
    42,  47,  48,  50, 
    48,  50,  51,  56, 
    47,  49,  50,  55, 
    41,  46,  47,  49, 
    42,  45,  50,  51, 
    39,  42,  47,  48, 
    42,  48,  50,  51, 
    43,  44,  49,  52, 
    46,  47,  49,  55, 
    44,  45,  50,  53, 
    47,  48,  50,  56, 
    41,  47,  49,  50, 
    44,  49,  50,  52, 
    45,  50,  51,  53, 
    49,  50,  52,  58, 
    44,  50,  52,  53, 
    50,  51,  53,  59, 
    50,  52,  53,  58, 
    45,  51,  53,  54, 
    51,  53,  54,  59, 
    50,  53,  58,  59, 
    49,  50,  55,  58, 
    47,  50,  55,  56, 
    50,  51,  56,  59, 
    51,  56,  57,  59, 
    50,  55,  56,  58, 
    51,  54,  59,  60, 
    48,  51,  56,  57, 
    51,  57,  59,  60, 
    52,  53,  58,  61, 
    53,  54,  59,  62, 
    50,  56,  58,  59, 
    53,  58,  59,  61, 
    54,  59,  60,  62, 
    53,  59,  61,  62, 
    54,  60,  62,  63 };

  for ( j = 0; j < node_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      node_xyz[i+j*3] = node_xyz_save[i+j*3];
    }
  }

  for ( j = 0; j < tetra_num; j++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      tetra_node[i+j*4] = tetra_node_save[i+j*4] - 1;
    }
  }

  return;
}
//****************************************************************************80

void tet_mesh_order4_example_size ( int *node_num, int *tetra_num )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_EXAMPLE_SIZE sizes an example linear tet mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *TETRA_NUM, the number of tetrahedrons.
//
{
  *node_num = 63;
  *tetra_num = 144;

  return;
}
//****************************************************************************80

void tet_mesh_order4_refine_compute ( int node_num1, int tetra_num1, 
  double node_xyz1[], int tetra_node1[], int node_num2, int tetra_num2,
   int edge_data[], double node_xyz2[], int tetra_node2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_REFINE_COMPUTE computes a refined order 4 tet mesh
//
//  Discussion:
//
//    A refined 4-node tet mesh can be derived from a given
//    4-node tet mesh by interpolating nodes at the midpoint of
//    every edge of the mesh.
//
//    The mesh is described indirectly, as the sum of individual
//    tetrahedrons.  A single physical edge may be a logical edge of
//    any number of tetrahedrons.  It is important, however, that a
//    new node be created exactly once for each edge, assigned an index,
//    and associated with every tetrahedron that shares this edge. 
//
//    This routine handles that problem.
//
//    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
//    data items, one item for every edge of every tetrahedron.  Each
//    data item records, for a given tetrahedron edge, the global indices
//    of the two endpoints, the local indices of the two endpoints,
//    and the index of the tetrahedron.
//
//    Through careful sorting, it is possible to arrange this data in
//    a way that allows the proper generation of the interpolated nodes.
//
//    Let us add the new nodes and temporarily assign them local indices
//    5 through X, based on the following ordering:
//
//      1, 2, 3, 4, (1+2), (1+3), (1+4), (2+3), (2+4), (3+4).
//
//    Then let us assign these nodes to eight subtetrahedrons as follows:
//
//      1, 5, 6, 7
//      2, 5, 8, 9
//      3, 6, 8, 9
//      4, 7, 9, X
//      5, 6, 7, 9
//      5, 6, 8, 9
//      6, 7, 9, X
//      6, 8, 9, X
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Anwei Liu, Barry Joe,
//    Quality Local Refinement of Tetrahedral Meshes Based
//    on 8-Subtetrahedron Subdivision,
//    Mathematics of Computation,
//    Volume 65, Number 215, July 1996, pages 1183-1200.
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes in the input mesh.
//
//    Input, int TETRA_NUM1, the number of tetrahedrons in the
//    input mesh.
//
//    Input, double NODE_XYZ1[3*NODE_NUM1], the coordinates of
//    the nodes that make up the input mesh.
//
//    Input, int TETRA_NODE1[4*TETRA_NUM], the indices of the nodes
//    in the input mesh.
//
//    Input, int NODE_NUM2, the number of nodes in the refined mesh.
//
//    Input, int TETRA_NUM2, the number of tetrahedrons in the
//    refined mesh.
//
//    Input, int EDGE_DATA[5*(6*TETRA_NUM1)], edge data.
//
//    Output, double NODE_XYZ2[3*NODE_NUM2], the coordinates of
//    the nodes that make up the refined mesh.
//
//    Output, int TETRA_NODE2[4*TETRA_NUM2], the indices of the nodes 
//    in the refined mesh.
//
{
  int dim_num = 3;
  int edge;
  int i;
  int j;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int node;
  int tetra_order = 4;
  int tetra1;
  int tetra2;
  int v;
  int v1;
  int v2;
//
//  Generate the index and coordinates of the new midside nodes, 
//  and update the tetradehron-node data.
//
  for ( j = 0; j < node_num1; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      node_xyz2[i+j*dim_num] = node_xyz1[i+j*dim_num];
    }
  }
  for ( j = 0; j < tetra_num2; j++ )
  {
    for ( i = 0; i < tetra_order; i++ )
    {
      tetra_node2[i+j*tetra_order] = -1;
    }
  }
//
//  The vertices of the input tetrahedron can be assigned now.
//
  for ( tetra1 = 0; tetra1 < tetra_num1; tetra1++ )
  {
    tetra_node2[0+(tetra1*8+0)*tetra_order] = tetra_node1[0+tetra1*tetra_order];
    tetra_node2[0+(tetra1*8+1)*tetra_order] = tetra_node1[1+tetra1*tetra_order];
    tetra_node2[0+(tetra1*8+2)*tetra_order] = tetra_node1[2+tetra1*tetra_order];
    tetra_node2[0+(tetra1*8+3)*tetra_order] = tetra_node1[3+tetra1*tetra_order];
  }
  node = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 6 * tetra_num1; edge++ )
  {
//
//  Read the data defining the edge.
//
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];
//
//  If this edge is new, create the coordinates and index.
//
    if ( n1 != n1_old || n2 != n2_old )
    {
      if ( node_num2 <= node )
      {
        cout << "\n";
        cout << "TET_MESH_ORDER4_REFINE_COMPUTE - Fatal error!\n";
        cout << "  Node index exceeds NODE_NUM2.\n";
        exit ( 1 );
      }

      for ( i = 0; i < dim_num; i++ )
      {
        node_xyz2[i+node*dim_num] = 
        ( node_xyz2[i+(n1-1)*dim_num] + node_xyz2[i+(n2-1)*dim_num] ) / 2.0;
      }
      node = node + 1;
      n1_old = n1;
      n2_old = n2;
    }
//
//  Assign the node to the tetrahedron.
//
    v1 = edge_data[2+edge*5];
    v2 = edge_data[3+edge*5];
    tetra1 = edge_data[4+edge*5];
//
//  We know the two vertices that bracket this new node.
//  This tells us whether it is new node number 5, 6, 7, 8, 9 or 10.
//  This tells us which of the new subtetrahedrons it belongs to,
//  and what position it occupies.
//
    if ( v1 == 1 && v2 == 2 )
    {
      tetra_node2[1+(tetra1*8+0)*tetra_order] = node;
      tetra_node2[1+(tetra1*8+1)*tetra_order] = node;
      tetra_node2[0+(tetra1*8+4)*tetra_order] = node;
      tetra_node2[0+(tetra1*8+5)*tetra_order] = node;
    }
    else if ( v1 == 1 && v2 == 3 )
    {
      tetra_node2[2+(tetra1*8+0)*tetra_order] = node;
      tetra_node2[1+(tetra1*8+2)*tetra_order] = node;
      tetra_node2[1+(tetra1*8+4)*tetra_order] = node;
      tetra_node2[1+(tetra1*8+5)*tetra_order] = node;
      tetra_node2[0+(tetra1*8+6)*tetra_order] = node;
      tetra_node2[0+(tetra1*8+7)*tetra_order] = node;
    }
    else if ( v1 == 1 && v2 == 4 )
    {
      tetra_node2[3+(tetra1*8+0)*tetra_order] = node;
      tetra_node2[1+(tetra1*8+3)*tetra_order] = node;
      tetra_node2[2+(tetra1*8+4)*tetra_order] = node;
      tetra_node2[1+(tetra1*8+6)*tetra_order] = node;
    }
    else if ( v1 == 2 && v2 == 3 )
    {
      tetra_node2[2+(tetra1*8+1)*tetra_order] = node;
      tetra_node2[2+(tetra1*8+2)*tetra_order] = node;
      tetra_node2[2+(tetra1*8+5)*tetra_order] = node;
      tetra_node2[1+(tetra1*8+7)*tetra_order] = node;
    }
    else if ( v1 == 2 && v2 == 4 )
    {
      tetra_node2[3+(tetra1*8+1)*tetra_order] = node;
      tetra_node2[3+(tetra1*8+2)*tetra_order] = node;
      tetra_node2[2+(tetra1*8+3)*tetra_order] = node;
      tetra_node2[3+(tetra1*8+4)*tetra_order] = node;
      tetra_node2[3+(tetra1*8+5)*tetra_order] = node;
      tetra_node2[2+(tetra1*8+6)*tetra_order] = node;
      tetra_node2[2+(tetra1*8+7)*tetra_order] = node;
    }
    else if ( v1 == 3 && v2 == 4 )
    {
      tetra_node2[3+(tetra1*8+3)*tetra_order] = node;
      tetra_node2[3+(tetra1*8+6)*tetra_order] = node;
      tetra_node2[3+(tetra1*8+7)*tetra_order] = node;
    }
  }

  return;
}
//****************************************************************************80

void tet_mesh_order4_refine_size ( int node_num1, int tetra_num1, 
  int tetra_node1[], int *node_num2, int *tetra_num2, int edge_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_REFINE_SIZE sizes a refined order 4 tet mesh.
//
//  Discussion:
//
//    A refined tet mesh can be derived from an existing one by interpolating 
//    nodes at the midpoint of every edge of the mesh.
//
//    The mesh is described indirectly, as the sum of individual
//    tetrahedrons.  A single physical edge may be a logical edge of
//    any number of tetrahedrons.  It is important, however, that a
//    new node be created exactly once for each edge, assigned an index,
//    and associated with every tetrahedron that shares this edge. 
//
//    This routine handles that problem.
//
//    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
//    data items, one item for every edge of every tetrahedron.  Each
//    data item records, for a given tetrahedron edge, the global indices
//    of the two endpoints, the local indices of the two endpoints,
//    and the index of the tetrahedron.
//
//    Through careful sorting, it is possible to arrange this data in
//    a way that allows the proper generation of the interpolated nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes in the original mesh.
//
//    Input, int TETRA_NUM1, the number of tetrahedrons in the
//    original mesh.
//
//    Input, int TETRA_NODE1[4*TETRA_NUM1], the indices of the nodes
//    in the original mesh.
//
//    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
//
//    Output, int *TETRA_NUM2, the number of tetrahedrons in the refined mesh.
//
//    Output, int EDGE_DATA[5*(6*TETRA_NUM1)], edge data.
//
{
  int a;
  int b;
  int edge;
  int i;
  int j;
  int k;
  int l;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int tetra;
  int tetra_order = 4;
//
//  Step 1.
//  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
//  construct the six edge relations:
//
//    (I,J,1,2,T)
//    (I,K,1,3,T)
//    (I,L,1,4,T)
//    (J,K,2,3,T)
//    (J,L,2,4,T)
//    (K,L,3,4,T)
//
//  In order to make matching easier, we reorder each pair of nodes
//  into ascending order.
//
  for ( tetra = 0; tetra < tetra_num1; tetra++ )
  {
    i = tetra_node1[0+tetra*tetra_order];
    j = tetra_node1[1+tetra*tetra_order];
    k = tetra_node1[2+tetra*tetra_order];
    l = tetra_node1[3+tetra*tetra_order];

    i4i4_sort_a ( i, j, &a, &b );

    edge_data[0+(6*tetra)*5] = a;
    edge_data[1+(6*tetra)*5] = b;
    edge_data[2+(6*tetra)*5] = 1;
    edge_data[3+(6*tetra)*5] = 2;
    edge_data[4+(6*tetra)*5] = tetra;

    i4i4_sort_a ( i, k, &a, &b );

    edge_data[0+(6*tetra+1)*5] = a;
    edge_data[1+(6*tetra+1)*5] = b;
    edge_data[2+(6*tetra+1)*5] = 1;
    edge_data[3+(6*tetra+1)*5] = 3;
    edge_data[4+(6*tetra+1)*5] = tetra;

    i4i4_sort_a ( i, l, &a, &b );

    edge_data[0+(6*tetra+2)*5] = a;
    edge_data[1+(6*tetra+2)*5] = b;
    edge_data[2+(6*tetra+2)*5] = 1;
    edge_data[3+(6*tetra+2)*5] = 4;
    edge_data[4+(6*tetra+2)*5] = tetra;

    i4i4_sort_a ( j, k, &a, &b );

    edge_data[0+(6*tetra+3)*5] = a;
    edge_data[1+(6*tetra+3)*5] = b;
    edge_data[2+(6*tetra+3)*5] = 2;
    edge_data[3+(6*tetra+3)*5] = 3;
    edge_data[4+(6*tetra+3)*5] = tetra;

    i4i4_sort_a ( j, l, &a, &b );

    edge_data[0+(6*tetra+4)*5] = a;
    edge_data[1+(6*tetra+4)*5] = b;
    edge_data[2+(6*tetra+4)*5] = 2;
    edge_data[3+(6*tetra+4)*5] = 4;
    edge_data[4+(6*tetra+4)*5] = tetra;

    i4i4_sort_a ( k, l, &a, &b );

    edge_data[0+(6*tetra+5)*5] = a;
    edge_data[1+(6*tetra+5)*5] = b;
    edge_data[2+(6*tetra+5)*5] = 3;
    edge_data[3+(6*tetra+5)*5] = 4;
    edge_data[4+(6*tetra+5)*5] = tetra;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:2; the routine we call here
//  sorts on the full column but that won't hurt us.
//
//  What we need is to find all cases where tetrahedrons share an edge.
//  By sorting the columns of the EDGE_DATA array, we will put shared edges
//  next to each other.
//
  i4col_sort_a ( 5, 6*tetra_num1, edge_data );
//
//  Step 3. All the tetrahedrons which share an edge show up as consecutive
//  columns with identical first two entries.  Figure out how many new
//  nodes there are, and allocate space for their coordinates.
//
  *node_num2 = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 6 * tetra_num1; edge++ )
  {
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];
    if ( n1 != n1_old || n2 != n2_old )
    {
      *node_num2 = *node_num2 + 1;
      n1_old = n1;
      n2_old = n2;
    }
  }

  *tetra_num2 = 8 * tetra_num1;

  return;
}
//****************************************************************************80

void tet_mesh_order4_to_order10_compute ( int tetra_num, int tetra_node1[], 
  int node_num1, double node_xyz1[], int edge_data[], int tetra_node2[], 
  int node_num2, double node_xyz2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_TO_ORDER10_COMPUTE computes a quadratic tet mesh from a linear one.
//
//  Discussion:
//
//    A quadratic (10 node) tet mesh can be derived from a linear
//    (4 node) tet mesh by interpolating nodes at the midpoint of
//    every edge of the mesh.
//
//    The mesh is described indirectly, as the sum of individual
//    tetrahedrons.  A single physical edge may be a logical edge of
//    any number of tetrahedrons.  It is important, however, that a
//    new node be created exactly once for each edge, assigned an index,
//    and associated with every tetrahedron that shares this edge. 
//
//    This routine handles that problem.
//
//    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
//    data items, one item for every edge of every tetrahedron.  Each
//    data item records, for a given tetrahedron edge, the global indices
//    of the two endpoints, the local indices of the two endpoints,
//    and the index of the tetrahedron.
//
//    Through careful sorting, it is possible to arrange this data in
//    a way that allows the proper generation of the interpolated nodes.
//
//    The node ordering for the quadratic tetrahedron is somewhat
//    arbitrary.  In the current scheme, the vertices are listed
//    first, followed by the 6 midside nodes.  Each midside node
//    may be identified by the two vertices that bracket it.  Thus,
//    the node ordering may be suggested by:
//
//      1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TETRA_NUM, the number of tetrahedrons in the
//    linear mesh.
//
//    Input, int TETRA_NODE1[4*TETRA_NUM], the indices of the nodes
//    in the linear mesh.
//
//    Input, int NODE_NUM1, the number of nodes for the linear mesh.
//
//    Input, double NODE_XYZ1[3*NODE_NUM1], the coordinates of
//    the nodes that make up the linear mesh.
//
//    Input, int EDGE_DATA[5*(6*TETRA_NUM)], edge data.
//
//    Output, int TETRA_NODE2[10*TETRA_NUM], the indices of the nodes 
//    in the quadratic mesh.
//
//    Input, int NODE_NUM2, the number of nodes for the quadratic mesh.
//
//    Output, double NODE_XYZ2[3*NODE_NUM2], the coordinates of
//    the nodes that make up the quadratic mesh.
//
{
  int dim_num = 3;
  int edge;
  int i;
  int j;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int node;
  int tetra;
  int tetra_order1 = 4;
  int tetra_order2 = 10;
  int v;
  int v1;
  int v2;
//
//  Generate the index and coordinates of the new midside nodes, 
//  and update the tetradehron-node data.
//
  for ( j = 0; j < node_num1; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      node_xyz2[i+j*dim_num] = node_xyz1[i+j*dim_num];
    }
  }
  for ( j = 0; j < tetra_num; j++ )
  {
    for ( i = 0; i < tetra_order1; i++ )
    {
      tetra_node2[i+j*tetra_order2] = tetra_node1[i+j*tetra_order1];
    }
  }
  node = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 6 * tetra_num; edge++ )
  {
//
//  Read the data defining the edge.
//
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];
//
//  If this edge is new, create the coordinates and index.
//
    if ( n1 != n1_old || n2 != n2_old )
    {
      if ( node_num2 <= node )
      {
        cout << "\n";
        cout << "TET_MESH_ORDER4_TO_ORDER10_COMPUTE - Fatal error!\n";
        cout << "  Node index exceeds NODE_NUM2.\n";
        exit ( 1 );
      }

      for ( i = 0; i < dim_num; i++ )
      {
        node_xyz2[i+node*dim_num] = 
        ( node_xyz2[i+(n1-1)*dim_num] + node_xyz2[i+(n2-1)*dim_num] ) / 2.0;
      }
      node = node + 1;
      n1_old = n1;
      n2_old = n2;
    }
//
//  Assign the node to the tetrahedron.
//
    v1 = edge_data[2+edge*5];
    v2 = edge_data[3+edge*5];
//
//  Here is where the local ordering of the nodes is effected:
//
    if ( v1 == 1 && v2 == 2 )
    {
      v = 5;
    }
    else if ( v1 == 1 && v2 == 3 )
    {
      v = 6;
    }
    else if ( v1 == 1 && v2 == 4 )
    {
      v = 7;
    }
    else if ( v1 == 2 && v2 == 3 )
    {
      v = 8;
    }
    else if ( v1 == 2 && v2 == 4 )
    {
      v = 9;
    }
    else if ( v1 == 3 && v2 == 4 )
    {
      v = 10;
    }

    tetra = edge_data[4+edge*5];

    tetra_node2[v-1+tetra*tetra_order2] = node;
  }

  return;
}
//****************************************************************************80

void tet_mesh_order4_to_order10_size ( int tetra_num, int tetra_node1[], 
  int node_num1, int edge_data[], int *node_num2 )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_TO_ORDER10_SIZE sizes a quadratic tet mesh from a linear one.
//
//  Discussion:
//
//    A quadratic (10 node) tet mesh can be derived from a linear
//    (4 node) tet mesh by interpolating nodes at the midpoint of
//    every edge of the mesh.
//
//    The mesh is described indirectly, as the sum of individual
//    tetrahedrons.  A single physical edge may be a logical edge of
//    any number of tetrahedrons.  It is important, however, that a
//    new node be created exactly once for each edge, assigned an index,
//    and associated with every tetrahedron that shares this edge. 
//
//    This routine handles that problem.
//
//    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
//    data items, one item for every edge of every tetrahedron.  Each
//    data item records, for a given tetrahedron edge, the global indices
//    of the two endpoints, the local indices of the two endpoints,
//    and the index of the tetrahedron.
//
//    Through careful sorting, it is possible to arrange this data in
//    a way that allows the proper generation of the interpolated nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TETRA_NUM, the number of tetrahedrons in the
//    linear mesh.
//
//    Input, int TETRA_NODE1[4*TETRA_NUM], the indices of the nodes
//    in the linear mesh.
//
//    Input, int NODE_NUM1, the number of nodes for the linear mesh.
//
//    Output, int EDGE_DATA[5*(6*TETRA_NUM)], edge data.
//
//    Output, int *NODE_NUM2, the number of nodes for the quadratic mesh.
//
{
  int a;
  int b;
  int edge;
  int i;
  int j;
  int k;
  int l;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int tetra;
  int tetra_order1 = 4;
//
//  Step 1.
//  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
//  construct the six edge relations:
//
//    (I,J,1,2,T)
//    (I,K,1,3,T)
//    (I,L,1,4,T)
//    (J,K,2,3,T)
//    (J,L,2,4,T)
//    (K,L,3,4,T)
//
//  In order to make matching easier, we reorder each pair of nodes
//  into ascending order.
//
  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    i = tetra_node1[0+tetra*tetra_order1];
    j = tetra_node1[1+tetra*tetra_order1];
    k = tetra_node1[2+tetra*tetra_order1];
    l = tetra_node1[3+tetra*tetra_order1];

    i4i4_sort_a ( i, j, &a, &b );

    edge_data[0+(6*tetra)*5] = a;
    edge_data[1+(6*tetra)*5] = b;
    edge_data[2+(6*tetra)*5] = 1;
    edge_data[3+(6*tetra)*5] = 2;
    edge_data[4+(6*tetra)*5] = tetra;

    i4i4_sort_a ( i, k, &a, &b );

    edge_data[0+(6*tetra+1)*5] = a;
    edge_data[1+(6*tetra+1)*5] = b;
    edge_data[2+(6*tetra+1)*5] = 1;
    edge_data[3+(6*tetra+1)*5] = 3;
    edge_data[4+(6*tetra+1)*5] = tetra;

    i4i4_sort_a ( i, l, &a, &b );

    edge_data[0+(6*tetra+2)*5] = a;
    edge_data[1+(6*tetra+2)*5] = b;
    edge_data[2+(6*tetra+2)*5] = 1;
    edge_data[3+(6*tetra+2)*5] = 4;
    edge_data[4+(6*tetra+2)*5] = tetra;

    i4i4_sort_a ( j, k, &a, &b );

    edge_data[0+(6*tetra+3)*5] = a;
    edge_data[1+(6*tetra+3)*5] = b;
    edge_data[2+(6*tetra+3)*5] = 2;
    edge_data[3+(6*tetra+3)*5] = 3;
    edge_data[4+(6*tetra+3)*5] = tetra;

    i4i4_sort_a ( j, l, &a, &b );

    edge_data[0+(6*tetra+4)*5] = a;
    edge_data[1+(6*tetra+4)*5] = b;
    edge_data[2+(6*tetra+4)*5] = 2;
    edge_data[3+(6*tetra+4)*5] = 4;
    edge_data[4+(6*tetra+4)*5] = tetra;

    i4i4_sort_a ( k, l, &a, &b );

    edge_data[0+(6*tetra+5)*5] = a;
    edge_data[1+(6*tetra+5)*5] = b;
    edge_data[2+(6*tetra+5)*5] = 3;
    edge_data[3+(6*tetra+5)*5] = 4;
    edge_data[4+(6*tetra+5)*5] = tetra;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:2; the routine we call here
//  sorts on the full column but that won't hurt us.
//
//  What we need is to find all cases where tetrahedrons share an edge.
//  By sorting the columns of the EDGE_DATA array, we will put shared edges
//  next to each other.
//
  i4col_sort_a ( 5, 6*tetra_num, edge_data );
//
//  Step 3. All the tetrahedrons which share an edge show up as consecutive
//  columns with identical first two entries.  Figure out how many new
//  nodes there are, and allocate space for their coordinates.
//
  *node_num2 = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 6 * tetra_num; edge++ )
  {
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];
    if ( n1 != n1_old || n2 != n2_old )
    {
      *node_num2 = *node_num2 + 1;
      n1_old = n1;
      n2_old = n2;
    }
  }

  return;
}
//****************************************************************************80

void tet_mesh_order10_adj_count ( int node_num, int tet_num, 
  int tet_node[], int *adj_num, int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER10_ADJ_COUNT counts the number of nodal adjacencies.
//
//  Discussion:
//
//    Assuming that the tet mesh is to be used in a finite element
//    computation, we declare that two distinct nodes are "adjacent" if and
//    only if they are both included in some tetrahedron.
//
//    It is the purpose of this routine to determine the number of
//    such adjacency relationships.
//
//    The initial count gets only the (I,J) relationships, for which
//    node I is strictly less than node J.  This value is doubled
//    to account for symmetry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TET_NUM, the number of tetrahedrons.
//
//    Input, int TET_NODE[10*TET_NUM], the indices of the nodes.
//
//    Output, int *ADJ_NUM, the total number of adjacency relationships,
//
//    Output, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
//
{
  int i;
  int j;
  int k;
  int l;
  int node;
  int *pair;
  int pair_num;
  int pair_unique_num;
//
//  Each order 10 tetrahedron defines 45 adjacency pairs.
//
  pair = new int[2*45*tet_num];

  k = 0;
  for ( i = 0; i < 9; i++ )
  {
    for ( j = i + 1; j < 10; j++ )
    {
      for ( l = 0; l < tet_num; l++ )
      {
        pair[0+(k*tet_num+l)*2] = tet_node[i+l*10];
        pair[1+(k*tet_num+l)*2] = tet_node[j+l*10];
      }
      k = k + 1;
    }
  }
//
//  Force the nodes of each pair to be listed in ascending order.
//
  pair_num = 45 * tet_num;

  i4col_sort2_a ( 2, pair_num, pair );
//
//  Rearrange the columns in ascending order.
//
  i4col_sort_a ( 2, pair_num, pair );
//
//  Get the number of unique columns.
//
  pair_unique_num = i4col_sorted_unique_count ( 2, pair_num, pair );
//
//  The number of adjacencies is TWICE this value, plus the number of nodes.
//
  *adj_num = 2 * pair_unique_num;
//
//  Now set up the ADJ_ROW counts.
//
  for ( node = 0; node < node_num; node++ )
  {
    adj_row[node] = 0;
  }

  for ( k = 0; k < pair_num; k++ )
  {
    if ( 0 < k )
    {
      if ( pair[0+(k-1)*2] == pair[0+k*2] &&
           pair[1+(k-1)*2] == pair[1+k*2] )
      {
        continue;
      }
    }
    i = pair[0+k*2];
    j = pair[1+k*2];

    adj_row[i-1] = adj_row[i-1] + 1;
    adj_row[j-1] = adj_row[j-1] + 1;
  }
//
//  We used ADJ_ROW to count the number of entries in each row.
//  Convert it to pointers into the ADJ array.
//
  for ( node = node_num-1; 0 <= node; node-- )
  {
    adj_row[node] = adj_row[node+1];
  }

  adj_row[0] = 1;
  for ( node = 1; node <= node_num; node++ )
  {
    adj_row[node] = adj_row[node-1] + adj_row[i];
  }

  delete [] pair;

  return;
}
//****************************************************************************80

int *tet_mesh_order10_adj_set ( int node_num, int tet_num, 
  int tet_node[], int adj_num, int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER10_ADJ_SET sets the nodal adjacency matrix.
//
//  Discussion:
//
//    A compressed format is used for the nodal adjacency matrix.
//
//    It is assumed that we know ADJ_NUM, the number of adjacency entries
//    and the ADJ_ROW array, which keeps track of the list of slots
//    in ADJ where we can store adjacency information for each row.
//
//    We essentially repeat the work of TET_MESH_ORDER4_ADJ_COUNT, but
//    now we have a place to store the adjacency information.
//
//    A copy of the ADJ_ROW array is useful, as we can use it to keep track
//    of the next available entry in ADJ for adjacencies associated with
//    a given row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TET_NUM, the number of tetrahedrons.
//
//    Input, int TET_NODE[10*TET_NUM], the indices of the nodes.
//
//    Input, int ADJ_NUM, the total number of adjacency relationships,
//
//    Input, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
//
//    Output, int TET_MESH_ORDER4_ADJ_SET[ADJ_NUM], 
//    the adjacency information.
//
{
  int *adj;
  int *adj_row_copy;
  int i;
  int j;
  int k;
  int l;
  int node;
  int *pair;
  int pair_num;
//
//  Each order 10 tetrahedron defines 45 adjacency pairs.
//
  pair = new int[2*45*tet_num];

  k = 0;
  for ( i = 0; i < 9; i++ )
  {
    for ( j = i + 1; j < 10; j++ )
    {
      for ( l = 0; l < tet_num; l++ )
      {
        pair[0+(k*tet_num+l)*2] = tet_node[i+l*10];
        pair[1+(k*tet_num+l)*2] = tet_node[j+l*10];
      }
      k = k + 1;
    }
  }
//
//  Force the nodes of each pair to be listed in ascending order.
//
  pair_num = 45 * tet_num;

  i4col_sort2_a ( 2, pair_num, pair );
//
//  Rearrange the columns in ascending order.
//
  i4col_sort_a ( 2, pair_num, pair );
//
//  Mark all entries of ADJ so we will know later if we missed one.
//
  adj = new int[adj_num];

  for ( i = 0; i < adj_num; i++ )
  {
    adj[i] = -1;
  }
//
//  Copy the ADJ_ROW array and use it to keep track of the next
//  free entry for each row.
//
  adj_row_copy = new int[node_num];

  for ( node = 0; node < node_num; node++ )
  {
    adj_row_copy[node] = adj_row[node];
  }
//
//  Now set up the ADJ_ROW counts.
//
  for ( k = 0; k < pair_num; k++ )
  {
    if ( 0 < k )
    {
      if ( pair[0+(k-1)*2] == pair[0+k*2] &&
           pair[1+(k-1)*2] == pair[1+k*2] )
      {
        continue;
      }
    }
    i = pair[0+k*2];
    j = pair[1+k*2];

    adj[adj_row_copy[i]] = j;
    adj_row_copy[i] = adj_row_copy[i] + 1;
    adj[adj_row_copy[j]] = i;
    adj_row_copy[j] = adj_row_copy[j] + 1;
  }
  delete [] adj_row_copy;
  delete [] pair;

  return adj;
}
//****************************************************************************80

void tet_mesh_order10_example_set ( int node_num, int tetra_num, 
  double node_xyz[], int tetra_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER10_EXAMPLE_SET sets an example quadratic tet mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Output, double NODE_XYZ[3*NODE_NUM], the node coordinates.
//
//    Output, int TETRA_NODE[10*TETRA_NUM], the nodes forming each tet.
//
{
  int i;
  int j;
  double node_xyz_save[3*27] = {
   0.0,  0.0,  0.0, 
   0.0,  0.0,  1.0, 
   0.0,  1.0,  0.0, 
   0.0,  1.0,  1.0, 
   1.0,  0.0,  0.0, 
   1.0,  0.0,  1.0, 
   1.0,  1.0,  0.0, 
   1.0,  1.0,  1.0, 
   0.0,  0.0,  0.5, 
   0.0,  0.5,  0.0, 
   0.0,  0.5,  0.5, 
   0.5,  0.0,  0.0, 
   0.0,  0.5,  1.0, 
   0.5,  0.0,  0.5, 
   0.5,  0.0,  1.0, 
   0.0,  1.0,  0.5, 
   0.5,  0.5,  0.0, 
   0.5,  1.0,  0.0, 
   0.5,  0.5,  0.5, 
   0.5,  0.5,  1.0, 
   0.5,  1.0,  0.5, 
   0.5,  1.0,  1.0, 
   1.0,  0.0,  0.5, 
   1.0,  0.5,  0.0, 
   1.0,  0.5,  0.5, 
   1.0,  0.5,  1.0, 
   1.0,  1.0,  0.5  };
  int tetra_node_save[10*6] = {
    4,   3,   5,   1,  16,  19,  17,  11,  10,  12, 
    4,   2,   5,   1,  13,  19,  14,  11,   9,  12, 
    4,   7,   3,   5,  21,  16,  18,  19,  24,  17, 
    4,   7,   8,   5,  21,  22,  27,  19,  24,  25, 
    4,   6,   2,   5,  20,  13,  15,  19,  23,  14, 
    4,   6,   8,   5,  20,  22,  26,  19,  23,  25 };

  for ( j = 0; j < node_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      node_xyz[i+j*3] = node_xyz_save[i+j*3];
    }
  }

  for ( j = 0; j < tetra_num; j++ )
  {
    for ( i = 0; i < 10; i++ )
    {
      tetra_node[i+j*10] = tetra_node_save[i+j*10] - 1;
    }
  }

  return;
}
//****************************************************************************80

void tet_mesh_order10_example_size ( int *node_num, int *tetra_num )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER10_EXAMPLE_SIZE sizes an example quadratic tet mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *TETRA_NUM, the number of tetrahedrons.
//
{
  *node_num = 27;
  *tetra_num = 6;

  return;
}
//****************************************************************************80

void tet_mesh_order10_to_order4_compute ( int tetra_num1, int tetra_node1[], 
  int tetra_num2, int tetra_node2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER10_TO_ORDER4_COMPUTE linearizes a quadratic tet mesh.
//
//  Discussion:
//
//    A quadratic tet mesh is assumed to consist of 10-node
//    tetrahedrons.
//
//    This routine rearranges the information so as to define a 4-node
//    tet mesh.
//
//    The same nodes are used, but there are 8 times as many
//    tetrahedrons.
//
//    The node ordering for the quadratic tetrahedron is somewhat
//    arbitrary.  In the current scheme, the vertices are listed
//    first, followed by the 6 midside nodes.  Each midside node
//    may be identified by the two vertices that bracket it.  Thus,
//    the node ordering may be suggested by:
//
//      1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Anwei Liu, Barry Joe,
//    Quality Local Refinement of Tetrahedral Meshes Based
//    on 8-Subtetrahedron Subdivision,
//    Mathematics of Computation,
//    Volume 65, Number 215, July 1996, pages 1183-1200.
//
//  Parameters:
//
//    Input, int TETRA_NUM1, the number of tetrahedrons in the quadratic
//    tet mesh.
//
//    Input, int TETRA_NODE1[10*TETRA_NUM1], the indices of the nodes
//    that made up the quadratic mesh.
//
//    Input, int TETRA_NUM2, the number of tetrahedrons in the linear
//    tet mesh.  TETRA_NUM2 = 8 * TETRA_NUM1.
//
//    Output, int TETRA_NODE2[4*TETRA_NUM2], the indices of the nodes
//    that make up the linear mesh.
//
{
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int n7;
  int n8;
  int n9;
  int nx;
  int tetra1;
  int tetra2;

  tetra2 = 0;

  for ( tetra1 = 0; tetra1 < tetra_num1; tetra1++ )
  {
    n1 = tetra_node1[0+tetra1*10];
    n2 = tetra_node1[1+tetra1*10];
    n3 = tetra_node1[2+tetra1*10];
    n4 = tetra_node1[3+tetra1*10];
    n5 = tetra_node1[4+tetra1*10];
    n6 = tetra_node1[5+tetra1*10];
    n7 = tetra_node1[6+tetra1*10];
    n8 = tetra_node1[7+tetra1*10];
    n9 = tetra_node1[8+tetra1*10];
    nx = tetra_node1[9+tetra1*10];

    tetra_node2[0+tetra2*4] = n1;
    tetra_node2[1+tetra2*4] = n5;
    tetra_node2[2+tetra2*4] = n6;
    tetra_node2[3+tetra2*4] = n7;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n2;
    tetra_node2[1+tetra2*4] = n5;
    tetra_node2[2+tetra2*4] = n8;
    tetra_node2[3+tetra2*4] = n9;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n3;
    tetra_node2[1+tetra2*4] = n6;
    tetra_node2[2+tetra2*4] = n8;
    tetra_node2[3+tetra2*4] = n9;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n4;
    tetra_node2[1+tetra2*4] = n7;
    tetra_node2[2+tetra2*4] = n9;
    tetra_node2[3+tetra2*4] = nx;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n5;
    tetra_node2[1+tetra2*4] = n6;
    tetra_node2[2+tetra2*4] = n7;
    tetra_node2[3+tetra2*4] = n9;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n5;
    tetra_node2[1+tetra2*4] = n6;
    tetra_node2[2+tetra2*4] = n8;
    tetra_node2[3+tetra2*4] = n9;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n6;
    tetra_node2[1+tetra2*4] = n7;
    tetra_node2[2+tetra2*4] = n9;
    tetra_node2[3+tetra2*4] = nx;
    tetra2 = tetra2 + 1;

    tetra_node2[0+tetra2*4] = n6;
    tetra_node2[1+tetra2*4] = n8;
    tetra_node2[2+tetra2*4] = n9;
    tetra_node2[3+tetra2*4] = nx;
    tetra2 = tetra2 + 1;
  }

  return;
}
//****************************************************************************80

void tet_mesh_order10_to_order4_size ( int node_num1, int tetra_num1,
  int *node_num2, int *tetra_num2 )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER10_TO_ORDER4_SIZE sizes a linear tet mesh from a quadratic one.
//
//  Discussion:
//
//    A linear (4 node) tet mesh can be derived from a quadratic
//    (10 node) tet mesh using the same set of nodes, but reassigning
//    the nodes of each quadratic tet among 8 linear subtets.
//
//    This routine returns the number of nodes and tetrahedra in the
//    linear mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Anwei Liu, Barry Joe,
//    Quality Local Refinement of Tetrahedral Meshes Based
//    on 8-Subtetrahedron Subdivision,
//    Mathematics of Computation,
//    Volume 65, Number 215, July 1996, pages 1183-1200.
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes in the quadratic mesh.
//
//    Input, int TETRA_NUM1, the number of tetrahedrons in the
//    quadratic mesh.
//
//    Output, int *NODE_NUM2, the number of nodes for the linear mesh.
//
//    Output, int *TETRA_NUM2, the number of tetrahedrons in the
//    linear mesh.
//
{
  *node_num2 = node_num1;
  *tetra_num2 = 8 * tetra_num1;

  return;
}
//****************************************************************************80

void tet_mesh_quad ( int node_num, double node_xyz[], int tetra_order, 
  int tetra_num, int tetra_node[], 
  void quad_fun ( int n, double xyz_vec[], double fvec[] ), 
  int quad_num, double quad_xyz[], double quad_w[], double *quad_value, 
  double *region_volume )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_QUAD approximates an integral over a tet mesh.
//
//  Discussion:
//
//    The routine will accept tetrahedral meshes of order higher than 4.
//    However, only the first four nodes (the vertices) of each
//    tetrahedron will be used.  This will still produce correct results
//    for higher order tet meshes, as long as the sides of each
//    tetrahedron are flat (linear).
//
//    We assume that the vertices of each tetrahedron are listed first
//    in the description of higher order tetrahedrons.
//
//    The approximation of the integral is made using a quadrature rule 
//    defined on the unit tetrahedron, and supplied by the user.  
//
//    The user also supplies the name of a subroutine, here called "QUAD_FUN", 
//    which evaluates the integrand at a set of points.  The form is:
//
//      void quad_fun ( int n, double xyz_vec[3*n], double f_vec[n] )
//
//    and it returns in each entry F_VEC(1:N), the value of the integrand
//    at XYZ_VEC(1:3,1:N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes in the tet mesh.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TETRA_ORDER, the order of tetrahedrons in the tet mesh.
//
//    Input, int TETRA_NUM, the number of tetrahedrons in the tet mesh.
//
//    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], indices of the nodes.
//
//    Input, void QUAD_FUN ( int N, double XYZ_VEC[3*N], F_VEC[N] ), the name 
//    of the routine that evaluates the integrand.
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_XYZ[3*QUAD_NUM], the abscissas of the 
//    quadrature rule, in the unit tetrahedron.
//
//    Input, double QUAD_W[QUAD_NUM], the weights of the 
//    quadrature rule.
//
//    Output, double *QUAD_VALUE, the estimate of the integral
//    of F(X,Y) over the region covered by the tet mesh.
//
//    Output, double *REGION_VOLUME, the volume of the region.
//
{
  int i;
  int j;
  int quad;
  double quad_f[quad_num];
  double quad2_xyz[3*quad_num];
  double temp;
  int tet;
  double tetra_volume;
  double tetra_xyz[3*4];

  *quad_value = 0.0;
  *region_volume = 0.0;

  for ( tet = 0; tet < tetra_num; tet++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        tetra_xyz[i+j*3] = node_xyz[i+(tetra_node[j+tet*4]-1)*3];
      }
    }

    tetra_volume = tetrahedron_volume ( tetra_xyz );

    tetrahedron_order4_reference_to_physical ( tetra_xyz, quad_num, 
      quad_xyz, quad2_xyz );

    quad_fun ( quad_num, quad2_xyz, quad_f );

    temp = 0.0;
    for ( quad = 0; quad < quad_num; quad++ )
    {
      temp = temp + quad_w[quad] * quad_f[quad];
    }
    *quad_value = *quad_value + tetra_volume * temp;

    *region_volume = *region_volume + tetra_volume;
  }

  return;
}
//****************************************************************************80

void tet_mesh_quality1 ( int node_num, double node_xyz[], 
  int tetra_order, int tetra_num, int tetra_node[], double *value_min, 
  double *value_mean, double *value_max, double *value_var )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_QUALITY1 returns a tet mesh quality factor.
//
//  Discussion:
//
//    The tet mesh quality measure is the minimum of the 
//    corresponding tetrahedron quality measure, over all tetrahedrons in the 
//    tet mesh.
//
//    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
//    tet mesh, but it simply ignores the extra nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TETRA_ORDER, the order of the tetrahedrons.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
//
//    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
//    the minimum, mean, maximum and variance of the quality measure.
//
{
# define DIM_NUM 3

  int i;
  int j;
  int node;
  int tetra;
  double tetrahedron[DIM_NUM*4];
  double *tetrahedron_quality;

  tetrahedron_quality = new double[tetra_num];

  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      node = tetra_node[j+tetra*tetra_order];
      for ( i = 0; i < DIM_NUM; i++ )
      {
        tetrahedron[i+j*DIM_NUM] = node_xyz[i+(node-1)*DIM_NUM];
      }
    }
    tetrahedron_quality[tetra] = tetrahedron_quality1_3d ( tetrahedron );
  }

  *value_max = r8vec_max ( tetra_num, tetrahedron_quality );
  *value_min = r8vec_min ( tetra_num, tetrahedron_quality );
  *value_mean = r8vec_mean ( tetra_num, tetrahedron_quality );
  *value_var = r8vec_variance ( tetra_num, tetrahedron_quality );

  delete [] tetrahedron_quality;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void tet_mesh_quality2 ( int node_num, double node_xyz[], int tetra_order, 
  int tetra_num, int tetra_node[], double *value_min, double *value_mean, 
  double *value_max, double *value_var )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_QUALITY2 returns a tet mesh quality factor.
//
//  Discussion:
//
//    The tet mesh quality measure is the minimum of the 
//    corresponding tetrahedron quality measure, over all tetrahedrons in the 
//    tet mesh.
//
//    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
//    tet mesh, but it simply ignores the extra nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TETRA_ORDER, the order of the tetrahedrons.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
//
//    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
//    the minimum, mean, maximum and variance of the quality measure.
//
{
# define DIM_NUM 3

  int i;
  int j;
  int node;
  int tetra;
  double tetrahedron[DIM_NUM*4];
  double *tetrahedron_quality;

  tetrahedron_quality = new double[tetra_num];

  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      node = tetra_node[j+tetra*tetra_order];
      for ( i = 0; i < DIM_NUM; i++ )
      {
        tetrahedron[i+j*DIM_NUM] = node_xyz[i+(node-1)*DIM_NUM];
      }
    }
    tetrahedron_quality[tetra] = tetrahedron_quality2_3d ( tetrahedron );
  }

  *value_max = r8vec_max ( tetra_num, tetrahedron_quality );
  *value_min = r8vec_min ( tetra_num, tetrahedron_quality );
  *value_mean = r8vec_mean ( tetra_num, tetrahedron_quality );
  *value_var = r8vec_variance ( tetra_num, tetrahedron_quality );

  delete [] tetrahedron_quality;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void tet_mesh_quality3 ( int node_num, double node_xyz[], int tetra_order, 
  int tetra_num, int tetra_node[], double *value_min, double *value_mean, 
  double *value_max, double *value_var )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_QUALITY3 returns a tet mesh quality factor.
//
//  Discussion:
//
//    The tet mesh quality measure is the minimum of the 
//    corresponding tetrahedron quality measure, over all tetrahedrons in the 
//    tet mesh.
//
//    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
//    tet mesh, but it simply ignores the extra nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TETRA_ORDER, the order of the tetrahedrons.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
//
//    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
//    the minimum, mean, maximum and variance of the quality measure.
//
{
# define DIM_NUM 3

  int i;
  int j;
  int node;
  int tetra;
  double tetrahedron[DIM_NUM*4];
  double *tetrahedron_quality;

  tetrahedron_quality = new double[tetra_num];

  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      node = tetra_node[j+tetra*tetra_order];
      for ( i = 0; i < DIM_NUM; i++ )
      {
        tetrahedron[i+j*DIM_NUM] = node_xyz[i+(node-1)*DIM_NUM];
      }
    }
    tetrahedron_quality[tetra] = tetrahedron_quality3_3d ( tetrahedron );
  }

  *value_max = r8vec_max ( tetra_num, tetrahedron_quality );
  *value_min = r8vec_min ( tetra_num, tetrahedron_quality );
  *value_mean = r8vec_mean ( tetra_num, tetrahedron_quality );
  *value_var = r8vec_variance ( tetra_num, tetrahedron_quality );

  delete [] tetrahedron_quality;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void tet_mesh_quality4 ( int node_num, double node_xyz[], int tetra_order, 
  int tetra_num, int tetra_node[], double *value_min, double *value_mean, 
  double *value_max, double *value_var )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_QUALITY4 returns a tet mesh quality factor.
//
//  Discussion:
//
//    The tet mesh quality measure is the minimum of the 
//    corresponding tetrahedron quality measure, over all tetrahedrons in the 
//    tet mesh.
//
//    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
//    tet mesh, but it simply ignores the extra nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TETRA_ORDER, the order of the tetrahedrons.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
//
//    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
//    the minimum, mean, maximum and variance of the quality measure.
//
{
# define DIM_NUM 3

  int i;
  int j;
  int node;
  int tetra;
  double tetrahedron[DIM_NUM*4];
  double *tetrahedron_quality;

  tetrahedron_quality = new double[tetra_num];

  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      node = tetra_node[j+tetra*tetra_order];
      for ( i = 0; i < DIM_NUM; i++ )
      {
        tetrahedron[i+j*DIM_NUM] = node_xyz[i+(node-1)*DIM_NUM];
      }
    }
    tetrahedron_quality[tetra] = tetrahedron_quality4_3d ( tetrahedron );
  }

  *value_max = r8vec_max ( tetra_num, tetrahedron_quality );
  *value_min = r8vec_min ( tetra_num, tetrahedron_quality );
  *value_mean = r8vec_mean ( tetra_num, tetrahedron_quality );
  *value_var = r8vec_variance ( tetra_num, tetrahedron_quality );

  delete [] tetrahedron_quality;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void tet_mesh_quality5 ( int node_num, double node_xyz[], int tetra_order, 
  int tetra_num, int tetra_node[], double *value_min, double *value_mean, 
  double *value_max, double *value_var )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_QUALITY5 returns a tet mesh quality factor.
//
//  Discussion:
//
//    The tet mesh quality measure is the ratio of the minimum
//    tetrahedron volume to the maximum tetrahedron volume.
//
//    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
//    tet mesh, but it simply ignores the extra nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TETRA_ORDER, the order of the tetrahedrons.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
//
//    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
//    the minimum, mean, maximum and variance of the quality measure.
//
{
# define DIM_NUM 3

  int i;
  int j;
  int node;
  double quality;
  int tetra;
  double tetrahedron[DIM_NUM*4];
  double *tetrahedron_quality;
  double volume_max;

  tetrahedron_quality = new double[tetra_num];

  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      node = tetra_node[j+tetra*tetra_order];
      for ( i = 0; i < DIM_NUM; i++ )
      {
        tetrahedron[i+j*DIM_NUM] = node_xyz[i+(node-1)*DIM_NUM];
      }
    }
    tetrahedron_quality[tetra] = tetrahedron_volume ( tetrahedron );
  }

  volume_max = r8vec_max ( tetra_num, tetrahedron_quality );

  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    tetrahedron_quality[tetra] = tetrahedron_quality[tetra] / volume_max;
  }

  *value_max = r8vec_max ( tetra_num, tetrahedron_quality );
  *value_min = r8vec_min ( tetra_num, tetrahedron_quality );
  *value_mean = r8vec_mean ( tetra_num, tetrahedron_quality );
  *value_var = r8vec_variance ( tetra_num, tetrahedron_quality );

  delete [] tetrahedron_quality;

  return;
# undef DIM_NUM
}
//****************************************************************************80

int tet_mesh_search_delaunay ( int node_num, double node_xyz[], int tet_order, 
  int tet_num, int tet_node[], int tet_neighbor[], double p[], int *face, 
  int *step_num )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_SEARCH_DELAUNAY searches a Delaunay tet mesh for a point.
//
//  Discussion:
//
//    The algorithm "walks" from one tetrahedron to its neighboring tetrahedron,
//    and so on, until a tetrahedron is found containing point P, or P is found
//    to be outside the convex hull.
//
//    The algorithm computes the barycentric coordinates of the point with
//    respect to the current tetrahedron.  If all 4 quantities are positive,
//    the point is contained in the tetrahedron.  If the I-th coordinate is
//    negative, then P lies on the far side of edge I, which is opposite
//    from vertex I.  This gives a hint as to where to search next.
//
//    For a Delaunay tet mesh, the search is guaranteed to terminate.
//    For other meshes, a continue may occur.
//
//    Note the surprising fact that, even for a Delaunay tet mesh of
//    a set of nodes, the nearest node to P need not be one of the
//    vertices of the tetrahedron containing P.
//
//    The code can be called for tet meshes of any order, but only
//    the first 4 nodes in each tetrahedron are considered.  Thus, if
//    higher order tetrahedrons are used, and the extra nodes are intended
//    to give the tetrahedron a polygonal shape, these will have no effect,
//    and the results obtained here might be misleading.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2009
//
//  Author:
//
//    John Burkardt.
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
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of 
//    the nodes.
//
//    Input, int TET_ORDER, the order of the tetrahedrons.
//
//    Input, int TET_NUM, the number of tetrahedrons.
//
//    Input, int TET_NODE[TET_ORDER*TET_NUM],
//    the nodes that make up each tetrahedron.
//
//    Input, int TET_NEIGHBOR[4*TET_NUM], the 
//    tetrahedron neighbor list.
//
//    Input, double P[3], the coordinates of a point.
//
//    Output, int *FACE, indicates the position of the point P in
//    face TET_INDEX:
//    0, the interior or boundary of the tetrahedron;
//    -1, outside the convex hull of the tet mesh, past face 1;
//    -2, outside the convex hull of the tet mesh, past face 2;
//    -3, outside the convex hull of the tet mesh, past face 3.
//    -4, outside the convex hull of the tet mesh, past face 4.
//
//    Output, int *STEP_NUM, the number of steps taken.
//
//    Output, int TET_MESH_SEARCH_DELAUNAY, the index of the tetrahedron 
//    where the search ended.  If a cycle occurred, then -1 is returned.
//
{
  double *alpha;
  int i;
  int j;
  int k;
  int tet_index;
  double tet_xyz[3*4];
  static int tet_index_save = -1;
//
//  If possible, start with the previous successful value of TET_INDEX.
//
  if ( tet_index_save < 1 || tet_num < tet_index_save )
  {
    tet_index = ( tet_num + 1 ) / 2;
  }
  else
  {
    tet_index = tet_index_save;
  }

  *step_num = -1;
  *face = 0;

  for ( ; ; )
  {
    *step_num = *step_num + 1;

    if ( tet_num < *step_num )
    {
      cerr << "\n";
      cerr << "TET_MESH_SEARCH_DELAUNAY - Fatal error!\n";
      cerr << "  The algorithm seems to be cycling.\n";
      tet_index = -1;
      *face = -1;
      exit ( 1 );
    }

    for ( j = 0; j < 4; j++ )
    {
      k = tet_node[j+tet_index*4];
      for ( i = 0; i < 3; i++ )
      {
        tet_xyz[i+j*3] = node_xyz[i+k*3];
      }
    }

    alpha = tetrahedron_barycentric ( tet_xyz, p );
//
//  If the barycentric coordinates are all positive, then the point
//  is inside the tetrahedron and we're done.
//
    if ( 0.0 <= alpha[0] && 0.0 <= alpha[1] && 0.0 <= alpha[2] && 0.0 <= alpha[3] )
    {
      break;
    }
//
//  At least one barycentric coordinate is negative.
//
//  If there is a negative barycentric coordinate for which there exists an
//  opposing tetrahedron neighbor closer to the point, move to that tetrahedron.
//
    if ( alpha[0] < 0.0 && 0 < tet_neighbor[0+tet_index*4] )
    {
      tet_index = tet_neighbor[0+tet_index*4];
      continue;
    }
    else if ( alpha[1] < 0.0 && 0 < tet_neighbor[1+tet_index*4] )
    {
      tet_index = tet_neighbor[1+tet_index*4];
      continue;
    }
    else if ( alpha[2] < 0.0 && 0 < tet_neighbor[2+tet_index*4] )
    {
      tet_index = tet_neighbor[2+tet_index*4];
      continue;
    }
    else if ( alpha[3] < 0.0 && 0 < tet_neighbor[3+tet_index*4] )
    {
      tet_index = tet_neighbor[3+tet_index*4];
      continue;
    }
//
//  All negative barycentric coordinates correspond to vertices opposite
//  faces on the convex hull.
//
//  Note the face and exit.
//
    if ( alpha[0] < 0.0 )
    {
      *face = -1;
      break;
    }
    else if ( alpha[1] < 0.0 )
    {
      *face = -2;
      break;
    }
    else if ( alpha[2] < 0.0 )
    {
      *face = -3;
      break;
    }
    else if ( alpha[3] < 0.0 )
    {
      *face = -4;
      break;
    }
  }

  tet_index_save = tet_index;

  return tet_index;
}
//****************************************************************************80

int tet_mesh_search_naive ( int node_num, double node_xyz[],
  int tet_order, int tet_num, int tet_node[], double p[], int *step_num )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_SEARCH_NAIVE naively searches a tet mesh.
//
//  Discussion:
//
//    The algorithm simply checks each tetrahedron to see if point P is
//    contained in it.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates 
//    of the nodes.
//
//    Input, int TET_ORDER, the order of the tetrahedrons.
//
//    Input, int TET_NUM, the number of tetrahedrons in
//    the mesh.
//
//    Input, int TET_NODE[TET_ORDER*TET_NUM], 
//    the nodes that make up each tetrahedron.
//
//    Input, double P[3], the coordinates of a point.
//
//    Output, int TET_MESH_ORDER4_SEARCH_NAIE, the index of the tetrahedron
//    where the search ended, or -1 if no tetrahedron was found containing
//    the point.
//
//    Output, int *STEP_NUM, the number of tetrahedrons examined.
{
  double *alpha;
  int i;
  int j;
  int tet;
  int tet_index;
  double tet_xyz[3*4];

  tet_index = -1;
  *step_num = 0;

  for ( tet = 0; tet < tet_num; tet++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        tet_xyz[i+j*3] = node_xyz[i+tet_node[j+tet*4]*3];
      }
    }
    alpha = tetrahedron_barycentric ( tet_xyz, p );

    if ( r8vec_is_nonnegative ( 4, alpha ) )
    {
      tet_index = tet;
      *step_num = tet;
      return tet_index;
    }

    delete [] alpha;
  }

  return tet_index;
}
//****************************************************************************80

double *tetrahedron_barycentric ( double tetra[3*4], double p[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_BARYCENTRIC returns the barycentric coordinates of a point.
//
//  Discussion:
//
//    The barycentric coordinates of a point P with respect to
//    a tetrahedron are a set of four values C(1:4), each associated
//    with a vertex of the tetrahedron.  The values must sum to 1.
//    If all the values are between 0 and 1, the point is contained
//    within the tetrahedron.
//
//    The barycentric coordinate of point X related to vertex A can be
//    interpreted as the ratio of the volume of the tetrahedron with 
//    vertex A replaced by vertex X to the volume of the original 
//    tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Input, double P[3], the point to be checked.
//
//    Output, double C[4], the barycentric coordinates of the point with
//    respect to the tetrahedron.
//
{
# define N 3
# define RHS_NUM 1

  double a[N*(N+RHS_NUM)];
  double *c;
  int info;
//
//  Set up the linear system
//
//    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
//    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
//    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
//
//  which is satisfied by the barycentric coordinates.
//

  a[0+0*N] = tetra[0+1*3] - tetra[0+0*3];
  a[1+0*N] = tetra[1+1*3] - tetra[1+0*3];
  a[2+0*N] = tetra[2+1*3] - tetra[2+0*3];

  a[0+1*N] = tetra[0+2*3] - tetra[0+0*3];
  a[1+1*N] = tetra[1+2*3] - tetra[1+0*3];
  a[2+1*N] = tetra[2+2*3] - tetra[2+0*3];

  a[0+2*N] = tetra[0+3*3] - tetra[0+0*3];
  a[1+2*N] = tetra[1+3*3] - tetra[1+0*3];
  a[2+2*N] = tetra[2+3*3] - tetra[2+0*3];

  a[0+3*N] = p[0]         - tetra[0+0*3];
  a[1+3*N] = p[1]         - tetra[1+0*3];
  a[2+3*N] = p[2]         - tetra[2+0*3];
//
//  Solve the linear system.
//
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TETRAHEDRON_BARYCENTRIC - Fatal error!\n";
    cout << "  The linear system is singular.\n";
    cout << "  The input data does not form a proper tetrahedron.\n";
    exit ( 1 );
  }

  c = new double[4];

  c[1] = a[0+3*N];
  c[2] = a[1+3*N];
  c[3] = a[2+3*N];

  c[0] = 1.0 - c[1] - c[2] - c[3];

  return c;
# undef N
# undef RHS_NUM
}
//****************************************************************************80

void tetrahedron_circumsphere_3d ( double tetra[3*4], double *r, double pc[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere of a tetrahedron in 3D.
//
//  Discussion:
//
//    The circumsphere, or circumscribed sphere, of a tetrahedron is the sphere that
//    passes through the four vertices.  The circumsphere is not necessarily
//    the smallest sphere that contains the tetrahedron.
//
//    Surprisingly, the diameter of the sphere can be found by solving
//    a 3 by 3 linear system.  This is because the vectors P2 - P1,
//    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
//    right triangle with the diameter through P1.  Hence, the dot product of
//    P2 - P1 with that diameter is equal to the square of the length
//    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
//    the diameter vector originating at P1, and hence the radius and
//    center.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Output, double *R, PC[3], the coordinates of the center of the
//    circumscribed sphere, and its radius.  If the linear system is
//    singular, then R = -1, PC[] = 0.
//
{
# define DIM_NUM 3
# define RHS_NUM 1

  double a[DIM_NUM*(DIM_NUM+RHS_NUM)];
  int info;
//
//  Set up the linear system.
//
  a[0+0*3] = tetra[0+1*3] - tetra[0+0*3];
  a[0+1*3] = tetra[1+1*3] - tetra[1+0*3];
  a[0+2*3] = tetra[2+1*3] - tetra[2+0*3];
  a[0+3*3] = pow ( tetra[0+1*3] - tetra[0+0*3], 2 ) 
           + pow ( tetra[1+1*3] - tetra[1+0*3], 2 ) 
           + pow ( tetra[2+1*3] - tetra[2+0*3], 2 );

  a[1+0*3] = tetra[0+2*3] - tetra[0+0*3];
  a[1+1*3] = tetra[1+2*3] - tetra[1+0*3];
  a[1+2*3] = tetra[2+2*3] - tetra[2+0*3];
  a[1+3*3] = pow ( tetra[0+2*3] - tetra[0+0*3], 2 ) 
           + pow ( tetra[1+2*3] - tetra[1+0*3], 2 ) 
           + pow ( tetra[2+2*3] - tetra[2+0*3], 2 );

  a[2+0*3] = tetra[0+3*3] - tetra[0+0*3];
  a[2+1*3] = tetra[1+3*3] - tetra[1+0*3];
  a[2+2*3] = tetra[2+3*3] - tetra[2+0*3];
  a[2+3*3] = pow ( tetra[0+3*3] - tetra[0+0*3], 2 ) 
           + pow ( tetra[1+3*3] - tetra[1+0*3], 2 ) 
           + pow ( tetra[2+3*3] - tetra[2+0*3], 2 );
//
//  Solve the linear system.
//
  info = r8mat_solve ( DIM_NUM, RHS_NUM, a );
//
//  If the system was singular, return a consolation prize.
//
  if ( info != 0 )
  {
    *r = -1.0;
    r8vec_zero ( DIM_NUM, pc );
    return;
  }
//
//  Compute the radius and center.
//
  *r = 0.5 * sqrt 
    ( a[0+3*3] * a[0+3*3] 
    + a[1+3*3] * a[1+3*3] 
    + a[2+3*3] * a[2+3*3] );

  pc[0] = tetra[0+0*3] + 0.5 * a[0+3*3];
  pc[1] = tetra[1+0*3] + 0.5 * a[1+3*3];
  pc[2] = tetra[2+0*3] + 0.5 * a[2+3*3];

  return;
# undef DIM_NUM
# undef RHS_NUM
}
//****************************************************************************80

double *tetrahedron_edge_length_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_EDGE_LENGTH_3D returns edge lengths of a tetrahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Output, double EDGE_LENGTH[6], the length of the edges.
//
{
# define DIM_NUM 3

  double *edge_length;
  int i;
  int j1;
  int j2;
  int k;
  double v[DIM_NUM];

  edge_length = new double[6];

  k = 0;
  for ( j1 = 0; j1 < 3; j1++ )
  {
    for ( j2 = j1 + 1; j2 < 4; j2++ )
    {
      for ( i = 0; i < DIM_NUM; i++ )
      {
        v[i] = tetra[i+j2*DIM_NUM] - tetra[i+j1*DIM_NUM];
      }
      edge_length[k] = r8vec_length ( DIM_NUM, v );
      k = k + 1;
    }
  }

  return edge_length;
# undef DIM_NUM
}
//****************************************************************************80

void tetrahedron_insphere_3d ( double tetra[3*4], double *r, double pc[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_INSPHERE_3D finds the insphere of a tetrahedron in 3D.
//
//  Discussion:
//
//    The insphere of a tetrahedron is the inscribed sphere, which touches
//    each face of the tetrahedron at a single point.
//
//    The points of contact are the centroids of the triangular faces
//    of the tetrahedron.  Therefore, the point of contact for a face
//    can be computed as the average of the vertices of that face.
//
//    The sphere can then be determined as the unique sphere through
//    the four given centroids.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Schneider, David Eberly,
//    Geometric Tools for Computer Graphics,
//    Elsevier, 2002,
//    ISBN: 1558605940,
//    LC: T385.G6974.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Output, double *R, PC[3], the radius and the center
//    of the sphere.
//
{
# define DIM_NUM 3

  double b[4*4];
  double gamma;
  int i;
  int j;
  double l123;
  double l124;
  double l134;
  double l234;
  double *n123;
  double *n124;
  double *n134;
  double *n234;
  double v21[DIM_NUM];
  double v31[DIM_NUM];
  double v41[DIM_NUM];
  double v32[DIM_NUM];
  double v42[DIM_NUM];
  double v43[DIM_NUM];

  for ( i = 0; i < DIM_NUM; i++ )
  {
    v21[i] = tetra[i+1*DIM_NUM] - tetra[i+0*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v31[i] = tetra[i+2*DIM_NUM] - tetra[i+0*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v41[i] = tetra[i+3*DIM_NUM] - tetra[i+0*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v32[i] = tetra[i+2*DIM_NUM] - tetra[i+1*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v42[i] = tetra[i+3*DIM_NUM] - tetra[i+1*DIM_NUM];
  }
  for ( i = 0; i < DIM_NUM; i++ )
  {
    v43[i] = tetra[i+3*DIM_NUM] - tetra[i+2*DIM_NUM];
  }

  n123 = r8vec_cross_3d ( v21, v31 );
  n124 = r8vec_cross_3d ( v41, v21 );
  n134 = r8vec_cross_3d ( v31, v41 );
  n234 = r8vec_cross_3d ( v42, v32 );

  l123 = r8vec_length ( DIM_NUM, n123 );
  l124 = r8vec_length ( DIM_NUM, n124 );
  l134 = r8vec_length ( DIM_NUM, n134 );
  l234 = r8vec_length ( DIM_NUM, n234 );

  delete [] n123;
  delete [] n124;
  delete [] n134;
  delete [] n234;

  for ( i = 0; i < DIM_NUM; i++ )
  {
    pc[i] = ( l234 * tetra[i+0*DIM_NUM]
            + l134 * tetra[i+1*DIM_NUM]
            + l124 * tetra[i+2*DIM_NUM]
            + l123 * tetra[i+3*DIM_NUM] )
            / ( l234 + l134 + l124 + l123 );
  }

  for ( j = 0; j < 4; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      b[i+j*4] = tetra[i+j*DIM_NUM];
    }
    b[3+j*4] = 1.0;
  }
  
  gamma = fabs ( r8mat_det_4d ( b ) );

  *r = gamma / ( l234 + l134 + l124 + l123 );

  return;
# undef DIM_NUM
}
//****************************************************************************80

void tetrahedron_order4_physical_to_reference ( double tetra[], int n, 
  double phy[], double ref[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE maps physical points to reference points.
//
//  Discussion:
//
//    Given the vertices of an order 4 physical tetrahedron and a point
//    (X,Y,Z) in the physical tetrahedron, the routine computes the value
//    of the corresponding image point (R,S,T) in reference space.
//
//    This routine may be appropriate for an order 10 tetrahedron,
//    if the mapping between reference and physical space is linear.  
//    This implies, in particular, that the edges of the image tetrahedron 
//    are straight, the faces are flat, and the "midside" nodes in the
//    physical tetrahedron are halfway along the sides of the physical 
//    tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.  
//    The vertices are assumed to be the images of
//    (0,0,0), (1,0,0), (0,1,0) and (0,0,1) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double PHY[3*N], the coordinates of physical points
//    to be transformed.
//
//    Output, double REF[3*N], the coordinates of the corresponding
//    points in the reference space.
//
{
  double a[3*3];
  double det;
  int i;
  int j;
//
//  Set up the matrix.
//
  for ( i = 0; i < 3; i++ )
  {
    a[i+0*3] = tetra[i+1*3] - tetra[i+0*3];
    a[i+1*3] = tetra[i+2*3] - tetra[i+0*3];
    a[i+2*3] = tetra[i+3*3] - tetra[i+0*3];
  }
//
//  Compute the determinant.
//
  det =  a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] ) 
       + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] ) 
       + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );
//
//  If the determinant is zero, bail out.
//
  if ( det == 0.0 )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        ref[i+j*3] = 0.0;
      }
    }
    return;
  }
//
//  Compute the solution.
//
  for ( j = 0; j < n; j++ )
  {
    ref[0+j*3] = (   ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] ) 
                     * ( phy[0+j*3] - tetra[0+0*3] ) 
                   - ( a[0+1*3] * a[2+2*3] - a[0+2*3] * a[2+1*3] ) 
                     * ( phy[1+j*3] - tetra[1+0*3] ) 
                   + ( a[0+1*3] * a[1+2*3] - a[0+2*3] * a[1+1*3] ) 
                     * ( phy[2+j*3] - tetra[2+0*3] ) 
                 ) / det;

    ref[1+j*3] = ( - ( a[1+0*3] * a[2+2*3] - a[1+2*3] * a[2+0*3] ) 
                     * ( phy[0+j*3] - tetra[0+0*3] ) 
                   + ( a[0+0*3] * a[2+2*3] - a[0+2*3] * a[2+0*3] ) 
                     * ( phy[1+j*3] - tetra[1+0*3] ) 
                   - ( a[0+0*3] * a[1+2*3] - a[0+2*3] * a[1+0*3] ) 
                     * ( phy[2+j*3] - tetra[2+0*3] ) 
                 ) / det;

    ref[2+j*3] = (   ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] ) 
                     * ( phy[0+j*3] - tetra[0+0*3] ) 
                   - ( a[0+0*3] * a[2+1*3] - a[0+1*3] * a[2+0*3] ) 
                     * ( phy[1+j*3] - tetra[1+0*3] ) 
                   + ( a[0+0*3] * a[1+1*3] - a[0+1*3] * a[1+0*3] ) 
                     * ( phy[2+j*3] - tetra[2+0*3] ) 
                 ) / det;
  }
  return;
}
//****************************************************************************80

void tetrahedron_order4_reference_to_physical ( double tetra[], int n, 
  double ref[], double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL maps reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 4 physical tetrahedron and a point
//    (R,S,T) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y,Z) in physical space.
//
//    This routine will also be correct for an order 10 tetrahedron,
//    if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image tetrahedron are straight, the faces are flat, and 
//    the "midside" nodes in the physical tetrahedron are
//    halfway along the edges of the physical tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0,0), (1,0,0),
//    (0,1,0) and (0,0,1) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double REF[3*N], points in the reference tetrahedron
//
//    Output, double PHY[3*N], corresponding points in the
//    physical tetrahedron.
//
{
  int i;
  int j;

  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*3] = tetra[i+0*3] * ( 1.0 - ref[0+j*3] - ref[1+j*3] - ref[2+j*3] ) 
                 + tetra[i+1*3] *       + ref[0+j*3]                
                 + tetra[i+2*3] *                    + ref[1+j*3]
                 + tetra[i+3*3] *                                 + ref[2+j*3];
    }
  }

  return;
}
//****************************************************************************80

double tetrahedron_quality1_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_QUALITY1_3D: "quality" of a tetrahedron in 3D.
//
//  Discussion:
//
//    The quality of a tetrahedron is 3.0 times the ratio of the radius of
//    the inscribed sphere divided by that of the circumscribed sphere.
//
//    An equilateral tetrahredron achieves the maximum possible quality of 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Output, double TETRAHEDRON_QUALITY1_3D, the quality of the tetrahedron.
//
{
# define DIM_NUM 3

  double pc[DIM_NUM];
  double quality;
  double r_in;
  double r_out;

  tetrahedron_circumsphere_3d ( tetra, &r_out, pc );

  tetrahedron_insphere_3d ( tetra, &r_in, pc );

  quality = 3.0 * r_in / r_out;

  return quality;
# undef DIM_NUM
}
//****************************************************************************80

double tetrahedron_quality2_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_QUALITY2_3D: "quality" of a tetrahedron in 3D.
//
//  Discussion:
//
//    The quality measure #2 of a tetrahedron is:
//
//      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
//
//    where
//
//      RIN = radius of the inscribed sphere;
//      LMAX = length of longest side of the tetrahedron.
//
//    An equilateral tetrahredron achieves the maximum possible quality of 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Desheng Wang,
//    The Optimal Centroidal Voronoi Tesselations and the Gersho's
//    Conjecture in the Three-Dimensional Space,
//    Computers and Mathematics with Applications,
//    Volume 49, 2005, pages 1355-1373.
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Output, double TETRAHEDRON_QUALITY2_3D, the quality of the tetrahedron.
//
{
# define DIM_NUM 3

  double *edge_length;
  double l_max;
  double pc[DIM_NUM];
  double quality2;
  double r_in;

  edge_length = tetrahedron_edge_length_3d ( tetra );

  l_max = r8vec_max ( 6, edge_length );

  tetrahedron_insphere_3d ( tetra, &r_in, pc );

  quality2 = 2.0 * sqrt ( 6.0 ) * r_in / l_max;

  delete [] edge_length;

  return quality2;
# undef DIM_NUM
}
//****************************************************************************80

double tetrahedron_quality3_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_QUALITY3_3D computes the mean ratio of a tetrahedron.
//
//  Discussion:
//
//    This routine computes QUALITY3, the eigenvalue or mean ratio of
//    a tetrahedron.
//
//      QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
//
//    This value may be used as a shape quality measure for the tetrahedron.
//
//    For an equilateral tetrahedron, the value of this quality measure
//    will be 1.  For any other tetrahedron, the value will be between
//    0 and 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2005
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
//    Input, double TETRA(3,4), the coordinates of the vertices.
//
//    Output, double TETRAHEDRON_QUALITY3_3D, the mean ratio of the tetrahedron.
//
{
# define DIM_NUM 3

  double ab[DIM_NUM];
  double ac[DIM_NUM];
  double ad[DIM_NUM];
  double bc[DIM_NUM];
  double bd[DIM_NUM];
  double cd[DIM_NUM];
  double denom;
  int i;
  double lab;
  double lac;
  double lad;
  double lbc;
  double lbd;
  double lcd;
  double quality3;
  double volume;
//
//  Compute the vectors representing the sides of the tetrahedron.
//
  for ( i = 0; i < DIM_NUM; i++ )
  {
    ab[i] = tetra[i+1*DIM_NUM] - tetra[i+0*DIM_NUM];
    ac[i] = tetra[i+2*DIM_NUM] - tetra[i+0*DIM_NUM];
    ad[i] = tetra[i+3*DIM_NUM] - tetra[i+0*DIM_NUM];
    bc[i] = tetra[i+2*DIM_NUM] - tetra[i+1*DIM_NUM];
    bd[i] = tetra[i+3*DIM_NUM] - tetra[i+1*DIM_NUM];
    cd[i] = tetra[i+3*DIM_NUM] - tetra[i+2*DIM_NUM];
  }
//
//  Compute the squares of the lengths of the sides.
//
  lab = pow ( ab[0], 2 ) + pow ( ab[1], 2 ) + pow ( ab[2], 2 );
  lac = pow ( ac[0], 2 ) + pow ( ac[1], 2 ) + pow ( ac[2], 2 );
  lad = pow ( ad[0], 2 ) + pow ( ad[1], 2 ) + pow ( ad[2], 2 );
  lbc = pow ( bc[0], 2 ) + pow ( bc[1], 2 ) + pow ( bc[2], 2 );
  lbd = pow ( bd[0], 2 ) + pow ( bd[1], 2 ) + pow ( bd[2], 2 );
  lcd = pow ( cd[0], 2 ) + pow ( cd[1], 2 ) + pow ( cd[2], 2 );
//
//  Compute the volume.
//
  volume = fabs ( 
      ab[0] * ( ac[1] * ad[2] - ac[2] * ad[1] ) 
    + ab[1] * ( ac[2] * ad[0] - ac[0] * ad[2] ) 
    + ab[2] * ( ac[0] * ad[1] - ac[1] * ad[0] ) ) / 6.0;

  denom = lab + lac + lad + lbc + lbd + lcd;

  if ( denom == 0.0 )
  {
    quality3 = 0.0;
  }
  else
  {
    quality3 = 12.0 * pow ( 3.0 * volume, 2.0 / 3.0 ) / denom;
  }

  return quality3;
# undef DIM_NUM
}
//****************************************************************************80

double tetrahedron_quality4_3d ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_QUALITY4_3D computes the minimum solid angle of a tetrahedron.
//
//  Discussion:
//
//    This routine computes a quality measure for a tetrahedron, based
//    on the sine of half the minimum of the four solid angles.
//
//  Modified:
//
//    17 August 2005
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
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Output, double QUALITY4, the value of the quality measure.
//
{
# define DIM_NUM 3

  double ab[DIM_NUM];
  double ac[DIM_NUM];
  double ad[DIM_NUM];
  double bc[DIM_NUM];
  double bd[DIM_NUM];
  double cd[DIM_NUM];
  double denom;
  int i;
  double l1;
  double l2;
  double l3;
  double lab;
  double lac;
  double lad;
  double lbc;
  double lbd;
  double lcd;
  double quality4;
  double volume;
//
//  Compute the vectors that represent the sides.
//
  for ( i = 0; i < DIM_NUM; i++ )
  {
    ab[i] = tetra[i+1*DIM_NUM] - tetra[i+0*DIM_NUM];
    ac[i] = tetra[i+2*DIM_NUM] - tetra[i+0*DIM_NUM];
    ad[i] = tetra[i+3*DIM_NUM] - tetra[i+0*DIM_NUM];
    bc[i] = tetra[i+2*DIM_NUM] - tetra[i+1*DIM_NUM];
    bd[i] = tetra[i+3*DIM_NUM] - tetra[i+1*DIM_NUM];
    cd[i] = tetra[i+3*DIM_NUM] - tetra[i+2*DIM_NUM];
  }
//
//  Compute the lengths of the sides.
//
  lab = r8vec_length ( DIM_NUM, ab );
  lac = r8vec_length ( DIM_NUM, ac );
  lad = r8vec_length ( DIM_NUM, ad );
  lbc = r8vec_length ( DIM_NUM, bc );
  lbd = r8vec_length ( DIM_NUM, bd );
  lcd = r8vec_length ( DIM_NUM, cd );
//
//  Compute the volume.
//
  volume = fabs ( 
      ab[0] * ( ac[1] * ad[2] - ac[2] * ad[1] ) 
    + ab[1] * ( ac[2] * ad[0] - ac[0] * ad[2] ) 
    + ab[2] * ( ac[0] * ad[1] - ac[1] * ad[0] ) ) / 6.0;

  quality4 = 1.0;

  l1 = lab + lac;
  l2 = lab + lad;
  l3 = lac + lad;

  denom = ( l1 + lbc ) * ( l1 - lbc ) 
        * ( l2 + lbd ) * ( l2 - lbd ) 
        * ( l3 + lcd ) * ( l3 - lcd );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = r8_min ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  l1 = lab + lbc;
  l2 = lab + lbd;
  l3 = lbc + lbd;

  denom = ( l1 + lac ) * ( l1 - lac ) 
        * ( l2 + lad ) * ( l2 - lad ) 
        * ( l3 + lcd ) * ( l3 - lcd );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = r8_min ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  l1 = lac + lbc;
  l2 = lac + lcd;
  l3 = lbc + lcd;

  denom = ( l1 + lab ) * ( l1 - lab ) 
        * ( l2 + lad ) * ( l2 - lad ) 
        * ( l3 + lbd ) * ( l3 - lbd );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = r8_min ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  l1 = lad + lbd;
  l2 = lad + lcd;
  l3 = lbd + lcd;

  denom = ( l1 + lab ) * ( l1 - lab ) 
        * ( l2 + lac ) * ( l2 - lac ) 
        * ( l3 + lbc ) * ( l3 - lbc );

  if ( denom <= 0.0 )
  {
    quality4 = 0.0;
  }
  else
  {
    quality4 = r8_min ( quality4, 12.0 * volume / sqrt ( denom ) );
  }

  quality4 = quality4 * 1.5 * sqrt ( 6.0 );

  return quality4;
# undef DIM_NUM
}
//****************************************************************************80

void tetrahedron_reference_sample ( int n, int *seed, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_REFERENCE_SAMPLE samples points in the reference tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the  number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double P[3*N], random points in the tetrahedron.
//
{
  double alpha;
  double beta;
  double gamma;
  int j;
  double r;

  for ( j = 0; j < n; j++ )
  {
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the tetrahedron's volume.
//
//  Imagine a plane, parallel to face 1, so that the volume between
//  vertex 1 and the plane is R percent of the full tetrahedron volume.
//
//  The plane will intersect sides 12, 13, and 14 at a fraction
//  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
//
    alpha = pow ( r, 1.0 / 3.0 );
//
//  Determine the coordinates of the points on sides 12, 13 and 14 intersected
//  by the plane, which form a triangle TR.
//
//  Now choose, uniformly at random, a point in this triangle.
//
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
    beta = sqrt ( r );
//
//  Determine the coordinates of the points on sides 2 and 3 intersected
//  by line L.
//
//  Now choose, uniformly at random, a point on the line L.
//
    gamma = r8_uniform_01 ( seed );

    p[0+j*3] = alpha * ( 1.0 - beta ) *         gamma;
    p[1+j*3] = alpha *         beta   * ( 1.0 - gamma );
    p[2+j*3] = alpha *         beta   *         gamma;
  }

  return;
}
//****************************************************************************80

void tetrahedron_sample ( double tetra[3*4], int n, int *seed, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_SAMPLE returns random points in a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double P[3*N], random points in the tetrahedron.
//
{
# define DIM_NUM 3

  double alpha;
  double beta;
  double gamma;
  int i;
  int j;
  int k;
  double *p12;
  double *p13;
  double r;
  double *t;

  p12 = new double[DIM_NUM];
  p13 = new double[DIM_NUM];
  t = new double[DIM_NUM*3];

  for ( k = 0; k < n; k++ )
  {
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the tetrahedron's volume.
//
//  Imagine a plane, parallel to face 1, so that the volume between
//  vertex 1 and the plane is R percent of the full tetrahedron volume.
//
//  The plane will intersect sides 12, 13, and 14 at a fraction
//  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
//
    alpha = pow ( r, 1.0 / 3.0 );
//
//  Determine the coordinates of the points on sides 12, 13 and 14 intersected
//  by the plane, which form a triangle TR.
//
    for ( i = 0; i < DIM_NUM; i++ )
    {
      for ( j = 0; j < 3; j++ )
      {
        t[i+j*3] = ( 1.0 - alpha ) * tetra[i+0*3] 
                 +         alpha   * tetra[i+(j+1)*3];
      }
    }
//
//  Now choose, uniformly at random, a point in this triangle.
//
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
    beta = sqrt ( r );
//
//  Determine the coordinates of the points on sides 2 and 3 intersected
//  by line L.
//
    for ( i = 0; i < DIM_NUM; i++ )
    {
      p12[i] = ( 1.0 - beta ) * t[i+0*3] 
             +         beta   * t[i+1*3];

      p13[i] = ( 1.0 - beta ) * t[i+0*3] 
             +         beta   * t[i+2*3];
    }
//
//  Now choose, uniformly at random, a point on the line L.
//
    gamma = r8_uniform_01 ( seed );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      p[i+k*3] = gamma * p12[i] + ( 1.0 - gamma ) * p13[i];
    }
  }

  delete [] p12;
  delete [] p13;
  delete [] t;

  return;
# undef DIM_NUM
}
//****************************************************************************80

double tetrahedron_volume ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Output, double TETRAHEDRON_VOLUME, the volume of the tetrahedron.
//
{
  double a[4*4];
  int i;
  int j;
  double volume;

  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < 4; j++ )
    { 
      a[i+j*4] = tetra[i+j*3];
    }
  }

  i = 3;
  for ( j = 0; j < 4; j++ )
  {
    a[i+j*4] = 1.0;
  }

  volume = fabs ( r8mat_det_4d ( a ) ) / 6.0;

  return volume;
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
