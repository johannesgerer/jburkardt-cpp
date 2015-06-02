# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "lagrange_nd.hpp"

//****************************************************************************80

int *comp_unrank_grlex ( int kc, int rank )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_UNRANK_GRLEX computes the composition of given grlex rank.
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
//    Input, int RANK, the rank of the composition.
//    1 <= RANK.
//
//    Output, int COMP_UNRANK_GRLEX[KC], the composition XC of the given rank.
//    For each I, 0 <= XC[I] <= NC, and 
//    sum ( 1 <= I <= KC ) XC[I] = NC.
//
{
  int i;
  int j;
  int ks;
  int nc;
  int nksub;
  int ns;
  int r;
  int rank1;
  int rank2;
  int *xc;
  int *xs;
//
//  Ensure that 1 <= KC.
//
  if ( kc < 1 )
  {
    cerr << "\n";
    cerr << "COMP_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  KC < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 1 <= RANK.
//
  if ( rank < 1 )
  {
    cerr << "\n";
    cerr << "COMP_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  RANK < 1\n";
    exit ( 1 );
  }
//
//  Determine the appropriate value of NC.
//  Do this by adding up the number of compositions of sum 0, 1, 2, 
//  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
//  gives you the rank of the composition within the set of compositions
//  of sum NC.  And that's the number you need in order to do the
//  unranking.
//
  rank1 = 1;
  nc = -1;
  for ( ; ; )
  {
    nc = nc + 1;
    r = i4_choose ( nc + kc - 1, nc );
    if ( rank < rank1 + r )
    {
      break;
    }
    rank1 = rank1 + r;
  }

  rank2 = rank - rank1;
//
//  Convert to KSUBSET format.
//  Apology: an unranking algorithm was available for KSUBSETS,
//  but not immediately for compositions.  One day we will come back
//  and simplify all this.
//
  ks = kc - 1;
  ns = nc + kc - 1;
  xs = new int[ks];

  nksub = i4_choose ( ns, ks );

  j = 1;

  for ( i = 1; i <= ks; i++ )
  {
    r = i4_choose ( ns - j, ks - i );

    while ( r <= rank2 && 0 < r )
    {
      rank2 = rank2 - r;
      j = j + 1;
      r = i4_choose ( ns - j, ks - i );
    }
    xs[i-1] = j;
    j = j + 1;
  }
//
//  Convert from KSUBSET format to COMP format.
//
  xc = new int[kc];
  xc[0] = xs[0] - 1;
  for ( i = 2; i < kc; i++ )
  {
    xc[i-1] = xs[i-1] - xs[i-2] - 1;
  }
  xc[kc-1] = ns - xs[ks-1];

  delete [] xs;

  return xc;
}
//****************************************************************************80

int i4_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE computes the binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, int I4_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = k;
  if ( n - k < mn )
  {
    mn = n - k;
  }

  if ( mn < 0 )
  {
    value = 0;
  }
  else if ( mn == 0 )
  {
    value = 1;
  }
  else
  {
    mx = k;
    if ( mx < n - k )
    {
      mx = n - k;
    }
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
    }
  }

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

void i4vec_concatenate ( int n1, int a[], int n2, int b[], int c[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_CONCATENATE concatenates two I4VEC's.
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
//    22 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, the number of entries in the first vector.
//
//    Input, int A[N1], the first vector.
//
//    Input, int N2, the number of entries in the second vector.
//
//    Input, int B[N2], the second vector.
//
//    Output, int C[N1+N2], the concatenated vector.
//
{
  int i;

  for ( i = 0; i < n1; i++ )
  {
    c[i] = a[i];
  }
  for ( i = 0; i < n2; i++ )
  {
    c[n1+i] = b[i];
  }

  return;
}
//****************************************************************************80

void i4vec_permute ( int n, int p[], int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PERMUTE permutes an I4VEC in place.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    This routine permutes an array of integer "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   1,   3,   4,   0,   2 )
//      A = (   1,   2,   3,   4,   5 )
//
//    Output:
//
//      A    = (   2,   4,   5,   1,   3 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.
//
//    Input/output, int A[N], the array to be permuted.
//
{
  int a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  perm_check0 ( n, p );
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is 0.
//  So temporarily add 1 to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1;
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "I4VEC_PERMUTE - Fatal error!\n";
          cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
//
//  Restore the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1;
  }

  return;
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

int *i4vec_sort_heap_index_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
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
//      i4vec_permute ( n, indx, a )
//
//    after which a(*) is sorted.
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
//    Input, int A[N], an array to be index-sorted.
//
//    Output, int I4VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  int aval;
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

double *interpolant_value ( int d, int r, int pn, int po[], double pc[], 
  int pe[], double pd[], int ni, double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    INTERPOLANT_VALUE evaluates a Lagrange interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int R, the maximum number of terms in a polynomial.
//
//    Input, int PN, the number of polynomials.
//
//    Input, int PO[PN], the "order" of the polynomials.
//
//    Input, double PC[PN*R], the coefficients of the polynomial.
//
//    Input, int PE[PN*R], the indices of the exponents of the polynomial.
//
//    Input, double PD[PN], the coefficient of each polynomial.  
//    For a Lagrange interpolant, this is the data value at each Lagrange point.
//
//    Input, int NI, the number of interpolant evaluation points.
//
//    Input, double XI[D*NI], the coordinates of the interpolation evaluation points.
//
//    Output, double INTERPOLANT_VALUE[NI], the value of the interpolant at XI.
//
{
  double *c;
  int *e;
  int i;
  int j;
  int oj;
  double *value;
  double *yi;

  yi = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    yi[i] = 0.0;
  }

  c = new double[r];
  e = new int[r];

  for ( j = 0; j < pn; j++ )
  {
    oj = po[j];
    for ( i = 0; i < oj; i++ )
    {
      c[i] = pc[j+i*pn];
      e[i] = pe[j+i*pn];
    }
    value = polynomial_value ( d, oj, c, e, ni, xi );
    for ( i = 0; i < ni; i++ )
    {
      yi[i] = yi[i] + pd[j] * value[i];
    }
    delete [] value;
  }

  delete [] c;
  delete [] e;

  return yi;
}
//****************************************************************************80

void lagrange_complete ( int d, int n, int r, int nd, double xd[], int po[], 
  double pc[], int pe[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_COMPLETE: Complete Lagrange polynomial basis from data.
//
//  Discussion:
//
//    This function represents algorithm 4.1 in the reference.
//
//    This function is given XD, a set of ND distinct data points in a 
//    D dimensional space, and returns information defining a set of 
//    ND Lagrange polynomials L(i)(X) with the property that:
//
//      L(i)(XD(j)) = delta(i,j)
//
//    In order for this computation to be carried out, it is necessary that
//    ND, the number of data points, is equal to R, the dimension of the 
//    space of polynomials in D dimensions and total degree N or less, that is:
//
//      ND = R = Choose ( N + D, N )
//
//    There will be ND polynomials returned.  Each polynomial can have
//    as many as R coefficients.
//
//    Each polynomial is given as a vector, with each entry corresponding
//    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
//
//      PO(i) is the order, that is, the number of nonzero coefficients;
//      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
//      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
//
//    The exponent codes are a compact way of recording the exponent vector
//    associated with each monomial.  If PE(i,j) = k, then the corresponding
//    vector of D exponents can be determined by:
//
//      E = mono_unrank_grlex ( D, k );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Tomas Sauer, Yuan Xu,
//    On multivariate Lagrange interpolation,
//    Mathematics of Computation,
//    Volume 64, Number 211, July 1995, pages 1147-1170.
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N, the maximum total degree.
//
//    Input, int R, the number of monomials in D dimensions 
//    of total degree N or less.
//
//    Input, int ND, the number of data points.
//    This function requires that the ND is equal to R.
//
//    Input, double XD[D*ND], the data points, which must be distinct.
//
//    Output, int PO[ND], the order (number of nonzero coefficients),
//    for the Lagrange basis polynomials.
//
//    Output, double PC[ND*R], the coefficients for the 
//    Lagrange basis polynomials.
//
//    Output, int PE[ND*R], the  exponent indices for the 
//    Lagrange basis polynomials.
//
{
  double *c;
  double *cj;
  double *ck;
  double d_max;
  double d_min;
  double d_tol;
  int *e;
  int *ej;
  int *ek;
  int i;
  int j;
  int k;
  int l;
  int o;;
  int oj;
  int ok;
  double *qc;
  int *qe;
  int *qo;
  double *value;
//
//  Verify that R is correct.
//
  if ( r != mono_upto_enum ( d, n ) )
  {
    cerr << "\n";
    cerr << "LAGRANGE_COMPLETE - Fatal error!\n";
    cerr << "  The value R is not correct.\n";
    exit ( 1 );
  }

  if ( r != nd )
  {
    cerr << "\n";
    cerr << "LAGRANGE_COMPLETE - Fatal error!\n";
    cerr << "  The value R = " << r << "\n";
    cerr << "  does not equal ND = " << nd << "\n";
    exit ( 1 );
  }
//
//  Verify that the points are sufficiently distinct.
//
  r8col_separation ( d, nd, xd, d_min, d_max );
  d_tol = sqrt ( r8_epsilon ( ) );

  if ( d_min < d_tol )
  {
    cerr << "\n";
    cerr << "LAGRANGE_COMPLETE - Fatal error!\n";
    cerr << "  Some points are too close!\n";
    cerr << "  Minimum data point separation is = " << d_min << "\n";
    exit ( 1 );
  }
//
//  Make some work space.
//
  c = new double[r];
  cj = new double[r];
  ck = new double[r];
  e = new int[r];
  ej = new int[r];
  ek = new int[r];
//
//  Initialize the polynomials Q, which span the space of
//  N-th degree polynomials.
//
//  Default option: 
//  * all ND-dimensional monomials of degree N or less.
//    in 2D, this might be 1, x, y, x^2, xy, y^2, ...
//
  qo = new int[r];
  qc = new double[r*r];
  qe = new int[r*r];

  for ( k = 0; k < r; k++ )
  {
    qo[k] = 1;
    qc[k+0*r] = 1.0;
    qe[k+0*r] = k + 1;
    for ( j = 1; j < r; j++ )
    {
      qc[k+j*r] = 0.0;
      qe[k+j*r] = 0;
    }
  }
//
//  Now set up the P polynomials.
//
  for ( k = 0; k < r; k++ )
  {
    po[k] = 0;
    for ( j = 0; j < r; j++ )
    {
      pc[k+j*r] = 0.0;
      pe[k+j*r] = 0;
    }
  }

  for ( k = 0; k < nd; k++ )
  {
//
//  Find the first polynomial Q(K:R)(X) which is nonzero at X(K).
//
    i = r + 1;

    for ( j = k; j < r; j++ )
    {
      o = qo[j];
      for ( l = 0; l < o; l++ )
      {
        c[l] = qc[j+l*r];
        e[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, o, c, e, 1, xd + k * d );

      if ( value[0] != 0.0 )
      {
        i = j;
        break;
      }
      delete [] value;
    }

    if ( i == r + 1 )
    {
      cerr << "\n";
      cerr << "LAGRANGE_COMPLETE - Fatal error!\n";
      cerr << "  I = R+1.\n";
      exit ( 1 );
    }
//
//  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
//
    o = qo[i];
    po[k] = qo[i];
    for ( l = 0; l < o; l++ )
    {
      pc[k+l*r] = qc[i+l*r] / value[0];
      pe[k+l*r] = qe[i+l*r];
    }

    delete [] value;
//
//  Modify P(1:k-1)(X).
//
    for ( j = 0; j < k; j++ )
    {
      oj = po[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = pc[j+l*r];
        ej[l] = pe[j+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*r];
        ek[l] = pe[k+l*r];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      po[j] = o;
      for ( l = 0; l < o; l++ )
      {
        pc[j+l*r] = c[l];
        pe[j+l*r] = e[l];
      }

      delete [] value;
    }
//
//  Modify Q(I:downto:K+1)
//
    for ( j = i; k < j; j-- )
    {
      oj = qo[j-1];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j-1+l*r];
        ej[l] = qe[j-1+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );
 
      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*r];
        ek[l] = pe[k+l*r];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
//
//  Modify Q(I+1:R)
//
    for ( j = i + 1; j < r; j++ )
    {
      oj = qo[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j+l*r];
        ej[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*r];
        ek[l] = pe[k+l*r];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
  }
//
//  Get rid of tiny coefficients.
//
  for ( i = 0; i < nd; i++ )
  {
    oj = po[i];
    for ( l = 0; l < oj; l++ )
    {
      cj[l] = pc[i+l*r];
      ej[l] = pe[i+l*r];
    }

    polynomial_compress ( oj, cj, ej, ok, ck, ek );

    po[i] = ok;
    for ( l = 0; l < ok; l++ )
    {
      pc[i+l*r] = ck[l];
      pe[i+l*r] = ek[l];
    }
  }
//
//  Free memory.
//
  delete [] c;
  delete [] cj;
  delete [] ck;
  delete [] e;
  delete [] ej;
  delete [] ek;
  delete [] qc;
  delete [] qe;
  delete [] qo;

  return;
}
//****************************************************************************80

void lagrange_complete2 ( int d, int n, int r, int nd, double xd[], int po[], 
  double pc[], int pe[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_COMPLETE2: Complete Lagrange polynomial basis from data.
//
//  Discussion:
//
//    This function represents algorithm 4.1 in the reference,
//    with the further modification that a form of "pivoting" is used
//    to select the next polynomial as the one with maximum absolute
//    value at the current node.
//
//    This function is given XD, a set of ND distinct data points in a 
//    D dimensional space, and returns information defining a set of 
//    ND Lagrange polynomials L(i)(X) with the property that:
//
//      L(i)(XD(j)) = delta(i,j)
//
//    In order for this computation to be carried out, it is necessary that
//    ND, the number of data points, is equal to R, the dimension of the 
//    space of polynomials in D dimensions and total degree N or less, that is:
//
//      ND = R = Choose ( N + D, N )
//
//    There will be ND polynomials returned.  Each polynomial can have
//    as many as R coefficients.
//
//    Each polynomial is given as a vector, with each entry corresponding
//    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
//
//      PO(i) is the order, that is, the number of nonzero coefficients;
//      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
//      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
//
//    The exponent codes are a compact way of recording the exponent vector
//    associated with each monomial.  If PE(i,j) = k, then the corresponding
//    vector of D exponents can be determined by:
//
//      E = mono_unrank_grlex ( D, k );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Tomas Sauer, Yuan Xu,
//    On multivariate Lagrange interpolation,
//    Mathematics of Computation,
//    Volume 64, Number 211, July 1995, pages 1147-1170.
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N, the maximum total degree.
//
//    Input, int R, the number of monomials in D dimensions 
//    of total degree N or less.
//
//    Input, int ND, the number of data points.
//    This function requires that the ND is equal to R.
//
//    Input, double XD[D*ND], the data points, which must be distinct.
//
//    Output, int PO[ND], the order (number of nonzero coefficients),
//    for the Lagrange basis polynomials.
//
//    Output, double PC[ND*R], the coefficients for the 
//    Lagrange basis polynomials.
//
//    Output, int PE[ND*R], the  exponent indices for the 
//    Lagrange basis polynomials.
//
{
  double *c;
  double *cj;
  double *ck;
  double d_max;
  double d_min;
  double d_tol;
  int *e;
  int *ej;
  int *ek;
  int i;
  int j;
  int k;
  int l;
  int o;;
  int oj;
  int ok;
  double *qc;
  int *qe;
  int *qo;
  double *value;
  double value_max;
//
//  Verify that R is correct.
//
  if ( r != mono_upto_enum ( d, n ) )
  {
    cerr << "\n";
    cerr << "LAGRANGE_COMPLETE2 - Fatal error!\n";
    cerr << "  The value R is not correct.\n";
    exit ( 1 );
  }

  if ( r != nd )
  {
    cerr << "\n";
    cerr << "LAGRANGE_COMPLETE2 - Fatal error!\n";
    cerr << "  The value R = " << r << "\n";
    cerr << "  does not equal ND = " << nd << "\n";
    exit ( 1 );
  }
//
//  Verify that the points are sufficiently distinct.
//
  r8col_separation ( d, nd, xd, d_min, d_max );
  d_tol = sqrt ( r8_epsilon ( ) );

  if ( d_min < d_tol )
  {
    cerr << "\n";
    cerr << "LAGRANGE_COMPLETE2 - Fatal error!\n";
    cerr << "  Some points are too close!\n";
    cerr << "  Minimum data point separation is = " << d_min << "\n";
    exit ( 1 );
  }
//
//  Make some work space.
//
  c = new double[r];
  cj = new double[r];
  ck = new double[r];
  e = new int[r];
  ej = new int[r];
  ek = new int[r];
//
//  Initialize the polynomials Q, which span the space of
//  N-th degree polynomials.
//
//  Default option: 
//  * all ND-dimensional monomials of degree N or less.
//    in 2D, this might be 1, x, y, x^2, xy, y^2, ...
//
  qo = new int[r];
  qc = new double[r*r];
  qe = new int[r*r];

  for ( k = 0; k < r; k++ )
  {
    qo[k] = 1;
    qc[k+0*r] = 1.0;
    qe[k+0*r] = k + 1;
    for ( j = 1; j < r; j++ )
    {
      qc[k+j*r] = 0.0;
      qe[k+j*r] = 0;
    }
  }
//
//  Now set up the P polynomials.
//
  for ( k = 0; k < r; k++ )
  {
    po[k] = 0;
    for ( j = 0; j < r; j++ )
    {
      pc[k+j*r] = 0.0;
      pe[k+j*r] = 0;
    }
  }

  for ( k = 0; k < nd; k++ )
  {
//
//  Find the first polynomial Q(K:R)(X) which is nonzero at X(K).
//
    i = r + 1;
    value_max = 0.0;

    for ( j = k; j < r; j++ )
    {
      o = qo[j];
      for ( l = 0; l < o; l++ )
      {
        c[l] = qc[j+l*r];
        e[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, o, c, e, 1, xd + k * d );

      if ( fabs ( value_max ) <= fabs ( value[0] ) )
      {
        i = j;
        value_max = value[0];
      }

      delete [] value;
    }

    if ( i == r + 1 )
    {
      cerr << "\n";
      cerr << "LAGRANGE_COMPLETE2 - Fatal error!\n";
      cerr << "  I = R+1.\n";
      exit ( 1 );
    }

    value[0] = value_max;
//
//  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
//
    o = qo[i];
    po[k] = qo[i];
    for ( l = 0; l < o; l++ )
    {
      pc[k+l*r] = qc[i+l*r] / value[0];
      pe[k+l*r] = qe[i+l*r];
    }

    delete [] value;
//
//  Modify P(1:k-1)(X).
//
    for ( j = 0; j < k; j++ )
    {
      oj = po[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = pc[j+l*r];
        ej[l] = pe[j+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*r];
        ek[l] = pe[k+l*r];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      po[j] = o;
      for ( l = 0; l < o; l++ )
      {
        pc[j+l*r] = c[l];
        pe[j+l*r] = e[l];
      }

      delete [] value;
    }
//
//  Modify Q(I:downto:K+1)
//
    for ( j = i; k < j; j-- )
    {
      oj = qo[j-1];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j-1+l*r];
        ej[l] = qe[j-1+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );
 
      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*r];
        ek[l] = pe[k+l*r];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
//
//  Modify Q(I+1:R)
//
    for ( j = i + 1; j < r; j++ )
    {
      oj = qo[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j+l*r];
        ej[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*r];
        ek[l] = pe[k+l*r];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
  }
//
//  Get rid of tiny coefficients.
//
  for ( i = 0; i < nd; i++ )
  {
    oj = po[i];
    for ( l = 0; l < oj; l++ )
    {
      cj[l] = pc[i+l*r];
      ej[l] = pe[i+l*r];
    }

    polynomial_compress ( oj, cj, ej, ok, ck, ek );

    po[i] = ok;
    for ( l = 0; l < ok; l++ )
    {
      pc[i+l*r] = ck[l];
      pe[i+l*r] = ek[l];
    }
  }
//
//  Free memory.
//
  delete [] c;
  delete [] cj;
  delete [] ck;
  delete [] e;
  delete [] ej;
  delete [] ek;
  delete [] qc;
  delete [] qe;
  delete [] qo;

  return;
}
//****************************************************************************80

void lagrange_partial ( int d, int n, int r, int nd, double xd[], int po[], 
  double pc[], int pe[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_PARTIAL: Partial Lagrange polynomial basis from data.
//
//  Discussion:
//
//    This function represents algorithm 4.1 in the reference,
//    modified for the case where the number of data points is less
//    than the dimension of the desired polynomial space.
//
//    This function is given XD, a set of ND distinct data points in a 
//    D dimensional space, and returns information defining a set of 
//    ND Lagrange polynomials L(i)(X) with the property that:
//
//      L(i)(XD(j)) = delta(i,j)
//
//    This function is used in cases where ND, the number of data points, 
//    is less than or equal to R, the dimension of the space of polynomials 
//    in D dimensions and total degree N or less, that is:
//
//      ND <= R = Choose ( N + D, N )
//
//    There will be ND polynomials returned.  Each polynomial can have
//    as many as R coefficients.
//
//    Each polynomial is given as a vector, with each entry corresponding
//    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
//
//      PO(i) is the order, that is, the number of nonzero coefficients;
//      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
//      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
//
//    The exponent codes are a compact way of recording the exponent vector
//    associated with each monomial.  If PE(i,j) = k, then the corresponding
//    vector of D exponents can be determined by:
//
//      E = mono_unrank_grlex ( D, k );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Tomas Sauer, Yuan Xu,
//    On multivariate Lagrange interpolation,
//    Mathematics of Computation,
//    Volume 64, Number 211, July 1995, pages 1147-1170.
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N, the maximum total degree.
//
//    Input, int R, the number of monomials in D dimensions 
//    of total degree N or less.
//
//    Input, int ND, the number of data points.
//    It must be the case that ND <= R.
//
//    Input, double XD[D*ND], the data points, which must be distinct.
//
//    Output, int PO[ND], the order (number of nonzero coefficients),
//    for the Lagrange basis polynomials.
//
//    Output, double PC[ND*R], the coefficients for the 
//    Lagrange basis polynomials.
//
//    Output, int PE[ND*R], the  exponent indices for the 
//    Lagrange basis polynomials.
//
{
  double *c;
  double *cj;
  double *ck;
  double d_max;
  double d_min;
  double d_tol;
  int *e;
  int *ej;
  int *ek;
  int i;
  int j;
  int k;
  int l;
  int o;
  int oj;
  int ok;
  double *qc;
  int *qe;
  int *qo;
  double *value;
//
//  Verify that R is correct.
//
  if ( r != mono_upto_enum ( d, n ) )
  {
    cerr << "\n";
    cerr << "LAGRANGE_PARTIAL - Fatal error!\n";
    cerr << "  The value R is not correct.\n";
    exit ( 1 );
  }

  if ( r < nd )
  {
    cerr << "\n";
    cerr << "LAGRANGE_PARTIAL - Fatal error!\n";
    cerr << "  The value R = " << r << "\n";
    cerr << "  is less than ND = " << nd << "\n";
    exit ( 1 );
  }
//
//  Verify that the points are sufficiently distinct.
//
  r8col_separation ( d, nd, xd, d_min, d_max );
  d_tol = sqrt ( r8_epsilon ( ) );

  if ( d_min < d_tol )
  {
    cerr << "\n";
    cerr << "LAGRANGE_PARTIAL - Fatal error!\n";
    cerr << "  Some points are too close!\n";
    cerr << "  Minimum data point separation is = " << d_min << "\n";
    exit ( 1 );
  }
//
//  Make some work space.
//
  c = new double[r];
  cj = new double[r];
  ck = new double[r];
  e = new int[r];
  ej = new int[r];
  ek = new int[r];
//
//  Initialize the polynomials Q, which span the space of
//  N-th degree polynomials.
//
//  Default option: 
//  * all ND-dimensional monomials of degree N or less.
//    in 2D, this might be 1, x, y, x^2, xy, y^2, ...
//
  qo = new int[r];
  qc = new double[r*r];
  qe = new int[r*r];

  for ( k = 0; k < r; k++ )
  {
    qo[k] = 1;
    qc[k+0*r] = 1.0;
    qe[k+0*r] = k + 1;
    for ( j = 1; j < r; j++ )
    {
      qc[k+j*r] = 0.0;
      qe[k+j*r] = 0;
    }
  }
//
//  Now set up the P polynomials.
//
  for ( k = 0; k < nd; k++ )
  {
    po[k] = 0;
    for ( j = 0; j < r; j++ )
    {
      pc[k+j*nd] = 0.0;
      pe[k+j*nd] = 0;
    }
  }
  for ( k = 0; k < nd; k++ )
  {
//
//  Find the first polynomial Q(K:R)(X) which is nonzero at X(K).
//
    i = r + 1;

    for ( j = k; j < r; j++ )
    {
      o = qo[j];
      for ( l = 0; l < o; l++ )
      {
        c[l] = qc[j+l*r];
        e[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, o, c, e, 1, xd + k * d );

      if ( value[0] != 0.0 )
      {
        i = j;
        break;
      }
      delete [] value;
    }

    if ( i == r + 1 )
    {
      cerr << "\n";
      cerr << "LAGRANGE_PARTIAL - Fatal error!\n";
      cerr << "  I = R+1.\n";
      exit ( 1 );
    }
//
//  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
//
    o = qo[i];
    po[k] = qo[i];
    for ( l = 0; l < o; l++ )
    {
      pc[k+l*nd] = qc[i+l*r] / value[0];
      pe[k+l*nd] = qe[i+l*r];
    }

    delete [] value;
//
//  Modify P(1:k-1)(X).
//
    for ( j = 0; j < k; j++ )
    {
      oj = po[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = pc[j+l*nd];
        ej[l] = pe[j+l*nd];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*nd];
        ek[l] = pe[k+l*nd];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      po[j] = o;
      for ( l = 0; l < o; l++ )
      {
        pc[j+l*nd] = c[l];
        pe[j+l*nd] = e[l];
      }

      delete [] value;
    }
//
//  Modify Q(I:downto:K+1)
//
    for ( j = i; k < j; j-- )
    {
      oj = qo[j-1];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j-1+l*r];
        ej[l] = qe[j-1+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );
 
      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*nd];
        ek[l] = pe[k+l*nd];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
//
//  Modify Q(I+1:R)
//
    for ( j = i + 1; j < r; j++ )
    {
      oj = qo[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j+l*r];
        ej[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*nd];
        ek[l] = pe[k+l*nd];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
  }
//
//  Get rid of tiny coefficients.
//
  for ( i = 0; i < nd; i++ )
  {
    oj = po[i];
    for ( l = 0; l < oj; l++ )
    {
      cj[l] = pc[i+l*nd];
      ej[l] = pe[i+l*nd];
    }

    polynomial_compress ( oj, cj, ej, ok, ck, ek );

    po[i] = ok;
    for ( l = 0; l < ok; l++ )
    {
      pc[i+l*nd] = ck[l];
      pe[i+l*nd] = ek[l];
    }
  }
//
//  Free memory.
//
  delete [] c;
  delete [] cj;
  delete [] ck;
  delete [] e;
  delete [] ej;
  delete [] ek;
  delete [] qc;
  delete [] qe;
  delete [] qo;

  return;
}
//****************************************************************************80

void lagrange_partial2 ( int d, int n, int r, int nd, double xd[], int po[], 
  double pc[], int pe[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_PARTIAL2: Partial Lagrange polynomial basis from data.
//
//  Discussion:
//
//    This function represents algorithm 4.1 in the reference,
//    modified for the case where the number of data points is less
//    than the dimension of the desired polynomial space,
//    with the further modification that a form of "pivoting" is used
//    to select the next polynomial as the one with maximum absolute
//    value at the current node.
//
//    This function is given XD, a set of ND distinct data points in a 
//    D dimensional space, and returns information defining a set of 
//    ND Lagrange polynomials L(i)(X) with the property that:
//
//      L(i)(XD(j)) = delta(i,j)
//
//    This function is used in cases where ND, the number of data points, 
//    is less than or equal to R, the dimension of the space of polynomials 
//    in D dimensions and total degree N or less, that is:
//
//      ND <= R = Choose ( N + D, N )
//
//    There will be ND polynomials returned.  Each polynomial can have
//    as many as R coefficients.
//
//    Each polynomial is given as a vector, with each entry corresponding
//    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
//
//      PO(i) is the order, that is, the number of nonzero coefficients;
//      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
//      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
//
//    The exponent codes are a compact way of recording the exponent vector
//    associated with each monomial.  If PE(i,j) = k, then the corresponding
//    vector of D exponents can be determined by:
//
//      E = mono_unrank_grlex ( D, k );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Tomas Sauer, Yuan Xu,
//    On multivariate Lagrange interpolation,
//    Mathematics of Computation,
//    Volume 64, Number 211, July 1995, pages 1147-1170.
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N, the maximum total degree.
//
//    Input, int R, the number of monomials in D dimensions 
//    of total degree N or less.
//
//    Input, int ND, the number of data points.
//    It must be the case that ND <= R.
//
//    Input, double XD[D*ND], the data points, which must be distinct.
//
//    Output, int PO[ND], the order (number of nonzero coefficients),
//    for the Lagrange basis polynomials.
//
//    Output, double PC[ND*R], the coefficients for the 
//    Lagrange basis polynomials.
//
//    Output, int PE[ND*R], the  exponent indices for the 
//    Lagrange basis polynomials.
//
{
  double *c;
  double *cj;
  double *ck;
  double d_max;
  double d_min;
  double d_tol;
  int *e;
  int *ej;
  int *ek;
  int i;
  int j;
  int k;
  int l;
  int o;
  int oj;
  int ok;
  double *qc;
  int *qe;
  int *qo;
  double *value;
  double value_max;
//
//  Verify that R is correct.
//
  if ( r != mono_upto_enum ( d, n ) )
  {
    cerr << "\n";
    cerr << "LAGRANGE_PARTIAL2 - Fatal error!\n";
    cerr << "  The value R is not correct.\n";
    exit ( 1 );
  }

  if ( r < nd )
  {
    cerr << "\n";
    cerr << "LAGRANGE_PARTIAL2 - Fatal error!\n";
    cerr << "  The value R = " << r << "\n";
    cerr << "  is less than ND = " << nd << "\n";
    exit ( 1 );
  }
//
//  Verify that the points are sufficiently distinct.
//
  r8col_separation ( d, nd, xd, d_min, d_max );
  d_tol = sqrt ( r8_epsilon ( ) );

  if ( d_min < d_tol )
  {
    cerr << "\n";
    cerr << "LAGRANGE_PARTIAL2 - Fatal error!\n";
    cerr << "  Some points are too close!\n";
    cerr << "  Minimum data point separation is = " << d_min << "\n";
    exit ( 1 );
  }
//
//  Make some work space.
//
  c = new double[r];
  cj = new double[r];
  ck = new double[r];
  e = new int[r];
  ej = new int[r];
  ek = new int[r];
//
//  Initialize the polynomials Q, which span the space of
//  N-th degree polynomials.
//
//  Default option: 
//  * all ND-dimensional monomials of degree N or less.
//    in 2D, this might be 1, x, y, x^2, xy, y^2, ...
//
  qo = new int[r];
  qc = new double[r*r];
  qe = new int[r*r];

  for ( k = 0; k < r; k++ )
  {
    qo[k] = 1;
    qc[k+0*r] = 1.0;
    qe[k+0*r] = k + 1;
    for ( j = 1; j < r; j++ )
    {
      qc[k+j*r] = 0.0;
      qe[k+j*r] = 0;
    }
  }
//
//  Now set up the P polynomials.
//
  for ( k = 0; k < nd; k++ )
  {
    po[k] = 0;
    for ( j = 0; j < r; j++ )
    {
      pc[k+j*nd] = 0.0;
      pe[k+j*nd] = 0;
    }
  }

  for ( k = 0; k < nd; k++ )
  {
//
//  Find the first polynomial Q(K:R)(X) which is nonzero at X(K).
//
    i = r + 1;
    value_max = 0.0;

    for ( j = k; j < r; j++ )
    {
      o = qo[j];
      for ( l = 0; l < o; l++ )
      {
        c[l] = qc[j+l*r];
        e[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, o, c, e, 1, xd + k * d );

      if ( fabs ( value_max ) <= fabs ( value[0] ) )
      {
        i = j;
        value_max = value[0];
      }

      delete [] value;
    }

    if ( i == r + 1 )
    {
      cerr << "\n";
      cerr << "LAGRANGE_PARTIAL2 - Fatal error!\n";
      cerr << "  I = R+1.\n";
      exit ( 1 );
    }

    value = new double[1];
    value[0] = value_max;
//
//  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
//
    o = qo[i];
    po[k] = qo[i];
    for ( l = 0; l < o; l++ )
    {
      pc[k+l*nd] = qc[i+l*r] / value[0];
      pe[k+l*nd] = qe[i+l*r];
    }

    delete [] value;
//
//  Modify P(1:k-1)(X).
//
    for ( j = 0; j < k; j++ )
    {
      oj = po[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = pc[j+l*nd];
        ej[l] = pe[j+l*nd];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*nd];
        ek[l] = pe[k+l*nd];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      po[j] = o;
      for ( l = 0; l < o; l++ )
      {
        pc[j+l*nd] = c[l];
        pe[j+l*nd] = e[l];
      }

      delete [] value;
    }
//
//  Modify Q(I:downto:K+1)
//
    for ( j = i; k < j; j-- )
    {
      oj = qo[j-1];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j-1+l*r];
        ej[l] = qe[j-1+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );
 
      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*nd];
        ek[l] = pe[k+l*nd];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
//
//  Modify Q(I+1:R)
//
    for ( j = i + 1; j < r; j++ )
    {
      oj = qo[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j+l*r];
        ej[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*nd];
        ek[l] = pe[k+l*nd];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
  }
//
//  Get rid of tiny coefficients.
//
  for ( i = 0; i < nd; i++ )
  {
    oj = po[i];
    for ( l = 0; l < oj; l++ )
    {
      cj[l] = pc[i+l*nd];
      ej[l] = pe[i+l*nd];
    }

    polynomial_compress ( oj, cj, ej, ok, ck, ek );

    po[i] = ok;
    for ( l = 0; l < ok; l++ )
    {
      pc[i+l*nd] = ck[l];
      pe[i+l*nd] = ek[l];
    }
  }
//
//  Free memory.
//
  delete [] c;
  delete [] cj;
  delete [] ck;
  delete [] e;
  delete [] ej;
  delete [] ek;
  delete [] qc;
  delete [] qe;
  delete [] qo;

  return;
}
//****************************************************************************80

void lagrange_partial3 ( int d, int n, int nd, double xd[], int option, 
  int po[], double **pc, int **pe, int &n2 )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_PARTIAL3: Partial Lagrange polynomial basis from data.
//
//  Discussion:
//
//    This function, together with lagrange_partial4(), is a representation
//    of algorithm 4.1 in the reference, modified:
//    * for the case where the number of data points is less
//      than the dimension of the desired polynomial space,
//    * so that a form of "pivoting" is used
//      to select the next polynomial as the one with maximum absolute
//      value at the current node;
//    * so that if the problem is not well posed, successively higher
//      values of N are tried.
//
//    This function is given XD, a set of ND distinct data points in a 
//    D dimensional space, and returns information defining a set of 
//    ND Lagrange polynomials L(i)(X) with the property that:
//
//      L(i)(XD(j)) = delta(i,j)
//
//    This function is used in cases where ND, the number of data points, 
//    is less than or equal to R, the dimension of the space of polynomials 
//    in D dimensions and total degree N or less, that is:
//
//      ND <= R = Choose ( N + D, N )
//
//    There will be ND polynomials returned.  Each polynomial can have
//    as many as R coefficients.
//
//    Each polynomial is given as a vector, with each entry corresponding
//    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
//
//      PO(i) is the order, that is, the number of nonzero coefficients;
//      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
//      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
//
//    The exponent codes are a compact way of recording the exponent vector
//    associated with each monomial.  If PE(i,j) = k, then the corresponding
//    vector of D exponents can be determined by:
//
//      E = mono_unrank_grlex ( D, k );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Tomas Sauer, Yuan Xu,
//    On multivariate Lagrange interpolation,
//    Mathematics of Computation,
//    Volume 64, Number 211, July 1995, pages 1147-1170.
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N, the maximum total degree.
//
//    Input, int ND, the number of data points.
//    It must be the case that ND <= R = the number of monomials 
//    of degree N in D dimensions.
//
//    Input, double XD[D*ND], the data points, which must be distinct.
//
//    Input, int OPTION, determines the initial basis:
//    0, monomials, 1, x, y, x^2, xy, y^2, x^3, ...
//    1, Legendre products, 1, y, x, (3y^2-1)/2, xy, (3^x^2-1), (5y^3-3)/2, ...
//
//    Output, int PO[ND], the order (number of nonzero coefficients) for the 
//    Lagrange basis polynomials.
//
//    Output, double **PC, the ND by R array of coefficients for the 
//    Lagrange basis polynomials.
//
//    Output, int **PE, the ND by R array of exponent indices for the 
//    Lagrange basis polynomials.
//
//    Output, int &N2, the adjusted value of N, which may have been
//    increased because the interpolation problem for N was not well posed.
//
{
  double d_max;
  double d_min;
  double d_tol;
  int j;
  int r;
  bool success;
  double tol;
//
//  Verify that the points are sufficiently distinct.
//
  r8col_separation ( d, nd, xd, d_min, d_max );
  d_tol = sqrt ( r8_epsilon ( ) );

  if ( d_min < d_tol )
  {
    cerr << "\n";
    cerr << "LAGRANGE_PARTIAL3 - Fatal error!\n";
    cerr << "  Some points are too close!\n";
    cerr << "  Minimum data point separation is = " << d_min << "\n";
    exit ( 1 );
  }
//
//  Search for the appropriate interpolation space.
//
  n2 = n;
  tol = 0.0001;

  for ( ; ; )
  {
    r = mono_upto_enum ( d, n2 );

    ( *pc ) = new double[nd*r];
    ( *pe ) = new int[nd*r]; 

    success = lagrange_partial4 ( d, n2, r, nd, xd, option, tol, po, *pc, *pe );

    if ( success )
    {
      return;
    }

    delete [] *pc;
    delete [] *pe;

    n2 = n2 + 1;
    cout << "LAGRANGE_PARTIAL3 - Increase N to " << n2 << "\n";
  }

  return;
}
//****************************************************************************80

bool lagrange_partial4 ( int d, int n, int r, int nd, double xd[], int option, 
  double tol, int po[], double pc[], int pe[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_PARTIAL4: Partial Lagrange polynomial basis from data.
//
//  Discussion:
//
//    This function, together with lagrange_partial3(), is a representation
//    of algorithm 4.1 in the reference, modified:
//    * for the case where the number of data points is less
//      than the dimension of the desired polynomial space,
//    * so that a form of "pivoting" is used
//      to select the next polynomial as the one with maximum absolute
//      value at the current node;
//    * so that if the problem is not well posed, successively higher
//      values of N are tried.
//
//    This function is given XD, a set of ND data points in a D dimensional
//    space, and returns information defining a set of ND Lagrange polynomials
//    L(i)(X) with the property that:
//
//      L(i)(XD(j)) = delta(i,j)
//
//    This function is used in cases where ND, the number of data points, 
//    is less than or equal to R, the dimension of the space of polynomials 
//    in D dimensions and total degree N or less, that is:
//
//      ND <= R = Choose ( N + D, N )
//
//    There will be ND polynomials returned.  Each polynomial can have
//    as many as R coefficients.
//
//    Each polynomial is given as a vector, with each entry corresponding
//    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
//
//      PO(i) is the order, that is, the number of nonzero coefficients;
//      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
//      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
//
//    The exponent codes are a compact way of recording the exponent vector
//    associated with each monomial.  If PE(i,j) = k, then the corresponding
//    vector of D exponents can be determined by:
//
//      E = mono_unrank_grlex ( D, k );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Tomas Sauer, Yuan Xu,
//    On multivariate Lagrange interpolation,
//    Mathematics of Computation,
//    Volume 64, Number 211, July 1995, pages 1147-1170.
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N, the maximum total degree.
//
//    Input, int R, the number of monomials in D dimensions 
//    of total degree N or less.
//
//    Input, int ND, the number of data points.
//    It must be the case that ND <= R.
//
//    Input, double XD[D*ND], the data points.
//
//    Input, int OPTION, determines the initial basis:
//    0, monomials, 1, x, y, x^2, xy, y^2, x^3, ...
//    1, Legendre products, 1, y, x, (3y^2-1)/2, xy, (3^x^2-1), (5y^3-3)/2, ...
//
//    Input, double TOL, a tolerance for the pivoting operation.
//    If no unused polynomial can be found with a value at least TOL
//    at the current point, the algorithm fails.
//
//    Output, int PO[ND], the order (number of nonzero coefficients) for the 
//    Lagrange basis polynomials.
//
//    Output, double PC[ND*R], the coefficients for the 
//    Lagrange basis polynomials.
//
//    Output, int PE[ND*R], the exponent indices for the 
//    Lagrange basis polynomials.
//
//    Output, bool LAGRANGE_PARTIAL4, is 0 if the algorithm failed
//    (in which case the other outputs are not useful),
//    and 1 if it was successful.
//
{
  double *c;
  double *cj;
  double *ck;
  int *e;
  int *ej;
  int *ek;
  int i;
  int j;
  int k;
  int l;
  int *lpp;
  int o;
  int oj;
  int ok;
  double *qc;
  int *qe;
  int *qo;
  bool success;
  char title[80];
  double *value;
  double value_max;

  success = true;
//
//  Verify that R is acceptable.
//
  if ( r < nd )
  {
    cerr << "\n";
    cerr << "LAGRANGE_PARTIAL4 - Fatal error!\n";
    cerr << "  The value R = " << r << "\n";
    cerr << "  is less than ND = " << nd << "\n";
    exit ( 1 );
  }
//
//  Make some work space.
//
  c = new double[r];
  cj = new double[r];
  ck = new double[r];
  e = new int[r];
  ej = new int[r];
  ek = new int[r];
//
//  Initialize the polynomials Q spanning the space of N-th degree polynomials.
//
  qo = new int[r];
  qc = new double[r*r];
  qe = new int[r*r];

  for ( j = 0; j < r; j++ )
  {
    qo[j] = 0;
    for ( i = 0; i < r; i++ )
    {
      qc[i+j*r] = 0.0;
      qe[i+j*r] = 0;
    }
  }
//
//  Option 0: First R D-dimensional monomials
//  Option 1: First R D-dimensional Legendre product polynomials.
//
  for ( k = 0; k < r; k++ )
  {
    if ( option == 0 )
    {
      o = 1;
      c[0] = 1.0;
      e[0] = k + 1;
    }
    else if ( option == 1 )
    {
      lpp = comp_unrank_grlex ( d, k + 1 );
      lpp_to_polynomial ( d, lpp, r, o, c, e );
      delete [] lpp;
    }

    qo[k] = o;
    for ( j = 0; j < o; j++ )
    {
      qc[k+j*r] = c[j];
      qe[k+j*r] = e[j];
    }
  }
//
//  Now set up the P polynomials.
//
  for ( k = 0; k < nd; k++ )
  {
    po[k] = 0;
    for ( j = 0; j < r; j++ )
    {
      pc[k+j*nd] = 0.0;
      pe[k+j*nd] = 0;
    }
  }

  for ( k = 0; k < nd; k++ )
  {
//
//  Find the polynomial Q(K:R)(X) which is most nonzero at X(K).
//
    i = r + 1;
    value_max = 0.0;

    for ( j = k; j < r; j++ )
    {
      o = qo[j];
      for ( l = 0; l < o; l++ )
      {
        c[l] = qc[j+l*r];
        e[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, o, c, e, 1, xd + k * d );

      if ( fabs ( value_max ) <= fabs ( value[0] ) )
      {
        i = j;
        value_max = value[0];
      }

      delete [] value;
    }
//
//  If the best nonzero value was too small or zero, fail.
//
    if ( fabs ( value_max ) < tol || i == r + 1 )
    {
      success = false;
      cout << "LAGRANGE_PARTIAL4 - Unacceptable VALUE_MAX = " << value_max << "\n";
      return success;
    }

    value = new double[1];
    value[0] = value_max;
//
//  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
//
    o = qo[i];
    po[k] = qo[i];
    for ( l = 0; l < o; l++ )
    {
      pc[k+l*nd] = qc[i+l*r] / value[0];
      pe[k+l*nd] = qe[i+l*r];
    }

    delete [] value;
//
//  Modify P(1:k-1)(X).
//
    for ( j = 0; j < k; j++ )
    {
      oj = po[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = pc[j+l*nd];
        ej[l] = pe[j+l*nd];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*nd];
        ek[l] = pe[k+l*nd];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      po[j] = o;
      for ( l = 0; l < o; l++ )
      {
        pc[j+l*nd] = c[l];
        pe[j+l*nd] = e[l];
      }

      delete [] value;
    }
//
//  Modify Q(I:downto:K+1)
//
    for ( j = i; k < j; j-- )
    {
      oj = qo[j-1];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j-1+l*r];
        ej[l] = qe[j-1+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );
 
      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*nd];
        ek[l] = pe[k+l*nd];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
//
//  Modify Q(I+1:R)
//
    for ( j = i + 1; j < r; j++ )
    {
      oj = qo[j];
      for ( l = 0; l < oj; l++ )
      {
        cj[l] = qc[j+l*r];
        ej[l] = qe[j+l*r];
      }

      value = polynomial_value ( d, oj, cj, ej, 1, xd + k * d );

      ok = po[k];
      for ( l = 0; l < ok; l++ )
      {
        ck[l] = pc[k+l*nd];
        ek[l] = pe[k+l*nd];
      }

      polynomial_axpy ( - value[0], ok, ck, ek, oj, cj, ej, o, c, e );

      delete [] value;

      qo[j] = o;
      for ( l = 0; l < o; l++ )
      {
        qc[j+l*r] = c[l];
        qe[j+l*r] = e[l];
      }
    }
  }
//
//  Get rid of tiny coefficients.
//
  for ( i = 0; i < nd; i++ )
  {
    oj = po[i];
    for ( l = 0; l < oj; l++ )
    {
      cj[l] = pc[i+l*nd];
      ej[l] = pe[i+l*nd];
    }

    polynomial_compress ( oj, cj, ej, ok, ck, ek );
    
    po[i] = ok;
    for ( l = 0; l < ok; l++ )
    {
      pc[i+l*nd] = ck[l];
      pe[i+l*nd] = ek[l];
    }
  }
//
//  Free memory.
//
  delete [] c;
  delete [] cj;
  delete [] ck;
  delete [] e;
  delete [] ej;
  delete [] ek;
  delete [] qc;
  delete [] qe;
  delete [] qo;

  return success;
}
//****************************************************************************80

void lp_coefficients ( int n, int &o, double c[], int f[] )

//****************************************************************************80
//
//  Purpose:
//
//    LP_COEFFICIENTS: coefficients of Legendre polynomials P(n,x).
//
//  First terms:
//
//     1
//     0     1
//    -1/2   0      3/2
//     0    -3/2    0     5/2
//     3/8   0    -30/8   0     35/8
//     0    15/8    0   -70/8    0     63/8
//    -5/16  0    105/16  0   -315/16   0    231/16
//     0   -35/16   0   315/16   0   -693/16   0    429/16
//
//     1.00000
//     0.00000  1.00000
//    -0.50000  0.00000  1.50000
//     0.00000 -1.50000  0.00000  2.5000
//     0.37500  0.00000 -3.75000  0.00000  4.37500
//     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
//    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
//     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Output, int &O, the number of coefficients.
//
//    Output, double C[(N+2)/2], the coefficients of the Legendre
//    polynomial of degree N.
//
//    Output, int F[(N+2)/2], the exponents.
//
{
  double *ctable;
  int i;
  int j;
  int k;

  ctable = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      ctable[i+j*(n+1)] = 0.0;
    }
  }

  ctable[0+0*(n+1)] = 1.0;

  if ( 0 < n )
  {
    ctable[1+1*(n+1)] = 1.0;

    for ( i = 2; i <= n; i++ )
    {
      for ( j = 0; j <= i-2; j++ )
      {
        ctable[i+j*(n+1)] =
            ( double ) (   - i + 1 ) * ctable[i-2+j*(n+1)] / ( double ) i;
      }
      for ( j = 1; j <= i; j++ )
      {
        ctable[i+j*(n+1)] = ctable[i+j*(n+1)]
          + ( double ) ( i + i - 1 ) * ctable[i-1+(j-1)*(n+1)] / ( double ) i;
      }
    }
  }
//
//  Extract the nonzero data from the alternating columns of the last row.
//
  o = ( n + 2 ) / 2;

  k = o;
  for ( j = n; 0 <= j; j = j - 2 )
  {
    k = k - 1;
    c[k] = ctable[n+j*(n+1)];
    f[k] = j;
  }

  delete [] ctable;

  return;
}
//****************************************************************************80

void lpp_to_polynomial ( int m, int l[], int o_max, int &o, double c[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    LPP_TO_POLYNOMIAL writes a Legendre Product Polynomial as a polynomial.
//
//  Discussion:
//
//    For example, if 
//      M = 3,
//      L = ( 1, 0, 2 ),
//    then
//      L(1,0,2)(X,Y,Z) 
//      = L(1)(X) * L(0)(Y) * L(2)(Z)
//      = X * 1 * ( 3Z^2-1)/2
//      = - 1/2 X + (3/2) X Z^2
//    so
//      O = 2 (2 nonzero terms)
//      C = -0.5
//           1.5
//      E = 4    <-- index in 3-space of exponent (1,0,0)
//          15   <-- index in 3-space of exponent (1,0,2)
//
//    The output value of O is no greater than
//      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int L[M], the index of each Legendre product 
//    polynomial factor.  0 <= L(*).
//
//    Input, int O_MAX, an upper limit on the size of the 
//    output arrays.
//      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2.
//
//    Output, int &O, the "order" of the polynomial product.
//
//    Output, double C[O], the coefficients of the polynomial product.
//
//    Output, int E[O], the indices of the exponents of the 
//    polynomial product.
//
{
  double *c1;
  double *c2;
  int *e1;
  int *e2;
  int *f2;
  int i;
  int i1;
  int i2;
  int j1;
  int j2;
  int o1;
  int o2;
  int *p;
  int *pp;

  c1 = new double[o_max];
  c2 = new double[o_max];
  e1 = new int[o_max];
  e2 = new int[o_max];
  f2 = new int[o_max];
  pp = new int[m];

  o1 = 1;
  c1[0] = 1.0;
  e1[0] = 1;
//
//  Implicate one factor at a time.
//
  for ( i = 0; i < m; i++ )
  {
    lp_coefficients ( l[i], o2, c2, f2 );
 
    o = 0;

    for ( j2 = 0; j2 < o2; j2++ )
    {
      for ( j1 = 0; j1 < o1; j1++ )
      {
        c[o] = c1[j1] * c2[j2];
        if ( 0 < i )
        {
          p = mono_unrank_grlex ( i, e1[j1] );
        }
        for ( i2 = 0; i2 < i; i2++ )
        {
          pp[i2] = p[i2];
        }
        pp[i] = f2[j2];
        e[o] = mono_rank_grlex ( i + 1, pp );
        o = o + 1;
        if ( 0 < i )
        {
          delete [] p;
        }
      }
    }

    polynomial_sort ( o, c, e );
    polynomial_compress ( o, c, e, o, c, e );

    o1 = o;
    for ( i1 = 0; i1 < o; i1++ )
    {
      c1[i1] = c[i1];
      e1[i1] = e[i1];
    }
  }

  delete [] c1;
  delete [] c2;
  delete [] e1;
  delete [] e2;
  delete [] f2;
  delete [] pp;

  return;
}
//****************************************************************************80

int mono_between_enum ( int d, int n1, int n2 )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_ENUM enumerates monomials in D dimensions of degrees in a range.
//
//  Discussion:
//
//    For D = 3, we have the following table:
//
//     N2 0  1  2  3  4  5  6   7   8
//   N1 +----------------------------
//    0 | 1  4 10 20 35 56 84 120 165
//    1 | 0  3  9 19 34 55 83 119 164
//    2 | 0  0  6 16 31 52 80 116 161
//    3 | 0  0  0 10 25 46 74 110 155
//    4 | 0  0  0  0 15 36 64 100 145
//    5 | 0  0  0  0  0 21 49  85 130
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N1, N2, the minimum and maximum degrees.
//    0 <= N1 <= N2.
//
//    Output, int MONO_BETWEEN_ENUM, the number of monomials 
//    in D variables, of total degree between N1 and N2 inclusive.
//
{
  int n0;
  int n1_copy;
  int value;

  n1_copy = i4_max ( n1, 0 );

  if ( n2 < n1_copy )
  {
    value = 0;
    return value;
  }

  if ( n1_copy == 0 )
  {
    value = i4_choose ( n2 + d, n2 );
  }
  else if ( n1_copy == n2 )
  {
    value = i4_choose ( n2 + d - 1, n2 );
  }
  else
  {
    n0 = n1_copy - 1;
    value = i4_choose ( n2 + d, n2 ) - i4_choose ( n0 + d, n0 );
  }

  return value;
}
//****************************************************************************80

void mono_between_next_grlex ( int d, int n1, int n2, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_BETWEEN_NEXT_GRLEX: grlex next monomial, degree between N1 and N2.
//
//  Discussion:
//
//    We consider all monomials in a D dimensional space, with total
//    degree N between N1 and N2, inclusive.
//
//    For example:
//
//    D = 3
//    N1 = 2
//    N2 = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     2        2
//    2 |  0     1     1        2
//    3 |  0     2     0        2
//    4 |  1     0     1        2
//    5 |  1     1     0        2
//    6 |  2     0     0        2
//      |
//    7 |  0     0     3        3
//    8 |  0     1     2        3
//    9 |  0     2     1        3
//   10 |  0     3     0        3
//   11 |  1     0     2        3
//   12 |  1     1     1        3
//   13 |  1     2     0        3
//   14 |  2     0     1        3
//   15 |  2     1     0        3
//   16 |  3     0     0        3
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
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N1, N2, the minimum and maximum degrees.
//    0 <= N1 <= N2.
//
//    Input/output, int X[D], the current monomial.
//    To start the sequence, set X = [ 0, 0, ..., 0, N1 ].
//    The last value in the sequence is X = [ N2, 0, ..., 0, 0 ].
//
{
  if ( n1 < 0 )
  {
    cerr << "\n";
    cerr << "MONO_BETWEEN_NEXT_GRLEX - Fatal error!\n";
    cerr << "  N1 < 0.\n";
    exit ( 1 );
  }

  if ( n2 < n1 )
  {
    cerr << "\n";
    cerr << "MONO_BETWEEN_NEXT_GRLEX - Fatal error!\n";
    cerr << "  N2 < N1.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( d, x ) < n1 )
  {
    cerr << "\n";
    cerr << "MONO_BETWEEN_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X sums to less than N1.\n";
    exit ( 1 );
  }

  if ( n2 < i4vec_sum ( d, x ) )
  {
    cerr << "\n";
    cerr << "MONO_BETWEEN_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X sums to more than N2.\n";
    exit ( 1 );
  }

  if ( n2 == 0 )
  {
    return;
  }

  if ( x[0] == n2 )
  {
    x[0] = 0;
    x[d-1] = n1;
  }
  else
  {
    mono_next_grlex ( d, x );
  }

  return;
}
//****************************************************************************80

void mono_next_grlex ( int d, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_NEXT_GRLEX returns the next monomial in grlex order.
//
//  Discussion:
//
//    Example:
//
//    D = 3
//
//    #  X(1)  X(2)  X(3)  Degree
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
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input/output, int X[D], the current monomial.
//    The first element is X = [ 0, 0, ..., 0, 0 ].
//
{
  int i;
  int im1;
  int j;
  int t;
//
//  Ensure that 1 <= D.
//
  if ( d < 1 )
  {
    cerr << "\n";
    cerr << "MONO_NEXT_GRLEX - Fatal error!\n";
    cerr << "  D < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 0 <= X(I).
//
  for ( i = 0; i < d; i++ )
  {
    if ( x[i] < 0 )
    {
      cerr << "\n";
      cerr << "MONO_NEXT_GRLEX - Fatal error!\n";
      cerr << "  X[I] < 0\n";
      exit ( 1 );
    }
  }
//
//  Find I, the index of the rightmost nonzero entry of X.
//
  i = 0;
  for ( j = d; 1 <= j; j-- )
  {
    if ( 0 < x[j-1] )
    {
      i = j;
      break;
    }
  }
//
//  set T = X(I)
//  set X(I) to zero,
//  increase X(I-1) by 1,
//  increment X(D) by T-1.
//
  if ( i == 0 )
  {
    x[d-1] = 1;
    return;
  }
  else if ( i == 1 )
  {
    t = x[0] + 1;
    im1 = d;
  }
  else if ( 1 < i )
  {
    t = x[i-1];
    im1 = i - 1;
  }

  x[i-1] = 0;
  x[im1-1] = x[im1-1] + 1;
  x[d-1] = x[d-1] + t - 1;

  return;
}
//****************************************************************************80

int mono_rank_grlex ( int m, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_RANK_GRLEX computes the graded lexicographic rank of a monomial.
//
//  Discussion:
//
//    The graded lexicographic ordering is used, over all monomials of
//    dimension M, for degree NM = 0, 1, 2, ...
//
//    For example, if M = 3, the ranking begins:
//
//    Rank  Sum    1  2  3
//    ----  ---   -- -- --
//       1    0    0  0  0
//
//       2    1    0  0  1
//       3    1    0  1  0
//       4    1    1  0  1
//
//       5    2    0  0  2
//       6    2    0  1  1
//       7    2    0  2  0
//       8    2    1  0  1
//       9    2    1  1  0
//      10    2    2  0  0
//
//      11    3    0  0  3
//      12    3    0  1  2
//      13    3    0  2  1
//      14    3    0  3  0
//      15    3    1  0  2
//      16    3    1  1  1
//      17    3    1  2  0
//      18    3    2  0  1
//      19    3    2  1  0
//      20    3    3  0  0
//
//      21    4    0  0  4
//      ..   ..   .. .. ..
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    1 <= D.
//
//    Input, int XC[M], the monomial.
//    For each 1 <= I <= M, we have 0 <= XC(I).
//
//    Output, int MONO_RANK_GRLEX, the rank.
//
{
  int i;
  int j;
  int ks;
  int n;
  int nm;
  int ns;
  int rank;
  int tim1;
  int *xs;
//
//  Ensure that 1 <= M.
//
  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "MONO_RANK_GRLEX - Fatal error!\n";
    cerr << "  M < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 0 <= X(I).
//
  for ( i = 0; i < m; i++ )
  {
    if ( x[i] < 0 )
    {
      cerr << "\n";
      cerr << "MONO_RANK_GRLEX - Fatal error!\n";
      cerr << "  X[I] < 0\n";
      exit ( 1 );
    }
  }
//
//  NM = sum ( X )
//
  nm = i4vec_sum ( m, x );
//
//  Convert to KSUBSET format.
//
  ns = nm + m - 1;
  ks = m - 1;
  xs = new int[ks];
  xs[0] = x[0] + 1;
  for ( i = 2; i < m; i++ )
  {
    xs[i-1] = xs[i-2] + x[i-1] + 1;
  }
//
//  Compute the rank.
//
  rank = 1;

  for ( i = 1; i <= ks; i++ )
  {
    if ( i == 1 )
    {
      tim1 = 0;
    }
    else
    {
      tim1 = xs[i-2];
    }

    if ( tim1 + 1 <= xs[i-1] - 1 )
    {
      for ( j = tim1 + 1; j <= xs[i-1] - 1; j++ )
      {
        rank = rank + i4_choose ( ns - j, ks - i );
      }
    }
  }

  for ( n = 0; n < nm; n++ )
  {
    rank = rank + i4_choose ( n + m - 1, n );
  }

  delete [] xs;

  return rank;
}
//****************************************************************************80

int mono_total_enum ( int d, int n )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_ENUM enumerates monomials in D dimensions of degree equal to N.
//
//  Discussion:
//
//    For D = 3, we have the following values:
//
//    N  VALUE
//
//    0    1
//    1    3
//    2    6
//    3   10
//    4   15
//    5   21
//
//    In particular, VALUE(3,3) = 10 because we have the 10 monomials:
//
//      x^3, x^2y, x^2z, xy^2, xyz, xz^3, y^3, y^2z, yz^2, z^3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N, the maximum degree.
//
//    Output, int MONO_TOTAL_ENUM, the number of monomials in D variables,
//    of total degree N.
//
{
  int value;

  value = i4_choose ( n + d - 1, n );

  return value;
}
//****************************************************************************80

void mono_total_next_grlex ( int d, int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_NEXT_GRLEX: grlex next monomial with total degree equal to N.
//
//  Discussion:
//
//    We consider all monomials in a D dimensional space, with total degree N.
//
//    For example:
//
//    D = 3
//    N = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     3        3
//    2 |  0     1     2        3
//    3 |  0     2     1        3
//    4 |  0     3     0        3
//    5 |  1     0     2        3
//    6 |  1     1     1        3
//    7 |  1     2     0        3
//    8 |  2     0     1        3
//    9 |  2     1     0        3
//   10 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N, the degree.
//    0 <= N.
//
//    Input/output, int X[D], the current monomial.
//    To start the sequence, set X = [ 0, 0, ..., 0, N ].
//    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
//
{
  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "MONO_TOTAL_NEXT_GRLEX - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( d, x ) != n )
  {
    cerr << "\n";
    cerr << "MONO_TOTAL_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X does not sum to N.\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    return;
  }

  if ( x[0] == n )
  {
    x[0] = 0;
    x[d-1] = n;
  }
  else
  {
    mono_next_grlex ( d, x );
  }

  return;
}
/******************************************************************************/

int *mono_unrank_grlex ( int d, int rank )

/******************************************************************************/
//
//  Purpose:
//
//    MONO_UNRANK_GRLEX computes the composition of given grlex rank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//    1 <= D.
//
//    Input, int RANK, the rank.
//    1 <= RANK.
//
//    Output, int MONO_UNRANK_GRLEX[D], the monomial X of the given rank.
//    For each I, 0 <= XC[I] <= NM, and 
//    sum ( 1 <= I <= D ) XC[I] = NM.
//
{
  int i;
  int j;
  int ks;
  int nksub;
  int nm;
  int ns;
  int r;
  int rank1;
  int rank2;
  int *x;
  int *xs;
//
//  Ensure that 1 <= D.
//
  if ( d < 1 )
  {
    cerr << "\n";
    cerr << "MONO_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  D < 1\n";
    cerr << "  D = " << d << "\n";
    exit ( 1 );
  }
//
//  Ensure that 1 <= RANK.
//
  if ( rank < 1 )
  {
    cerr << "\n";
    cerr << "MONO_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  RANK < 1\n";
    cerr << "  RANK = " << rank << "\n";
    exit ( 1 );
  }
//
//  Special case D == 1.
//
  if ( d == 1 )
  {
    x = new int[d];
    x[0] = rank - 1;
    return x;
  }
//
//  Determine the appropriate value of NM.
//  Do this by adding up the number of compositions of sum 0, 1, 2, 
//  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
//  gives you the rank of the composition within the set of compositions
//  of sum NM.  And that's the number you need in order to do the
//  unranking.
//
  rank1 = 1;
  nm = -1;
  for ( ; ; )
  {
    nm = nm + 1;
    r = i4_choose ( nm + d - 1, nm );
    if ( rank < rank1 + r )
    {
      break;
    }
    rank1 = rank1 + r;
  }

  rank2 = rank - rank1;
//
//  Convert to KSUBSET format.
//  Apology: an unranking algorithm was available for KSUBSETS,
//  but not immediately for compositions.  One day we will come back
//  and simplify all this.
//
  ks = d - 1;
  ns = nm + d - 1;
  xs = new int[ks];

  nksub = i4_choose ( ns, ks );

  j = 1;

  for ( i = 1; i <= ks; i++ )
  {
    r = i4_choose ( ns - j, ks - i );

    while ( r <= rank2 && 0 < r )
    {
      rank2 = rank2 - r;
      j = j + 1;
      r = i4_choose ( ns - j, ks - i );
    }
    xs[i-1] = j;
    j = j + 1;
  }
//
//  Convert from KSUBSET format to COMP format.
//
  x = new int[d];
  x[0] = xs[0] - 1;
  for ( i = 2; i < d; i++ )
  {
    x[i-1] = xs[i-1] - xs[i-2] - 1;
  }
  x[d-1] = ns - xs[ks-1];

  delete [] xs;

  return x;
}
//****************************************************************************80

int mono_upto_enum ( int d, int n )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_ENUM enumerates monomials in D dimensions of degree up to N.
//
//  Discussion:
//
//    For D = 2, we have the following values:
//
//    N  VALUE
//
//    0    1
//    1    3
//    2    6
//    3   10
//    4   15
//    5   21
//
//    In particular, VALUE(2,3) = 10 because we have the 10 monomials:
//
//      1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int N, the maximum degree.
//
//    Output, int MONO_UPTO_ENUM, the number of monomials in
//    D variables, of total degree N or less.
//
{
  int value;

  value = i4_choose ( n + d, n );

  return value;
}
//****************************************************************************80

double *mono_value ( int d, int nx, int f[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_VALUE evaluates a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int NX, the number of evaluation points.
//
//    Input, int F[D], the exponents of the monomial.
//
//    Input, double X[D*NX], the coordinates of the evaluation points.
//
//    Output, double MONO_VALUE[NX], the value of the monomial at X.
//
{
  int i;
  int j;
  double *v;

  v = new double[nx];

  for ( j = 0; j < nx; j++ )
  {
    v[j] = 1.0;
    for ( i = 0; i < d; i++ )
    {
      v[j] = v[j] * pow ( x[i+j*d], f[i] );
    }
  }

  return v;
}
//****************************************************************************80

void perm_check0 ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK0 checks a 0-based permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from 0 to
//    to N-1 occurs among the N entries of the permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
{
  int ierror;
  int location;
  int value;

  for ( value = 0; value < n; value++ )
  {
    ierror = 1;

    for ( location = 0; location < n; location++ )
    {
      if ( p[location] == value )
      {
        ierror = 0;
        break;
      }
    }

    if ( ierror != 0 )
    {
      cerr << "\n";
      cerr << "PERM_CHECK0 - Fatal error!\n";
      cerr << "  Permutation is missing value " << value << "\n";
      exit ( 1 );
    }

  }

  return;
}
//****************************************************************************80

void polynomial_axpy ( double s, int o1, double c1[], int e1[], int o2, 
  double c2[], int e2[], int &o, double c[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_AXPY adds a multiple of one polynomial to another.
//
//  Discussion:
//
//    P(X) = S * P1(X) + P2(X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double S, the multiplier of polynomial 1.
//
//    Input, int O1, the "order" of polynomial 1.
//
//    Input, double C1[O1], the coefficients of polynomial 1.
//
//    Input, int E1[O1], the indices of the exponents of 
//    polynomial 1.
//
//    Input, int O2, the "order" of polynomial 2.
//
//    Input, double C2[O2], the coefficients of polynomial 2.
//
//    Input, int E2[O2], the indices of the exponents of 
//    polynomial 2.
//
//    Output, int &O, the "order" of the polynomial sum.
//
//    Output, double C[O], the coefficients of the polynomial sum.
//
//    Output, int E[O], the indices of the exponents of 
//    the polynomial sum.
//
{
  double *c3;
  int *e3;
  int i;
  int o3;
  double *sc1;

  o3 = o1 + o2;

  c3 = new double[o3];
  e3 = new int[o3];
  sc1 = new double[o1];

  for ( i = 0; i < o1; i++ )
  {
    sc1[i] = s * c1[i];
  }
  r8vec_concatenate ( o1, sc1, o2, c2, c3 );
  i4vec_concatenate ( o1, e1, o2, e2, e3 );

  polynomial_sort ( o3, c3, e3 );
  polynomial_compress ( o3, c3, e3, o, c, e );

  delete [] c3;
  delete [] e3;
  delete [] sc1;

  return;
}
//****************************************************************************80

void polynomial_compress ( int o1, double c1[], int e1[], int &o2, double c2[], 
  int e2[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_COMPRESS compresses a polynomial.
//
//  Discussion:
//
//    The function polynomial_sort ( ) should be called first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int O1, the "order" of the polynomial.
//
//    Input, double C1[O1], the coefficients of the polynomial.
//
//    Input, int E1[O1], the indices of the exponents of 
//    the polynomial.
//
//    Output, int &O2, the "order" of the polynomial.
//
//    Output, double C2[O2], the coefficients of the polynomial.
//
//    Output, int E2[O2], the indices of the exponents of 
//    the polynomial.
//
{
  int get;
  int put;
  const double r8_epsilon_sqrt = 0.1490116119384766E-07;

  get = 0;
  put = 0;

  while ( get < o1 )
  {
    get = get + 1;

    if ( fabs ( c1[get-1] ) <= r8_epsilon_sqrt )
    {
      continue;
    }

    if ( 0 == put )
    {
      put = put + 1;
      c2[put-1] = c1[get-1];
      e2[put-1] = e1[get-1];
    }
    else
    {
      if ( e2[put-1] == e1[get-1] )
      {
        c2[put-1] = c2[put-1] + c1[get-1];
      }
      else
      {
        put = put + 1;
        c2[put-1] = c1[get-1];
        e2[put-1] = e1[get-1];
       }
    }
  }
 
  o2 = put;

  return;
}
//****************************************************************************80

void polynomial_print ( int d, int o, double c[], int e[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_PRINT prints a polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int O, the "order" of the polynomial, that is,
//    simply the number of terms.
//
//    Input, double C[O], the coefficients.
//
//    Input, int E[O], the indices of the exponents.
//
//    Input, string TITLE, a title.
//
{
  int *f;
  int i;
  int j;

  cout << title << "\n";

  if ( o == 0 )
  {
    cout << "      0.\n";
  }
  else
  {
    for ( j = 0; j < o; j++ )
    {
      cout << "    ";
      if ( c[j] < 0.0 )
      {
        cout << "- ";
      }
      else
      {
        cout << "+ ";
      }
      cout << fabs ( c[j] ) << " * x^(";

      f = mono_unrank_grlex ( d, e[j] );
      for ( i = 0; i < d; i++ )
      {
        cout << f[i];
        if ( i < d - 1 )
        {
          cout << ",";
        }
        else
        {
          cout << ")";
        }
      }
      delete [] f;

      if ( j == o - 1 )
      {
        cout << ".";
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

void polynomial_sort ( int o, double c[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_SORT sorts the information in a polynomial.
//
//  Discussion
//
//    The coefficients C and exponents E are rearranged so that 
//    the elements of E are in ascending order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int O, the "order" of the polynomial.
//
//    Input/output, double C[O], the coefficients of the polynomial.
//
//    Input/output, int E[O], the indices of the exponents of 
//    the polynomial.
//
{
  int *indx;

  indx = i4vec_sort_heap_index_a ( o, e );

  i4vec_permute ( o, indx, e );
  r8vec_permute ( o, indx, c );

  delete [] indx;

  return;
}
//****************************************************************************80

double *polynomial_value ( int d, int o, double c[], int e[], int nx, 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_VALUE evaluates a polynomial.
//
//  Discussion:
//
//    The polynomial is evaluated term by term, and no attempt is made to
//    use an approach such as Horner's method to speed up the process.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int D, the spatial dimension.
//
//    Input, int O, the "order" of the polynomial.
//
//    Input, double C[O], the coefficients of the polynomial.
//
//    Input, int E(O), the indices of the exponents 
//    of the polynomial.
//
//    Input, int NX, the number of evaluation points.
//
//    Input, double X[D*NX], the coordinates of the evaluation points.
//
//    Output, double POLYNOMIAL_VALUE[NX], the value of the polynomial at X.
//
{
  int *f;
  int j;
  int k;
  double *p;
  double *v;

  p = new double[nx];

  for ( k = 0; k < nx; k++ )
  {
    p[k] = 0.0;
  }

  for ( j = 0; j < o; j++ )
  {
    f = mono_unrank_grlex ( d, e[j] );
    v = mono_value ( d, nx, f, x );
    for ( k = 0; k < nx; k++ )
    {
      p[k] = p[k] + c[j] * v[k];
    }
    delete [] f;
    delete [] v;
  }

  return p;
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

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

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
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
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
  double value;

  value = 1.0E+30;

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

void r8col_separation ( int m, int n, double a[], double &d_min, double &d_max )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_SEPARATION returns the "separation" of an R8COL.
//
//  Discussion:
//
//    D_MIN is the minimum distance between two columns,
//    D_MAX is the maximum distance between two columns.
//
//    The distances are measured using the Loo norm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns 
//    in the array.  If N < 2, it does not make sense to call this routine.
//
//    Input, double A[M*N], the array whose variances are desired.
//
//    Output, double &D_MIN, &D_MAX, the minimum and maximum distances.
//
{
  double d;
  int i;
  int j1;
  int j2;

  d_min = r8_huge ( );
  d_max = 0.0;

  for ( j1 = 0; j1 < n; j1++ )
  {
    for ( j2 = j1 + 1; j2 < n; j2++ )
    {
      d = 0.0;
      for ( i = 0; i < m; i++ )
      {
        d = r8_max ( d, r8_abs ( a[i+j1*m] - a[i+j2*m] ) );
      }
      d_min = r8_min ( d_min, d );
      d_max = r8_max ( d_max, d );
    }
  }

  return;
}
//****************************************************************************80

double r8mat_is_identity ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
//
//  Discussion:
//
//    An R8MAT is a matrix of real ( kind = 8 ) values.
//
//    The routine returns the Frobenius norm of A - I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Output, double R8MAT_IS_IDENTITY, the Frobenius norm
//    of the difference matrix A - I, which would be exactly zero
//    if A were the identity matrix.
//
{
  double error_frobenius;
  int i;
  int j;
  double t;

  error_frobenius = 0.0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        t = a[i+j*n] - 1.0;
      }
      else
      {
        t = a[i+j*n];
      }
      error_frobenius = error_frobenius + t * t;
    }
  }
  error_frobenius = sqrt ( error_frobenius );

  return error_frobenius;
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
//    20 August 2010
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

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
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
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

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

void r8vec_concatenate ( int n1, double a[], int n2, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CONCATENATE concatenates two R8VEC's.
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
//    22 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, the number of entries in the first vector.
//
//    Input, double A[N1], the first vector.
//
//    Input, int N2, the number of entries in the second vector.
//
//    Input, double B[N2], the second vector.
//
//    Output, double C[N1+N2], the concatenated vector.
//
{
  int i;

  for ( i = 0; i < n1; i++ )
  {
    c[i] = a[i];
  }
  for ( i = 0; i < n2; i++ )
  {
    c[n1+i] = b[i];
  }

  return;
}
//****************************************************************************80

void r8vec_permute ( int n, int p[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PERMUTE permutes an R8VEC in place.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   1,   3,   4,   0,   2 )
//      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
//
//    Output:
//
//      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.
//
//    Input/output, double A[N], the array to be permuted.
//
{
  double a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  perm_check0 ( n, p );
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is 0.
//  So temporarily add 1 to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1;
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "R8VEC_PERMUTE - Fatal error!\n";
          cerr << "  A permutation index is out of range.\n";
          cerr << "  P(" << iput << ") = " << iget << "\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
//
//  Restore the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1;
  }
  return;
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
