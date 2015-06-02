# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
int *abscissa_level_closed_nd ( int level_max, int dim_num, int test_num, 
  int test_val[] );
double cc_abscissa ( int order, int i );
double *cc_weights ( int n );
int choose ( int n, int k );
void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
string i4_to_string ( int i4, string format );
int i4vec_product ( int n, int a[] );
int index_to_level_closed ( int dim_num, int t[], int order, int level_max );
void level_to_order_closed ( int dim_num, int level[], int order[] );
int *multigrid_index0 ( int dim_num, int order_1d[], int order_nd );
void multigrid_scale_closed ( int dim_num, int order_nd, int level_max, 
  int level_1d[], int grid_index[] );
double *product_weights_cc ( int dim_num, int order_1d[], int order_nd );
double r8_epsilon ( );
double r8_huge ( );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title );
int s_len_trim ( string s );
void sparse_grid_cc ( int dim_num, int level_max, int point_num, 
  double grid_weight[], double grid_point[] );
int *sparse_grid_cc_index ( int dim_num, int level_max, int point_num );
int sparse_grid_cfn_size ( int dim_num, int level_max );
void sparse_grid_cc_weights ( int dim_num, int level_max, int point_num, 
  int grid_index[], double grid_weight[] );
void timestamp ( );
void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_CC_DATASET.
//
//  Discussion:
//
//    This program computes a sparse grid quadrature rule based on 1D
//    Clenshaw-Curtis rules and writes it to a file.
//
//    The user specifies:
//    * the spatial dimension of the quadrature region,
//    * the level that defines the Smolyak grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
{
  int dim;
  int dim_num;
  int level_max;
  int level_min;
  int point;
  int point_num;
  double *r;
  string r_filename;
  double *w;
  string w_filename;
  double weight_sum;
  double *x;
  string x_filename;

  timestamp ( );
  cout << "\n";
  cout << "SPARSE_GRID_CC_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Compute the abscissas and weights of a quadrature rule\n";
  cout << "  associated with a sparse grid derived from a Smolyak\n";
  cout << "  construction based on 1D Clenshaw-Curtis rules.\n";
  cout << "\n";
  cout << "  Inputs to the program include:\n";
  cout << "\n";
  cout << "    DIM_NUM, the spatial dimension.\n";
  cout << "    (typically in the range of 2 to 10)\n";
  cout << "\n";
  cout << "    LEVEL_MAX, the level of the sparse grid.\n";
  cout << "    (typically in the range of 0, 1, 2, 3, ...\n";
  cout << "\n";
  cout << "  Output from the program includes:\n";
  cout << "\n";
  cout << "    * A printed table of the abscissas and weights.\n";
  cout << "\n";
  cout << "    * A set of 3 files that define the quadrature rule.\n";
  cout << "\n";
  cout << "    (1) cc_d?_level?_r.txt, the ranges;\n";
  cout << "    (2) cc_d?_level?_w.txt, the weights;\n";
  cout << "    (3) cc_d?_level?_x.txt, the abscissas.\n";
//
//  Get the spatial dimension.
//
  if ( 1 < argc )
  {
    dim_num = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of DIM_NUM (1 or greater)\n";
    cin >> dim_num;
  }

  cout << "\n";
  cout << "  Spatial dimension requested is = " << dim_num << "\n";
//
//  Get the level.
//
  if ( 2 < argc )
  {
    level_max = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of LEVEL_MAX (0 or greater).\n";
    cin >> level_max;
  }

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  cout << "\n";
  cout << "  LEVEL_MIN is = " << level_min << "\n";
  cout << "  LEVEL_MAX is = " << level_max << "\n";
// 
//  How many distinct points will there be?
//
  point_num = sparse_grid_cfn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  The number of distinct abscissas in the\n";
  cout << "  quadrature rule is determined from the spatial\n";
  cout << "  dimension DIM_NUM and the level LEVEL_MAX.\n";
  cout << "  For the given input, this value will be = " << point_num << "\n";
//
//  Allocate memory.
//
  r = new double[dim_num*2];
  w = new double[point_num];
  x = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  for ( dim = 0; dim < dim_num; dim++ )
  {
    r[dim+0*dim_num] = -1.0;
    r[dim+1*dim_num] = +1.0;
  }

  sparse_grid_cc ( dim_num, level_max, point_num, w, x );

  r8mat_transpose_print_some ( dim_num, point_num, x, 1, 1, dim_num, 
    10, "  First 10 grid points:" );

  r8vec_print_some ( point_num, w, 1, 10, "  First 10 grid weights:" );

  weight_sum = 0.0;
  for ( point = 0; point < point_num; point++ )
  {
    weight_sum = weight_sum + w[point];
  }

  cout << "\n";
  cout << "  Weights sum to   " << weight_sum << "\n";
  cout << "  Correct value is " << pow ( 2.0, dim_num ) << "\n";
//
//  Construct appropriate file names.
//
  r_filename = "cc_d" + i4_to_string ( dim_num, "%d" ) 
    + "_level" + i4_to_string ( level_max, "%d" ) + "_r.txt";
  w_filename = "cc_d" + i4_to_string ( dim_num, "%d" ) 
    + "_level" + i4_to_string ( level_max, "%d" ) + "_w.txt";
  x_filename = "cc_d" + i4_to_string ( dim_num, "%d" ) 
    + "_level" + i4_to_string ( level_max, "%d" ) + "_x.txt";
//
//  Write the rule to files.
//
  cout << "\n";
  cout << "  Creating R file = \"" << r_filename << "\".\n";

  r8mat_write ( r_filename, dim_num, 2, r );

  cout << "  Creating W file = \"" << w_filename << "\".\n";

  r8mat_write ( w_filename, 1, point_num, w );

  cout << "  Creating X file = \"" << x_filename << "\".\n";

  r8mat_write ( x_filename, dim_num, point_num, x );
//
//  Terminate.
//
  delete [] r;
  delete [] w;
  delete [] x;

  cout << "\n";
  cout << "SPARSE_GRID_CC_DATASET\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int *abscissa_level_closed_nd ( int level_max, int dim_num, int test_num, 
  int test_val[] )

//****************************************************************************80
//
//  Purpose:
//
//    ABSCISSA_LEVEL_CLOSED_ND: first level at which an abscissa is generated.
//
//  Discussion:
//
//    We need this routine because the sparse grid is generated as a sum of 
//    product grids, and many points in the sparse grid will belong to several
//    of these product grids, and we need to do something special the very 
//    first time we encounter such a point - namely, count it.  So this routine 
//    determines, for any point in the full product grid, the first level 
//    at which that point would be included.
//
//
//    We assume an underlying product grid.  In each dimension, this product
//    grid has order 2^LEVEL_MAX + 1.
//
//    We will say a sparse grid has total level LEVEL if each point in the
//    grid has a total level of LEVEL or less.
//
//    The "level" of a point is determined as the sum of the levels of the
//    point in each spatial dimension.
//
//    The level of a point in a single spatial dimension I is determined as
//    the level, between 0 and LEVEL_MAX, at which the point's I'th index
//    would have been generated.
//
//
//    This description is terse and perhaps unenlightening.  Keep in mind
//    that the product grid is the product of 1D grids,
//    that the 1D grids are built up by levels, having
//    orders (total number of points ) 1, 3, 5, 9, 17, 33 and so on,
//    and that these 1D grids are nested, so that each point in a 1D grid
//    has a first level at which it appears.
//
//    Our procedure for generating the points of a sparse grid, then, is
//    to choose a value LEVEL_MAX, to generate the full product grid,
//    but then only to keep those points on the full product grid whose
//    LEVEL is less than or equal to LEVEL_MAX.  
//
//
//    Note that this routine is really just testing out the idea of
//    determining the level.  Our true desire is to be able to start
//    with a value LEVEL, and determine, in a straightforward manner,
//    all the points that are generated exactly at that level, or
//    all the points that are generated up to and including that level.
//
//    This allows us to generate the new points to be added to one sparse
//    grid to get the next, or to generate a particular sparse grid at once.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int LEVEL_MAX, controls the size of the final sparse grid.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int TEST_NUM, the number of points to be tested.
//
//    Input, int TEST_VAL[DIM_NUM*TEST_NUM], the indices of the points 
//    to be tested.  Normally, each index would be between 0 and 2^LEVEL_MAX.
//
//    Output, int ABSCISSA_LEVEL_ND[TEST_NUM], the value of LEVEL at which the
//    point would first be generated, assuming that a standard sequence of
//    nested grids is used.
//
{
  int dim;
  int j;
  int order;
  int t;
  int *test_level;

  test_level = new int[test_num];

  if ( level_max == 0 )
  {
    for ( j = 0; j < test_num; j++ )
    {
      test_level[j] = 0;
    }
    return test_level;
  }

  order = i4_power ( 2, level_max ) + 1;

  for ( j = 0; j < test_num; j++ )
  {
    test_level[j] = index_to_level_closed ( dim_num, test_val+j*dim_num, 
      order, level_max );
  }

  return test_level;
}
//****************************************************************************80

double cc_abscissa ( int order, int i )

//****************************************************************************80
//
//  Purpose:
//
//    CC_ABSCISSA returns the I-th abscissa of the Clenshaw Curtis rule.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to
//    right.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//
//    Input, int I, the index of the desired abscissa.  1 <= I <= ORDER.
//
//    Output, double CC_ABSCISSA, the value of the I-th 
//    abscissa in the rule of order ORDER.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( order < 1 )
  {
    cout << "\n";
    cout << "CC_ABSCISSA - Fatal error!\n";
    cout << "  Input value of ORDER < 1.\n";
    cout << "  Input value of ORDER = " << order << "\n";
    exit ( 1 );
  }

  if ( i < 1 || order < i )
  {
    cout << "\n";
    cout << "CC_ABSCISSA - Fatal error!\n";
    cout << "  1 <= I <= ORDER is required.\n";
    cout << "  I = " << i << "\n";
    cout << "  ORDER = " << order << "\n";
    exit ( 1 );
  }

  if ( order == 1 )
  {
    value = 0.0;
    return value;
  }

  value = cos ( ( double ) ( order - i ) * pi 
              / ( double ) ( order - 1 ) );

  if ( 2 * i - 1 == order )
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

double *cc_weights ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CC_WEIGHTS computes Clenshaw Curtis weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Clenshaw, Alan Curtis,
//    A Method for Numerical Integration on an Automatic Computer,
//    Numerische Mathematik,
//    Volume 2, Number 1, December 1960, pages 197-205.
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Output, double CC_WEIGHTS[N], the weights of the rule.
//
{
  double b;
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;
  double *w;

  w = new double[n];

  if ( n == 1 )
  {
    w[0] = 2.0;
    return w;
  }

  for ( i = 1; i <= n; i++ )
  {
    theta = ( double ) ( i - 1 ) * pi / ( double ) ( n - 1 );

    w[i-1] = 1.0;

    for ( j = 1; j <= ( n - 1 ) / 2; j++ )
    {
      if ( 2 * j == ( n - 1 ) )
      {
        b = 1.0;
      }
      else
      {
        b = 2.0;
      }

      w[i-1] = w[i-1] - b * cos ( 2.0 * ( double ) ( j ) * theta ) 
           / ( double ) ( 4 * j * j - 1 );
    }
  }

  w[0] = w[0] / ( double ) ( n - 1 );
  for ( i = 1; i < n-1; i++ )
  {
    w[i] = 2.0 * w[i] / ( double ) ( n - 1 );
  }
  w[n-1] = w[n-1] / ( double ) ( n - 1 );

  return w;
}
//****************************************************************************80

int choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    CHOOSE computes the binomial coefficient C(N,K).
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
//    12 November 2007
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
//    Output, int CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

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
    mx = i4_max ( k, n - k );
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
    }
  }

  return value;
}
//****************************************************************************80

void comp_next ( int n, int k, int a[], bool *more, int *h, int *t )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed 
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//    This routine originally used a SAVE statement to maintain the
//    variables H and T.  I have decided (based on an wasting an
//    entire morning trying to track down a problem) that it is safer
//    to pass these variables as arguments, even though the user should
//    never alter them.  This allows this routine to safely shuffle
//    between several ongoing calculations.
//
//
//    There are 28 compositions of 6 into three parts.  This routine will
//    produce those compositions in the following order:
//
//     I         A
//     -     ---------
//     1     6   0   0
//     2     5   1   0
//     3     4   2   0
//     4     3   3   0
//     5     2   4   0
//     6     1   5   0
//     7     0   6   0
//     8     5   0   1
//     9     4   1   1
//    10     3   2   1
//    11     2   3   1
//    12     1   4   1
//    13     0   5   1
//    14     4   0   2
//    15     3   1   2
//    16     2   2   2
//    17     1   3   2
//    18     0   4   2
//    19     3   0   3
//    20     2   1   3
//    21     1   2   3
//    22     0   3   3
//    23     2   0   4
//    24     1   1   4
//    25     0   2   4
//    26     1   0   5
//    27     0   1   5
//    28     0   0   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
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
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool *MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;

  if ( !( *more ) )
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < *t )
    {
      *h = 0;
    }
    *h = *h + 1;
    *t = a[*h-1];
    a[*h-1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;
  }

  *more = ( a[k-1] != n );

  return;
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
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
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
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
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

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of integer division.
//
//  Formula:
//
//    If 
//      NREM = I4_MODP ( I, J ) 
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//  Comments:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Examples:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
// 
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is 
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cout << "\n";
    cout << "I4_MODP - Fatal error!\n";
    cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

string i4_to_string ( int i4, string format )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  char i4_char[80];
  string i4_string;

  sprintf ( i4_char, format.c_str ( ), i4 );

  i4_string = string ( i4_char );

  return i4_string;
}
//****************************************************************************80

int i4vec_product ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRODUCT multiplies the entries of an integer vector.
//
//  Example:
//
//    A = ( 1, 2, 3, 4 )
//
//    I4VEC_PRODUCT = 24
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector
//
//    Output, int I4VEC_PRODUCT, the product of the entries of A.
//
{
  int i;
  int product;

  product = 1;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
//****************************************************************************80

int index_to_level_closed ( int dim_num, int t[], int order, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_TO_LEVEL_CLOSED determines the level of a point given its index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int T[DIM_NUM], the grid indices of a point in a 1D closed rule.
//    0 <= T[I] <= ORDER.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, int LEVEL_MAX, the level with respect to which the
//    index applies.
//
//    Output, int INDEX_TO_LEVEL_CLOSED, the first level on which
//    the point associated with the given index will appear.
//
{
  int dim;
  int level;
  int s;
  int value;

  value = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    s = t[dim];

    s = i4_modp ( s, order );

    if ( s == 0 )
    {
      level = 0;
    }
    else
    {
      level = level_max;

      while ( ( s % 2 ) == 0 )
      {
        s = s / 2;
        level = level - 1;
      }
    }

    if ( level == 0 )
    {
      level = 1;
    }
    else if ( level == 1 )
    {
      level = 0;
    }
    value = value + level;
  }
  return value;
}
//****************************************************************************80

void level_to_order_closed ( int dim_num, int level[], int order[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEVEL_TO_ORDER_CLOSED converts a level to an order for closed rules.
//
//  Discussion:
//
//    Sparse grids can naturally be nested.  A natural scheme is to use
//    a series of one-dimensional rules arranged in a series of "levels"
//    whose order roughly doubles with each step.
//
//    The arrangement described here works naturally for the Clenshaw Curtis
//    and Newton Cotes closed rules.  
//
//    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single 
//    point at the center, and for all values afterwards, we use the 
//    relationship
//
//      ORDER = 2^LEVEL + 1
//
//    The following table shows how the growth will occur:
//
//    Level    Order
//
//    0          1
//    1          3 =  2 + 1
//    2          5 =  4 + 1
//    3          9 =  8 + 1
//    4         17 = 16 + 1
//    5         33 = 32 + 1
//
//    For the Clenshaw Curtis and Newton Cotes Closed rules, the point growth
//    is nested.  If we have ORDER points on a particular LEVEL, the next
//    level includes all these old points, plus ORDER-1 new points, formed
//    in the gaps between successive pairs of old points.
//
//    Level    Order = New + Old
//
//    0          1   =  1  +  0
//    1          3   =  2  +  1
//    2          5   =  2  +  3
//    3          9   =  4  +  5
//    4         17   =  8  +  9
//    5         33   = 16  + 17
//
//    In this routine, we assume that a vector of levels is given,
//    and the corresponding orders are desired.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL[DIM_NUM], the nesting level.
//
//    Output, int ORDER[DIM_NUM], the order (number of points) 
//    of the rule.
//
{
  int dim;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( level[dim] < 0 )
    {
      order[dim] = -1;
    }
    else if ( level[dim] == 0 )
    {
      order[dim] = 1;
    }
    else
    {
      order[dim] = i4_power ( 2, level[dim] ) + 1 ;
    }
  }
  return;
}
//****************************************************************************80

int *multigrid_index0 ( int dim_num, int order_1d[], int order_nd )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_INDEX0 returns an indexed multidimensional grid.
//
//  Discussion:
//
//    For dimension DIM, the second index of INDX may vary from 
//    0 to ORDER_1D[DIM]-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int ORDER_1D[DIM_NUM], the order of the
//    rule in each dimension.
//
//    Input, int ORDER_ND, the product of the entries of ORDER_1D.
//
//    Output, int INDX[DIM_NUM*ORDER_ND], the indices of the points in
//    the grid.  The second dimension of this array is equal to the
//    product of the entries of ORDER_1D.
//
{
  int *a;
  int dim;
  bool more;
  int p;
  int *indx;

  indx = new int[dim_num*order_nd];
  a = new int[dim_num];
  more = false;
  p = 0;

  for ( ; ; )
  {
    vec_colex_next2 ( dim_num, order_1d, a, &more );

    if ( !more )
    {
      break;
    }

    for ( dim = 0; dim < dim_num; dim++ )
    {
      indx[dim+p*dim_num] = a[dim];
    }
    p = p + 1;
  }

  delete [] a;

  return indx;
}
//****************************************************************************80

void multigrid_scale_closed ( int dim_num, int order_nd, int level_max, 
  int level_1d[], int grid_index[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_SCALE_CLOSED renumbers a grid as a subgrid on a higher level.
//
//  Discussion:
//
//    This routine takes a grid associated with a given value of
//    LEVEL, and multiplies all the indices by a power of 2, so that
//    the indices reflect the position of the same points, but in
//    a grid of level LEVEL_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int ORDER_ND, the number of points in the grid.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, int LEVEL_1D[DIM_NUM], the level in each dimension.
//
//    Input/output, int GRID_INDEX[DIM_NUM*POINT_NUM], the index
//    values for each grid point.  On input, these indices are based in
//    the level for which the grid was generated; on output, the
//    indices are appropriate for the grid as a subgrid of a grid
//    of level LEVEL_MAX.
//
{
  int dim;
  int factor;
  int order;
  int order_max;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( level_1d[dim] == 0 )
    {
      if ( 0 == level_max )
      {
        order_max = 1;
      }
      else
      {
        order_max = i4_power ( 2, level_max ) + 1;
      }
      for ( order = 0; order < order_nd; order++ )
      {
        grid_index[dim+order*dim_num] = ( order_max - 1 ) / 2;
      }
    }
    else
    {
      factor = i4_power ( 2, level_max - level_1d[dim] );
      for ( order = 0; order < order_nd; order++ )
      {
        grid_index[dim+order*dim_num] = grid_index[dim+order*dim_num] * factor;
      }
    }
  }

  return;
}
//****************************************************************************80

double *product_weights_cc ( int dim_num, int order_1d[], int order_nd )

//****************************************************************************80
//
//  Purpose:
//
//    PRODUCT_WEIGHTS_CC computes weights for a Clenshaw Curtis product rule.
//
//  Discussion:
//
//    This routine computes the weights for a quadrature rule which is
//    a product of closed rules of varying order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
//
//    Input, int ORDER_ND, the order of the product rule.
//
//    Output, double PRODUCT_WEIGHTS_CC[DIM_NUM*ORDER_ND], 
//    the product rule weights.
//
{
  int dim;
  int order;
  double *w_1d;
  double *w_nd;

  w_nd = new double[order_nd];

  for ( order = 0; order < order_nd; order++ )
  {
    w_nd[order] = 1.0;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    w_1d = cc_weights ( order_1d[dim] );

    r8vec_direct_product2 ( dim, order_1d[dim], w_1d, dim_num, 
      order_nd, w_nd );

    delete [] w_1d; 
  }

  return w_nd;
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
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file with no header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    To explain what is going on here, suppose we had to construct
//    a multidimensional quadrature rule as the product of K rules
//    for 1D quadrature.
//
//    The product rule will be represented as a list of points and weights.
//
//    The J-th item in the product rule will be associated with
//      item J1 of 1D rule 1,
//      item J2 of 1D rule 2, 
//      ..., 
//      item JK of 1D rule K.
//
//    In particular, 
//      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
//    and
//      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
//
//    So we can construct the quadrature rule if we can properly
//    distribute the information in the 1D quadrature rules.
//
//    This routine carries out that task for the weights W.
//
//    Another way to do this would be to compute, one by one, the
//    set of all possible indices (J1,J2,...,JK), and then index
//    the appropriate information.  An advantage of the method shown
//    here is that you can process the K-th set of information and
//    then discard it.
//
//  Example:
//
//    Rule 1: 
//      Order = 4
//      W(1:4) = ( 2, 3, 5, 7 )
//
//    Rule 2:
//      Order = 3
//      W(1:3) = ( 11, 13, 17 )
//
//    Rule 3:
//      Order = 2
//      W(1:2) = ( 19, 23 )
//
//    Product Rule:
//      Order = 24
//      W(1:24) =
//        ( 2 * 11 * 19 )
//        ( 3 * 11 * 19 )
//        ( 4 * 11 * 19 )
//        ( 7 * 11 * 19 )
//        ( 2 * 13 * 19 )
//        ( 3 * 13 * 19 )
//        ( 5 * 13 * 19 )
//        ( 7 * 13 * 19 )
//        ( 2 * 17 * 19 )
//        ( 3 * 17 * 19 )
//        ( 5 * 17 * 19 )
//        ( 7 * 17 * 19 )
//        ( 2 * 11 * 23 )
//        ( 3 * 11 * 23 )
//        ( 5 * 11 * 23 )
//        ( 7 * 11 * 23 )
//        ( 2 * 13 * 23 )
//        ( 3 * 13 * 23 )
//        ( 5 * 13 * 23 )
//        ( 7 * 13 * 23 )
//        ( 2 * 17 * 23 )
//        ( 3 * 17 * 23 )
//        ( 5 * 17 * 23 )
//        ( 7 * 17 * 23 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR_INDEX, the index of the factor being processed.
//    The first factor processed must be factor 0.
//
//    Input, int FACTOR_ORDER, the order of the factor.
//
//    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values for
//    factor FACTOR_INDEX.
//
//    Input, int FACTOR_NUM, the number of factors.
//
//    Input, int POINT_NUM, the number of elements in the direct product.
//
//    Input/output, double W[POINT_NUM], the elements of the
//    direct product, which are built up gradually.  
//
//  Local Parameters:
//
//    Local, integer START, the first location of a block of values to set.
//
//    Local, integer CONTIG, the number of consecutive values to set.
//
//    Local, integer SKIP, the distance from the current value of START
//    to the next location of a block of values to set.
//
//    Local, integer REP, the number of blocks of values to set.
//
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 1.0;
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( j = 0; j < factor_order; j++ )
  {
    start = 0 + j * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( i = start; i < start + contig; i++ )
      {
        w[i] = w[i] * factor_value[j];
      }
      start = start + skip;
    }
  }

  contig = contig * factor_order;

  return;
}
//****************************************************************************80

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_SOME prints "some" of an R8VEC.
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
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, integer I_LO, I_HI, the first and last indices to print.
//    The routine expects 1 <= I_LO <= I_HI <= N.
//
//    Input, string TITLE, an optional title.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    cout << "  " << setw(8)  << i       
         << "  " << setw(14) << a[i-1]  << "\n";
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
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

void sparse_grid_cc ( int dim_num, int level_max, int point_num, 
  double grid_weight[], double grid_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CC computes a sparse grid of Clenshaw Curtis points.
//
//  Discussion:
//
//    This program computes a quadrature rule and writes it to a file.
//
//    The quadrature rule is associated with a sparse grid derived from
//    a Smolyak construction using a closed 1D quadrature rule. 
//
//    The user specifies:
//    * the spatial dimension of the quadrature region,
//    * the level that defines the Smolyak grid.
//    * the closed 1D quadrature rule (Clenshaw-Curtis or Newton-Cotes Closed).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, controls the size of the final sparse grid.
//
//    Input, int POINT_NUM, the number of points in the grid, as determined
//    by SPARSE_GRID_CC_SIZE.
//
//    Output, double GRID_WEIGHTS[POINT_NUM], the weights.
//
//    Output, double GRID_POINTS[DIM_NUM*POINT_NUM], the points.
//
{
  int dim;
  int *grid_index;
  int order_max;
  int point;
//
//  Determine the index vector, relative to the full product grid,
//  that identifies the points in the sparse grid.
//
  grid_index = sparse_grid_cc_index ( dim_num, level_max, point_num );
//
//  Compute the physical coordinates of the abscissas.
//
  if ( 0 == level_max )
  {
    order_max = 1;
  }
  else
  {
    order_max = i4_power ( 2, level_max ) + 1;
  }

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      grid_point[dim+point*dim_num] = 
        cc_abscissa ( order_max, grid_index[dim+point*dim_num] + 1 );
    }
  }
//
//  Gather the weights.
//
  sparse_grid_cc_weights ( dim_num, level_max, point_num, grid_index, 
    grid_weight );

  delete [] grid_index;

  return;
}
//****************************************************************************80

int *sparse_grid_cc_index ( int dim_num, int level_max, int point_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CC_INDEX indexes the points forming a sparse grid.
//
//  Discussion:
//
//    The points forming the sparse grid are guaranteed to be a subset
//    of a certain product grid.  The product grid is formed by DIM_NUM
//    copies of a 1D rule of fixed order.  The orders of the 1D rule,
//    (called ORDER_1D) and the order of the product grid, (called ORDER)
//    are determined from the value LEVEL_MAX.
//
//    Thus, any point in the product grid can be identified by its grid index,
//    a set of DIM_NUM indices, each between 1 and ORDER_1D.
//
//    This routine creates the GRID_INDEX array, listing (uniquely) the
//    points of the sparse grid.  
//
//    An assumption has been made that the 1D rule is closed (includes
//    the interval endpoints) and nested (points that are part of a rule
//    of a given level will be part of every rule of higher level).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, int POINT_NUM, the total number of points in the grids.
//
//    Output, int SPARSE_GRID_CC_INDEX[DIM_NUM*POINT_NUM], a list of point 
//    indices, representing a subset of the product grid of level LEVEL_MAX,
//    representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
{
  int dim;
  int factor;
  int *grid_index;
  int *grid_index2;
  int *grid_level;
  int h;
  int j;
  int level;
  int *level_1d;
  bool more;
  int *order_1d;
  int order_nd;
  int point;
  int point_num2;
  int t;

  grid_index = new int[dim_num*point_num];
//
//  The outer loop generates LEVELs from 0 to LEVEL_MAX.
//
  point_num2 = 0;

  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( level = 0; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
//  that adds up to LEVEL.
//
    more = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );
//
//  Transform each 1D level to a corresponding 1D order.
//
      level_to_order_closed ( dim_num, level_1d, order_1d );
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//
      grid_index2 = multigrid_index0 ( dim_num, order_1d, order_nd );
//
//  Adjust these grid indices to reflect LEVEL_MAX.
//
      multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, 
        grid_index2 );
//
//  Determine the first level of appearance of each of the points.
//
      grid_level = abscissa_level_closed_nd ( level_max, dim_num, order_nd, 
        grid_index2 );
//
//  Only keep those points which first appear on this level.
//
      for ( point = 0; point < order_nd; point++ )
      {
        if ( grid_level[point] == level )
        {
          for ( dim = 0; dim < dim_num; dim++ )
          {
            grid_index[dim+point_num2*dim_num] =
              grid_index2[dim+point*dim_num];
          }
          point_num2 = point_num2 + 1;
        }
      }

      delete [] grid_index2;
      delete [] grid_level;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] order_1d;

  return grid_index;
}
//****************************************************************************80

int sparse_grid_cfn_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CFN_SIZE sizes a sparse grid using Closed Fully Nested rules.
//
//  Discussion:
//
//    The grid is defined as the sum of the product rules whose LEVEL
//    satisfies:
//
//      0 <= LEVEL <= LEVEL_MAX.
//
//    This calculation is much faster than a previous method.  It simply
//    computes the number of new points that are added at each level in the
//    1D rule, and then counts the new points at a given DIM_NUM dimensional
//    level vector as the product of the new points added in each dimension.
//
//    This approach will work for nested families, and may be extensible
//    to other families, and to mixed rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int SPARSE_GRID_CC_SIZE, the number of points in the grid.
//
{
  int dim;
  int h;
  int i;
  int j;
  int level;
  int *level_1d;
  bool more;
  int *new_1d;
  int point_num;
  int t;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Construct the vector that counts the new points in the 1D rule.
//
  new_1d = new int[level_max+1];

  new_1d[0] = 1;
  new_1d[1] = 2;

  j = 1;
  for ( i = 2; i <= level_max; i++ )
  {
    j = j * 2;
    new_1d[i] = j;
  }
//
//  Count the number of points by counting the number of new points 
//  associated with each level vector.
//
  level_1d = new int[dim_num];

  point_num = 0;

  for ( level = 0; level <= level_max; level++ )
  {
    more = false;
    h = 0;
    t = 0;

    for ( ; ;)
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );

      v = 1;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        v = v * new_1d[level_1d[dim]];
      }

      point_num = point_num + v;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] new_1d;

  return point_num;
}
//****************************************************************************80

void sparse_grid_cc_weights ( int dim_num, int level_max, int point_num, 
  int grid_index[], double grid_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CC_WEIGHTS gathers the weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, int POINT_NUM, the total number of points in the grids.
//
//    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of point indices,
//    representing a subset of the product grid of level LEVEL_MAX,
//    representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
//    Output, double GRID_WEIGHT[POINT_NUM], the weights
//    associated with the sparse grid points.
//
{
  bool all_equal;
  int coeff;
  int dim;
  int *grid_index2;
  double *grid_weight2;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int order_nd;
  int *order_1d;
  int point;
  int point2;
  int t;

  if ( level_max == 0 )
  {
    for ( point = 0; point < point_num; point++ )
    { 
      grid_weight[point] = pow ( 2.0, dim_num );
    }
    return;
  }

  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( point = 0; point < point_num; point++ )
  { 
    grid_weight[point] = 0.0;
  }

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
//  that adds up to LEVEL.
//
    more = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );
//
//  Transform each 1D level to a corresponding 1D order.
//
      level_to_order_closed ( dim_num, level_1d, order_1d );
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  Generate the indices of the points corresponding to the grid.
//
      grid_index2 = multigrid_index0 ( dim_num, order_1d, order_nd );
//
//  Compute the weights for this grid.
//
      grid_weight2 = product_weights_cc ( dim_num, order_1d, order_nd );
//
//  Adjust the grid indices to reflect LEVEL_MAX.
//
      multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, 
        grid_index2 );
//
//  Now determine the coefficient.
//
      coeff = i4_power ( -1, level_max - level ) 
        * choose ( dim_num - 1, level_max - level );

      for ( point2 = 0; point2 < order_nd; point2++ )
      {
        for ( point = 0; point < point_num; point++ )
        {
          all_equal = true;
          for ( dim = 0; dim < dim_num; dim++ )
          {
            if ( grid_index2[dim+point2*dim_num] !=
                 grid_index[dim+point*dim_num] )
            {
              all_equal = false;
              break;
            }
          }
          if ( all_equal )
          {
            grid_weight[point] = grid_weight[point] 
              + ( double ) ( coeff ) * grid_weight2[point2];
            break;
          }
        }
      }

      delete [] grid_index2;
      delete [] grid_weight2;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] order_1d;

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
//****************************************************************************80

void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_COLEX_NEXT2 generates vectors in colex order.
//
//  Discussion:
//
//    The vectors are produced in colexical order, starting with
//
//    (0,        0,        ...,0),
//    (1,        0,        ...,0),
//     ...
//    (BASE(1)-1,0,        ...,0)
//
//    (0,        1,        ...,0)
//    (1,        1,        ...,0)
//    ...
//    (BASE(1)-1,1,        ...,0)
//
//    (0,        2,        ...,0)
//    (1,        2,        ...,0)
//    ...
//    (BASE(1)-1,BASE(2)-1,...,BASE(DIM_NUM)-1).
//
//  Examples:
//
//    DIM_NUM = 2,
//    BASE = { 3, 3 }
//
//    0   0
//    1   0
//    2   0
//    0   1
//    1   1
//    2   1
//    0   2
//    1   2
//    2   2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int BASE[DIM_NUM], the bases to be used in each dimension.
//    In dimension I, entries will range from 0 to BASE[I]-1.
//
//    Output, int A[DIM_NUM], the next vector.
//
//    Input/output, bool *MORE.  Set this variable false before
//    the first call.  On return, MORE is TRUE if another vector has
//    been computed.  If MORE is returned FALSE, ignore the output 
//    vector and stop calling the routine.
//
{
  int i;

  if ( !( *more ) )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = 0;
    }
    *more = true;
  }
  else
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = a[i] + 1;

      if ( a[i] < base[i] )
      {
        return;
      }
      a[i] = 0;
    }
    *more = false;
  }

  return;
}

