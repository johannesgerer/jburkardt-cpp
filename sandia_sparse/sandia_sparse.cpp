# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "sandia_sparse.hpp"

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

int *abscissa_level_open_nd ( int level_max, int dim_num, int test_num, 
  int test_val[] )

//****************************************************************************80
//
//  Purpose:
//
//    ABSCISSA_LEVEL_OPEN_ND: first level at which given abscissa is generated.
//
//  Discussion:
//
//    We assume an underlying product grid.  In each dimension, this product
//    grid has order 2**(LEVEL_MAX+1) - 1.
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
//    orders (total number of points ) 1, 3, 7, 15, 31 and so on,
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
//  Modified:
//
//    19 April 2007
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
//    to be tested.  Normally, each index would be between 0 and 2**LEVEL_MAX.
//
//    Output, int ABSCISSA_OPEN_LEVEL_ND[TEST_NUM], the value of LEVEL at which the
//    point would first be generated, assuming that a standard sequence of
//    nested grids is used.
//
{
  int dim;
  int j;
  int level;
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
    test_level[j] = index_to_level_open ( dim_num, test_val+j*dim_num,
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
//    31 March 2008
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
  }
  else if ( 2 * ( order - i ) == order - 1 )
  {
    value = 0.0;
  }
  else
  {
    value = cos ( ( double ) ( order - i ) * pi 
                / ( double ) ( order - 1 ) );
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
//    09 November 2007
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
//    variables H and T.  I have decided that it is safer
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

double f1_abscissa ( int order, int i )

//****************************************************************************80
//
//  Purpose:
//
//    F1_ABSCISSA returns the I-th abscissa for the Fejer type 1 rule.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to
//    right.
//
//    This rule is defined on [-1,1].
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the Fejer type 1 rule.
//    1 <= ORDER.
//
//    Input, int I, the index of the desired abscissa.  
//    1 <= I <= ORDER.
//
//    Output, double F1_ABSCISSA, the value of the I-th 
//    abscissa in the Fejer type 1 rule of order ORDER.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( order < 1 )
  {
    value = - r8_huge ( );
    return value;
  }

  if ( i < 1 || order < i )
  {
    cout << "\n";
    cout << "F1_ABSCISSA - Fatal error!\n";
    cout << "  1 <= I <= ORDER is required.\n";
    cout << "  I = " << i << "\n";
    cout << "  ORDER = " << order << "\n";
    exit ( 1 );
  }

  if ( order == 1 )
  {
    value = 0.0;
  }
  else if ( 2 * ( 2 * order + 1 - 2 * i ) == 2 * order )
  {
    value = 0.0;
  }
  else
  {
    value = cos ( ( double ) ( 2 * order + 1 - 2 * i ) * pi 
                / ( double ) ( 2 * order ) );
  }

  return value;
}
//****************************************************************************80

double *f1_weights ( int order )

//****************************************************************************80
//
//  Purpose:
//
//    F1_WEIGHTS computes weights for a Fejer type 1 rule.
//
//  Modified:
//
//    28 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Output, double F1_WEIGHTS[ORDER], the weights.
//
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double *theta;
  double *w;

  if ( order < 1 )
  {
    cout << "\n";
    cout << "F1_WEIGHTS - Fatal error!\n";
    cout << "  ORDER < 1.\n";
    exit ( 1 );
  }

  w = new double[order];

  if ( order == 1 )
  {
    w[0] = 2.0;
    return w;
  }

  theta = new double[order];

  for ( i = 1; i <= order; i++ )
  {
    theta[i-1] = ( double ) ( 2 * ( order + 1 - i ) - 1 ) * pi 
               / ( double ) ( 2 * order );
  }

  for ( i = 1; i <= order; i++ )
  {
    w[i-1] = 1.0;
    for ( j = 1; j <= ( order / 2 ); j++ )
    {
      w[i-1] = w[i-1] - 2.0 * cos ( 2.0 * ( double ) ( j ) * theta[i-1] ) 
        / ( double ) ( 4 * j * j - 1 );
    }
  }

  for ( i = 0; i < order; i++ )
  {
    w[i] = 2.0 * w[i] / ( double ) ( order );
  }

  delete [] theta;

  return w;
}
//****************************************************************************80

double f2_abscissa ( int order, int i )

//****************************************************************************80
//
//  Purpose:
//
//    F2_ABSCISSA returns the I-th abscissa for the Fejer type 2 rule.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to
//    right.
//
//    This rule is defined on [-1,1].
//
//  Modified:
//
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the Fejer type 2 rule.
//    1 <= ORDER.
//
//    Input, int I, the index of the desired abscissa.  1 <= I <= ORDER.
//
//    Output, double F2_ABSCISSA, the value of the I-th 
//    abscissa in the Fejer type 2 rule of order ORDER.
//
{
  double pi = 3.141592653589793;
  double value;

  if ( order < 1 )
  {
    value = - r8_huge ( );
    return value;
  }

  if ( i < 1 || order < i )
  {
    cout << "\n";
    cout << "F2_ABSCISSA - Fatal error!\n";
    cout << "  1 <= I <= ORDER is required.\n";
    cout << "  I = " << i << "\n";
    cout << "  ORDER = " << order << "\n";
    exit ( 1 );
  }

  if ( order == 1 )
  {
    value = 0.0;
  }
  else if ( 2 * ( order + 1 - i ) == order + 1 )
  {
    value = 0.0;
  }
  else
  {
    value = cos ( ( double ) ( order + 1 - i ) * pi 
                / ( double ) ( order + 1 ) );
  }

  return value;
}
//****************************************************************************80

double *f2_weights ( int order )

//****************************************************************************80
//
//  Purpose:
//
//    F2_WEIGHTS computes weights for a Fejer type 2 rule.
//
//  Modified:
//
//    28 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Output, double F2_WEIGHTS[ORDER], the weights.
//
{
  int i;
  int j;
  double p;
  double pi = 3.141592653589793;
  double *theta;
  double *w;

  if ( order < 1 )
  {
    cout << "\n";
    cout << "F2_WEIGHTS - Fatal error!\n";
    cout << "  ORDER < 1.\n";
    exit ( 1 );
  }

  w = new double[order];

  if ( order == 1 )
  {
    w[0] = 2.0;
    return w;
  }
  else if ( order == 2 )
  {
    w[0] = 1.0;
    w[1] = 1.0;
    return w;
  }

  theta = new double[order];

  for ( i = 1; i <= order; i++ )
  {
    theta[i-1] = ( double ) ( order + 1 - i ) * pi 
               / ( double ) ( order + 1 );
  }

  for ( i = 1; i <= order; i++ )
  {
    w[i-1] = 1.0;

    for ( j = 1; j <= ( ( order - 1 ) / 2 ); j++ )
    {
      w[i-1] = w[i-1] - 2.0 * cos ( 2.0 * ( double ) ( j ) * theta[i-1] ) 
        / ( double ) ( 4 * j * j - 1 );
    }

    if ( 2 < order )
    {
      p = 2.0 * ( double ) ( ( ( order + 1 ) / 2 ) ) - 1.0;
      w[i-1] = w[i-1] - cos ( ( p + 1.0 ) * theta[i-1] ) / p;
    }

  }

  for ( i = 0; i < order; i++ )
  {
    w[i] = 2.0 * w[i] / ( double ) ( order + 1 );
  }
  delete [] theta;

  return w;
}
//****************************************************************************80

void gh_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    GH_ABSCISSA sets abscissas for multidimensional Gauss-Hermite quadrature.
//
//  Discussion:
//
//    The "nesting" as it occurs for Gauss-Hermite sparse grids simply
//    involves the use of a specified set of permissible orders for the
//    rule.  
//
//    The X array lists the (complete) Gauss-Legendre abscissas for rules 
//    of order 1, 3, 7, 15, 31, 63 or 127, in order. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the index of the abscissa
//    from the rule, for each dimension and point.
//
//    Input, int GRID_BASE[DIM_NUM], the number of points used in the 
//    rule for a given dimension.
//
//    Output, double GRID_POINT[DIM_NUM], the grid points of abscissas.
//
{
  int dim;
  int level;
  int point;
  int pointer;
  int skip[8] = { 0, 1, 4, 11, 26, 57, 120, 247 };
  double x[247] = {
    0.0E+00, 
   -0.122474487139158904909864203735E+01, 
    0.0E+00, 
    0.122474487139158904909864203735E+01, 
   -0.265196135683523349244708200652E+01, 
   -0.167355162876747144503180139830E+01, 
   -0.816287882858964663038710959027E+00, 
    0.0E+00, 
    0.816287882858964663038710959027E+00, 
    0.167355162876747144503180139830E+01, 
    0.265196135683523349244708200652E+01, 
   -0.449999070730939155366438053053E+01, 
   -0.366995037340445253472922383312E+01, 
   -0.296716692790560324848896036355E+01, 
   -0.232573248617385774545404479449E+01, 
   -0.171999257518648893241583152515E+01, 
   -0.113611558521092066631913490556E+01, 
   -0.565069583255575748526020337198E+00, 
    0.0E+00, 
    0.565069583255575748526020337198E+00, 
    0.113611558521092066631913490556E+01, 
    0.171999257518648893241583152515E+01, 
    0.232573248617385774545404479449E+01, 
    0.296716692790560324848896036355E+01, 
    0.366995037340445253472922383312E+01, 
    0.449999070730939155366438053053E+01, 
   -6.9956801237185402753248521473232E+00, 
   -6.2750787049428601427036567812530E+00, 
   -5.6739614446185883296332558789276E+00, 
   -5.1335955771123807045862968913996E+00, 
   -4.6315595063128599420667997654336E+00, 
   -4.1562717558181451724831352315314E+00, 
   -3.7007434032314694224497164589673E+00, 
   -3.2603207323135408104645401509648E+00, 
   -2.8316804533902054557015640151425E+00, 
   -2.4123177054804201051740184582119E+00, 
   -2.0002585489356389657975562598571E+00, 
   -1.5938858604721398261388419455550E+00, 
   -1.1918269983500464260821358649242E+00, 
   -0.79287697691530893968593032998830E+00, 
   -0.39594273647142311094670041663436E+00, 
    0.0000000000000000000000000000000E+00, 
    0.39594273647142311094670041663436E+00, 
    0.79287697691530893968593032998830E+00, 
    1.1918269983500464260821358649242E+00, 
    1.5938858604721398261388419455550E+00, 
    2.0002585489356389657975562598571E+00, 
    2.4123177054804201051740184582119E+00, 
    2.8316804533902054557015640151425E+00, 
    3.2603207323135408104645401509648E+00, 
    3.7007434032314694224497164589673E+00, 
    4.1562717558181451724831352315314E+00, 
    4.6315595063128599420667997654336E+00, 
    5.1335955771123807045862968913996E+00, 
    5.6739614446185883296332558789276E+00, 
    6.2750787049428601427036567812530E+00, 
    6.9956801237185402753248521473232E+00, 
   -10.435499877854168053468115427285E+00, 
   -9.8028759912974963635223935286507E+00, 
   -9.2792019543050391319404745506496E+00, 
   -8.8118581437284546442526628275570E+00, 
   -8.3807683451863219343010651043788E+00, 
   -7.9755950801420373181541806298501E+00, 
   -7.5901395198641066762479783194468E+00, 
   -7.2203167078889678461161324222529E+00, 
   -6.8632544331795368527353285876066E+00, 
   -6.5168348106821160605273395854042E+00, 
   -6.1794379922705969862418461787263E+00, 
   -5.8497884000810673462526582961482E+00, 
   -5.5268572526403031425047575122840E+00, 
   -5.2097979830408354861575136416263E+00, 
   -4.8979018644975742350745099214868E+00, 
   -4.5905665744435190229271294569091E+00, 
   -4.2872733352824404031727616199454E+00, 
   -3.9875699104197157485227052068068E+00, 
   -3.6910577000963465117322810559754E+00, 
   -3.3973817713303911852755941806287E+00, 
   -3.1062230279282566329138616746036E+00, 
   -2.8172919672837977750747135657355E+00, 
   -2.5303236304712010926855221718499E+00, 
   -2.2450734604812066298995918179330E+00, 
   -1.9613138583081485293922008411321E+00, 
   -1.6788312791720137520802800622638E+00, 
   -1.3974237486049625107570752063702E+00, 
   -1.1168987050996462690510970277840E+00, 
   -0.83707109558947615977737795461293E+00, 
   -0.55776166427908221668763665253822E+00, 
   -0.27879538567115223986687628627202E+00, 
    0.00000000000000000000000000000000E+00, 
    0.27879538567115223986687628627202E+00, 
    0.55776166427908221668763665253822E+00, 
    0.83707109558947615977737795461293E+00, 
    1.1168987050996462690510970277840E+00, 
    1.3974237486049625107570752063702E+00, 
    1.6788312791720137520802800622638E+00, 
    1.9613138583081485293922008411321E+00, 
    2.2450734604812066298995918179330E+00, 
    2.5303236304712010926855221718499E+00, 
    2.8172919672837977750747135657355E+00, 
    3.1062230279282566329138616746036E+00, 
    3.3973817713303911852755941806287E+00, 
    3.6910577000963465117322810559754E+00, 
    3.9875699104197157485227052068068E+00, 
    4.2872733352824404031727616199454E+00, 
    4.5905665744435190229271294569091E+00, 
    4.8979018644975742350745099214868E+00, 
    5.2097979830408354861575136416263E+00, 
    5.5268572526403031425047575122840E+00, 
    5.8497884000810673462526582961482E+00, 
    6.1794379922705969862418461787263E+00, 
    6.5168348106821160605273395854042E+00, 
    6.8632544331795368527353285876066E+00, 
    7.2203167078889678461161324222529E+00, 
    7.5901395198641066762479783194468E+00, 
    7.9755950801420373181541806298501E+00, 
    8.3807683451863219343010651043788E+00, 
    8.8118581437284546442526628275570E+00, 
    9.2792019543050391319404745506496E+00, 
    9.8028759912974963635223935286507E+00, 
    10.435499877854168053468115427285E+00, 
   -15.228338148167350978246954433464E+00, 
   -14.669595158833972632746354112896E+00, 
   -14.209085995284870755168244250887E+00, 
   -13.799722290211676634645246746673E+00, 
   -13.423518590070950062438258321855E+00, 
   -13.071208660474601901583995439649E+00, 
   -12.737235652415686338138003924072E+00, 
   -12.417939378869715805445879624069E+00, 
   -12.110749020947747600132123508132E+00, 
   -11.813772198267727195134584136191E+00, 
   -11.525565112572696599167888588564E+00, 
   -11.244994583785543445194384194300E+00, 
   -10.971150569840247423423040263881E+00, 
   -10.703288201027481347670940744690E+00, 
   -10.440787957772772867742591798027E+00, 
   -10.183127473450343888624126450357E+00, 
   -9.9298610495114250736847004273684E+00, 
   -9.6806044412474728038150712732737E+00, 
   -9.4350233389881650135019598506287E+00, 
   -9.1928244988460305715774195052527E+00, 
   -8.9537488108565404323807890169970E+00, 
   -8.7175658087076307363833999548548E+00, 
   -8.4840692689832473326097180339984E+00, 
   -8.2530736454457156579694124243888E+00, 
   -8.0244111514703375578594739796798E+00, 
   -7.7979293513870105420829120455591E+00, 
   -7.5734891556083454022834960763301E+00, 
   -7.3509631392269052701961258043733E+00, 
   -7.1302341220350710668064025713431E+00, 
   -6.9111939615465713197465633109366E+00, 
   -6.6937425208758294190074417381666E+00, 
   -6.4777867811645365448144903821487E+00, 
   -6.2632400742737354345609723857092E+00, 
   -6.0500214161419845694465474482388E+00, 
   -5.8380549248774187386601690807757E+00, 
   -5.6272693105464816659423455794909E+00, 
   -5.4175974259243240722848425872924E+00, 
   -5.2089758693153983587570258372239E+00, 
   -5.0013446320386360038520809107373E+00, 
   -4.7946467843764925009748509930857E+00, 
   -4.5888281947698372951606485031212E+00, 
   -4.3838372778464736294253744407459E+00, 
   -4.1796247675352031349421189892408E+00, 
   -3.9761435120673355916035814195920E+00, 
   -3.7733482881250526721004678400057E+00, 
   -3.5711956317782180447199756485249E+00, 
   -3.3696436841717397896643629240035E+00, 
   -3.1686520501953630191857798261495E+00, 
   -2.9681816685955910267761649521505E+00, 
   -2.7681946921824058801226545958892E+00, 
   -2.5686543769473501723144013022363E+00, 
   -2.3695249790490401080012474645702E+00, 
   -2.1707716587411506879498498083695E+00, 
   -1.9723603904195020079324743227565E+00, 
   -1.7742578780516791584676442103681E+00, 
   -1.5764314753267801315519597621879E+00, 
   -1.3788491099261778091441557053728E+00, 
   -1.1814792113700685848678583598423E+00, 
   -0.98429064194027277726568984213773E+00, 
   -0.78725263021825034151596831878971E+00, 
   -0.59033470680942102142230439346102E+00, 
   -0.39350664185130136568037826200185E+00, 
   -0.19673838392423251964272239737078E+00, 
    0.0000000000000000000000000000000E+00, 
    0.19673838392423251964272239737078E+00, 
    0.39350664185130136568037826200185E+00, 
    0.59033470680942102142230439346102E+00, 
    0.78725263021825034151596831878971E+00, 
    0.98429064194027277726568984213773E+00, 
    1.1814792113700685848678583598423E+00, 
    1.3788491099261778091441557053728E+00, 
    1.5764314753267801315519597621879E+00, 
    1.7742578780516791584676442103681E+00, 
    1.9723603904195020079324743227565E+00, 
    2.1707716587411506879498498083695E+00, 
    2.3695249790490401080012474645702E+00, 
    2.5686543769473501723144013022363E+00, 
    2.7681946921824058801226545958892E+00, 
    2.9681816685955910267761649521505E+00, 
    3.1686520501953630191857798261495E+00, 
    3.3696436841717397896643629240035E+00, 
    3.5711956317782180447199756485249E+00, 
    3.7733482881250526721004678400057E+00, 
    3.9761435120673355916035814195920E+00, 
    4.1796247675352031349421189892408E+00, 
    4.3838372778464736294253744407459E+00, 
    4.5888281947698372951606485031212E+00, 
    4.7946467843764925009748509930857E+00, 
    5.0013446320386360038520809107373E+00, 
    5.2089758693153983587570258372239E+00, 
    5.4175974259243240722848425872924E+00, 
    5.6272693105464816659423455794909E+00, 
    5.8380549248774187386601690807757E+00, 
    6.0500214161419845694465474482388E+00, 
    6.2632400742737354345609723857092E+00, 
    6.4777867811645365448144903821487E+00, 
    6.6937425208758294190074417381666E+00, 
    6.9111939615465713197465633109366E+00, 
    7.1302341220350710668064025713431E+00, 
    7.3509631392269052701961258043733E+00, 
    7.5734891556083454022834960763301E+00, 
    7.7979293513870105420829120455591E+00, 
    8.0244111514703375578594739796798E+00, 
    8.2530736454457156579694124243888E+00, 
    8.4840692689832473326097180339984E+00, 
    8.7175658087076307363833999548548E+00, 
    8.9537488108565404323807890169970E+00, 
    9.1928244988460305715774195052527E+00, 
    9.4350233389881650135019598506287E+00, 
    9.6806044412474728038150712732737E+00, 
    9.9298610495114250736847004273684E+00, 
    10.183127473450343888624126450357E+00, 
    10.440787957772772867742591798027E+00, 
    10.703288201027481347670940744690E+00, 
    10.971150569840247423423040263881E+00, 
    11.244994583785543445194384194300E+00, 
    11.525565112572696599167888588564E+00, 
    11.813772198267727195134584136191E+00, 
    12.110749020947747600132123508132E+00, 
    12.417939378869715805445879624069E+00, 
    12.737235652415686338138003924072E+00, 
    13.071208660474601901583995439649E+00, 
    13.423518590070950062438258321855E+00, 
    13.799722290211676634645246746673E+00, 
    14.209085995284870755168244250887E+00, 
    14.669595158833972632746354112896E+00, 
    15.228338148167350978246954433464E+00  
  };

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( grid_base[dim] < 0 )
    {
      cout << "\n";
      cout << "GH_ABSCISSA - Fatal error!\n";
      cout << "  Some base values are less than 0.\n";
      exit ( 1 );
    }
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 63 < grid_base[dim] )
    {
      cout << "\n";
      cout << "GH_ABSCISSA - Fatal error!\n";
      cout << "  Some base values are greater than 63.\n";
      exit ( 1 );
    }
  }

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level = i4_log_2 ( grid_base[dim] + 1 );

      pointer = skip[level] + ( grid_index[dim+point*dim_num] + grid_base[dim] );

      grid_point[dim+point*dim_num] = x[pointer];
    }
  }

  return;
}
//****************************************************************************80

double *gh_weights ( int order )

//****************************************************************************80
//
//  Purpose:
//
//    GH_WEIGHTS returns weights for certain Gauss-Hermite quadrature rules.
//
//  Discussion:
//
//    The allowed orders are 1, 3, 7, 15, 31, 63 and 127.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2007
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
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
//
//    Output, double WEIGHT[ORDER], the weights.
//    The weights are positive, symmetric and should sum to SQRT(PI).
//
{
  double *weight;

  weight = new double[order];

  if ( order == 1 )
  {
    weight[1-1] = 1.77245385090551602729816748334E+00;
  }
  else if ( order == 3 )
  {
    weight[1-1] = 0.295408975150919337883027913890E+00;
    weight[2-1] = 0.118163590060367735153211165556E+01;
    weight[3-1] = 0.295408975150919337883027913890E+00;
  }
  else if ( order == 7 )
  {
    weight[1-1] = 0.971781245099519154149424255939E-03;
    weight[2-1] = 0.545155828191270305921785688417E-01;
    weight[3-1] = 0.425607252610127800520317466666E+00;
    weight[4-1] = 0.810264617556807326764876563813E+00;
    weight[5-1] = 0.425607252610127800520317466666E+00;
    weight[6-1] = 0.545155828191270305921785688417E-01;
    weight[7-1] = 0.971781245099519154149424255939E-03;
  }
  else if ( order == 15 )
  {
    weight[1-1] =  0.152247580425351702016062666965E-08;
    weight[2-1] =  0.105911554771106663577520791055E-05;
    weight[3-1] =  0.100004441232499868127296736177E-03;
    weight[4-1] =  0.277806884291277589607887049229E-02;
    weight[5-1] =  0.307800338725460822286814158758E-01;
    weight[6-1] =  0.158488915795935746883839384960E+00;
    weight[7-1] =  0.412028687498898627025891079568E+00;
    weight[8-1] =  0.564100308726417532852625797340E+00;
    weight[9-1] =  0.412028687498898627025891079568E+00;
    weight[10-1] = 0.158488915795935746883839384960E+00;
    weight[11-1] = 0.307800338725460822286814158758E-01;
    weight[12-1] = 0.277806884291277589607887049229E-02;
    weight[13-1] = 0.100004441232499868127296736177E-03;
    weight[14-1] = 0.105911554771106663577520791055E-05;
    weight[15-1] = 0.152247580425351702016062666965E-08;
  }
  else if ( order == 31 )
  {
    weight[  1-1] =   0.46189683944498305857470556847735E-21;
    weight[  2-1] =   0.51106090079112519643027197715274E-17;
    weight[  3-1] =   0.58995564987355133075257722133966E-14;
    weight[  4-1] =   0.18603735214463569590294465062239E-11;
    weight[  5-1] =   0.23524920032013205739850619940094E-09;
    weight[  6-1] =   0.14611988344865057576066495091513E-07;
    weight[  7-1] =   0.50437125589241034841778074689627E-06;
    weight[  8-1] =   0.10498602757642934202945441341697E-04;
    weight[  9-1] =   0.13952090395003623854995664958146E-03;
    weight[ 10-1] =   0.12336833073030489880608311394968E-02;
    weight[ 11-1] =   0.74827999140119116765002499116934E-02;
    weight[ 12-1] =   0.31847230731201222775249585776902E-01;
    weight[ 13-1] =   0.96717948160569462991143316029341E-01;
    weight[ 14-1] =   0.21213278866810461318136114862419E+00;
    weight[ 15-1] =   0.33877265789305344906000174083214E+00;
    weight[ 16-1] =   0.39577855609737786462923720809676E+00;
    weight[ 17-1] =   0.33877265789305344906000174083214E+00;
    weight[ 18-1] =   0.21213278866810461318136114862419E+00;
    weight[ 19-1] =   0.96717948160569462991143316029341E-01;
    weight[ 20-1] =   0.31847230731201222775249585776902E-01;
    weight[ 21-1] =   0.74827999140119116765002499116934E-02;
    weight[ 22-1] =   0.12336833073030489880608311394968E-02;
    weight[ 23-1] =   0.13952090395003623854995664958146E-03;
    weight[ 24-1] =   0.10498602757642934202945441341697E-04;
    weight[ 25-1] =   0.50437125589241034841778074689627E-06;
    weight[ 26-1] =   0.14611988344865057576066495091513E-07;
    weight[ 27-1] =   0.23524920032013205739850619940094E-09;
    weight[ 28-1] =   0.18603735214463569590294465062239E-11;
    weight[ 29-1] =   0.58995564987355133075257722133966E-14;
    weight[ 30-1] =   0.51106090079112519643027197715274E-17;
    weight[ 31-1] =   0.46189683944498305857470556847735E-21;
  }
  else if ( order == 63 )
  {
    weight[  1-1] =   0.37099206434787551197827130470031E-47;
    weight[  2-1] =   0.10400778615192299534481914814892E-41;
    weight[  3-1] =   0.19796804708258311251124226474396E-37;
    weight[  4-1] =   0.84687478191640015120141181138947E-34;
    weight[  5-1] =   0.13071305930779945903630127634063E-30;
    weight[  6-1] =   0.93437837175367456929765381518998E-28;
    weight[  7-1] =   0.36027426635173044862245783257252E-25;
    weight[  8-1] =   0.82963863115951789374753323156164E-23;
    weight[  9-1] =   0.12266629909105281472971700203949E-20;
    weight[ 10-1] =   0.12288435628797061539461585325494E-18;
    weight[ 11-1] =   0.86925536958188009075932426691516E-17;
    weight[ 12-1] =   0.44857058689176221240330804981619E-15;
    weight[ 13-1] =   0.17335817955735154599902643794700E-13;
    weight[ 14-1] =   0.51265062385038307838565047455223E-12;
    weight[ 15-1] =   0.11808921844532942490513037158404E-10;
    weight[ 16-1] =   0.21508698297808025739828859845140E-09;
    weight[ 17-1] =   0.31371929535285447801497640621672E-08;
    weight[ 18-1] =   0.37041625984781705796752840204084E-07;
    weight[ 19-1] =   0.35734732949879669663960738150956E-06;
    weight[ 20-1] =   0.28393114498380927832990899215541E-05;
    weight[ 21-1] =   0.18709113003730498008961134765721E-04;
    weight[ 22-1] =   0.10284880800653635546698378640623E-03;
    weight[ 23-1] =   0.47411702610173128107201781718693E-03;
    weight[ 24-1] =   0.18409222622384813438539657470055E-02;
    weight[ 25-1] =   0.60436044551187631655712178246467E-02;
    weight[ 26-1] =   0.16829299199599730926458559757600E-01;
    weight[ 27-1] =   0.39858264027692992170237391875317E-01;
    weight[ 28-1] =   0.80467087993950415219587554532823E-01;
    weight[ 29-1] =   0.13871950817615293377792092082674E+00;
    weight[ 30-1] =   0.20448695346833761570957197160475E+00;
    weight[ 31-1] =   0.25799889943058042204920467417642E+00;
    weight[ 32-1] =   0.27876694884838411919175686949858E+00;
    weight[ 33-1] =   0.25799889943058042204920467417642E+00;
    weight[ 34-1] =   0.20448695346833761570957197160475E+00;
    weight[ 35-1] =   0.13871950817615293377792092082674E+00;
    weight[ 36-1] =   0.80467087993950415219587554532823E-01;
    weight[ 37-1] =   0.39858264027692992170237391875317E-01;
    weight[ 38-1] =   0.16829299199599730926458559757600E-01;
    weight[ 39-1] =   0.60436044551187631655712178246467E-02;
    weight[ 40-1] =   0.18409222622384813438539657470055E-02;
    weight[ 41-1] =   0.47411702610173128107201781718693E-03;
    weight[ 42-1] =   0.10284880800653635546698378640623E-03;
    weight[ 43-1] =   0.18709113003730498008961134765721E-04;
    weight[ 44-1] =   0.28393114498380927832990899215541E-05;
    weight[ 45-1] =   0.35734732949879669663960738150956E-06;
    weight[ 46-1] =   0.37041625984781705796752840204084E-07;
    weight[ 47-1] =   0.31371929535285447801497640621672E-08;
    weight[ 48-1] =   0.21508698297808025739828859845140E-09;
    weight[ 49-1] =   0.11808921844532942490513037158404E-10;
    weight[ 50-1] =   0.51265062385038307838565047455223E-12;
    weight[ 51-1] =   0.17335817955735154599902643794700E-13;
    weight[ 52-1] =   0.44857058689176221240330804981619E-15;
    weight[ 53-1] =   0.86925536958188009075932426691516E-17;
    weight[ 54-1] =   0.12288435628797061539461585325494E-18;
    weight[ 55-1] =   0.12266629909105281472971700203949E-20;
    weight[ 56-1] =   0.82963863115951789374753323156164E-23;
    weight[ 57-1] =   0.36027426635173044862245783257252E-25;
    weight[ 58-1] =   0.93437837175367456929765381518998E-28;
    weight[ 59-1] =   0.13071305930779945903630127634063E-30;
    weight[ 60-1] =   0.84687478191640015120141181138947E-34;
    weight[ 61-1] =   0.19796804708258311251124226474396E-37;
    weight[ 62-1] =   0.10400778615192299534481914814892E-41;
    weight[ 63-1] =   0.37099206434787551197827130470031E-47;
  }
  else if ( order == 127 )
  {
    weight[  1-1] =   0.12504497577050595552677230002883E-100;
    weight[  2-1] =   0.17272798059419131415318615789672E-93;
    weight[  3-1] =   0.89321681571986548608031150791499E-88;
    weight[  4-1] =   0.77306185240893578449625186483810E-83;
    weight[  5-1] =   0.20143957652648255497735460506196E-78;
    weight[  6-1] =   0.21503714733610239701351039429345E-74;
    weight[  7-1] =   0.11341924208594594813715533569504E-70;
    weight[  8-1] =   0.33489139011795051950683388483136E-67;
    weight[  9-1] =   0.60486548964016681064424451668405E-64;
    weight[ 10-1] =   0.71375092946352177824971347343892E-61;
    weight[ 11-1] =   0.57884563374885556636801095624030E-58;
    weight[ 12-1] =   0.33581166223858230300409326551248E-55;
    weight[ 13-1] =   0.14394641949253923568603163698953E-52;
    weight[ 14-1] =   0.46821808383216117724080263903889E-50;
    weight[ 15-1] =   0.11817054440684264071348471955361E-47;
    weight[ 16-1] =   0.23581659156008927203181682045005E-45;
    weight[ 17-1] =   0.37814427940797540210712758405540E-43;
    weight[ 18-1] =   0.49411031115771638145610738414006E-41;
    weight[ 19-1] =   0.53255303775425059266087298458297E-39;
    weight[ 20-1] =   0.47854390680131484999315199332765E-37;
    weight[ 21-1] =   0.36191883445952356128627543209554E-35;
    weight[ 22-1] =   0.23232083386343554805352497446119E-33;
    weight[ 23-1] =   0.12753331411008716683688974281454E-31;
    weight[ 24-1] =   0.60277753850758742112436095241270E-30;
    weight[ 25-1] =   0.24679773241777200207460855084439E-28;
    weight[ 26-1] =   0.88019567691698482573264198727415E-27;
    weight[ 27-1] =   0.27482489212040561315005725890593E-25;
    weight[ 28-1] =   0.75468218903085486125222816438456E-24;
    weight[ 29-1] =   0.18303134636280466270545996891835E-22;
    weight[ 30-1] =   0.39355990860860813085582448449811E-21;
    weight[ 31-1] =   0.75293161638581191068419292570042E-20;
    weight[ 32-1] =   0.12857997786722855037584105682618E-18;
    weight[ 33-1] =   0.19659326888445857792541925311450E-17;
    weight[ 34-1] =   0.26986511907214101894995783364250E-16;
    weight[ 35-1] =   0.33344414303198856330118301113874E-15;
    weight[ 36-1] =   0.37173303125150639885726463109574E-14;
    weight[ 37-1] =   0.37473954472839737091885387788983E-13;
    weight[ 38-1] =   0.34230094493397259538669512076007E-12;
    weight[ 39-1] =   0.28385303724993373166810860630552E-11;
    weight[ 40-1] =   0.21406920290454669208938772802828E-10;
    weight[ 41-1] =   0.14706331273431716244229273183839E-09;
    weight[ 42-1] =   0.92173940967434659264335883218167E-09;
    weight[ 43-1] =   0.52781663936972714041837056042506E-08;
    weight[ 44-1] =   0.27650497044951117835905283127679E-07;
    weight[ 45-1] =   0.13267855842539464770913063113371E-06;
    weight[ 46-1] =   0.58380944276113062188573331195042E-06;
    weight[ 47-1] =   0.23581561724775629112332165335800E-05;
    weight[ 48-1] =   0.87524468034280444703919485644809E-05;
    weight[ 49-1] =   0.29876790535909012274846532159647E-04;
    weight[ 50-1] =   0.93874435720072545206729594267039E-04;
    weight[ 51-1] =   0.27170762627931172053444716883938E-03;
    weight[ 52-1] =   0.72493929742498358979684249380921E-03;
    weight[ 53-1] =   0.17841208326763432884316727108264E-02;
    weight[ 54-1] =   0.40524855186046131499765636276283E-02;
    weight[ 55-1] =   0.85000263041544110385806705526917E-02;
    weight[ 56-1] =   0.16471142241609687824005585301760E-01;
    weight[ 57-1] =   0.29499296248213632269675010319119E-01;
    weight[ 58-1] =   0.48847387114300011006959603975676E-01;
    weight[ 59-1] =   0.74807989768583731416517226905270E-01;
    weight[ 60-1] =   0.10598520508090929403834368934301E+00;
    weight[ 61-1] =   0.13893945309051540832066283010510E+00;
    weight[ 62-1] =   0.16856236074207929740526975049765E+00;
    weight[ 63-1] =   0.18927849580120432177170145550076E+00;
    weight[ 64-1] =   0.19673340688823289786163676995151E+00;
    weight[ 65-1] =   0.18927849580120432177170145550076E+00;
    weight[ 66-1] =   0.16856236074207929740526975049765E+00;
    weight[ 67-1] =   0.13893945309051540832066283010510E+00;
    weight[ 68-1] =   0.10598520508090929403834368934301E+00;
    weight[ 69-1] =   0.74807989768583731416517226905270E-01;
    weight[ 70-1] =   0.48847387114300011006959603975676E-01;
    weight[ 71-1] =   0.29499296248213632269675010319119E-01;
    weight[ 72-1] =   0.16471142241609687824005585301760E-01;
    weight[ 73-1] =   0.85000263041544110385806705526917E-02;
    weight[ 74-1] =   0.40524855186046131499765636276283E-02;
    weight[ 75-1] =   0.17841208326763432884316727108264E-02;
    weight[ 76-1] =   0.72493929742498358979684249380921E-03;
    weight[ 77-1] =   0.27170762627931172053444716883938E-03;
    weight[ 78-1] =   0.93874435720072545206729594267039E-04;
    weight[ 79-1] =   0.29876790535909012274846532159647E-04;
    weight[ 80-1] =   0.87524468034280444703919485644809E-05;
    weight[ 81-1] =   0.23581561724775629112332165335800E-05;
    weight[ 82-1] =   0.58380944276113062188573331195042E-06;
    weight[ 83-1] =   0.13267855842539464770913063113371E-06;
    weight[ 84-1] =   0.27650497044951117835905283127679E-07;
    weight[ 85-1] =   0.52781663936972714041837056042506E-08;
    weight[ 86-1] =   0.92173940967434659264335883218167E-09;
    weight[ 87-1] =   0.14706331273431716244229273183839E-09;
    weight[ 88-1] =   0.21406920290454669208938772802828E-10;
    weight[ 89-1] =   0.28385303724993373166810860630552E-11;
    weight[ 90-1] =   0.34230094493397259538669512076007E-12;
    weight[ 91-1] =   0.37473954472839737091885387788983E-13;
    weight[ 92-1] =   0.37173303125150639885726463109574E-14;
    weight[ 93-1] =   0.33344414303198856330118301113874E-15;
    weight[ 94-1] =   0.26986511907214101894995783364250E-16;
    weight[ 95-1] =   0.19659326888445857792541925311450E-17;
    weight[ 96-1] =   0.12857997786722855037584105682618E-18;
    weight[ 97-1] =   0.75293161638581191068419292570042E-20;
    weight[ 98-1] =   0.39355990860860813085582448449811E-21;
    weight[ 99-1] =   0.18303134636280466270545996891835E-22;
    weight[100-1] =   0.75468218903085486125222816438456E-24;
    weight[101-1] =   0.27482489212040561315005725890593E-25;
    weight[102-1] =   0.88019567691698482573264198727415E-27;
    weight[103-1] =   0.24679773241777200207460855084439E-28;
    weight[104-1] =   0.60277753850758742112436095241270E-30;
    weight[105-1] =   0.12753331411008716683688974281454E-31;
    weight[106-1] =   0.23232083386343554805352497446119E-33;
    weight[107-1] =   0.36191883445952356128627543209554E-35;
    weight[108-1] =   0.47854390680131484999315199332765E-37;
    weight[109-1] =   0.53255303775425059266087298458297E-39;
    weight[110-1] =   0.49411031115771638145610738414006E-41;
    weight[111-1] =   0.37814427940797540210712758405540E-43;
    weight[112-1] =   0.23581659156008927203181682045005E-45;
    weight[113-1] =   0.11817054440684264071348471955361E-47;
    weight[114-1] =   0.46821808383216117724080263903889E-50;
    weight[115-1] =   0.14394641949253923568603163698953E-52;
    weight[116-1] =   0.33581166223858230300409326551248E-55;
    weight[117-1] =   0.57884563374885556636801095624030E-58;
    weight[118-1] =   0.71375092946352177824971347343892E-61;
    weight[119-1] =   0.60486548964016681064424451668405E-64;
    weight[120-1] =   0.33489139011795051950683388483136E-67;
    weight[121-1] =   0.11341924208594594813715533569504E-70;
    weight[122-1] =   0.21503714733610239701351039429345E-74;
    weight[123-1] =   0.20143957652648255497735460506196E-78;
    weight[124-1] =   0.77306185240893578449625186483810E-83;
    weight[125-1] =   0.89321681571986548608031150791499E-88;
    weight[126-1] =   0.17272798059419131415318615789672E-93;
    weight[127-1] =   0.12504497577050595552677230002883E-100;
  }
  else
  {
    cout << "\n";
    cout << "GH_WEIGHTS - Fatal error!\n";
    cout << "  Illegal value of ORDER = " << order << "\n";
    cout << "  Legal values are 1, 3, 7, 15, 31, 63 and 127.\n";
    exit ( 1 );
  }
  return weight;
}
//****************************************************************************80

void gl_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    GL_ABSCISSA sets abscissas for multidimensional Gauss-Legendre quadrature.
//
//  Discussion:
//
//    The "nesting" as it occurs for Gauss-Legendre sparse grids simply
//    involves the use of a specified set of permissible orders for the
//    rule.  
//
//    The X array lists the (complete) Gauss-Legendre abscissas for rules 
//    of order 1, 3, 7, 15, 31, 63 or 127, in order. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the index of the abscissa
//    from the Gauss-Legendre rule, for each dimension and point.
//
//    Input, int GRID_BASE[DIM_NUM], the number of points used in the 
//    Gauss-Legendre rule for a given dimension.
//
//    Output, double GRID_POINT[DIM_NUM], the grid points of
//    Gauss-Legendre abscissas.
//
{
  int dim;
  int level;
  int point;
  int pointer;
  int skip[8] = { 0, 1, 4, 11, 26, 57, 120, 247 };
  double x[247] = {
       0.0E+00, 
     - 0.774596669241483377035853079956E+00, 
       0.0E+00, 
       0.774596669241483377035853079956E+00, 
     - 0.949107912342758524526189684048E+00, 
     - 0.741531185599394439863864773281E+00, 
     - 0.405845151377397166906606412077E+00, 
       0.0E+00, 
       0.405845151377397166906606412077E+00, 
       0.741531185599394439863864773281E+00, 
       0.949107912342758524526189684048E+00, 
      - 0.987992518020485428489565718587E+00, 
      - 0.937273392400705904307758947710E+00, 
      - 0.848206583410427216200648320774E+00, 
      - 0.724417731360170047416186054614E+00, 
      - 0.570972172608538847537226737254E+00, 
      - 0.394151347077563369897207370981E+00, 
      - 0.201194093997434522300628303395E+00, 
        0.0E+00, 
        0.201194093997434522300628303395E+00, 
       0.394151347077563369897207370981E+00, 
       0.570972172608538847537226737254E+00, 
       0.724417731360170047416186054614E+00, 
       0.848206583410427216200648320774E+00, 
       0.937273392400705904307758947710E+00, 
       0.987992518020485428489565718587E+00, 
      -0.99708748181947707454263838179654,    
      -0.98468590966515248400211329970113,    
      -0.96250392509294966178905249675943,    
      -0.93075699789664816495694576311725,    
      -0.88976002994827104337419200908023,    
      -0.83992032014626734008690453594388,    
      -0.78173314841662494040636002019484,    
      -0.71577678458685328390597086536649,    
      -0.64270672292426034618441820323250,    
      -0.56324916140714926272094492359516,    
      -0.47819378204490248044059403935649,    
      -0.38838590160823294306135146128752,    
      -0.29471806998170161661790389767170,    
      -0.19812119933557062877241299603283,    
      -0.99555312152341520325174790118941E-01,
       0.00000000000000000000000000000000,
       0.99555312152341520325174790118941E-01,
       0.19812119933557062877241299603283,    
       0.29471806998170161661790389767170,    
       0.38838590160823294306135146128752,    
       0.47819378204490248044059403935649,    
       0.56324916140714926272094492359516,
       0.64270672292426034618441820323250,    
       0.71577678458685328390597086536649,    
       0.78173314841662494040636002019484,    
       0.83992032014626734008690453594388,    
       0.88976002994827104337419200908023,    
       0.93075699789664816495694576311725,    
       0.96250392509294966178905249675943,    
       0.98468590966515248400211329970113,   
       0.99708748181947707454263838179654,    
      -0.99928298402912378050701628988630E+00,     
      -0.99622401277797010860209018267357E+00,     
      -0.99072854689218946681089469460884E+00,     
      -0.98280881059372723486251140727639E+00,     
      -0.97248403469757002280196067864927E+00,     
      -0.95977944975894192707035416626398E+00,     
      -0.94472613404100980296637531962798E+00,     
      -0.92736092062184320544703138132518E+00,     
      -0.90772630277853155803695313291596E+00,     
      -0.88587032850785342629029845731337E+00,     
      -0.86184648236412371953961183943106E+00,     
      -0.83571355431950284347180776961571E+00,     
      -0.80753549577345676005146598636324E+00,     
      -0.77738126299037233556333018991104E+00,     
      -0.74532464831784741782932166103759E+00,     
      -0.71144409958484580785143153770401E+00,     
      -0.67582252811498609013110331596954E+00,     
      -0.63854710582136538500030695387338E+00,     
      -0.59970905187762523573900892686880E+00,     
      -0.55940340948628501326769780007005E+00,     
      -0.51772881329003324812447758452632E+00,     
      -0.47478724799480439992221230985149E+00,     
      -0.43068379879511160066208893391863E+00,     
      -0.38552639421224789247761502227440E+00,     
      -0.33942554197458440246883443159432E+00,     
      -0.29249405858625144003615715555067E+00,     
      -0.24484679324595336274840459392483E+00,     
      -0.19660034679150668455762745706572E+00,     
      -0.14787278635787196856983909655297E+00,     
      -0.98783356446945279529703669453922E-01, 
      -0.49452187116159627234233818051808E-01, 
       0.00000000000000000000000000000000E+00,     
       0.49452187116159627234233818051808E-01, 
       0.98783356446945279529703669453922E-01, 
       0.14787278635787196856983909655297E+00,     
       0.19660034679150668455762745706572E+00,     
       0.24484679324595336274840459392483E+00,     
       0.29249405858625144003615715555067E+00,     
       0.33942554197458440246883443159432E+00,     
       0.38552639421224789247761502227440E+00,     
       0.43068379879511160066208893391863E+00,     
       0.47478724799480439992221230985149E+00,     
       0.51772881329003324812447758452632E+00,     
       0.55940340948628501326769780007005E+00,     
       0.59970905187762523573900892686880E+00,     
       0.63854710582136538500030695387338E+00,     
       0.67582252811498609013110331596954E+00,     
       0.71144409958484580785143153770401E+00,     
       0.74532464831784741782932166103759E+00,     
       0.77738126299037233556333018991104E+00,     
       0.80753549577345676005146598636324E+00,     
       0.83571355431950284347180776961571E+00,     
       0.86184648236412371953961183943106E+00,     
       0.88587032850785342629029845731337E+00,     
       0.90772630277853155803695313291596E+00,     
       0.92736092062184320544703138132518E+00,     
       0.94472613404100980296637531962798E+00,     
       0.95977944975894192707035416626398E+00,     
       0.97248403469757002280196067864927E+00,     
       0.98280881059372723486251140727639E+00,     
       0.99072854689218946681089469460884E+00,     
       0.99622401277797010860209018267357E+00,  
       0.99928298402912378050701628988630E+00
  };

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( grid_base[dim] < 0 )
    {
      cout << "\n";
      cout << "GL_ABSCISSA - Fatal error!\n";
      cout << "  Some base values are less than 0.\n";
      exit ( 1 );
    }
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 63 < grid_base[dim] )
    {
      cout << "\n";
      cout << "GL_ABSCISSA - Fatal error!\n";
      cout << "  Some base values are greater than 63.\n";
      exit ( 1 );
    }
  }

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level = i4_log_2 ( grid_base[dim] + 1 );

      pointer = skip[level] + ( grid_index[dim+point*dim_num] + grid_base[dim] );

      grid_point[dim+point*dim_num] = x[pointer];
    }
  }

  return;
}
//****************************************************************************80

double *gl_weights ( int order )

//****************************************************************************80
//
//  Purpose:
//
//    GL_WEIGHTS returns weights for certain Gauss-Legendre quadrature rules.
//
//  Discussion:
//
//    The allowed orders are 1, 3, 7, 15, 31, 63 and 127.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2007
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
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
//
//    Output, double WEIGHT[ORDER], the weights.
//    The weights are positive, symmetric and should sum to 2.
//
{
  double *weight;

  weight = new double[order];

  if ( order == 1 )
  {
    weight[1-1] = 2.0E+00;
  }
  else if ( order == 3 )
  {
    weight[1-1] = 5.0E+00 / 9.0E+00;
    weight[2-1] = 8.0E+00 / 9.0E+00;
    weight[3-1] = 5.0E+00 / 9.0E+00;
  }
  else if ( order == 7 )
  {
    weight[1-1] = 0.129484966168869693270611432679E+00;
    weight[2-1] = 0.279705391489276667901467771424E+00;
    weight[3-1] = 0.381830050505118944950369775489E+00;
    weight[4-1] = 0.417959183673469387755102040816E+00;
    weight[5-1] = 0.381830050505118944950369775489E+00;
    weight[6-1] = 0.279705391489276667901467771424E+00;
    weight[7-1] = 0.129484966168869693270611432679E+00;
  }
  else if ( order == 15 )
  {
    weight[1-1] =  0.307532419961172683546283935772E-01;
    weight[2-1] =  0.703660474881081247092674164507E-01;
    weight[3-1] =  0.107159220467171935011869546686E+00;
    weight[4-1] =  0.139570677926154314447804794511E+00;
    weight[5-1] =  0.166269205816993933553200860481E+00;
    weight[6-1] =  0.186161000015562211026800561866E+00;
    weight[7-1] =  0.198431485327111576456118326444E+00;
    weight[8-1] =  0.202578241925561272880620199968E+00;
    weight[9-1] =  0.198431485327111576456118326444E+00;
    weight[10-1] = 0.186161000015562211026800561866E+00;
    weight[11-1] = 0.166269205816993933553200860481E+00;
    weight[12-1] = 0.139570677926154314447804794511E+00;
    weight[13-1] = 0.107159220467171935011869546686E+00;
    weight[14-1] = 0.703660474881081247092674164507E-01;
    weight[15-1] = 0.307532419961172683546283935772E-01;
  }
  else if ( order == 31 )
  {
    weight[ 1-1] =   0.74708315792487746093913218970494E-02;
    weight[ 2-1] =   0.17318620790310582463552990782414E-01;
    weight[ 3-1] =   0.27009019184979421800608642617676E-01;
    weight[ 4-1] =   0.36432273912385464024392008749009E-01;
    weight[ 5-1] =   0.45493707527201102902315857856518E-01;
    weight[ 6-1] =   0.54103082424916853711666259085477E-01;
    weight[ 7-1] =   0.62174786561028426910343543686657E-01;
    weight[ 8-1] =   0.69628583235410366167756126255124E-01;
    weight[ 9-1] =   0.76390386598776616426357674901331E-01;
    weight[10-1] =   0.82392991761589263903823367431962E-01;
    weight[11-1] =   0.87576740608477876126198069695333E-01;
    weight[12-1] =   0.91890113893641478215362871607150E-01;
    weight[13-1] =   0.95290242912319512807204197487597E-01;
    weight[14-1] =   0.97743335386328725093474010978997E-01;
    weight[15-1] =   0.99225011226672307874875514428615E-01;
    weight[16-1] =   0.99720544793426451427533833734349E-01;
    weight[17-1] =   0.99225011226672307874875514428615E-01;
    weight[18-1] =   0.97743335386328725093474010978997E-01;
    weight[19-1] =   0.95290242912319512807204197487597E-01;
    weight[20-1] =   0.91890113893641478215362871607150E-01;
    weight[21-1] =   0.87576740608477876126198069695333E-01;
    weight[22-1] =   0.82392991761589263903823367431962E-01;
    weight[23-1] =   0.76390386598776616426357674901331E-01;
    weight[24-1] =   0.69628583235410366167756126255124E-01;
    weight[25-1] =   0.62174786561028426910343543686657E-01;
    weight[26-1] =   0.54103082424916853711666259085477E-01;
    weight[27-1] =   0.45493707527201102902315857856518E-01;
    weight[28-1] =   0.36432273912385464024392008749009E-01;
    weight[29-1] =   0.27009019184979421800608642617676E-01;
    weight[30-1] =   0.17318620790310582463552990782414E-01;
    weight[31-1] =   0.74708315792487746093913218970494E-02;
  }
  else if ( order == 63 )
  {
    weight[ 1-1] =   0.18398745955770837880499331680577E-02;
    weight[ 2-1] =   0.42785083468637618661951422543371E-02;
    weight[ 3-1] =   0.67102917659601362519069109850892E-02;
    weight[ 4-1] =   0.91259686763266563540586445877022E-02;
    weight[ 5-1] =   0.11519376076880041750750606118707E-01;
    weight[ 6-1] =   0.13884612616115610824866086365937E-01;
    weight[ 7-1] =   0.16215878410338338882283672974995E-01;
    weight[ 8-1] =   0.18507464460161270409260545805144E-01;
    weight[ 9-1] =   0.20753761258039090775341953421471E-01;
    weight[10-1] =   0.22949271004889933148942319561770E-01;
    weight[11-1] =   0.25088620553344986618630138068443E-01;
    weight[12-1] =   0.27166574359097933225189839439413E-01;
    weight[13-1] =   0.29178047208280526945551502154029E-01;
    weight[14-1] =   0.31118116622219817508215988557189E-01;
    weight[15-1] =   0.32982034883779341765683179672459E-01;
    weight[16-1] =   0.34765240645355877697180504642788E-01;
    weight[17-1] =   0.36463370085457289630452409787542E-01;
    weight[18-1] =   0.38072267584349556763638324927889E-01;
    weight[19-1] =   0.39587995891544093984807928149202E-01;
    weight[20-1] =   0.41006845759666398635110037009072E-01;
    weight[21-1] =   0.42325345020815822982505485403028E-01;
    weight[22-1] =   0.43540267083027590798964315704401E-01;
    weight[23-1] =   0.44648638825941395370332669516813E-01;
    weight[24-1] =   0.45647747876292608685885992608542E-01;
    weight[25-1] =   0.46535149245383696510395418746953E-01;
    weight[26-1] =   0.47308671312268919080604988338844E-01;
    weight[27-1] =   0.47966421137995131411052756195132E-01;
    weight[28-1] =   0.48506789097883847864090099145802E-01;
    weight[29-1] =   0.48928452820511989944709361549215E-01;
    weight[30-1] =   0.49230380423747560785043116988145E-01;
    weight[31-1] =   0.49411833039918178967039646116705E-01;
    weight[32-1] =   0.49472366623931020888669360420926E-01;
    weight[33-1] =   0.49411833039918178967039646116705E-01;
    weight[34-1] =   0.49230380423747560785043116988145E-01;
    weight[35-1] =   0.48928452820511989944709361549215E-01;
    weight[36-1] =   0.48506789097883847864090099145802E-01;
    weight[37-1] =   0.47966421137995131411052756195132E-01;
    weight[38-1] =   0.47308671312268919080604988338844E-01;
    weight[39-1] =   0.46535149245383696510395418746953E-01;
    weight[40-1] =   0.45647747876292608685885992608542E-01;
    weight[41-1] =   0.44648638825941395370332669516813E-01;
    weight[42-1] =   0.43540267083027590798964315704401E-01;
    weight[43-1] =   0.42325345020815822982505485403028E-01;
    weight[44-1] =   0.41006845759666398635110037009072E-01;
    weight[45-1] =   0.39587995891544093984807928149202E-01;
    weight[46-1] =   0.38072267584349556763638324927889E-01;
    weight[47-1] =   0.36463370085457289630452409787542E-01;
    weight[48-1] =   0.34765240645355877697180504642788E-01;
    weight[49-1] =   0.32982034883779341765683179672459E-01;
    weight[50-1] =   0.31118116622219817508215988557189E-01;
    weight[51-1] =   0.29178047208280526945551502154029E-01;
    weight[52-1] =   0.27166574359097933225189839439413E-01;
    weight[53-1] =   0.25088620553344986618630138068443E-01;
    weight[54-1] =   0.22949271004889933148942319561770E-01;
    weight[55-1] =   0.20753761258039090775341953421471E-01;
    weight[56-1] =   0.18507464460161270409260545805144E-01;
    weight[57-1] =   0.16215878410338338882283672974995E-01;
    weight[58-1] =   0.13884612616115610824866086365937E-01;
    weight[59-1] =   0.11519376076880041750750606118707E-01;
    weight[60-1] =   0.91259686763266563540586445877022E-02;
    weight[61-1] =   0.67102917659601362519069109850892E-02;
    weight[62-1] =   0.42785083468637618661951422543371E-02;
    weight[63-1] =   0.18398745955770837880499331680577E-02;
  }
  else if ( order == 127 )
  {
    weight[  1-1] =   0.45645726109586654495731936146574E-03;
    weight[  2-1] =   0.10622766869538486959954760554099E-02;
    weight[  3-1] =   0.16683488125171936761028811985672E-02;
    weight[  4-1] =   0.22734860707492547802810838362671E-02;
    weight[  5-1] =   0.28772587656289004082883197417581E-02;
    weight[  6-1] =   0.34792893810051465908910894094105E-02;
    weight[  7-1] =   0.40792095178254605327114733456293E-02;
    weight[  8-1] =   0.46766539777779034772638165662478E-02;
    weight[  9-1] =   0.52712596565634400891303815906251E-02;
    weight[ 10-1] =   0.58626653903523901033648343751367E-02;
    weight[ 11-1] =   0.64505120486899171845442463868748E-02;
    weight[ 12-1] =   0.70344427036681608755685893032552E-02;
    weight[ 13-1] =   0.76141028256526859356393930849227E-02;
    weight[ 14-1] =   0.81891404887415730817235884718726E-02;
    weight[ 15-1] =   0.87592065795403145773316804234385E-02;
    weight[ 16-1] =   0.93239550065309714787536985834029E-02;
    weight[ 17-1] =   0.98830429087554914716648010899606E-02;
    weight[ 18-1] =   0.10436130863141005225673171997668E-01;
    weight[ 19-1] =   0.10982883090068975788799657376065E-01;
    weight[ 20-1] =   0.11522967656921087154811609734510E-01;
    weight[ 21-1] =   0.12056056679400848183529562144697E-01;
    weight[ 22-1] =   0.12581826520465013101514365424172E-01;
    weight[ 23-1] =   0.13099957986718627426172681912499E-01;
    weight[ 24-1] =   0.13610136522139249906034237533759E-01;
    weight[ 25-1] =   0.14112052399003395774044161633613E-01;
    weight[ 26-1] =   0.14605400905893418351737288078952E-01;
    weight[ 27-1] =   0.15089882532666922992635733981431E-01;
    weight[ 28-1] =   0.15565203152273955098532590262975E-01;
    weight[ 29-1] =   0.16031074199309941802254151842763E-01;
    weight[ 30-1] =   0.16487212845194879399346060358146E-01;
    weight[ 31-1] =   0.16933342169871654545878815295200E-01;
    weight[ 32-1] =   0.17369191329918731922164721250350E-01;
    weight[ 33-1] =   0.17794495722974774231027912900351E-01;
    weight[ 34-1] =   0.18208997148375106468721469154479E-01;
    weight[ 35-1] =   0.18612443963902310429440419898958E-01;
    weight[ 36-1] =   0.19004591238555646611148901044533E-01;
    weight[ 37-1] =   0.19385200901246454628112623489471E-01;
    weight[ 38-1] =   0.19754041885329183081815217323169E-01;
    weight[ 39-1] =   0.20110890268880247225644623956287E-01;
    weight[ 40-1] =   0.20455529410639508279497065713301E-01;
    weight[ 41-1] =   0.20787750081531811812652137291250E-01;
    weight[ 42-1] =   0.21107350591688713643523847921658E-01;
    weight[ 43-1] =   0.21414136912893259295449693233545E-01;
    weight[ 44-1] =   0.21707922796373466052301324695331E-01;
    weight[ 45-1] =   0.21988529885872983756478409758807E-01;
    weight[ 46-1] =   0.22255787825930280235631416460158E-01;
    weight[ 47-1] =   0.22509534365300608085694429903050E-01;
    weight[ 48-1] =   0.22749615455457959852242553240982E-01;
    weight[ 49-1] =   0.22975885344117206754377437838947E-01;
    weight[ 50-1] =   0.23188206663719640249922582981729E-01;
    weight[ 51-1] =   0.23386450514828194170722043496950E-01;
    weight[ 52-1] =   0.23570496544381716050033676844306E-01;
    weight[ 53-1] =   0.23740233018760777777714726703424E-01;
    weight[ 54-1] =   0.23895556891620665983864481754172E-01;
    weight[ 55-1] =   0.24036373866450369675132086026456E-01;
    weight[ 56-1] =   0.24162598453819584716522917710986E-01;
    weight[ 57-1] =   0.24274154023278979833195063936748E-01;
    weight[ 58-1] =   0.24370972849882214952813561907241E-01;
    weight[ 59-1] =   0.24452996155301467956140198471529E-01;
    weight[ 60-1] =   0.24520174143511508275183033290175E-01;
    weight[ 61-1] =   0.24572466031020653286354137335186E-01;
    weight[ 62-1] =   0.24609840071630254092545634003360E-01;
    weight[ 63-1] =   0.24632273575707679066033370218017E-01;
    weight[ 64-1] =   0.24639752923961094419579417477503E-01;
    weight[ 65-1] =   0.24632273575707679066033370218017E-01;
    weight[ 66-1] =   0.24609840071630254092545634003360E-01;
    weight[ 67-1] =   0.24572466031020653286354137335186E-01;
    weight[ 68-1] =   0.24520174143511508275183033290175E-01;
    weight[ 69-1] =   0.24452996155301467956140198471529E-01;
    weight[ 70-1] =   0.24370972849882214952813561907241E-01;
    weight[ 71-1] =   0.24274154023278979833195063936748E-01;
    weight[ 72-1] =   0.24162598453819584716522917710986E-01;
    weight[ 73-1] =   0.24036373866450369675132086026456E-01;
    weight[ 74-1] =   0.23895556891620665983864481754172E-01;
    weight[ 75-1] =   0.23740233018760777777714726703424E-01;
    weight[ 76-1] =   0.23570496544381716050033676844306E-01;
    weight[ 77-1] =   0.23386450514828194170722043496950E-01;
    weight[ 78-1] =   0.23188206663719640249922582981729E-01;
    weight[ 79-1] =   0.22975885344117206754377437838947E-01;
    weight[ 80-1] =   0.22749615455457959852242553240982E-01;
    weight[ 81-1] =   0.22509534365300608085694429903050E-01;
    weight[ 82-1] =   0.22255787825930280235631416460158E-01;
    weight[ 83-1] =   0.21988529885872983756478409758807E-01;
    weight[ 84-1] =   0.21707922796373466052301324695331E-01;
    weight[ 85-1] =   0.21414136912893259295449693233545E-01;
    weight[ 86-1] =   0.21107350591688713643523847921658E-01;
    weight[ 87-1] =   0.20787750081531811812652137291250E-01;
    weight[ 88-1] =   0.20455529410639508279497065713301E-01;
    weight[ 89-1] =   0.20110890268880247225644623956287E-01;
    weight[ 90-1] =   0.19754041885329183081815217323169E-01;
    weight[ 91-1] =   0.19385200901246454628112623489471E-01;
    weight[ 92-1] =   0.19004591238555646611148901044533E-01;
    weight[ 93-1] =   0.18612443963902310429440419898958E-01;
    weight[ 94-1] =   0.18208997148375106468721469154479E-01;
    weight[ 95-1] =   0.17794495722974774231027912900351E-01;
    weight[ 96-1] =   0.17369191329918731922164721250350E-01;
    weight[ 97-1] =   0.16933342169871654545878815295200E-01;
    weight[ 98-1] =   0.16487212845194879399346060358146E-01;
    weight[ 99-1] =   0.16031074199309941802254151842763E-01;
    weight[100-1] =   0.15565203152273955098532590262975E-01;
    weight[101-1] =   0.15089882532666922992635733981431E-01;
    weight[102-1] =   0.14605400905893418351737288078952E-01;
    weight[103-1] =   0.14112052399003395774044161633613E-01;
    weight[104-1] =   0.13610136522139249906034237533759E-01;
    weight[105-1] =   0.13099957986718627426172681912499E-01;
    weight[106-1] =   0.12581826520465013101514365424172E-01;
    weight[107-1] =   0.12056056679400848183529562144697E-01;
    weight[108-1] =   0.11522967656921087154811609734510E-01;
    weight[109-1] =   0.10982883090068975788799657376065E-01;
    weight[110-1] =   0.10436130863141005225673171997668E-01;
    weight[111-1] =   0.98830429087554914716648010899606E-02;
    weight[112-1] =   0.93239550065309714787536985834029E-02;
    weight[113-1] =   0.87592065795403145773316804234385E-02;
    weight[114-1] =   0.81891404887415730817235884718726E-02;
    weight[115-1] =   0.76141028256526859356393930849227E-02;
    weight[116-1] =   0.70344427036681608755685893032552E-02;
    weight[117-1] =   0.64505120486899171845442463868748E-02;
    weight[118-1] =   0.58626653903523901033648343751367E-02;
    weight[119-1] =   0.52712596565634400891303815906251E-02;
    weight[120-1] =   0.46766539777779034772638165662478E-02;
    weight[121-1] =   0.40792095178254605327114733456293E-02;
    weight[122-1] =   0.34792893810051465908910894094105E-02;
    weight[123-1] =   0.28772587656289004082883197417581E-02;
    weight[124-1] =   0.22734860707492547802810838362671E-02;
    weight[125-1] =   0.16683488125171936761028811985672E-02;
    weight[126-1] =   0.10622766869538486959954760554099E-02;
    weight[127-1] =   0.45645726109586654495731936146574E-03;
  }
  else
  {
    cout << "\n";
    cout << "GL_WEIGHTS - Fatal error!\n";
    cout << "  Illegal value of ORDER = " << order << "\n";
    cout << "  Legal values are 1, 3, 7, 15, 31, 63 and 127.\n";
    exit ( 1 );
  }
  return weight;
}
//****************************************************************************80

double gp_abscissa ( int level, int index )

//****************************************************************************80
//
//  Purpose:
//
//    GP_ABSCISSA returns the I-th abscissa for a Gauss-Patterson rule.
//
//  Discussion:
//
//    The rule is specified by its level.
//
//    The number of points in the rule, known as the order, is
//    related to the level by the formula:
//
//      ORDER = 2^(LEVEL+1)-1.
//
//    Only rules of order 1, 3, 7, 15, 31, 63 and 127 are allowed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Prem Kythe, Michael Schaeferkotter,
//    Handbook of Computational Methods for Integration,
//    Chapman and Hall, 2004,
//    ISBN: 1-58488-428-2,
//    LC: QA299.3.K98.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int LEVEL, the level of the rule.
//    0 <= LEVEL <= 6.
//
//    Input, int INDEX, the index of the point in the rule.
//
//    Output, double GP_ABSCISSA, the value of the INDEX-th 
//    abscissa in the rule of level LEVEL.
//
{
  int order;
  static double x_001[1] = {
     0.0 };
  static double x_003[3] = {
    -0.77459666924148337704, 
     0.0, 
     0.77459666924148337704 };
  static double x_007[7] = {
    -0.96049126870802028342, 
    -0.77459666924148337704, 
    -0.43424374934680255800, 
     0.0, 
     0.43424374934680255800, 
     0.77459666924148337704, 
     0.96049126870802028342 };
  static double x_015[15] = {
    -0.99383196321275502221, 
    -0.96049126870802028342, 
    -0.88845923287225699889, 
    -0.77459666924148337704, 
    -0.62110294673722640294, 
    -0.43424374934680255800, 
    -0.22338668642896688163, 
     0.0, 
     0.22338668642896688163, 
     0.43424374934680255800, 
     0.62110294673722640294, 
     0.77459666924148337704, 
     0.88845923287225699889, 
     0.96049126870802028342, 
     0.99383196321275502221 };
  static double x_031[31] = {
    -0.99909812496766759766, 
    -0.99383196321275502221, 
    -0.98153114955374010687, 
    -0.96049126870802028342, 
    -0.92965485742974005667, 
    -0.88845923287225699889, 
    -0.83672593816886873550, 
    -0.77459666924148337704, 
    -0.70249620649152707861, 
    -0.62110294673722640294, 
    -0.53131974364437562397, 
    -0.43424374934680255800, 
    -0.33113539325797683309, 
    -0.22338668642896688163, 
    -0.11248894313318662575, 
     0.0, 
     0.11248894313318662575, 
     0.22338668642896688163, 
     0.33113539325797683309, 
     0.43424374934680255800, 
     0.53131974364437562397, 
     0.62110294673722640294, 
     0.70249620649152707861, 
     0.77459666924148337704, 
     0.83672593816886873550, 
     0.88845923287225699889, 
     0.92965485742974005667, 
     0.96049126870802028342, 
     0.98153114955374010687, 
     0.99383196321275502221, 
     0.99909812496766759766 };
  static double x_063[63] = {
    -0.99987288812035761194, 
    -0.99909812496766759766, 
    -0.99720625937222195908, 
    -0.99383196321275502221, 
    -0.98868475754742947994, 
    -0.98153114955374010687, 
    -0.97218287474858179658, 
    -0.96049126870802028342, 
    -0.94634285837340290515, 
    -0.92965485742974005667, 
    -0.91037115695700429250, 
    -0.88845923287225699889, 
    -0.86390793819369047715, 
    -0.83672593816886873550, 
    -0.80694053195021761186, 
    -0.77459666924148337704, 
    -0.73975604435269475868, 
    -0.70249620649152707861, 
    -0.66290966002478059546, 
    -0.62110294673722640294, 
    -0.57719571005204581484, 
    -0.53131974364437562397, 
    -0.48361802694584102756, 
    -0.43424374934680255800, 
    -0.38335932419873034692, 
    -0.33113539325797683309, 
    -0.27774982202182431507, 
    -0.22338668642896688163, 
    -0.16823525155220746498, 
    -0.11248894313318662575, 
    -0.056344313046592789972, 
     0.0, 
     0.056344313046592789972, 
     0.11248894313318662575, 
     0.16823525155220746498, 
     0.22338668642896688163, 
     0.27774982202182431507, 
     0.33113539325797683309, 
     0.38335932419873034692, 
     0.43424374934680255800, 
     0.48361802694584102756, 
     0.53131974364437562397, 
     0.57719571005204581484, 
     0.62110294673722640294, 
     0.66290966002478059546, 
     0.70249620649152707861, 
     0.73975604435269475868, 
     0.77459666924148337704, 
     0.80694053195021761186, 
     0.83672593816886873550, 
     0.86390793819369047715, 
     0.88845923287225699889, 
     0.91037115695700429250, 
     0.92965485742974005667, 
     0.94634285837340290515, 
     0.96049126870802028342, 
     0.97218287474858179658, 
     0.98153114955374010687, 
     0.98868475754742947994, 
     0.99383196321275502221, 
     0.99720625937222195908, 
     0.99909812496766759766, 
     0.99987288812035761194 };
  static double x_127[127] = {
    -0.99998243035489159858, 
    -0.99987288812035761194, 
    -0.99959879967191068325, 
    -0.99909812496766759766, 
    -0.99831663531840739253, 
    -0.99720625937222195908, 
    -0.99572410469840718851, 
    -0.99383196321275502221, 
    -0.99149572117810613240, 
    -0.98868475754742947994, 
    -0.98537149959852037111, 
    -0.98153114955374010687, 
    -0.97714151463970571416, 
    -0.97218287474858179658, 
    -0.96663785155841656709, 
    -0.96049126870802028342, 
    -0.95373000642576113641, 
    -0.94634285837340290515, 
    -0.93832039777959288365, 
    -0.92965485742974005667, 
    -0.92034002547001242073, 
    -0.91037115695700429250, 
    -0.89974489977694003664, 
    -0.88845923287225699889, 
    -0.87651341448470526974, 
    -0.86390793819369047715, 
    -0.85064449476835027976, 
    -0.83672593816886873550, 
    -0.82215625436498040737, 
    -0.80694053195021761186, 
    -0.79108493379984836143, 
    -0.77459666924148337704, 
    -0.75748396638051363793, 
    -0.73975604435269475868, 
    -0.72142308537009891548, 
    -0.70249620649152707861, 
    -0.68298743109107922809, 
    -0.66290966002478059546, 
    -0.64227664250975951377, 
    -0.62110294673722640294, 
    -0.59940393024224289297, 
    -0.57719571005204581484, 
    -0.55449513263193254887, 
    -0.53131974364437562397, 
    -0.50768775753371660215, 
    -0.48361802694584102756, 
    -0.45913001198983233287, 
    -0.43424374934680255800, 
    -0.40897982122988867241, 
    -0.38335932419873034692, 
    -0.35740383783153215238, 
    -0.33113539325797683309, 
    -0.30457644155671404334, 
    -0.27774982202182431507, 
    -0.25067873030348317661, 
    -0.22338668642896688163, 
    -0.19589750271110015392, 
    -0.16823525155220746498, 
    -0.14042423315256017459, 
    -0.11248894313318662575, 
    -0.084454040083710883710, 
    -0.056344313046592789972, 
    -0.028184648949745694339, 
     0.0, 
     0.028184648949745694339, 
     0.056344313046592789972, 
     0.084454040083710883710, 
     0.11248894313318662575, 
     0.14042423315256017459, 
     0.16823525155220746498, 
     0.19589750271110015392, 
     0.22338668642896688163, 
     0.25067873030348317661, 
     0.27774982202182431507, 
     0.30457644155671404334, 
     0.33113539325797683309, 
     0.35740383783153215238, 
     0.38335932419873034692, 
     0.40897982122988867241, 
     0.43424374934680255800, 
     0.45913001198983233287, 
     0.48361802694584102756, 
     0.50768775753371660215, 
     0.53131974364437562397, 
     0.55449513263193254887, 
     0.57719571005204581484, 
     0.59940393024224289297, 
     0.62110294673722640294, 
     0.64227664250975951377, 
     0.66290966002478059546, 
     0.68298743109107922809, 
     0.70249620649152707861, 
     0.72142308537009891548, 
     0.73975604435269475868, 
     0.75748396638051363793, 
     0.77459666924148337704, 
     0.79108493379984836143, 
     0.80694053195021761186, 
     0.82215625436498040737, 
     0.83672593816886873550, 
     0.85064449476835027976, 
     0.86390793819369047715, 
     0.87651341448470526974, 
     0.88845923287225699889, 
     0.89974489977694003664, 
     0.91037115695700429250, 
     0.92034002547001242073, 
     0.92965485742974005667, 
     0.93832039777959288365, 
     0.94634285837340290515, 
     0.95373000642576113641, 
     0.96049126870802028342, 
     0.96663785155841656709, 
     0.97218287474858179658, 
     0.97714151463970571416, 
     0.98153114955374010687, 
     0.98537149959852037111, 
     0.98868475754742947994, 
     0.99149572117810613240, 
     0.99383196321275502221, 
     0.99572410469840718851, 
     0.99720625937222195908, 
     0.99831663531840739253, 
     0.99909812496766759766, 
     0.99959879967191068325, 
     0.99987288812035761194, 
     0.99998243035489159858 };
  static double x_255[255] = {
    -0.99999759637974846462, 
    -0.99998243035489159858, 
    -0.99994399620705437576, 
    -0.99987288812035761194, 
    -0.99976049092443204733, 
    -0.99959879967191068325, 
    -0.99938033802502358193, 
    -0.99909812496766759766, 
    -0.99874561446809511470, 
    -0.99831663531840739253, 
    -0.99780535449595727456, 
    -0.99720625937222195908, 
    -0.99651414591489027385, 
    -0.99572410469840718851, 
    -0.99483150280062100052, 
    -0.99383196321275502221, 
    -0.99272134428278861533, 
    -0.99149572117810613240, 
    -0.99015137040077015918, 
    -0.98868475754742947994, 
    -0.98709252795403406719, 
    -0.98537149959852037111, 
    -0.98351865757863272876, 
    -0.98153114955374010687, 
    -0.97940628167086268381, 
    -0.97714151463970571416, 
    -0.97473445975240266776, 
    -0.97218287474858179658, 
    -0.96948465950245923177, 
    -0.96663785155841656709, 
    -0.96364062156981213252, 
    -0.96049126870802028342, 
    -0.95718821610986096274, 
    -0.95373000642576113641, 
    -0.95011529752129487656, 
    -0.94634285837340290515, 
    -0.94241156519108305981, 
    -0.93832039777959288365, 
    -0.93406843615772578800, 
    -0.92965485742974005667, 
    -0.92507893290707565236, 
    -0.92034002547001242073, 
    -0.91543758715576504064, 
    -0.91037115695700429250, 
    -0.90514035881326159519, 
    -0.89974489977694003664, 
    -0.89418456833555902286, 
    -0.88845923287225699889, 
    -0.88256884024734190684, 
    -0.87651341448470526974, 
    -0.87029305554811390585, 
    -0.86390793819369047715, 
    -0.85735831088623215653, 
    -0.85064449476835027976, 
    -0.84376688267270860104, 
    -0.83672593816886873550, 
    -0.82952219463740140018, 
    -0.82215625436498040737, 
    -0.81462878765513741344, 
    -0.80694053195021761186, 
    -0.79909229096084140180, 
    -0.79108493379984836143, 
    -0.78291939411828301639, 
    -0.77459666924148337704, 
    -0.76611781930376009072, 
    -0.75748396638051363793, 
    -0.74869629361693660282, 
    -0.73975604435269475868, 
    -0.73066452124218126133, 
    -0.72142308537009891548, 
    -0.71203315536225203459, 
    -0.70249620649152707861, 
    -0.69281376977911470289, 
    -0.68298743109107922809, 
    -0.67301883023041847920, 
    -0.66290966002478059546, 
    -0.65266166541001749610, 
    -0.64227664250975951377, 
    -0.63175643771119423041, 
    -0.62110294673722640294, 
    -0.61031811371518640016, 
    -0.59940393024224289297, 
    -0.58836243444766254143, 
    -0.57719571005204581484, 
    -0.56590588542365442262, 
    -0.55449513263193254887, 
    -0.54296566649831149049, 
    -0.53131974364437562397, 
    -0.51955966153745702199, 
    -0.50768775753371660215, 
    -0.49570640791876146017, 
    -0.48361802694584102756, 
    -0.47142506587165887693, 
    -0.45913001198983233287, 
    -0.44673538766202847374, 
    -0.43424374934680255800, 
    -0.42165768662616330006, 
    -0.40897982122988867241, 
    -0.39621280605761593918, 
    -0.38335932419873034692, 
    -0.37042208795007823014, 
    -0.35740383783153215238, 
    -0.34430734159943802278, 
    -0.33113539325797683309, 
    -0.31789081206847668318, 
    -0.30457644155671404334, 
    -0.29119514851824668196, 
    -0.27774982202182431507, 
    -0.26424337241092676194, 
    -0.25067873030348317661, 
    -0.23705884558982972721, 
    -0.22338668642896688163, 
    -0.20966523824318119477, 
    -0.19589750271110015392, 
    -0.18208649675925219825, 
    -0.16823525155220746498, 
    -0.15434681148137810869, 
    -0.14042423315256017459, 
    -0.12647058437230196685, 
    -0.11248894313318662575, 
    -0.098482396598119202090, 
    -0.084454040083710883710, 
    -0.070406976042855179063, 
    -0.056344313046592789972, 
    -0.042269164765363603212, 
    -0.028184648949745694339, 
    -0.014093886410782462614, 
    0.0, 
    0.014093886410782462614, 
    0.028184648949745694339, 
    0.042269164765363603212, 
    0.056344313046592789972, 
    0.070406976042855179063, 
    0.084454040083710883710, 
    0.098482396598119202090, 
    0.11248894313318662575, 
    0.12647058437230196685, 
    0.14042423315256017459, 
    0.15434681148137810869, 
    0.16823525155220746498, 
    0.18208649675925219825, 
    0.19589750271110015392, 
    0.20966523824318119477, 
    0.22338668642896688163, 
    0.23705884558982972721, 
    0.25067873030348317661, 
    0.26424337241092676194, 
    0.27774982202182431507, 
    0.29119514851824668196, 
    0.30457644155671404334, 
    0.31789081206847668318, 
    0.33113539325797683309, 
    0.34430734159943802278, 
    0.35740383783153215238, 
    0.37042208795007823014, 
    0.38335932419873034692, 
    0.39621280605761593918, 
    0.40897982122988867241, 
    0.42165768662616330006, 
    0.43424374934680255800, 
    0.44673538766202847374, 
    0.45913001198983233287, 
    0.47142506587165887693, 
    0.48361802694584102756, 
    0.49570640791876146017, 
    0.50768775753371660215, 
    0.51955966153745702199, 
    0.53131974364437562397, 
    0.54296566649831149049, 
    0.55449513263193254887, 
    0.56590588542365442262, 
    0.57719571005204581484, 
    0.58836243444766254143, 
    0.59940393024224289297, 
    0.61031811371518640016, 
    0.62110294673722640294, 
    0.63175643771119423041, 
    0.64227664250975951377, 
    0.65266166541001749610, 
    0.66290966002478059546, 
    0.67301883023041847920, 
    0.68298743109107922809, 
    0.69281376977911470289, 
    0.70249620649152707861, 
    0.71203315536225203459, 
    0.72142308537009891548, 
    0.73066452124218126133, 
    0.73975604435269475868, 
    0.74869629361693660282, 
    0.75748396638051363793, 
    0.76611781930376009072, 
    0.77459666924148337704, 
    0.78291939411828301639, 
    0.79108493379984836143, 
    0.79909229096084140180, 
    0.80694053195021761186, 
    0.81462878765513741344, 
    0.82215625436498040737, 
    0.82952219463740140018, 
    0.83672593816886873550, 
    0.84376688267270860104, 
    0.85064449476835027976, 
    0.85735831088623215653, 
    0.86390793819369047715, 
    0.87029305554811390585, 
    0.87651341448470526974, 
    0.88256884024734190684, 
    0.88845923287225699889, 
    0.89418456833555902286, 
    0.89974489977694003664, 
    0.90514035881326159519, 
    0.91037115695700429250, 
    0.91543758715576504064, 
    0.92034002547001242073, 
    0.92507893290707565236, 
    0.92965485742974005667, 
    0.93406843615772578800, 
    0.93832039777959288365, 
    0.94241156519108305981, 
    0.94634285837340290515, 
    0.95011529752129487656, 
    0.95373000642576113641, 
    0.95718821610986096274, 
    0.96049126870802028342, 
    0.96364062156981213252, 
    0.96663785155841656709, 
    0.96948465950245923177, 
    0.97218287474858179658, 
    0.97473445975240266776, 
    0.97714151463970571416, 
    0.97940628167086268381, 
    0.98153114955374010687, 
    0.98351865757863272876, 
    0.98537149959852037111, 
    0.98709252795403406719, 
    0.98868475754742947994, 
    0.99015137040077015918, 
    0.99149572117810613240, 
    0.99272134428278861533, 
    0.99383196321275502221, 
    0.99483150280062100052, 
    0.99572410469840718851, 
    0.99651414591489027385, 
    0.99720625937222195908, 
    0.99780535449595727456, 
    0.99831663531840739253, 
    0.99874561446809511470, 
    0.99909812496766759766, 
    0.99938033802502358193, 
    0.99959879967191068325, 
    0.99976049092443204733, 
    0.99987288812035761194, 
    0.99994399620705437576, 
    0.99998243035489159858, 
    0.99999759637974846462 };
  double value;

  order = i4_power ( 2, level + 1 ) - 1;

  if ( order < 1 )
  {
    value = - r8_huge ( );
  }
  else if ( index < 1 || order < index )
  {
    value = - r8_huge ( );
  }
  else if ( order == 1 )
  {
    value = x_001[index-1];
  }
  else if ( order == 3 )
  {
    value = x_003[index-1];
  }
  else if ( order == 7 )
  {
    value = x_007[index-1];
  }
  else if ( order == 15 )
  {
    value = x_015[index-1];
  }
  else if ( order == 31 )
  {
    value = x_031[index-1];
  }
  else if ( order == 63 )
  {
    value = x_063[index-1];
  }
  else if ( order == 127 )
  {
    value = x_127[index-1];
  }
  else if ( order == 255 )
  {
    value = x_255[index-1];
  }
  else
  {
    value = - r8_huge ( );
  }
  return value;
}
//****************************************************************************80

double *gp_weights ( int order )

//****************************************************************************80
//
//  Purpose:
//
//    GP_WEIGHTS sets weights for a Gauss-Patterson rule.
//
//  Discussion:
//
//    The zeroth rule, of order 1, is the standard Gauss-Legendre rule.
//
//    The first rule, of order 3, is the standard Gauss-Legendre rule.
//
//    The second rule, of order 7, includes the abscissas of the previous
//    rule.
//
//    Each subsequent rule is nested in a similar way.  Rules are available
//    of orders 1, 3, 7, 15, 31, 63, 127 and 255
//
//  Modified:
//
//    23 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Prem Kythe, Michael Schaeferkotter,
//    Handbook of Computational Methods for Integration,
//    Chapman and Hall, 2004,
//    ISBN: 1-58488-428-2,
//    LC: QA299.3.K98.
//
//    Thomas Patterson,
//    The Optimal Addition of Points to Quadrature Formulae,
//    Mathematics of Computation,
//    Volume 22, Number 104, October 1968, pages 847-856.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    ORDER must be 1, 3, 7, 15, 31, 63, 127 or 255.
//
//    Output, double W[ORDER], the weights of the rule.
//    The weights are positive, symmetric and should sum to 2.
//
{
  double *w;
  static double w_001[1] =
  {
    2.0
  };
  static double w_003[3] = 
  {
    0.555555555555555555556,
    0.888888888888888888889,
    0.555555555555555555556
  };
  static double w_007[7] =
  {
    0.104656226026467265194,
    0.268488089868333440729,
    0.401397414775962222905,
    0.450916538658474142345,
    0.401397414775962222905,
    0.268488089868333440729,
    0.104656226026467265194
  };
  static double w_015[15] =
  {
    0.0170017196299402603390,
    0.0516032829970797396969,
    0.0929271953151245376859,
    0.134415255243784220360,
    0.171511909136391380787,
    0.200628529376989021034,
    0.219156858401587496404,
    0.225510499798206687386,
    0.219156858401587496404,
    0.200628529376989021034,
    0.171511909136391380787,
    0.134415255243784220360,
    0.0929271953151245376859,
    0.0516032829970797396969,
    0.0170017196299402603390
  };
  static double w_031[31] =
  {
    0.00254478079156187441540,
    0.00843456573932110624631,
    0.0164460498543878109338,
    0.0258075980961766535646,
    0.0359571033071293220968,
    0.0464628932617579865414,
    0.0569795094941233574122,
    0.0672077542959907035404,
    0.0768796204990035310427,
    0.0857559200499903511542,
    0.0936271099812644736167,
    0.100314278611795578771,
    0.105669893580234809744,
    0.109578421055924638237,
    0.111956873020953456880,
    0.112755256720768691607,
    0.111956873020953456880,
    0.109578421055924638237,
    0.105669893580234809744,
    0.100314278611795578771,
    0.0936271099812644736167,
    0.0857559200499903511542,
    0.0768796204990035310427,
    0.0672077542959907035404,
    0.0569795094941233574122,
    0.0464628932617579865414,
    0.0359571033071293220968,
    0.0258075980961766535646,
    0.0164460498543878109338,
    0.00843456573932110624631,
    0.00254478079156187441540
  };
  static double w_063[63] =
  {
    0.000363221481845530659694,
    0.00126515655623006801137,
    0.00257904979468568827243,
    0.00421763044155885483908,
    0.00611550682211724633968,
    0.00822300795723592966926,
    0.0104982469096213218983,
    0.0129038001003512656260,
    0.0154067504665594978021,
    0.0179785515681282703329,
    0.0205942339159127111492,
    0.0232314466399102694433,
    0.0258696793272147469108,
    0.0284897547458335486125,
    0.0310735511116879648799,
    0.0336038771482077305417,
    0.0360644327807825726401,
    0.0384398102494555320386,
    0.0407155101169443189339,
    0.0428779600250077344929,
    0.0449145316536321974143,
    0.0468135549906280124026,
    0.0485643304066731987159,
    0.0501571393058995374137,
    0.0515832539520484587768,
    0.0528349467901165198621,
    0.0539054993352660639269,
    0.0547892105279628650322,
    0.0554814043565593639878,
    0.0559784365104763194076,
    0.0562776998312543012726,
    0.0563776283603847173877,
    0.0562776998312543012726,
    0.0559784365104763194076,
    0.0554814043565593639878,
    0.0547892105279628650322,
    0.0539054993352660639269,
    0.0528349467901165198621,
    0.0515832539520484587768,
    0.0501571393058995374137,
    0.0485643304066731987159,
    0.0468135549906280124026,
    0.0449145316536321974143,
    0.0428779600250077344929,
    0.0407155101169443189339,
    0.0384398102494555320386,
    0.0360644327807825726401,
    0.0336038771482077305417,
    0.0310735511116879648799,
    0.0284897547458335486125,
    0.0258696793272147469108,
    0.0232314466399102694433,
    0.0205942339159127111492,
    0.0179785515681282703329,
    0.0154067504665594978021,
    0.0129038001003512656260,
    0.0104982469096213218983,
    0.00822300795723592966926,
    0.00611550682211724633968,
    0.00421763044155885483908,
    0.00257904979468568827243,
    0.00126515655623006801137,
    0.000363221481845530659694
  };
  static double w_127[127] =
  {
    0.0000505360952078625176247,
    0.000180739564445388357820,
    0.000377746646326984660274,
    0.000632607319362633544219,
    0.000938369848542381500794,
    0.00128952408261041739210,
    0.00168114286542146990631,
    0.00210881524572663287933,
    0.00256876494379402037313,
    0.00305775341017553113613,
    0.00357289278351729964938,
    0.00411150397865469304717,
    0.00467105037211432174741,
    0.00524912345480885912513,
    0.00584344987583563950756,
    0.00645190005017573692280,
    0.00707248999543355546805,
    0.00770337523327974184817,
    0.00834283875396815770558,
    0.00898927578406413572328,
    0.00964117772970253669530,
    0.0102971169579563555237,
    0.0109557333878379016480,
    0.0116157233199551347270,
    0.0122758305600827700870,
    0.0129348396636073734547,
    0.0135915710097655467896,
    0.0142448773729167743063,
    0.0148936416648151820348,
    0.0155367755558439824399,
    0.0161732187295777199419,
    0.0168019385741038652709,
    0.0174219301594641737472,
    0.0180322163903912863201,
    0.0186318482561387901863,
    0.0192199051247277660193,
    0.0197954950480974994880,
    0.0203577550584721594669,
    0.0209058514458120238522,
    0.0214389800125038672465,
    0.0219563663053178249393,
    0.0224572658268160987071,
    0.0229409642293877487608,
    0.0234067774953140062013,
    0.0238540521060385400804,
    0.0242821652033365993580,
    0.0246905247444876769091,
    0.0250785696529497687068,
    0.0254457699654647658126,
    0.0257916269760242293884,
    0.0261156733767060976805,
    0.0264174733950582599310,
    0.0266966229274503599062,
    0.0269527496676330319634,
    0.0271855132296247918192,
    0.0273946052639814325161,
    0.0275797495664818730349,
    0.0277407021782796819939,
    0.0278772514766137016085,
    0.0279892182552381597038,
    0.0280764557938172466068,
    0.0281388499156271506363,
    0.0281763190330166021307,
    0.0281888141801923586938,
    0.0281763190330166021307,
    0.0281388499156271506363,
    0.0280764557938172466068,
    0.0279892182552381597038,
    0.0278772514766137016085,
    0.0277407021782796819939,
    0.0275797495664818730349,
    0.0273946052639814325161,
    0.0271855132296247918192,
    0.0269527496676330319634,
    0.0266966229274503599062,
    0.0264174733950582599310,
    0.0261156733767060976805,
    0.0257916269760242293884,
    0.0254457699654647658126,
    0.0250785696529497687068,
    0.0246905247444876769091,
    0.0242821652033365993580,
    0.0238540521060385400804,
    0.0234067774953140062013,
    0.0229409642293877487608,
    0.0224572658268160987071,
    0.0219563663053178249393,
    0.0214389800125038672465,
    0.0209058514458120238522,
    0.0203577550584721594669,
    0.0197954950480974994880,
    0.0192199051247277660193,
    0.0186318482561387901863,
    0.0180322163903912863201,
    0.0174219301594641737472,
    0.0168019385741038652709,
    0.0161732187295777199419,
    0.0155367755558439824399,
    0.0148936416648151820348,
    0.0142448773729167743063,
    0.0135915710097655467896,
    0.0129348396636073734547,
    0.0122758305600827700870,
    0.0116157233199551347270,
    0.0109557333878379016480,
    0.0102971169579563555237,
    0.00964117772970253669530,
    0.00898927578406413572328,
    0.00834283875396815770558,
    0.00770337523327974184817,
    0.00707248999543355546805,
    0.00645190005017573692280,
    0.00584344987583563950756,
    0.00524912345480885912513,
    0.00467105037211432174741,
    0.00411150397865469304717,
    0.00357289278351729964938,
    0.00305775341017553113613,
    0.00256876494379402037313,
    0.00210881524572663287933,
    0.00168114286542146990631,
    0.00128952408261041739210,
    0.000938369848542381500794,
    0.000632607319362633544219,
    0.000377746646326984660274,
    0.000180739564445388357820,
    0.0000505360952078625176247
  };
  static double w_255[255] =
  {
    0.69379364324108267170E-05,
    0.25157870384280661489E-04,
    0.53275293669780613125E-04,
    0.90372734658751149261E-04,
    0.13575491094922871973E-03,
    0.18887326450650491366E-03,
    0.24921240048299729402E-03,
    0.31630366082226447689E-03,
    0.38974528447328229322E-03,
    0.46918492424785040975E-03,
    0.55429531493037471492E-03,
    0.64476204130572477933E-03,
    0.74028280424450333046E-03,
    0.84057143271072246365E-03,
    0.94536151685852538246E-03,
    0.10544076228633167722E-02,
    0.11674841174299594077E-02,
    0.12843824718970101768E-02,
    0.14049079956551446427E-02,
    0.15288767050877655684E-02,
    0.16561127281544526052E-02,
    0.17864463917586498247E-02,
    0.19197129710138724125E-02,
    0.20557519893273465236E-02,
    0.21944069253638388388E-02,
    0.23355251860571608737E-02,
    0.24789582266575679307E-02,
    0.26245617274044295626E-02,
    0.27721957645934509940E-02,
    0.29217249379178197538E-02,
    0.30730184347025783234E-02,
    0.32259500250878684614E-02,
    0.33803979910869203823E-02,
    0.35362449977167777340E-02,
    0.36933779170256508183E-02,
    0.38516876166398709241E-02,
    0.40110687240750233989E-02,
    0.41714193769840788528E-02,
    0.43326409680929828545E-02,
    0.44946378920320678616E-02,
    0.46573172997568547773E-02,
    0.48205888648512683476E-02,
    0.49843645647655386012E-02,
    0.51485584789781777618E-02,
    0.53130866051870565663E-02,
    0.54778666939189508240E-02,
    0.56428181013844441585E-02,
    0.58078616599775673635E-02,
    0.59729195655081658049E-02,
    0.61379152800413850435E-02,
    0.63027734490857587172E-02,
    0.64674198318036867274E-02,
    0.66317812429018878941E-02,
    0.67957855048827733948E-02,
    0.69593614093904229394E-02,
    0.71224386864583871532E-02,
    0.72849479805538070639E-02,
    0.74468208324075910174E-02,
    0.76079896657190565832E-02,
    0.77683877779219912200E-02,
    0.79279493342948491103E-02,
    0.80866093647888599710E-02,
    0.82443037630328680306E-02,
    0.84009692870519326354E-02,
    0.85565435613076896192E-02,
    0.87109650797320868736E-02,
    0.88641732094824942641E-02,
    0.90161081951956431600E-02,
    0.91667111635607884067E-02,
    0.93159241280693950932E-02,
    0.94636899938300652943E-02,
    0.96099525623638830097E-02,
    0.97546565363174114611E-02,
    0.98977475240487497440E-02,
    0.10039172044056840798E-01,
    0.10178877529236079733E-01,
    0.10316812330947621682E-01,
    0.10452925722906011926E-01,
    0.10587167904885197931E-01,
    0.10719490006251933623E-01,
    0.10849844089337314099E-01,
    0.10978183152658912470E-01,
    0.11104461134006926537E-01,
    0.11228632913408049354E-01,
    0.11350654315980596602E-01,
    0.11470482114693874380E-01,
    0.11588074033043952568E-01,
    0.11703388747657003101E-01,
    0.11816385890830235763E-01,
    0.11927026053019270040E-01,
    0.12035270785279562630E-01,
    0.12141082601668299679E-01,
    0.12244424981611985899E-01,
    0.12345262372243838455E-01,
    0.12443560190714035263E-01,
    0.12539284826474884353E-01,
    0.12632403643542078765E-01,
    0.12722884982732382906E-01,
    0.12810698163877361967E-01,
    0.12895813488012114694E-01,
    0.12978202239537399286E-01,
    0.13057836688353048840E-01,
    0.13134690091960152836E-01,
    0.13208736697529129966E-01,
    0.13279951743930530650E-01,
    0.13348311463725179953E-01,
    0.13413793085110098513E-01,
    0.13476374833816515982E-01,
    0.13536035934956213614E-01,
    0.13592756614812395910E-01,
    0.13646518102571291428E-01,
    0.13697302631990716258E-01,
    0.13745093443001896632E-01,
    0.13789874783240936517E-01,
    0.13831631909506428676E-01,
    0.13870351089139840997E-01,
    0.13906019601325461264E-01,
    0.13938625738306850804E-01,
    0.13968158806516938516E-01,
    0.13994609127619079852E-01,
    0.14017968039456608810E-01,
    0.14038227896908623303E-01,
    0.14055382072649964277E-01,
    0.14069424957813575318E-01,
    0.14080351962553661325E-01,
    0.14088159516508301065E-01,
    0.14092845069160408355E-01,
    0.14094407090096179347E-01,
    0.14092845069160408355E-01,
    0.14088159516508301065E-01,
    0.14080351962553661325E-01,
    0.14069424957813575318E-01,
    0.14055382072649964277E-01,
    0.14038227896908623303E-01,
    0.14017968039456608810E-01,
    0.13994609127619079852E-01,
    0.13968158806516938516E-01,
    0.13938625738306850804E-01,
    0.13906019601325461264E-01,
    0.13870351089139840997E-01,
    0.13831631909506428676E-01,
    0.13789874783240936517E-01,
    0.13745093443001896632E-01,
    0.13697302631990716258E-01,
    0.13646518102571291428E-01,
    0.13592756614812395910E-01,
    0.13536035934956213614E-01,
    0.13476374833816515982E-01,
    0.13413793085110098513E-01,
    0.13348311463725179953E-01,
    0.13279951743930530650E-01,
    0.13208736697529129966E-01,
    0.13134690091960152836E-01,
    0.13057836688353048840E-01,
    0.12978202239537399286E-01,
    0.12895813488012114694E-01,
    0.12810698163877361967E-01,
    0.12722884982732382906E-01,
    0.12632403643542078765E-01,
    0.12539284826474884353E-01,
    0.12443560190714035263E-01,
    0.12345262372243838455E-01,
    0.12244424981611985899E-01,
    0.12141082601668299679E-01,
    0.12035270785279562630E-01,
    0.11927026053019270040E-01,
    0.11816385890830235763E-01,
    0.11703388747657003101E-01,
    0.11588074033043952568E-01,
    0.11470482114693874380E-01,
    0.11350654315980596602E-01,
    0.11228632913408049354E-01,
    0.11104461134006926537E-01,
    0.10978183152658912470E-01,
    0.10849844089337314099E-01,
    0.10719490006251933623E-01,
    0.10587167904885197931E-01,
    0.10452925722906011926E-01,
    0.10316812330947621682E-01,
    0.10178877529236079733E-01,
    0.10039172044056840798E-01,
    0.98977475240487497440E-02,
    0.97546565363174114611E-02,
    0.96099525623638830097E-02,
    0.94636899938300652943E-02,
    0.93159241280693950932E-02,
    0.91667111635607884067E-02,
    0.90161081951956431600E-02,
    0.88641732094824942641E-02,
    0.87109650797320868736E-02,
    0.85565435613076896192E-02,
    0.84009692870519326354E-02,
    0.82443037630328680306E-02,
    0.80866093647888599710E-02,
    0.79279493342948491103E-02,
    0.77683877779219912200E-02,
    0.76079896657190565832E-02,
    0.74468208324075910174E-02,
    0.72849479805538070639E-02,
    0.71224386864583871532E-02,
    0.69593614093904229394E-02,
    0.67957855048827733948E-02,
    0.66317812429018878941E-02,
    0.64674198318036867274E-02,
    0.63027734490857587172E-02,
    0.61379152800413850435E-02,
    0.59729195655081658049E-02,
    0.58078616599775673635E-02,
    0.56428181013844441585E-02,
    0.54778666939189508240E-02,
    0.53130866051870565663E-02,
    0.51485584789781777618E-02,
    0.49843645647655386012E-02,
    0.48205888648512683476E-02,
    0.46573172997568547773E-02,
    0.44946378920320678616E-02,
    0.43326409680929828545E-02,
    0.41714193769840788528E-02,
    0.40110687240750233989E-02,
    0.38516876166398709241E-02,
    0.36933779170256508183E-02,
    0.35362449977167777340E-02,
    0.33803979910869203823E-02,
    0.32259500250878684614E-02,
    0.30730184347025783234E-02,
    0.29217249379178197538E-02,
    0.27721957645934509940E-02,
    0.26245617274044295626E-02,
    0.24789582266575679307E-02,
    0.23355251860571608737E-02,
    0.21944069253638388388E-02,
    0.20557519893273465236E-02,
    0.19197129710138724125E-02,
    0.17864463917586498247E-02,
    0.16561127281544526052E-02,
    0.15288767050877655684E-02,
    0.14049079956551446427E-02,
    0.12843824718970101768E-02,
    0.11674841174299594077E-02,
    0.10544076228633167722E-02,
    0.94536151685852538246E-03,
    0.84057143271072246365E-03,
    0.74028280424450333046E-03,
    0.64476204130572477933E-03,
    0.55429531493037471492E-03,
    0.46918492424785040975E-03,
    0.38974528447328229322E-03,
    0.31630366082226447689E-03,
    0.24921240048299729402E-03,
    0.18887326450650491366E-03,
    0.13575491094922871973E-03,
    0.90372734658751149261E-04,
    0.53275293669780613125E-04,
    0.25157870384280661489E-04,
    0.69379364324108267170E-05
  };

  w = new double[order];

  if ( order == 1 )
  {
    r8vec_copy ( order, w_001, w );
  }
  else if ( order == 3 )
  {
    r8vec_copy ( order, w_003, w );
  }
  else if ( order == 7 )
  {
    r8vec_copy ( order, w_007, w );
  }
  else if ( order == 15 )
  {
    r8vec_copy ( order, w_015, w );
  }
  else if ( order == 31 )
  {
    r8vec_copy ( order, w_031, w );
  }
  else if ( order == 63 )
  {
    r8vec_copy ( order, w_063, w );
  }
  else if ( order == 127 )
  {
    r8vec_copy ( order, w_127, w );
  }
  else if ( order == 255 )
  {
    r8vec_copy ( order, w_255, w );
  }
  else
  {
    std::cerr << "\n";
    std::cerr << "GP_WEIGHTS - Fatal error!\n";
    std::cerr << "  Unexpected value of ORDER = " << order << ".\n";
    std::exit ( 1 );
  }

  return w;
}
//****************************************************************************80

int i4_log_2 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
//
//  Example:
//
//        I  I4_LOG_10
//    -----  --------
//        0    0
//        1    0
//        2    1
//        3    1
//        4    2
//        5    2
//        7    2
//        8    3
//        9    3
//     1000    9
//     1024   10
//
//  Discussion:
//
//    I4_LOG_2 ( I ) + 1 is the number of binary digits in I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number whose logarithm base 2 is desired.
//
//    Output, int I4_LOG_2, the integer part of the logarithm base 2 of
//    the absolute value of X.
//
{
  int i_abs;
  int two_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    two_pow = 2;

    i_abs = abs ( i );

    while ( two_pow <= i_abs )
    {
      value = value + 1;
      two_pow = two_pow * 2;
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

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If 
//      NREM = I4_MODP ( I, J ) 
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
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
//    26 May 1999
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
//    01 April 2004
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
//    09 November 2007
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

int *index_level_own ( int level, int level_max, int dim_num, int point_num, 
  int grid_index[], int grid_base[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_LEVEL_OWN: determine first level at which given index is generated.
//
//  Discussion:
//
//    We are constructing a sparse grid of OWN points.  The grid
//    is built up of product grids, with a characteristic LEVEL.  
//
//    We are concerned with identifying points in this product grid which
//    have actually been generated previously, on a lower value of LEVEL.
//
//    This routine determines the lowest value of LEVEL at which each of
//    the input points would be generated.
//
//    In 1D, given LEVEL, the number of points is ORDER = 2**(LEVEL+1) + 1,
//    (except that LEVEL = 0 implies ORDER = 1), the BASE is (ORDER-1)/2, 
//    and the point INDEX values range from -BASE to +BASE.
//
//    The values of INDEX and BASE allow us to determine the abstract
//    properties of the point.  In particular, if INDEX is 0, the corresponding
//    abscissa is 0, the special "nested" value we need to take care of.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2007
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
//    Input, int LEVEL, the level at which these points were 
//    generated.  LEVEL_MIN <= LEVEL <= LEVEL_MAX.
//
//    Input, int LEVEL_MAX, the maximum level.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points to be tested.
//
//    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the indices of the 
//    points to be tested.
//
//    Input, int GRID_BASE[DIM_NUM], the "base", which is essentially
//    the denominator of the index.
//
//    Output, int INDEX_LEVEL_OWN[POINT_NUM], the value of LEVEL at 
//    which the point would first be generated.  This will be the same as
//    the input value of LEVEL, unless the point has an INDEX of 0 and
//    a corresponding BASE that is NOT zero.
//
{
  int dim;
  int *grid_level;
  int level_min;
  int point;
  
  grid_level = new int[point_num];
  
  if ( dim_num == 1 )
  {
    level_min = level_max;
  }
  else
  {
    level_min = 0;
  }
//
//  If a point has a DIM-th component whose INDEX is 0, then the 
//  value of LEVEL at which this point would first be generated is
//  less than LEVEL, unless the DIM-th component of GRID_BASE is 0.
//
  for ( point = 0; point < point_num; point++ )
  {
    grid_level[point] = i4_max ( level, level_min );

    for ( dim = 0; dim < dim_num; dim++ )
    {
      if ( grid_index[dim+point*dim_num] == 0 )
      {
        grid_level[point] = i4_max ( grid_level[point] - grid_base[dim], level_min );
      }
    }
  }

  return grid_level;
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

int index_to_level_open ( int dim_num, int t[], int order, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_TO_LEVEL_OPEN determines the level of a point given its index.
//
//  Modified:
//
//    19 April 2007
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
//    Input, int T[DIM_NUM], the grid index of a point.
//
//    Input, int ORDER, the order of the rule.
//
//    Input, int LEVEL_MAX, the level with respect to which the
//    index applies.
//
//    Output, int INDEX_TO_LEVEL_OPEN, the first level on which
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

void level_to_order_open ( int dim_num, int level[], int order[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEVEL_TO_ORDER_OPEN converts a level to an order for open rules.
//
//  Discussion:
//
//    Sparse grids can naturally be nested.  A natural scheme is to use
//    a series of one-dimensional rules arranged in a series of "levels"
//    whose order roughly doubles with each step.
//
//    The arrangement described here works naturally for the Fejer Type 1,
//    Fejer Type 2, Newton Cotes Open, Newton Cotes Half Open,
//    and Gauss-Patterson rules.  It also can be used, partially, to describe
//    the growth of Gauss-Legendre rules.
//
//    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single 
//    point at the center, and for all values afterwards, we use the relationship
//
//      ORDER = 2**(LEVEL+1) - 1.
//
//    The following table shows how the growth will occur:
//
//    Level    Order
//
//    0          1
//    1          3 =  4 - 1
//    2          7 =  8 - 1
//    3         15 = 16 - 1
//    4         31 = 32 - 1
//    5         63 = 64 - 1
//
//    For the Fejer Type 1, Fejer Type 2, Newton Cotes Open, 
//    Newton Cotes Open Half, and Gauss-Patterson rules, the point growth is
//    nested.  If we have ORDER points on a particular LEVEL, the next level 
//    includes all these old points, plus ORDER+1 new points, formed in the 
//    gaps between successive pairs of old points plus an extra point at each 
//    end.
//
//    Level    Order = New + Old
//
//    0          1   =  1  +  0
//    1          3   =  2  +  1
//    2          7   =  4  +  3
//    3         15   =  8  +  7
//    4         31   = 16  + 15
//    5         63   = 32  + 31
//
//    If we use a series of Gauss-Legendre rules, then there is almost no 
//    nesting, except that the central point is shared.  If we insist on 
//    producing a comparable series of such points, then the "nesting" behavior
//    is as follows:
//
//    Level    Order = New + Old
//
//    0          1   =  1  +  0
//    1          3   =  2  +  1
//    2          7   =  6  +  1
//    3         15   = 14  +  1
//    4         31   = 30  +  1
//    5         63   = 62  +  1
//
//    Moreover, if we consider ALL the points used in such a set of "nested" 
//    Gauss-Legendre rules, then we must sum the "NEW" column, and we see that
//    we get roughly twice as many points as for the truly nested rules.
//
//    In this routine, we assume that a vector of levels is given,
//    and the corresponding orders are desired.
//
//  Modified:
//
//    19 April 2007
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
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1 ;
    }
  }
  return;
}
//****************************************************************************80

void levels_index ( int dim_num, int level_max, int rule, int point_num, 
  int grid_index[], int grid_base[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX indexes a sparse grid.
//
//  Discussion:
//
//    The sparse grid is the logical sum of product grids with total LEVEL 
//    between LEVEL_MIN and LEVEL_MAX.
//
//    The necessary dimensions of GRID_INDEX can be determined by 
//    calling LEVELS_INDEX_SIZE first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 March 2008
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
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Input, int POINT_NUM, the total number of points 
//    in the grids.
//
//    Output, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of 
//    point indices, representing a subset of the product grid of level 
//    LEVEL_MAX, representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
//    Output, int GRID_BASE[DIM_NUM*POINT_NUM], a list of 
//    the orders of the rules associated with each point and dimension.
//
{
  if ( rule == 1 )
  {
   levels_index_cfn ( dim_num, level_max, point_num, grid_index, grid_base );
  }
  else if ( 2 <= rule && rule <= 4 )
  {
   levels_index_ofn ( dim_num, level_max, point_num, grid_index, grid_base );
  }
  else if ( 5 <= rule && rule <= 6 )
  {
    levels_index_own ( dim_num, level_max, point_num, grid_index, grid_base );
  }
  else if ( 7 == rule )
  {
    levels_index_onn ( dim_num, level_max, point_num, grid_index, grid_base );
  }
  else
  {
    cout << "\n";
    cout << "LEVELS_INDEX - Fatal error!\n";
    cout << "  Unrecognized rule number = " << rule << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void levels_index_cfn ( int dim_num, int level_max, int point_num, 
  int grid_index[], int grid_base[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_CFN indexes a sparse grid made from CFN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of CLOSED FULLY NESTED 1D quadrature rules.
//
//    CFN rules include Clenshaw Curtis rules.
//
//    The sparse grid is the logical sum of product grids with total LEVEL 
//    between LEVEL_MIN and LEVEL_MAX.
//
//    The necessary dimensions of GRID_INDEX can be determined by 
//    calling LEVELS_INDEX_SIZE_CFN first.
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
//    Output, int LEVELS_INDEX_CFN[DIM_NUM*POINT_NUM], a list of point 
//    indices, representing a subset of the product grid of level LEVEL_MAX,
//    representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
{
  int dim;
  int factor;
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
      grid_index2 = multigrid_index_cfn ( dim_num, order_1d, order_nd );
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
          if ( point_num <= point_num2 )
          {
            cout << "\n";
            cout << "LEVELS_INDEX_CFN - Fatal error!\n";
            cout << "  Exceeding maximum point index POINT_NUM = "
                 << point_num << "\n";
            exit ( 1 );
          }

          for ( dim = 0; dim < dim_num; dim++ )
          {
            grid_base[dim+point_num2*dim_num] = order_1d[dim];
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

  if ( point_num2 < point_num )
  {
    cout << "\n";
    cout << "LEVELS_INDEX_CFN - Fatal error!\n";
    cout << "  Set fewer points than POINT_NUM = " << point_num << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void levels_index_ofn ( int dim_num, int level_max, int point_num, 
  int grid_index[], int grid_base[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_OFN indexes a sparse grid made from OFN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of OPEN FULLY NESTED 1D quadrature rules.
//
//    OFN rules include Fejer 1, Fejer 2, and Gauss Patterson rules.
//
//    The sparse grid is the logical sum of product grids with total LEVEL 
//    between LEVEL_MIN and LEVEL_MAX.
//
//    The necessary dimensions of GRID_INDEX can be determined by 
//    calling LEVELS_INDEX_SIZE_OFN first.
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
//    Output, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of
//    point indices, representing a subset of the product grid of level 
//    LEVEL_MAX, representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
//    Output, int GRID_BASE[DIM_NUM*POINT_NUM], a list of 
//    the orders of the rules associated with each point and dimension.
//
{
  int dim;
  int factor;
  int *grid_index2;
  int h;
  int level;
  int *level_1d;
  bool more;
  int *order_1d;
  int order_nd;
  int point;
  int point_num2;
  int t;
  bool test;
//
//  The outer loop generates LEVELs from 0 to LEVEL_MAX.
//
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  point_num2 = 0;

  for ( level = 0; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition that adds up to LEVEL.
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
      level_to_order_open ( dim_num, level_1d, order_1d );
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//
      grid_index2 = multigrid_index_ofn ( dim_num, order_1d, order_nd );
//
//  Only keep those points which first appear on this level.
//  If you keep a point, it is necessary to rescale each of its components
//  so that we save the coordinates as they apply on the final grid.
//
      for ( point = 0; point < order_nd; point++ )
      {
        test = true;
        for ( dim = 0; dim < dim_num; dim++ )
        {
          if ( grid_index2[dim+point*dim_num] % 2 == 0 )
          {
            test = false;
          }
        }

        if ( test )
        {

          if ( point_num <= point_num2 )
          {
            cout << "LEVELS_INDEX_OFN - Fatal error!\n";
            cout << "  Exceeding maximum point index POINT_NUM = " 
                 << point_num << "\n";
            exit ( 1 );
          }

          for ( dim = 0; dim < dim_num; dim++ )
          {
            grid_base[dim+point_num2*dim_num] = order_1d[dim];

            grid_index[dim+point_num2*dim_num] =
              i4_power ( 2, level_max - level_1d[dim] )
              * grid_index2[dim+point*dim_num];
          }
          point_num2 = point_num2 + 1;
        }
      }

      delete [] grid_index2;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] order_1d;

  if ( point_num2 < point_num )
  {
    cout << "\n";
    cout << "LEVELS_INDEX_OFN - Fatal error!\n";
    cout << "  Set fewer points than POINT_NUM = " << point_num << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void levels_index_onn ( int dim_num, int level_max, int point_num, 
  int grid_index [], int grid_base[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_ONN indexes a sparse grid made from ONN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of OPEN NON NESTED 1D quadrature rules.
//
//    ONN rules include Gauss Laguerre.
//
//    The sparse grid is the logical sum of product grids with total LEVEL 
//    between LEVEL_MIN and LEVEL_MAX.
//
//    The necessary dimensions of GRID_INDEX can be determined by 
//    calling LEVELS_INDEX_SIZE_ONN first.
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
//    Output, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of 
//    point indices, representing a subset of the product grid of level 
//    LEVEL_MAX, representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
//    Output, int GRID_BASE[DIM_NUM*POINT_NUM], a list of 
//    the orders of the rules associated with each point and dimension.
//
{
  int dim;
  int factor;
  int *grid_base2;
  int *grid_index2;
  int h;
  int j;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int *order_1d;
  int order_nd;
  int point;
  int point_num2;
  int t;
//
//  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
//
  point_num2 = 0;

  level_min = i4_max ( 0, level_max + 1 - dim_num );
  
  grid_base2 = new int[dim_num];
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

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
      level_to_order_open ( dim_num, level_1d, order_1d );

      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_base2[dim] = order_1d[dim];
      }
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//
      grid_index2 = multigrid_index_onn ( dim_num, order_1d, order_nd );
//
//  Only keep those points which first appear on this level.
//
      for ( point = 0; point < order_nd; point++ )
      {
        if ( point_num <= point_num2 )
        {
          cout << "\n";
          cout << "LEVELS_INDEX_ONN - Fatal error!\n";
          cout << "  Exceeding maximum point index POINT_NUM = " 
               << point_num << "\n";
          exit ( 1 );
        }
        for ( dim = 0; dim < dim_num; dim++ )
        {
          grid_index[dim+point_num2*dim_num] = grid_index2[dim+point*dim_num];
          grid_base[dim+point_num2*dim_num] = grid_base2[dim];
        }
        point_num2 = point_num2 + 1;
      }

      delete [] grid_index2;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] grid_base2;
  delete [] level_1d;
  delete [] order_1d;

  if ( point_num2 < point_num )
  {
    cout << "\n";
    cout << "LEVELS_INDEX_ONN - Fatal error!\n";
    cout << "  Set fewer points than POINT_NUM = " << point_num << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void levels_index_own ( int dim_num, int level_max, int point_num, 
  int grid_index [], int grid_base[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_OWN indexes a sparse grid made from OWN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of OPEN WEAKLY NESTED 1D quadrature rules.
//
//    OWN rules include Gauss Hermite and Gauss Legendre.
//
//    The sparse grid is the logical sum of product grids with total LEVEL 
//    between LEVEL_MIN and LEVEL_MAX.
//
//    The necessary dimensions of GRID_INDEX can be determined by 
//    calling LEVELS_INDEX_SIZE_OWN first.
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
//    Output, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of 
//    point indices, representing a subset of the product grid of level 
//    LEVEL_MAX, representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
//    Output, int GRID_BASE[DIM_NUM*POINT_NUM], a list of 
//    the orders of the rules associated with each point and dimension.
//
{
  int dim;
  int factor;
  int *grid_base2;
  int *grid_index2;
  int *grid_level;
  int h;
  int j;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int *order_1d;
  int order_nd;
  int point;
  int point_num2;
  int t;
//
//  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
//
  point_num2 = 0;

  if ( dim_num == 1 )
  {
    level_min = level_max;
  }
  else
  {
    level_min = 0;
  }
  
  grid_base2 = new int[dim_num];
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

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
      level_to_order_open ( dim_num, level_1d, order_1d );

      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_base2[dim] = ( order_1d[dim] - 1 ) / 2;
      }
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//
      grid_index2 = multigrid_index_own ( dim_num, order_1d, order_nd );
//
//  Determine the first level of appearance of each of the points.
//  This allows us to flag certain points as being repeats of points
//  generated on a grid of lower level.  
//
//  This is SLIGHTLY tricky.
//
      grid_level = index_level_own ( level, level_max, dim_num, order_nd, 
        grid_index2, grid_base2 );
//
//  Only keep those points which first appear on this level.
//
      for ( point = 0; point < order_nd; point++ )
      {
        if ( grid_level[point] == level )
        {
          if ( point_num <= point_num2 )
          {
            cout << "\n";
            cout << "LEVELS_INDEX_OWN - Fatal error!\n";
            cout << "  Exceeding maximum point index POINT_NUM = " 
                 << point_num << "\n";
            exit ( 1 );
          }
          for ( dim = 0; dim < dim_num; dim++ )
          {
            grid_index[dim+point_num2*dim_num] =
              grid_index2[dim+point*dim_num];
            grid_base[dim+point_num2*dim_num] = grid_base2[dim];
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
  delete [] grid_base2;
  delete [] level_1d;
  delete [] order_1d;

  if ( point_num2 < point_num )
  {
    cout << "\n";
    cout << "LEVELS_INDEX_OWN - Fatal error!\n";
    cout << "  Set fewer points than POINT_NUM = " << point_num << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

int levels_index_size ( int dim_num, int level_max, int rule )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_SIZE sizes a sparse grid.
//
//  Discussion:
//
//    The sparse grid is the logical sum of product grids with total LEVEL 
//    between LEVEL_MIN and LEVEL_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 December 2009
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
//    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
//
//    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, integer ( kind = 4 ) RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Output, int  LEVELS_INDEX_SIZE, the total number of unique 
//    points in the grids.
//
{
  int point_num;

  if ( rule == 1 )
  {
    point_num = sparse_grid_cc_size ( dim_num, level_max );
  }
  else if ( 2 <= rule && rule <= 4 )
  {
    point_num = sparse_grid_ofn_size ( dim_num, level_max );
  }
  else if ( 5 <= rule && rule <= 6 )
  {
    point_num = levels_index_size_own ( dim_num, level_max );
  }
  else if ( 7 == rule )
  {
    point_num = levels_index_size_onn ( dim_num, level_max );
  }
  else
  {
    point_num = -1;
    cout << "\n";
    cout << "LEVELS_INDEX_SIZE - Fatal error!\n";
    cout << "  Unrecognized value of RULE = " << rule << "\n";
    exit ( 1 );
  }

  return point_num;
}
//****************************************************************************80

int levels_index_size_cfn ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_SIZE_CFN sizes a sparse grid made from CFN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of CLOSED FULLY NESTED 1D quadrature rules.
//
//    CFN rules include Clenshaw Curtis rules.
//
//    The sparse grid is the logical sum of product grids with total LEVEL 
//    between LEVEL_MIN and LEVEL_MAX.
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
//    Output, int LEVELS_INDEX_SIZE_CFN, the number of points in the grid.
//
{
  int dim;
  int factor;
  int *grid_index;
  int *grid_level;
  int h;
  int j;
  int level;
  int *level_1d;
  bool more;
  int *order_1d;
  int order_max;
  int order_nd;
  int point;
  int point_num;
  int t;
//
//  Special case.
//
  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  The outer loop generates LEVELs from 0 to LEVEL_MAX.
//
  point_num = 0;

  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( level = 0; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition that adds up to LEVEL.
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
      grid_index = multigrid_index_cfn ( dim_num, order_1d, order_nd );
//
//  Adjust these grid indices to reflect LEVEL_MAX.
//
      multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, 
        grid_index );
//
//  Determine the first level of appearance of each of the points.
//
      grid_level = abscissa_level_closed_nd ( level_max, dim_num, order_nd, 
        grid_index );
//
//  Only keep those points which first appear on this level.
//
      for ( point = 0; point < order_nd; point++ )
      {
        if ( grid_level[point] == level )
        {
          point_num = point_num + 1;
        }
      }

      delete [] grid_index;
      delete [] grid_level;

      if ( !more )
      {
        break;
      }
    }
  }

  delete [] level_1d;
  delete [] order_1d;

  return point_num;
}
//****************************************************************************80

int levels_index_size_onn ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_SIZE_ONN sizes a sparse grid made from ONN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of OPEN NON-NESTED 1D quadrature rules.
//
//    ONN rules include Gauss Laguerre.
//
//    The sparse grid is the logical sum of product grids with total LEVEL 
//    between LEVEL_MIN and LEVEL_MAX.
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
//    Output, int LEVELS_INDEX_SIZE_ONN, the number of points in the grid.
//
{
  int dim;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int num;
  int *order_1d;
  int order_nd;
  int point_num;
  int t;
//
//  Special case.
//
  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  The outer loop generates LEVELs from 0 to LEVEL_MAX.
//
  point_num = 0;

  level_min = i4_max ( 0, level_max + 1 - dim_num );
  
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition that adds up to LEVEL.
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
      level_to_order_open ( dim_num, level_1d, order_1d );

      point_num = point_num + i4vec_product ( dim_num, order_1d );

      if ( !more )
      {
        break;
      }
    }
  }

  delete [] level_1d;
  delete [] order_1d;

  return point_num;
}
//****************************************************************************80

int levels_index_size_own ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_SIZE_OWN sizes a sparse grid made from OWN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of OPEN WEAKLY NESTED 1D quadrature rules.
//
//    OWN rules include Gauss Hermite and Gauss Legendre.
//
//    The sparse grid is the logical sum of product grids with total LEVEL 
//    between LEVEL_MIN and LEVEL_MAX.
//
//    Oddly enough, in order to count the number of points, we will
//    behave as though LEVEL_MIN was zero.  This is because our computation
//    concentrates on throwing away all points generated at lower levels,
//    but, in fact, if we start at a nonzero level, we need to include
//    on that level all the points that would have been generated on lower
//    levels.
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
//    Output, int LEVELS_INDEX_SIZE_OWN, the number of points in the grid.
//
{
  int dim;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int num;
  int *order_1d;
  int order_nd;
  int point_num;
  int t;
//
//  Special case.
//
  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
//
//  The normal definition of LEVEL_MIN:
//
//   level_min = max ( 0, level_max + 1 - dim_num )
//
//  Our somewhat artificial temporary local definition of LEVEL_MIN:
//
  if ( dim_num == 1 )
  {
    level_min = level_max;
    point_num = 1;
  }
  else
  {
    level_min = 0;
    point_num = 0;
  }
  
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition that adds up to LEVEL.
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
      level_to_order_open ( dim_num, level_1d, order_1d );

      for ( dim = 0; dim < dim_num; dim++ )
      {
//
//  If we can reduce the level in this dimension by 1 and
//  still not go below LEVEL_MIN.
//
        if ( 1 < order_1d[dim] )
        {
          order_1d[dim] = order_1d[dim] - 1;
        }
      }
      point_num = point_num + i4vec_product ( dim_num, order_1d );

      if ( !more )
      {
        break;
      }
    }
  }

  delete [] level_1d;
  delete [] order_1d;

  return point_num;
}
//****************************************************************************80

void lg_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    LG_ABSCISSA sets abscissas for multidimensional Gauss-Laguerre quadrature.
//
//  Discussion:
//
//    The "nesting" as it occurs for Gauss-Laguerre sparse grids simply
//    involves the use of a specified set of permissible orders for the
//    rule.  
//
//    The X array lists the (complete) Gauss-Legendre abscissas for rules 
//    of order 1, 3, 7, 15, 31, 63 or 127, in order. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the index of the abscissa
//    from the rule, for each dimension and point.
//
//    Input, int GRID_BASE[DIM_NUM], the number of points used in the 
//    rule for a given dimension.
//
//    Output, double GRID_POINT[DIM_NUM], the grid points of abscissas.
//
{
  int dim;
  int level;
  int point;
  int pointer;
  int skip[8] = { 0, 1, 4, 11, 26, 57, 120, 247 };
  double x[247] = {
    1.0E+00, 
    0.415774556783479083311533873128E+00, 
    0.229428036027904171982205036136E+01, 
    0.628994508293747919686641576551E+01, 
    0.193043676560362413838247885004E+00, 
    0.102666489533919195034519944317E+01, 
    0.256787674495074620690778622666E+01, 
    0.490035308452648456810171437810E+01, 
    0.818215344456286079108182755123E+01, 
    0.127341802917978137580126424582E+02, 
    0.193957278622625403117125820576E+02, 
    0.933078120172818047629030383672E-01, 
    0.492691740301883908960101791412E+00, 
    0.121559541207094946372992716488E+01, 
    0.226994952620374320247421741375E+01, 
    0.366762272175143727724905959436E+01, 
    0.542533662741355316534358132596E+01, 
    0.756591622661306786049739555812E+01, 
    0.101202285680191127347927394568E+02, 
    0.131302824821757235640991204176E+02, 
    0.166544077083299578225202408430E+02, 
    0.207764788994487667729157175676E+02, 
    0.256238942267287801445868285977E+02, 
    0.314075191697539385152432196202E+02, 
    0.385306833064860094162515167595E+02, 
    0.480260855726857943465734308508E+02, 
    0.45901947621108290743496080275224E-01, 
    0.24198016382477204890408974151714E+00, 
    0.59525389422235073707330165005414E+00, 
    1.1066894995329987162111308789792E+00, 
    1.7775956928747727211593727482675E+00, 
    2.6097034152566806503893375925315E+00, 
    3.6051968023400442698805817554243E+00, 
    4.7667470844717611313629127271123E+00, 
    6.0975545671817409269925429328463E+00, 
    7.6014009492331374229360106942867E+00, 
    9.2827143134708894182536695297710E+00, 
    11.146649755619291358993815629587E+00, 
    13.199189576244998522464925028637E+00, 
    15.447268315549310075809325891801E+00, 
    17.898929826644757646725793817752E+00, 
    20.563526336715822170743048968779E+00, 
    23.451973482011858591050255575933E+00, 
    26.577081352118260459975876986478E+00, 
    29.953990872346445506951917840024E+00, 
    33.600759532902202735410313885784E+00, 
    37.539164407330440882887902558001E+00, 
    41.795830870182219981347945853330E+00, 
    46.403866806411123136029227604386E+00, 
    51.405314476797755161861461088395E+00, 
    56.854992868715843620511922055660E+00, 
    62.826855908786321453677523304806E+00, 
    69.425277191080345623322251656443E+00, 
    76.807047763862732837609972285484E+00, 
    85.230358607545669169387065607043E+00, 
    95.188939891525629981308606853957E+00, 
    107.95224382757871475002440117666E+00, 
    0.22768893732576153785994330248562E-01, 
    0.11998325242727824715771416426383E+00, 
    0.29494185444770149577427738517405E+00, 
    0.54779087896237725363865073775856E+00, 
    0.87869061179931901673895567052285E+00, 
    1.2878464335919706302309207788611E+00, 
    1.7755123815388553763979463268728E+00, 
    2.3419925567085989256055628337716E+00, 
    2.9876423223246473939976731053629E+00, 
    3.7128695992018000346299637413422E+00, 
    4.5181363349503584391105568561550E+00, 
    5.4039601781825946286902599782736E+00, 
    6.3709163787865330220392250891777E+00, 
    7.4196399339311711154888493199004E+00, 
    8.5508280008403328312589048722235E+00, 
    9.7652425999245366807004592977996E+00, 
    11.063713635140661736220550410604E+00, 
    12.447142262356492749798687569289E+00, 
    13.916504641057818562912967008183E+00, 
    15.472856110036296424777143607779E+00, 
    17.117335833863588753116900303886E+00, 
    18.851171974154856850873483787506E+00, 
    20.675687448056515660377265667433E+00, 
    22.592306346311528381292277759986E+00, 
    24.602561094972638883700642760037E+00, 
    26.708100458737343969779087998829E+00, 
    28.910698500451382640177718103234E+00, 
    31.212264631175912885477773820802E+00, 
    33.614854909101154836598842888345E+00, 
    36.120684774484823056306328740825E+00, 
    38.732143442933582145626041607663E+00, 
    41.451810222318741191114726181363E+00, 
    44.282473071479233839358857134636E+00, 
    47.227149784295686898935095231536E+00, 
    50.289112264240695761749021839419E+00, 
    53.471914456788652808348280619542E+00, 
    56.779424636342062213099781057119E+00, 
    60.215862909019862886417550114424E+00, 
    63.785845004235974631701139601836E+00, 
    67.494433702293885830374325695045E+00, 
    71.347199604295266286654803376075E+00, 
    75.350293425653234254290504744279E+00, 
    79.510532629986309149555391354778E+00, 
    83.835506080872257843339817658508E+00, 
    88.333701570354369086112766326498E+00, 
    93.014662728558547405303399037100E+00, 
    97.889184147578140043386727677112E+00, 
    102.96955690741381650783952746778E+00, 
    108.26988161961595392226350967206E+00, 
    113.80647350287462738934485955901E+00, 
    119.59839538830458666962452963285E+00, 
    125.66817255856119431291196303280E+00, 
    132.04277272091165746585590583045E+00, 
    138.75498418103789078167590567526E+00, 
    145.84541318313540358283994248439E+00, 
    153.36548459497863623710815962660E+00, 
    161.38215194813761243562172669592E+00, 
    169.98570600665839438795175301156E+00, 
    179.30366247401580910251827858515E+00, 
    189.52789596532475473668721332981E+00, 
    200.97521159924656741628671841018E+00, 
    214.25368536638788642698056296400E+00, 
    230.93465747089703971246562985079E+00, 
    0.11339635298518611691893169631306E-01, 
    0.59749753435726620281348237057387E-01, 
    0.14685098690746167612388223687431E+00, 
    0.27267590735859553131378008278900E+00, 
    0.43724600644192665554577035869932E+00, 
    0.64058688222566929533576416399983E+00, 
    0.88272968639058364481487653650042E+00, 
    1.1637114160166537661560584700951E+00, 
    1.4835750152834613891313584861012E+00, 
    1.8423694351613565380686320809853E+00, 
    2.2401496839579024244513315656522E+00, 
    2.6769768780141303692167869961238E+00, 
    3.1529182957082825565771508308846E+00, 
    3.6680474360304752540226339926515E+00, 
    4.2224440823301888455977876667425E+00, 
    4.8161943715870502475665535087286E+00, 
    5.4493908694559416755862178908416E+00, 
    6.1221326512997254193944584763155E+00, 
    6.8345253894122668112237994973336E+00, 
    7.5866814466367472174205986836847E+00, 
    8.3787199765932725254842120659452E+00, 
    9.2107670307426558777922506102445E+00, 
    10.082955672528643809166439353647E+00, 
    10.995426098858125429803147358780E+00, 
    11.948325769197725997610605127857E+00, 
    12.941809542585531053723381098192E+00, 
    13.976039822878506520014405668679E+00, 
    15.051186712579523631574796365435E+00, 
    16.167428175612852922977395051768E+00, 
    17.324950209443673446561163712616E+00, 
    18.523947026965688560811711309349E+00, 
    19.764621248611504104071669386884E+00, 
    21.047184105173183606877044020054E+00, 
    22.371855651855542817648123918101E+00, 
    23.738864994122497183652313788712E+00, 
    25.148450525937368234077278385644E+00, 
    26.600860181041749607253384279755E+00, 
    28.096351697964619201753961292129E+00, 
    29.635192899504178910610227138642E+00, 
    31.217661987479759144214467152615E+00, 
    32.844047853610430460522951341338E+00, 
    34.514650407441149149105635947422E+00, 
    36.229780922306804019615388508885E+00, 
    37.989762400399956435968780140278E+00, 
    39.794929958089961778396437141707E+00, 
    41.645631232730180705153990897484E+00, 
    43.542226812286859549950892993822E+00, 
    45.485090689228791137996151336673E+00, 
    47.474610740231964719468766599146E+00, 
    49.511189233379087716728884584381E+00, 
    51.595243364671244443182771266934E+00, 
    53.727205825819316758288140069145E+00, 
    55.907525405447553305830605991732E+00, 
    58.136667626022439197077526025660E+00, 
    60.415115419018590295707192053805E+00, 
    62.743369841051809700207126742685E+00, 
    65.121950833949996311956025417139E+00, 
    67.551398031997886314411872443149E+00, 
    70.032271619884584511229871192030E+00, 
    72.565153245206849090888669416801E+00, 
    75.150646989739935299354362325096E+00, 
    77.789380404085816000647405462136E+00, 
    80.482005610750729205803962926758E+00, 
    83.229200481195914886796120019048E+00, 
    86.031669892953582966798238732643E+00, 
    88.890147073512051099652518544282E+00, 
    91.805395038358177994971250170499E+00, 
    94.778208131331583205387031034825E+00, 
    97.809413676305116411054110115424E+00, 
    100.89987375017285940371939762172E+00, 
    104.05048708821598934704076845022E+00, 
    107.26219113414600428423116401414E+00, 
    110.53596424851500530602771351277E+00, 
    113.87282809075839485348376187652E+00, 
    117.27385019192517774095477886379E+00, 
    120.74014673718880106173978002719E+00, 
    124.27288557955698354259506446928E+00, 
    127.87328950885942645093841745425E+00, 
    131.54263980314366921809377742137E+00, 
    135.28228009311836970132738106369E+00, 
    139.09362057432970013964422086977E+00, 
    142.97814260643601776808227753574E+00, 
    146.93740374437366549441080969072E+00, 
    150.97304325252187127492511437460E+00, 
    155.08678816034612572229641420609E+00, 
    159.28045992663288235401956989889E+00, 
    163.55598178957571104015967182053E+00, 
    167.91538689194360134245547184721E+00, 
    172.36082728473812536838156191681E+00, 
    176.89458392960192176311674993508E+00, 
    181.51907784036813069227528834025E+00, 
    186.23688252828112373861202530357E+00, 
    191.05073794450929196790836610789E+00, 
    195.96356614879879837839002542988E+00, 
    200.97848897600025153696475526130E+00, 
    206.09884802468871112127283042753E+00, 
    211.32822735671655260572377256981E+00, 
    216.67047937658230323477089465777E+00, 
    222.12975445929687246267304963754E+00, 
    227.71053502072232419089132431317E+00, 
    233.41767488282602453367775322563E+00, 
    239.25644498830308620018749667089E+00, 
    245.23258677871567172531254018984E+00, 
    251.35237488718128030005500991754E+00, 
    257.62269123792061413076191882313E+00, 
    264.05111322908240551754377241831E+00, 
    270.64601945722796749299111718606E+00, 
    277.41671750163651071798388218104E+00, 
    284.37359974220870326674402873120E+00, 
    291.52833521346495719581282021650E+00, 
    298.89410837028248600878895615414E+00, 
    306.48591978262611320418112423947E+00, 
    314.32096986471177487400007507615E+00, 
    322.41915589128679683349440361344E+00, 
    330.80372663802405651933847334878E+00, 
    339.50216127832433747735367595958E+00, 
    348.54737559472697355480761787441E+00, 
    357.97942028029845454049007443090E+00, 
    367.84794520076004578858341422871E+00, 
    378.21590623135532818332979188889E+00, 
    389.16539141251004101579475325153E+00, 
    400.80729331451702589996361286427E+00, 
    413.29853681779384418008260081859E+00, 
    426.87579153663675538288509017051E+00, 
    441.93085485310841412460309271842E+00, 
    459.21804639888429981971267313224E+00, 
    480.69378263388373859704269229304E+00   
  };

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( grid_base[dim] < 1 )
    {
      cout << "\n";
      cout << "LG_ABSCISSA - Fatal error!\n";
      cout << "  Some base values are less than 1.\n";
      exit ( 1 );
    }
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 127 < grid_base[dim] )
    {
      cout << "\n";
      cout << "LG_ABSCISSA - Fatal error!\n";
      cout << "  Some base values are greater than 127.\n";
      exit ( 1 );
    }
  }

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level = i4_log_2 ( grid_base[dim] + 1 ) - 1;

      pointer = skip[level] + grid_index[dim+point*dim_num];

     if ( pointer < 1 || 247 < pointer )
     {
        cout << "\n";
        cout << "LG_ABSCISSA - Fatal error!\n";
        cout << "  POINTER out of bounds.\n";
        cout << "  POINTER    = " << pointer << "\n";
        cout << "  POINT      = " << point << "\n";
        cout << "  DIM        = " << dim << "\n";
        cout << "  GRID_BASE  = " << grid_base[dim] << "\n";
        cout << "  LEVEL      = " << level << "\n";
        cout << "  GRID_INDEX = " << grid_index[dim+point*dim_num] << "\n";
        exit ( 1 );
      }

      grid_point[dim+point*dim_num] = x[pointer-1];
    }
  }

  return;
}
//****************************************************************************80

double *lg_weights ( int order )

//****************************************************************************80
//
//  Purpose:
//
//    LG_WEIGHTS returns weights for certain Gauss-Laguerre quadrature rules.
//
//  Discussion:
//
//    The allowed orders are 1, 3, 7, 15, 31, 63 and 127.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 October 2007
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
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
//
//    Output, double WEIGHT[ORDER], the weights.
//    The weights are positive, symmetric and should sum to 1.
//
{
  double *weight;

  weight = new double[order];

  if ( order == 1 )
  {
    weight[1-1] = 1.0E+00;
  }
  else if ( order == 3 )
  {
    weight[1-1] =  0.711093009929173015449590191143E+00;
    weight[2-1] =  0.278517733569240848801444888457E+00;
    weight[3-1] =  0.103892565015861357489649204007E-01;
  }
  else if ( order == 7 )
  {
    weight[1-1] =  0.409318951701273902130432880018E+00;
    weight[2-1] =  0.421831277861719779929281005417E+00;
    weight[3-1] =  0.147126348657505278395374184637E+00;
    weight[4-1] =  0.206335144687169398657056149642E-01;
    weight[5-1] =  0.107401014328074552213195962843E-02;
    weight[6-1] =  0.158654643485642012687326223234E-04;
    weight[7-1] =  0.317031547899558056227132215385E-07;
  }
  else if ( order == 15 )
  {
    weight[1-1] =  0.218234885940086889856413236448E+00;
    weight[2-1] =  0.342210177922883329638948956807E+00;
    weight[3-1] =  0.263027577941680097414812275022E+00;
    weight[4-1] =  0.126425818105930535843030549378E+00;
    weight[5-1] =  0.402068649210009148415854789871E-01;
    weight[6-1] =  0.856387780361183836391575987649E-02;
    weight[7-1] =  0.121243614721425207621920522467E-02;
    weight[8-1] =  0.111674392344251941992578595518E-03;
    weight[9-1] =  0.645992676202290092465319025312E-05;
    weight[10-1] = 0.222631690709627263033182809179E-06;
    weight[11-1] = 0.422743038497936500735127949331E-08;
    weight[12-1] = 0.392189726704108929038460981949E-10;
    weight[13-1] = 0.145651526407312640633273963455E-12;
    weight[14-1] = 0.148302705111330133546164737187E-15;
    weight[15-1] = 0.160059490621113323104997812370E-19;
  }
  else if ( order == 31 )
  {
    weight[  1-1] =   0.11252789550372583820847728082801E+00;
    weight[  2-1] =   0.21552760818089123795222505285045E+00;
    weight[  3-1] =   0.23830825164569654731905788089234E+00;
    weight[  4-1] =   0.19538830929790229249915303390711E+00;
    weight[  5-1] =   0.12698283289306190143635272904602E+00;
    weight[  6-1] =   0.67186168923899300670929441993508E-01;
    weight[  7-1] =   0.29303224993879487404888669311974E-01;
    weight[  8-1] =   0.10597569915295736089529380314433E-01;
    weight[  9-1] =   0.31851272582386980320974842433019E-02;
    weight[ 10-1] =   0.79549548307940382922092149012477E-03;
    weight[ 11-1] =   0.16480052126636687317862967116412E-03;
    weight[ 12-1] =   0.28229237864310816393860971468993E-04;
    weight[ 13-1] =   0.39802902551008580387116174900106E-05;
    weight[ 14-1] =   0.45931839841801061673729694510289E-06;
    weight[ 15-1] =   0.43075545187731100930131457465897E-07;
    weight[ 16-1] =   0.32551249938271570855175749257884E-08;
    weight[ 17-1] =   0.19620246675410594996247151593142E-09;
    weight[ 18-1] =   0.93190499086617587129534716431331E-11;
    weight[ 19-1] =   0.34377541819411620520312597898311E-12;
    weight[ 20-1] =   0.96795247130446716997405035776206E-14;
    weight[ 21-1] =   0.20368066110115247398010624219291E-15;
    weight[ 22-1] =   0.31212687280713526831765358632585E-17;
    weight[ 23-1] =   0.33729581704161052453395678308350E-19;
    weight[ 24-1] =   0.24672796386616696011038363242541E-21;
    weight[ 25-1] =   0.11582201904525643634834564576593E-23;
    weight[ 26-1] =   0.32472922591425422434798022809020E-26;
    weight[ 27-1] =   0.49143017308057432740820076259666E-29;
    weight[ 28-1] =   0.34500071104808394132223135953806E-32;
    weight[ 29-1] =   0.87663710117162041472932760732881E-36;
    weight[ 30-1] =   0.50363643921161490411297172316582E-40;
    weight[ 31-1] =   0.19909984582531456482439549080330E-45;
  }
  else if ( order == 63 )
  {
    weight[  1-1] =   0.57118633213868979811587283390476E-01;
    weight[  2-1] =   0.12067476090640395283319932036351E+00;
    weight[  3-1] =   0.15925001096581873723870561096472E+00;
    weight[  4-1] =   0.16875178327560799234596192963585E+00;
    weight[  5-1] =   0.15366641977668956696193711310131E+00;
    weight[  6-1] =   0.12368770614716481641086652261948E+00;
    weight[  7-1] =   0.89275098854848671545279150057422E-01;
    weight[  8-1] =   0.58258485446105944957571825725160E-01;
    weight[  9-1] =   0.34546657545992580874717085812508E-01;
    weight[ 10-1] =   0.18675685985714656798286552591203E-01;
    weight[ 11-1] =   0.92233449044093536528490075241649E-02;
    weight[ 12-1] =   0.41671250684839592762582663470209E-02;
    weight[ 13-1] =   0.17238120299900582715386728541955E-02;
    weight[ 14-1] =   0.65320845029716311169340559359043E-03;
    weight[ 15-1] =   0.22677644670909586952405173207471E-03;
    weight[ 16-1] =   0.72127674154810668410750270234861E-04;
    weight[ 17-1] =   0.21011261180466484598811536851241E-04;
    weight[ 18-1] =   0.56035500893357212749181536071292E-05;
    weight[ 19-1] =   0.13673642785604888017836641282292E-05;
    weight[ 20-1] =   0.30507263930195817240736097189550E-06;
    weight[ 21-1] =   0.62180061839309763559981775409241E-07;
    weight[ 22-1] =   0.11566529551931711260022448996296E-07;
    weight[ 23-1] =   0.19614588267565478081534781863335E-08;
    weight[ 24-1] =   0.30286171195709411244334756404054E-09;
    weight[ 25-1] =   0.42521344539400686769012963452599E-10;
    weight[ 26-1] =   0.54202220578073819334698791381873E-11;
    weight[ 27-1] =   0.62627306838597672554166850420603E-12;
    weight[ 28-1] =   0.65474443156573322992307089591924E-13;
    weight[ 29-1] =   0.61815575808729181846302500000047E-14;
    weight[ 30-1] =   0.52592721363507381404263991342633E-15;
    weight[ 31-1] =   0.40230920092646484015391506025408E-16;
    weight[ 32-1] =   0.27600740511819536505013824207729E-17;
    weight[ 33-1] =   0.16936946756968296053322009855265E-18;
    weight[ 34-1] =   0.92689146872177087314963772462726E-20;
    weight[ 35-1] =   0.45093739060365632939780140603959E-21;
    weight[ 36-1] =   0.19435162876132376573629962695374E-22;
    weight[ 37-1] =   0.73926270895169207037999639194513E-24;
    weight[ 38-1] =   0.24714364154434632615980126000066E-25;
    weight[ 39-1] =   0.72288649446741597655145390616476E-27;
    weight[ 40-1] =   0.18407617292614039362985209905608E-28;
    weight[ 41-1] =   0.40583498566841960105759537058880E-30;
    weight[ 42-1] =   0.77000496416438368114463925286343E-32;
    weight[ 43-1] =   0.12488505764999334328843314866038E-33;
    weight[ 44-1] =   0.17185000226767010697663950619912E-35;
    weight[ 45-1] =   0.19896372636672396938013975755522E-37;
    weight[ 46-1] =   0.19199671378804058267713164416870E-39;
    weight[ 47-1] =   0.15278588285522166920459714708240E-41;
    weight[ 48-1] =   0.99054752688842142955854138884590E-44;
    weight[ 49-1] =   0.51597523673029211884228858692990E-46;
    weight[ 50-1] =   0.21249846664084111245693912887783E-48;
    weight[ 51-1] =   0.67903852766852910591172042494884E-51;
    weight[ 52-1] =   0.16466654148296177467908300517887E-53;
    weight[ 53-1] =   0.29509065402691055027053659375033E-56;
    weight[ 54-1] =   0.37838420647571051984882241014675E-59;
    weight[ 55-1] =   0.33358130068542431878174667995217E-62;
    weight[ 56-1] =   0.19223461022273880981363303073329E-65;
    weight[ 57-1] =   0.67812696961083016872779388922288E-69;
    weight[ 58-1] =   0.13404752802440604607620468935693E-72;
    weight[ 59-1] =   0.13109745101805029757648048223928E-76;
    weight[ 60-1] =   0.52624863881401787388694579143866E-81;
    weight[ 61-1] =   0.63780013856587414257760666006511E-86;
    weight[ 62-1] =   0.12997078942372924566347473916943E-91;
    weight[ 63-1] =   0.10008511496968754063443740168421E-98;
  }
  else if ( order == 127 )
  {
    weight[  1-1] =   0.28773246692000124355770010301506E-01;
    weight[  2-1] =   0.63817468175134649363480949265236E-01;
    weight[  3-1] =   0.91919669721570571389864194652717E-01;
    weight[  4-1] =   0.11054167914413766381245463002967E+00;
    weight[  5-1] =   0.11879771633375850188328329422643E+00;
    weight[  6-1] =   0.11737818530052695148804451630074E+00;
    weight[  7-1] =   0.10819305984180551488335145581193E+00;
    weight[  8-1] =   0.93827075290489628080377261401107E-01;
    weight[  9-1] =   0.76966450960588843995822485928431E-01;
    weight[ 10-1] =   0.59934903912939714332570730063476E-01;
    weight[ 11-1] =   0.44417742073889001371708316272923E-01;
    weight[ 12-1] =   0.31385080966252320983009372215062E-01;
    weight[ 13-1] =   0.21172316041924506411370709025015E-01;
    weight[ 14-1] =   0.13650145364230541652171185564626E-01;
    weight[ 15-1] =   0.84172852710599172279366657385445E-02;
    weight[ 16-1] =   0.49674990059882760515912858620175E-02;
    weight[ 17-1] =   0.28069903895001884631961957446400E-02;
    weight[ 18-1] =   0.15192951003941952460445341057817E-02;
    weight[ 19-1] =   0.78789028751796084086217287140548E-03;
    weight[ 20-1] =   0.39156751064868450584507324648999E-03;
    weight[ 21-1] =   0.18652434268825860550093566260060E-03;
    weight[ 22-1] =   0.85173160415576621908809828160247E-04;
    weight[ 23-1] =   0.37285639197853037712145321577724E-04;
    weight[ 24-1] =   0.15648416791712993947447805296768E-04;
    weight[ 25-1] =   0.62964340695224829035692735524979E-05;
    weight[ 26-1] =   0.24288929711328724574541379938222E-05;
    weight[ 27-1] =   0.89824607890051007201922871545035E-06;
    weight[ 28-1] =   0.31844174740760353710742966328091E-06;
    weight[ 29-1] =   0.10821272905566839211861807542741E-06;
    weight[ 30-1] =   0.35245076750635536015902779085340E-07;
    weight[ 31-1] =   0.11001224365719347407063839761738E-07;
    weight[ 32-1] =   0.32904079616717932125329343003261E-08;
    weight[ 33-1] =   0.94289145237889976419772700772988E-09;
    weight[ 34-1] =   0.25882578904668318184050195309296E-09;
    weight[ 35-1] =   0.68047437103370762630942259017560E-10;
    weight[ 36-1] =   0.17131398805120837835399564475632E-10;
    weight[ 37-1] =   0.41291744524052865469443922304935E-11;
    weight[ 38-1] =   0.95264189718807273220707664873469E-12;
    weight[ 39-1] =   0.21032604432442425932962942047474E-12;
    weight[ 40-1] =   0.44427151938729352860940434285789E-13;
    weight[ 41-1] =   0.89760500362833703323319846405449E-14;
    weight[ 42-1] =   0.17341511407769287074627948346848E-14;
    weight[ 43-1] =   0.32028099548988356631494379835210E-15;
    weight[ 44-1] =   0.56531388950793682022660742095189E-16;
    weight[ 45-1] =   0.95329672799026591234588044025896E-17;
    weight[ 46-1] =   0.15353453477310142565288509437552E-17;
    weight[ 47-1] =   0.23608962179467365686057842132176E-18;
    weight[ 48-1] =   0.34648742794456611332193876653230E-19;
    weight[ 49-1] =   0.48515241897086461320126957663545E-20;
    weight[ 50-1] =   0.64786228633519813428137373790678E-21;
    weight[ 51-1] =   0.82476020965403242936448553126316E-22;
    weight[ 52-1] =   0.10005361880214719793491658282977E-22;
    weight[ 53-1] =   0.11561395116207304954233181263632E-23;
    weight[ 54-1] =   0.12719342731167922655612134264961E-24;
    weight[ 55-1] =   0.13316584714165372967340004160814E-25;
    weight[ 56-1] =   0.13261218454678944033646108509198E-26;
    weight[ 57-1] =   0.12554995447643949807286074138324E-27;
    weight[ 58-1] =   0.11294412178579462703240913107219E-28;
    weight[ 59-1] =   0.96491020279562119228500608131696E-30;
    weight[ 60-1] =   0.78241846768302099396733076955632E-31;
    weight[ 61-1] =   0.60181503542219626658249939076636E-32;
    weight[ 62-1] =   0.43882482704961741551510518054138E-33;
    weight[ 63-1] =   0.30314137647517256304035802501863E-34;
    weight[ 64-1] =   0.19826016543944539545224676057020E-35;
    weight[ 65-1] =   0.12267623373665926559013654872402E-36;
    weight[ 66-1] =   0.71763931692508888943812834967620E-38;
    weight[ 67-1] =   0.39659378833836963584113716149270E-39;
    weight[ 68-1] =   0.20688970553868040099581951696677E-40;
    weight[ 69-1] =   0.10179587017979517245268418427523E-41;
    weight[ 70-1] =   0.47200827745986374625714293679649E-43;
    weight[ 71-1] =   0.20606828985553374825744353490744E-44;
    weight[ 72-1] =   0.84627575907305987245899032156188E-46;
    weight[ 73-1] =   0.32661123687088798658026998931647E-47;
    weight[ 74-1] =   0.11833939207883162380564134612682E-48;
    weight[ 75-1] =   0.40211209123895013807243250164050E-50;
    weight[ 76-1] =   0.12799824394111125389430292847476E-51;
    weight[ 77-1] =   0.38123877747548846504399051365162E-53;
    weight[ 78-1] =   0.10612057542701156767898551949650E-54;
    weight[ 79-1] =   0.27571446947200403594113572720812E-56;
    weight[ 80-1] =   0.66772544240928492881306904862856E-58;
    weight[ 81-1] =   0.15052438383868234954068178600268E-59;
    weight[ 82-1] =   0.31538986800113758526689068500772E-61;
    weight[ 83-1] =   0.61326614299483180785237418887960E-63;
    weight[ 84-1] =   0.11048510030324810567549119229368E-64;
    weight[ 85-1] =   0.18410563538091348076979665543900E-66;
    weight[ 86-1] =   0.28323926570052832195543883237652E-68;
    weight[ 87-1] =   0.40154409843763655508670978777418E-70;
    weight[ 88-1] =   0.52351530215683708779772201956106E-72;
    weight[ 89-1] =   0.62634476665005100555787696642851E-74;
    weight[ 90-1] =   0.68612210535666530365348093803922E-76;
    weight[ 91-1] =   0.68651298840956019297134099761855E-78;
    weight[ 92-1] =   0.62581388433728084867318704240915E-80;
    weight[ 93-1] =   0.51833271237514904046803469968027E-82;
    weight[ 94-1] =   0.38893621571918443533108973497673E-84;
    weight[ 95-1] =   0.26357711379476932781525533730623E-86;
    weight[ 96-1] =   0.16078851293917979699005509638883E-88;
    weight[ 97-1] =   0.87978042070968939637972577886624E-91;
    weight[ 98-1] =   0.43013405077495109903408697802188E-93;
    weight[ 99-1] =   0.18713435881342838527144321803729E-95;
    weight[100-1] =   0.72125744708060471675805761366523E-98;
    weight[101-1] =   0.24508746062177874383231742333023E-100;
    weight[102-1] =   0.73042094619470875777647865078327E-103;
    weight[103-1] =   0.18983290818383463537886818579820E-105;
    weight[104-1] =   0.42757400244246684123093264825902E-108;
    weight[105-1] =   0.82894681420515755691423485228897E-111;
    weight[106-1] =   0.13729432219324400013067050156048E-113;
    weight[107-1] =   0.19265464126404973222043166489406E-116;
    weight[108-1] =   0.22693344503301354826140809941334E-119;
    weight[109-1] =   0.22209290603717355061909071271535E-122;
    weight[110-1] =   0.17851087685544512662856555121755E-125;
    weight[111-1] =   0.11630931990387164467431190485525E-128;
    weight[112-1] =   0.60524443584652392290952805077893E-132;
    weight[113-1] =   0.24729569115063528647628375096400E-135;
    weight[114-1] =   0.77789065006489410364997205809045E-139;
    weight[115-1] =   0.18409738662712607039570678274636E-142;
    weight[116-1] =   0.31900921131079114970179071968597E-146;
    weight[117-1] =   0.39179487139174199737617666077555E-150;
    weight[118-1] =   0.32782158394188697053774429820559E-154;
    weight[119-1] =   0.17793590713138888062819640128739E-158;
    weight[120-1] =   0.58882353408932623157467835381214E-163;
    weight[121-1] =   0.10957236509071169877747203273886E-167;
    weight[122-1] =   0.10281621114867000898285076975760E-172;
    weight[123-1] =   0.41704725557697758145816510853967E-178;
    weight[124-1] =   0.58002877720316101774638319601971E-184;
    weight[125-1] =   0.18873507745825517106171619101120E-190;
    weight[126-1] =   0.69106601826730911682786705950895E-198;
    weight[127-1] =   0.43506813201105855628383313334402E-207;
  }
  else
  {
    cout << "\n";
    cout << "LG_WEIGHTS - Fatal error!\n";
    cout << "  Illegal value of ORDER = " << order << "\n";
    cout << "  Legal values are 1, 3, 7, 15, 31, 63 and 127.\n";
    exit ( 1 );
  }
  return weight;
}
//****************************************************************************80

double monomial_integral_hermite ( int dim_num, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_INTEGRAL_HERMITE integrates a Hermite mononomial.
//
//  Discussion:
//
//    H(d,n) = Integral ( -Infinity < x < Infinity ) 
//      x1^n1 * x2^n2...*xd^nd * exp(-x1^2-x2^2...-xd^2 ) dx
//
//    H(d,n) is 0 if any n(i) odd.
//
//    H(d,n) = product ( 1 <= i <= d ) 
//      ( (n(i)-1)!! * sqrt(pi) / 2^(n(i)/2) for all n(i) even.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2007
//
//  Author:
//
//    John Burkardt
//
// Parameters:
//
//    Input, int DIM_NUM, the dimension of the integral.
//
//    Input, int EXPON[DIM_NUM], the order of the integral.  
//    0 <= EXPON(1:DIM_NUM).
//
//    Output, double HERMITE_INTEGRAL, the value of the integral.
//
{ 
  int dim;
  double pi = 3.141592653589793;
  double value;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( expon[dim] < 0 )
    {
      value = - r8_huge ( );
      return value;
    }
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( ( expon[dim] % 2 ) == 1 )
    {
      value = 0.0;
      return value;
    }
  }

  value = 1.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    value = value * r8_factorial2 ( expon[dim] - 1 ) * sqrt ( pi ) 
      / ( double ) i4_power ( 2, expon[dim] / 2 );
  }

  return value;
}
//****************************************************************************80

double monomial_integral_laguerre ( int dim_num, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_INTEGRAL_LAGUERRE integrates a Laguerre monomial.
//
//  Discussion:
//
//    L(1,n) = Integral ( 0 <= x < Infinity ) x^n exp ( -x ) dx
//           = n!
//
//    L(d,n) = Integral ( 0 <= x(i) < Infinity ) 
//             x1^n1 * x2^n2...*xd^nd * exp(-x1-x2...-xd ) dx
//           = Product ( n(i)! ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 October 2007
//
//  Author:
//
//    John Burkardt
//
// Parameters:
//
//    Input, int DIM_NUM, the dimension of the integral.
//
//    Input, int EXPON[DIM_NUM], the order of the integral.  
//    0 <= EXPON(1:DIM_NUM).
//
//    Output, double MONOMIAL_INTEGRAL_LAGUERRE, the value of the integral.
//
{ 
  int dim;
  double value;

  value = 1.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    value = value * r8_factorial ( expon[dim] );
  }

  return value;
}
//****************************************************************************80

double monomial_integral_legendre ( int dim_num, int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_INTEGRAL_LEGENDRE integrates a Legendre monomial.
//
//  Discussion:
//
//    This routine returns the exact value of a multidimensional Legendre 
//    type integral:
//
//      integral ( -1 <= x(1:n) <= +1 ) f(x) dx
//
//    where f(x) is a monomial of the form:
//
//      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
//
//    and the exponents are nonnegative integers.  Note that the combination 
//    0^0 is treated as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Output, double MONOMIAL_INTEGRAL_LEGENDRE, the value of the integral 
//    of the monomial.
//
{
  int dim;
  double value;

  value = 0.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( ( expon[dim] % 2 ) == 1 )
    {
      return value;
    }
  }

  value = 1.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    value = value * 2.0 / ( double ) ( expon[dim] + 1 );
  }

  return value;
}
//****************************************************************************80

double monomial_quadrature ( int dim_num, int expon[], int point_num, 
  double weight[], double x[], int rule )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
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
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Input, int POINT_NUM, the number of points in the rule.
//
//    Input, double WEIGHT[POINT_NUM], the quadrature weights.
//
//    Input, double X[DIM_NUM*POINT_NUM], the quadrature points.
//
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Output, double MONOMIAL_QUADRATURE, the quadrature error.
//
{
  double exact;
  int point;
  double quad;
  double quad_error;
  double scale;
  double *value;
//
//  Get the exact value of the integral of the unscaled monomial.
//
  if ( 1 <= rule && rule <= 5 )
  {
    exact = monomial_integral_legendre ( dim_num, expon );
  }
  else if ( rule == 6 )
  {
    exact = monomial_integral_hermite ( dim_num, expon );
  }
  else if ( rule == 7 )
  {
    exact = monomial_integral_laguerre ( dim_num, expon );
  }
  else
  {
    cout << "\n";
    cout << "MONOMIAL_QUADRATURE - Fatal error!\n";
    cout << "  Unrecognized value of RULE.\n";
    exit ( 1 );
  }
//
//  Evaluate the monomial at the quadrature points.
//
  value = monomial_value ( dim_num, point_num, x, expon );
//
//  Compute the quadrature sum.
//
  quad = 0.0;
  for ( point = 0; point < point_num; point++ )
  {
    quad = quad + weight[point] * value[point];
  }
//
//  Absolute error if EXACT = 0, relative error otherwise:
//
  if ( exact == 0.0 )
  {
    quad_error = r8_abs ( quad - exact );
  }
  else
  {
    quad_error = r8_abs ( quad - exact ) / r8_abs ( exact );
  }

  delete [] value;

  return quad_error;
}
//****************************************************************************80

double *monomial_value ( int dim_num, int point_num, double x[], int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
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
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points at which the
//    monomial is to be evaluated.
//
//    Input, double X[DIM_NUM*POINT_NUM], the point coordinates.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Output, double MONOMIAL_VALUE[POINT_NUM], the value of the monomial.
//
{
  int dim;
  int point;
  double *value;

  value = new double[point_num];

  for ( point = 0; point < point_num; point++ )
  {
    value[point] = 1.0;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 0 != expon[dim] )
    {
      for ( point = 0; point < point_num; point++ )
      {
        value[point] = value[point] * pow ( x[dim+point*dim_num], expon[dim] );
      }
    }
  }

  return value;
}
//****************************************************************************80

int *multigrid_index_cfn ( int dim_num, int order_1d[], int order_nd )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_INDEX_CFN indexes a sparse grid based on CFN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of CLOSED FULLY NESTED 1D quadrature rules.
//
//    CFN rules include Clenshaw Curtis rules.
//
//    For dimension DIM, the second index of INDX may vary from 
//    0 to ORDER_1D(DIM)-1.
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

int *multigrid_index_ofn ( int dim_num, int order_1d[], int order_nd )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_INDEX_OFN indexes a sparse grid based on OFN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of OPEN FULLY NESTED 1D quadrature rules.
//
//    OFN rules include Fejer 1, Fejer 2, and Gauss Patterson rules.
//
//    For dimension DIM, the second index of INDX may vary from 
//    1 to ORDER_1D(DIM).
//
//  Modified:
//
//    26 March 2008
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
//    Input, int DIM_NUM, the spatial dimension of the points.
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
      indx[dim+p*dim_num] = a[dim] + 1;
    }
    p = p + 1;
  }

  delete [] a;

  return indx;
}
//****************************************************************************80

int *multigrid_index_onn ( int dim_num, int order_1d[], int order_nd )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_INDEX_ONN indexes a sparse grid based on ONN 1D rules.
//
//  Discussion:
//
//    The sparse grid is presumed to have been created from products
//    of OPEN NON-NESTED 1D quadrature rules.
//
//    ONN rules include Gauss Laguerre.
//
//    For dimension DIM, the number of points is ORDER_1D(DIM).
//
//    We index the points as
//      1, 2, 3, ..., ORDER_1D(DIM).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 October 2007
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
//    Input, int DIM_NUM, the spatial dimension of the points.
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
//
//  The values of A(DIM) are between 0 and ORDER_1D(DIM)-1 = N - 1 = 2 * M.
//  Subtracting M sets the range to -M to +M, as we wish.
//
    for ( dim = 0; dim < dim_num; dim++ )
    {
      indx[dim+p*dim_num] = a[dim] + 1;
    }
    p = p + 1;
  }

  delete [] a;

  return indx;
}
//****************************************************************************80

int *multigrid_index_own ( int dim_num, int order_1d[], int order_nd )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_INDEX_OWN returns an indexed multidimensional grid.
//
//  Discussion:
//
//    For dimension DIM, the number of points is ORDER_1D[DIM].
//
//    We assume that ORDER_1D[DIM] is an odd number,
//      ORDER_1D[DIM] = N = 2 * M + 1
//    so that the points have coordinates
//      -M/M, -(M-1)/M, ..., -1/M, 0/M, 1/M, 2/M, 3/M, ..., (M-1)/M, M/M.
//    and we index them as
//      -M,   -(M-1),        -1,   0,   1,   2,   3,   ...,  M-1,    M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2007
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
//    Input, int DIM_NUM, the spatial dimension of the points.
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
//
//  The values of A(DIM) are between 0 and ORDER_1D(DIM)-1 = N - 1 = 2 * M.
//  Subtracting M sets the range to -M to +M, as we wish.
//
    for ( dim = 0; dim < dim_num; dim++ )
    {
      indx[dim+p*dim_num] = a[dim] - ( order_1d[dim] - 1 ) / 2;
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

void multigrid_scale_open ( int dim_num, int order_nd, int level_max, 
  int level_1d[], int grid_index[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_SCALE_OPEN renumbers a grid as a subgrid on a higher level.
//
//  Discussion:
//
//    This routine takes a grid associated with a given value of
//    LEVEL, and multiplies all the indices by a power of 2, so that
//    the indices reflect the position of the same points, but in
//    a grid of level LEVEL_MAX.
//
//    For an open grid, going from one level to the next, a set of indices
//    will be rescaled by 2*INDEX-1.
//
//  Modified:
//
//    08 June 2007
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

  for ( dim = 0; dim < dim_num; dim++ )
  {
    factor = i4_power ( 2, level_max - level_1d[dim] );

    for ( order = 0; order < order_nd; order++ )
    {
      grid_index[dim+order*dim_num] = grid_index[dim+order*dim_num] * factor;
    }
  }

  return;
}
//****************************************************************************80

double *product_weights ( int dim_num, int order_1d[], int order_nd, int rule )

//****************************************************************************80
//
//  Purpose:
//
//    PRODUCT_WEIGHTS computes the weights of a product rule.
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
//    09 November 2007
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
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
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
    if ( rule == 1 )
    {
      w_1d = cc_weights ( order_1d[dim] );
    }
    else if ( rule == 2 )
    {
      w_1d = f1_weights ( order_1d[dim] );
    }
    else if ( rule == 3 )
    {
      w_1d = f2_weights ( order_1d[dim] );
    }
    else if ( rule == 4 )
    {
      w_1d = gp_weights ( order_1d[dim] );
    }
    else if ( rule == 5 )
    {
      w_1d = gl_weights ( order_1d[dim] );
    }
    else if ( rule == 6 )
    {
      w_1d = gh_weights ( order_1d[dim] );
    }
    else if ( rule == 7 )
    {
      w_1d = lg_weights ( order_1d[dim] );
    }
    else
    {
      cout << "\n";
      cout << "PRODUCT_WEIGHTS - Fatal error!\n";
      cout << "  Unrecognized rule number = " << rule << "\n";
      exit ( 1 );
    }

    r8vec_direct_product2 ( dim, order_1d[dim], w_1d, dim_num, 
      order_nd, w_nd );

    delete [] w_1d; 
  }

  return w_nd;
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
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
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
//    24 March 2008
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
//    Output, double R8_CHOOSE, the number of combinations of N
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
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = N! = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
//****************************************************************************80

double r8_factorial2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL2 computes the double factorial function.
//
//  Discussion:
//
//    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//  Example:
//
//     N    N!!
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the double factorial 
//    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
//
//    Output, double R8_FACTORIAL2, the value of N!!.
//
{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
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

double r8_mop ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOP returns the I-th power of -1 as an R8 value.
//
//  Discussion:
//
//    An R8 is an double value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, double R8_MOP, the I-th power of -1.
//
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
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
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
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
//***************************************************************************80 

void sparse_grid ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] )
  
//***************************************************************************80 
// 
//  Purpose:
//
//    SPARSE_GRID computes a sparse grid. 
// 
//  Discussion: 
//
//    A Smolyak construction is used to create a multidimensional sparse grid.  
// 
//    The user specifies: 
//    * the spatial dimension of the quadrature region, 
//    * the level that defines the Smolyak grid. 
//    * the 1D quadrature rule. 
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
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
//    Input, int LEVEL_MAX, controls the size of the final 
//    sparse grid. 
// 
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Input, int POINT_NUM, the number of points in the grid, 
//    as determined by LEVELS_INDEX_SIZE.
//
//    Output, double GRID_WEIGHT[POINT_NUM], the weights. 
// 
//    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points. 
// 
{
  if ( rule == 1 )
  {
    sparse_grid_cfn ( dim_num, level_max, rule, point_num, grid_weight, 
      grid_point );
  }
  else if ( 2 <= rule && rule <= 4 )
  {
    sparse_grid_ofn ( dim_num, level_max, rule, point_num, grid_weight, 
      grid_point );
  }
  else if ( 5 <= rule && rule <= 6 )
  {
    sparse_grid_own ( dim_num, level_max, rule, point_num, grid_weight, 
      grid_point );
  }
  else if ( 7 == rule )
  {
    sparse_grid_onn ( dim_num, level_max, rule, point_num, grid_weight, 
      grid_point );
  }
  else
  {
    cout << "\n";
    cout << "SPARSE_GRID - Fatal error!\n";
    cout << "  Illegal input rule index = " << rule << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int sparse_grid_cc_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CC_SIZE sizes a sparse grid using Clenshaw Curtis rules.
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
  int j;
  int l;
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
  for ( l = 2; l <= level_max; l++ )
  {
    j = j * 2;
    new_1d[l] = j;
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

void sparse_grid_cfn ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_CFN computes a sparse grid based on a CFN 1D rule. 
//
//  Discussion:
//
//    The 1D quadrature rule is assumed to be Closed Fully Nested.
//
//    Closed Fully Nested rules include Clenshaw Curtis rules.
//
//    A Smolyak construction is used to create a multidimensional sparse grid.  
// 
//    The user specifies: 
//    * the spatial dimension of the quadrature region, 
//    * the level that defines the Smolyak grid. 
//    * the quadrature rule.
//    * the number of points. 
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
//    Input, int LEVEL_MAX, controls the size of the final sparse grid.
//
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Input, int POINT_NUM, the number of points in the grid, as determined
//    by SPARSE_GRID_SIZE_CFN.
//
//    Output, double GRID_WEIGHTS[POINT_NUM], the weights.
//
//    Output, double GRID_POINTS[DIM_NUM*POINT_NUM], the points.
//
{
  int dim;
  int *grid_base;
  int *grid_index;
  int order_max;
  int point;

  if ( rule != 1 )
  {
    cout << "\n";
    cout << "SPARSE_GRID_CFN - Fatal error!\n";
    cout << "  Illegal input rule index = " << rule << "\n";
    exit ( 1 );
  }
//
//  Determine the index vector, relative to the full product grid,
//  that identifies the points in the sparse grid.
//
  grid_index = new int[dim_num*point_num];
  grid_base = new int[dim_num*point_num];

  levels_index_cfn ( dim_num, level_max, point_num, grid_index, grid_base );
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
      if ( rule == 1 )
      {
        grid_point[dim+point*dim_num] = 
          cc_abscissa ( order_max, grid_index[dim+point*dim_num] + 1 );
      }
    }
  }
//
//  Gather the weights.
//
  sparse_grid_weights_cfn ( dim_num, level_max, rule, point_num, grid_index, 
    grid_weight );

  delete [] grid_base;
  delete [] grid_index;

  return;
}
//***************************************************************************80 

void sparse_grid_ofn ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] )
  
//***************************************************************************80 
// 
//  Purpose:
//
//    SPARSE_GRID_OFN computes a sparse grid based on an OFN 1D rule. 
// 
//  Discussion: 
//
//    The 1D quadrature rule is assumed to be Open Fully Nested.
//
//    Open Fully Nested rules include Fejer 1, Fejer 2, and Gauss Patterson rules.
//
//    A Smolyak construction is used to create a multidimensional sparse grid.  
// 
//    The user specifies: 
//    * the spatial dimension of the quadrature region, 
//    * the level that defines the Smolyak grid. 
//    * the 1D quadrature rule. 
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2008
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
//    Input, int LEVEL_MAX, controls the size of the final 
//    sparse grid. 
// 
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Input, int POINT_NUM, the number of points in the grid, 
//    as determined by LEVELS_INDEX_SIZE.
//
//    Output, double GRID_WEIGHT[POINT_NUM], the weights. 
// 
//    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points. 
// 
{
  int dim;
  int *grid_base;
  int *grid_index;
  int order_max;
  int point;

  if ( rule < 2 || 4 < rule )
  {
    cout << "\n";
    cout << "SPARSE_GRID_OFN - Fatal error!\n";
    cout << "  Illegal input rule index = " << rule << "\n";
    exit ( 1 );
  }
// 
//  Determine the index vector, relative to the full product grid, 
//  that identifies the points in the sparse grid. 
//
  grid_base = new int[dim_num*point_num];
  grid_index = new int[dim_num*point_num];

  levels_index_ofn ( dim_num, level_max, point_num, grid_index, grid_base );
// 
//  Compute the physical coordinates of the abscissas. 
//
  order_max = i4_power ( 2, level_max + 1 ) - 1;

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      if ( rule == 2 )
      {
        grid_point[dim+point*dim_num] = 
          f1_abscissa ( order_max, grid_index[dim+point*dim_num] );
      }
      else if ( rule == 3 )
      {
        grid_point[dim+point*dim_num] =
          f2_abscissa ( order_max, grid_index[dim+point*dim_num] );
      }
      else if ( rule == 4 )
      {
        grid_point[dim+point*dim_num] = 
          gp_abscissa ( order_max, grid_index[dim+point*dim_num] );
      }
    }
  }
// 
//  Gather the weights. 
//
  sparse_grid_weights_ofn ( dim_num, level_max, rule, point_num, 
    grid_index, grid_weight );

  delete [] grid_base;
  delete [] grid_index;

  return;
}
//****************************************************************************80

int sparse_grid_ofn_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_OFN_SIZE sizes a sparse grid using Open Fully Nested rules.
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
//    23 December 2009
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
  int j;
  int l;
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
  for ( l = 1; l <= level_max; l++ )
  {
    new_1d[l] = 2 * new_1d[l-1];
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

void sparse_grid_onn ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_ONN computes a sparse grid based on a ONN 1D rule.
//
//  Discussion:
//
//    The 1D quadrature rule is assumed to be Open Non-Nested.
//    Such rules include Gauss Laguerre rules.
//
//    A Smolyak construction is used to create a multidimensional sparse grid.  
//
//    The user specifies: 
//    * the spatial dimension of the quadrature region, 
//    * the level that defines the Smolyak grid. 
//    * the quadrature rule;
//    * the number of points in the rule.
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
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Input, int POINT_NUM, the number of points in the grid, as determined
//    by SPARSE_GRID_SIZE_ONN.
//
//    Output, double GRID_WEIGHT[POINT_NUM], the weights.
//
//    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points.
//
{
  double coeff;
  int dim;
  int *grid_base2;
  int *grid_index2;
  double *grid_weight2;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int *order_1d;
  int order_nd;
  int order_max;
  int point;
  int point_num2;
  int t;

  for ( point = 0; point < point_num; point++ )
  {
    grid_weight[point] = 0.0;
  }
//
//  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
//
  point_num2 = 0;

  level_min = i4_max ( 0, level_max + 1 - dim_num );
  
  grid_base2 = new int[dim_num];
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

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
//  The relationship is the same as for other OPEN rules.
//
      level_to_order_open ( dim_num, level_1d, order_1d );
      
      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_base2[dim] = order_1d[dim];
      }
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  Compute the weights for this product grid.
//
      grid_weight2 = product_weights ( dim_num, order_1d, order_nd, rule );
//
//  Now determine the coefficient of the weight.
//
      coeff = r8_mop ( level_max - level ) 
        * r8_choose ( dim_num - 1, level_max - level );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
//
      grid_index2 = multigrid_index_onn ( dim_num, order_1d, order_nd );

      for ( point = 0; point < order_nd; point++ )
      {

        if ( point_num <= point_num2 )
        {
          cout << "\n";
          cout << "SPARSE_GRID_ONN - Fatal error!\n";
          cout << "  Exceeding maximum point index POINT_NUM = " 
               << point_num << "\n";
          exit ( 1 );
        }

        lg_abscissa ( dim_num, 1, grid_index2+point*dim_num, 
          grid_base2, grid_point+point_num2*dim_num );

        grid_weight[point_num2] = coeff * grid_weight2[point];
	  
        point_num2 = point_num2 + 1;
      }

      delete [] grid_index2;
      delete [] grid_weight2;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] grid_base2;
  delete [] level_1d;
  delete [] order_1d;
  
  if ( point_num2 < point_num )
  {
    cout << "\n";
    cout << "SPARSE_GRID_ONN - Fatal error!\n";
    cout << "  Set fewer points than POINT_NUM = " << point_num << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void sparse_grid_own ( int dim_num, int level_max, int rule, int point_num, 
  double grid_weight[], double grid_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_OWN computes a sparse grid based on an OWN 1D rule.
//
//  Discussion:
//
//    The 1D quadrature rule is assumed to be Open Weakly Nested.
//    Such rules include Gauss Hermite and Gauss Legendre rules.
//
//    A Smolyak construction is used to create a multidimensional sparse grid.  
// 
//    The user specifies: 
//    * the spatial dimension of the quadrature region, 
//    * the level that defines the Smolyak grid,
//    * the rule;
//    * the number of points.
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
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Input, int POINT_NUM, the number of points in the grid, as determined
//    by LEVELS_INDEX_SIZE_OWN.
//
//    Output, double GRID_WEIGHT[POINT_NUM], the weights.
//
//    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points.
//
{
  double coeff;
  int dim;
  int *grid_base2;
  int *grid_index2;
  int *grid_level;
  double *grid_point_temp;
  double *grid_weight2;
  int h;
  int level;
  int *level_1d;
  int level_min;
  int level_min2;
  bool more;
  int *order_1d;
  int order_nd;
  int order_max;
  int point;
  int point_num2;
  int point2;
  int point3;
  int t;

  for ( point = 0; point < point_num; point++ )
  {
    grid_weight[point] = 0.0;
  }
//
//  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
//
  point_num2 = 0;

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  if ( dim_num == 1 )
  {
    level_min2 = level_min;
  }
  else
  {
    level_min2 = 0;
  }
  
  grid_base2 = new int[dim_num];
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( level = level_min2; level <= level_max; level++ )
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
//  The relationship is the same as for other OPEN rules.
//  The GL rule differs from the other OPEN rules only in the nesting behavior.
//
      level_to_order_open ( dim_num, level_1d, order_1d );
      
      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_base2[dim] = ( order_1d[dim] - 1 ) / 2;
      }
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  Compute the weights for this product grid.
//
      grid_weight2 = product_weights ( dim_num, order_1d, order_nd, rule );
//
//  Now determine the coefficient of the weight.
//
      coeff = r8_mop ( level_max - level ) 
        * r8_choose ( dim_num - 1, level_max - level );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
//
      grid_index2 = multigrid_index_own ( dim_num, order_1d, order_nd );
//
//  Determine the first level of appearance of each of the points.
//  This allows us to flag certain points as being repeats of points
//  generated on a grid of lower level.  
//
//  This is SLIGHTLY tricky.
//
      grid_level = index_level_own ( level, level_max, dim_num, order_nd, 
        grid_index2, grid_base2 );
//
//  Only keep those points which first appear on this level.
//
      for ( point = 0; point < order_nd; point++ )
      {
//
//  Either a "new" point (increase count, create point, create weight)
//
        if ( grid_level[point] == level )
        {

          if ( point_num <= point_num2 )
          {
            cout << "\n";
            cout << "SPARSE_GRID_OWN - Fatal error!\n";
            cout << "  Exceeding maximum point index POINT_NUM = "
                 << point_num << "\n";
            exit ( 1 );
          }

          if ( rule == 5 )
          {
            gl_abscissa ( dim_num, 1, grid_index2+point*dim_num, 
              grid_base2, grid_point+point_num2*dim_num );
          }
          else if ( rule == 6 )
          {
            gh_abscissa ( dim_num, 1, grid_index2+point*dim_num, 
              grid_base2, grid_point+point_num2*dim_num );
          }
          else
          {
            cout << "\n";
            cout << "SPARSE_GRID_OWN - Fatal error!\n";
            cout << "  Unrecognized rule number = " << rule << "\n";
            exit ( 1 );
          }

          if ( level_min <= level )
          {
            grid_weight[point_num2] = coeff * grid_weight2[point];
          }
	  
          point_num2 = point_num2 + 1;
        }
//
//  or an already existing point (create point temporarily, find match,
//  add weight to matched point's weight).
//
        else
        {
          if ( level_min <= level )
          {
            grid_point_temp = new double[dim_num];

            if ( rule == 5 )
            {
              gl_abscissa ( dim_num, 1, grid_index2+point*dim_num, 
                grid_base2, grid_point_temp );
            }
            else if ( rule == 6 )
            {
              gh_abscissa ( dim_num, 1, grid_index2+point*dim_num, 
                grid_base2, grid_point_temp );
            }
            else
            {
              cout << "\n";
              cout << "SPARSE_GRID_OWN - Fatal error!\n";
              cout << "  Unrecognized rule number = " << rule << "\n";
              exit ( 1 );
            }

            for ( point2 = 0; point2 < point_num2; point2++ )
            {
              point3 = point2;
              for ( dim = 0; dim < dim_num; dim++ )
              {
                if ( grid_point[dim+point2*dim_num] != grid_point_temp[dim] )
                {
                  point3 = -1;
                  break;
                }
              }
              if ( point3 == point2 )
              {
                break;
              }
            }
          
            if ( point3 == -1 )
            {
              cout << "\n";
              cout << "SPARSE_GRID_OWN - Fatal error!\n";
              cout << "  Could not match point.\n";
              exit ( 1 );
            }
          
            grid_weight[point3] = grid_weight[point3] + 
              ( double ) ( coeff ) * grid_weight2[point];       
          } 
        }
      }

      delete [] grid_index2;
      delete [] grid_level;
      delete [] grid_weight2;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] grid_base2;
  delete [] level_1d;
  delete [] order_1d;

  if ( point_num2 < point_num )
  {
    cout << "\n";
    cout << "SPARSE_GRID_OWN - Fatal error!\n";
    cout << "  Set fewer points than POINT_NUM = " << point_num << "\n";
    exit ( 1 );
  }
  
  return;
}
//****************************************************************************80

void sparse_grid_weights_cfn ( int dim_num, int level_max, int rule, 
  int point_num, int grid_index[], double grid_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_WEIGHTS_CFN computes sparse grid weights based on a CFN 1D rule.
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
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
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
  double coeff;
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
      grid_weight[point] = i4_power ( 2, dim_num );
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
      grid_index2 = multigrid_index_cfn ( dim_num, order_1d, order_nd );
//
//  Compute the weights for this grid.
//
      grid_weight2 = product_weights ( dim_num, order_1d, order_nd, rule );
//
//  Adjust the grid indices to reflect LEVEL_MAX.
//
      multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, 
        grid_index2 );
//
//  Now determine the coefficient.
//
      coeff = r8_mop ( level_max - level ) 
        * r8_choose ( dim_num - 1, level_max - level );

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
              + coeff * grid_weight2[point2];
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

void sparse_grid_weights_ofn ( int dim_num, int level_max, int rule, 
  int point_num, int grid_index[], double grid_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_WEIGHTS_OFN computes sparse grid weights based on a OFN 1D rule.
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
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
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
  double coeff;
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
      grid_weight[point] = i4_power ( 2, dim_num );
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
      level_to_order_open ( dim_num, level_1d, order_1d );
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  Generate the indices of the points corresponding to the grid.
//
      grid_index2 = multigrid_index_ofn ( dim_num, order_1d, order_nd );
//
//  Compute the weights for this grid.
//
      grid_weight2 = product_weights ( dim_num, order_1d, order_nd, rule );
//
//  Adjust the grid indices to reflect LEVEL_MAX.
//
      multigrid_scale_open ( dim_num, order_nd, level_max, level_1d, 
        grid_index2 );
//
//  Now determine the coefficient.
//
      coeff = r8_mop ( level_max - level ) 
        * r8_choose ( dim_num - 1, level_max - level );

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
              + coeff * grid_weight2[point2];
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
//  Example:
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
//    25 May 2007
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
