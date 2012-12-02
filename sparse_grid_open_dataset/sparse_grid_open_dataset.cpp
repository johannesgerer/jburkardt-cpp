# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <string>

using namespace std;

int main ( int argc, char *argv[] );
int *abscissa_level_open_nd ( int level_max, int dim_num, int test_num, 
  int test_val[] );
void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
double f2_abscissa ( int order, int i );
double *f2_weights ( int order );
double gp_abscissa ( int order, int i );
double *gp_weights ( int order );
int i4_choose ( int n, int k );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
string i4_to_string ( int i4, string format );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
int i4vec_product ( int n, int a[] );
int index_to_level_open ( int dim_num, int t[], int order, int level_max );
void level_to_order_open ( int dim_num, int level[], int order[] );
int *multigrid_index1 ( int dim_num, int order_1d[], int order_nd );
void multigrid_scale_open ( int dim_num, int order_nd, int level_max, 
  int level_1d[], int grid_index[] );
double nco_abscissa ( int order, int i );
double *nco_weights ( int order );
double *product_weights_open ( int dim_num, int order_1d[], int order_nd, 
  int rule );
double r8_epsilon ( );
double r8_huge ( );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title );
double r8vec_sum ( int n, double a[] );
int s_len_trim ( string s );
int sparse_grid_ofn_size ( int dim_num, int level_max );
int *spgrid_open_index ( int dim_num, int level_max, int point_num );
double *spgrid_open_weights ( int dim_num, int level_max, int point_num, 
  int grid_index[], int rule );
void timestamp ( );
double ts_abscissa ( int order, int i );
double *ts_weights ( int order );
void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_OPEN_DATASET.
//
//  Discussion:
//
//    This program computes a quadrature rule and writes it to a file.
//
//    The quadrature rule is associated with a sparse grid derived from
//    a Smolyak construction using an open 1D quadrature rule. 
//
//    The user specifies:
//    * the spatial dimension of the quadrature region,
//    * the level that defines the Smolyak grid.
//    * the open 1D quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2009
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
  int *grid_index;
  double *grid_point;
  double *grid_region;
  double *grid_weight;
  double h;
  int level_max;
  int m;
  int n;
  int order_max;
  int point;
  int point_num;
  string r_filename;
  int rule;
  string w_filename;
  double weight_sum;
  string x_filename;

  timestamp ( );
  cout << "\n";
  cout << "SPARSE_GRID_OPEN_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Compute the abscissas and weights of a quadrature rule\n";
  cout << "  associated with a sparse grid derived from a Smolyak\n";
  cout << "  construction based on an open quadrature rule.\n";
  cout << "\n";
  cout << "  Inputs to the program include:\n";
  cout << "\n";
  cout << "    DIM_NUM, the spatial dimension.\n";
  cout << "    (typically in the range of 2 to 10)\n";
  cout << "\n";
  cout << "    LEVEL_MAX, the \"level\" of the sparse grid.\n";
  cout << "    (typically in the range of 0, 1, 2, 3, ...\n";
  cout << "\n";
  cout << "    RULE, the 1D quadrature rule\n";
  cout << "    2: Fejer Type 2 (\"F2\").\n";
  cout << "    3: Gauss-Patterson (\"GP\");\n";
  cout << "    4: Newton-Cotes Open (\"NCO\").\n";
  cout << "    5: Tanh-Sinh (\"TS\").\n";
  cout << "\n";
  cout << "  Output from the program includes:\n";
  cout << "\n";
  cout << "    A printed table of the abscissas and weights.\n";
  cout << "\n";
  cout << "    A set of files defining the quadrature rules.\n";
  cout << "\n";
  cout << "    \"***_d?_level?_x.txt\", a file of the abscissas;\n";
  cout << "    \"***_d?_level?_w.txt\", a file of the weights;\n";
  cout << "    \"***_d?_level?_r.txt\", a file of the ranges.\n";
//
//  Get the spatial dimension:
//
  if ( 1 < argc ) 
  {
    dim_num = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "SPARSE_GRID_OPEN_DATASET:\n";
    cout << "  Enter the value of DIM_NUM.\n";

    cin >> dim_num;
  }
  cout << "\n";
  cout << "  Spatial dimension requested is = " << dim_num << "\n";
//
//  Get the product file root name:
//
  if ( 2 < argc ) 
  {
    level_max = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "SPARSE_GRID_OPEN_DATASET:\n";
    cout << "  Enter the value of LEVEL_MAX.\n";
    cin >> level_max;
  }
  cout << "\n";
  cout << "  The sparse grid level is = " << level_max << "\n";
//
//  Get the rule index:
//
  if ( 3 < argc ) 
  {
    rule = atoi ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "SPARSE_GRID_OPEN_DATASET:\n";
    cout << "  Enter the value of RULE.\n";
    cout << "  2 = F2   = Fejer Type 2 Rule,\n";
    cout << "  3 = GP   = Gauss-Patterson,\n";
    cout << "  4 = NCO  = Newton-Cotes Open,\n";
    cout << "  5 = TS   = Tanh-Sinh.\n";
    cin >> rule;
  }
  cout << "\n";
  cout << "  The 1D quadrature rule index = " << rule << "\n";

  if ( rule == 2 )
  {
    cout << "  F2:   Fejer Type 2 Rule.\n";
  }
  else if ( rule == 3 )
  {
    cout << "  GP:   Gauss-Patterson Rule.\n";
  }
  else if ( rule == 4 )
  {
    cout << "  NCO:  Newton-Cotes Open Rule.\n";
  }
  else if ( rule == 5 )
  {
    cout << "  TS:   Tanh-Sinh Rule.\n";
  }
  else
  {
    cout << "\n";
    cout << "SPARSE_GRID_OPEN_DATASET - Fatal error!\n";
    cout << "  Illegal value of RULE.\n";
    exit ( 1 );
  }
//
//  How many distinct points will there be?
//
  point_num = sparse_grid_ofn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  The number of distinct abscissas in the\n";
  cout << "  quadrature rule is determined from the spatial\n";
  cout << "  dimension DIM_NUM and the level LEVEL_MAX.\n";
  cout << "  For the given input, this value will be = " << point_num << "\n";

  grid_point = new double[dim_num*point_num];
//
//  Determine the index vector, relative to the full product grid,
//  that identifies the points in the sparse grid.
//
  grid_index = spgrid_open_index ( dim_num, level_max, point_num );

  i4mat_transpose_print_some ( dim_num, point_num, grid_index, 1, 1,
    dim_num, 10, "  First 10 entries of grid index:" );
//
//  Compute the physical coordinates of the abscissas.
//
  order_max = i4_power ( 2, level_max + 1 ) - 1;

  if ( rule == 5 )
  {
    m = level_max - 3;
    n = ( ( order_max + 1 ) / 2 ) - 1;
    h = 4.0 / ( double ) ( order_max + 1 );

    cout << "  M = " << m 
         << "  ORDER_MAX = " << order_max
         << "  N = " << n
         << "  H = " << h << "\n";
  }

  if ( rule == 2 )
  {
    for ( point = 0; point < point_num; point++ )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_point[dim+point*dim_num] = 
          f2_abscissa ( order_max, grid_index[dim+point*dim_num] );
      }
    }
  }
  else if ( rule == 3 )
  {
    for ( point = 0; point < point_num; point++ )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_point[dim+point*dim_num] = 
          gp_abscissa ( order_max, grid_index[dim+point*dim_num] );
      }
    }
  }
  else if ( rule == 4 )
  {
    for ( point = 0; point < point_num; point++ )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_point[dim+point*dim_num] = 
          nco_abscissa ( order_max, grid_index[dim+point*dim_num] );
      }
    }
  }
  else if ( rule == 5 )
  {
    for ( point = 0; point < point_num; point++ )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_point[dim+point*dim_num] = 
          ts_abscissa ( order_max, grid_index[dim+point*dim_num] );
      }
    }
  }
  r8mat_transpose_print_some ( dim_num, point_num, grid_point, 1, 1,
    dim_num, 10, "  First 10 entries of grid point:" );
//
//  Gather the weights.
//
  grid_weight = spgrid_open_weights ( dim_num, level_max, point_num, 
    grid_index, rule );

  r8vec_print_some ( point_num, grid_weight, 1, 10, 
    "  First 10 grid weights:" );

  weight_sum = r8vec_sum ( point_num, grid_weight );

  cout << "\n";
  cout << "  Weights sum to   " 
       << setprecision(16) << setw(24) << weight_sum << "\n";
  cout << "  Correct value is " 
       << setprecision(16) << setw(24) << pow ( 2.0, dim_num ) << "\n";
//
//  Write the rule to files.
//
  if ( rule == 2 )
  {
    r_filename = "f2_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_r.txt";
    w_filename = "f2_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_w.txt";
    x_filename = "f2_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_x.txt";
  }
  else if ( rule == 3 )
  {
    r_filename = "gp_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_r.txt";
    w_filename = "gp_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_w.txt";
    x_filename = "gp_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_x.txt";
  }
  else if ( rule == 4 )
  {
    r_filename = "nco_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_r.txt";
    w_filename = "nco_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_w.txt";
    x_filename = "nco_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_x.txt";
  }
  else if ( rule == 5 )
  {
    r_filename = "ts_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_r.txt";
    w_filename = "ts_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_w.txt";
    x_filename = "ts_d" + i4_to_string ( dim_num, "%d" )
      + "_level" + i4_to_string ( level_max, "%d" ) + "_x.txt";
  }

  cout << "\n";
  cout << "  Creating X file = \"" << x_filename << "\".\n";

  r8mat_write ( x_filename, dim_num, point_num, grid_point );

  cout << "  Creating W file = \"" << w_filename << "\".\n";

  r8mat_write ( w_filename, 1, point_num, grid_weight );
 
  grid_region = new double[dim_num*2];

  for ( dim = 0; dim < dim_num; dim++ )
  {
    grid_region[dim+0*dim_num] = -1.0;
    grid_region[dim+1*dim_num] = +1.0;
  }
  cout << "  Creating R file = \"" << r_filename << "\".\n";

  r8mat_write ( r_filename, dim_num, 2, grid_region );

  delete [] grid_index;
  delete [] grid_point;
  delete [] grid_region;
  delete [] grid_weight;

  cout << "\n";
  cout << "SPARSE_GRID_OPEN_DATASET:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
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
    exit ( 1 );
  }

  if ( order == 1 )
  {
    value = 0.0;
    return value;
  }

  value = cos ( ( double ) ( order + 1 - i ) * pi 
              / ( double ) ( order + 1 ) );

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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
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
//    22 May 2007
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
//    Input, int N, K, are the values of N and K.
//
//    Output, int I4_CHOOSE, the number of combinations of N
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
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of i1 and i2.
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
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
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
//    I4VEC_PRODUCT multiplies the entries of an I4VEC.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_PRODUCT = 24
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

int index_to_level_open ( int dim_num, int t[], int order, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_TO_LEVEL_OPEN determines the level of a point given its index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
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
//    The arrangement described here works naturally for the 
//    Fejer Type 2, Newton Cotes Open,
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
//    For the Fejer Type 2, Newton Cotes Open, 
//    and Gauss-Patterson rules, the point growth is
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
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

int *multigrid_index1 ( int dim_num, int order_1d[], int order_nd )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_INDEX1 returns an indexed multidimensional grid.
//
//  Discussion:
//
//    For dimension DIM, the second index of INDX may vary from 
//    1 to ORDER_1D[DIM].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
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

double nco_abscissa ( int order, int i )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_ABSCISSA returns the I-th abscissa for the Newton Cotes open rule.
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
//    25 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    1 <= ORDER.
//
//    Input, int I, the index of the desired abscissa.  
//    1 <= I <= ORDER.
//
//    Output, double NCO_ABSCISSA, the value of the I-th 
//    abscissa in the Newton Cotes open rule of order ORDER.
//
{
  double value;
  double x_max = +1.0;
  double x_min = -1.0;

  if ( order < 1 )
  {
    value = - r8_huge ( );
    return value;
  }

  if ( i < 1 || order < i )
  {
    cout << "\n";
    cout << "NCO_ABSCISSA - Fatal error!\n";
    cout << "  1 <= I <= ORDER is required.\n";
    exit ( 1 );
  }

  if ( order == 1 )
  {
    value = ( x_min + x_max ) / 2.0;
    return value;
  }

  value = ( ( double ) ( order - i + 1 ) * x_min
          + ( double ) (         i     ) * x_max )
          / ( double ) ( order     + 1 );

  return value;
}
//****************************************************************************80

double *nco_weights ( int order )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_WEIGHTS computes weights for a Newton-Cotes Open rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order.
//
//    Output, double W[ORDER], the weights.
//
{
  double *diftab;
  int i;
  int j;
  int k;
  double *w;
  double x_max = +1.0;
  double x_min = -1.0;
  double *x;
  double yvala;
  double yvalb;

  diftab = new double[order];
  w = new double[order];
  x = new double[order];

  for ( i = 1; i <= order; i++ )
  {
    x[i-1] = ( ( double ) ( order + 1 - i ) * x_min   
             + ( double ) (             i ) * x_max ) 
           /   ( double ) ( order + 1     );
  }

  for ( i = 1; i <= order; i++ )
  {
//
//  Compute the Lagrange basis polynomial which is 1 at X(I),
//  and zero at the other nodes.
//
    for ( j = 0; j < order; j++ )
    {
      diftab[j] = 0.0;
    }
    diftab[i-1] = 1.0;

    for ( j = 2; j <= order; j++ )
    {
      for ( k = j; k <= order; k++ )
      {
        diftab[order+j-k-1] = ( diftab[order+j-k-1-1] - diftab[order+j-k-1] ) 
          / ( x[order+1-k-1] - x[order+j-k-1] );
      }
    }

    for ( j = 1; j < order; j++ )
    {
      for ( k = 1; k <= order - j; k++ )
      {
        diftab[order-k-1] = diftab[order-k-1] - x[order-k-j] * diftab[order-k];
      }
    }
//
//  Evaluate the antiderivative of the polynomial at the left and
//  right endpoints.
//
    yvala = diftab[order-1] / ( double ) ( order );
    for ( j = order-1; 1 <= j; j-- )
    {
      yvala = yvala * x_min + diftab[j-1] / ( double ) ( j );
    }
    yvala = yvala * x_min;

    yvalb = diftab[order-1] / ( double ) ( order );
    for ( j = order-1; 1 <= j; j-- )
    {
      yvalb = yvalb * x_max + diftab[j-1] / ( double ) ( j );
    }
    yvalb = yvalb * x_max;

    w[i-1] = yvalb - yvala;
  }

  delete [] diftab;
  delete [] x;

  return w;
}
//****************************************************************************80

double *product_weights_open ( int dim_num, int order_1d[], int order_nd, 
  int rule )

//****************************************************************************80
//
//  Purpose:
//
//    PRODUCT_WEIGHTS_OPEN: weights for an open product rule.
//
//  Discussion:
//
//    This routine computes the weights for a quadrature rule which is
//    a product of 1D rules of varying order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2009
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
//    Input, int RULE, the 1D quadrature rule being used.
//    2, Fejer Type 2 Rule;
//    3, Gauss-Patterson Rule,
//    4, Newton-Cotes Open Rule,
//    5, Tanh-Sinh Rule.
//
//    Output, double PRODUCT_WEIGHTS_OPEN[DIM_NUM*ORDER_ND], the product 
//    rule weights.
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
    if ( rule == 2 )
    {
      w_1d = f2_weights ( order_1d[dim] );
    }
    else if ( rule == 3 )
    {
      w_1d = gp_weights ( order_1d[dim] );
    }
    else if ( rule == 4 )
    {
      w_1d = nco_weights ( order_1d[dim] );
    }
    else if ( rule == 5 )
    {
      w_1d = ts_weights ( order_1d[dim] );
    }

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
//    01 July 2004
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
  double value;

  value = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + value )  )
  {
    value = value / 2.0;
  }

  value = 2.0 * value;

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
//    R8MAT_WRITE writes an R8MAT file.
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
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    cout << "  " << setw(8)  << i       
         << "  " << setw(14) << a[i-1]  << "\n";
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
//    An R8VEC is a double precision vector.
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

int *spgrid_open_index ( int dim_num, int level_max, int point_num )

//****************************************************************************80
//
//  Purpose:
//
//    LEVELS_OPEN_INDEX computes open grids with 0 <= LEVEL <= LEVEL_MAX.
//
//  Discussion:
//
//    The necessary dimensions of GRID_INDEX can be
//    determined by calling SPGRID_OPEN_SIZE first.
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
//    Output, int LEVELS_MAX_INDEX[DIM_NUM*POINT_NUM], a list of point indices,
//    representing a subset of the product grid of level LEVEL_MAX,
//    representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
{
  int dim;
  int factor;
  int *grid_index;
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
  grid_index = new int[dim_num*point_num];
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
      grid_index2 = multigrid_index1 ( dim_num, order_1d, order_nd );
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
          for ( dim = 0; dim < dim_num; dim++ )
          {
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

  return grid_index;
}
//****************************************************************************80

double *spgrid_open_weights ( int dim_num, int level_max, int point_num, 
  int grid_index[], int rule )

//****************************************************************************80
//
//  Purpose:
//
//    SPGRID_OPEN_WEIGHTS gathers the weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 July 2008
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
//    Input, int RULE, the 1D quadrature rule being used.
//    2, Fejer Type 2 Rule;
//    3, Gauss-Patterson Rule,
//    4, Newton-Cotes Open Rule,
//    5, Newton-Cotes Open Half Rule.
//
//    Output, double SPGRID_OPEN_WEIGHTS[POINT_NUM], the weights
//    associated with the sparse grid points.
//
{
  bool all_equal;
  int coeff;
  int dim;
  int *grid_index2;
  double *grid_weight;
  double *grid_weight2;
  int h;
  int level;
  int *level_1d;
  int level_min;
  int match;
  bool more;
  int order_nd;
  int *order_1d;
  int point;
  int point2;
  int t;

  grid_weight = new double[point_num];

  if ( level_max == 0 )
  {
    for ( point = 0; point < point_num; point++ )
    {
      grid_weight[point] = i4_power ( 2, dim_num );
    }
    return grid_weight;
  }

  for ( point = 0; point < point_num; point++ )
  {
    grid_weight[point] = 0.0;
  }

  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

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
      grid_index2 = multigrid_index1 ( dim_num, order_1d, order_nd );
//
//  Compute the weights for this grid.
//
      grid_weight2 = product_weights_open ( dim_num, order_1d, order_nd, rule );
//
//  Adjust the grid indices to reflect LEVEL_MAX.
//
      multigrid_scale_open ( dim_num, order_nd, level_max, level_1d, 
        grid_index2 );
//
//  Now determine the coefficient.
//
      coeff = i4_power ( -1, level_max - level ) 
        * i4_choose ( dim_num - 1, level_max - level );

      for ( point2 = 0; point2 < order_nd; point2++ )
      {
        match = -1;

        for ( point = 0; point < point_num; point++ )
        {
          all_equal = true;
          for ( dim = 0; dim < dim_num; dim++ )
          {
            if ( grid_index2[dim+point2*dim_num] 
              != grid_index[dim+point*dim_num] )
            {
              all_equal = false;
              break;
            }
          }
          if ( all_equal )
          {
            grid_weight[point] = grid_weight[point] 
              + ( double ) ( coeff ) * grid_weight2[point2];
            match = point;
            break;
          }
        }

        if ( match == -1 )
        {
          cout << "\n";
          cout << "SPGRID_OPEN_WEIGHTS - Fatal error!\n";
          cout << "  Could not match grid index.\n";
          cout << "  Point index = " << point2 << "\n";
          cout << "\n";
          cout << "  LEVEL = " << level << "\n";
          cout << "\n";
          cout << "  LEVEL_1D:\n";
          for ( dim = 0; dim < dim_num; dim++ )
          {
            cout << setw(6) << level_1d[dim];
          }
          cout << "\n";
          cout << "\n";
          cout << "  ORDER_1D:\n";
          for ( dim = 0; dim < dim_num; dim++ )
          {
            cout << setw(6) << order_1d[dim];
          }
          cout << "\n";
          cout << "\n";
          cout << "  GRID_INDEX2\n";
          for ( dim = 0; dim < dim_num; dim++ )
          {
            cout << setw(6) << grid_index2[dim+point2*dim_num];
          }
          cout << "\n";
          exit ( 1 );
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

  return grid_weight;
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

double ts_abscissa ( int order, int i )

//****************************************************************************80
//
//  Purpose:
//
//    TS_ABSCISSA returns the I-th abscissa for the tanh-sinh rule.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to
//    right.
//
//    This rule is defined on [-1,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//
//    Input, int I, the index of the desired abscissa. 
//    1 <= I <= ORDER.
//
//    Output, double TS_ABSCISSA, the value of the I-th abscissa 
//    in the rule of order ORDER.
//
{
  double ct;
  double ct2;
  double h;
  double pi = 3.141592653589793;
  double st;
  double t;
  double value;

  if ( order < 1 )
  {
    value = - r8_huge ( );
  }
  else if ( i < 1 || order < i )
  {
    value = - r8_huge ( );
  }
  else if ( order == 1 )
  {
    value = 0.0;
  }
  else if ( 2 * i - order - 1 == 0 )
  {
    value = 0.0;
  }
  else
  {
    h = 4.0 / ( double ) ( order + 1 );

    t = ( double ) ( 2 * i - order - 1 ) * h / 2.0;

    ct = cosh ( t );
    st = sinh ( t );
    ct2 = cosh ( 0.5 * pi * st );

    value = tanh ( 0.5 * pi * st );
  }

  return value;
}
//****************************************************************************80

double *ts_weights ( int order )

//****************************************************************************80
//
//  Purpose:
//
//    TS_WEIGHTS computes weights for a tanh-sinh rule.
//
//  Discussion:
//
//    In the 1D case, a sequence of rules is used of increasing order.
//    For low order, the weights do not sum to 2, but with increasing 
//    order, the sum quickly converges to 2.
//
//    However, for sparse grid applications, the lowest order rules are
//    involved in every grid, so it seems it might be useful to force
//    the weights to sum to 2 immediately.  This addresses only one very
//    obvious defect of the lower order rules.  I am not sure what to do
//    about the fact the none of the rules have a definable precision,
//    and the family of rules has not precision but asymptotic accuracy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//
//    Output, double W[ORDER], the weights of the rule.
//
{
  double ct;
  double ct2;
  double h;
  int i;
  double pi = 3.141592653589793;
  double st;
  double t;
  double *w;
  double w_sum;

  if ( order < 1 )
  {
    cerr << "\n";
    cerr << "TS_WEIGHTS - Fatal error!\n";
    cerr << "  ORDER < 1.\n";
    exit ( 1 );
  }

  w = new double[order];

  h = 4.0 / ( double ) ( order + 1 );

  for ( i = 0; i < order; i++ )
  {
    t = ( double ) ( 2 * i - order + 1 ) * h / 2.0;

    ct = cosh ( t );
    st = sinh ( t );
    ct2 = cosh ( 0.5 * pi * st );;

    w[i] = 0.5 * pi * h * ct / ct2 / ct2;
  }
//
//  Normalize the weights so that they sum to 2.0.
//
  w_sum = 0.0;
  for ( i = 0; i < order; i++ )
  {
    w_sum = w_sum + w[i];
  }
  for ( i = 0; i < order; i++ )
  {
    w[i] = 2.0 * w[i] / w_sum;
  }

  return w;
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
