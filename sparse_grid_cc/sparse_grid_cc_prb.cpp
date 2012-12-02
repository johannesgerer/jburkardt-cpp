# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <string>

using namespace std;

# include "sparse_grid_cc.H"

using namespace webbur;

int main ( );

void test01 ( int dim_min, int dim_max, int level_max_min, int level_max_max );
void test015 ( int dim_min, int dim_max, int level_max_min, int level_max_max );
void test02 ( int dim_num, int level_max );
void test03 ( int dim_num, int level_max );
void test04 ( int dim_num, int level_max );
void test05 ( int dim_num, int level_max, int degree_max );
void test06 ( int dim_num, int level_max );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_CC_PRB.
//
//  Discussion:
//
//    SPARSE_GRID_CC_PRB tests the SPARSE_GRID_CC routines.
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
{
  int dim_max;
  int dim_min;
  int dim_num;
  int level_max;
  int level_max_max;
  int level_max_min;

  webbur::timestamp ( );
  cout << "\n";
  cout << "SPARSE_GRID_CC_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPARSE_GRID_CC library.\n";
//
//  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
//
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 10;
  test01 ( dim_min, dim_max, level_max_min, level_max_max );

  cout << "\n";
  webbur::timestamp ( );

  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 10;
  test01 ( dim_min, dim_max, level_max_min, level_max_max );

  cout << "\n";
  webbur::timestamp ( );

  dim_min = 100;
  dim_max = 100;
  level_max_min = 0;
  level_max_max = 4;
  test01 ( dim_min, dim_max, level_max_min, level_max_max );

  cout << "\n";
  webbur::timestamp ( );
//
//  Count number of points in sparse rule from DIM_MIN to DIM_MAX, LEVEL_MAX_MAX.
//
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 10;
  test015 ( dim_min, dim_max, level_max_min, level_max_max );

  cout << "\n";
  webbur::timestamp ( );

  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 10;
  test015 ( dim_min, dim_max, level_max_min, level_max_max );

  cout << "\n";
  webbur::timestamp ( );

  dim_min = 100;
  dim_max = 100;
  level_max_min = 0;
  level_max_max = 4;
  test015 ( dim_min, dim_max, level_max_min, level_max_max );

  cout << "\n";
  webbur::timestamp ( );
//
//  Compute abstract grid indices of sparse grid points as selected from product grid
//  for DIMENSION, LEVEL_MAX.
//
  test02 ( 2, 3 );
  test02 ( 2, 4 );
  test02 ( 3, 0 );
  test02 ( 3, 2 );
  test02 ( 6, 2 );
//
//  Compute sparse Clenshaw-Curtis rule for DIMENSION, LEVEL_MAX.
//
  test03 ( 2, 3 );
  test03 ( 3, 0 );
  test03 ( 3, 1 );
//
//  Test sum of weights for DIMENSION, LEVEL_MAX.
//
  test04 ( 2, 4 );
  test04 ( 3, 0 );
  test04 ( 3, 1 );
  test04 ( 3, 6 );
  test04 ( 10, 3 );
//
//  Test monomial exactness for DIMENSION, LEVEL_MAX, DEGREE_MAX.
// 
  test05 ( 2, 0, 3 );
  test05 ( 2, 1, 5 );
  test05 ( 2, 2, 7 );
  test05 ( 2, 3, 9 );
  test05 ( 2, 4, 11 );
  test05 ( 2, 5, 13 );
  
  test05 ( 3, 0, 2 );
  test05 ( 3, 1, 4 );
  test05 ( 3, 2, 6 );
  test05 ( 3, 3, 8 );
//
//  Show how to write a rule to a file.
//
  dim_num = 2;
  level_max = 3;

  test06 ( dim_num, level_max );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPARSE_GRID_CC_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );
  
  return 0;
}
//****************************************************************************80

void test01 ( int dim_min, int dim_max, int level_max_min, int level_max_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests SPARSE_GRID_CFN_SIZE.
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
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
{
  int dim_num;
  int level_max;
  int point_num;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  SPARSE_GRID_CFN_SIZE returns the number of distinct\n";
  cout << "  points in a sparse grid of Closed Fully Nested rules.\n";
  cout << "\n";
  cout << "  Each sparse grid is of spatial dimension DIM,\n";
  cout << "  and is made up of all product grids of levels up to LEVEL_MAX.\n";
  cout << "\n";
  cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    cout << "  " << setw(8) << dim_num;
  }
  cout << "\n";
  cout << "\n";
  cout << "   LEVEL_MAX\n";
  cout << "\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    cout << "    " << setw(4) << level_max;
    for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
    {
      point_num = webbur::sparse_grid_cfn_size ( dim_num, level_max );
      cout << "  " << setw(8) << point_num;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test015 ( int dim_min, int dim_max, int level_max_min, int level_max_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST015 tests SPARSE_GRID_CCS_SIZE.
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
//  Parameters:
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
{
  int dim_num;
  int level_max;
  int point_num;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  SPARSE_GRID_CCS_SIZE returns the number of distinct\n";
  cout << "  points in a Clenshaw Curtis Slow-Growth sparse grid.\n";
  cout << "\n";
  cout << "  Each sparse grid is of spatial dimension DIM,\n";
  cout << "  and is made up of all product grids of levels up to LEVEL_MAX.\n";
  cout << "\n";
  cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    cout << "  " << setw(8) << dim_num;
  }
  cout << "\n";
  cout << "\n";
  cout << "   LEVEL_MAX\n";
  cout << "\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    cout << "    " << setw(4) << level_max;
    for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
    {
      point_num = webbur::sparse_grid_ccs_size ( dim_num, level_max );
      cout << "  " << setw(8) << point_num;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests SPARSE_GRID_CC_INDEX.
//
//  Discussion:
//
//    The routine computes the indices of the unique points used in a sparse 
//    multidimensional grid whose size is controlled by a parameter LEVEL_MAX.
//
//    Once these indices are returned, they can be converted into
//    Clenshaw Curtis points.
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
//    Input, int LEVEL_MAX, the level.
//
{
  int dim;
  int *grid_index;
  int point;
  int point_num;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  SPARSE_GRID_CC_INDEX returns all grid indexes\n";
  cout << "  whose level value satisfies\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "  Here, LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  and the order of the rule is 2**LEVEL + 1.\n";

  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  point_num = webbur::sparse_grid_cfn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  Number of unique points in the grid = " << point_num << "\n";
//
//  Compute the orders and points.
//
  grid_index = webbur::sparse_grid_cc_index ( dim_num, level_max, point_num );
//
//  Now we're done.  Print the merged grid data.
//
  cout << "\n";
  cout << "  Grid index:\n";
  cout << "\n";
  for ( point = 0; point < point_num; point++ )
  {
    cout << "  " << setw(4) << point << "  ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(6) << grid_index[dim+point*dim_num];
    }
    cout << "\n";
  }

  delete [] grid_index;

  return;
}
//****************************************************************************80

void test03 ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 call SPARSE_GRID_CC to create a Clenshaw Curtis grid.
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
//    Input, int LEVEL_MAX, the level.
//
{
  int dim;
  double *grid_point;
  double *grid_weight;
  int point;
  int point_num;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  SPARSE_GRID_CC makes a sparse Clenshaw Curtis grid.\n";
  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
//
//  Determine the number of points.
//
  point_num = webbur::sparse_grid_cfn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  Number of unique points in the grid = " << point_num << "\n";
//
//  Allocate space for the weights and points.
//
  grid_weight = new double[point_num];
  grid_point = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  webbur::sparse_grid_cc ( dim_num, level_max, point_num, grid_weight, grid_point );
//
//  Print them out.
//
  cout << "\n";
  cout << "  Grid weights:\n";
  cout << "\n";
  for ( point = 0; point < point_num; point++ )
  {
    cout << "  " << setw(4) << point
         << fixed << setw(10) << grid_weight[point] << "\n";
  }
  
  cout << "\n";
  cout << "  Grid points:\n";
  cout << "\n";
  for ( point = 0; point < point_num; point++ )
  {
    cout << "  " << setw(4) << point;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << fixed << setw(10) << grid_point[dim+point*dim_num];
    }
    cout << "\n";
  }

  delete [] grid_point;
  delete [] grid_weight;

  return;
}
//****************************************************************************80

void test04 ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 sums the weights and compares them to 2^DIM_NUM.
//
//  Discussion:
//
//    This routine gets the sparse grid indices and determines the 
//    corresponding sparse grid abscissas.
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
//    Input, int LEVEL_MAX, the level.
//
{
  double *grid_point;
  double *grid_weight;
  int point;
  int point_num;
  double weight_sum;
  double weight_sum_error;
  double weight_sum_exact;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  Compute the weights of a Clenshaw Curtis sparse grid .\n";
  cout << "\n";
  cout << "  As a simple test, sum these weights.\n";
  cout << "  They should sum to exactly 2^DIM_NUM.\n";
  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
//
//  Determine the number of points.
//
  point_num = webbur::sparse_grid_cfn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  Number of unique points in the grid = " << point_num << "\n";
//
//  Allocate space for the weights and points.
//
  grid_weight = new double[point_num];
  grid_point = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  webbur::sparse_grid_cc ( dim_num, level_max, point_num, grid_weight, grid_point );
//
//  Sum the weights.
//
  weight_sum = 0.0;
  for ( point = 0; point < point_num; point++ )
  {
    weight_sum = weight_sum + grid_weight[point];  }
  
  weight_sum_exact = pow ( 2.0, dim_num );
  
  weight_sum_error = fabs ( weight_sum - weight_sum_exact );
  
  cout << "\n";
  cout << "    Weight sum     Exact sum    Difference\n";
  cout << "\n";
  cout << "  " << setw(12) << weight_sum
       << "  " << setw(12) << weight_sum_exact
       << "  " << setw(12) << weight_sum_error << "\n";

  delete [] grid_point;
  delete [] grid_weight;

  return;
}
//****************************************************************************80

void test05 ( int dim_num, int level_max, int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests a Clenshaw Curtis sparse grid rule for monomial exactness.
//
//  Discussion:
//
//    This test is going to check EVERY monomial of total degree DEGREE_MAX
//    or less.  Even for a moderately high dimension of DIM_NUM = 10, you
//    do NOT want to use a large value of DEGREE_MAX, since there are
//
//      1         monomials of total degree 0,
//      DIM_NUM   monomials of total degree 1,
//      DIM_NUM^2 monomials of total degree 2,
//      DIM_NUM^3 monomials of total degree 3, and so on.
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
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the level.
//
//    Input, int DEGREE_MAX, the maximum monomial total degree to check.
//
{
  int degree;
  int dim;
  bool error;
  int *expon;
  double *grid_point;
  double *grid_weight;
  int h;
  int last;
  bool more;
  int point;
  int point_num;
  double quad_error;
  int t;
  double volume;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Check the exactness of a Clenshaw Curtis sparse grid quadrature rule,\n";
  cout << "  applied to all monomials of orders 0 to DEGREE_MAX.\n";
  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "\n";
  cout << "  The maximum total degree to be checked is DEGREE_MAX = " << degree_max << "\n";
  cout << "\n";
  cout << "  We expect this rule to be accurate up to and including total degree " 
	   << 2 * level_max + 1 << "\n";
//
//  Determine the number of points in the rule.
//
  point_num = webbur::sparse_grid_cfn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  Number of unique points in the grid = " << point_num << "\n";
//
//  Allocate space for the weights and points.
//
  grid_weight = new double[point_num];
  grid_point = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  webbur::sparse_grid_cc ( dim_num, level_max, point_num, grid_weight, 
  grid_point );
//
//  Rescale the weights, and translate the abscissas.
//
  volume = pow ( 2.0, dim_num );

  for ( point = 0; point < point_num; point++ )
  {
    grid_weight[point] = grid_weight[point] / volume;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    for ( point = 0; point < point_num; point++ )
    {
      grid_point[dim+point*dim_num] = ( grid_point[dim+point*dim_num] + 1.0 ) 
      / 2.0; 
    }
  }
//
//  Explore the monomials.
//
  expon = new int[dim_num];

  cout << "\n";
  cout << "      Error      Total   Monomial\n";
  cout << "                 Degree  Exponents\n";
  cout << "\n";

  for ( degree = 0; degree <= degree_max; degree++ )
  {
    more = false;

    for ( ; ; )
    {
      webbur::comp_next ( degree, dim_num, expon, &more, &h, &t );

      quad_error = webbur::monomial_quadrature ( dim_num, expon, point_num, 
        grid_weight, grid_point );

      cout << "  " << setw(12) << quad_error
           << "     " << setw(2)  << degree
           << "      ";

      for ( dim = 0; dim < dim_num; dim++ )
      {
        cout << setw(2) << expon[dim];
      }
      cout << "\n";

      if ( !more )
      {
        break;
      }
    }
    cout << "\n";
  }

  delete [] expon;
  delete [] grid_point;
  delete [] grid_weight;

  return;
}
//***************************************************************************80

void test06 ( int dim_num, int level_max )

//***************************************************************************80
//
//  Purpose:
//
//    TEST06 creates a sparse Clenshaw-Curtis grid and writes it to a file.
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
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the level.
//
{
  int dim;
  int point;
  int point_num;
  double *r;
  string r_filename;
  double *w;
  string w_filename;
  double *x;
  string x_filename;

  cout << "\n";
  cout << "TEST06:\n";
  cout << "  Call SPARSE_GRID_CC to make a sparse Clenshaw-Curtis grid.\n";
  cout << "  Write the data to a set of quadrature files.\n";
    
  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
//
//  Determine the number of points.
//
  point_num = webbur::sparse_grid_cfn_size ( dim_num, level_max );
//
//  Allocate space for the weights and points.
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

  webbur::sparse_grid_cc ( dim_num, level_max, point_num, w, x );
//
//  Write the data out.
//
  r_filename = "cc_d" + i4_to_string ( dim_num, "%d" ) + "_level" 
    + i4_to_string ( level_max, "%d" ) + "_r.txt";
  w_filename = "cc_d" + i4_to_string ( dim_num, "%d" ) + "_level" 
    + i4_to_string ( level_max, "%d" ) + "_w.txt";
  x_filename = "cc_d" + i4_to_string ( dim_num, "%d" ) + "_level" 
    + i4_to_string ( level_max, "%d" ) + "_x.txt";

  webbur::r8mat_write ( r_filename, dim_num, 2,         r );
  webbur::r8mat_write ( w_filename, 1,       point_num, w );
  webbur::r8mat_write ( x_filename, dim_num, point_num, x );

  cout << "\n";
  cout << "  R data written to \"" << r_filename << "\".\n";
  cout << "  W data written to \"" << w_filename << "\".\n";
  cout << "  X data written to \"" << x_filename << "\".\n";

  delete [] r;
  delete [] w;
  delete [] x;

  return;
}

