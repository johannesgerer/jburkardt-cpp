# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "sandia_sparse.hpp"

int main ( );
void levels_index_size_test ( int rule, int dim_min, int dim_max, 
  int level_max_min, int level_max_max );
void levels_index_test ( int rule, int dim_num, int level_max );
void sparse_grid_compute_test ( int rule, int dim_num, int level_max );
void sparse_grid_weight_test ( int rule, int dim_num, int level_max );
void sparse_grid_monomial_test ( int rule, int dim_num, int level_max, 
  int degree_max );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SANDIA_SPARSE_PRB.
//
//  Discussion:
//
//    SANDIA_SPARSE_PRB calls a set of problems for SANDIA_SPARSE.
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
{
  int degree_max;
  int dim_max;
  int dim_min;
  int dim_num;
  int level;
  int level_max;
  int level_max_max;
  int level_max_min;
  int rule;
  int rule_max = 7;

  timestamp ( );

  cout << "\n";
  cout << "SANDIA_SPARSE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SANDIA_SPARSE library.\n";
//
//  Test LEVELS_INDEX_SIZE for one example each of CFN, OFN, OWN and ONN rules.
//
  if ( true )
  {
    rule = 1;

    dim_min = 1;
    dim_max = 1;
    level_max_min = 0;
    level_max_max = 10;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 1;
    dim_max = 6;
    level_max_min = 0;
    level_max_max = 6;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 5;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 100;
    dim_max = 100;
    level_max_min = 0;
    level_max_max = 2;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    rule = 2;

    dim_min = 1;
    dim_max = 1;
    level_max_min = 0;
    level_max_max = 10;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 1;
    dim_max = 6;
    level_max_min = 0;
    level_max_max = 6;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 5;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 100;
    dim_max = 100;
    level_max_min = 0;
    level_max_max = 2;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    rule = 5;

    dim_min = 1;
    dim_max = 1;
    level_max_min = 0;
    level_max_max = 10;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 1;
    dim_max = 6;
    level_max_min = 0;
    level_max_max = 6;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 5;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 100;
    dim_max = 100;
    level_max_min = 0;
    level_max_max = 2;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    rule = 7;

    dim_min = 1;
    dim_max = 1;
    level_max_min = 0;
    level_max_max = 10;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 1;
    dim_max = 6;
    level_max_min = 0;
    level_max_max = 6;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 5;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );

    dim_min = 100;
    dim_max = 100;
    level_max_min = 0;
    level_max_max = 2;

    levels_index_size_test ( rule, dim_min, dim_max, level_max_min, 
      level_max_max );
  }
//
//  Test LEVELS_INDEX for one example each of CFN, OFN, OWN and ONN rules.
//
  if ( true )
  {
    rule = 1;

    dim_num = 2;
    level_max = 1;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 2;
    level_max = 3;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 3;
    level_max = 0;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 3;
    level_max = 2;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 6;
    level_max = 2;
    levels_index_test ( rule, dim_num, level_max );

    rule = 2;

    dim_num = 2;
    level_max = 1;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 2;
    level_max = 3;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 3;
    level_max = 0;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 3;
    level_max = 2;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 6;
    level_max = 2;
    levels_index_test ( rule, dim_num, level_max );

    rule = 5;

    dim_num = 2;
    level_max = 1;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 2;
    level_max = 3;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 3;
    level_max = 0;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 3;
    level_max = 2;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 6;
    level_max = 2;
    levels_index_test ( rule, dim_num, level_max );

    rule = 7;

    dim_num = 2;
    level_max = 1;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 2;
    level_max = 3;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 3;
    level_max = 0;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 3;
    level_max = 2;
    levels_index_test ( rule, dim_num, level_max );

    dim_num = 6;
    level_max = 2;
    levels_index_test ( rule, dim_num, level_max );
  }
//
//  Test SPARSE_GRID by having it compute a few sparse grids based on each rule.
//
  if ( true )
  {
    for ( rule = 1; rule <= rule_max; rule++ )
    {
      dim_num = 2;
      level_max = 1;
      sparse_grid_compute_test ( rule, dim_num, level_max );

      dim_num = 2;
      level_max = 2;
      sparse_grid_compute_test ( rule, dim_num, level_max );

      dim_num = 3;
      level_max = 1;
      sparse_grid_compute_test ( rule, dim_num, level_max );
    }
  }
//
//  Test SPARSE_GRID by having it compute a few sparse grids based on each rule,
//  and comparing the sum of the quadrature weights to the expected sum.
//
  if ( true )
  {
    for ( rule = 1; rule <= rule_max; rule++ )
    {
      dim_num = 2;
      level_max = 4;
      sparse_grid_weight_test ( rule, dim_num, level_max );

      dim_num = 3;
      level_max = 0;
      sparse_grid_weight_test ( rule, dim_num, level_max );

      dim_num = 3;
      level_max = 1;
      sparse_grid_weight_test ( rule, dim_num, level_max );

      dim_num = 3;
      level_max = 6;
      sparse_grid_weight_test ( rule, dim_num, level_max );

      dim_num = 10;
      level_max = 3;
      sparse_grid_weight_test ( rule, dim_num, level_max );
    }
  }
//
//  Test SPARSE_GRID by having it compute a few sparse grids based on each rule,
//  and comparing estimated and exact monomial integrals.
//
  if ( true )
  {
    for ( rule = 1; rule <= rule_max; rule++ )
    {
      dim_num = 2;
      level_max = 0;
      degree_max = 3;
      sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max );

      dim_num = 2;
      level_max = 1;
      degree_max = 5;
      sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max );

      dim_num = 2;
      level_max = 2;
      degree_max = 7;
      sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max );

      dim_num = 2;
      level_max = 3;
      degree_max = 9;
      sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max );

      dim_num = 3;
      level_max = 0;
      degree_max = 2;
      sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max );

      dim_num = 3;
      level_max = 1;
      degree_max = 4;
      sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max );

      dim_num = 3;
      level_max = 2;
      degree_max = 6;
      sparse_grid_monomial_test ( rule, dim_num, level_max, degree_max );
    }
  }
//
//  All done.
//
  cout << "\n";
  cout << "SANDIA_SPARSE_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//***************************************************************************80

void levels_index_size_test ( int rule, int dim_min, int dim_max, 
  int level_max_min, int level_max_max )

//***************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_SIZE_TEST tests LEVELS_INDEX_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
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
//    Input, int DIM_MIN, the minimum spatial dimension.
//
//    Input, int DIM_MAX, the maximum spatial dimension.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX.
//
{
  int dim_num;
  int level_max;

  cout << "\n";
  cout << "LEVELS_INDEX_SIZE_TEST\n";
  cout << "  LEVELS_INDEX_SIZE returns the number of distinct\n";
  cout << "  points in a sparse grid derived from a 1D rule.\n";
  cout << "\n";
  cout << "  We are looking at rules like rule " << rule << "\n";
  cout << "\n";
  cout << "  Each sparse grid is of spatial dimension DIM,\n";
  cout << "  and is made up of product grids such that\n";
  cout << "  LEVEL_MIN <= LEVEL <= LEVEL_MAX.\n";

  cout << "\n";
  cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
  {
    cout << "  " << setw(8) << dim_num;
  }
  cout << "\n";
  cout << "   LEVEL_MAX\n";
  cout << "   ---------\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    cout << "    " << setw(4) << level_max;
    for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
    {
      cout << "  " << setw(8) << levels_index_size ( dim_num, level_max, rule );
    }
    cout << "\n";
  }

  return;
}
//***************************************************************************80

void levels_index_test ( int rule, int dim_num, int level_max )

//***************************************************************************80
//
//  Purpose:
//
//    LEVELS_INDEX_TEST tests LEVELS_INDEX.
//
//  Discussion:
//
//    The routine computes the indices of the unique points used in a sparse 
//    multidimensional grid whose size is controlled by a parameter LEVEL_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
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
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the level.
//
{
  int dim;
  int *grid_base;
  int *grid_index;
  int level_min;
  int point;
  int point_num;

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  cout << "\n";
  cout << "LEVELS_INDEX_TEST\n";
  cout << "  LEVELS_INDEX returns all grid indexes\n";
  cout << "  whose level value satisfies\n";
  cout << "    LEVEL_MIN <= LEVEL <= LEVEL_MAX.\n";
  cout << "  Here, LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  and the order of the rule is 2^LEVEL + 1.\n";
  cout << "\n";
  cout << "  We are looking at rules like rule " << rule << "\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  LEVEL_MIN =                 " << level_min << "\n";
  cout << "  LEVEL_MAX =                 " << level_max << "\n";

  point_num = levels_index_size ( dim_num, level_max, rule );

  cout << "  Unique points in the grid = " << point_num << "\n";
//
//  Compute the orders and points.
//
  grid_base = new int[dim_num*point_num];
  grid_index = new int[dim_num*point_num];

  levels_index ( dim_num, level_max, rule, point_num, grid_index, 
    grid_base );
//
//  Now we're done.  Print the merged grid data.
//
  cout << "\n";
  cout << " Point     Grid indices:\n";
  cout << "           Grid bases:\n";
  cout << "\n";
  for ( point = 0; point < point_num; point++ )
  {
    cout << "  " << setw(4) << point << "  ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(6) << grid_index[dim+point*dim_num];
    }
    cout << "\n";
    cout << "        ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(6) << grid_base[dim+point*dim_num];
    }
    cout << "\n";
  }

  delete [] grid_base;
  delete [] grid_index;

  return;
}
//***************************************************************************80

void sparse_grid_compute_test ( int rule, int dim_num, int level_max )

//***************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_COMPUTE_TEST computes and prints a sparse grid rule.
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
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the level.
//
{
  int dim;
  double *grid_point;
  double *grid_weight;
  int level_min;
  int point;
  int point_num;

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  cout << "\n";
  cout << "SPARSE_GRID_COMPUTE_TEST:\n";
  cout << "  SPARSE_GRID can make a sparse grid.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  LEVEL_MIN =                 " << level_min << "\n";
  cout << "  LEVEL_MAX =                 " << level_max << "\n";
  cout << "  1D quadrature index RULE =  " << rule << "\n";
//
//  Determine the number of points.
//
  point_num = levels_index_size ( dim_num, level_max, rule );

  cout << "  Unique points in the grid = " << point_num << "\n";
//
//  Allocate space for the weights and points.
//
  grid_weight = new double[point_num];
  grid_point = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  sparse_grid ( dim_num, level_max, rule, point_num, grid_weight, 
    grid_point );
//
//  Print them out.
//
  cout << "\n";
  cout << "  Grid weights:\n";
  cout << "\n";
  for ( point = 0; point < point_num; point++ )
  {
    cout << "  " << setw(4) << point 
         << "  " << setw(14) << grid_weight[point] << "\n";
  }
  
  cout << "\n";
  cout << "  Grid points:\n";
  cout << "\n";
  for ( point = 0; point < point_num; point++ )
  {
    cout << "  " << point << "  ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(14) << grid_point[dim+point*dim_num];
    }
    cout << "\n";
  }

  delete [] grid_point;
  delete [] grid_weight;

  return;
}
//***************************************************************************80

void sparse_grid_weight_test ( int rule, int dim_num, int level_max )

//***************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_WEIGHT_TEST checks the sum of the quadrature weights.
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
//    31 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
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
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the level.
//
{
  double *grid_point;
  double *grid_weight;
  int level_min;
  double pi = 3.141592653589793;
  int point;
  int point_num;
  double weight_sum;
  double weight_sum_error;
  double weight_sum_exact;

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  if ( 1 <= rule && rule <= 5 )
  {
    weight_sum_exact = pow ( 2.0, dim_num );
  }
  else if ( 6 == rule )
  {
    weight_sum_exact = sqrt ( pow ( pi, dim_num ) );
  }
  else if ( 7 == rule )
  {
    weight_sum_exact = 1.0;
  }

  cout << "\n";
  cout << "SPARSE_GRID_WEIGHT_TEST:\n";
  cout << "  Compute the weights of a sparse grid.\n";
  cout << "\n";
  cout << "  As a simple test, sum these weights.\n";
  cout << "  They should sum to exactly " << weight_sum_exact << "\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  LEVEL_MIN =                 " << level_min << "\n";
  cout << "  LEVEL_MAX =                 " << level_max << "\n";
  cout << "  1D quadrature index RULE =  " << rule << "\n";
//
//  Determine the number of points.
//
  point_num = levels_index_size ( dim_num, level_max, rule );

  cout << "  Unique points in the grid = " << point_num << "\n";
//
//  Allocate space for the weights and points.
//
  grid_weight = new double[point_num];
  grid_point = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  sparse_grid ( dim_num, level_max, rule, point_num, grid_weight, 
    grid_point );
//
//  Sum the weights.
//
  weight_sum = 0.0;
  for ( point = 0; point < point_num; point++ )
  {
    weight_sum = weight_sum + grid_weight[point];
  }
   
  weight_sum_error = r8_abs ( weight_sum - weight_sum_exact );
  
  cout << "\n";
  cout << "    Weight sum  Expected sum    Difference\n";
  cout << "\n";
  cout << "  " << setw(14) << weight_sum
       << "  " << setw(14) << weight_sum_exact
       << "  " << setw(14) << weight_sum_error << "\n";

  delete [] grid_point;
  delete [] grid_weight;

  return;
}
//***************************************************************************80

void sparse_grid_monomial_test ( int rule, int dim_num, int level_max, 
  int degree_max )

//***************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_MONOMIAL_TEST tests monomial exactness of the sparse grid rules.
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
//    Input, int RULE, the index of the rule.
//    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
//    2, "F1", Fejer 1 Open Fully Nested rule.
//    3, "F2", Fejer 2 Open Fully Nested rule.
//    4, "GP", Gauss Patterson Open Fully Nested rule.
//    5, "GL", Gauss Legendre Open Weakly Nested rule.
//    6, "GH", Gauss Hermite Open Weakly Nested rule.
//    7, "LG", Gauss Laguerre Open Non Nested rule.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the level.
//
//    Input, int DEGREE_MAX, the maximum monomial total 
//    degree to check.
//
{
  int degree;
  int dim;
  int *expon;
  double *grid_point;
  double *grid_weight;
  int h;
  int last;
  int level_min;
  bool more;
  int point;
  int point_num;
  double quad_error;
  int t;
  double volume;

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  cout << "\n";
  cout << "SPARSE_GRID_MONOMIAL_TEST\n";
  cout << "  Check the exactness of a sparse grid quadrature rule,\n";
  cout << "  applied to all monomials of orders 0 to DEGREE_MAX.\n";
  cout << "\n";
  cout << "  For cases where the dimension is greater than 1,\n";
  cout << "  many sparse grid of this level have accuracy through\n";
  cout << "  monomials of total degree   " << 2 * level_max + 1 << "\n";;
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";;
  cout << "  LEVEL_MIN =                 " << level_min << "\n";;
  cout << "  LEVEL_MAX =                 " << level_max << "\n";;
  cout << "  1D quadrature index RULE =  " << rule << "\n";;
  cout << "  Check up to DEGREE_MAX =    " << degree_max << "\n";
//
//  Determine the number of points in the rule.
//
  point_num = levels_index_size ( dim_num, level_max, rule );

  cout << "  Unique points in the grid = " << point_num << "\n";
//
//  Allocate space for the weights and points.
//
  grid_weight = new double[point_num];
  grid_point = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  sparse_grid ( dim_num, level_max, rule, point_num, grid_weight, grid_point );
//
//  Compare exact and estimated values of the integrals of various monomials.
//
  expon = new int[dim_num];

  cout << "\n";
  cout << "      Error      Total   Monomial\n";
  cout << "                 Degree  Exponents\n";
  cout << "\n";

  for ( degree = 0; degree <= degree_max; degree++ )
  {
    more = 0;

    for ( ; ; )
    {
      comp_next ( degree, dim_num, expon, &more, &h, &t );

      quad_error = monomial_quadrature ( dim_num, expon, point_num, 
        grid_weight, grid_point, rule );

      cout << "  " << setw(14) << quad_error
           << "  " << setw(2) << degree
           << "  ";

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

    if ( 1 < dim_num )
    {
      cout << "\n";
    }
  }

  delete [] expon;
  delete [] grid_point;
  delete [] grid_weight;

  return;
}
