# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "sparse_grid_open.hpp"

int main ( );

void test01 ( int dim_min, int dim_max, int level_max_min, int level_max_max );
void test011 ( int dim_min, int dim_max, int level_max_min, int level_max_max );
void test012 ( int dim_min, int dim_max, int level_max_min, int level_max_max );
void test013 ( int dim_min, int dim_max, int level_max_min, int level_max_max );
void test015 ( int dim_min, int dim_max, int level_max_min, int level_max_max );
void test02 ( int dim_num, int level_max );
void test03 ( int dim_num, int level_max );
void test04 ( int dim_num, int level_max );
void test05 ( int dim_num, int level_max );
void test06 ( int dim_num, int level_max );
void test08 ( int dim_num, int level_max );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_OPEN_PRB.
//
//  Discussion:
//
//    SPARSE_GRID_OPEN_PRB tests the SPARSE_GRID_OPEN library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2010
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
  int level_max_max;
  int level_max_min;

  timestamp ( );
  cout << "\n";
  cout << "SPARSE_GRID_OPEN_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPARSE_GRID_OPEN library.\n";
//
//  Point counts for Open Fully Nested "OFN" rules.
//
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 10;
  test01 ( dim_min, dim_max, level_max_min, level_max_max );

  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 10;
  test01 ( dim_min, dim_max, level_max_min, level_max_max );
//
//  Point counts for Open Non Nested "ONN" rules.
//
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 10;
  test011 ( dim_min, dim_max, level_max_min, level_max_max );

  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 10;
  test011 ( dim_min, dim_max, level_max_min, level_max_max );
//
//  Point counts for Open Weakly Nested "OWN" rules.
//
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 10;
  test012 ( dim_min, dim_max, level_max_min, level_max_max );

  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 10;
  test012 ( dim_min, dim_max, level_max_min, level_max_max );
//
//  Point counts for Fejer Type 2 Slow rules.
//
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 10;
  test013 ( dim_min, dim_max, level_max_min, level_max_max );

  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 10;
  test013 ( dim_min, dim_max, level_max_min, level_max_max );
//
//  Point counts for Gauss-Patterson-Slow rules.
//
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 10;
  test015 ( dim_min, dim_max, level_max_min, level_max_max );

  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 10;
  test015 ( dim_min, dim_max, level_max_min, level_max_max );

  test02 ( 2, 2 );
  test02 ( 2, 3 );
  test02 ( 2, 4 );
  test02 ( 3, 2 );
  test02 ( 6, 2 );

  test04 ( 2, 3 );

  test05 ( 2, 3 );

  test06 ( 2, 4 );

  test08 ( 2, 4 );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPARSE_GRID_OPEN_PRB\n";
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
//    TEST01 tests SPARSE_GRID_OFN_SIZE.
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
  cout << "  SPARSE_GRID_OFN_SIZE returns the number of \n";
  cout << "  distinct points in a sparse grid made up of all \n";
  cout << "  product grids formed from open fully nested\n";
  cout << "  quadrature rules.\n";
  cout << "\n";
  cout << "  The sparse grid is the sum of all product grids\n";
  cout << "  of order LEVEL, with\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "\n";
  cout << "  LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  the order of the 1D rule is 2^(LEVEL+1) - 1,\n";
  cout << "  the region is [-1,1]^DIM_NUM.\n";
  cout << "\n";
  cout << "  For this kind of rule, there is complete nesting,\n";
  cout << "  that is, a sparse grid of a given level includes\n";
  cout << "  ALL the points on grids of lower levels.\n";
  cout << "\n";
  cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    cout << "  " << setw(10) << dim_num;
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
      point_num = sparse_grid_ofn_size ( dim_num, level_max );
      cout << "  " << setw(10) << point_num;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test011 ( int dim_min, int dim_max, int level_max_min, int level_max_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST011 tests SPARSE_GRID_ONN_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2010
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
  cout << "TEST011\n";
  cout << "  SPARSE_GRID_ONN_SIZE returns the number of \n";
  cout << "  distinct points in a sparse grid made up of all \n";
  cout << "  product grids formed from open non nested \n";
  cout << "  quadrature rules.\n";
  cout << "\n";
  cout << "  The sparse grid is the sum of all product grids\n";
  cout << "  of order LEVEL, with\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "\n";
  cout << "  LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  the order of the 1D rule is 2^(LEVEL+1) - 1,\n";
  cout << "  the region is [-1,1]^DIM_NUM.\n";
  cout << "\n";
  cout << "  For this kind of rule, there is no nesting.\n";
  cout << "\n";
  cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    cout << "  " << setw(10) << dim_num;
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
      point_num = sparse_grid_onn_size ( dim_num, level_max );
      cout << "  " << setw(10) << point_num;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test012 ( int dim_min, int dim_max, int level_max_min, int level_max_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST012 tests SPARSE_GRID_OWN_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 January 2010
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
  cout << "TEST012\n";
  cout << "  SPARSE_GRID_OWN_SIZE returns the number of \n";
  cout << "  distinct points in a sparse grid made up of all \n";
  cout << "  product grids formed from open weakly nested \n";
  cout << "  quadrature rules.\n";
  cout << "\n";
  cout << "  The sparse grid is the sum of all product grids\n";
  cout << "  of order LEVEL, with\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "\n";
  cout << "  LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  the order of the 1D rule is 2^(LEVEL+1) - 1,\n";
  cout << "  the region is [-1,1]^DIM_NUM.\n";
  cout << "\n";
  cout << "  For this kind of rule, there is weak nesting,\n";
  cout << "  that is, 0.0 is the only point any two rules have in common.\n";
  cout << "\n";
  cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    cout << "  " << setw(10) << dim_num;
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
      point_num = sparse_grid_own_size ( dim_num, level_max );
      cout << "  " << setw(10) << point_num;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test013 ( int dim_min, int dim_max, int level_max_min, int level_max_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST013 tests SPARSE_GRID_F2S_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 December 2009
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
  cout << "TEST013\n";
  cout << "  SPARSE_GRID_F2S_SIZE returns the number of \n";
  cout << "  distinct points in a sparse grid made up of all \n";
  cout << "  product grids formed from Fejer Type 2 Slow \n";
  cout << "  quadrature rules.\n";
  cout << "\n";
  cout << "  The sparse grid is the sum of all product grids\n";
  cout << "  of order LEVEL, with\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "\n";
  cout << "  LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  the order of the 1D rule is 2^(LEVEL+1) - 1,\n";
  cout << "  the region is [-1,1]^DIM_NUM.\n";
  cout << "\n";
  cout << "  For this kind of rule, there is complete nesting,\n";
  cout << "  that is, a sparse grid of a given level includes\n";
  cout << "  ALL the points on grids of lower levels.\n";
  cout << "\n";
  cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    cout << "  " << setw(10) << dim_num;
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
      point_num = sparse_grid_f2s_size ( dim_num, level_max );
      cout << "  " << setw(10) << point_num;
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
//    TEST015 tests SPARSE_GRID_GPS_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2009
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
  cout << "TEST015\n";
  cout << "  SPARSE_GRID_GPS_SIZE returns the number of \n";
  cout << "  distinct points in a sparse grid made up of all \n";
  cout << "  product grids formed from Gauss-Patterson-Slow \n";
  cout << "  quadrature rules.\n";
  cout << "\n";
  cout << "  The sparse grid is the sum of all product grids\n";
  cout << "  of order LEVEL, with\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "\n";
  cout << "  LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  the order of the 1D rule is 2^(LEVEL+1) - 1,\n";
  cout << "  the region is [-1,1]^DIM_NUM.\n";
  cout << "\n";
  cout << "  For this kind of rule, there is complete nesting,\n";
  cout << "  that is, a sparse grid of a given level includes\n";
  cout << "  ALL the points on grids of lower levels.\n";
  cout << "\n";
  cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    cout << "  " << setw(10) << dim_num;
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
      point_num = sparse_grid_gps_size ( dim_num, level_max );
      cout << "  " << setw(10) << point_num;
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
//    TEST02 tests LEVELS_OPEN_INDEX.
//
//  Discussion:
//
//    The routine under study computes the indices of the unique points
//    used in a sparse multidimensional grid whose size is controlled
//    by a parameter LEVEL.
//
//    Once these indices are returned, they can be converted into
//    Fejer type 2 points, Gauss-Patterson, 
//    Newton-Cotes Open, scaled to the appropriate 
//    interval.
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
{
  int dim;
  int *grid_index;
  int grid_num;
  int j;
  int level;
  int point_num;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  LEVELS_OPEN_INDEX returns all grid indexes\n";
  cout << "  whose level value satisfies\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "  Here, LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  and the order of the rule is 2^(LEVEL+1) - 1.\n";

  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  point_num = sparse_grid_ofn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  Number of unique points in the grid = " << point_num << "\n";
//
//  Compute the orders and points.
//
  grid_index = levels_open_index ( dim_num, level_max, point_num );
//
//  Now we're done.  Print the merged grid data.
//
  cout << "\n";
  cout << "  Grid index:\n";
  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j << "  ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(6) << grid_index[dim+j*dim_num];
    }
    cout << "\n";
  }

  delete [] grid_index;

  return;
}
//****************************************************************************80

void test04 ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests LEVELS_OPEN_INDEX to create a Fejer Type 2 grid.
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
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int *grid_index;
  int grid_num;
  double *grid_point;
  int j;
  int level;
  int order_max;
  int point_num;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  Make a sparse Fejer Type 2 grid.\n";
  cout << "\n";
  cout << "  LEVELS_OPEN_INDEX returns all grid indexes\n";
  cout << "  whose level value satisfies\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "  Here, LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  and the order of the rule is 2^(LEVEL+1) - 1.\n";
  cout << "\n";
  cout << "  Now we demonstrate how to convert grid indices\n";
  cout << "  into physical grid points.  In this case, we\n";
  cout << "  want points on [-1,+1]^DIM_NUM.\n";
  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  point_num = sparse_grid_ofn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  Number of unique points in the grid = " << point_num << "\n";
//
//  Compute the grid indices.
//
  grid_index = levels_open_index ( dim_num, level_max, point_num );
//
//  Print the grid indices.
//
  cout << "\n";
  cout << "  Grid index:\n";
  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(6) << grid_index[dim+j*dim_num];
    }
    cout << "\n";
  }
//
//  Convert index information to physical information.
//
  order_max = i4_power ( 2, level_max + 1 ) - 1;

  grid_point = new double[dim_num*point_num];

  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      grid_point[dim+j*dim_num] = 
        f2_abscissa ( order_max, grid_index[dim+j*dim_num] );
    }
  }

  cout << "\n";
  cout << "  Grid points:\n";
  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(10) << grid_point[dim+j*dim_num];
    }
    cout << "\n";
  }

  delete [] grid_index;
  delete [] grid_point;

  return;
}
//****************************************************************************80

void test05 ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests LEVELS_OPEN_INDEX to create a Gauss Patterson grid.
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
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int *grid_index;
  int grid_num;
  double *grid_point;
  int j;
  int level;
  int order_max;
  int point_num;

  cout << "\n";
  cout << "TEST05:\n";
  cout << "  Make a sparse Gauss Patterson grid.\n";
  cout << "\n";
  cout << "  LEVELS_OPEN_INDEX returns all grid indexes\n";
  cout << "  whose level value satisfies\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "  Here, LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  and the order of the rule is 2^(LEVEL+1) - 1.\n";
  cout << "\n";
  cout << "  Now we demonstrate how to convert grid indices\n";
  cout << "  into physical grid points.  In this case, we\n";
  cout << "  want points on [-1,+1]^DIM_NUM.\n";
  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  point_num = sparse_grid_ofn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  Number of unique points in the grid = " << point_num << "\n";
//
//  Compute the grid indices.
//
  grid_index = levels_open_index ( dim_num, level_max, point_num );
//
//  Print the grid indices.
//
  cout << "\n";
  cout << "  Grid index:\n";
  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(6) << grid_index[dim+j*dim_num];
    }
    cout << "\n";
  }
//
//  Convert index information to physical information.
//  Note that GP_ABSCISSA expects the LEVEL value, not the ORDER!
//
  grid_point = new double[dim_num*point_num];

  order_max = i4_power ( 2, level_max + 1 ) - 1;

  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      grid_point[dim+j*dim_num] = 
        gp_abscissa ( order_max, grid_index[dim+j*dim_num] );
    }
  }

  cout << "\n";
  cout << "  Grid points:\n";
  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(10) << grid_point[dim+j*dim_num];
    }
    cout << "\n";
  }

  delete [] grid_index;
  delete [] grid_point;

  return;
}
//****************************************************************************80

void test06 ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests levels_open_INDEX to create a Newton Cotes Open grid.
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
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int *grid_index;
  int grid_num;
  double *grid_point;
  int j;
  int level;
  int order_max;
  int point_num;

  cout << "\n";
  cout << "TEST06:\n";
  cout << "  Make a sparse Newton Cotes Open grid.\n";
  cout << "\n";
  cout << "  LEVELS_OPEN_INDEX returns all grid indexes\n";
  cout << "  whose level value satisfies\n";
  cout << "    0 <= LEVEL <= LEVEL_MAX.\n";
  cout << "  Here, LEVEL is the sum of the levels of the 1D rules,\n";
  cout << "  and the order of the rule is 2^(LEVEL+1) - 1.\n";
  cout << "\n";
  cout << "  Now we demonstrate how to convert grid indices\n";
  cout << "  into physical grid points.  In this case, we\n";
  cout << "  want points on [0,+1]^DIM_NUM.\n";
  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  point_num = sparse_grid_ofn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  Number of unique points in the grid = " << point_num << "\n";
//
//  Compute the grid indices.
//
  grid_index = levels_open_index ( dim_num, level_max, point_num );
//
//  Print the grid indices.
//
  cout << "\n";
  cout << "  Grid index:\n";
  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(6) << grid_index[dim+j*dim_num];
    }
    cout << "\n";
  }
//
//  Convert index information to physical information.
//
  order_max = i4_power ( 2, level_max + 1 ) - 1;

  grid_point = new double[dim_num*point_num];

  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      grid_point[dim+j*dim_num] = 
        nco_abscissa ( order_max, grid_index[dim+j*dim_num] );
    }
  }

  cout << "\n";
  cout << "  Grid points:\n";
  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(10) << grid_point[dim+j*dim_num];
    }
    cout << "\n";
  }

  delete [] grid_index;
  delete [] grid_point;

  return;
}
//****************************************************************************80

void test08 ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 creates and writes sparse grid files of all types.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2009
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  string file_name;
  int *grid_index;
  int grid_num;
  double *grid_point;
  int j;
  int level;
  int order_max;
  int point_num;

  cout << "\n";
  cout << "TEST08:\n";
  cout << "  Make sparse grids and write to files.\n";
  cout << "\n";
  cout << "  LEVEL_MAX = " << level_max << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  point_num = sparse_grid_ofn_size ( dim_num, level_max );

  cout << "\n";
  cout << "  Number of unique points in the grid = " << point_num << "\n";
//
//  Compute the grid indices.
//
  grid_index = levels_open_index ( dim_num, level_max, point_num );
//
//  Print the grid indices.
//
  cout << "\n";
  cout << "  Grid index:\n";
  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << "  " << setw(4) << j;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << setw(6) << grid_index[dim+j*dim_num];
    }
    cout << "\n";
  }
//
//  Convert index information to physical information.
//
  order_max = i4_power ( 2, level_max + 1 ) - 1;

  grid_point = new double[dim_num*point_num];
//
//  Create F2 data and write to file.
//
  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      grid_point[dim+j*dim_num] = 
        f2_abscissa ( order_max, grid_index[dim+j*dim_num] );
    }
  }
  file_name = "f2_d" + i4_to_string ( dim_num, "%d" )
    + "_level" + i4_to_string ( level_max, "%d" ) + ".txt";

  r8mat_write ( file_name, dim_num, point_num, grid_point );

  cout << "  Wrote file \"" << file_name << "\".\n";
//
//  Create GP data and write to file.
//  Note that GP_ABSCISSA expects the value of LEVEL, not ORDER!
//
  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      grid_point[dim+j*dim_num] = 
        gp_abscissa ( order_max, grid_index[dim+j*dim_num] );
    }
  }
  file_name = "gp_d" + i4_to_string ( dim_num, "%d" )
    + "_level" + i4_to_string ( level_max, "%d" ) + ".txt";

  r8mat_write ( file_name, dim_num, point_num, grid_point );

  cout << "  Wrote file \"" << file_name << "\".\n";
//
//  Create NCO data and write to file.
//
  for ( j = 0; j < point_num; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      grid_point[dim+j*dim_num] = 
        nco_abscissa ( order_max, grid_index[dim+j*dim_num] );
    }
  }
  file_name = "nco_d" + i4_to_string ( dim_num, "%d" )
    + "_level" + i4_to_string ( level_max, "%d" ) + ".txt";
  
  r8mat_write ( file_name, dim_num, point_num, grid_point );

  cout << "  Wrote file \"" << file_name << "\".\n";

  delete [] grid_index;
  delete [] grid_point;

  return;
}

