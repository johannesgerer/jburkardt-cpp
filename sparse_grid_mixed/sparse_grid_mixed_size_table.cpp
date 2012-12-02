# include "sandia_rules.hpp"
# include "sparse_grid_mixed.hpp"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

int main ( );

void sparse_grid_mixed_size_tabulate ( int rule_1d, double alpha_1d,
  double beta_1d, int dim_min, int dim_max, int level_max_min,
  int level_max_max );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_MIXED_SIZE_TABLE.
//
//  Discussion:
//
//    SPARSE_GRID_MIXED_SIZE_TABLE_PRB makes a point count table.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2011
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
  double alpha_1d;
  double beta_1d;
  int dim_max;
  int dim_min;
  int level_max_max;
  int level_max_min;
  double rule_1d;

  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_SIZE_TABLE\n";
  std::cout << "  C++ version\n";
  std::cout << "  Tabulate sparse grid size for various rules, levels, dimensions.\n";
//
//  "Rule 0"
//  A special input that indicates we want a table of the number of polynomials
//  of degree DEGREE or less in a space of dimension DIM.
//
  rule_1d = 0;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );

  rule_1d = 0;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );
//
//  Rule 1.
//  Clenshaw-Curtis Exponential Growth
//
  rule_1d = 1;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );

  rule_1d = 1;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );

  rule_1d = 1;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 100;
  dim_max = 100;
  level_max_min = 0;
  level_max_max = 2;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );
//
//  Rule 3.
//  Gauss Patterson Exponential Growth
//
  rule_1d = 3;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );

  rule_1d = 3;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );
//
//  Rule 11
//  Clenshaw-Curtis Slow Exponential Growth
//
  rule_1d = 11;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );

  rule_1d = 11;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );
//
//  Rule 13
//  Gauss Patterson Slow Exponential Growth
//
  rule_1d = 13;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );

  rule_1d = 13;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );
//
//  Rule 16
//  Gauss Patterson Moderate Exponential Growth
//
  rule_1d = 16;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );

  rule_1d = 16;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );
//
//  Rule 17
//
  rule_1d = 17;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );

  rule_1d = 17;
  alpha_1d = 0.0;
  beta_1d = 0.0;
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;

  sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min,
    dim_max, level_max_min, level_max_max );

  std::cout << "\n";
  webbur::timestamp ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_SIZE_TABLE\n";
  std::cout << "  Normal end of execution.\n";

  std::cout << "\n";
  webbur::timestamp ( );

  return 0;
}
//***************************************************************************80

void sparse_grid_mixed_size_tabulate ( int rule_1d, double alpha_1d,
  double beta_1d, int dim_min, int dim_max, int level_max_min,
  int level_max_max )

//***************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_MIXED_SIZE_TABULATE tests SPARSE_GRID_MIXED_SIZE.
//
//  Discussion:
//
//    We do NOT consider mixed rules.  Instead, we are looking at sparse grid
//    rules for which all dimensions use the same 1D rule family.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int RULE_1D, the 1D rule.
//
//    Input, double ALPHA_1D, BETA_1D, the optional parameters.
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
  double *alpha;
  double *beta;
  int dim;
  int dim_num;
  int level_max;
  int point_num;
  int *rule;
  double tol;

  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_SIZE_TABULATE\n";
  std::cout << "  SPARSE_GRID_MIXED_SIZE returns the number of distinct\n";
  std::cout << "  points in a sparse grid.\n";

  if ( rule_1d == 0 )
  {
    std::cout << "\n";
    std::cout << "  Here we report the total number of polynomials of\n";
    std::cout << "  degree DEGREE or less in the given DIM dimensional space.\n";
    std::cout << "\n";

    std::cout << "   DIM: ";
    for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
    {
      std::cout << "  " << std::setw(8) << dim_num;
    }
    std::cout<< "\n";
    std::cout << "\n";
    std::cout << "   DEGREE\n";
    std::cout << "\n";

    for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
    {
      std::cout << "    " << std::setw(4) << level_max;

      for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
      {
        point_num = webbur::i4_choose ( dim_num + level_max, dim_num );

        std::cout << "  " << std::setw(8) << point_num;
      }
      std::cout << "\n";
    }
  }
  else
  {
    std::cout << "\n";
    std::cout << "  We use the same rule in all dimensions, and count the points\n";
    std::cout << "  for a range of dimensions and levels.\n";
    std::cout << "\n";
    std::cout << "  1D rule index is " << rule_1d << "\n";
    std::cout << "  ALPHA parameter is " << alpha_1d << "\n";
    std::cout << "  BETA parameter is  " << beta_1d << "\n";
    std::cout << "\n";

    tol = std::sqrt ( webbur::r8_epsilon ( ) );

    std::cout << "   DIM: ";
    for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
    {
      std::cout << "  " << std::setw(8) << dim_num;
    }
    std::cout<< "\n";
    std::cout << "\n";
    std::cout << "   LEVEL_MAX\n";
    std::cout << "\n";

    for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
    {
      std::cout << "    " << std::setw(4) << level_max;

      for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
      {
        rule = new int[dim_num];
        alpha = new double[dim_num];
        beta = new double[dim_num];

        for ( dim = 0; dim < dim_num; dim++ )
        {
          rule[dim] = rule_1d;
          alpha[dim] = alpha_1d;
          beta[dim] = beta_1d;
        }

        point_num = webbur::sparse_grid_mixed_size ( dim_num, level_max, rule,
          alpha, beta, tol );

        std::cout << "  " << std::setw(8) << point_num;

        delete [] alpha;
        delete [] beta;
        delete [] rule;
      }
      std::cout << "\n";
    }
  }
  return;
}
