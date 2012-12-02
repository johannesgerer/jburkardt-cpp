# include "sandia_rules.hpp"
# include "sparse_grid_mixed.hpp"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

int main ( );

void sparse_grid_mixed_weight_tests ( double tol );
void sparse_grid_mixed_weight_test ( int dim_num, int level_max_min,
  int level_max_max, int rule[], double alpha[], double beta[], double tol );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_MIXED_WEIGHT_PRB.
//
//  Discussion:
//
//    SPARSE_GRID_MIXED_WEIGHT_PRB tests SPARSE_GRID_MIXED_WEIGHT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2009
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
  double tol;

  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_WEIGHT_PRB\n";
  std::cout << "  C++ version\n";

  sparse_grid_mixed_weight_tests ( tol );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_WEIGHT_PRB\n";
  std::cout << "  Normal end of execution.\n";

  std::cout << "\n";
  webbur::timestamp ( );

  return 0;
}
//****************************************************************************80

void sparse_grid_mixed_weight_tests ( double tol )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_MIXED_WEIGHT_TESTS calls SPARSE_GRID_MIXED_WEIGHT_TEST.
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
//  Parameters:
//
//    Input, double TOL, a tolerance for point equality.
//    A value of sqrt ( eps ) is reasonable, and will allow the code to
//    consolidate points which are equal, or very nearly so.  A value of
//    -1.0, on the other hand, will force the code to use every point, regardless
//    of duplication.
//
{
  double *alpha;
  double *beta;
  int dim_num;
  int level_max_max;
  int level_max_min;
  int *order_1d;
  int order_nd;
  int *rule;

  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_WEIGHT_TESTS\n";
  std::cout << "  Call SPARSE_GRID_MIXED_WEIGHT_TEST with various arguments.\n";
  std::cout << "  All tests will use a point equality tolerance of " << tol << "\n";

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 1;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 3;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 4;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 7;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 1.5;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 8;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.5;
  beta[0] = 0.0;
  beta[1] = 1.5;
  rule[0] = 2;
  rule[1] = 9;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 2.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 6;
  rule[1] = 4;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 3;
  level_max_min = 0;
  level_max_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  beta[2] = 0.0;
  rule[0] = 1;
  rule[1] = 2;
  rule[2] = 5;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 3.
//
  dim_num = 2;
  level_max_min = 0;
  level_max_max = 4;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 3;
  rule[1] = 3;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 13.
//
  dim_num = 2;
  level_max_min = 0;
  level_max_max = 4;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 13;
  rule[1] = 13;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 16.
//
  dim_num = 2;
  level_max_min = 0;
  level_max_max = 4;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 16;
  rule[1] = 16;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 17.
//
  dim_num = 2;
  level_max_min = 0;
  level_max_max = 4;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 17;
  rule[1] = 17;
  sparse_grid_mixed_weight_test ( dim_num, level_max_min, level_max_max, rule,
    alpha, beta, tol );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  return;
}
//***************************************************************************80

void sparse_grid_mixed_weight_test ( int dim_num, int level_max_min,
  int level_max_max, int rule[], double alpha[], double beta[], double tol )

//***************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_MIXED_WEIGHT_TEST checks the sum of the quadrature weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX_MIN, LEVEL_MAX_MAX, the minimum and
//    maximum values of LEVEL_MAX.
//
//    Input, int RULE[DIM_NUM], the rule in each dimension.
//     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
//     2, "F2",  Fejer Type 2, Open Fully Nested.
//     3, "GP",  Gauss Patterson, Open Fully Nested.
//     4, "GL",  Gauss Legendre, Open Weakly Nested.
//     5, "GH",  Gauss Hermite, Open Weakly Nested.
//     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
//     7, "LG",  Gauss Laguerre, Open Non Nested.
//     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
//     9, "GJ",  Gauss Jacobi, Open Non Nested.
//    10, "GW",  Golub Welsch, (presumed) Open Non Nested.
//    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested.
//    12, "F2_SE", Fejer Type 2 Slow Exponential, Closed Fully Nested.
//    13, "GP_SE", Gauss Patterson Slow Exponential, Closed Fully Nested.
//    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested.
//    15, "F2_ME", Fejer Type 2 Moderate Exponential, Closed Fully Nested.
//    16, "GP_ME", Gauss Patterson Moderate Exponential, Closed Fully Nested.
//    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
//
//    Input, double ALPHA[DIM_NUM], BETA[DIM_NUM], parameters used for
//    Generalized Gauss Hermite, Generalized Gauss Laguerre, and Gauss Jacobi rules.
//
//    Input, double TOL, a tolerance for point equality.
//
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  int dim;
  int level_max;
  double pi = 3.141592653589793;
  int point;
  int point_num;
  int point_total_num;
  int *sparse_unique_index;
  double *sparse_weight;
  double value1;
  double value2;
  double weight_sum;
  double weight_sum_error;
  double weight_sum_exact;

  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_WEIGHT_TEST\n";
  std::cout << "  Compute the weights of a sparse grid.\n";
  std::cout << "\n";
  std::cout << "  Each sparse grid is of spatial dimension DIM_NUM,\n";
  std::cout << "  and is made up of product grids of levels up to LEVEL_MAX.\n";
  std::cout << "\n";
  std::cout << " Dimension      Rule     Alpha          Beta\n";
  std::cout << "\n";

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 2 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 3 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 4 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 5 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 6 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(14) << alpha[dim] << "\n";
    }
    else if ( rule[dim] == 7 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 8 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(14) << alpha[dim] << "\n";
    }
    else if ( rule[dim] == 9 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(14) << alpha[dim]
           << "  " << std::setw(14) << beta[dim] << "\n";
    }
    else if ( rule[dim] == 10 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 11 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 12 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 13 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 14 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 15 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 16 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else if ( rule[dim] == 17 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim] << "\n";
    }
    else
    {
      std::cout << "\n";
      std::cout << "SPARSE_GRID_MIXED_WEIGHT_TEST - Fatal error!\n";
      std::cout << "  Unexpected value of RULE = " << rule[dim] << "\n";
      exit ( 1 );
    }
  }

  weight_sum_exact = 1.0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 2 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 3 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 4 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 5 )
    {
      weight_sum_exact = weight_sum_exact * sqrt ( pi );
    }
    else if ( rule[dim] == 6 )
    {
      weight_sum_exact = weight_sum_exact * webbur::r8_gamma ( 0.5 * ( alpha[dim] + 1.0 ) );
    }
    else if ( rule[dim] == 7 )
    {
      weight_sum_exact = weight_sum_exact * 1.0;
    }
    else if ( rule[dim] == 8 )
    {
      weight_sum_exact = weight_sum_exact * webbur::r8_gamma ( alpha[dim] + 1.0 );
    }
    else if ( rule[dim] == 9 )
    {
      arg1 = - alpha[dim];
      arg2 = 1.0;
      arg3 = beta[dim] + 2.0;
      arg4 = - 1.0;
      value1 = webbur::r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );
      arg1 = - beta[dim];
      arg2 = 1.0;
      arg3 = alpha[dim] + 2.0;
      arg4 = - 1.0;
      value2 = webbur::r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );
      weight_sum_exact = weight_sum_exact * (
        value1 / ( beta[dim] + 1.0 ) + value2 / ( alpha[dim] + 1.0 ) );
    }
    else if ( rule[dim] == 10 )
    {
      std::cerr << "\n";
      std::cerr << "SPARSE_GRID_MIXED_WEIGHT_TEST - Fatal error!\n";
      std::cerr << "  Unexpected value of RULE = 10.\n";
      exit ( 1 );
    }
    else if ( rule[dim] == 11 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 12 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 13 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 14 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 15 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 16 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else if ( rule[dim] == 17 )
    {
      weight_sum_exact = weight_sum_exact * 2.0;
    }
    else
    {
      std::cerr << "\n";
      std::cerr << "SPARSE_GRID_MIXED_WEIGHT_TEST - Fatal error!\n";
      std::cerr << "  Unexpected value of RULE[" << dim << "] = " << rule[dim] << ".\n";
      exit ( 1 );
    }
  }

  std::cout << "\n";
  std::cout << "  As a simple test, sum these weights.\n";
  std::cout << "  They should sum to exactly " << weight_sum_exact << "\n";

  std::cout << "\n";
  std::cout << "     Level      Weight sum  Expected sum    Difference\n";
  std::cout << "\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    point_total_num = webbur::sparse_grid_mixed_size_total ( dim_num, level_max, rule );

    point_num = webbur::sparse_grid_mixed_size ( dim_num, level_max, rule, alpha,
      beta, tol );

    sparse_unique_index = new int[point_total_num];

    webbur::sparse_grid_mixed_unique_index ( dim_num, level_max, rule, alpha, beta,
      tol, point_num, point_total_num, sparse_unique_index );

    sparse_weight = new double[point_num];

    webbur::sparse_grid_mixed_weight ( dim_num, level_max, rule, alpha, beta,
      point_num, point_total_num, sparse_unique_index, sparse_weight );

    weight_sum = webbur::r8vec_sum ( point_num, sparse_weight );

    weight_sum_error = webbur::r8_abs ( weight_sum - weight_sum_exact );

    std::cout << "  " << std::setw(8)  << level_max
         << "  " << std::setw(14) << weight_sum
         << "  " << std::setw(14) << weight_sum_exact
         << "  " << std::setw(14) << weight_sum_error << "\n";

    delete [] sparse_unique_index;
    delete [] sparse_weight;
  }
  return;
}
