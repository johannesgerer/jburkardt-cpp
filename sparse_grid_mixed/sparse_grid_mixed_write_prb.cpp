# include "sandia_rules.hpp"
# include "sparse_grid_mixed.hpp"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

int main ( );

void sparse_grid_mixed_write_tests ( double tol );
void sparse_grid_mixed_write_test ( int dim_num, int level_max, int rule[],
  double alpha[], double beta[], double tol, std::string file_name );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_MIXED_WRITE_PRB.
//
//  Discussion:
//
//    SPARSE_GRID_MIXED_WRITE_PRB tests SPARSE_GRID_MIXED_WRITE.
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
  std::cout << "SPARSE_GRID_MIXED_WRITE_PRB\n";
  std::cout << "  C++ version\n";

  sparse_grid_mixed_write_tests ( tol );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_WRITE_PRB\n";
  std::cout << "  Normal end of execution.\n";

  std::cout << "\n";
  webbur::timestamp ( );

  return 0;
}
//****************************************************************************80

void sparse_grid_mixed_write_tests ( double tol )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_MIXED_WRITE_TESTS calls SPARSE_GRID_MIXED_WRITE_TEST.
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
  std::string file_name;
  int level_max;
  int level_max_max;
  int level_max_min;
  int *order_1d;
  int order_nd;
  int *rule;

  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_WRITE_TESTS\n";
  std::cout << "  Call SPARSE_GRID_MIXED_WRITE_TEST with various arguments.\n";
  std::cout << "  All tests will use a point equality tolerance of " << tol << "\n";

  dim_num = 2;
  level_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 1;
  file_name = "sparse_grid_mixed_d2_l2_ccxcc";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max = 3;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 3;
  file_name = "sparse_grid_mixed_d2_l3_ccxgp";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 4;
  file_name = "sparse_grid_mixed_d2_l2_ccxgl";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 7;
  file_name = "sparse_grid_mixed_d2_l2_ccxlg";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 1.5;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 1;
  rule[1] = 8;
  file_name = "sparse_grid_mixed_d2_l2_ccxglg";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.5;
  beta[0] = 0.0;
  beta[1] = 1.5;
  rule[0] = 2;
  rule[1] = 9;
  file_name = "sparse_grid_mixed_d2_l2_f2xgj";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 2;
  level_max = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 2.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 6;
  rule[1] = 4;
  file_name = "sparse_grid_mixed_d2_l2_gghxgl";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  dim_num = 3;
  level_max = 2;
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
  file_name = "sparse_grid_mixed_d3_l2_ccxf2xgh";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 3.
//
  dim_num = 2;
  level_max = 4;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 3;
  rule[1] = 3;
  file_name = "sparse_grid_mixed_d2_l4_gpxgp";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 13.
//
  dim_num = 2;
  level_max = 4;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 13;
  rule[1] = 13;
  file_name = "sparse_grid_mixed_d2_l4_gpsexgpse";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 16.
//
  dim_num = 2;
  level_max = 4;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 16;
  rule[1] = 16;
  file_name = "sparse_grid_mixed_d2_l4_gpmexgpme";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 17.
//
  dim_num = 2;
  level_max = 4;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  rule[0] = 17;
  rule[1] = 17;
  file_name = "sparse_grid_mixed_d2_l4_ccnxccn";
  sparse_grid_mixed_write_test ( dim_num, level_max, rule, alpha,
    beta, tol, file_name );
  delete [] alpha;
  delete [] beta;
  delete [] rule;

  return;
}
//***************************************************************************80

void sparse_grid_mixed_write_test ( int dim_num, int level_max, int rule[],
  double alpha[], double beta[], double tol, std::string file_name )

//***************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_MIXED_WRITE_TEST tests SPARSE_GRID_MIXED_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer DIM_NUM, the spatial dimension.
//
//    Input, integer LEVEL_MAX, the level that defines the grid.
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
//    Generalized Gauss Hermite, Generalized Gauss Laguerre,
//    and Gauss Jacobi rules.
//
//    Input, double TOL, a tolerance for point equality.
//
//    Input, string FILE_NAME, the main name of the
//    output files.
//
{
  int point_num;
  int point_total_num;
  int *sparse_index;
  int *sparse_order;
  double *sparse_point;
  int *sparse_unique_index;
  double *sparse_weight;

  std::cout << "\n";
  std::cout << "SPARSE_GRID_MIXED_WRITE_TEST\n";
  std::cout << "  SPARSE_GRID_MIXED_WRITE writes a sparse grid rule\n";
  std::cout << "  to X, W and R files.\n";
//
//  Compute necessary data.
//
  point_total_num = webbur::sparse_grid_mixed_size_total ( dim_num, level_max, rule );

  point_num = webbur::sparse_grid_mixed_size ( dim_num, level_max, rule, alpha,
    beta, tol );

  sparse_unique_index = new int[point_total_num];

  webbur::sparse_grid_mixed_unique_index ( dim_num, level_max, rule, alpha, beta,
    tol, point_num, point_total_num, sparse_unique_index );

  sparse_order = new int[dim_num*point_num];
  sparse_index = new int[dim_num*point_num];

  webbur::sparse_grid_mixed_index ( dim_num, level_max, rule, point_num,
    point_total_num, sparse_unique_index, sparse_order, sparse_index );
//
//  Compute points and weights.
//
  sparse_point = new double [ dim_num * point_num ];

  webbur::sparse_grid_mixed_point ( dim_num, level_max, rule, alpha, beta,
    point_num, sparse_order, sparse_index, sparse_point );

  sparse_weight = new double[point_num];

  webbur::sparse_grid_mixed_weight ( dim_num, level_max, rule, alpha, beta,
    point_num, point_total_num, sparse_unique_index, sparse_weight );
//
//  Write points and weights to files.
//
  webbur::sparse_grid_mixed_write ( dim_num, rule, alpha, beta, point_num,
    sparse_weight, sparse_point, file_name );

  delete [] sparse_index;
  delete [] sparse_order;
  delete [] sparse_point;
  delete [] sparse_unique_index;
  delete [] sparse_weight;

  return;
}
