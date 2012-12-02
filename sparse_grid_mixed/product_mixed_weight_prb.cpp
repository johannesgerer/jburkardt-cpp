# include "sandia_rules.hpp"
# include "sparse_grid_mixed.hpp"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

int main ( );

void product_mixed_weight_tests ( );
void product_mixed_weight_test ( int dim_num, int order_1d[], int order_nd,
  int rule[], double alpha[], double beta[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PRODUCT_MIXED_WEIGHT_PRB.
//
//  Discussion:
//
//    PRODUCT_MIXED_WEIGHT_PRB tests PRODUCT_MIXED_WEIGHT.
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
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "PRODUCT_MIXED_WEIGHT_PRB\n";
  std::cout << "  C++ version\n";

  product_mixed_weight_tests ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "PRODUCT_MIXED_WEIGHT_PRB\n";
  std::cout << "  Normal end of execution.\n";

  std::cout << "\n";
  webbur::timestamp ( );

  return 0;
}
//****************************************************************************80

void product_mixed_weight_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    PRODUCT_MIXED_WEIGHT_TESTS calls PRODUCT_MIXED_WEIGHT_TEST with various arguments.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *alpha;
  double *beta;
  int dim_num;
  int *order_1d;
  int order_nd;
  int *rule;

  std::cout << "\n";
  std::cout << "PRODUCT_MIXED_WEIGHT_TESTS\n";
  std::cout << "  Call PRODUCT_MIXED_WEIGHT_TEST with various arguments.\n";

  dim_num = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  order_1d[0] = 3;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 1;
  rule[1] = 1;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;

  dim_num = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  order_1d[0] = 3;
  order_1d[1] = 7;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 1;
  rule[1] = 5;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;

  dim_num = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  order_1d[0] = 3;
  order_1d[1] = 3;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 3;
  rule[1] = 7;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;

  dim_num = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 1.5;
  beta[0] = 0.0;
  beta[1] = 0.0;
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 1;
  rule[1] = 8;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;

  dim_num = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.5;
  beta[0] = 0.0;
  beta[1] = 1.5;
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 2;
  rule[1] = 9;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;

  dim_num = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 2.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  order_1d[0] = 7;
  order_1d[1] = 7;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 6;
  rule[1] = 4;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;

  dim_num = 3;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  beta[2] = 0.0;
  order_1d[0] = 2;
  order_1d[1] = 3;
  order_1d[2] = 3;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 1;
  rule[1] = 3;
  rule[2] = 5;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;
//
//  Dimension 2, Rule 13
//
  dim_num = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  order_1d[0] = 15;
  order_1d[1] = 15;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 13;
  rule[1] = 13;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;
//
//  Dimension 2, Rule 16
//
  dim_num = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  order_1d[0] = 15;
  order_1d[1] = 15;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 16;
  rule[1] = 16;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;
//
//  Dimension 2, Rule 17
//
  dim_num = 2;
  alpha = new double[dim_num];
  beta = new double[dim_num];
  order_1d = new int[dim_num];
  rule = new int[dim_num];
  alpha[0] = 0.0;
  alpha[1] = 0.0;
  beta[0] = 0.0;
  beta[1] = 0.0;
  order_1d[0] = 15;
  order_1d[1] = 15;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule[0] = 17;
  rule[1] = 17;
  product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta );
  delete [] alpha;
  delete [] beta;
  delete [] order_1d;
  delete [] rule;

  return;
}
//***************************************************************************80

void product_mixed_weight_test ( int dim_num, int order_1d[], int order_nd,
  int rule[], double alpha[], double beta[] )

//***************************************************************************80
//
//  Purpose:
//
//    PRODUCT_MIXED_WEIGHT_TEST computes the weights of a mixed factor product rule.
//
//  Discussion:
//
//    This routine gets the sparse grid indices and determines the
//    corresponding sparse grid weights.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 February 2010
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
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  int dim;
  double pi = 3.141592653589793;
  double value1;
  double value2;
  double *weight;
  double weight_sum;
  double weight_sum_error;
  double weight_sum_exact;

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
      std::cerr << "PRODUCT_MIXED_WEIGHT_TEST - Fatal error!\n";
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
      std::cerr << "PRODUCT_MIXED_WEIGHT_TEST - Fatal error!\n";
      std::cerr << "  Unexpected value of RULE[" << dim << "] = " << rule[dim] << ".\n";
      exit ( 1 );
    }
  }

  std::cout << "\n";
  std::cout << "PRODUCT_MIXED_WEIGHT_TEST:\n";
  std::cout << "  Compute the weights of a mixed factor product grid.\n";
  std::cout << "\n";
  std::cout << "  As a simple test, sum these weights.\n";
  std::cout << "  They should sum to exactly " << weight_sum_exact << "\n";
  std::cout << "\n";
  std::cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  std::cout << "\n";
  std::cout << " Dimension      Rule     Order        Alpha          Beta\n";
  std::cout << "\n";

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 2 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 3 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 4 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 5 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 6 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim]
           << "  " << std::setw(14) << alpha[dim] << "\n";
    }
    else if ( rule[dim] == 7 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 8 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim]
           << "  " << std::setw(14) << alpha[dim] << "\n";
    }
    else if ( rule[dim] == 9 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim]
           << "  " << std::setw(14) << alpha[dim]
           << "  " << std::setw(14) << beta[dim] << "\n";
    }
    else if ( rule[dim] == 10 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 11 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 12 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 13 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 14 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 15 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 16 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 17 )
    {
      std::cout << "  " << std::setw(8) << dim
           << "  " << std::setw(8) << rule[dim]
           << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else
    {
      std::cout << "\n";
      std::cout << "PRODUCT_MIXED_WEIGHT_TEST - Fatal error!\n";
      std::cout << "  Unexpected value of RULE = " << rule[dim] << "\n";
      exit ( 1 );
    }
  }
//
//  Compute the weights and points.
//
  weight = new double[order_nd];

  webbur::product_mixed_weight ( dim_num, order_1d, order_nd, rule, alpha, beta,
    weight );
//
//  Sum the weights.
//
  weight_sum = webbur::r8vec_sum ( order_nd, weight );

  weight_sum_error = webbur::r8_abs ( weight_sum - weight_sum_exact );

  std::cout << "\n";
  std::cout << "    Weight sum  Expected sum    Difference\n";
  std::cout << "\n";
  std::cout << "  " << std::setw(14) << weight_sum
       << "  " << std::setw(14) << weight_sum_exact
       << "  " << std::setw(14) << weight_sum_error << "\n";

  delete [] weight;

  return;
}
