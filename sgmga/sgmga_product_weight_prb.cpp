# include "sandia_rules.hpp"
# include "sgmga.hpp"

# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>

int main ( );

void sgmga_product_weight_tests ( );
void sgmga_product_weight_test ( int dim_num, int order_1d[], 
  int order_nd, int rule[], int np[], double p[],
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ) );

typedef void ( *GWPointer ) ( int order, int np, double p[], double w[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SGMGA_PRODUCT_WEIGHT_PRB.
//
//  Discussion:
//
//    SGMGA_PRODUCT_WEIGHT_PRB tests the SGMGA_PRODUCT_WEIGHT function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
{
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SGMGA_PRODUCT_WEIGHT_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the SGMGA_PRODUCT_WEIGHT function.\n";
//
//  Make sure the individual product rule weights are computed correctly.
//
  sgmga_product_weight_tests ( );
//
//  That's all.
//
  std::cout << "\n";
  std::cout << "SGMGA_PRODUCT_WEIGHT_PRB\n";
  std::cout << "  Normal end of execution.\n";

  std::cout << "\n";
  webbur::timestamp ( );
  
  return 0;
}
//****************************************************************************80

void sgmga_product_weight_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_PRODUCT_WEIGHT_TESTS calls SGMGA_PRODUCT_WEIGHT_TEST.
//
//  Discussion:
//
//    To test Golub Welsch rules for a spatial dimension DIM, we can
//    set RULE[DIM] = 10, and set the corresponding entry of 
//    GW_COMPUTE_WEIGHTS to the name of a function that we know is already
//    available, such as "webbur::clenshaw_curtis_compute_weights".
//
//    Note that, for ALL the tests, we set every entry of the GW_COMPUTE_WEIGHTS
//    array.  However, a particular entry is only inspected if the corresponding
//    entry of RULE is 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  GWPointer *gw_compute_weights;
  int *np;
  int np_sum;
  int *order_1d;
  int order_nd;
  double *p;
  int *rule;

  std::cout << "\n";
  std::cout << "SGMGA_PRODUCT_WEIGHT_TESTS\n";
  std::cout << "  Call SGMGA_PRODUCT_WEIGHT_TEST with various arguments.\n";

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 3;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 1;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = webbur::clenshaw_curtis_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] np;
  delete [] order_1d;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 3;
  order_1d[1] = 7;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 5;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = webbur::hermite_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] np;
  delete [] order_1d;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 3;
  order_1d[1] = 3;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule = new int[dim_num];
  rule[0] = 3;
  rule[1] = 7;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::patterson_lookup_weights_np;
  gw_compute_weights[1] = webbur::laguerre_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] np;
  delete [] order_1d;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 8;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 1;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  p[0] = 1.5;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = webbur::gen_laguerre_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] np;
  delete [] order_1d;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule = new int[dim_num];
  rule[0] = 2;
  rule[1] = 9;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  p[0] = 0.5;
  p[1] = 1.5;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::fejer2_compute_weights_np;
  gw_compute_weights[1] = webbur::jacobi_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] np;
  delete [] order_1d;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 7;
  order_1d[1] = 9;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule = new int[dim_num];
  rule[0] = 6;
  rule[1] = 10;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  p[0] = 2.0;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::gen_hermite_compute_weights_np;
  gw_compute_weights[1] = webbur::legendre_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] np;
  delete [] order_1d;
  delete [] p;
  delete [] rule;

  dim_num = 3;
  order_1d = new int[dim_num];
  order_1d[0] = 2;
  order_1d[1] = 3;
  order_1d[2] = 3;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 4;
  rule[2] = 5;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = webbur::legendre_compute_weights_np;
  gw_compute_weights[2] = webbur::hermite_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] np;
  delete [] order_1d;
  delete [] p;
  delete [] rule;
//
//  Repeat, treating  rules #2 and #3 as Golub Welsch rules.
//
  dim_num = 3;
  order_1d = new int[dim_num];
  order_1d[0] = 2;
  order_1d[1] = 3;
  order_1d[2] = 3;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 11;
  rule[2] = 11;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_compute_weights_np;
  gw_compute_weights[1] = webbur::legendre_compute_weights_np;
  gw_compute_weights[2] = webbur::hermite_compute_weights_np;
  sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] np;
  delete [] order_1d;
  delete [] p;
  delete [] rule;

  return;
}
//***************************************************************************80

void sgmga_product_weight_test ( int dim_num, int order_1d[], 
  int order_nd, int rule[], int np[], double p[],
  void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ) )

//***************************************************************************80
//
//  Purpose:
//
//    SGMGA_PRODUCT_WEIGHT_TEST: weights of a mixed factor product rule.
//
//  Discussion:
//
//    This routine computes a sparse grid and compares the sum of the weights
//    to the expected exact value.
//
//    The routine cannot produce a result for rules that include one or more
//    component rules of type 10, that is, Golub-Welsch rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2010
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
//    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
//    11, "UO",  User supplied Open, presumably Non Nested.
//    12, "UC",  User supplied Closed, presumably Non Nested.
//
//    Input, int NP[RULE_NUM], the number of parameters used by each rule.
//
//    Input, double P[sum(NP[*])], the parameters needed by each rule.
//
//    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int np, double p[], double w[] ),
//    an array of pointers to functions which return the 1D quadrature weights 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
{
  double alpha;
  double beta;
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  int dim;
  int i;
  int p_index;
  double pi = 3.141592653589793;
  double value1;
  double value2;
  double *weight;
  double weight_sum;
  double weight_sum_error;
  double weight_sum_exact;
//
//  Determine the integral of 1 over the multidimensional weighted region.
//
  p_index = 0;

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
      weight_sum_exact = weight_sum_exact * std::sqrt ( pi );
    }
    else if ( rule[dim] == 6 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;

      weight_sum_exact = weight_sum_exact 
        * webbur::r8_gamma ( 0.5 * ( alpha + 1.0 ) );
    }
    else if ( rule[dim] == 7 )
    {
      weight_sum_exact = weight_sum_exact * 1.0;
    }
    else if ( rule[dim] == 8 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;

      weight_sum_exact = weight_sum_exact * webbur::r8_gamma ( alpha + 1.0 );
    }
    else if ( rule[dim] == 9 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      beta = p[p_index];
      p_index = p_index + 1;
      arg1 = - alpha;
      arg2 = 1.0;
      arg3 = beta + 2.0;
      arg4 = - 1.0;
      value1 = webbur::r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );
      arg1 = - beta;
      arg2 = 1.0;
      arg3 = alpha + 2.0;
      arg4 = - 1.0;
      value2 = webbur::r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );
      weight_sum_exact = weight_sum_exact * ( 
        value1 / ( beta + 1.0 ) + value2 / ( alpha + 1.0 ) );
    }
    else if ( rule[dim] == 10 )
    {
      weight_sum_exact = weight_sum_exact * std::sqrt ( pi );
    }
    else if ( rule[dim] == 11 )
    {
      for ( i = 0; i < np[dim]; i++ )
      {
        alpha = p[p_index];
        p_index = p_index + 1;
      } 
      weight_sum_exact = 0.0;
    }
    else if ( rule[dim] == 12 )
    {
      for ( i = 0; i < np[dim]; i++ )
      {
        alpha = p[p_index];
        p_index = p_index + 1;
      } 
      weight_sum_exact = 0.0;
    }
    else
    {
      std::cerr << "\n";
      std::cerr << "SGMGA_PRODUCT_WEIGHT_TEST - Fatal error!\n";
      std::cerr << "  Unexpected value of RULE[" << dim << "] = " 
           << rule[dim] << ".\n";
      std::exit ( 1 );
    }
  }

  std::cout << "\n";
  std::cout << "SGMGA_PRODUCT_WEIGHT_TEST:\n";
  std::cout << "  Compute the weights of a mixed factor product grid.\n";
  if ( weight_sum_exact != 0.0 )
  {
    std::cout << "\n";
    std::cout << "  As a simple test, sum these weights.\n";
    std::cout << "  They should sum to exactly " << weight_sum_exact << "\n";
  }
  std::cout << "\n";
  std::cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  std::cout << "\n";
  std::cout << " Dimension      Rule    Growth        Parameters\n";
  std::cout << "\n";

  p_index = 0;
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
      alpha = p[p_index];
      p_index = p_index + 1;
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << order_1d[dim] 
                << "  " << std::setw(14) << alpha << "\n";
    }
    else if ( rule[dim] == 7 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << order_1d[dim] << "\n";
    }
    else if ( rule[dim] == 8 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << order_1d[dim] 
                << "  " << std::setw(14) << alpha << "\n";
    }
    else if ( rule[dim] == 9 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      beta = p[p_index];
      p_index = p_index + 1;
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << order_1d[dim] 
                << "  " << std::setw(14) << alpha 
                << "  " << std::setw(14) << beta << "\n";
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
                << "  " << std::setw(8) << order_1d[dim];
      for ( i = 0; i < np[dim]; i++ )
      {
        alpha = p[p_index];
        p_index = p_index + 1;
        std::cout << "  " << std::setw(14) << alpha;
      } 
      std::cout << "\n";
    }
    else if ( rule[dim] == 12 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << order_1d[dim];
      for ( i = 0; i < np[dim]; i++ )
      {
        alpha = p[p_index];
        p_index = p_index + 1;
        std::cout << "  " << std::setw(14) << alpha;
      } 
      std::cout << "\n";
    }
    else
    {
      std::cerr << "\n";
      std::cerr << "SGMGA_PRODUCT_WEIGHT_TEST - Fatal error!\n";
      std::cerr << "  Cannot perform test for rule = " << rule[dim] << "\n";
      std::exit ( 1 );
    }
  }
//
//  Compute the weights.
//
  weight = new double[order_nd];

  webbur::sgmga_product_weight ( dim_num, order_1d, order_nd, rule, 
    np, p, gw_compute_weights, weight );
//
//  Sum the weights to get the approximation to the integral of 1.
//
  weight_sum = webbur::r8vec_sum ( order_nd, weight );
//
//  Compare the exact and estimated integrals.
//
  weight_sum_error = webbur::r8_abs ( weight_sum - weight_sum_exact );

  if ( weight_sum_exact != 0.0 )
  {
    std::cout << "\n";
    std::cout << "    Weight sum  Expected sum    Difference\n";
    std::cout << "\n";
    std::cout << "  " << std::setw(14) << weight_sum
              << "  " << std::setw(14) << weight_sum_exact
              << "  " << std::setw(14) << weight_sum_error << "\n";
  }
  else
  {
    std::cout << "\n";
    std::cout << "    Weight sum\n";
    std::cout << "\n";
    std::cout << "  " << std::setw(14) << weight_sum << "\n";
  }
  delete [] weight;

  return;
}

