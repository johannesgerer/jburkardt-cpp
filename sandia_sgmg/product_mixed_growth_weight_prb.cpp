# include "sandia_rules.hpp"
# include "sandia_rules2.hpp"
# include "sandia_sgmg.hpp"

# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>
//
//  Two global variables needed to support the "parameter" function.
//
double *P;
int *NP;

int main ( );

void product_mixed_growth_weight_tests ( );

void product_mixed_growth_weight_test 
( 
  int dim_num,
  int order_1d[], 
  int order_nd,
  void ( *gw_compute_weights[] ) ( int order, int dim, double w[] )
);

typedef void ( *GWPointer ) ( int order, int dim, double w[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PRODUCT_MIXED_GROWTH_WEIGHT_PRB.
//
//  Discussion:
//
//    PRODUCT_MIXED_GROWTH_WEIGHT tests PRODUCT_MIXED_GROWTH_WEIGHT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2012
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
  std::cout << "PRODUCT_MIXED_GROWTH_WEIGHT_PRB\n";
  std::cout << "  C++ version\n";
//
//  Make sure the individual product rule weights are computed correctly.
//
  product_mixed_growth_weight_tests ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "PRODUCT_MIXED_GROWTH_WEIGHT\n";
  std::cout << "  Normal end of execution.\n";

  std::cout << "\n";
  webbur::timestamp ( );
  
  return 0;
}
namespace webbur 
{
//****************************************************************************80

double parameter ( int dim, int offset )

//****************************************************************************80
//
//  Purpose:
//
//    PARAMETER is a user-supplied routine to retrieve parameters.
//
//  Discussion:
//
//    The declaration for this function is in SANDIA_RULES.H
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM, the spatial dimension.
//
//    Input, int OFFSET, the offset of the parameter within the 
//    spatial dimension.
//
//    Output, double PARAMETER, the value of the OFFSET-th parameter
//    associated with the DIM-th dimension.
//
{
  int i;
  int j;
  double value;

  j = 0;
  for ( i = 0; i < dim; i++ )
  {
    j = j + NP[i];
  }
  value = P[j+offset];

  return value;
}
}
//****************************************************************************80

void product_mixed_growth_weight_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    PRODUCT_MIXED_GROWTH_WEIGHT_TESTS calls PRODUCT_MIXED_GROWTH_WEIGHT_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  GWPointer *gw_compute_weights;
  int np_sum;
  int *order_1d;
  int order_nd;

  std::cout << "\n";
  std::cout << "PRODUCT_MIXED_GROWTH_WEIGHT_TESTS\n";
  std::cout << "  Call PRODUCT_MIXED_GROWTH_WEIGHT_TEST with various arguments.\n";

  dim_num = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  order_1d = new int[dim_num];
  order_1d[0] = 3;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::clenshaw_curtis_weights;
  product_mixed_growth_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  order_1d = new int[dim_num];
  order_1d[0] = 3;
  order_1d[1] = 7;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::hermite_weights;
  product_mixed_growth_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  order_1d = new int[dim_num];
  order_1d[0] = 3;
  order_1d[1] = 3;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::patterson_weights;
  gw_compute_weights[1] = webbur::laguerre_weights;
  product_mixed_growth_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 1;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 1.5;
  order_1d = new int[dim_num];
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::gen_laguerre_weights;
  product_mixed_growth_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 0.5;
  P[1] = 1.5;
  order_1d = new int[dim_num];
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::fejer2_weights;
  gw_compute_weights[1] = webbur::jacobi_weights;
  product_mixed_growth_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 2.0;
  order_1d = new int[dim_num];
  order_1d[0] = 7;
  order_1d[1] = 9;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::gen_hermite_weights;
  gw_compute_weights[1] = webbur::hermite_genz_keister_weights;
  product_mixed_growth_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 3;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  NP[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  order_1d = new int[dim_num];
  order_1d[0] = 2;
  order_1d[1] = 3;
  order_1d[2] = 3;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::legendre_weights;
  gw_compute_weights[2] = webbur::hermite_weights;
  product_mixed_growth_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  return;
}
//***************************************************************************80

void product_mixed_growth_weight_test 
( 
  int dim_num,
  int order_1d[], 
  int order_nd, 
  void ( *gw_compute_weights[] ) ( int order, int dim, double w[] )
)

//***************************************************************************80
//
//  Purpose:
//
//    PRODUCT_MIXED_GROWTH_WEIGHT_TEST: weights of a mixed factor product rule.
//
//  Discussion:
//
//    This routine computes a sparse grid and compares the sum of the weights
//    to the expected exact value.
//
//    The routine cannot produce a result for rules that include one or more
//    component rules of type 11 or 12, that is, User Supplied rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2012
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
//    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int dim, double w[] ),
//    an array of pointers to functions which return the 1D quadrature weights 
//    associated with each spatial dimension for which a User Supplied rule 
//    is used.
//
{
  int dim;
  double *weight;
  double weight_sum;
//
//  Determine the integral of 1 over the multidimensional weighted region.
//
  std::cout << "\n";
  std::cout << "PRODUCT_MIXED_GROWTH_WEIGHT_TEST:\n";
  std::cout << "  Compute the weights of a mixed factor product grid.\n";
  std::cout << "\n";
  std::cout << "  As a simple test, sum these weights.\n";
  std::cout << "\n";
  std::cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  std::cout << "\n";
  std::cout << " Dimension     Order\n";
  std::cout << "\n";

  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(8) << dim
              << "  " << std::setw(8) << order_1d[dim] << "\n";
  }
//
//  Compute the weights and points of the sparse grid mixed rule.
//
  weight = new double[order_nd];

  webbur::product_mixed_growth_weight ( dim_num, order_1d, order_nd,
    gw_compute_weights, weight );
//
//  Sum the weights to get the SGM approximation to the integral of 1.
//
  weight_sum = webbur::r8vec_sum ( order_nd, weight );

  std::cout << "\n";
  std::cout << "    Weight sum\n";
  std::cout << "\n";
  std::cout << "  " << std::setw(14) << weight_sum << "\n";

  delete [] weight;

  return;
}
