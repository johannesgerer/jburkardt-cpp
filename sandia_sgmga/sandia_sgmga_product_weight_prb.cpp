# include "sandia_rules.hpp"
# include "sandia_rules2.hpp"
# include "sandia_sgmga.hpp"

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

void sandia_sgmga_product_weight_tests ( );

void sandia_sgmga_product_weight_test 
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
//    MAIN is the main program for SANDIA_SGMGA_PRODUCT_WEIGHT_PRB.
//
//  Discussion:
//
//    SANDIA_SGMGA_PRODUCT_WEIGHT_PRB tests the SANDIA_SGMGA_PRODUCT_WEIGHT function.
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
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_PRODUCT_WEIGHT_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the SANDIA_SGMGA_PRODUCT_WEIGHT function.\n";
//
//  Make sure the individual product rule weights are computed correctly.
//
  sandia_sgmga_product_weight_tests ( );
//
//  That's all.
//
  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_PRODUCT_WEIGHT_PRB\n";
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

void sandia_sgmga_product_weight_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_PRODUCT_WEIGHT_TESTS calls SANDIA_SGMGA_PRODUCT_WEIGHT_TEST.
//
//  Discussion:
//
//    To test user supplied rules for a spatial dimension DIM, we can
//    set RULE[DIM] = 11 or 12, and set the corresponding entry of 
//    GW_COMPUTE_WEIGHTS to the name of a function that we know is already
//    available, such as "webbur::clenshaw_curtis_compute_weights".
//
//    Note that, for ALL the tests, we set every entry of the GW_COMPUTE_WEIGHTS
//    array.  However, a particular entry is only inspected if the corresponding
//    entry of RULE is 11 or 12.
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
  std::cout << "SANDIA_SGMGA_PRODUCT_WEIGHT_TESTS\n";
  std::cout << "  Call SANDIA_SGMGA_PRODUCT_WEIGHT_TEST with various arguments.\n";

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 3;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::clenshaw_curtis_weights;
  sandia_sgmga_product_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 3;
  order_1d[1] = 7;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::hermite_weights;
  sandia_sgmga_product_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 3;
  order_1d[1] = 3;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::patterson_weights;
  gw_compute_weights[1] = webbur::laguerre_weights;
  sandia_sgmga_product_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 1;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 1.5;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::gen_laguerre_weights;
  sandia_sgmga_product_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 5;
  order_1d[1] = 5;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 0.5;
  P[1] = 1.5;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::fejer2_weights;
  gw_compute_weights[1] = webbur::jacobi_weights;
  sandia_sgmga_product_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;
//
//  GEN HERMITE x HERMITE GENZ-KEISTER
//
  dim_num = 2;
  order_1d = new int[dim_num];
  order_1d[0] = 7;
  order_1d[1] = 9;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 2.0;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::gen_hermite_weights;
  gw_compute_weights[1] = webbur::hermite_genz_keister_weights;
  sandia_sgmga_product_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  dim_num = 3;
  order_1d = new int[dim_num];
  order_1d[0] = 2;
  order_1d[1] = 3;
  order_1d[2] = 3;
  order_nd = webbur::i4vec_product ( dim_num, order_1d );
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  NP[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::legendre_weights;
  gw_compute_weights[2] = webbur::hermite_weights;
  sandia_sgmga_product_weight_test ( dim_num, order_1d, order_nd, gw_compute_weights );
  delete [] gw_compute_weights;
  delete [] NP;
  delete [] order_1d;
  delete [] P;

  return;
}
//***************************************************************************80

void sandia_sgmga_product_weight_test 
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
//    SANDIA_SGMGA_PRODUCT_WEIGHT_TEST: weights of a mixed factor product rule.
//
//  Discussion:
//
//    This routine computes a sparse grid and compares the sum of the weights
//    to the expected exact value.
//
//    The routine cannot produce a result for rules that include one or more
//    component rules of type 11 or 12, that is, user supplied rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 November 2011
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
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
{
  int dim;
  double *weight;
  double weight_sum;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_PRODUCT_WEIGHT_TEST:\n";
  std::cout << "  Compute the weights of a mixed factor product grid.\n";
  std::cout << "\n";
  std::cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";

  std::cout << "\n";
  std::cout << " Dimension      Order\n";
  std::cout << "\n";

  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(8) << dim
              << "  " << std::setw(8) << order_1d[dim] << "\n";
  }
//
//  Compute the weights.
//
  weight = new double[order_nd];

  webbur::sandia_sgmga_product_weight ( dim_num, order_1d, order_nd, 
    gw_compute_weights, weight );
//
//  Sum the weights to get the approximation to the integral of 1.
//
  weight_sum = webbur::r8vec_sum ( order_nd, weight );

  std::cout << "\n";
  std::cout << "    Weight sum\n";
  std::cout << "\n";
  std::cout << "  " << std::setw(14) << weight_sum << "\n";

  delete [] weight;

  return;
}

