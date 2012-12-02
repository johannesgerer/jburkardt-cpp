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

void sandia_sgmga_weight_tests ( );

void sandia_sgmga_weight_test 
(
  int dim_num,
  double importance[], 
  double level_weight[],
  int level_max_min,
  int level_max_max,
  int growth, 
  void ( *gw_compute_points[] ) ( int order, int dim, double w[] ),
  void ( *gw_compute_weights[] ) ( int order, int dim, double w[] ),
  double tol, 
  int ( *gw_compute_order[] ) ( int level, int growth ) 
);

typedef void ( *GWPointer ) ( int order, int dim, double w[] );
typedef int ( *GWPointer2 ) ( int level, int growth );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SANDIA_SGMGA_WEIGHT_PRB.
//
//  Discussion:
//
//    SANDIA_SGMGA_WEIGHT_PRB tests the SANDIA_SGMGA_WEIGHT function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_WEIGHT_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the SANDIA_SGMGA_WEIGHT function.\n";
//
//  Generate the weights for sparse grid rules.
//
  sandia_sgmga_weight_tests ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_WEIGHT_PRB\n";
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

void sandia_sgmga_weight_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_WEIGHT_TESTS calls SANDIA_SGMGA_WEIGHT_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Local Parameters:
//
//    Local, double TOL, a tolerance for point equality.
//    A value of sqrt ( eps ) is reasonable, and will allow the code to
//    consolidate points which are equal, or very nearly so.  A value of
//    -1.0, on the other hand, will force the code to use every point, 
//    regardless of duplication.
//
{
  int dim;
  int dim_num;
  int growth;
  GWPointer2 *gw_compute_order;
  GWPointer *gw_compute_points;
  GWPointer *gw_compute_weights;
  double *importance;
  int level_max_max;
  int level_max_min;
  double *level_weight;
  int np_sum;
  int *order_1d;
  int order_nd;
  double tol;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_WEIGHT_TESTS\n";
  std::cout << "  Call SANDIA_SGMGA_WEIGHT_TEST with various arguments.\n";
//
//  Set the point equality tolerance.
//
  tol = std::sqrt ( webbur::r8_epsilon ( ) );
  std::cout << "\n";
  std::cout << "  All tests will use a point equality tolerance of " << tol << "\n";
//
//  Test #1
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::clenshaw_curtis_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::clenshaw_curtis_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #2
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::clenshaw_curtis_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::clenshaw_curtis_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #3
//
  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  NP[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::clenshaw_curtis_points;
  gw_compute_points[2] = webbur::clenshaw_curtis_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[2] = webbur::clenshaw_curtis_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  gw_compute_order[2] = webbur::level_to_order_exp_cc;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #4
//
  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  NP[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::clenshaw_curtis_points;
  gw_compute_points[2] = webbur::clenshaw_curtis_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[2] = webbur::clenshaw_curtis_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  gw_compute_order[2] = webbur::level_to_order_exp_cc;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #5
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::patterson_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::patterson_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_gp;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #6
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::legendre_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::legendre_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_wn;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #7
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::laguerre_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::laguerre_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_nn;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #8
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 1;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 1.5;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::gen_laguerre_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::gen_laguerre_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_nn;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #9
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 2;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 0.5;
  P[1] = 1.5;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::fejer2_points;
  gw_compute_points[1] = webbur::jacobi_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::fejer2_weights;
  gw_compute_weights[1] = webbur::jacobi_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_f2;
  gw_compute_order[1] = webbur::level_to_order_linear_nn;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #10
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  NP = new int[dim_num];
  NP[0] = 1;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 2.0;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::gen_hermite_points;
  gw_compute_points[1] = webbur::hermite_genz_keister_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::gen_hermite_weights;
  gw_compute_weights[1] = webbur::hermite_genz_keister_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_linear_wn;
  gw_compute_order[1] = webbur::level_to_order_exp_hgk;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  LEVEL_MAX 0 to 5
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 5;
  growth = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::clenshaw_curtis_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::clenshaw_curtis_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #11
//
  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  NP[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::legendre_points;
  gw_compute_points[2] = webbur::hermite_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::legendre_weights;
  gw_compute_weights[2] = webbur::hermite_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_wn;
  gw_compute_order[2] = webbur::level_to_order_linear_wn;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Test #12
//
  dim_num = 3;
  importance = new double[dim_num];
  importance[0] = 1.0;
  importance[1] = 0.0;
  importance[2] = 1.0;
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 3;
  growth = 2;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  NP[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::clenshaw_curtis_points;
  gw_compute_points[2] = webbur::clenshaw_curtis_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[1] = webbur::clenshaw_curtis_weights;
  gw_compute_weights[2] = webbur::clenshaw_curtis_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  gw_compute_order[2] = webbur::level_to_order_exp_cc;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  A 2D case, with a rule of type 11 in dimension 1, and type 12 in dimension 2.
//  Full exponential growth is requested, so:
//    Dimension 1: 1, 3, 7
//    Dimension 2: 1, 3, 5
//
  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 3;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  growth = 2;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::nco_points;
  gw_compute_points[1] = webbur::ncc_points;
  gw_compute_weights = new GWPointer[dim_num];
  gw_compute_weights[0] = webbur::nco_weights;
  gw_compute_weights[1] = webbur::ncc_weights;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_f2;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  sandia_sgmga_weight_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, gw_compute_weights, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] gw_compute_weights;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

  return;
}
//***************************************************************************80

void sandia_sgmga_weight_test 
(
  int dim_num,
  double importance[],
  double level_weight[],
  int level_max_min,
  int level_max_max,
  int growth,
  void ( *gw_compute_points[] ) ( int order, int dim, double w[] ),
  void ( *gw_compute_weights[] ) ( int order, int dim, double w[] ),
  double tol, 
  int ( *gw_compute_order[] ) ( int level, int growth )
)

//***************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_WEIGHT_TEST checks the sum of the quadrature weights.
//
//  Discussion:
//
//    If any component rule is of Golub-Welsch type, we cannot compute
//    the exact weight sum, which we set, instead, to zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double IMPORTANCE[DIM_NUM], the anisotropic importance.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the anisotropic weights.
//
//    Input, int LEVEL_MAX_MIN, LEVEL_MAX_MAX, the minimum and
//    maximum values of LEVEL_MAX.
//
//    Input, int GROWTH, the growth rule. 
//    0, slow;
//    1, moderate;
//    2, full.
//
//    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int dim, double x[] ),
//    an array of pointers to functions which return the 1D quadrature points 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, void ( *GW_COMPUTE_WEIGHTS[] ) ( int order, int dim, double w[] ),
//    an array of pointers to functions which return the 1D quadrature weights 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, double TOL, a tolerance for point equality.
//
//    Input, int ( *GW_COMPUTE_ORDER[] ) ( int level, int growth ),
//    an array of pointers to functions which return the order of the
//    1D quadrature rule of a given level and growth rule.
//
{
  int dim;
  int i;
  int level_max;
  int point;
  int point_num;
  int point_total_num;
  int *sparse_unique_index;
  double *sparse_weight;
  double weight_sum;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_WEIGHT_TEST\n";
  std::cout << "  Compute the weights of a sparse grid.\n";
  std::cout << "\n";
  std::cout << "  Each sparse grid is of spatial dimension DIM_NUM,\n";
  std::cout << "  and is made up of product grids of levels up to LEVEL_MAX.\n";
  std::cout << "\n";
  std::cout << "  IMPORTANCE:  ";
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(14) << importance[dim];
  }
  std::cout << "\n";
  std::cout << "  LEVEL_WEIGHT:";
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(14) << level_weight[dim];
  }

  std::cout << "\n";
  std::cout << "     Level      Weight sum\n";
  std::cout << "\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    point_total_num = webbur::sandia_sgmga_size_total ( dim_num, level_weight,
      level_max, growth, gw_compute_order );

    point_num = webbur::sandia_sgmga_size ( dim_num, level_weight, level_max, 
      gw_compute_points, tol, growth, gw_compute_order );

    sparse_unique_index = new int[point_total_num];

    webbur::sandia_sgmga_unique_index ( dim_num, level_weight, level_max, 
      gw_compute_points, tol, point_num, point_total_num, 
      growth, gw_compute_order, sparse_unique_index );

    sparse_weight = new double[point_num];

    webbur::sandia_sgmga_weight ( dim_num, level_weight, level_max, 
      gw_compute_weights, point_num, point_total_num, sparse_unique_index, 
      growth, gw_compute_order, sparse_weight );

    weight_sum = webbur::r8vec_sum ( point_num, sparse_weight );

    std::cout << "  " << std::setw(8)  << level_max
              << "  " << std::setw(14) << weight_sum << "\n";

    delete [] sparse_unique_index;
    delete [] sparse_weight;
  }
  return;
}
