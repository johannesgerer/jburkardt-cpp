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

void sandia_sgmga_size_tests ( );

void sandia_sgmga_size_test 
( 
  int dim_num, 
  double importance[], 
  double level_weight[],
  int level_max_min, 
  int level_max_max, 
  int growth, 
  void ( *gw_compute_points[] ) ( int order, int dim, double w[] ),
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
//    MAIN is the main program for SANDIA_SGMGA_SIZE_PRB.
//
//  Discussion:
//
//    SANDIA_SGMGA_SIZE_PRB tests SANDIA_SGMGA_SIZE and SANDIA_SGMGA_SIZE_TOTAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2010
//
//  Author:
//
//    John Burkardt
//
{
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_SIZE_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test SANDIA_SGMGA_SIZE and SANDIA_SGMGA_SIZE_TOTAL.\n";

  sandia_sgmga_size_tests ( );

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_SIZE_PRB\n";
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

void sandia_sgmga_size_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_SIZE_TESTS calls SANDIA_SGMGA_SIZE_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2012
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
  double *importance;
  int level_max_max;
  int level_max_min;
  double *level_weight;
  int np_sum;
  int *order_1d;
  int order_nd;
  double tol;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_SIZE_TESTS\n";
  std::cout << "  Call SANDIA_SGMGA_SIZE_TEST with various arguments.\n";
//
//  Set the point equality tolerance.
//
  tol = std::sqrt ( webbur::r8_epsilon ( ) );
  std::cout << "\n";
  std::cout << "  Point equality tolerance = " << tol << "\n";

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 5;
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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  gw_compute_order[2] = webbur::level_to_order_exp_cc;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

  dim_num = 3;
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
  NP[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::clenshaw_curtis_points;
  gw_compute_points[2] = webbur::clenshaw_curtis_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  gw_compute_order[2] = webbur::level_to_order_exp_cc;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_gp;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_wn;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_nn;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_nn;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_f2;
  gw_compute_order[1] = webbur::level_to_order_linear_nn;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  GEN HERMITE x HERMITE GENZ-KEISTER
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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_linear_wn;
  gw_compute_order[1] = webbur::level_to_order_exp_hgk;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_wn;
  gw_compute_order[2] = webbur::level_to_order_linear_wn;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;
//
//  Try a case involving a dimension of "0" importance.
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
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  gw_compute_order[2] = webbur::level_to_order_exp_cc;
  sandia_sgmga_size_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max, growth, gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] NP;
  delete [] P;

  return;
}
//***************************************************************************80

void sandia_sgmga_size_test 
( 
  int dim_num,
  double importance[], 
  double level_weight[],
  int level_max_min,
  int level_max_max,
  int growth, 
  void ( *gw_compute_points[] ) ( int order, int dim, double w[] ),
  double tol,
  int ( *gw_compute_order[] ) ( int level, int growth )
)

//***************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_SIZE_TEST tests SGMG_SIZE, SGMG_SIZE_TOTAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double IMPORTANCE[DIM_NUM], the importance for each dimension.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the weights for each dimension.
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
//    Input, double TOL, a tolerance for point equality.
//
//    Input, int ( *GW_COMPUTE_ORDER[] ) ( int level, int growth ),
//    an array of pointers to functions which return the order of the
//    1D quadrature rule of a given level and growth rule.
//
{
  int dim;
  int level_max;
  int p_index;
  int point_num;
  int point_total_num;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_SIZE_TEST\n";
  std::cout << "  SANDIA_SGMGA_SIZE_TOTAL counts the total number of points,\n";
  std::cout << "  including duplications, in an SANDIA_SGMGA sparse grid.\n";
  std::cout << "  SANDIA_SGMGA_SIZE counts the total number of points,\n";
  std::cout << "  excluding duplications, in an SANDIA_SGMGA sparse grid.\n";
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
  std::cout << "   DIM_NUM LEVEL_MAX POINT_NUM POINT_NUM\n";
  std::cout << "                        Unique     Total\n";
  std::cout << "\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    point_total_num = webbur::sandia_sgmga_size_total ( dim_num, level_weight,
      level_max, growth, gw_compute_order );

    point_num = webbur::sandia_sgmga_size ( dim_num, level_weight, level_max, 
      gw_compute_points, tol, growth, gw_compute_order );

    std::cout << "  " << std::setw(8) << dim_num
              << "  " << std::setw(8) << level_max
              << "  " << std::setw(8) << point_num
              << "  " << std::setw(8) << point_total_num << "\n";
  }

  return;
}
