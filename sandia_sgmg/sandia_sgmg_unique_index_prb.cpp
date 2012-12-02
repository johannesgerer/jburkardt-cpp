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

void sgmg_unique_index_tests ( double tol );

void sgmg_unique_index_test 
(
  int dim_num,
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
//    MAIN tests SGMG_UNIQUE_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2012
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
  std::cout << "SGMG_UNIQUE_INDEX_PRB\n";
  std::cout << "  C++ version\n";
//
//  Check that we can determine a UNIQUE_INDEX mapping, so that we can
//  generate all the points, and then select unique representatives, and
//  then match each point to its representative.
//
  tol = std::sqrt ( webbur::r8_epsilon ( ) );
  sgmg_unique_index_tests ( tol );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SGMG_UNIQUE_INDEX_PRB\n";
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

void sgmg_unique_index_tests ( double tol )

//****************************************************************************80
//
//  Purpose:
//
//    SGMG_UNIQUE_INDEX_TESTS calls SGMG_UNIQUE_INDEX_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2011
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
//    -1.0, on the other hand, will force the code to use every point, 
//    regardless of duplication.
//
{
  int dim_num;
  int growth;
  GWPointer2 *gw_compute_order;
  GWPointer *gw_compute_points;
  int level_max_max;
  int level_max_min;
  int np_sum;

  std::cout << "\n";
  std::cout << "SGMG_UNIQUE_INDEX_TESTS\n";
  std::cout << "  Call SGMG_UNIQUE_INDEX_TEST with various arguments.\n";
  std::cout << "\n";
  std::cout << "  All tests will use a point equality tolerance of " << tol << "\n";

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  growth = 2;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::clenshaw_curtis_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_cc;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  growth = 2;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::patterson_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_exp_gp;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::legendre_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_wn;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::laguerre_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_nn;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::gen_laguerre_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_nn;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 1;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 1.5;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::fejer2_points;
  gw_compute_points[1] = webbur::jacobi_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_f2;
  gw_compute_order[1] = webbur::level_to_order_linear_nn;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 2;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 0.5;
  P[1] = 1.5;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::gen_hermite_points;
  gw_compute_points[1] = webbur::hermite_genz_keister_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_linear_wn;
  gw_compute_order[1] = webbur::level_to_order_exp_hgk;
  NP = new int[dim_num];
  NP[0] = 1;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  P[0] = 2.0;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;

  dim_num = 3;
  level_max_min = 0;
  level_max_max = 2;
  growth = 1;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_points;
  gw_compute_points[1] = webbur::legendre_points;
  gw_compute_points[2] = webbur::hermite_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_cc;
  gw_compute_order[1] = webbur::level_to_order_linear_wn;
  gw_compute_order[2] = webbur::level_to_order_linear_wn;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  NP[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;
//
//  Dimension 2, Level 4, Rule 3, slow exponential growth.
//
  dim_num = 2;
  level_max_min = 4;
  level_max_max = 4;
  growth = 0;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::patterson_points;
  gw_compute_points[1] = webbur::patterson_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_gp;
  gw_compute_order[1] = webbur::level_to_order_exp_gp;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;
//
//  Dimension 2, Level 4, Rule 3, moderate exponential growth.
//
  dim_num = 2;
  level_max_min = 4;
  level_max_max = 4;
  growth = 1;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::patterson_points;
  gw_compute_points[1] = webbur::patterson_points;
  gw_compute_order = new GWPointer2[dim_num];
  gw_compute_order[0] = webbur::level_to_order_exp_gp;
  gw_compute_order[1] = webbur::level_to_order_exp_gp;
  NP = new int[dim_num];
  NP[0] = 0;
  NP[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, NP );
  P = new double[np_sum];
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, growth, 
    gw_compute_points, tol, gw_compute_order );
  delete [] gw_compute_order;
  delete [] gw_compute_points;
  delete [] NP;
  delete [] P;

  return;
}
//***************************************************************************80

void sgmg_unique_index_test 
(
  int dim_num,
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
//    SGMG_UNIQUE_INDEX_TEST tests SGMG_UNIQUE_INDEX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2011
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
//    Input, int GROWTH, the growth rule. 
//    0, slow;
//    1, moderate;
//    2, full.
//
//    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int dim, double x[] ),
//    an array of pointers to functions which return the 1D quadrature points 
//    associated with each spatial dimension.
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
  int level_min;
  int point;
  int point_num;
  int point_total_num;
  int *sparse_unique_index;

  std::cout << "\n";
  std::cout << "SGMG_UNIQUE_INDEX_TEST\n";
  std::cout << "  SGMG_UNIQUE_INDEX returns a mapping between\n";
  std::cout << "  the nonunique and unique points in a sparse grid.\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    point_total_num = webbur::sgmg_size_total ( dim_num, level_max, growth,
      gw_compute_order );

    point_num = webbur::sgmg_size ( dim_num, level_max, gw_compute_points, tol, 
      growth, gw_compute_order );

    std::cout << "\n";
    std::cout << " LEVEL_MIN LEVEL_MAX POINT_NUM POINT_NUM\n";
    std::cout << "                        Unique     Total\n";

    level_min = webbur::i4_max ( 0, level_max + 1 - dim_num );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << level_min
              << "  " << std::setw(8) << level_max
              << "  " << std::setw(8) << point_num
              << "  " << std::setw(8) << point_total_num << "\n";

    sparse_unique_index = new int[point_total_num];

    webbur::sgmg_unique_index ( dim_num, level_max, gw_compute_points, tol, 
      point_num, point_total_num, growth, gw_compute_order, sparse_unique_index );

    std::cout << "\n";
    std::cout << "     POINT    UNIQUE\n";
    std::cout << "\n";
    for ( point = 0; point < point_total_num; point++ )
    {
      std::cout << "  " << std::setw(8) << point
                << "  " << std::setw(8) << sparse_unique_index[point] << "\n";
    }
    delete [] sparse_unique_index;
  }
  return;
}

