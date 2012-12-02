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

void sandia_sgmga_size_tabulate 
( 
  int gw_compute_order_1d ( int level, int growth ),
  int growth_1d,
  int pd_1d, 
  double p_1d[],
  int dim_min,
  int dim_max,
  int level_max_min,
  int level_max_max, 
  void gw_compute_points_1d ( int order, int dim, double x[] )
 );

typedef void ( *GWPointer ) ( int order, int dim, double w[] );
typedef int ( *GWPointer2 ) ( int level, int growth );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SANDIA_SGMGA_SIZE_TABLE.
//
//  Discussion:
//
//    SANDIA_SGMGA_INDEX_PRB tests the SANDIA_SGMGA_INDEX function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  int dim_max;
  int dim_min;
  int growth_1d;
  int level_max_max;
  int level_max_min;
  int np_1d;
  double *p_1d;

  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_SIZE_TABLE\n";
  std::cout << "  C++ version\n";
  std::cout << "  Make tables of point counts.\n";
//
//  Clenshaw-Curtis Grid, slow growth.
//
  growth_1d = 0;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::clenshaw_curtis_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Clenshaw-Curtis Grid, slow growth.
//
  growth_1d = 0;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::clenshaw_curtis_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Clenshaw-Curtis Grid, full growth.
//
  growth_1d = 2;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::clenshaw_curtis_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Clenshaw-Curtis Grid, full growth.
//
  growth_1d = 2;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::clenshaw_curtis_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Clenshaw-Curtis Grid, exponential growth.
//
  growth_1d = 2;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 100;
  dim_max = 100;
  level_max_min = 0;
  level_max_max = 2;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::clenshaw_curtis_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss-Patterson Grid, Slow growth.
//
  growth_1d = 0;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_exp_gp, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::patterson_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss-Patterson Grid, Moderate growth.
//
  growth_1d = 1;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_exp_gp, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::patterson_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss-Patterson Grid, Full growth.
//
  growth_1d = 2;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 6;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_exp_gp, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::patterson_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss Legendre Grid, Slow growth.
//
  growth_1d = 0;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_linear_wn, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::legendre_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss Legendre Grid, Moderate growth.
//
  growth_1d = 1;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_linear_wn, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::legendre_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss Laguerre Grid, Moderate growth.
//
  growth_1d = 1;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_linear_nn, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::laguerre_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Hermite Genz Keister, Slow growth.
//
  growth_1d = 0;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sandia_sgmga_size_tabulate ( webbur::level_to_order_exp_hgk, growth_1d, np_1d, 
    p_1d, dim_min, dim_max, level_max_min, level_max_max, webbur::hermite_genz_keister_points );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_SIZE_TABLE\n";
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

void sandia_sgmga_size_tabulate 
(
  int gw_compute_order_1d ( int level, int growth ), 
  int growth_1d,
  int np_1d, 
  double p_1d[],
  int dim_min,
  int dim_max,
  int level_max_min,
  int level_max_max,
  void gw_compute_points_1d ( int order, int dim, double x[] ) 
)

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_SIZE_TABULATE tests SANDIA_SGMGA_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int GW_COMPUTE_ORDER_1D ( int level, int growth ),
//    a function which returns the order of the
//    1D quadrature rule of a given level and growth rule.
//
//    Input, int GROWTH, the growth rule. 
//    0, slow;
//    1, moderate;
//    2, full.
//
//    Input, int NP_1D, the number of parameters in the 1D rule.
//
//    Input, double P_1D[NP_1D], the parameters.
//
//    Input, int DIM_MIN, the minimum spatial dimension to consider.
//
//    Input, int DIM_MAX, the maximum spatial dimension to consider.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX to consider.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX to consider.
//
//    Input, GW_COMPUTE_POINTS_1D ( int order, int dim, double x[] ),
//    a function which return the 1D quadrature points.
//
{
  int dim;
  int dim_num;
  int growth;
  GWPointer2 *gw_compute_order;
  GWPointer *gw_compute_points;
  int i;
  int level_max;
  double *level_weight;
  int np_sum;
  int point_num;
  double tol;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_SIZE_TABULATE\n";
  std::cout << "  SANDIA_SGMGA_SIZE returns the number of distinct\n";
  std::cout << "  points in a sparse grid.\n";
  std::cout << "\n";
  std::cout << "  We use the same rule in all dimensions, and count the points,\n";
  std::cout << "  for a range of dimensions and levels.\n";
  std::cout << "\n";
  std::cout << "  Growth rule = " << growth_1d << "\n";
  std::cout << "\n";

  tol = std::sqrt ( webbur::r8_epsilon ( ) );

  std::cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++)
  {
    std::cout << "  " << std::setw(8) << dim_num;
  }
  std::cout << "\n";
  std::cout << "\n";
  std::cout << "   LEVEL_MAX\n";
  std::cout << "\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    std::cout << "    " << std::setw(4) << level_max;
    for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
    {
      level_weight = new double[dim_num];
      growth = growth_1d;
      NP = new int[dim_num];
      np_sum = dim_num * np_1d;
      P = new double[np_sum];
      gw_compute_points = new GWPointer[dim_num];
      gw_compute_order = new GWPointer2[dim_num];

      for ( dim = 0; dim < dim_num; dim++ )
      {
        level_weight[dim] = 1.0;
        NP[dim] = np_1d;
        for ( i = 0; i < np_1d; i++ )
        {
          P[i+dim*np_1d] = p_1d[i];
        }
        gw_compute_points[dim] = gw_compute_points_1d;
        gw_compute_order[dim] = gw_compute_order_1d;
      }
      point_num = webbur::sandia_sgmga_size ( dim_num, level_weight, level_max, 
        gw_compute_points, tol, growth, gw_compute_order );

      std::cout << "  " << std::setw(8) << point_num;

      delete [] gw_compute_order;
      delete [] gw_compute_points;
      delete [] level_weight;
      delete [] NP;
      delete [] P;
    }
    std::cout << "\n";
  }

  return;
}
