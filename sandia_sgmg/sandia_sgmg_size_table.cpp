# include "sandia_rules.hpp"
# include "sandia_rules2.hpp"
# include "sandia_sgmg.hpp"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
//
//  Two global variables needed to support the "parameter" function.
//
double *P;
int *NP;

int main ( );
void sgmg_size_tabulate 
( 
  int gw_compute_order_1d ( int level, int growth ), 
  int growth_1d,
  int np_1d,
  double p_1d[], 
  void gw_compute_points_1d ( int order, int dim, double w[] ),
  int dim_min,
  int dim_max,
  int level_max_min,
  int level_max_max 
);

typedef void ( *GWPointer ) ( int order, int dim, double w[] );
typedef int ( *GWPointer2 ) ( int level, int growth );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN makes a point count table.
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
  std::cout << "SGMG_SIZE_TABLE\n";
  std::cout << "  C++ version\n";
//
//  Clenshaw-Curtis, Slow Exponential Growth.
//
  if ( true )
  {
    growth_1d = 0;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 1;
    dim_max = 5;
    level_max_min = 0;
    level_max_max = 7;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
      p_1d, webbur::clenshaw_curtis_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Clenshaw-Curtis, Slow Exponential Growth.
//
  if ( true )
  {
    growth_1d = 0;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 5;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
      p_1d, webbur::clenshaw_curtis_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Clenshaw-Curtis Grid, Full exponential growth.
//
  if ( true )
  {
    growth_1d = 2;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 1;
    dim_max = 5;
    level_max_min = 0;
    level_max_max = 7;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
      p_1d, webbur::clenshaw_curtis_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Clenshaw-Curtis Grid, Full exponential growth.
//
  if ( true )
  {
    growth_1d = 2;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 5;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
      p_1d, webbur::clenshaw_curtis_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Clenshaw-Curtis Grid, Full exponential growth.
//
  if ( true )
  {
    growth_1d = 2;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 100;
    dim_max = 100;
    level_max_min = 0;
    level_max_max = 2;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_cc, growth_1d, np_1d, 
      p_1d, webbur::clenshaw_curtis_points, 
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Gauss-Patterson, Slow Exponential Growth.
//
  if ( true )
  {
    growth_1d = 0;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 1;
    dim_max = 5;
    level_max_min = 0;
    level_max_max = 7;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_gp, growth_1d, np_1d, 
      p_1d, webbur::patterson_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Gauss-Patterson, Slow Exponential Growth.
//
  if ( true )
  {;
    growth_1d = 0;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 5;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_gp, growth_1d, np_1d, 
      p_1d, webbur::patterson_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Gauss-Patterson, Moderate Exponential Growth.
//
  if ( true )
  {
    growth_1d = 1;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 1;
    dim_max = 5;
    level_max_min = 0;
    level_max_max = 7;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_gp, growth_1d, np_1d, 
      p_1d, webbur::patterson_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Gauss-Patterson, Moderate Exponential Growth.
//
  if ( true )
  {
    growth_1d = 1;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 5;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_gp, growth_1d, np_1d, 
      p_1d, webbur::patterson_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Gauss-Patterson, Full Exponential Growth
//
  if ( true )
  {
    growth_1d = 2;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 1;
    dim_max = 5;
    level_max_min = 0;
    level_max_max = 7;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_gp, growth_1d, np_1d, 
      p_1d, webbur::patterson_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Gauss-Patterson, Full Exponential Growth.
//
  if ( true )
  {
    growth_1d = 2;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 6;
    dim_max = 10;
    level_max_min = 0;
    level_max_max = 5;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_exp_gp, growth_1d, np_1d, 
      p_1d, webbur::patterson_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Gauss-Legendre, slow linear odd growth.
//
  if ( true )
  {
    growth_1d = 0;
    np_1d = 0;
    p_1d = new double[np_1d];
    dim_min = 1;
    dim_max = 5;
    level_max_min = 0;
    level_max_max = 7;
    ctime = webbur::cpu_time ( );
    sgmg_size_tabulate ( webbur::level_to_order_linear_wn, growth_1d, np_1d, 
      p_1d, webbur::legendre_points,
      dim_min, dim_max, level_max_min, level_max_max );
    ctime = webbur::cpu_time ( ) - ctime;
    std::cout << "\n";
    std::cout << "  CPU Time = " << ctime << "\n";
    delete [] p_1d;
  }
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SGMG_SIZE_TABLE\n";
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
//***************************************************************************80

void sgmg_size_tabulate 
( 
  int gw_compute_order_1d ( int level, int growth ), 
  int growth_1d,
  int np_1d,
  double p_1d[], 
  void gw_compute_points_1d ( int order, int dim, double w[] ),
  int dim_min,
  int dim_max,
  int level_max_min,
  int level_max_max 
)

//***************************************************************************80
//
//  Purpose:
//
//    SGMG_SIZE_TABULATE tests SGMG_SIZE.
//
//  Discussion:
//
//    We do NOT consider mixed rules.  Instead, we are looking at sparse grid
//    rules for which all dimensions use the same 1D rule family.
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
//    Input, double P_1D[NP_1D], the parameters for the 1D rule.
//
//    Input, void GW_COMPUTE_POINTS_1D ( int order, int dim, double w[] ),
//    the function to be used to compute points for the 1D rule.
//
//    Input, int DIM_MIN, the minimum spatial dimension.
//
//    Input, int DIM_MAX, the maximum spatial dimension.
//
//    Input, int LEVEL_MAX_MIN, the minimum value of LEVEL_MAX.
//
//    Input, int LEVEL_MAX_MAX, the maximum value of LEVEL_MAX.
//
{
  int dim;
  int dim_num;
  int growth;
  GWPointer2 *gw_compute_order;
  GWPointer *gw_compute_points;
  int i;
  int j;
  int level_max;
  int *np;
  int np_sum;
  double *p;
  int point_num;
  double tol;

  std::cout << "\n";
  std::cout << "SGMG_SIZE_TABULATE\n";
  std::cout << "  SGMG_SIZE returns the number of distinct\n";
  std::cout << "  points in a sparse grid.\n";
  std::cout << "\n";
  std::cout << "  We use the same rule in all dimensions, and count the points\n";
  std::cout << "  for a range of dimensions and levels.\n";
  std::cout << "\n";

  tol = std::sqrt ( webbur::r8_epsilon ( ) );

  std::cout << "   DIM: ";
  for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
  {
    std::cout << "  " << std::setw(8) << dim_num;
  }
  std::cout<< "\n";
  std::cout << "\n";
  std::cout << "   LEVEL_MAX\n";
  std::cout << "\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    std::cout << "    " << std::setw(4) << level_max;

    for ( dim_num = dim_min; dim_num <= dim_max; dim_num++ )
    {
      growth = growth_1d;
      NP = new int[dim_num];
      np_sum = dim_num * np_1d;
      P = new double[np_sum];
      gw_compute_points = new GWPointer[dim_num];
      gw_compute_order = new GWPointer2[dim_num];

      j = 0;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        NP[dim] = np_1d;
        for ( i = 0; i < np_sum; i++ )
        {
          P[j] = p_1d[i];
        }
        gw_compute_points[dim] = gw_compute_points_1d;
        gw_compute_order[dim] = gw_compute_order_1d;
      }

      point_num = webbur:: sgmg_size ( dim_num, level_max, gw_compute_points, 
        tol, growth, gw_compute_order );

      std::cout << "  " << std::setw(6) << point_num;

      delete [] gw_compute_order;
      delete [] gw_compute_points;
      delete [] NP;
      delete [] P;
    }
    std::cout << "\n";
  }
  return;
}
