# include "sandia_rules.hpp"
# include "sgmg.hpp"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

int main ( );
void sgmg_size_tabulate ( int rule_1d, int growth_1d,
  int np_1d, double p_1d[], 
  void gw_compute_points_1d ( int order, int np, double p[], double w[] ),
  int dim_min, int dim_max, int level_max_min, int level_max_max );

typedef void ( *GWPointer ) ( int order, int np, double p[], double w[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SGMG_SIZE_TABLE.
//
//  Discussion:
//
//    SGMG_SIZE_TABLE_PRB makes a point count table.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2010
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
  double rule_1d;

  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SGMG_SIZE_TABLE\n";
  std::cout << "  C++ version\n";
//
//  Clenshaw-Curtis (1), slow exponential growth (4).
//
  rule_1d = 1;
  growth_1d = 4;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::clenshaw_curtis_compute_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Clenshaw-Curtis (1), slow exponential growth (4).
//
  rule_1d = 1;
  growth_1d = 4;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::clenshaw_curtis_compute_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Clenshaw-Curtis Grid (1), exponential growth (6).
//
  rule_1d = 1;
  growth_1d = 6;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::clenshaw_curtis_compute_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Clenshaw-Curtis Grid (1), exponential growth (6).
//
  rule_1d = 1;
  growth_1d = 6;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::clenshaw_curtis_compute_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Clenshaw-Curtis Grid (1), exponential growth (6).
//
  rule_1d = 1;
  growth_1d = 6;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 100;
  dim_max = 100;
  level_max_min = 0;
  level_max_max = 2;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::clenshaw_curtis_compute_points_np, 
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss-Patterson (3), slow exponential growth (4).
//
  rule_1d = 3;
  growth_1d = 4;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::patterson_lookup_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss-Patterson (3), slow exponential growth (4).
//
  rule_1d = 3;
  growth_1d = 4;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::patterson_lookup_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss-Patterson (3), Moderate Exponential Growth (5).
//
  rule_1d = 3;
  growth_1d = 5;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::patterson_lookup_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss-Patterson (3), Moderate Exponential Growth (5).
//
  rule_1d = 3;
  growth_1d = 5;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::patterson_lookup_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss-Patterson (3), exponential growth (6).
//
  rule_1d = 3;
  growth_1d = 6;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::patterson_lookup_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss-Patterson (3), exponential growth (6).
//
  rule_1d = 3;
  growth_1d = 6;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 6;
  dim_max = 10;
  level_max_min = 0;
  level_max_max = 5;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::patterson_lookup_points_np,
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU Time = " << ctime << "\n";
  delete [] p_1d;
//
//  Gauss Legendre Grid (4), slow linear odd growth (2)
//
  rule_1d = 4;
  growth_1d = 2;
  np_1d = 0;
  p_1d = new double[np_1d];
  dim_min = 1;
  dim_max = 5;
  level_max_min = 0;
  level_max_max = 7;
  ctime = webbur::cpu_time ( );
  sgmg_size_tabulate ( rule_1d, growth_1d, np_1d, 
    p_1d, webbur::legendre_compute_points_np, 
    dim_min, dim_max, level_max_min, level_max_max );
  ctime = webbur::cpu_time ( ) - ctime;
  std::cout << "\n";
  std::cout << "  CPU_TIME = " << ctime << "\n";
  delete [] p_1d;
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
//***************************************************************************80

void sgmg_size_tabulate ( int rule_1d, int growth_1d,
  int np_1d, double p_1d[], 
  void gw_compute_points_1d ( int order, int np, double p[], double w[] ),
  int dim_min, int dim_max, int level_max_min, int level_max_max )

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
//    22 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int RULE_1D, the 1D rule.
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
//    Input, int GROWTH_1D, the 1D growth rule. 
//    0, "DF", default growth associated with this quadrature rule;
//    1, "SL", slow linear, L+1;
//    2  "SO", slow linear odd, O=1+2((L+1)/2)
//    3, "ML", moderate linear, 2L+1;
//    4, "SE", slow exponential;
//    5, "ME", moderate exponential;
//    6, "FE", full exponential.
//
//    Input, int NP_1D, the number of parameters in the 1D rule.
//
//    Input, double P_1D[NP_1D], the parameters for the 1D rule.
//
//    Input, void gw_compute_points_1d ( int order, int np, double p[], double w[] ),
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
  int *growth;
  GWPointer *gw_compute_points;
  int i;
  int j;
  int level_max;
  int *np;
  int np_sum;
  double *p;
  int point_num;
  int *rule;
  double tol;

  std::cout << "\n";
  std::cout << "SGMG_SIZE_TABULATE\n";
  std::cout << "  SGMG_SIZE returns the number of distinct\n";
  std::cout << "  points in a sparse grid.\n";
  std::cout << "\n";
  std::cout << "  We use the same rule in all dimensions, and count the points\n";
  std::cout << "  for a range of dimensions and levels.\n";
  std::cout << "\n";
  std::cout << "  1D rule index is " << rule_1d << "\n";
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
      rule = new int[dim_num];
      growth = new int[dim_num];
      np = new int[dim_num];
      np_sum = dim_num * np_1d;
      p = new double[np_sum];
      gw_compute_points = new GWPointer[dim_num];

      j = 0;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        rule[dim] = rule_1d;
        growth[dim] = growth_1d;
        np[dim] = np_1d;
        for ( i = 0; i < np_sum; i++ )
        {
          p[j] = p_1d[i];
        }
        gw_compute_points[dim] = gw_compute_points_1d;
      }

      point_num = webbur:: sgmg_size ( dim_num, level_max, rule, 
        np, p, gw_compute_points, tol, growth );

      std::cout << "  " << std::setw(6) << point_num;

      delete [] growth;
      delete [] gw_compute_points;
      delete [] np;
      delete [] p;
      delete [] rule;
    }
    std::cout << "\n";
  }
  return;
}
