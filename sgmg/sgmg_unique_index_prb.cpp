# include "sandia_rules.hpp"
# include "sgmg.hpp"

# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>

int main ( );

void sgmg_unique_index_tests ( double tol );
void sgmg_unique_index_test ( int dim_num, int level_max_min, 
  int level_max_max, int rule[], int growth[], int np[], double p[],
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double w[] ),
  double tol );

typedef void ( *GWPointer ) ( int order, int np, double p[], double w[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SGMG_UNIQUE_INDEX_PRB.
//
//  Discussion:
//
//    SGMG_UNIQUE_INDEX_PRB tests SGMG_UNIQUE_INDEX.
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
//    20 June 2010
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
  int *growth;
  GWPointer *gw_compute_points;
  int level_max_max;
  int level_max_min;
  int *np;
  int np_sum;
  int *order_1d;
  int order_nd;
  double *p;
  int *rule;

  std::cout << "\n";
  std::cout << "SGMG_UNIQUE_INDEX_TESTS\n";
  std::cout << "  Call SGMG_UNIQUE_INDEX_TEST with various arguments.\n";
  std::cout << "\n";
  std::cout << "  All tests will use a point equality tolerance of " << tol << "\n";

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 1;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 6;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::clenshaw_curtis_compute_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 3;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 6;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::patterson_lookup_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 4;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 3;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::legendre_compute_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 7;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 3;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::laguerre_compute_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 1;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  p[0] = 1.5;
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 8;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 3;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::gen_laguerre_compute_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 2;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  p[0] = 0.5;
  p[1] = 1.5;
  rule = new int[dim_num];
  rule[0] = 2;
  rule[1] = 9;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 3;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::fejer2_compute_points_np;
  gw_compute_points[1] = webbur::jacobi_compute_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  level_max_min = 0;
  level_max_max = 2;
  np = new int[dim_num];
  np[0] = 1;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  p[0] = 2.0;
  rule = new int[dim_num];
  rule[0] = 6;
  rule[1] = 10;
  growth = new int[dim_num];
  growth[0] = 3;
  growth[1] = 4;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::gen_hermite_compute_points_np;
  gw_compute_points[1] = webbur::hermite_genz_keister_lookup_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 3;
  level_max_min = 0;
  level_max_max = 2;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 4;
  rule[2] = 5;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 3;
  growth[2] = 3;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::legendre_compute_points_np;
  gw_compute_points[2] = webbur::hermite_compute_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;
//
//  Repeat, treating  rules #2 and #3 as Golub Welsch rules.
//
  dim_num = 3;
  level_max_min = 0;
  level_max_max = 2;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 11;
  rule[2] = 11;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 3;
  growth[2] = 3;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::legendre_compute_points_np;
  gw_compute_points[2] = webbur::hermite_compute_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 3, exponential growth.
//
  dim_num = 2;
  level_max_min = 4;
  level_max_max = 4;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 3;
  rule[1] = 3;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 6;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::patterson_lookup_points_np;
  gw_compute_points[1] = webbur::patterson_lookup_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 3, slow exponential growth.
//
  dim_num = 2;
  level_max_min = 4;
  level_max_max = 4;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 3;
  rule[1] = 3;
  growth = new int[dim_num];
  growth[0] = 4;
  growth[1] = 4;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::patterson_lookup_points_np;
  gw_compute_points[1] = webbur::patterson_lookup_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;
//
//  Dimension 2, Level 4, Rule 3, moderate exponential growth.
//
  dim_num = 2;
  level_max_min = 4;
  level_max_max = 4;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 3;
  rule[1] = 3;
  growth = new int[dim_num];
  growth[0] = 5;
  growth[1] = 5;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::patterson_lookup_points_np;
  gw_compute_points[1] = webbur::patterson_lookup_points_np;
  sgmg_unique_index_test ( dim_num, level_max_min, level_max_max, 
    rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] np;
  delete [] p;
  delete [] rule;

  return;
}
//***************************************************************************80

void sgmg_unique_index_test ( int dim_num, int level_max_min, 
  int level_max_max, int rule[], int growth[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double w[] ),
  double tol )

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
//    20 June 2010
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
//    Input, int GROWTH[DIM_NUM], the growth rule in each dimension. 
//    0, "DF", default growth associated with this quadrature rule;
//    1, "SL", slow linear, L+1;
//    2  "SO", slow linear odd, O=1+2((L+1)/2)
//    3, "ML", moderate linear, 2L+1;
//    4, "SE", slow exponential;
//    5, "ME", moderate exponential;
//    6, "FE", full exponential.
//
//    Input, int NP[RULE_NUM], the number of parameters used by each rule.
//
//    Input, double P[sum(NP[*])], the parameters needed by each rule.
//
//    Input, void ( *GW_COMPUTE_POINTS[] ) ( int order, int np, double p[], double x[] ),
//    an array of pointers to functions which return the 1D quadrature points 
//    associated with each spatial dimension for which a Golub Welsch rule 
//    is used.
//
//    Input, double TOL, a tolerance for point equality.
//
{
  double alpha;
  double beta;
  int dim;
  int i;
  int level_max;
  int level_min;
  int p_index;
  int point;
  int point_num;
  int point_total_num;
  int *sparse_unique_index;

  std::cout << "\n";
  std::cout << "SGMG_UNIQUE_INDEX_TEST\n";
  std::cout << "  SGMG_UNIQUE_INDEX returns a mapping between\n";
  std::cout << "  the nonunique and unique points in a sparse grid.\n";
  std::cout << "\n";
  std::cout << " Dimension      Rule  Growth rate      Parameters\n";
  std::cout << "\n";

  p_index = 0;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( rule[dim] == 1 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim] 
                << "  " << std::setw(8) << growth[dim]  << "\n";
    }
    else if ( rule[dim] == 2 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << growth[dim]  << "\n";
    }
    else if ( rule[dim] == 3 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << growth[dim]  << "\n";
    }
    else if ( rule[dim] == 4 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << growth[dim]  << "\n";
    }
    else if ( rule[dim] == 5 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << growth[dim]  << "\n";
    }
    else if ( rule[dim] == 6 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << growth[dim]
                << "  " << std::setw(14) << alpha << "\n";
    }
    else if ( rule[dim] == 7 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << growth[dim]  << "\n";
    }
    else if ( rule[dim] == 8 )
    {
      alpha = p[p_index];
      p_index = p_index + 1;
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << growth[dim]
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
                << "  " << std::setw(8) << growth[dim]
                << "  " << std::setw(14) << alpha 
                << "  " << std::setw(14) << beta << "\n";
    }
    else if ( rule[dim] == 10 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << growth[dim]  << "\n";
    }
    else if ( rule[dim] == 11 )
    {
      std::cout << "  " << std::setw(8) << dim
                << "  " << std::setw(8) << rule[dim]
                << "  " << std::setw(8) << growth[dim];
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
                << "  " << std::setw(8) << growth[dim];
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
      std::cerr << "SGMG_UNIQUE_INDEX_TEST - Fatal error!\n";
      std::cerr << "  Unexpected value of RULE = " << rule[dim] << "\n";
      std::exit ( 1 );
    }
  }

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    point_total_num = webbur::sgmg_size_total ( dim_num, 
      level_max, rule, growth );

    point_num = webbur::sgmg_size ( dim_num, level_max, 
      rule, np, p, gw_compute_points, tol, growth );

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

    webbur::sgmg_unique_index ( dim_num, level_max, rule, 
      np, p, gw_compute_points, tol, point_num, point_total_num, 
      growth, sparse_unique_index );

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

