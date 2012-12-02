# include "sandia_rules.hpp"
# include "sgmga.hpp"

# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>

int main ( );

void sgmga_unique_index_tests ( );
void sgmga_unique_index_test ( int dim_num, double importance[],
  double level_weight[], int level_max_min, int level_max_max, int rule[], 
  int growth[], int np[], double p[],
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double w[] ),
  double tol );

typedef void ( *GWPointer ) ( int order, int np, double p[], double w[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SGMGA_UNIQUE_INDEX_PRB.
//
//  Discussion:
//
//    SGMGA_UNIQUE_INDEX_PRB tests the SGMGA_UNIQUE_INDEX routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 November 2009
//
//  Author:
//
//    John Burkardt
//
{
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SGMGA_UNIQUE_INDEX_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the routines in the SGMGA library.\n";
//
//  Check that we can determine a UNIQUE_INDEX mapping, so that we can
//  generate all the points, and then select unique representatives, and
//  then match each point to its representative.
//
  sgmga_unique_index_tests ( );
//
//  That's all.
//
  std::cout << "\n";
  std::cout << "SGMGA_UNIQUE_INDEX_PRB\n";
  std::cout << "  Normal end of execution.\n";

  std::cout << "\n";
  webbur::timestamp ( );
  
  return 0;
}
//****************************************************************************80

void sgmga_unique_index_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_UNIQUE_INDEX_TESTS calls SGMGA_UNIQUE_INDEX_TEST.
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
  int *growth;
  GWPointer *gw_compute_points;
  double *importance;
  int level_max_max;
  int level_max_min;
  double *level_weight;
  int *np;
  int np_sum;
  int *order_1d;
  int order_nd;
  double *p;
  int *rule;
  double tol;

  std::cout << "\n";
  std::cout << "SGMGA_UNIQUE_INDEX_TESTS\n";
  std::cout << "  Call SGMGA_UNIQUE_INDEX_TEST with various arguments\n";
//
//  Set the point equality tolerance.
//
  tol = std::sqrt ( webbur::r8_epsilon ( ) );
  std::cout << "\n";
  std::cout << "  All tests will use a point equality tolerance of " << tol << "\n";

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  rule[1] = 1;
  rule[2] = 1;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 6;
  growth[2] = 6;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::clenshaw_curtis_compute_points_np;
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  rule[1] = 1;
  rule[2] = 1;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 6;
  growth[2] = 6;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::clenshaw_curtis_compute_points_np;
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  gw_compute_points[1] = webbur::legendre_compute_points_np;
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;
//
//  Repeat, treating  rules #2 and #3 as Golub Welsch rules.
//
  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
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
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;
//
//  Try a case which includes a dimension of "0 importance".
//
  dim_num = 3;
  importance = new double[dim_num];
  importance[0] = 1.0;
  importance[1] = 0.0;
  importance[2] = 1.0;
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 3;
  np = new int[dim_num];
  np[0] = 0;
  np[1] = 0;
  np[2] = 0;
  np_sum = webbur::i4vec_sum ( dim_num, np );
  p = new double[np_sum];
  rule = new int[dim_num];
  rule[0] = 1;
  rule[1] = 1;
  rule[2] = 1;
  growth = new int[dim_num];
  growth[0] = 6;
  growth[1] = 6;
  growth[2] = 6;
  gw_compute_points = new GWPointer[dim_num];
  gw_compute_points[0] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[1] = webbur::clenshaw_curtis_compute_points_np;
  gw_compute_points[2] = webbur::clenshaw_curtis_compute_points_np;
  sgmga_unique_index_test ( dim_num, importance, level_weight, level_max_min,
    level_max_max, rule, growth, np, p, gw_compute_points, tol );
  delete [] growth;
  delete [] gw_compute_points;
  delete [] importance;
  delete [] level_weight;
  delete [] np;
  delete [] p;
  delete [] rule;

  return;
}
//***************************************************************************80

void sgmga_unique_index_test ( int dim_num, double importance[],
  double level_weight[], int level_max_min, int level_max_max, int rule[], 
  int growth[], int np[], double p[], 
  void ( *gw_compute_points[] ) ( int order, int np, double p[], double w[] ),
  double tol )

//***************************************************************************80
//
//  Purpose:
//
//    SGMGA_UNIQUE_INDEX_TEST tests SGMGA_UNIQUE_INDEX.
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
//    Input, double IMPORTANCE[DIM_NUM], the importance for each dimension.
//
//    Input, double LEVEL_WEIGHT[DIM_NUM], the weights for each dimension.
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
//    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
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
  int p_index;
  int point;
  int point_num;
  int point_total_num;
  int *sparse_unique_index;

  std::cout << "\n";
  std::cout << "SGMGA_UNIQUE_INDEX_TEST\n";
  std::cout << "  SGMGA_UNIQUE_INDEX returns a mapping between\n";
  std::cout << "  the nonunique and unique points in a sparse grid.\n";
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
      std::cerr << "SGMGA_UNIQUE_INDEX_TEST - Fatal error!\n";
      std::cerr << "  Unexpected value of RULE = " << rule[dim] << "\n";
      std::exit ( 1 );
    }
  }

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    point_total_num = webbur::sgmga_size_total ( dim_num, level_weight,
      level_max, rule, growth );

    point_num = webbur::sgmga_size ( dim_num, level_weight, level_max,
      rule, np, p, gw_compute_points, tol, growth );

    std::cout << "\n";
    std::cout << " LEVEL_MAX POINT_NUM POINT_NUM\n";
    std::cout << "              Unique     Total\n";

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << level_max
              << "  " << std::setw(8) << point_num
              << "  " << std::setw(8) << point_total_num << "\n";

    sparse_unique_index = new int[point_total_num];

    webbur::sgmga_unique_index ( dim_num, level_weight, level_max, rule, 
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

