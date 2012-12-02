# include "sandia_rules.hpp"
# include "sgmga.hpp"

# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>

int main ( );

void sgmga_vcn_coef_tests ( );
void sgmga_vcn_coef_test ( int dim_num, double importance[], 
  double level_weight[], int level_max_min, int level_max_max );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SGMGA_VCN_COEF_PRB.
//
//  Discussion:
//
//    SGMGA_VCN_COEF_PRB tests the SGMGA_VCN_COEF function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2009
//
//  Author:
//
//    John Burkardt
//
{
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SGMGA_VCN_COEF_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the SGMGA_VCN_COEF function.\n";
//
//  Compute examples of the combinatorial coefficent.
//
  sgmga_vcn_coef_tests ( );
//
//  That's all.
//
  std::cout << "\n";
  std::cout << "SGMGA_VCN_COEF_PRB\n";
  std::cout << "  Normal end of execution.\n";

  std::cout << "\n";
  webbur::timestamp ( );
  
  return 0;
}
//****************************************************************************80

void sgmga_vcn_coef_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_VCN_COEF_TESTS calls SGMGA_VCN_COEF_TEST.
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
  int dim;
  int dim_num;
  double *importance;
  int level_max;
  int level_max_max;
  int level_max_min;
  double *level_weight;

  std::cout << "\n";
  std::cout << "SGMGA_VCN_COEF_TESTS\n";
  std::cout << "  calls SGMGA_VCN_COEF_TEST.\n";

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 4;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  delete [] importance;
  delete [] level_weight;

  dim_num = 2;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 4;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  delete [] importance;
  delete [] level_weight;

  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = 1.0;
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 4;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  delete [] importance;
  delete [] level_weight;

  dim_num = 3;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 4;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  delete [] importance;
  delete [] level_weight;

  dim_num = 4;
  importance = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    importance[dim] = ( double ) ( dim + 1 );
  }
  level_weight = new double[dim_num];
  webbur::sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max_min = 0;
  level_max_max = 3;
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  delete [] importance;
  delete [] level_weight;
//
//  Try a case with a dimension of "0 importance".
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
  sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, 
    level_max_max );
  delete [] importance;
  delete [] level_weight;

  return;
}
//****************************************************************************80

void sgmga_vcn_coef_test ( int dim_num, double importance[], 
  double level_weight[], int level_max_min, int level_max_max )

//****************************************************************************80
//
//  Purpose:
//
//    SGMGA_VCN_COEF_TEST tests SGMGA_VCN_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  double coef1;
  double coef1_sum;
  double coef2;
  double coef2_sum;
  int dim;
  int i;
  int *level_1d;
  int *level_1d_max;
  int *level_1d_min;
  int level_max;
  double level_weight_min_pos;
  double level_weight_norm;
  bool more_grids;
  double q;
  double q_max;
  double q_min;

  level_1d = new int[dim_num];
  level_1d_max = new int[dim_num];
  level_1d_min = new int[dim_num];

  std::cout << "\n";
  std::cout << "SGMGA_VCN_COEF_TEST\n";
  std::cout << "  For anisotropic problems, a \"combinatorial coefficent\"\n";
  std::cout << "  must be computed for each component product grid.\n";
  std::cout << "  SGMGA_VCN_COEF_NAIVE does this in a simple, inefficient way.\n";
  std::cout << "  SGMGA_VCN_COEF tries to be more efficient.\n";
  std::cout << "  Here, we simply compare COEF1 and COEF2, the same\n";
  std::cout << "  coefficient computed by the naive and efficient ways.\n";
  std::cout << "\n";
  std::cout << "  IMPORTANCE:\n";
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(14) << importance[dim];
  }
  std::cout << "\n";
  std::cout << "  LEVEL_WEIGHT:\n";
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(14) << level_weight[dim];
  }
  std::cout << "\n";

  for ( level_max = level_max_min; level_max <= level_max_max; level_max++ )
  {
    i = 0;
    coef1_sum = 0.0;
    coef2_sum = 0.0;
//
//  Initialization.
//
    level_weight_min_pos = webbur::r8vec_min_pos ( dim_num, level_weight );
    q_min = ( double ) ( level_max ) * level_weight_min_pos 
      - webbur::r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max ) * level_weight_min_pos;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level_1d_min[dim] = 0;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      if ( 0.0 < level_weight[dim] )
      {
        level_1d_max[dim] = webbur::r8_floor ( q_max / level_weight[dim] ) + 1;
        if ( q_max <= ( level_1d_max[dim] - 1 ) * level_weight[dim] )
        {
          level_1d_max[dim] = level_1d_max[dim] - 1;
        }
      }
      else
      {
        level_1d_max[dim] = 0;
      }
    }
    more_grids = false;


    std::cout << "\n";
    std::cout << "     I               Q       Coef1       Coef2   X\n";
    std::cout << "   MIN" << "  " << std::setw(14) << q_min
              << "                        ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(2) << level_1d_min[dim];
    }
    std::cout << "\n";
//
//  Seek all vectors LEVEL_1D which satisfy the constraint:
//
//    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
//      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT[I] * LEVEL_1D[I]
//      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
//
    for ( ; ; )
    {
      webbur::sgmga_vcn_ordered_naive ( dim_num, level_weight, level_1d_max, 
        level_1d, q_min, q_max, &more_grids );

      if ( !more_grids )
      {
        break;
      }
//
//  Compute the combinatorial coefficient.
//
      coef1 = webbur::sgmga_vcn_coef_naive ( dim_num, level_weight, level_1d_max,
        level_1d, q_min, q_max );

      coef2 = webbur::sgmga_vcn_coef ( dim_num, level_weight, level_1d, q_max );

      i = i + 1;

      q = 0.0;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        q = q + level_weight[dim] * ( double ) level_1d[dim];
      }

      coef1_sum = coef1_sum + coef1;
      coef2_sum = coef2_sum + coef2;

      std::cout << "  " << std::setw(4)  << i
                << "  " << std::setw(14) << q
                << "  " << std::setw(10) << coef1
                << "  " << std::setw(10) << coef2;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        std::cout << "  " << std::setw(2) << level_1d[dim];
      }
      std::cout << "\n";
    }
    std::cout << "   MAX" << "  " << std::setw(14) << q_max
              << "                        ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(2) << level_1d_max[dim];
    }
    std::cout << "\n";
    std::cout << "   SUM                " 
              << "  " << std::setw(10) << coef1_sum 
              << "  " << std::setw(10) << coef2_sum << "\n";
  }
  delete [] level_1d;
  delete [] level_1d_max;
  delete [] level_1d_min;

  return;
}

