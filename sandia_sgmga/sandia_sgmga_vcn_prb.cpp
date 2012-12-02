# include "sandia_rules.hpp"
# include "sandia_rules2.hpp"
# include "sandia_sgmga.hpp"

# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>
# include <ctime>
//
//  Two global variables needed to support the "parameter" function.
//
double *P;
int *NP;

int main ( );
void sandia_sgmga_vcn_tests ( );
void sandia_sgmga_vcn_test ( int dim_num, double importance[], 
  double level_weight[], double q_min, double q_max );
void sandia_sgmga_vcn_timing_tests ( );
void sandia_sgmga_vcn_timing_test ( int dim_num, double importance[], 
  double level_weight[], double q_min, double q_max );
void sandia_sgmga_vcn_ordered_tests ( );
void sandia_sgmga_vcn_ordered_test ( int dim_num, double importance[], 
  double level_weight[], double q_min, double q_max );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SANDIA_SGMGA_VCN_PRB.
//
//  Discussion:
//
//    SANDIA_SGMGA_VCN_PRB tests SANDIA_SGMGA_VCN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_VCN_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the SANDIA_SGMGA_VCN and SANDIA_SGMGA_VCN_ORDERED functions.\n";

  sandia_sgmga_vcn_tests ( );

  sandia_sgmga_vcn_timing_tests ( );

  sandia_sgmga_vcn_ordered_tests ( );

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_VCN_PRB\n";
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

void sandia_sgmga_vcn_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_VCN_TESTS calls SANDIA_SGMGA_VCN_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int dim_num;
  int dim_num_array[12] = {
    2, 2, 2, 2, 2, 
    3, 3, 3, 3, 3, 
    4, 4 };
  double *importance;
  int level_max;
  int level_max_array[12] = {
    0, 1, 2, 3, 4, 
    0, 1, 2, 3, 4, 
    2, 3 };
  double *level_weight;
  double q_max;
  double q_min;
  int test;
  int test_num = 12;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_VCN_TESTS\n";
  std::cout << "  calls SANDIA_SGMGA_VCN_TEST.\n";
//
//  Isotropic examples.
//
  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_array[test];
    importance = new double[dim_num];
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = 1.0;
    }
    level_weight = new double[dim_num];
    webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = level_max_array[test];
    q_min = ( double ) ( level_max ) - webbur::r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sandia_sgmga_vcn_test ( dim_num, importance, level_weight, q_min, q_max );

    delete [] importance;
    delete [] level_weight;
  }
//
//  Zero weight example.
//
  dim_num = 3;
  importance = new double[dim_num];
  importance[0] = 1.0;
  importance[1] = 0.0;
  importance[2] = 1.0;
  level_weight = new double[dim_num];
  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
  level_max = 2;
  q_min = ( double ) ( level_max ) - webbur::r8vec_sum ( dim_num, level_weight );
  q_max = ( double ) ( level_max );

  sandia_sgmga_vcn_test ( dim_num, importance, level_weight, q_min, q_max );

  delete [] importance;
  delete [] level_weight;
//
//  Anisotropic examples.
//
  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_array[test];
    importance = new double[dim_num];
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = ( double ) ( dim + 1 );
    }
    level_weight = new double[dim_num];
    webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = level_max_array[test];
    q_min = ( double ) ( level_max ) - webbur::r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sandia_sgmga_vcn_test ( dim_num, importance, level_weight, q_min, q_max );

    delete [] importance;
    delete [] level_weight;
  }

  return;
}
//****************************************************************************80

void sandia_sgmga_vcn_test ( int dim_num, double importance[], double level_weight[], 
  double q_min, double q_max )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_VCN_TEST tests SANDIA_SGMGA_VCN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2011
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int i;
  int *level_1d;
  int *level_1d_max;
  int *level_1d_min;
  double level_weight_norm;
  bool more_grids;
  double q;
  int test;

  level_1d = new int[dim_num];
  level_1d_max = new int[dim_num];
  level_1d_min = new int[dim_num];

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_VCN_TEST\n";
  std::cout << "  Consider vectors 0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),\n";
  std::cout << "  Set Q = sum ( LEVEL_WEIGHT(1:N) * LEVEL_1D(1:N) )\n";
  std::cout << "  Accept vectors for which Q_MIN < Q <= Q_MAX\n";
  std::cout << "  No particular order is imposed on the LEVEL_1D values.\n";
  std::cout << "  SANDIA_SGMGA_VCN_NAIVE uses a naive approach;\n";
  std::cout << "  SANDIA_SGMGA_VCN tries to be more efficient.\n";
  std::cout << "  Here, we just compare the results.\n";

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
  for ( dim = 0; dim < dim_num; dim++ )
  {
    level_1d_min[dim] = 0;
  }

  more_grids = false;

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
  std::cout << "\n";
  std::cout << "  SANDIA_SGMGA_VCN_NAIVE\n";
  std::cout << "     I               Q   X\n";
  std::cout << "   MIN" << "  " << std::setw(14) << q_min;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_min[dim];
  }
  std::cout << "\n";

  i = 0;

  for ( ; ; )
  {
    webbur::sandia_sgmga_vcn_naive ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }

    q = 0.0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      q = q + level_weight[dim] * ( double ) level_1d[dim];
    }
    i = i + 1;
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(14) << q;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(2) << level_1d[dim];
    }
    std::cout << "\n";
  }
  std::cout << "   MAX" << "  " << std::setw(14) << q_max;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_max[dim];
  }
  std::cout << "\n";

  std::cout << "\n";
  std::cout << "  SANDIA_SGMGA_VCN\n";
  std::cout << "     I               Q   X\n";
  std::cout << "   MIN" << "  " << std::setw(14) << q_min;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_min[dim];
  }
  std::cout << "\n";

  i = 0;

  for ( ; ; )
  {
    webbur::sandia_sgmga_vcn ( dim_num, level_weight,level_1d, q_min, q_max,
      &more_grids );

    if ( !more_grids )
    {
      break;
    }

    q = 0.0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      q = q + level_weight[dim] * ( double ) level_1d[dim];
    }
    i = i + 1;
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(14) << q;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(2) << level_1d[dim];
    }
    std::cout << "\n";
  }
  std::cout << "   MAX" << "  " << std::setw(14) << q_max;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_max[dim];
  }
  std::cout << "\n";

  delete [] level_1d;
  delete [] level_1d_max;
  delete [] level_1d_min;

  return;
}
//****************************************************************************80

void sandia_sgmga_vcn_timing_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_VCN_TIMING_TESTS calls SANDIA_SGMGA_VCN_TIMING_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2010
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
  double *level_weight;
  double q_max;
  double q_min;
  int test;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_VCN_TIMING_TESTS\n";
  std::cout << "  calls SANDIA_SGMGA_VCN_TIMING_TEST.\n";
//
//  Isotropic examples.
//
  dim_num = 2;

  for ( test = 0; test < 2; test++ )
  {
    dim_num = dim_num * 2;
    importance = new double[dim_num];
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = 1.0;
    }
    level_weight = new double[dim_num];
    webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = 2;
    q_min = ( double ) ( level_max ) - webbur::r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sandia_sgmga_vcn_timing_test ( dim_num, importance, level_weight, q_min, q_max );

    delete [] importance;
    delete [] level_weight;
  }
//
//  Anisotropic examples.
//
  dim_num = 2;

  for ( test = 0; test < 2; test++ )
  {
    dim_num = dim_num * 2;
    importance = new double[dim_num];
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = ( double ) ( dim + 1 );
    }
    level_weight = new double[dim_num];
    webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = 2;
    q_min = ( double ) ( level_max ) - webbur::r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sandia_sgmga_vcn_timing_test ( dim_num, importance, level_weight, q_min, q_max );

    delete [] importance;
    delete [] level_weight;
  }

  return;
}
//****************************************************************************80

void sandia_sgmga_vcn_timing_test ( int dim_num, double importance[], 
  double level_weight[], double q_min, double q_max )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_VCN_TIMING_TEST times SANDIA_SGMGA_VCN and SANDIA_SGMGA_VCN_NAIVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int i;
  int *level_1d;
  int *level_1d_max;
  int *level_1d_min;
  double level_weight_norm;
  bool more_grids;
  double q;
  clock_t t1;
  clock_t t2;
  int test;

  level_1d = new int[dim_num];
  level_1d_max = new int[dim_num];
  level_1d_min = new int[dim_num];

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_VCN_TIMING_TEST\n";
  std::cout << "  Consider vectors 0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),\n";
  std::cout << "  Set Q = sum ( LEVEL_WEIGHT(1:N) * LEVEL_1D(1:N) )\n";
  std::cout << "  Accept vectors for which Q_MIN < Q <= Q_MAX\n";
  std::cout << "  No particular order is imposed on the LEVEL_1D values.\n";
  std::cout << "  SANDIA_SGMGA_VCN_NAIVE uses a naive approach;\n";
  std::cout << "  SANDIA_SGMGA_VCN tries to be more efficient.\n";
  std::cout << "  Here, we compare the timings.\n";

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
  for ( dim = 0; dim < dim_num; dim++ )
  {
    level_1d_min[dim] = 0;
  }

  more_grids = false;

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
  std::cout << "\n";
  std::cout << "  SANDIA_SGMGA_VCN_NAIVE\n";
  std::cout << "     I               Q   X\n";
  std::cout << "   MIN" << "  " << std::setw(14) << q_min;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_min[dim];
  }
  std::cout << "\n";

  i = 0;
  t1 = std::clock ( );
  for ( ; ; )
  {
    webbur::sandia_sgmga_vcn_naive ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }
  }
  t2 = std::clock ( );
  std::cout << "   MAX" << "  " << std::setw(14) << q_max;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_max[dim];
  }
  std::cout << "\n";
  std::cout << "  TIME" << "  " << std::setw(14)
            << ( double ) ( t2 - t1 ) / CLOCKS_PER_SEC << "\n";

  std::cout << "\n";
  std::cout << "  SANDIA_SGMGA_VCN\n";
  std::cout << "     I               Q   X\n";
  std::cout << "   MIN" << "  " << std::setw(14) << q_min;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_min[dim];
  }
  std::cout << "\n";

  i = 0;
  t1 = std::clock ( );
  for ( ; ; )
  {
    webbur::sandia_sgmga_vcn ( dim_num, level_weight,level_1d, q_min, q_max,
      &more_grids );

    if ( !more_grids )
    {
      break;
    }
  }
  t2 = std::clock ( );
  std::cout << "   MAX" << "  " << std::setw(14) << q_max;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_max[dim];
  }
  std::cout << "\n";
  std::cout << "  TIME" << "  " << std::setw(14) 
            << ( double ) ( t2 - t1 ) / CLOCKS_PER_SEC << "\n";

  delete [] level_1d;
  delete [] level_1d_max;
  delete [] level_1d_min;

  return;
}
//****************************************************************************80

void sandia_sgmga_vcn_ordered_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_VCN_ORDERED_TESTS calls SANDIA_SGMGA_VCN_ORDERED_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int dim_num;
  int dim_num_array[12] = {
    2, 2, 2, 2, 2, 
    3, 3, 3, 3, 3, 
    4, 4 };
  double *importance;
  int level_max;
  int level_max_array[12] = {
    0, 1, 2, 3, 4, 
    0, 1, 2, 3, 4, 
    2, 3 };
  double *level_weight;
  double q_max;
  double q_min;
  int test;
  int test_num = 12;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_VCN_ORDERED_TESTS\n";
  std::cout << "  calls SANDIA_SGMGA_VCN_ORDERED_TEST.\n";
//
//  Isotropic examples.
//
  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_array[test];
    importance = new double[dim_num];
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = 1.0;
    }
    level_weight = new double[dim_num];
    webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = level_max_array[test];
    q_min = ( double ) ( level_max ) - webbur::r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sandia_sgmga_vcn_ordered_test ( dim_num, importance, level_weight, q_min, q_max );

    delete [] importance;
    delete [] level_weight;
  }
//
//  Anisotropic examples.
//
  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_array[test];
    importance = new double[dim_num];
    for ( dim = 0; dim < dim_num; dim++ )
    {
      importance[dim] = ( double ) ( dim + 1 );
    }
    level_weight = new double[dim_num];
    webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );
    level_max = level_max_array[test];
    q_min = ( double ) ( level_max ) - webbur::r8vec_sum ( dim_num, level_weight );
    q_max = ( double ) ( level_max );

    sandia_sgmga_vcn_ordered_test ( dim_num, importance, level_weight, q_min, q_max );

    delete [] importance;
    delete [] level_weight;
  }

  return;
}
//****************************************************************************80

void sandia_sgmga_vcn_ordered_test ( int dim_num, double importance[], 
  double level_weight[], double q_min, double q_max )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_VCN_ORDERED_TEST tests SANDIA_SGMGA_VCN_ORDERED and SANDIA_SGMGA_VCN_ORDERED_NAIVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int i;
  int *level_1d;
  int *level_1d_max;
  int *level_1d_min;
  double level_weight_norm;
  bool more_grids;
  double q;
  int test;

  level_1d = new int[dim_num];
  level_1d_max = new int[dim_num];
  level_1d_min = new int[dim_num];

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_VCN_ORDERED_TEST\n";

  std::cout << "  Consider vectors 0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),\n";
  std::cout << "  Set Q = sum ( LEVEL_WEIGHT(1:N) * LEVEL_1D(1:N) )\n";
  std::cout << "  Accept only vectors for which Q_MIN < Q <= Q_MAX\n";
  std::cout << "  The solutions are weakly ordered by the value of Q.\n";
  std::cout << "  SANDIA_SGMGA_VCN_ORDERED_NAIVE calls SANDIA_SGMGA_VCN_NAIVE;\n";
  std::cout << "  SANDIA_SGMGA_VCN_ORDERED calls SANDIA_SGMGA_VCN.\n";

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
  for ( dim = 0; dim < dim_num; dim++ )
  {
    level_1d_min[dim] = 0;
  }

  more_grids = false;

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
  std::cout << "\n";
  std::cout << "  SANDIA_SGMGA_VCN_ORDERED_NAIVE:\n";
  std::cout << "     I               Q   X\n";
  std::cout << "   MIN" << "  " << std::setw(14) << q_min;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_min[dim];
  }
  std::cout << "\n";

  i = 0;

  for ( ; ; )
  {
    webbur::sandia_sgmga_vcn_ordered_naive ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }

    q = 0.0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      q = q + level_weight[dim] * ( double ) level_1d[dim];
    }
    i = i + 1;
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(14) << q;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(2) << level_1d[dim];
    }
    std::cout << "\n";
  }
  std::cout << "   MAX" << "  " << std::setw(14) << q_max;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_max[dim];
  }
  std::cout << "\n";

  std::cout << "\n";
  std::cout << "  SANDIA_SGMGA_VCN_ORDERED:\n";
  std::cout << "     I               Q   X\n";
  std::cout << "   MIN" << "  " << std::setw(14) << q_min;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_min[dim];
  }
  std::cout << "\n";

  i = 0;

  for ( ; ; )
  {
    webbur::sandia_sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, level_1d, 
      q_min, q_max, &more_grids );

    if ( !more_grids )
    {
      break;
    }

    q = 0.0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      q = q + level_weight[dim] * ( double ) level_1d[dim];
    }
    i = i + 1;
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(14) << q;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(2) << level_1d[dim];
    }
    std::cout << "\n";
  }
  std::cout << "   MAX" << "  " << std::setw(14) << q_max;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(2) << level_1d_max[dim];
  }
  std::cout << "\n";

  delete [] level_1d;
  delete [] level_1d_max;
  delete [] level_1d_min;

  return;
}


