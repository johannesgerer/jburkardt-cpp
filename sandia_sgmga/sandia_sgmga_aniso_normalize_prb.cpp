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
void sandia_sgmga_aniso_balance_tests ( );
void sandia_sgmga_aniso_balance_test ( double alpha_max, int dim_num, 
  double level_weight[] );
void sandia_sgmga_aniso_normalize_tests ( );
void sandia_sgmga_aniso_normalize_test ( int dim_num, double level_weight[] );
void sandia_sgmga_importance_to_aniso_tests ( );
void sandia_sgmga_importance_to_aniso_test ( int dim_num, double importance[], 
  double level_weight[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SANDIA_SGMGA_ANISO_NORMALIZE_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  webbur::timestamp ( );

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_ANISO_NORMALIZE_PRB:\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the SANDIA_SGMGA_ANISO_NORMALIZE and\n";
  std::cout << "  SANDIA_SGMGA_IMPORTANCE_TO_ANISO functions.\n";

  sandia_sgmga_aniso_balance_tests ( );

  sandia_sgmga_aniso_normalize_tests ( );

  sandia_sgmga_importance_to_aniso_tests ( );

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_ANISO_NORMALIZE_PRB:\n";
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

void sandia_sgmga_aniso_balance_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_ANISO_BALANCE_TESTS call SANDIA_SGMGA_ANISO_BALANCE_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 November 2009
//
//  Author:
//
//    John Burkardt
//
{
  double alpha_max;
  int dim;
  int dim_num;
  double *level_weight;
  int seed;
  int test;
  int test_num;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_ANISO_BALANCE_TESTS\n";
  std::cout << "  Call SANDIA_SGMGA_ANISO_BALANCE_TEST with various arguments.\n";

  test_num = 5;
  dim_num = 5;
  level_weight = new double[dim_num];

  alpha_max = 10.0;
  seed = 123456789;
  for ( test = 1; test <= test_num; test++ )
  {
    webbur::r8vec_uniform_01 ( dim_num, &seed, level_weight );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level_weight[dim] = 10.0 * level_weight[dim];
    }
    sandia_sgmga_aniso_balance_test ( alpha_max, dim_num, level_weight );
  }

  alpha_max = 5.0;
  seed = 123456789;
  for ( test = 1; test <= test_num; test++ )
  {
    webbur::r8vec_uniform_01 ( dim_num, &seed, level_weight );
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level_weight[dim] = 10.0 * level_weight[dim];
    }
    sandia_sgmga_aniso_balance_test ( alpha_max, dim_num, level_weight );
  }

  delete [] level_weight;

  return;
}
//****************************************************************************80

void sandia_sgmga_aniso_balance_test ( double alpha_max, int dim_num, 
  double level_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_ANISO_BALANCE_TEST calls SANDIA_SGMGA_ANISO_BALANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  double *level_weight2;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_ANISO_BALANCE_TEST\n";
  std::cout << "  ALPHA_MAX = " << alpha_max << "\n";
  std::cout << "  Input weight sum: "
    << webbur::r8vec_sum ( dim_num, level_weight ) << "\n";
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(12) << level_weight[dim];
  }
  std::cout << "\n";

  level_weight2 = webbur::sandia_sgmga_aniso_balance ( alpha_max, dim_num, 
    level_weight );

  std::cout << "  Output weight sum: "
    << webbur::r8vec_sum ( dim_num, level_weight2 ) << "\n";
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(12) << level_weight2[dim];
  }
  std::cout << "\n";

  delete [] level_weight2;

  return;
}
//****************************************************************************80

void sandia_sgmga_aniso_normalize_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_ANISO_NORMALIZE_TESTS call SANDIA_SGMGA_ANISO_NORMALIZE_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 November 2009
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  double *level_weight;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_ANISO_NORMALIZE_TESTS\n";
  std::cout << "  Call SANDIA_SGMGA_ANISO_NORMALIZE_TEST with various arguments.\n";

  dim_num = 2;
  level_weight = new double[dim_num];
  level_weight[0] = 1.0;
  level_weight[1] = 1.0;
  sandia_sgmga_aniso_normalize_test ( dim_num, level_weight );
  delete [] level_weight;

  dim_num = 2;
  level_weight = new double[dim_num];
  level_weight[0] = 10.0;
  level_weight[1] = 10.0;
  sandia_sgmga_aniso_normalize_test ( dim_num, level_weight );
  delete [] level_weight;

  dim_num = 2;
  level_weight = new double[dim_num];
  level_weight[0] = 10.0;
  level_weight[1] = 2.0;
  sandia_sgmga_aniso_normalize_test ( dim_num, level_weight );
  delete [] level_weight;

  dim_num = 2;
  level_weight = new double[dim_num];
  level_weight[0] = 1.0;
  level_weight[1] = 2.0;
  sandia_sgmga_aniso_normalize_test ( dim_num, level_weight );
  delete [] level_weight;

  dim_num = 3;
  level_weight = new double[dim_num];
  level_weight[0] = 1.0;
  level_weight[1] = 2.0;
  level_weight[2] = 3.0;
  sandia_sgmga_aniso_normalize_test ( dim_num, level_weight );
  delete [] level_weight;
//
//  Try a case in which one variable has 0 weight.
//
  dim_num = 3;
  level_weight = new double[dim_num];
  level_weight[0] = 2.0;
  level_weight[1] = 0.0;
  level_weight[2] = 1.5;
  sandia_sgmga_aniso_normalize_test ( dim_num, level_weight );
  delete [] level_weight;

  dim_num = 4;
  level_weight = new double[dim_num];
  level_weight[0] = 1.0;
  level_weight[1] = 2.0;
  level_weight[2] = 3.0;
  level_weight[3] = 4.0;
  sandia_sgmga_aniso_normalize_test ( dim_num, level_weight );
  delete [] level_weight;

  return;
}
//****************************************************************************80

void sandia_sgmga_aniso_normalize_test ( int dim_num, double level_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_ANISO_NORMALIZE_TEST calls SANDIA_SGMGA_ANISO_NORMALIZE.
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
  int option;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_ANISO_NORMALIZE_TEST\n";
  std::cout << "  Input weight sum: "
    << webbur::r8vec_sum ( dim_num, level_weight ) << "\n";
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(12) << level_weight[dim];
  }
  std::cout << "\n";

  for ( option = 0; option <= 2; option++ )
  {
    webbur::sandia_sgmga_aniso_normalize ( option, dim_num, level_weight );

    std::cout << "  For OPTION = " << option
      << "  Normalized weight sum: "
      << webbur::r8vec_sum ( dim_num, level_weight ) << "\n";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      std::cout << "  " << std::setw(12) << level_weight[dim];
    }
    std::cout << "\n";
  }

  return;
}
//****************************************************************************80

void sandia_sgmga_importance_to_aniso_tests ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_IMPORTANCE_TO_ANISO_TESTS call SANDIA_SGMGA_IMPORTANCE_TO_ANISO_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 September 2009
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  double *importance;
  double *level_weight;

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_IMPORTANCE_TO_ANISO_TESTS\n";
  std::cout << "  Call SANDIA_SGMGA_IMPORTANCE_TO_ANISO_TEST with various arguments.\n";

  dim_num = 2;
  importance = new double[dim_num];
  level_weight = new double[dim_num];
  importance[0] = 1.0;
  importance[1] = 1.0;
  sandia_sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  delete [] importance;
  delete [] level_weight;

  dim_num = 2;
  importance = new double[dim_num];
  level_weight = new double[dim_num];
  importance[0] = 10.0;
  importance[1] = 10.0;
  sandia_sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  delete [] importance;
  delete [] level_weight;

  dim_num = 2;
  importance = new double[dim_num];
  level_weight = new double[dim_num];
  importance[0] = 10.0;
  importance[1] = 2.0;
  sandia_sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  delete [] importance;
  delete [] level_weight;

  dim_num = 2;
  importance = new double[dim_num];
  level_weight = new double[dim_num];
  importance[0] = 1.0;
  importance[1] = 2.0;
  sandia_sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  delete [] importance;
  delete [] level_weight;

  dim_num = 3;
  importance = new double[dim_num];
  level_weight = new double[dim_num];
  importance[0] = 1.0;
  importance[1] = 2.0;
  importance[2] = 3.0;
  sandia_sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  delete [] importance;
  delete [] level_weight;
//
//  Try a case in which one variable has 0 importance.
//
  dim_num = 3;
  importance = new double[dim_num];
  level_weight = new double[dim_num];
  importance[0] = 2.0;
  importance[1] = 0.0;
  importance[2] = 1.5;
  sandia_sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  delete [] importance;
  delete [] level_weight;

  dim_num = 4;
  importance = new double[dim_num];
  level_weight = new double[dim_num];
  importance[0] = 1.0;
  importance[1] = 2.0;
  importance[2] = 3.0;
  importance[3] = 4.0;
  sandia_sgmga_importance_to_aniso_test ( dim_num, importance, level_weight );
  delete [] importance;
  delete [] level_weight;

  return;
}
//****************************************************************************80

void sandia_sgmga_importance_to_aniso_test ( int dim_num, double importance[], 
  double level_weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGA_IMPORTANCE_TO_ANISO_TEST calls SANDIA_SGMGA_IMPORTANCE_TO_ANISO.
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

  std::cout << "\n";
  std::cout << "SANDIA_SGMGA_IMPORTANCE_TO_ANISO_TEST\n";
  std::cout << "  Importances:\n";
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(12) << importance[dim];
  }
  std::cout << "\n";

  webbur::sandia_sgmga_importance_to_aniso ( dim_num, importance, level_weight );

  std::cout << "  Anisotropic coefficients:\n";
  for ( dim = 0; dim < dim_num; dim++ )
  {
    std::cout << "  " << std::setw(12) << level_weight[dim];
  }
  std::cout << "\n";

  return;
}
