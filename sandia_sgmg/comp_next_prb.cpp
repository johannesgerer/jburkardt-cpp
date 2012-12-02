# include "sandia_rules.hpp"
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
void comp_next_test ( int dim_num, int level_max );

typedef void ( *GWPointer ) ( int order, int dim, double w[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for COMP_NEXT_PRB.
//
//  Discussion:
//
//     COMP_NEXT_PRB tests COMP_NEXT.
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
  int dim_num;
  int dim_num_array[12] = {
    2, 2, 2, 2, 2, 
    3, 3, 3, 3, 3, 
    4, 4 };
  int level_max;
  int level_max_array[12] = {
    0, 1, 2, 3, 4, 
    0, 1, 2, 3, 4, 
    2, 3 };
  int test;
  int test_num = 12;

  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << " COMP_NEXT_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "\n";
  std::cout << "  Test COMP_NEXT.\n";
//
//  Check that COMP_NEXT generates compositions correctly.
//
  for ( test = 0; test < test_num; test++ )
  {
    dim_num = dim_num_array[test];
    level_max = level_max_array[test];
    comp_next_test ( dim_num, level_max );
  }
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << " COMP_NEXT_PRB\n";
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

void comp_next_test ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT_TEST tests COMP_NEXT, which computes 1D level vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum level.
//
{
  int dim;
  int h;
  int i;
  int level;
  int *level_1d;
  int level_min;
  bool more_grids;
  int t;

  level_1d = new int[dim_num];
  level_min = webbur::i4_max ( 0, level_max + 1 - dim_num );

  std::cout << "\n";
  std::cout << "COMP_NEXT_TEST\n";
  std::cout << "  COMP_NEXT generates, one at a time, vectors\n";
  std::cout << "  LEVEL_1D(1:DIM_NUM) whose components add up to LEVEL.\n";
  std::cout << "\n";
  std::cout << "  We call with:\n";
  std::cout << "  DIM_NUM = " << dim_num << "\n";
  std::cout << "  " << level_min << " = LEVEL_MIN <= LEVEL <= LEVEL_MAX = "
            << level_max << "\n";
  std::cout << "\n";
  std::cout << "     LEVEL     INDEX  LEVEL_1D Vector\n";
//
//  The outer loop generates values of LEVEL from LEVEL_MIN to LEVEL_MAX.
//
  for ( level = level_min; level <= level_max; level++ )
  {
    std::cout << "\n";
//
//  The inner loop generates vectors LEVEL_1D(1:DIM_NUM) whose components 
//  add up to LEVEL.
//
    more_grids = false;
    h = 0;
    t = 0;
    i = 0;

    for ( ; ; )
    {
      webbur::comp_next ( level, dim_num, level_1d, &more_grids, &h, &t );

      i = i + 1;
      std::cout << "  " << std::setw(8) << level
                << "  " << std::setw(8) << i;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        std::cout << "  " << std::setw(8) << level_1d[dim];
      }
      std::cout << "\n";

      if ( !more_grids )
      {
        break;
      }
    }
  }
  delete [] level_1d;

  return;
}
