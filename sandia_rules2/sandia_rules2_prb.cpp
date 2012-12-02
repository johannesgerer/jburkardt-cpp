# include "sandia_rules.hpp"
# include "sandia_rules2.hpp"

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
void test185 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SANDIA_RULES2_PRB.
//
//  Discussion:
//
//    SANDIA_RULES2_PRB calls a set of tests for the SANDIA_RULES2 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 April 2010
//
//  Author:
//
//    John Burkardt
//
{
  int r;

  webbur::timestamp ( );

  std::cout << "\n";
  std::cout << "SANDIA_RULES2_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the SANDIA_RULES2 library.\n";

  test185 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SANDIA_RULES2_PRB\n";
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

void test185 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST185 tests JACOBI_POINTS and JACOBI_WEIGHTS.
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
{
# define TEST_NUM 3

  double alpha;
  double alpha_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.5 };
  int dim;
  int i;
  int order;
  int order_max = 10;
  int test1;
  int test2;
  double w_diff;
  double *w1;
  double *w2;
  double x_diff;
  double *x1;
  double *x2;

  std::cout << "\n";
  std::cout << "TEST185\n";
  std::cout << "  JACOBI_POINTS and JACOBI_WEIGHTS compute a Gauss-Jacobi rule\n";
  std::cout << "  which is appropriate for integrands of the form\n";
  std::cout << "\n";
  std::cout << "    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.\n";
  std::cout << "\n";
  std::cout << "  For technical reasons, the parameters ALPHA and BETA are to\n";
  std::cout << "  supplied by a function called PARAMETER.\n";
  std::cout << "\n";
  std::cout << "   Order       ALPHA        BETA   ||X1-X2||   ||W1-W2||\n";
  std::cout << "\n";

  for ( test1 = 0; test1 < TEST_NUM; test1++ )
  {
    alpha = alpha_test[test1];

    for ( test2 = 0; test2 < TEST_NUM; test2++ )
    {
      beta = beta_test[test2];

      dim = 0;
      NP = new int[1];
      NP[0] = 2;
      P = new double[2];
      P[0] = alpha;
      P[1] = beta;

      for ( order = 1; order <= order_max; order++ )
      {
        w1 = new double[order];
        w2 = new double[order];
        x1 = new double[order];
        x2 = new double[order];

        webbur::jacobi_compute ( order, alpha, beta, x1, w1 );

        webbur::jacobi_points ( order, dim, x2 );
        webbur::jacobi_weights ( order, dim, w2 );

        x_diff = webbur::r8vec_diff_norm_li ( order, x1, x2 );
        w_diff = webbur::r8vec_diff_norm_li ( order, w1, w2 );

        std::cout << "  " << std::setw(6) << order
                  << "  " << std::setw(10) << alpha
                  << "  " << std::setw(10) << beta
                  << "  " << std::setw(10) << x_diff
                  << "  " << std::setw(10) << w_diff << "\n";
        delete [] w1;
        delete [] w2;
        delete [] x1;
        delete [] x2;
      }
      delete [] NP;
      delete [] P;
    }
  }
  return;
# undef TEST_NUM
}
