# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "laguerre_test_int.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LAGUERRE_TEST_INT_PRB.
//
//  Discussion:
//
//    LAGUERRE_TEST_INT_PRB tests the LAGUERRE_TEST_INT library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LAGUERRE_TEST_INT_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LAGUERRE_TEST_INT library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LAGUERRE_TEST_INT_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests P00_PROBLEM_NUM and P00_TITLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  int problem;
  int problem_num;
  string title;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  P00_PROBLEM_NUM returns the number of problems.\n";
  cout << "  P00_TITLE returns the title of a problem.\n";

  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "  P00_PROBLEM_NUM: number of problems is " << problem_num << "\n";
  cout << "\n";
  cout << "   Problem       Title\n";
  cout << "\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "  " << setw(8) << problem
         << "  \"" << title << "\".\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests P00_ALPHA and P00_EXACT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2007
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double exact;
  int problem;
  int problem_num;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  P00_ALPHA returns the lower limit of integration.\n";
  cout << "  P00_EXACT returns the \"exact\" integral.\n";

  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "   Problem       ALPHA           EXACT\n";
  cout << "\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    alpha = p00_alpha ( problem );

    exact = p00_exact ( problem );

    cout << "  " << setw(8) << problem
         << "  " << setw(14) << alpha
         << "  " << setprecision(16) << setw(24) << exact << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests P00_GAUSS_LAGUERRE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  int order;
  int order_log;
  int problem;
  int problem_num;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  P00_GAUSS_LAGUERRE applies a Gauss-Laguerre rule\n";
  cout << "  to estimate an integral on [ALPHA,+oo).\n";

  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "                              Exact\n";
  cout << "   Problem     Order          Estimate        Error\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 1;

    cout << "\n";
    cout << "  " << setw(8) << problem
         << "  " << "        "
         << "  " << setprecision(6) << setw(14) << exact << "\n";

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_gauss_laguerre ( problem, order );

      error = r8_abs ( exact - estimate );

      cout << "  " << "        "
           << "  " << setw(8) << order
           << "  " << setprecision(6) << setw(14) << estimate
           << "  " << setprecision(6) << setw(14) << error << "\n";

      order = order * 2;
    }
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests P00_EXP_TRANSFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  int order;
  int order_log;
  int problem;
  int problem_num;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  P00_EXP_TRANSFORM applies an exponential transform\n";
  cout << "  to estimate an integral on [ALPHA,+oo)\n";
  cout << "  as a transformed integral on (0,exp(-ALPHA)],\n";
  cout << "  and applying a Gauss-Legendre rule.\n";

  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "                              Exact\n";
  cout << "   Problem     Order          Estimate        Error\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 1;

    cout << "\n";
    cout << "  " << setw(8) << problem
         << "  " << "        "
         << "  " << setprecision(6) << setw(14) << exact << "\n";

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_exp_transform ( problem, order );

      error = r8_abs ( exact - estimate );

      cout << "  " << "        "
           << "  " << setw(8) << order
           << "  " << setprecision(6) << setw(14) << estimate
           << "  " << setprecision(6) << setw(14) << error << "\n";

      order = order * 2;
    }
  }
  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests P00_RAT_TRANSFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  int order;
  int order_log;
  int problem;
  int problem_num;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  P00_RAT_TRANSFORM applies a rational transform\n";
  cout << "  to estimate an integral on [ALPHA,+oo)\n";
  cout << "  as a transformed integral on (0,1/(1+ALPHA)],\n";
  cout << "  and applying a Gauss-Legendre rule.\n";

  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "                              Exact\n";
  cout << "   Problem     Order          Estimate        Error\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 1;

    cout << "\n";
    cout << "  " << setw(8) << problem
         << "  " << "        "
         << "  " << setprecision(6) << setw(14) << exact << "\n";

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_rat_transform ( problem, order );

      error = r8_abs ( exact - estimate );

      cout << "  " << "        "
           << "  " << setw(8) << order
           << "  " << setprecision(6) << setw(14) << estimate
           << "  " << setprecision(6) << setw(14) << error << "\n";

      order = order * 2;
    }
  }
  return;
}
