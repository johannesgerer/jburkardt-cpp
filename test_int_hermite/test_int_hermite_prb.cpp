# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "test_int_hermite.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_INT_HERMITE_PRB.
//
//  Discussion:
//
//    TEST_INT_HERMITE_PRB demonstrates the use of the TEST_INT_HERMITE
//    integration test functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TEST_INT_HERMITE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_INT_HERMITE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_INT_HERMITE_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  exit ( 0 );
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
//    26 May 2009 
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
//    TEST02 tests P00_EXACT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
{
  double exact;
  int m;
  int problem;
  int problem_num;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  P00_EXACT returns the \"exact\" integral.\n";

  problem_num = p00_problem_num ( );

  m = 4;
  p06_param ( 'S', 'M', &m );

  cout << "\n";
  cout << "   Problem       EXACT\n";
  cout << "\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    cout << "  " << setw(8) << problem
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
//    TEST03 tests P00_GAUSS_HERMITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  int m;
  int order;
  int order_log;
  int problem;
  int problem_num;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  P00_GAUSS_HERMITE applies a Gauss-Hermite rule\n";
  cout << "  to estimate an integral on (-oo,+oo).\n";

  problem_num = p00_problem_num ( );

  m = 4;
  p06_param ( 'S', 'M', &m );

  cout << "\n";
  cout << 
    "   Problem     Order          Estimate        Exact          Error\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 1;

    cout << "\n";

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_gauss_hermite ( problem, order );

      error = r8_abs ( exact - estimate );

      cout << "  " << setw(8) << problem
           << "  " << setw(8) << order
           << "  " << setprecision(6) << setw(14) << estimate
           << "  " << setprecision(6) << setw(14) << exact
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
//    TEST04 tests P00_TURING.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  double h;
  int m;
  int n;
  int order_log;
  int problem;
  int problem_num;
  int test;
  double tol;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  P00_TURING applies a Turing procedure\n";
  cout << "  to estimate an integral on (-oo,+oo).\n";

  problem_num = p00_problem_num ( );

  m = 4;
  p06_param ( 'S', 'M', &m );

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      tol = 1.0E-4;
    }
    else if ( test == 2 )
    {
      tol = 1.0E-07;
    }
    cout << "\n";
    cout << "  Using a tolerance of TOL = " << tol << "\n";
    cout << "\n";
    cout << 
      "   Problem     Order          Estimate        Exact          Error\n";

    for ( problem = 1; problem <= problem_num; problem++ )
    {
      exact = p00_exact ( problem );

      h = 1.0;

      cout << "\n";

      for ( order_log = 0; order_log <= 6; order_log++ )
      {
        estimate = p00_turing ( problem, h, tol, &n );

        error = r8_abs ( exact - estimate );

        cout << "  " << setw(8) << problem
             << "  " << setw(10) << h
             << "  " << setw(8) << n
             << "  " << setprecision(6) << setw(14) << estimate
             << "  " << setprecision(6) << setw(14) << exact
             << "  " << setprecision(6) << setw(14) << error << "\n";

        h = h / 2.0;
      }
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
//    TEST05 tests P00_GAUSS_HERMITE against the polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 2009
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  int m;
  int order;
  int order_log;
  int problem;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  P00_GAUSS_HERMITE applies a Gauss-Hermite rule to\n";
  cout << "  estimate the integral x^m exp(-x*x) over (-oo,+oo).\n";

  problem = 6;

  cout << "\n";
  cout << "         M     Order      Estimate        Exact           Error\n";

  for ( m = 0; m <= 6; m++ )
  {
    p06_param ( 'S', 'M', &m );

    exact = p00_exact ( problem );

    cout << "\n";

    for ( order = 1; order <= 3 + ( m / 2 ); order++ )
    {
      estimate = p00_gauss_hermite ( problem, order );

      error = r8_abs ( exact - estimate );

      cout << "  " << setw(8) << m
           << "  " << setw(8) << order
           << "  " << setw(14) << estimate
           << "  " << setw(14) << exact
           << "  " << setw(14) << error << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests P00_MONTE_CARLO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double estimate;
  double exact;
  int m;
  int order;
  int order_log;
  int problem;
  int problem_num;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  P00_MONTE_CARLO uses a weighted form of the Monte Carlo method\n";
  cout << "  to estimate a Hermite integral on (-oo,+oo).\n";

  problem_num = p00_problem_num ( );

  m = 4;
  p06_param ( 'S', 'M', &m );

  cout << "\n";
  cout << 
    "   Problem     Order          Estimate        Exact          Error\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    exact = p00_exact ( problem );

    order = 128;

    cout << "\n";

    for ( order_log = 0; order_log <= 6; order_log++ )
    {
      estimate = p00_monte_carlo ( problem, order );

      error = r8_abs ( exact - estimate );

      cout << "  " << setw(8) << problem
           << "  " << setw(8) << order
           << "  " << setprecision(6) << setw(14) << estimate
           << "  " << setprecision(6) << setw(14) << exact
           << "  " << setprecision(6) << setw(14) << error << "\n";

      order = order * 4;
    }
  }
  return;
}
