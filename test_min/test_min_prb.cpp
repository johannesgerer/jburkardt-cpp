# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "test_min.hpp"

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
//    MAIN is the main program for TEST_MIN_PRB.
//
//  Discussion:
//
//    TEST_MIN_PRB calls the TEST_MIN tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp (  );

  cout << "\n";
  cout << "TEST_MIN_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_MIN library.\n";

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
  cout << "TEST_MIN_PRB\n";
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
//    TEST01 prints the title of each problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  int problem_num;
  int problem;
  string title;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For each problem, print the title.\n";
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "   Problem Title\n";
  cout << "\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "  " << setw(2) << problem
         << "  "            << title << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 evaluates the objective function at each starting point.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  double f_sol;
  double f_start;
  int know;
  int problem_num;
  int problem;
  string title;
  double x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For each problem, evaluate the function\n";
  cout << "  at the starting point and the solution.\n";
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "\n";
    cout << "  Problem " << problem << "\n";
    cout << "  " << title << "\n";
    cout << "\n";

    x = p00_start ( problem );

    f_start = p00_f ( problem, x );

    cout << "    F(X_START) = " << f_start << "\n";

    p00_sol ( problem, &know, &x );

    if ( 0 < know )
    {
      f_sol = p00_f ( problem, x );
      cout << "    F(X_SOL) = " << f_sol << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 compares the exact and approximate first derivatives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  double f1;
  double f1_dif;
  int problem_num;
  int problem;
  string title;
  double x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For each problem, compare the exact and\n";
  cout << "  approximate gradients at the starting point.\n";
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "\n";
    cout << "  Problem " << problem << "\n";
    cout << "  " << title << "\n";

    x = p00_start ( problem );

    f1 = p00_f1 ( problem, x );

    f1_dif = p00_f1_dif ( problem, x );

    cout << "\n";
    cout << "  X\n";
    cout << "  " << x << "\n";
    cout << "  F'(X) (exact)\n";
    cout << "  " << f1 << "\n";
    cout << "  F'(X) (difference)\n";
    cout << "  " << f1_dif << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 compares the exact and approximate second derivatives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  double f2;
  double f2_dif;
  int problem_num;
  int problem;
  string title;
  double x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For each problem, compare the exact and\n";
  cout << "  approximate second derivatives at the starting point.\n";
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "\n";
    cout << "  Problem " << problem << "\n";
    cout << "  " << title << "\n";

    x = p00_start ( problem );

    cout << "\n";
    cout << "  X:\n";
    cout << "  " << x << "\n";

    f2 = p00_f2 ( problem, x );

    cout << "  F\"(X) (exact):\n";
    cout << "  " << f2 << "\n";

    f2_dif = p00_f2_dif ( problem, x );

    cout << "  F\"(X) (difference):\n";
    cout << "  " << f2_dif << "\n";
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 carries out a simple bisection method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double fa;
  double fb;
  double fc;
  double fd;
  double fe;
  int i;
  int max_step = 10;
  int problem_num;
  int problem;
  string title;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For each problem, take a few steps of \n";
  cout << "  the bisection method.\n";
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "\n";
    cout << "  Problem " << problem << "\n";
    cout << "  " << title << "\n";

    p00_interval ( problem, &a, &c );
    b = 0.5 * ( a + c );
    fa = p00_f ( problem, a );
    fc = p00_f ( problem, c );
    fb = p00_f ( problem, b );

    i = 0;
    cout << "\n";
    cout << "  " << i << "\n";
    cout << "  X:"
         << "  " << setw(10) << a
         << "  " << setw(10) << b
         << "  " << setw(10) << c << "\n";
    cout << "  F:"
         << "  " << setw(10) << fa
         << "  " << setw(10) << fb
         << "  " << setw(10) << fc << "\n";

    for ( i = 1; i <= max_step; i++ )
    {
      d = 0.5 * ( a + b );
      fd = p00_f ( problem, d );

      e = 0.5 * ( b + c );
      fe = p00_f ( problem, e );

      if ( fd <= fb )
      {
        c = b;
        fc = fb;
        b = d;
        fb = fd;
      }
      else if ( fe <= fb )
      {
        a = b;
        fa = fb;
        b = e;
        fb = fe;
      }
      else
      {
        a = d;
        fa = fd;
        c = e;
        fc = fe;
      }

      cout << "  " << i << "\n";
      cout << "  X:"
           << "  " << setw(10) << a
           << "  " << setw(10) << b
           << "  " << setw(10) << c << "\n";
      cout << "  F:"
           << "  " << setw(10) << fa
           << "  " << setw(10) << fb
           << "  " << setw(10) << fc << "\n";
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
//    TEST06 carries out a version of Brent's derivative-free minimizer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2009
//
//  Author:
//
//    John Burkardt
//
{
  double fa;
  double fb;
  double fx;
  int problem_num;
  int problem;
  string title;
  double tol = 0.000001;
  double x;
  double a;
  double b;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  For each problem, use Brent's method.\n";
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "\n";
    cout << "  Problem " << problem << "\n";
    cout << "  " << title << "\n";

    p00_interval ( problem, &a, &b );

    fa = p00_f ( problem, a );
    fb = p00_f ( problem, b );

    cout << "\n";
    cout << "  Initial interval [A,B]:\n";
    cout << "\n";
    cout << "   A,       B:"
         << "  " << setw(16) << a
         << "  " << "                "
         << "  " << setw(16) << b << "\n";
    cout << "  FA,      FB:"
         << "  " << setw(16) << fa
         << "  " << "                "
         << "  " << setw(16) << fb << "\n";

    x = p00_fmin ( &a, &b, problem, tol );

    fa = p00_f ( problem, a );
    fb = p00_f ( problem, b );
    fx = p00_f ( problem, x );

    cout << "\n";
    cout << "  Final interval [A,X*,B]:\n";
    cout << "\n";
    cout << "   A,  X*,  B:"
         << "  " << setw(16) << a
         << "  " << setw(16) << x
         << "  " << setw(16) << b << "\n";
    cout << "  FA, FX*, FB:"
         << "  " << setw(16) << fa
         << "  " << setw(16) << fx
         << "  " << setw(16) << fb << "\n";
  }

  return;
}
