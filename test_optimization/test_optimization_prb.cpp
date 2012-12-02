# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_optimization.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_OPTIMIZATION_PRB.
//
//  Discussion:
//
//    TEST_OPTIMIZATION_PRB calls the TEST_OPTIMIZATION tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TEST_OPTIMIZATION_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_OPTIMIZATION library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_OPTIMIZATION_PRB\n";
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
//    TEST01 simply prints the title of each problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 February 2012
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
  cout << "  For each problem, print the title.\n";
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "  Problem  Title\n";
  cout << "\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );
    cout << "  " << setw(7) << problem 
         << "  " << title << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 samples the function at 1,000 points and prints the minimum.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double *f;
  double f_min;
  int know;
  int m = 2;
  int n = 1000;
  int problem;
  int problem_num;
  int seed;
  string title;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For each problem, using dimension M = 2\n";
  cout << "  sample the function at N = 1000 points,\n";
  cout << "  and print the minimum and maximum.\n";

  seed = 123456789;
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "  Problem     Minimum  Sample Minimum  Sample Maximum\n";
  cout << "\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    know = 0;
    x = p00_sol ( problem, m, know );
    if ( know != 0 )
    {
      f = p00_f ( problem, m, 1, x );
      f_min = f[0];
      delete [] x;
      delete [] f;
    }

    a = new double[m];
    b = new double[m];
    p00_ab ( problem, m, a, b );
    x = r8col_uniform_new ( m, n, a, b, &seed );
    f = p00_f ( problem, m, n, x );
    if ( know != 0 )
    {
      cout << "  " << setw(7) << problem
           << "  " << setw(14) << f_min
           << "  " << setw(14) << r8vec_min ( n, f )
           << "  " << setw(14) << r8vec_max ( n, f ) << "\n";
    }
    else
    {
      cout << "  " << setw(7) << problem
           << "  " << "              "
           << "  " << setw(14) << r8vec_min ( n, f )
           << "  " << setw(14) << r8vec_max ( n, f ) << "\n";
    }
    delete [] a;
    delete [] b;
    delete [] x;
    delete [] f;
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tries Compass Search on each problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double delta_init;
  double delta_tol;
  double *f;
  double *fv;
  double fx;
  int i;
  int k;
  int k_max;
  int know;
  int m = 2;
  int n = 1000;
  int problem;
  int problem_num;
  double s;
  int seed;
  string title;
  double *x;
  double *x0;

  delta_tol = 0.000001;
  k_max = 20000;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For each problem, using dimension M = 2\n";
  cout << "  try compass search.\n";
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    seed = 123456789;

    a = new double[m];
    b = new double[m];
    p00_ab ( problem, m, a, b );
    x0 = r8col_uniform_new ( m, 1, a, b, &seed );
    fv = p00_f ( problem, m, 1, x0 );
    s = 0.0;
    for ( i = 0; i < m; i++ )
    {
      s = s + x0[i] * x0[i];
    }
    delta_init = 0.3 * sqrt ( s ) / ( double ) ( m );
    delta_init = r8_max ( delta_init, 1000.0 * delta_tol );
    cout << "\n";
    cout << "  Problem " << setw(2) << problem
         << "  DELTA_INIT = " << setw(14) << delta_init << "\n";
    cout << "  Initial:" << setw(14) << x0[0] 
         << "  " << setw(14) << x0[1]
         << "  " << setw(14) << fv[0] << "\n";
    x = p00_compass_search ( problem, m, x0, delta_tol, delta_init, 
      k_max, &fx, &k );
    cout << "  Final:  " << setw(14) << x[0] 
         << "  " << setw(14) << x[1]
         << "  " << setw(14) << fx 
         << "  Steps = " << k << "\n";
    delete [] fv;
    delete [] x;

    know = 0;
    for ( ; ; )
    {
      x = p00_sol ( problem, m, know );
      if ( know == 0 )
      {
        break;
      }
      fv = p00_f ( problem, m, 1, x );
      cout << "  Exact:"
           << "  " << setw(14) << x[0]
           << "  " << setw(14) << x[1]
           << "  " << setw(14) << fv[0] << "\n";
      delete [] fv;
      delete [] x;
    }

    delete [] a;
    delete [] b;
    delete [] x0;
  }
  return;
}